#include "GrainsGPU.hh"
#include "ContactForceModelFactory.hh"
#include "Grains.hh"
#include "PostProcessingWriterFactory.hh"
#include "RigidBodyFactory.hh"
#include "TimeIntegratorFactory.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
GrainsGPU<T>::GrainsGPU()
{
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
GrainsGPU<T>::~GrainsGPU()
{
}

// -------------------------------------------------------------------------------------------------
// Setup the GPU and its parameters
template <typename T>
void GrainsGPU<T>::setupGPUDevice()
{
    using GP = GrainsParameters<T>;

    // Check available devices first
    int deviceCount = 0;
    cudaErrCheck(cudaGetDeviceCount(&deviceCount));
    GAssert(deviceCount > 0, "No CUDA devices found!");

    // Set the device to the first one
    uint device = 0;
    cudaErrCheck(cudaSetDevice(device));

    // Get the device properties
    cudaDeviceProp prop;
    cudaErrCheck(cudaGetDeviceProperties(&prop, device));
    GoutWI(3, "GPU Device:");
    GoutWI(6, "Device: ", prop.name);
    GoutWI(6, "Compute capability: ", prop.major, ".", prop.minor);
    GoutWI(6, "Number of SMs: ", prop.multiProcessorCount);
    GoutWI(6, "Total global memory: ", prop.totalGlobalMem / (1024 * 1024), " MB");
    GoutWI(6, "Max threads per block: ", prop.maxThreadsPerBlock);
    GoutWI(6,
           "Max threads dimensions: ",
           prop.maxThreadsDim[0],
           "x",
           prop.maxThreadsDim[1],
           "x",
           prop.maxThreadsDim[2]);

    // Since this function is invoked, it means a GPU simulation is requested.
    GP::m_isGPU = true;
    GP::m_GPU   = prop;
}

// -------------------------------------------------------------------------------------------------
// Initializes the simulation using the XML input
template <typename T>
void GrainsGPU<T>::initialize(DOMElement* rootElement)
{
    // Read using the base class Grains<T> with GPU context ready
    Grains<T>::initialize(rootElement);

    // Reading different blocks of the input XML
    // Note that most of the reading is done at Grains<T>::initialize
    // Herein, we initialize specifically for GPU
    Gout(std::string(80, '='));
    Gout("Setting up the simulation on device ...");
    Gout(std::string(80, '='));
    setupGPUDevice();
    Construction(rootElement);
    Forces(rootElement);
    AdditionalFeatures(rootElement);
    Gout(std::string(80, '='));
    Gout("Setting up the simulation on device completed!");
    Gout(std::string(80, '='));
}

// -------------------------------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
template <typename T>
void GrainsGPU<T>::simulate()
{
    using GP = GrainsParameters<T>;
    auto& SS = GP::m_simulationState;

    Gout(std::string(80, '='));
    Gout("Starting the simulation on GPU");
    Gout(std::string(80, '='));

    auto& timer = GP::m_simTimer;
    timer.start(SimStage::Total);

    // Insertion: particle insertion on host + copy to device
    timer.start(SimStage::Insertion);
    Grains<T>::m_components->insertParticles(Grains<T>::m_insertion);
    Grains<T>::m_components->updateSubBodyPositions();
    cout << "Copying the inserted particles to the device ..." << endl;
    Grains<T>::m_components->copyTo(m_d_components);
    cout << "Copying completed!" << endl;
    timer.stop(SimStage::Insertion);  // auto-syncs (IsGPU == true)

    cout << "\nTime \t TO \tend \tParticles \tIn \tOut" << endl;

    // Pre-compute forces on the initial configuration so the first moveParticles call (and the
    // KDK first half-kick) has valid accelerations.
    if(GP::m_isLeapFrog)
    {
        timer.start(SimStage::DetectCollisions);
        m_d_components->detectCollisions();
        timer.stop(SimStage::DetectCollisions);  // auto-syncs

        timer.start(SimStage::ComputeContactForces);
        m_d_components->computeContactForces(m_d_contactForce);
        timer.stop(SimStage::ComputeContactForces);  // auto-syncs
    }

    // Write initial state (t = tStart) before advancing
    SS.time = GP::m_tStart;
    timer.start(SimStage::PostProcess);
    Grains<T>::postProcess(Grains<T>::m_components);
    timer.stop(SimStage::PostProcess);

    uint stepCount = 0;
    for(SS.time = GP::m_tStart + GP::m_dt; SS.time <= GP::m_tEnd; SS.time += GP::m_dt)
    {
        stepCount++;
        // Output time
        if(GP::m_verbosityFrequency > 0 && (stepCount % GP::m_verbosityFrequency == 0))
        {
            ostringstream oss;
            oss.width(10);
            oss << left << SS.time;
            std::cout << oss.str() << "  \t" << GP::m_tEnd << std::endl;
        }

        if(GP::m_isLeapFrog)
        {
            // KDK Step 1: half-kick + drift using f_n (from pre-loop or previous step).
            timer.start(SimStage::MoveParticles);
            m_d_components->moveParticles(m_d_timeIntegrator);
            timer.stop(SimStage::MoveParticles);  // auto-syncs
            // Detect collisions and compute forces at x_{n+1}.
            timer.start(SimStage::DetectCollisions);
            m_d_components->detectCollisions();
            timer.stop(SimStage::DetectCollisions);  // auto-syncs
            timer.start(SimStage::ComputeContactForces);
            m_d_components->computeContactForces(m_d_contactForce);
            timer.stop(SimStage::ComputeContactForces);  // auto-syncs
            // KDK Step 3: second half-kick using f_{n+1}.
            timer.start(SimStage::AdvanceVelocity);
            m_d_components->advanceVelocity(m_d_timeIntegrator);
            timer.stop(SimStage::AdvanceVelocity);  // auto-syncs
        }
        else
        {
            // Single-pass scheme: compute forces at x_n, then advance.
            timer.start(SimStage::DetectCollisions);
            m_d_components->detectCollisions();
            timer.stop(SimStage::DetectCollisions);  // auto-syncs
            timer.start(SimStage::ComputeContactForces);
            m_d_components->computeContactForces(m_d_contactForce);
            timer.stop(SimStage::ComputeContactForces);  // auto-syncs
            timer.start(SimStage::MoveParticles);
            m_d_components->moveParticles(m_d_timeIntegrator);
            timer.stop(SimStage::MoveParticles);  // auto-syncs
        }

        // Post-Processing
        timer.start(SimStage::PostProcess);
        Grains<T>::postProcess(m_d_components);
        timer.stop(SimStage::PostProcess);
    }
    cudaDeviceSynchronize();

    timer.stop(SimStage::Total);
    if(timer.isEnabled())
    {
        if(GP::m_cdmTimer.isEnabled())
            GP::m_cdmTimer.printSummary();
        if(GP::m_fmTimer.isEnabled())
            GP::m_fmTimer.printSummary();
        timer.printSummary();
    }
}

/* ============================================================================================== */
/* Low-Level Methods                                                                              */
/* ============================================================================================== */
// Constructs the simulation -- Reads the Construction part of the XML input to
// set the parameters
template <typename T>
void GrainsGPU<T>::Construction(DOMElement* rootElement)
{
    using GP = GrainsParameters<T>;

    // ---------------------------------------------------------------------------------------------
    // Particles
    GoutWI(3, "Copying rigid bodies to device ...");
    RigidBodyFactory<T>::copyHostToDevice(Grains<T>::m_rigidBodyList, m_d_rigidBodyList);
    GoutWI(3, "Copying rigid bodies to device completed!");

    // ---------------------------------------------------------------------------------------------
    // Contact force models
    // It is a GPU simulation, and we have already read contact force models
    // on the host. We allocate memory on device and copy the models over.
    GoutWI(3, "Copying contact force models to device ...");
    m_d_contactForce.reserve(GP::m_numContactPairs);
    ContactForceModelFactory<T>::copyHostToDevice(Grains<T>::m_contactForce, m_d_contactForce);
    GoutWI(3, "Copying contact force models to device completed!");

    // ---------------------------------------------------------------------------------------------
    // Temporal setting and time integration
    // It is a GPU simulation, and we have already read time integration on the
    // host. We allocate memory on device and copy the scheme over.
    GoutWI(3, "Copying time integration scheme to device ...");
    TimeIntegratorFactory<T>::copyHostToDevice(Grains<T>::m_timeIntegrator, m_d_timeIntegrator);
    GoutWI(3, "Copying time integration scheme to device completed!");

    // ---------------------------------------------------------------------------------------------
    // Setting up the component managers
    // Collect all initial state from the CPU manager as HOST buffers.
    // The DEVICE constructor uploads them to device and creates CDM + ForceModule.
    auto&                                         hm = *Grains<T>::m_components;
    GrainsMemBuffer<uint, MemType::HOST>          hostBodyTags;
    GrainsMemBuffer<Vector3<T>, MemType::HOST>    hostPos;
    GrainsMemBuffer<Quaternion<T>, MemType::HOST> hostOri;
    GrainsMemBuffer<Vector3<T>, MemType::HOST>    hostLocalPos;
    GrainsMemBuffer<Quaternion<T>, MemType::HOST> hostLocalQuat;
    hm.getBodyTag(hostBodyTags);
    hm.getPosition(hostPos);
    hm.getQuaternion(hostOri);
    hm.getLocalPos(hostLocalPos);
    hm.getLocalQuat(hostLocalQuat);

    m_d_components
        = std::make_unique<ComponentManager<T, MemType::DEVICE>>(&m_d_rigidBodyList,
                                                                 std::move(hostBodyTags),
                                                                 std::move(hostPos),
                                                                 std::move(hostOri),
                                                                 std::move(hostLocalPos),
                                                                 std::move(hostLocalQuat),
                                                                 hm.getNumberOfObstacles(),
                                                                 hm.getNumberOfParticles(),
                                                                 hm.getNumberOfComposites(),
                                                                 hm.getNumberOfSubBodies());
}

// -------------------------------------------------------------------------------------------------
// External force definition
template <typename T>
void GrainsGPU<T>::Forces(DOMElement* rootElement)
{
}

// -------------------------------------------------------------------------------------------------
// Additional features of the simulation: insertion, post-processing
template <typename T>
void GrainsGPU<T>::AdditionalFeatures(DOMElement* rootElement)
{
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class GrainsGPU<float>;
template class GrainsGPU<double>;