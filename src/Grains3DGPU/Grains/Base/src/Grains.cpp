#include "Grains.hh"
#include "ComponentManager.hh"
#include "ContactForceModelFactory.hh"
#include "PostProcessingWriterFactory.hh"
#include "RigidBodyFactory.hh"
#include "TimeIntegratorFactory.hh"
#include "VectorMath.hh"

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
// Default constructor
template <typename T>
Grains<T>::Grains()
{
    Gout(std::string(80, '='));
    Gout("Starting Grains3D ...");
    Gout(std::string(80, '='));
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
Grains<T>::~Grains()
{
}

// -------------------------------------------------------------------------------------------------
// Initializes the simulation using the XML input
template <typename T>
void Grains<T>::initialize(DOMElement* rootElement)
{
    // Reading different blocks of the input XML
    Gout(std::string(80, '='));
    Gout("Reading the input file ...");
    Gout(std::string(80, '='));
    Construction(rootElement);
    Forces(rootElement);
    AdditionalFeatures(rootElement);
    Gout(std::string(80, '='));
    Gout("Reading the input file completed!");
    Gout(std::string(80, '='));

    // Post-processing start
    for(auto& pp : m_postProcessor)
        pp->PostProcessing_start();
}

// -------------------------------------------------------------------------------------------------
// Performs post-processing
template <typename T>
template <MemType MT>
void Grains<T>::postProcess(const std::unique_ptr<ComponentManager<T, MT>>& cm)
{
    using GP = GrainsParameters<T>;
    auto& SS = GP::m_simulationState;

    if(GP::m_tSave.empty())
        return;

    if(GP::m_tSave.front() - SS.time < 0.01 * GP::m_dt)
    {
        GP::m_tSave.pop();
        if constexpr(MT == MemType::DEVICE)
        {
            cm->copyTo_PostProcessing(m_components);
            for(auto& pp : m_postProcessor)
                pp->PostProcessing(m_rigidBodyList, m_components, SS.time);
        }
        else
        {
            for(auto& pp : m_postProcessor)
                pp->PostProcessing(m_rigidBodyList, cm, SS.time);
        }
    }
    // In case we get past the saveTime, we need to remove it from the queue
    if(!GP::m_tSave.empty() && SS.time > GP::m_tSave.front())
        GP::m_tSave.pop();
}

// -------------------------------------------------------------------------------------------------
// Performs tasks after time-stepping
template <typename T>
void Grains<T>::finalize()
{
    for(auto& pp : m_postProcessor)
        pp->PostProcessing_end();
}

/* ============================================================================================== */
/* Low-Level Methods                                                                              */
/* ============================================================================================== */
// Sets up rigid bodies, body tags, local transforms, and ComponentManager from XML
template <typename T>
void Grains<T>::setupComponents(DOMNode* root, DOMElement* rootElement)
{
    using GP = GrainsParameters<T>;
    auto& LC = GP::m_collisionDetection.linkedCellParameters;

    DOMNode*        particles = ReaderXML::getNode(root, "Particles");
    DOMNode*        obstacles = ReaderXML::getNode(root, "Obstacles");
    ParticleData<T> obstacleData, particleData;
    uint            numObstacles = 0, numParticles = 0;

    GoutWI(6, "Reading rigid bodies ...");
    RigidBodyFactory<T>::create(obstacles,
                                particles,
                                obstacleData,
                                particleData,
                                numObstacles,
                                numParticles);
    GoutWI(6, "Reading rigid bodies completed!");

    // Count composites and sub-bodies
    uint nComposites = 0;
    uint nSubBodies  = 0;
    for(uint i = 0; i < particleData.numTemplates; ++i)
    {
        if(particleData.isComposite[i])
        {
            nComposites += particleData.numEachRef[i];
            nSubBodies += particleData.numEachRef[i] * particleData.subBodiesPerTemplate[i];
        }
    }

    // Allocate rigid body list and per-component initial state buffers
    const uint totalNumComponents = numObstacles + numParticles;
    GAssert(totalNumComponents > 0, "No components found in the simulation!");
    m_rigidBodyList.initialize(totalNumComponents);
    GrainsMemBuffer<Vector3<T>>    initialPosition(totalNumComponents);
    GrainsMemBuffer<Quaternion<T>> initialOrientation(totalNumComponents);
    GrainsMemBuffer<uint>          initialBodyTags(totalNumComponents);
    GrainsMemBuffer<Vector3<T>>    initialLocalPos(totalNumComponents);
    GrainsMemBuffer<Quaternion<T>> initialLocalQuat(totalNumComponents);

    uint offset = 0;

    // Expand obstacles (subBodiesPerTemplate = 1, numEachRef = 1 always)
    for(uint i = 0; i < obstacleData.numTemplates; ++i)
    {
        m_rigidBodyList[offset]    = new RigidBody<T>(*obstacleData.refRB[i]);
        initialBodyTags[offset]    = makeStandaloneBodyTag(i);
        initialPosition[offset]    = obstacleData.refInitialPos[i];
        initialOrientation[offset] = obstacleData.refInitialOri[i];
        initialLocalPos[offset]    = Vector3<T>(T(0), T(0), T(0));
        initialLocalQuat[offset]   = Quaternion<T>(T(0), T(0), T(0), T(1));
        T r                        = obstacleData.refRB[i]->getCircumscribedRadius();
        if(r > LC.maxObstacleRadius)
            LC.maxObstacleRadius = r;
        ++offset;
    }

    // Expand particles
    T    maxParticleRadius = T(0);
    uint compositeIdx      = 1;  // 1-based
    uint protoBase         = 0;  // flat index into refRB/refLocalPos/refLocalQuat
    for(uint i = 0; i < particleData.numTemplates; ++i)
    {
        uint S_i = particleData.subBodiesPerTemplate[i];
        uint M_i = particleData.numEachRef[i];
        bool isc = (particleData.isComposite[i] != 0u);

        for(uint j = 0; j < M_i; ++j)
        {
            for(uint k = 0; k < S_i; ++k)
            {
                uint rbIdx                 = protoBase + k;
                m_rigidBodyList[offset]    = new RigidBody<T>(*particleData.refRB[rbIdx]);
                initialPosition[offset]    = particleData.refInitialPos[i];
                initialOrientation[offset] = particleData.refInitialOri[i];
                initialLocalPos[offset]    = particleData.refLocalPos[rbIdx];
                initialLocalQuat[offset]   = particleData.refLocalQuat[rbIdx];
                // shapeId = numObstacles + protoBase + k (prototype index, not slot index)
                // All instances of the same particle prototype share the same shapeId,
                // keeping the value well within the 10-bit field regardless of particle count.
                const uint sid = numObstacles + protoBase + k;
                GAssert(sid < 1024u, "shapeId overflow: too many unique particle prototypes");
                initialBodyTags[offset]
                    = isc ? makeSubBodyTag(sid, compositeIdx, k) : makeStandaloneBodyTag(sid);
                T r = particleData.refRB[rbIdx]->getCircumscribedRadius();
                if(r > maxParticleRadius)
                    maxParticleRadius = r;
                ++offset;
            }
            if(isc)
                ++compositeIdx;
        }
        protoBase += S_i;
    }

    // Store radii in LC params for use in Construction()
    LC.minCellSize = T(2) * maxParticleRadius;

    // Compute maxNumCellsPerObstacle now that minCellSize, maxObstacleRadius, and
    // cellSizeFactor are all known (cellSizeFactor set by CD XML parsing before this call)
    if(LC.minCellSize > T(0))
    {
        const uint adj            = (LC.maxObstacleRadius > T(0))
                                        ? static_cast<uint>(std::ceil((T(2) * LC.maxObstacleRadius)
                                                           / (LC.cellSizeFactor * LC.minCellSize)))
                                        : 1u;
        LC.maxNumCellsPerObstacle = adj * adj * adj;
    }

    // Free prototype rigid bodies
    for(uint i = 0; i < obstacleData.refRB.getSize(); ++i)
    {
        delete obstacleData.refRB[i];
        obstacleData.refRB[i] = nullptr;
    }
    for(uint i = 0; i < particleData.refRB.getSize(); ++i)
    {
        delete particleData.refRB[i];
        particleData.refRB[i] = nullptr;
    }

    // Contact Force Models
    // NOTE: Read after rigid bodies so that m_materialMap (populated by RigidBodyFactory)
    // is available to compute m_numContactPairs and resolve material IDs.
    {
        uint numMaterials     = GP::m_materialMap.size();
        GP::m_numContactPairs = numMaterials * (numMaterials + 1) / 2;
        DOMNode* contacts     = ReaderXML::getNode(root, "ContactForceModels");
        if(contacts)
        {
            GoutWI(6, "Reading contact force models ...");
            ContactForceModelFactory<T>::create(rootElement, m_contactForce);
            GoutWI(6, "Reading contact force models completed!");
            if(ReaderXML::getNodeAttr_String(contacts, "EnableTimings") == "true")
                GP::m_fmTimer.enable(GP::m_isGPU);
            if(ReaderXML::getNodeAttr_String(contacts, "Compaction") == "false")
            {
                GP::m_useCompaction = false;
                GoutWI(9, "Compaction disabled!");
            }
        }
    }

    // Construct ComponentManager. Contact forces MUST be parsed before this call.
    m_components
        = std::make_unique<ComponentManager<T, MemType::HOST>>(&m_rigidBodyList,
                                                               std::move(initialBodyTags),
                                                               std::move(initialPosition),
                                                               std::move(initialOrientation),
                                                               std::move(initialLocalPos),
                                                               std::move(initialLocalQuat),
                                                               numObstacles,
                                                               numParticles,
                                                               nComposites,
                                                               nSubBodies);
}

// -------------------------------------------------------------------------------------------------
// Constructs the simulation -- Reads the Construction part of the XML input to
// set the parameters
template <typename T>
void Grains<T>::Construction(DOMElement* rootElement)
{
    using GP = GrainsParameters<T>;
    auto& CD = GP::m_collisionDetection;
    auto& LC = CD.linkedCellParameters;

    // Output message
    GoutWI(3, "Construction");
    // ---------------------------------------------------------------------------------------------
    // Checking if Construction node is available
    DOMNode* root = ReaderXML::getNode(rootElement, "Construction");
    GAssert(root, "Construction node is mandatory!");

    // ---------------------------------------------------------------------------------------------
    // Domain size: origin, max coordinates and periodicity
    DOMNode* nOrigin = ReaderXML::getNode(root, "Origin");
    if(nOrigin)
        GP::m_origin.setValue(T(ReaderXML::getNodeAttr_Double(nOrigin, "X")),
                              T(ReaderXML::getNodeAttr_Double(nOrigin, "Y")),
                              T(ReaderXML::getNodeAttr_Double(nOrigin, "Z")));
    else
        GP::m_origin.setValue(T(0), T(0), T(0));

    DOMNode* nDomain = ReaderXML::getNode(root, "MaxCoordinate");
    GP::m_maxCoordinate.setValue(T(ReaderXML::getNodeAttr_Double(nDomain, "X")),
                                 T(ReaderXML::getNodeAttr_Double(nDomain, "Y")),
                                 T(ReaderXML::getNodeAttr_Double(nDomain, "Z")));

    // if the simulation is periodic
    DOMNode* nPeriodicity = ReaderXML::getNode(root, "Periodicity");
    if(nPeriodicity)
    {
        int PX = ReaderXML::getNodeAttr_Int(nPeriodicity, "PX");
        int PY = ReaderXML::getNodeAttr_Int(nPeriodicity, "PY");
        int PZ = ReaderXML::getNodeAttr_Int(nPeriodicity, "PZ");
        GAssert(PX * PY * PZ == 0, "Periodicity is not implemented!");
        GP::m_isPeriodic = false;
    }

    // ---------------------------------------------------------------------------------------------
    // Setting up collision detection
    GoutWI(6, "Reading collision detection ...");
    LC.minCorner                = GP::m_origin;
    LC.maxCorner                = GP::m_maxCoordinate;
    DOMNode* collisionDetection = ReaderXML::getNode(root, "CollisionDetection");
    GAssert(collisionDetection, "CollisionDetection node is mandatory!");
    if(ReaderXML::getNodeAttr_String(collisionDetection, "EnableTimings") == "true")
    {
        GoutWI(9, "Collision detection timings enabled!");
        GP::m_cdmTimer.enable(GP::m_isGPU);
    }
    {
        const std::string relTransVal
            = ReaderXML::getNodeAttr_String(collisionDetection, "UseRelativeTransformations");
        if(relTransVal == "true")
        {
            GoutWI(9, "Collision detection using relative transformations enabled!");
            CD.useRelativeTransformations = true;
        }
        else if(relTransVal == "false")
        {
            GoutWI(9, "Collision detection using relative transformations disabled!");
            CD.useRelativeTransformations = false;
        }
    }

    // Neighbor list
    DOMNode* nNeighborList = ReaderXML::getNode(collisionDetection, "NeighborList");
    GAssert(nNeighborList, "NeighborList node is mandatory!");
    std::string neighborListType = ReaderXML::getNodeAttr_String(nNeighborList, "Type");
    if(neighborListType == "BruteForce")
        CD.neighborListType = NeighborListType::NSQ;
    else if(neighborListType == "LinkedCell")
        CD.neighborListType = NeighborListType::LINKEDCELL;
    else
        GAbort("Unknown NeighborList type! Aborting Grains!");
    GoutWI(9, "NeighborList: " + neighborListType);

    // Linked cell
    if(CD.neighborListType == NeighborListType::LINKEDCELL)
    {
        DOMNode* nLinkedCell = ReaderXML::getNode(collisionDetection, "LinkedCell");
        GAssert(nLinkedCell, "LinkedCell node is mandatory when using LinkedCell neighbor list!");
        std::string linkedCellType = ReaderXML::getNodeAttr_String(nLinkedCell, "Type");
        if(linkedCellType == "Host")
            LC.type = LinkedCellType::HOST;
        else if(linkedCellType == "Device_SortBased")
            LC.type = LinkedCellType::SORTBASED;
        else if(linkedCellType == "Device_Atomic")
            LC.type = LinkedCellType::ATOMIC;
        else if(linkedCellType == "Device_AtomicFixed")
            LC.type = LinkedCellType::ATOMICFIXED;
        else
            GAbort("Unknown LinkedCell type! Aborting Grains!");

        // Cell size factor
        LC.cellSizeFactor = T(ReaderXML::getNodeAttr_Double(nLinkedCell, "CellSizeFactor"));

        // Update and sort frequency
        LC.updateFrequency = ReaderXML::getNodeAttr_Int(nLinkedCell, "UpdatingFrequency");
        LC.sortFrequency   = ReaderXML::getNodeAttr_Int(nLinkedCell, "SortingFrequency");

        GoutWI(9,
               "LinkedCell: " + linkedCellType + +", cell size factor "
                   + std::to_string(LC.cellSizeFactor) + ", updating frequency "
                   + std::to_string(LC.updateFrequency) + ", sorting frequency "
                   + std::to_string(LC.sortFrequency) + " ...");
    }
    else  // BruteForce - still need basic LinkedCell params for insertion checks
    {
        LC.type           = LinkedCellType::HOST;
        LC.cellSizeFactor = T(1.0);
    }

    // Bounding volume
    DOMNode* nBoundingVolume = ReaderXML::getNode(collisionDetection, "BoundingVolume");
    if(nBoundingVolume)
    {
        std::string boundingVolumeType = ReaderXML::getNodeAttr_String(nBoundingVolume, "Type");
        if(boundingVolumeType == "OFF")
            CD.boundingVolumeType = BoundingVolumeType::OFF;
        else if(boundingVolumeType == "OBB")
            CD.boundingVolumeType = BoundingVolumeType::OBB;
        else if(boundingVolumeType == "OBC")
            CD.boundingVolumeType = BoundingVolumeType::OBC;
        else
            GAbort("Unknown bounding volume type! Aborting Grains!");
        GoutWI(9, "BoundingVolume: " + boundingVolumeType);
    }

    // Narrow phase detection
    DOMNode* nNarrowPhase = ReaderXML::getNode(collisionDetection, "NarrowPhase");
    if(nNarrowPhase)
    {
        std::string narrowPhaseType = ReaderXML::getNodeAttr_String(nNarrowPhase, "Type");
        if(narrowPhaseType == "GJK")
            CD.narrowPhaseType = NarrowPhaseType::GJK;
        else if(narrowPhaseType == "GJK_SV")
            CD.narrowPhaseType = NarrowPhaseType::GJK_SV;
        else
            GAbort("Unknown narrow phase type! Aborting Grains!");
        // GJK warm-start acceleration: "true" to enable (default: false)
        std::string gjkAcc = ReaderXML::getNodeAttr_String(nNarrowPhase, "Acceleration");
        CD.gjkAcceleration = (gjkAcc == "true");
        // Pre-built ShapeData for vtable-free GJK support evaluation
        std::string prebuilt = ReaderXML::getNodeAttr_String(nNarrowPhase, "PrebuiltShapes");
        CD.usePrebuiltShapes = (prebuilt == "true");
        GoutWI(9,
               "NarrowPhase: " + narrowPhaseType
                   + ", Acceleration: " + (CD.gjkAcceleration ? "true" : "false")
                   + ", PrebuiltShapes: " + (CD.usePrebuiltShapes ? "true" : "false"));
    }
    GoutWI(6, "Reading collision detection completed!");

    // ---------------------------------------------------------------------------------------------
    // Temporal setting and time integration
    DOMNode* tempSetting = ReaderXML::getNode(root, "TemporalSetting");
    if(tempSetting)
    {
        DOMNode* nTime  = ReaderXML::getNode(tempSetting, "TimeInterval");
        T        tStart = ReaderXML::getNodeAttr_Double(nTime, "Start");
        T        tEnd   = ReaderXML::getNodeAttr_Double(nTime, "End");
        T        tStep  = ReaderXML::getNodeAttr_Double(nTime, "dt");
        GP::m_tStart    = tStart;
        GP::m_dt        = tStep;
        GP::m_tEnd      = tEnd + 0.01 * tStep;  // Adding a small tolerance
        DOMNode* nTI    = ReaderXML::getNode(tempSetting, "TimeIntegration");
        if(nTI)
        {
            GoutWI(6, "Reading time integration model ...");
            TimeIntegratorFactory<T>::create(nTI, GP::m_dt, m_timeIntegrator);
            GP::m_isLeapFrog
                = (m_timeIntegrator[0]->getTimeIntegratorType() == SECONDORDERLEAPFROG);
            GoutWI(6, "Reading time integration model completed!");
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Components setup (reads rigid bodies, then contact forces, then builds ComponentManager)
    setupComponents(root, rootElement);
}

// -------------------------------------------------------------------------------------------------
// External force definition
template <typename T>
void Grains<T>::Forces(DOMElement* rootElement)
{
    assert(rootElement != NULL);
    DOMNode* root = ReaderXML::getNode(rootElement, "Forces");

    // Output message
    GoutWI(3, "Forces");

    // Read the forces
    if(root)
    {
        // Gravity
        DOMNode* nGravity = ReaderXML::getNode(root, "Gravity");
        GAssert(nGravity, "Gravity node is mandatory!");
        GrainsParameters<T>::m_gravity[X] = T(ReaderXML::getNodeAttr_Double(nGravity, "GX"));
        GrainsParameters<T>::m_gravity[Y] = T(ReaderXML::getNodeAttr_Double(nGravity, "GY"));
        GrainsParameters<T>::m_gravity[Z] = T(ReaderXML::getNodeAttr_Double(nGravity, "GZ"));
        GoutWI(6, "Gravity =", Vector3ToString(GrainsParameters<T>::m_gravity));
    }
}

// -------------------------------------------------------------------------------------------------
// Additional features of the simulation: insertion, post-processing
template <typename T>
void Grains<T>::AdditionalFeatures(DOMElement* rootElement)
{
    using GP = GrainsParameters<T>;

    // Output message
    GoutWI(3, "Simulation");
    // ---------------------------------------------------------------------------------------------
    // Checking if Simulation node is available
    assert(rootElement != NULL);
    DOMNode* root = ReaderXML::getNode(rootElement, "Simulation");
    GAssert(root, "Simulation node is mandatory!");

    // ---------------------------------------------------------------------------------------------
    // Verbosity
    DOMNode* nVerbosity = ReaderXML::getNode(root, "Verbosity");
    if(nVerbosity)
    {
        GP::m_verbosityFrequency = ReaderXML::getNodeAttr_Int(nVerbosity, "Frequency");
        GoutWI(6, "Verbosity frequency set to", GP::m_verbosityFrequency);
    }

    // ---------------------------------------------------------------------------------------------
    // Simulation-level timings
    DOMNode* nTimings = ReaderXML::getNode(root, "Timings");
    if(nTimings && ReaderXML::getNodeAttr_String(nTimings, "Enable") == "true")
    {
        GP::m_simTimer.enable(GP::m_isGPU);
        GoutWI(6, "Simulation timings enabled.");
    }

    // ---------------------------------------------------------------------------------------------
    // Insertion policies
    DOMNode* nInsertion = ReaderXML::getNode(root, "ParticleInsertion");
    GoutWI(6, "Reading insertion policies ...");
    if(nInsertion)
        m_insertion = std::make_unique<Insertion<T>>(nInsertion);
    else
    {
        GoutWI(9, "No policy found, setting the insertion policy to default");
        m_insertion = std::make_unique<Insertion<T>>();
    }
    GoutWI(6, "Reading insertion policies completed.");

    // ---------------------------------------------------------------------------------------------
    // Post-processing writers
    DOMNode* nPostProcessing = ReaderXML::getNode(root, "PostProcessing");
    if(nPostProcessing)
    {
        GoutWI(3, "Post-processing");
        // Post-processing save time
        DOMNode* nTime = ReaderXML::getNode(nPostProcessing, "TimeSave");
        T        tStart, tEnd;
        if(ReaderXML::hasNodeAttr(nTime, "Start"))
            tStart = ReaderXML::getNodeAttr_Double(nTime, "Start");
        else
            tStart = GrainsParameters<T>::m_tStart;
        if(ReaderXML::hasNodeAttr(nTime, "End"))
            tEnd = ReaderXML::getNodeAttr_Double(nTime, "End");
        else
            tEnd = GrainsParameters<T>::m_tEnd;
        T tStep = ReaderXML::getNodeAttr_Double(nTime, "dt");
        for(T t = tStart; t <= tEnd; t += tStep)
            GrainsParameters<T>::m_tSave.push(t);
        // Save for tEnd as well
        GrainsParameters<T>::m_tSave.push(tEnd);

        // Post-processing writers
        DOMNode* nWriters = ReaderXML::getNode(nPostProcessing, "Writers");
        if(nWriters)
        {
            GoutWI(6, "Reading the post processing writers ...");
            DOMNodeList* allPPW = ReaderXML::getNodes(nWriters);
            for(uint i = 0; i < allPPW->getLength(); ++i)
            {
                DOMNode* nPPW = allPPW->item(i);
                m_postProcessor.push_back(std::unique_ptr<PostProcessingWriter<T>>(
                    PostProcessingWriterFactory<T>::create(nPPW)));
            }
            GoutWI(6, "Reading the post processing writers completed!");
        }
    }
    else
        GoutWI(6, "No postprocessing writer!");
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Grains<float>;
template class Grains<double>;
template void Grains<float>::postProcess<MemType::HOST>(
    const std::unique_ptr<ComponentManager<float, MemType::HOST>>&);
template void Grains<float>::postProcess<MemType::DEVICE>(
    const std::unique_ptr<ComponentManager<float, MemType::DEVICE>>&);
template void Grains<double>::postProcess<MemType::HOST>(
    const std::unique_ptr<ComponentManager<double, MemType::HOST>>&);
template void Grains<double>::postProcess<MemType::DEVICE>(
    const std::unique_ptr<ComponentManager<double, MemType::DEVICE>>&);