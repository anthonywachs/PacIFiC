#ifndef _COLLISIONDETECTIONBENCHMARK_HH_
#define _COLLISIONDETECTIONBENCHMARK_HH_

#include <algorithm>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "BenchmarkConfig.hh"
#include "Box.hh"
#include "CSVWriter.hh"
#include "ComponentManagerCommon.hh"
#include "ComponentManagerGPU_Kernels.hh"
#include "GJK.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "Insertion.hh"
#include "InsertionWindow.hh"
#include "Kinematics.hh"
#include "LinkedCell.hh"
#include "LinkedCellFactory.hh"
#include "NeighborList.hh"
#include "NeighborListFactory.hh"
#include "Quaternion.hh"
#include "QuaternionMath.hh"
#include "RigidBody.hh"
#include "RigidBodyFactory.hh"
#include "Sphere.hh"
#include "StepTimer.hh"
#include "Superquadric.hh"
#include "Transform3.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief GJK Collision Detection Benchmark

    This class benchmarks GJK collision detection with neighbor list generation.
    Directly uses detectCollisionsComponents_common without ComponentManager objects.
    Results are exported to CSV for post-processing.

    @author A.Yazdani - 2026 - Collision Detection Performance Validation */
// =================================================================================================
template <typename T>
class CollisionDetectionBenchmark
{
private:
    BenchmarkConfig            m_config;
    std::unique_ptr<CSVWriter> m_csvWriter;

    // Storage for particles
    std::unique_ptr<Convex<T>> m_shapeTemplate;

    // Host-side storage
    GrainsMemBuffer<RigidBody<T>*, MemType::HOST>  m_rigidBodies;
    GrainsMemBuffer<Vector3<T>, MemType::HOST>     m_positions;
    GrainsMemBuffer<Quaternion<T>, MemType::HOST>  m_quaternions;
    GrainsMemBuffer<Transform3<T>, MemType::HOST>  m_transforms;
    GrainsMemBuffer<Kinematics<T>, MemType::HOST>  m_kinematics;
    GrainsMemBuffer<Vector3<T>, MemType::HOST>     m_relativePositions;
    GrainsMemBuffer<Quaternion<T>, MemType::HOST>  m_relativeQuaternions;
    GrainsMemBuffer<Transform3<T>, MemType::HOST>  m_relativeTransforms;
    GrainsMemBuffer<ContactInfo<T>, MemType::HOST> m_contactInfo;

    // GPU-side storage
    GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>  m_d_rigidBodies;
    GrainsMemBuffer<Vector3<T>, MemType::DEVICE>     m_d_positions;
    GrainsMemBuffer<Quaternion<T>, MemType::DEVICE>  m_d_quaternions;
    GrainsMemBuffer<Transform3<T>, MemType::DEVICE>  m_d_transforms;
    GrainsMemBuffer<Vector3<T>, MemType::DEVICE>     m_d_relativePositions;
    GrainsMemBuffer<Quaternion<T>, MemType::DEVICE>  m_d_relativeQuaternions;
    GrainsMemBuffer<Transform3<T>, MemType::DEVICE>  m_d_relativeTransforms;
    GrainsMemBuffer<ContactInfo<T>, MemType::DEVICE> m_d_contactInfo;

public:
    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with configuration */
    CollisionDetectionBenchmark(const BenchmarkConfig& config,
                                const std::string&     csvFilename,
                                bool                   appendToCSV = false)
        : m_config(config)
    {
        m_csvWriter = std::make_unique<CSVWriter>(csvFilename, appendToCSV);

        // Initialize CSV header only if not appending
        if(!appendToCSV)
        {
            std::vector<std::string> columns = {"TrialID",
                                                "Platform",
                                                "Precision",
                                                "ParticleCount",
                                                "ShapeType",
                                                "ParticleSize",
                                                "AspectRatio",
                                                "GJKAlgo",
                                                "GJKRepresentation",
                                                "UseRelativeTransform",
                                                "NeighborListTime_ms",
                                                "RelativeTransformTime_ms",
                                                "GJKTime_ms",
                                                "TotalTime_ms",
                                                "PairCount"};
            m_csvWriter->writeHeader(columns);
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    ~CollisionDetectionBenchmark()
    {
        // Delete all rigid bodies
        for(uint i = 0; i < m_config.numParticles; ++i)
        {
            if(m_rigidBodies[i] != nullptr)
            {
                delete m_rigidBodies[i];
                m_rigidBodies[i] = nullptr;
            }
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Run benchmark with current configuration */
    void runBenchmark()
    {
        // Create shape template
        auto size = Vector3<T>(static_cast<T>(m_config.particleSize[X]),
                               static_cast<T>(m_config.particleSize[Y]),
                               static_cast<T>(m_config.particleSize[Z]));
        switch(m_config.shapeType)
        {
        case ParticleShapeType::BOX:
            m_shapeTemplate = std::make_unique<Box<T>>(size[X], size[Y], size[Z]);
            break;
        case ParticleShapeType::SPHERE:
            m_shapeTemplate = std::make_unique<Sphere<T>>(size[X]);
            break;
        case ParticleShapeType::SUPERQUADRIC:
            m_shapeTemplate
                = std::make_unique<Superquadric<T>>(size[X], size[Y], size[Z], T(2.0), T(2.0));
            break;
        default:
            GAbort("Invalid shape type!");
        }

        // Allocate memory
        uint nTotal = m_config.numParticles;
        m_rigidBodies.initialize(nTotal);
        m_positions.initialize(nTotal);
        m_quaternions.initialize(nTotal);
        m_kinematics.initialize(nTotal);
        m_transforms.initialize(nTotal);
        if(m_config.platform == PLATFORM::GPU || m_config.platform == PLATFORM::BOTH)
        {
            m_d_rigidBodies.initialize(nTotal);
            m_d_positions.initialize(nTotal);
            m_d_quaternions.initialize(nTotal);
            m_d_transforms.initialize(nTotal);
        }

        // Create rigid bodies
        for(uint i = 0; i < nTotal; ++i)
        {
            Convex<T>* cvx   = m_shapeTemplate->clone();
            m_rigidBodies[i] = new RigidBody<T>(cvx, T(0), 0, 1);
            m_positions[i]   = Vector3<T>(T(0), T(0), T(0));
            m_quaternions[i] = Quaternion<T>(T(1), T(0), T(0), T(0));
        }

        // Create CollisionDetectionParameters from config
        CollisionDetectionParameters<T> cdParams;
        cdParams.neighborListType               = NeighborListType::LINKEDCELL;
        cdParams.linkedCellParameters.type      = LinkedCellType::ATOMIC;
        cdParams.linkedCellParameters.minCorner = Vector3<T>(static_cast<T>(m_config.domainMin[X]),
                                                             static_cast<T>(m_config.domainMin[Y]),
                                                             static_cast<T>(m_config.domainMin[Z]));
        cdParams.linkedCellParameters.maxCorner = Vector3<T>(static_cast<T>(m_config.domainMax[X]),
                                                             static_cast<T>(m_config.domainMax[Y]),
                                                             static_cast<T>(m_config.domainMax[Z]));
        cdParams.linkedCellParameters.minCellSize
            = 2. * m_shapeTemplate->computeCircumscribedRadius();
        cdParams.linkedCellParameters.cellSizeFactor  = 1.0;
        cdParams.linkedCellParameters.updateFrequency = 0;
        cdParams.linkedCellParameters.sortFrequency   = 0;

        // Insert particles
        insertParticles(cdParams.linkedCellParameters, m_config.randomSeed);

        // Estimate maximum number of pairs for buffer allocation
        uint maxPairs = cdParams.linkedCellParameters.initialNumberOfPairsPerParticle
                        * m_config.numParticles;  // Conservative estimate

        // Allocate collision detection buffers for CPU
        m_relativePositions.initialize(maxPairs);
        m_relativeQuaternions.initialize(maxPairs);
        m_relativeTransforms.initialize(maxPairs);
        m_contactInfo.initialize(maxPairs);

        for(uint i = 0; i < nTotal; ++i)
        {
            m_transforms[i] = Transform3<T>(m_quaternions[i], m_positions[i]);
        }

        // Allocate collision detection buffers for GPU if needed
        if(m_config.platform == PLATFORM::GPU || m_config.platform == PLATFORM::BOTH)
        {
            // Copy all data to device
            RigidBodyFactory<T>::copyHostToDevice(m_rigidBodies, m_d_rigidBodies);
            m_d_positions.copyFrom(m_positions);
            m_d_quaternions.copyFrom(m_quaternions);
            m_d_transforms.copyFrom(m_transforms);

            m_d_relativePositions.initialize(maxPairs);
            m_d_relativeQuaternions.initialize(maxPairs);
            m_d_relativeTransforms.initialize(maxPairs);
            m_d_contactInfo.initialize(maxPairs);
        }

        // Create neighbor lists once per platform (they only depend on positions)
        std::unique_ptr<NeighborList<T, MemType::HOST>>   cpuNeighborList;
        std::unique_ptr<NeighborList<T, MemType::DEVICE>> gpuNeighborList;
        const uint2*                                      cpuPairList  = nullptr;
        const uint2*                                      gpuPairList  = nullptr;
        uint                                              cpuPairCount = 0;
        uint                                              gpuPairCount = 0;

        if(m_config.platform == PLATFORM::CPU || m_config.platform == PLATFORM::BOTH)
        {
            cpuNeighborList = NeighborListFactory<T, MemType::HOST>::create(&m_rigidBodies,
                                                                            m_positions,
                                                                            m_quaternions,
                                                                            cdParams,
                                                                            0,
                                                                            m_config.numParticles);
            cpuNeighborList->updateNeighborList(m_positions, 0, m_config.numParticles);
            cpuPairCount = cpuNeighborList->getSize();
            cpuPairList  = cpuNeighborList->getData();
        }

        if(m_config.platform == PLATFORM::GPU || m_config.platform == PLATFORM::BOTH)
        {
            gpuNeighborList
                = NeighborListFactory<T, MemType::DEVICE>::create(&m_d_rigidBodies,
                                                                  m_d_positions,
                                                                  m_d_quaternions,
                                                                  cdParams,
                                                                  0,
                                                                  m_config.numParticles);
            gpuNeighborList->updateNeighborList(m_d_positions, 0, m_config.numParticles);
            gpuPairCount = gpuNeighborList->getSize();
            gpuPairList  = gpuNeighborList->getData();
        }

        // Run all configurations, with all trials for each config before moving to next config
        std::vector<GJKRepresentationType> representations
            = {GJKRepresentationType::QUATERNION, GJKRepresentationType::TRANSFORM};
        std::vector<GJKVariantType> gjkVariants
            = {GJKVariantType::JOHNSON, GJKVariantType::SIGNEDVOLUME};
        std::vector<bool> transformModes = {true, false};  // relative vs global

        for(auto representation : representations)
        {
            for(auto gjkVariant : gjkVariants)
            {
                for(bool useRelativeTransform : transformModes)
                {
                    // Run all trials for this configuration on CPU
                    if(m_config.platform == PLATFORM::CPU || m_config.platform == PLATFORM::BOTH)
                    {
                        for(uint trial = 0; trial < m_config.numTrials; ++trial)
                        {
                            runSingleConfiguration<MemType::HOST>(cpuPairList,
                                                                  cpuPairCount,
                                                                  trial,
                                                                  representation,
                                                                  gjkVariant,
                                                                  useRelativeTransform);
                        }
                        if(m_config.validateContacts)
                        {
                            writeContactsToFile(cpuNeighborList.get(), &m_contactInfo);
                        }
                    }
                    // Run all trials for this configuration on GPU
                    if(m_config.platform == PLATFORM::GPU || m_config.platform == PLATFORM::BOTH)
                    {
                        for(uint trial = 0; trial < m_config.numTrials; ++trial)
                        {
                            runSingleConfiguration<MemType::DEVICE>(gpuPairList,
                                                                    gpuPairCount,
                                                                    trial,
                                                                    representation,
                                                                    gjkVariant,
                                                                    useRelativeTransform);
                        }
                    }
                    if(m_config.validateContacts)
                    {
                        writeContactsToFile(gpuNeighborList.get(), &m_d_contactInfo);
                    }
                }
            }
        }
    }

private:
    // ---------------------------------------------------------------------------------------------
    /** @brief Insert particles using Insertion class */
    void insertParticles(const LinkedCellParameters<T>& linkedCellParams, uint seed)
    {
        // Use Insertion class for particle placement
        // Create insertion window covering the domain for positions
        std::vector<InsertionWindow<T>> positionWindows;
        Vector3<T>                      domainMin(static_cast<T>(m_config.domainMin[X]),
                             static_cast<T>(m_config.domainMin[Y]),
                             static_cast<T>(m_config.domainMin[Z]));
        Vector3<T>                      domainMax(static_cast<T>(m_config.domainMax[X]),
                             static_cast<T>(m_config.domainMax[Y]),
                             static_cast<T>(m_config.domainMax[Z]));
        positionWindows.emplace_back(domainMin, domainMax, seed);

        // Create insertion windows for random orientations (Euler angles: 0 to 2*pi)
        std::vector<InsertionWindow<T>> orientationWindows;
        Vector3<T>                      angleMin(T(0), T(0), T(0));
        Vector3<T>                      angleMax(T(2.0 * M_PI), T(2.0 * M_PI), T(2.0 * M_PI));
        orientationWindows.emplace_back(angleMin, angleMax, seed + 1);

        // Create Insertion object
        auto insertion = std::make_unique<Insertion<T>>();
        insertion->setPositionInsertionInfo(positionWindows);
        insertion->setOrientationInsertionInfo(orientationWindows);

        // Perform insertion
        insertion->insert(&m_rigidBodies,
                          m_positions,
                          m_quaternions,
                          m_kinematics,
                          linkedCellParams,
                          0,
                          m_config.numParticles);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Run a single GJK configuration with timing */
    template <MemType M>
    void runSingleConfiguration(const uint2*          pairList,
                                uint                  pairCount,
                                uint                  trialID,
                                GJKRepresentationType representation,
                                GJKVariantType        gjkVariant,
                                bool                  useRelativeTransform)
    {
        // Select appropriate buffers based on memory type
        GrainsMemBuffer<RigidBody<T>*, M>*  rbBuffer;
        GrainsMemBuffer<Vector3<T>, M>*     posBuffer;
        GrainsMemBuffer<Quaternion<T>, M>*  quatBuffer;
        GrainsMemBuffer<Transform3<T>, M>*  transformBuffer;
        GrainsMemBuffer<Vector3<T>, M>*     relPosBuffer;
        GrainsMemBuffer<Quaternion<T>, M>*  relQuatBuffer;
        GrainsMemBuffer<Transform3<T>, M>*  relTransformBuffer;
        GrainsMemBuffer<ContactInfo<T>, M>* contactBuffer;

        if constexpr(M == MemType::HOST)
        {
            rbBuffer           = &m_rigidBodies;
            posBuffer          = &m_positions;
            quatBuffer         = &m_quaternions;
            transformBuffer    = &m_transforms;
            relPosBuffer       = &m_relativePositions;
            relQuatBuffer      = &m_relativeQuaternions;
            relTransformBuffer = &m_relativeTransforms;
            contactBuffer      = &m_contactInfo;
        }
        else
        {
            rbBuffer           = &m_d_rigidBodies;
            posBuffer          = &m_d_positions;
            quatBuffer         = &m_d_quaternions;
            transformBuffer    = &m_d_transforms;
            relPosBuffer       = &m_d_relativePositions;
            relQuatBuffer      = &m_d_relativeQuaternions;
            relTransformBuffer = &m_d_relativeTransforms;
            contactBuffer      = &m_d_contactInfo;
        }

        // Resize contact buffer if needed
        contactBuffer->resize(pairCount);
        relPosBuffer->resize(pairCount);
        relQuatBuffer->resize(pairCount);
        relTransformBuffer->resize(pairCount);

        // Kernel launch parameters
        uint numThreads, numBlocks;
        computeOptimalThreadsAndBlocks(pairCount,
                                       GrainsParameters<T>::m_GPU,
                                       numBlocks,
                                       numThreads);

        // Step 2: Compute relative transformations (if using relative mode)
        StepTimer relTransformTimer;
        double    relTransformTime = 0.0;

        if(useRelativeTransform && pairCount > 0)
        {
            relTransformTimer.start();

            if(representation == GJKRepresentationType::TRANSFORM)
            {
                // Use Transform3 representation
                if constexpr(M == MemType::HOST)
                {
                    // CPU version
                    for(uint pairID = 0; pairID < pairCount; ++pairID)
                    {
                        computeRelativeTransformations_common<T>(pairList,
                                                                 transformBuffer->getData(),
                                                                 relTransformBuffer->getData(),
                                                                 pairID);
                    }
                }
                else
                {
                    // GPU version - launch kernel
                    computeRelativeTransformations_Kernel<<<numBlocks, numThreads>>>(
                        pairList,
                        transformBuffer->getData(),
                        relTransformBuffer->getData(),
                        pairCount);
                    cudaErrCheck(cudaDeviceSynchronize());
                }
            }
            else  // Quaternion representation
            {
                // Use Quaternion representation
                if constexpr(M == MemType::HOST)
                {
                    // CPU version
                    for(uint pairID = 0; pairID < pairCount; ++pairID)
                    {
                        computeRelativeTransformations_common<T>(pairList,
                                                                 posBuffer->getData(),
                                                                 quatBuffer->getData(),
                                                                 relPosBuffer->getData(),
                                                                 relQuatBuffer->getData(),
                                                                 pairID);
                    }
                }
                else
                {
                    // GPU version - launch kernel
                    computeRelativeTransformations_Kernel<<<numBlocks, numThreads>>>(
                        pairList,
                        posBuffer->getData(),
                        quatBuffer->getData(),
                        relPosBuffer->getData(),
                        relQuatBuffer->getData(),
                        pairCount);
                    cudaErrCheck(cudaDeviceSynchronize());
                }
            }
            relTransformTimer.stop();
            relTransformTime = relTransformTimer.getElapsedMilliseconds();
        }

        // Step 3: Run GJK collision detection
        StepTimer gjkTimer;
        gjkTimer.start();

        if(pairCount > 0 && useRelativeTransform)
        {
            // Use relative transformations
            if(representation == GJKRepresentationType::TRANSFORM)
            {
                // Transform3 representation
                if constexpr(M == MemType::HOST)
                {
                    // CPU version
                    if(gjkVariant == GJKVariantType::JOHNSON)
                    {
                        for(uint pairID = 0; pairID < pairCount; ++pairID)
                        {
                            detectCollisionsComponents_common<T, GJKType::JOHNSON, false>(
                                pairList,
                                rbBuffer->getData(),
                                relTransformBuffer->getData(),
                                contactBuffer->getData(),
                                pairID);
                        }
                    }
                    else  // SIGNEDVOLUME
                    {
                        for(uint pairID = 0; pairID < pairCount; ++pairID)
                        {
                            detectCollisionsComponents_common<T, GJKType::SIGNEDVOLUME, false>(
                                pairList,
                                rbBuffer->getData(),
                                relTransformBuffer->getData(),
                                contactBuffer->getData(),
                                pairID);
                        }
                    }
                }
                else
                {
                    // GPU version - launch kernel
                    detectCollisionsComponents_Kernel<<<numBlocks, numThreads>>>(
                        pairList,
                        rbBuffer->getData(),
                        relTransformBuffer->getData(),
                        contactBuffer->getData(),
                        pairCount);
                    cudaErrCheck(cudaDeviceSynchronize());
                }
            }
            else  // Quaternion representation
            {
                // Quaternion representation
                if constexpr(M == MemType::HOST)
                {
                    // CPU version
                    if(gjkVariant == GJKVariantType::JOHNSON)
                    {
                        for(uint pairID = 0; pairID < pairCount; ++pairID)
                        {
                            detectCollisionsComponents_common<T, GJKType::JOHNSON, false>(
                                pairList,
                                rbBuffer->getData(),
                                relPosBuffer->getData(),
                                relQuatBuffer->getData(),
                                contactBuffer->getData(),
                                pairID);
                        }
                    }
                    else  // SIGNEDVOLUME
                    {
                        for(uint pairID = 0; pairID < pairCount; ++pairID)
                        {
                            detectCollisionsComponents_common<T, GJKType::SIGNEDVOLUME, false>(
                                pairList,
                                rbBuffer->getData(),
                                relPosBuffer->getData(),
                                relQuatBuffer->getData(),
                                contactBuffer->getData(),
                                pairID);
                        }
                    }
                }
                else
                {
                    // GPU version - launch kernel
                    detectCollisionsComponents_Kernel<<<numBlocks, numThreads>>>(
                        pairList,
                        rbBuffer->getData(),
                        relPosBuffer->getData(),
                        relQuatBuffer->getData(),
                        contactBuffer->getData(),
                        pairCount);
                    cudaErrCheck(cudaDeviceSynchronize());
                }
            }
        }
        else if(pairCount > 0)
        {
            // Use global coordinates directly
            if(representation == GJKRepresentationType::TRANSFORM)
            {
                // Transform3 representation
                if constexpr(M == MemType::HOST)
                {
                    // CPU version
                    if(gjkVariant == GJKVariantType::JOHNSON)
                    {
                        for(uint pairID = 0; pairID < pairCount; ++pairID)
                        {
                            detectCollisionsComponentsGlobal_common<T, GJKType::JOHNSON, false>(
                                pairList,
                                rbBuffer->getData(),
                                transformBuffer->getData(),
                                contactBuffer->getData(),
                                pairID);
                        }
                    }
                    else  // SIGNEDVOLUME
                    {
                        for(uint pairID = 0; pairID < pairCount; ++pairID)
                        {
                            detectCollisionsComponentsGlobal_common<T,
                                                                    GJKType::SIGNEDVOLUME,
                                                                    false>(
                                pairList,
                                rbBuffer->getData(),
                                transformBuffer->getData(),
                                contactBuffer->getData(),
                                pairID);
                        }
                    }
                }
                else
                {
                    // GPU version - launch kernel
                    detectCollisionsComponentsGlobal_Kernel<<<numBlocks, numThreads>>>(
                        pairList,
                        rbBuffer->getData(),
                        transformBuffer->getData(),
                        contactBuffer->getData(),
                        pairCount);
                    cudaErrCheck(cudaDeviceSynchronize());
                }
            }
            else  // Quaternion representation
            {
                // Quaternion representation
                if constexpr(M == MemType::HOST)
                {
                    // CPU version
                    if(gjkVariant == GJKVariantType::JOHNSON)
                    {
                        for(uint pairID = 0; pairID < pairCount; ++pairID)
                        {
                            detectCollisionsComponentsGlobal_common<T, GJKType::JOHNSON, false>(
                                pairList,
                                rbBuffer->getData(),
                                posBuffer->getData(),
                                quatBuffer->getData(),
                                contactBuffer->getData(),
                                pairID);
                        }
                    }
                    else  // SIGNEDVOLUME
                    {
                        for(uint pairID = 0; pairID < pairCount; ++pairID)
                        {
                            detectCollisionsComponentsGlobal_common<T,
                                                                    GJKType::SIGNEDVOLUME,
                                                                    false>(pairList,
                                                                           rbBuffer->getData(),
                                                                           posBuffer->getData(),
                                                                           quatBuffer->getData(),
                                                                           contactBuffer->getData(),
                                                                           pairID);
                        }
                    }
                }
                else
                {
                    // GPU version - launch kernel
                    detectCollisionsComponentsGlobal_Kernel<<<numBlocks, numThreads>>>(
                        pairList,
                        rbBuffer->getData(),
                        posBuffer->getData(),
                        quatBuffer->getData(),
                        contactBuffer->getData(),
                        pairCount);
                    cudaErrCheck(cudaDeviceSynchronize());
                }
            }
        }

        gjkTimer.stop();
        double gjkTime = gjkTimer.getElapsedMilliseconds();

        // Step 4: Calculate total time
        double totalTime = relTransformTime + gjkTime;

        // Step 5: Write results to CSV (neighbor list time reported separately)
        m_csvWriter->writeRow(
            trialID,
            M == MemType::HOST ? "CPU" : "GPU",
            sizeof(T) == sizeof(float) ? "Single" : "Double",
            m_config.numParticles,
            m_config.shapeType == ParticleShapeType::BOX            ? "Box"
            : m_config.shapeType == ParticleShapeType::SPHERE       ? "Sphere"
            : m_config.shapeType == ParticleShapeType::SUPERQUADRIC ? "Superquadric"
                                                                    : "Unknown",
            m_config.particleSize[X],
            m_config.aspectRatio,
            gjkVariant == GJKVariantType::JOHNSON ? "Johnson" : "SignedVolume",
            representation == GJKRepresentationType::TRANSFORM ? "Transform" : "Quaternion",
            useRelativeTransform ? 1 : 0,
            0.0,  // Neighbor list time measured separately, once per platform
            relTransformTime,
            gjkTime,
            totalTime,
            pairCount);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Write contact information to file
        @param neighborList Neighbor list containing pair info
        @param contactBuffer Contact info buffer */
    template <MemType M>
    void writeContactsToFile(const NeighborList<T, M>*                 neighborList,
                             const GrainsMemBuffer<ContactInfo<T>, M>* contactBuffer)
    {
        if(neighborList->getSize() == 0)
        {
            Gout("No contacts to write for platform "
                 + std::string((M == MemType::HOST) ? "CPU" : "GPU"));
            return;
        }

        // Get pair list
        uint                                  numContacts = neighborList->getSize();
        GrainsMemBuffer<uint2, MemType::HOST> pairList(numContacts);
        neighborList->getBuffer().copyTo(pairList);

        // Get contact info
        GrainsMemBuffer<ContactInfo<T>, MemType::HOST> contactInfo(numContacts);
        contactBuffer->copyTo(contactInfo);

        // Create filename
        std::string   filename = "data/contacts_" + std::to_string(std::rand()) + ".txt";
        std::ofstream outFile(filename);

        if(!outFile.is_open())
        {
            GAbort(("Failed to open file: " + filename).c_str());
        }

        // Write header
        outFile << "# PairID Body1 Body2 ContactPoint_X ContactPoint_Y ContactPoint_Z "
                << "ContactVector_X ContactVector_Y ContactVector_Z OverlapDistance\n";

        // Create indices for sorting
        std::vector<uint> sortedIndices(numContacts);
        for(uint i = 0; i < numContacts; ++i)
        {
            sortedIndices[i] = i;
        }

        // Sort indices by Body1 (pair.x) ascending, then Body2 (pair.y) ascending
        std::sort(sortedIndices.begin(), sortedIndices.end(), [&pairList](uint a, uint b) {
            if(pairList[a].x != pairList[b].x)
                return pairList[a].x < pairList[b].x;
            return pairList[a].y < pairList[b].y;
        });

        // Write contact data in sorted order
        for(uint idx = 0; idx < numContacts; ++idx)
        {
            uint        i       = sortedIndices[idx];
            const auto& contact = contactInfo[i];
            const auto& pair    = pairList[i];

            Vector3<T> contactPt  = contact.getContactPoint();
            Vector3<T> contactVec = contact.getContactVector();
            T          overlap    = contact.getOverlapDistance();

            outFile << i << " " << pair.x << " " << pair.y << " " << contactPt[X] << " "
                    << contactPt[Y] << " " << contactPt[Z] << " " << contactVec[X] << " "
                    << contactVec[Y] << " " << contactVec[Z] << " " << overlap << "\n";
        }

        outFile.close();
        Gout("Wrote " + std::to_string(numContacts) + " contacts to " + filename);
    }
};

#endif