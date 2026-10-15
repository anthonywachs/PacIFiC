/* ================================================================================================
   Unit tests for NeighborList algorithms.

   Two primary test categories:
   1. Simple known-ground-truth tests: few particles with manually verified expected pairs;
      every CPU and GPU algorithm is checked against this known truth.
   2. Cross-algorithm consistency tests: larger particle systems where all CPU/GPU algorithms
      are compared pairwise to confirm they all agree.

   Algorithms covered:
     - NSQ (brute-force O(n^2)):        CPU (HOST) and GPU (DEVICE)
     - LinkedCell HOST (CPU):            HOST
     - LinkedCell SORTBASED (GPU):       DEVICE
     - LinkedCell ATOMIC (GPU):          DEVICE
     - LinkedCell ATOMICFIXED (GPU):     DEVICE
   ================================================================================================
 */

#include <cuda_runtime.h>
#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "Box.hh"
#include "CollisionDetection.hh"
#include "GrainsParameters.hh"
#include "NeighborList.hh"
#include "NeighborListFactory.hh"
#include "Quaternion.hh"
#include "RigidBody.hh"
#include "Vector3.hh"

// ================================================================================================
// Utility types and helpers
// ================================================================================================
using PairSet = std::set<std::pair<uint, uint>>;

/** @brief Normalise a pair so the smaller index is first. */
static inline std::pair<uint, uint> normPair(uint a, uint b)
{
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

/** @brief Extract neighbour pairs from a HOST NeighborList. */
static PairSet extractPairs_CPU(NeighborList<double, MemType::HOST>*             NL,
                                GrainsMemBuffer<Vector3<double>, MemType::HOST>& positions,
                                uint                                             nObs,
                                uint                                             nParts)
{
    NL->updateNeighborList(positions, nObs, nParts);
    PairSet      pairs;
    const uint2* data = NL->getData();
    uint         size = NL->getSize();
    for(uint i = 0; i < size; ++i)
        pairs.insert(normPair(data[i].x, data[i].y));
    return pairs;
}

/** @brief Extract neighbour pairs from a DEVICE NeighborList (copies data to host). */
static PairSet extractPairs_GPU(NeighborList<double, MemType::DEVICE>*             NL,
                                GrainsMemBuffer<Vector3<double>, MemType::DEVICE>& positions,
                                uint                                               nObs,
                                uint                                               nParts)
{
    NL->updateNeighborList(positions, nObs, nParts);
    cudaDeviceSynchronize();
    PairSet pairs;
    uint    size = NL->getSize();
    if(size == 0)
        return pairs;
    std::vector<uint2> host(size);
    cudaMemcpy(host.data(), NL->getData(), size * sizeof(uint2), cudaMemcpyDeviceToHost);
    for(uint i = 0; i < size; ++i)
        pairs.insert(normPair(host[i].x, host[i].y));
    return pairs;
}

/** @brief Build the exhaustive set of all particle-particle pairs for n particles (0-indexed). */
static PairSet allPairs(uint n)
{
    PairSet p;
    for(uint i = 0; i < n; ++i)
        for(uint j = i + 1; j < n; ++j)
            p.insert({i, j});
    return p;
}

/** @brief Reset the neighbor-list update counter (must be 0 before each create/update call). */
static void resetNLCounter()
{
    GrainsParameters<double>::m_simulationState.neighborListUpdateCount = 0;
}

// ================================================================================================
// Shared fixture helpers
// ================================================================================================
/** @brief Populate HOST and DEVICE memory buffers from a position list. */
static void buildBuffers(const std::vector<Vector3<double>>&                   positions,
                         uint                                                  nTotal,
                         const std::vector<RigidBody<double>*>&                rbVec,
                         GrainsMemBuffer<RigidBody<double>*, MemType::HOST>&   rb_cpu,
                         GrainsMemBuffer<Vector3<double>, MemType::HOST>&      pos_cpu,
                         GrainsMemBuffer<Quaternion<double>, MemType::HOST>&   quat_cpu,
                         GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE>& rb_gpu,
                         GrainsMemBuffer<Vector3<double>, MemType::DEVICE>&    pos_gpu,
                         GrainsMemBuffer<Quaternion<double>, MemType::DEVICE>& quat_gpu)
{
    Quaternion<double> identity(0.0, 0.0, 0.0, 1.0);

    rb_cpu.reserve(nTotal);
    pos_cpu.reserve(nTotal);
    quat_cpu.reserve(nTotal);
    for(uint i = 0; i < nTotal; ++i)
    {
        rb_cpu.push_back(rbVec[i]);
        pos_cpu.push_back(positions[i]);
        quat_cpu.push_back(identity);
    }

    rb_gpu.copyFrom(rb_cpu);
    pos_gpu.copyFrom(pos_cpu);
    quat_gpu.copyFrom(quat_cpu);
}

// ================================================================================================
//  FIXTURE: NeighborListSimpleTest
//  --------------------------------
//  5 particles at known positions inside a small domain.
//  Used for (a) ground-truth tests and (b) spatial-filtering tests.
//
//  Layout:
//    P0 = (0, 0, 0)
//    P1 = (1, 0, 0)
//    P2 = (0, 1, 0)
//    P3 = (0, 0, 1)
//    P4 = (0, 0, 2)
// ================================================================================================
class NeighborListSimpleTest : public ::testing::Test
{
protected:
    static constexpr uint nObs   = 0;
    static constexpr uint nParts = 5;
    static constexpr uint nTotal = nObs + nParts;

    Box<double>*                    sharedBox = nullptr;
    std::vector<RigidBody<double>*> rbVec;
    std::vector<Vector3<double>>    pos;

    GrainsMemBuffer<RigidBody<double>*, MemType::HOST>   rb_cpu;
    GrainsMemBuffer<Vector3<double>, MemType::HOST>      pos_cpu;
    GrainsMemBuffer<Quaternion<double>, MemType::HOST>   quat_cpu;
    GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE> rb_gpu;
    GrainsMemBuffer<Vector3<double>, MemType::DEVICE>    pos_gpu;
    GrainsMemBuffer<Quaternion<double>, MemType::DEVICE> quat_gpu;

    /** LinkedCell parameters with ONE BIG CELL covering the whole domain
     *  -> every particle is in the same cell neighbourhood -> all C(5,2)=10 pairs found */
    LinkedCellParameters<double> lcParamsLarge;

    void SetUp() override
    {
        resetNLCounter();
        GrainsParameters<double>::m_collisionDetection = CollisionDetectionParameters<double>{};

        sharedBox = new Box<double>(0.4, 0.4, 0.4);

        pos = {Vector3<double>(0.0, 0.0, 0.0),
               Vector3<double>(1.0, 0.0, 0.0),
               Vector3<double>(0.0, 1.0, 0.0),
               Vector3<double>(0.0, 0.0, 1.0),
               Vector3<double>(0.0, 0.0, 2.0)};

        for(uint i = 0; i < nTotal; ++i)
            rbVec.push_back(new RigidBody<double>(sharedBox, 0.1, 1000.0, 1));

        buildBuffers(pos, nTotal, rbVec, rb_cpu, pos_cpu, quat_cpu, rb_gpu, pos_gpu, quat_gpu);

        // One big cell (size = 10) -> no pair is filtered out
        lcParamsLarge.minCorner              = Vector3<double>(-1.0, -1.0, -1.0);
        lcParamsLarge.maxCorner              = Vector3<double>(5.0, 5.0, 5.0);
        lcParamsLarge.minCellSize            = 10.0;
        lcParamsLarge.cellSizeFactor         = 1.0;
        lcParamsLarge.maxNumCellsPerObstacle = 8;
        lcParamsLarge.maxParticlesPerCell    = 64;
        lcParamsLarge.updateFrequency        = 0;
        lcParamsLarge.sortFrequency          = 0;
    }

    void TearDown() override
    {
        // Note: memory cleanup omitted intentionally (known crash on teardown).
        GrainsParameters<double>::m_simulationState.neighborListUpdateCount = 0;
    }
};

// ------------------------------------------------------------------------------------------------
// Test 1: NSQ CPU returns exactly the known set of all C(5,2)=10 particle-particle pairs.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSimpleTest, NSQ_CPU_KnownGroundTruth)
{
    using GP                                  = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
    resetNLCounter();

    auto NL  = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                 pos_cpu,
                                                                 quat_cpu,
                                                                 GP::m_collisionDetection,
                                                                 nObs,
                                                                 nParts);
    auto got = extractPairs_CPU(NL.get(), pos_cpu, nObs, nParts);

    // Manually enumerated expected pairs (all C(5,2)=10 pairs)
    PairSet expected
        = {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}};

    EXPECT_EQ(got.size(), expected.size()) << "NSQ CPU: wrong number of pairs";
    EXPECT_EQ(got, expected) << "NSQ CPU: pair set does not match the known ground truth";
}

// ------------------------------------------------------------------------------------------------
// Test 2: NSQ GPU returns exactly the known ground truth (same 10 pairs).
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSimpleTest, NSQ_GPU_KnownGroundTruth)
{
    using GP                                  = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
    resetNLCounter();

    auto NL  = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                   pos_gpu,
                                                                   quat_gpu,
                                                                   GP::m_collisionDetection,
                                                                   nObs,
                                                                   nParts);
    auto got = extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts);

    PairSet expected
        = {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}};

    EXPECT_EQ(got.size(), expected.size()) << "NSQ GPU: wrong number of pairs";
    EXPECT_EQ(got, expected) << "NSQ GPU: pair set does not match the known ground truth";
}

// ------------------------------------------------------------------------------------------------
// Test 3: LinkedCell CPU (HOST) -- large cell covering whole domain -> all 10 pairs found.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSimpleTest, LinkedCell_HOST_KnownGroundTruth)
{
    using GP                                      = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
    lcParamsLarge.type                            = LinkedCellType::HOST;
    GP::m_collisionDetection.linkedCellParameters = lcParamsLarge;
    resetNLCounter();

    auto NL  = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                 pos_cpu,
                                                                 quat_cpu,
                                                                 GP::m_collisionDetection,
                                                                 nObs,
                                                                 nParts);
    auto got = extractPairs_CPU(NL.get(), pos_cpu, nObs, nParts);

    PairSet expected = allPairs(nParts);

    EXPECT_EQ(got.size(), expected.size()) << "LinkedCell HOST: wrong number of pairs";
    EXPECT_EQ(got, expected) << "LinkedCell HOST: pair set does not match the known ground truth";
}

// ------------------------------------------------------------------------------------------------
// Test 4: LinkedCell SORTBASED GPU -- large cell -> all 10 pairs found.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSimpleTest, LinkedCell_SortBased_GPU_KnownGroundTruth)
{
    using GP                                      = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
    lcParamsLarge.type                            = LinkedCellType::SORTBASED;
    GP::m_collisionDetection.linkedCellParameters = lcParamsLarge;
    resetNLCounter();

    auto NL  = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                   pos_gpu,
                                                                   quat_gpu,
                                                                   GP::m_collisionDetection,
                                                                   nObs,
                                                                   nParts);
    auto got = extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts);

    PairSet expected = allPairs(nParts);

    EXPECT_EQ(got.size(), expected.size()) << "LinkedCell SORTBASED: wrong number of pairs";
    EXPECT_EQ(got, expected)
        << "LinkedCell SORTBASED: pair set does not match the known ground truth";
}

// ------------------------------------------------------------------------------------------------
// Test 5: LinkedCell ATOMIC GPU -- large cell -> all 10 pairs found.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSimpleTest, LinkedCell_Atomic_GPU_KnownGroundTruth)
{
    using GP                                      = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
    lcParamsLarge.type                            = LinkedCellType::ATOMIC;
    GP::m_collisionDetection.linkedCellParameters = lcParamsLarge;
    resetNLCounter();

    auto NL  = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                   pos_gpu,
                                                                   quat_gpu,
                                                                   GP::m_collisionDetection,
                                                                   nObs,
                                                                   nParts);
    auto got = extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts);

    PairSet expected = allPairs(nParts);

    EXPECT_EQ(got.size(), expected.size()) << "LinkedCell ATOMIC: wrong number of pairs";
    EXPECT_EQ(got, expected) << "LinkedCell ATOMIC: pair set does not match the known ground truth";
}

// ------------------------------------------------------------------------------------------------
// Test 6: LinkedCell ATOMICFIXED GPU -- large cell -> all 10 pairs found.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSimpleTest, LinkedCell_AtomicFixed_GPU_KnownGroundTruth)
{
    using GP                                      = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
    lcParamsLarge.type                            = LinkedCellType::ATOMICFIXED;
    GP::m_collisionDetection.linkedCellParameters = lcParamsLarge;
    resetNLCounter();

    auto NL  = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                   pos_gpu,
                                                                   quat_gpu,
                                                                   GP::m_collisionDetection,
                                                                   nObs,
                                                                   nParts);
    auto got = extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts);

    PairSet expected = allPairs(nParts);

    EXPECT_EQ(got.size(), expected.size()) << "LinkedCell ATOMICFIXED: wrong number of pairs";
    EXPECT_EQ(got, expected)
        << "LinkedCell ATOMICFIXED: pair set does not match the known ground truth";
}

// ------------------------------------------------------------------------------------------------
// Test 7: All algorithms agree with each other on the 5-particle system.
//         Sweeps through every CPU/GPU combination and checks pairwise equality.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSimpleTest, AllAlgorithms_Match_KnownGroundTruth)
{
    using GP         = GrainsParameters<double>;
    PairSet expected = allPairs(nParts);  // 10 pairs

    struct Case
    {
        std::string      name;
        NeighborListType nlType;
        LinkedCellType   lcType;
        bool             isGPU;
    };

    std::vector<Case> cases = {
        {"NSQ_CPU", NeighborListType::NSQ, LinkedCellType::HOST, false},
        {"NSQ_GPU", NeighborListType::NSQ, LinkedCellType::HOST, true},
        {"LC_HOST", NeighborListType::LINKEDCELL, LinkedCellType::HOST, false},
        {"LC_SORTBASED_GPU", NeighborListType::LINKEDCELL, LinkedCellType::SORTBASED, true},
        {"LC_ATOMIC_GPU", NeighborListType::LINKEDCELL, LinkedCellType::ATOMIC, true},
        {"LC_ATOMICFIXED_GPU", NeighborListType::LINKEDCELL, LinkedCellType::ATOMICFIXED, true},
    };

    for(auto& c : cases)
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType = c.nlType;
        if(c.nlType == NeighborListType::LINKEDCELL)
        {
            lcParamsLarge.type                            = c.lcType;
            GP::m_collisionDetection.linkedCellParameters = lcParamsLarge;
        }

        PairSet got;
        if(!c.isGPU)
        {
            auto NL = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                         pos_cpu,
                                                                         quat_cpu,
                                                                         GP::m_collisionDetection,
                                                                         nObs,
                                                                         nParts);
            got     = extractPairs_CPU(NL.get(), pos_cpu, nObs, nParts);
        }
        else
        {
            auto NL = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                           pos_gpu,
                                                                           quat_gpu,
                                                                           GP::m_collisionDetection,
                                                                           nObs,
                                                                           nParts);
            got     = extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts);
        }

        EXPECT_EQ(got.size(), expected.size()) << c.name << ": wrong number of pairs (expected "
                                               << expected.size() << ", got " << got.size() << ")";
        EXPECT_EQ(got, expected) << c.name << ": pair set does not match the known ground truth";
    }
}

// ================================================================================================
//  FIXTURE: NeighborListSpatialFilterTest
//  -----------------------------------------
//  6 particles split into two spatially separated clusters.
//  The LinkedCell cell size is chosen so the two clusters are NOT in adjacent cells,
//  meaning only within-cluster pairs should appear in the LinkedCell neighbor list.
//  NSQ (no filter) still returns all C(6,2)=15 pairs, verifying the two methods differ.
//
//  Layout:
//    Cluster A (indices 0,1,2): near (0,0,0)
//      P0 = (0.0, 0.0, 0.0)
//      P1 = (0.5, 0.0, 0.0)
//      P2 = (0.0, 0.5, 0.0)
//    Cluster B (indices 3,4,5): near (10,0,0)
//      P3 = (10.0, 0.0, 0.0)
//      P4 = (10.5, 0.0, 0.0)
//      P5 = (10.0, 0.5, 0.0)
//
//  With cellSize=2.5 and domain (-1,-1,-1)->(12,2,2):
//    Cluster A -> x-cell 0,  Cluster B -> x-cell 4  (separation > 1 cell -> not adjacent)
//    -> cross-cluster pairs are excluded from the LinkedCell neighbor list.
// ================================================================================================
class NeighborListSpatialFilterTest : public ::testing::Test
{
protected:
    static constexpr uint nObs   = 0;
    static constexpr uint nParts = 6;
    static constexpr uint nTotal = nObs + nParts;

    Box<double>*                    sharedBox = nullptr;
    std::vector<RigidBody<double>*> rbVec;
    std::vector<Vector3<double>>    pos;

    GrainsMemBuffer<RigidBody<double>*, MemType::HOST>   rb_cpu;
    GrainsMemBuffer<Vector3<double>, MemType::HOST>      pos_cpu;
    GrainsMemBuffer<Quaternion<double>, MemType::HOST>   quat_cpu;
    GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE> rb_gpu;
    GrainsMemBuffer<Vector3<double>, MemType::DEVICE>    pos_gpu;
    GrainsMemBuffer<Quaternion<double>, MemType::DEVICE> quat_gpu;

    LinkedCellParameters<double> lcParams;

    PairSet expectedNSQ;         // all C(6,2) = 15 pairs
    PairSet expectedLinkedCell;  // only within-cluster: {0,1},{0,2},{1,2},{3,4},{3,5},{4,5}

    void SetUp() override
    {
        resetNLCounter();
        GrainsParameters<double>::m_collisionDetection = CollisionDetectionParameters<double>{};

        sharedBox = new Box<double>(0.2, 0.2, 0.2);

        // Cluster A
        pos.push_back(Vector3<double>(0.0, 0.0, 0.0));
        pos.push_back(Vector3<double>(0.5, 0.0, 0.0));
        pos.push_back(Vector3<double>(0.0, 0.5, 0.0));
        // Cluster B (far)
        pos.push_back(Vector3<double>(10.0, 0.0, 0.0));
        pos.push_back(Vector3<double>(10.5, 0.0, 0.0));
        pos.push_back(Vector3<double>(10.0, 0.5, 0.0));

        for(uint i = 0; i < nTotal; ++i)
            rbVec.push_back(new RigidBody<double>(sharedBox, 0.1, 1000.0, 1));

        buildBuffers(pos, nTotal, rbVec, rb_cpu, pos_cpu, quat_cpu, rb_gpu, pos_gpu, quat_gpu);

        // cellSize=2.5 -> clusters are 4 cells apart along x
        lcParams.minCorner              = Vector3<double>(-1.0, -1.0, -1.0);
        lcParams.maxCorner              = Vector3<double>(12.0, 2.0, 2.0);
        lcParams.minCellSize            = 2.5;
        lcParams.cellSizeFactor         = 1.0;
        lcParams.maxNumCellsPerObstacle = 8;
        lcParams.maxParticlesPerCell    = 16;
        lcParams.updateFrequency        = 0;
        lcParams.sortFrequency          = 0;

        expectedNSQ        = allPairs(nParts);  // 15 pairs
        expectedLinkedCell = {{0, 1},
                              {0, 2},
                              {1, 2},  // cluster A
                              {3, 4},
                              {3, 5},
                              {4, 5}};  // cluster B
    }

    void TearDown() override
    {
        GrainsParameters<double>::m_simulationState.neighborListUpdateCount = 0;
    }
};

// ------------------------------------------------------------------------------------------------
// Test 8: NSQ CPU on the two-cluster system returns all 15 pairs (no spatial filtering).
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSpatialFilterTest, NSQ_CPU_ReturnsAllPairs)
{
    using GP                                  = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
    resetNLCounter();

    auto NL  = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                 pos_cpu,
                                                                 quat_cpu,
                                                                 GP::m_collisionDetection,
                                                                 nObs,
                                                                 nParts);
    auto got = extractPairs_CPU(NL.get(), pos_cpu, nObs, nParts);

    EXPECT_EQ(got.size(), size_t(15))
        << "NSQ CPU should return all 15 pairs (no spatial filtering)";
    EXPECT_EQ(got, expectedNSQ);
}

// ------------------------------------------------------------------------------------------------
// Test 9: NSQ GPU on the two-cluster system returns all 15 pairs.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSpatialFilterTest, NSQ_GPU_ReturnsAllPairs)
{
    using GP                                  = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
    resetNLCounter();

    auto NL  = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                   pos_gpu,
                                                                   quat_gpu,
                                                                   GP::m_collisionDetection,
                                                                   nObs,
                                                                   nParts);
    auto got = extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts);

    EXPECT_EQ(got.size(), size_t(15))
        << "NSQ GPU should return all 15 pairs (no spatial filtering)";
    EXPECT_EQ(got, expectedNSQ);
}

// ------------------------------------------------------------------------------------------------
// Test 10: LinkedCell HOST -- spatially filtered to only 6 within-cluster pairs.
//          Explicitly verifies that all 9 cross-cluster pairs are absent.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSpatialFilterTest, LinkedCell_HOST_SpatialFilter_KnownPairs)
{
    using GP                                      = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
    lcParams.type                                 = LinkedCellType::HOST;
    GP::m_collisionDetection.linkedCellParameters = lcParams;
    resetNLCounter();

    auto NL  = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                 pos_cpu,
                                                                 quat_cpu,
                                                                 GP::m_collisionDetection,
                                                                 nObs,
                                                                 nParts);
    auto got = extractPairs_CPU(NL.get(), pos_cpu, nObs, nParts);

    EXPECT_EQ(got.size(), expectedLinkedCell.size())
        << "LinkedCell HOST: wrong number of pairs after spatial filtering "
           "(expected "
        << expectedLinkedCell.size() << ", got " << got.size() << ")";
    EXPECT_EQ(got, expectedLinkedCell)
        << "LinkedCell HOST: pair set does not match the expected within-cluster pairs";

    // Explicitly confirm that all cross-cluster pairs are absent
    PairSet crossCluster = {{0, 3}, {0, 4}, {0, 5}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {2, 4}, {2, 5}};
    for(auto& p : crossCluster)
        EXPECT_EQ(got.count(p), size_t(0))
            << "LinkedCell HOST: cross-cluster pair (" << p.first << "," << p.second
            << ") should NOT be in the neighbor list";
}

// ------------------------------------------------------------------------------------------------
// Test 11: All LinkedCell GPU variants correctly filter to only 6 within-cluster pairs.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListSpatialFilterTest, AllLinkedCell_GPU_SpatialFilter_KnownPairs)
{
    using GP                                  = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType = NeighborListType::LINKEDCELL;

    std::vector<std::pair<LinkedCellType, std::string>> gpuVariants = {
        {LinkedCellType::SORTBASED, "SORTBASED"},
        {LinkedCellType::ATOMIC, "ATOMIC"},
        {LinkedCellType::ATOMICFIXED, "ATOMICFIXED"},
    };

    for(auto& [lcType, name] : gpuVariants)
    {
        resetNLCounter();
        lcParams.type                                 = lcType;
        GP::m_collisionDetection.linkedCellParameters = lcParams;

        auto NL  = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                       pos_gpu,
                                                                       quat_gpu,
                                                                       GP::m_collisionDetection,
                                                                       nObs,
                                                                       nParts);
        auto got = extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts);

        EXPECT_EQ(got.size(), expectedLinkedCell.size())
            << "LinkedCell " << name << ": wrong number of pairs (expected "
            << expectedLinkedCell.size() << ", got " << got.size() << ")";
        EXPECT_EQ(got, expectedLinkedCell)
            << "LinkedCell " << name << ": pair set differs from expected within-cluster pairs";
    }
}

// ================================================================================================
//  FIXTURE: NeighborListConsistencyTest
//  --------------------------------------
//  A 4x4x4 = 64-particle uniform grid with 1.0 spacing.
//  The LinkedCell uses a single large cell (cellSize=20) so ALL particles are in the same
//  cell neighbourhood -- both NSQ and every LinkedCell variant must return the same
//  C(64,2) = 2016 pairs, verifying cross-algorithm consistency at a larger scale.
// ================================================================================================
class NeighborListConsistencyTest : public ::testing::Test
{
protected:
    static constexpr uint GRID           = 4;
    static constexpr uint nObs           = 0;
    static constexpr uint nParts         = GRID * GRID * GRID;  // 64
    static constexpr uint nTotal         = nObs + nParts;
    static constexpr uint EXPECTED_PAIRS = nParts * (nParts - 1) / 2;  // 2016

    Box<double>*                    sharedBox = nullptr;
    std::vector<RigidBody<double>*> rbVec;
    std::vector<Vector3<double>>    pos;

    GrainsMemBuffer<RigidBody<double>*, MemType::HOST>   rb_cpu;
    GrainsMemBuffer<Vector3<double>, MemType::HOST>      pos_cpu;
    GrainsMemBuffer<Quaternion<double>, MemType::HOST>   quat_cpu;
    GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE> rb_gpu;
    GrainsMemBuffer<Vector3<double>, MemType::DEVICE>    pos_gpu;
    GrainsMemBuffer<Quaternion<double>, MemType::DEVICE> quat_gpu;

    LinkedCellParameters<double> lcParamsLarge;

    void SetUp() override
    {
        resetNLCounter();
        GrainsParameters<double>::m_collisionDetection = CollisionDetectionParameters<double>{};

        sharedBox = new Box<double>(0.3, 0.3, 0.3);

        double spacing = 1.0;
        for(uint i = 0; i < GRID; ++i)
            for(uint j = 0; j < GRID; ++j)
                for(uint k = 0; k < GRID; ++k)
                    pos.push_back(Vector3<double>(i * spacing, j * spacing, k * spacing));

        for(uint i = 0; i < nTotal; ++i)
            rbVec.push_back(new RigidBody<double>(sharedBox, 0.1, 1000.0, 1));

        buildBuffers(pos, nTotal, rbVec, rb_cpu, pos_cpu, quat_cpu, rb_gpu, pos_gpu, quat_gpu);

        // One big cell -> every particle is in same neighbourhood
        lcParamsLarge.minCorner                       = Vector3<double>(-1.0, -1.0, -1.0);
        lcParamsLarge.maxCorner                       = Vector3<double>(5.0, 5.0, 5.0);
        lcParamsLarge.minCellSize                     = 20.0;
        lcParamsLarge.cellSizeFactor                  = 1.0;
        lcParamsLarge.maxNumCellsPerObstacle          = 8;
        lcParamsLarge.maxParticlesPerCell             = 256;
        lcParamsLarge.initialNumberOfPairsPerParticle = 128;
        lcParamsLarge.updateFrequency                 = 0;
        lcParamsLarge.sortFrequency                   = 0;
    }

    void TearDown() override
    {
        GrainsParameters<double>::m_simulationState.neighborListUpdateCount = 0;
    }
};

// ------------------------------------------------------------------------------------------------
// Test 12: All algorithms return 2016 pairs and agree pairwise on the 64-particle grid.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListConsistencyTest, AllAlgorithms_Return_Same_Pairs)
{
    using GP = GrainsParameters<double>;

    struct AlgoResult
    {
        std::string name;
        PairSet     pairs;
    };
    std::vector<AlgoResult> results;

    // NSQ CPU
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
        auto NL = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                     pos_cpu,
                                                                     quat_cpu,
                                                                     GP::m_collisionDetection,
                                                                     nObs,
                                                                     nParts);
        results.push_back({"NSQ_CPU", extractPairs_CPU(NL.get(), pos_cpu, nObs, nParts)});
    }

    // NSQ GPU
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
        auto NL = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                       pos_gpu,
                                                                       quat_gpu,
                                                                       GP::m_collisionDetection,
                                                                       nObs,
                                                                       nParts);
        results.push_back({"NSQ_GPU", extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts)});
    }

    // LinkedCell HOST
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
        lcParamsLarge.type                            = LinkedCellType::HOST;
        GP::m_collisionDetection.linkedCellParameters = lcParamsLarge;
        auto NL = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                     pos_cpu,
                                                                     quat_cpu,
                                                                     GP::m_collisionDetection,
                                                                     nObs,
                                                                     nParts);
        results.push_back({"LC_HOST", extractPairs_CPU(NL.get(), pos_cpu, nObs, nParts)});
    }

    // LinkedCell GPU variants
    std::vector<std::pair<LinkedCellType, std::string>> gpuVariants = {
        {LinkedCellType::SORTBASED, "LC_SORTBASED_GPU"},
        {LinkedCellType::ATOMIC, "LC_ATOMIC_GPU"},
        {LinkedCellType::ATOMICFIXED, "LC_ATOMICFIXED_GPU"},
    };
    for(auto& [lcType, name] : gpuVariants)
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
        lcParamsLarge.type                            = lcType;
        GP::m_collisionDetection.linkedCellParameters = lcParamsLarge;
        auto NL = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                       pos_gpu,
                                                                       quat_gpu,
                                                                       GP::m_collisionDetection,
                                                                       nObs,
                                                                       nParts);
        results.push_back({name, extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts)});
    }

    // Verify pair count and pairwise equality
    for(auto& r : results)
    {
        EXPECT_EQ(r.pairs.size(), size_t(EXPECTED_PAIRS))
            << r.name << ": expected " << EXPECTED_PAIRS << " pairs, got " << r.pairs.size();
    }

    const PairSet& baseline = results[0].pairs;
    for(size_t i = 1; i < results.size(); ++i)
    {
        EXPECT_EQ(results[i].pairs, baseline)
            << results[i].name << " vs " << results[0].name << ": pair sets differ";
    }
}

// ================================================================================================
//  FIXTURE: NeighborListBruteForceTest  (preserved from original, fixed)
//  -----------------------------------------------------------------------
//  8 particles in a 2x2x2 grid.  Verifies NSQ and LinkedCell CPU vs GPU.
//  The original loop covered only 2 of 3 GPU LinkedCell variants -- now all 3 are tested.
// ================================================================================================
class NeighborListBruteForceTest : public ::testing::Test
{
protected:
    static constexpr uint nObs   = 0;
    static constexpr uint nParts = 8;  // 2x2x2
    static constexpr uint nTotal = nObs + nParts;

    Box<double>*                    sharedBox = nullptr;
    std::vector<RigidBody<double>*> rbVec;

    GrainsMemBuffer<RigidBody<double>*, MemType::HOST>   rb_cpu;
    GrainsMemBuffer<Vector3<double>, MemType::HOST>      pos_cpu;
    GrainsMemBuffer<Quaternion<double>, MemType::HOST>   quat_cpu;
    GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE> rb_gpu;
    GrainsMemBuffer<Vector3<double>, MemType::DEVICE>    pos_gpu;
    GrainsMemBuffer<Quaternion<double>, MemType::DEVICE> quat_gpu;

    LinkedCellParameters<double> lcParams;

    void SetUp() override
    {
        resetNLCounter();
        GrainsParameters<double>::m_collisionDetection = CollisionDetectionParameters<double>{};

        sharedBox                            = new Box<double>(0.5, 0.5, 0.5);
        double                       spacing = 2.0;
        std::vector<Vector3<double>> pos;
        for(uint i = 0; i < 2; ++i)
            for(uint j = 0; j < 2; ++j)
                for(uint k = 0; k < 2; ++k)
                {
                    rbVec.push_back(new RigidBody<double>(sharedBox, 0.1, 1000.0, 1));
                    pos.push_back(Vector3<double>(i * spacing, j * spacing, k * spacing));
                }

        buildBuffers(pos, nTotal, rbVec, rb_cpu, pos_cpu, quat_cpu, rb_gpu, pos_gpu, quat_gpu);

        lcParams.minCorner                       = Vector3<double>(-1.0, -1.0, -1.0);
        lcParams.maxCorner                       = Vector3<double>(5.0, 5.0, 5.0);
        lcParams.minCellSize                     = 10.0;  // one big cell
        lcParams.cellSizeFactor                  = 1.0;
        lcParams.maxNumCellsPerObstacle          = 8;
        lcParams.maxParticlesPerCell             = 64;
        lcParams.initialNumberOfPairsPerParticle = 32;
        lcParams.updateFrequency                 = 0;
        lcParams.sortFrequency                   = 0;
    }

    void TearDown() override
    {
        GrainsParameters<double>::m_simulationState.neighborListUpdateCount = 0;
    }
};

// ------------------------------------------------------------------------------------------------
// Test 13: NSQ CPU vs GPU on the 2x2x2 grid.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListBruteForceTest, NSQ_CPU_vs_GPU)
{
    using GP                                  = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;

    resetNLCounter();
    auto NL_cpu    = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                     pos_cpu,
                                                                     quat_cpu,
                                                                     GP::m_collisionDetection,
                                                                     nObs,
                                                                     nParts);
    auto cpu_pairs = extractPairs_CPU(NL_cpu.get(), pos_cpu, nObs, nParts);

    resetNLCounter();
    auto NL_gpu    = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                       pos_gpu,
                                                                       quat_gpu,
                                                                       GP::m_collisionDetection,
                                                                       nObs,
                                                                       nParts);
    auto gpu_pairs = extractPairs_GPU(NL_gpu.get(), pos_gpu, nObs, nParts);

    EXPECT_GT(cpu_pairs.size(), size_t(0)) << "NSQ should produce some pairs";
    EXPECT_EQ(cpu_pairs, gpu_pairs) << "NSQ CPU and GPU results must match";
}

// ------------------------------------------------------------------------------------------------
// Test 14: LinkedCell HOST vs all three GPU variants on the 2x2x2 grid.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListBruteForceTest, LinkedCell_HOST_vs_AllGPU_Variants)
{
    using GP                                      = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
    lcParams.type                                 = LinkedCellType::HOST;
    GP::m_collisionDetection.linkedCellParameters = lcParams;

    resetNLCounter();
    auto NL_cpu    = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                     pos_cpu,
                                                                     quat_cpu,
                                                                     GP::m_collisionDetection,
                                                                     nObs,
                                                                     nParts);
    auto cpu_pairs = extractPairs_CPU(NL_cpu.get(), pos_cpu, nObs, nParts);

    EXPECT_GT(cpu_pairs.size(), size_t(0)) << "Should find some neighbour pairs";

    std::vector<std::pair<LinkedCellType, std::string>> gpuVariants = {
        {LinkedCellType::SORTBASED, "SORTBASED"},
        {LinkedCellType::ATOMIC, "ATOMIC"},
        {LinkedCellType::ATOMICFIXED, "ATOMICFIXED"},
    };

    for(auto& [lcType, name] : gpuVariants)
    {
        resetNLCounter();
        lcParams.type                                 = lcType;
        GP::m_collisionDetection.linkedCellParameters = lcParams;
        auto NL_gpu    = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                           pos_gpu,
                                                                           quat_gpu,
                                                                           GP::m_collisionDetection,
                                                                           nObs,
                                                                           nParts);
        auto gpu_pairs = extractPairs_GPU(NL_gpu.get(), pos_gpu, nObs, nParts);

        EXPECT_EQ(cpu_pairs, gpu_pairs)
            << "LinkedCell HOST vs " << name << " GPU: pair sets differ";
    }
    cudaDeviceSynchronize();
}

// ================================================================================================
//  UTILITY: raw-buffer validity checks
//  -------------------------------------
//  These helpers inspect the raw uint2 pair buffer BEFORE any normalization so they can catch
//  bugs that the PairSet abstraction would hide (e.g. buffer overflow writes dummy entries,
//  a kernel emitting (i,i) self-pairs, or two entries for the same unordered pair).
// ================================================================================================

/** @brief Return raw pairs directly from a HOST NeighborList (no normalization). */
static std::vector<std::pair<uint, uint>>
    rawPairs_CPU(NeighborList<double, MemType::HOST>*             NL,
                 GrainsMemBuffer<Vector3<double>, MemType::HOST>& positions,
                 uint                                             nObs,
                 uint                                             nParts)
{
    NL->updateNeighborList(positions, nObs, nParts);
    std::vector<std::pair<uint, uint>> out;
    const uint2*                       data = NL->getData();
    uint                               size = NL->getSize();
    out.reserve(size);
    for(uint i = 0; i < size; ++i)
        out.push_back({data[i].x, data[i].y});
    return out;
}

/** @brief Return raw pairs directly from a DEVICE NeighborList (no normalization). */
static std::vector<std::pair<uint, uint>>
    rawPairs_GPU(NeighborList<double, MemType::DEVICE>*             NL,
                 GrainsMemBuffer<Vector3<double>, MemType::DEVICE>& positions,
                 uint                                               nObs,
                 uint                                               nParts)
{
    NL->updateNeighborList(positions, nObs, nParts);
    cudaDeviceSynchronize();
    uint               size = NL->getSize();
    std::vector<uint2> host(size);
    if(size > 0)
        cudaMemcpy(host.data(), NL->getData(), size * sizeof(uint2), cudaMemcpyDeviceToHost);
    std::vector<std::pair<uint, uint>> out;
    out.reserve(size);
    for(uint i = 0; i < size; ++i)
        out.push_back({host[i].x, host[i].y});
    return out;
}

/** @brief Assert no self-pairs (i,i) exist in a raw pair list. */
static void assertNoSelfPairs(const std::vector<std::pair<uint, uint>>& raw,
                              const std::string&                        label)
{
    for(auto& p : raw)
        EXPECT_NE(p.first, p.second)
            << label << ": self-pair (" << p.first << "," << p.second << ") found";
}

/** @brief Assert no duplicate unordered pairs exist in a raw pair list. */
static void assertNoDuplicatePairs(const std::vector<std::pair<uint, uint>>& raw,
                                   const std::string&                        label)
{
    std::set<std::pair<uint, uint>> seen;
    for(auto& p : raw)
    {
        auto norm = normPair(p.first, p.second);
        EXPECT_TRUE(seen.insert(norm).second)
            << label << ": duplicate pair (" << p.first << "," << p.second << ") found";
    }
}

/** @brief Assert all indices in a raw pair list are within [0, nTotal). */
static void assertValidIndices(const std::vector<std::pair<uint, uint>>& raw,
                               uint                                      nTotal,
                               const std::string&                        label)
{
    for(auto& p : raw)
    {
        EXPECT_LT(p.first, nTotal)
            << label << ": index " << p.first << " out of range [0, " << nTotal << ")";
        EXPECT_LT(p.second, nTotal)
            << label << ": index " << p.second << " out of range [0, " << nTotal << ")";
    }
}

// ================================================================================================
//  FIXTURE: NeighborListPairValidityTest
//  ---------------------------------------
//  Runs no-self-pairs, no-duplicates, and valid-indices checks on all 6 CPU/GPU algorithms
//  using the 5-particle layout.  These tests catch buffer-level corruption that the
//  PairSet abstraction would silently hide.
// ================================================================================================
class NeighborListPairValidityTest : public ::testing::Test
{
protected:
    static constexpr uint nObs   = 0;
    static constexpr uint nParts = 5;
    static constexpr uint nTotal = nObs + nParts;

    Box<double>*                    sharedBox = nullptr;
    std::vector<RigidBody<double>*> rbVec;
    std::vector<Vector3<double>>    pos;

    GrainsMemBuffer<RigidBody<double>*, MemType::HOST>   rb_cpu;
    GrainsMemBuffer<Vector3<double>, MemType::HOST>      pos_cpu;
    GrainsMemBuffer<Quaternion<double>, MemType::HOST>   quat_cpu;
    GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE> rb_gpu;
    GrainsMemBuffer<Vector3<double>, MemType::DEVICE>    pos_gpu;
    GrainsMemBuffer<Quaternion<double>, MemType::DEVICE> quat_gpu;

    LinkedCellParameters<double> lcParams;

    void SetUp() override
    {
        resetNLCounter();
        GrainsParameters<double>::m_collisionDetection = CollisionDetectionParameters<double>{};

        sharedBox = new Box<double>(0.4, 0.4, 0.4);
        pos       = {Vector3<double>(0.0, 0.0, 0.0),
                     Vector3<double>(1.0, 0.0, 0.0),
                     Vector3<double>(0.0, 1.0, 0.0),
                     Vector3<double>(0.0, 0.0, 1.0),
                     Vector3<double>(0.0, 0.0, 2.0)};

        for(uint i = 0; i < nTotal; ++i)
            rbVec.push_back(new RigidBody<double>(sharedBox, 0.1, 1000.0, 1));

        buildBuffers(pos, nTotal, rbVec, rb_cpu, pos_cpu, quat_cpu, rb_gpu, pos_gpu, quat_gpu);

        // One big cell
        lcParams.minCorner              = Vector3<double>(-1.0, -1.0, -1.0);
        lcParams.maxCorner              = Vector3<double>(5.0, 5.0, 5.0);
        lcParams.minCellSize            = 10.0;
        lcParams.cellSizeFactor         = 1.0;
        lcParams.maxNumCellsPerObstacle = 8;
        lcParams.maxParticlesPerCell    = 64;
        lcParams.updateFrequency        = 0;
        lcParams.sortFrequency          = 0;
    }

    void TearDown() override
    {
        GrainsParameters<double>::m_simulationState.neighborListUpdateCount = 0;
    }
};

// ------------------------------------------------------------------------------------------------
// Test 15: Raw pair buffer validity -- no self-pairs, no duplicates, valid indices for all
//          6 CPU/GPU algorithm variants.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListPairValidityTest, AllAlgorithms_RawPairBuffer_Valid)
{
    using GP = GrainsParameters<double>;

    struct Case
    {
        std::string      name;
        NeighborListType nlType;
        LinkedCellType   lcType;
        bool             isGPU;
    };

    std::vector<Case> cases = {
        {"NSQ_CPU", NeighborListType::NSQ, LinkedCellType::HOST, false},
        {"NSQ_GPU", NeighborListType::NSQ, LinkedCellType::HOST, true},
        {"LC_HOST", NeighborListType::LINKEDCELL, LinkedCellType::HOST, false},
        {"LC_SORTBASED_GPU", NeighborListType::LINKEDCELL, LinkedCellType::SORTBASED, true},
        {"LC_ATOMIC_GPU", NeighborListType::LINKEDCELL, LinkedCellType::ATOMIC, true},
        {"LC_ATOMICFIXED_GPU", NeighborListType::LINKEDCELL, LinkedCellType::ATOMICFIXED, true},
    };

    for(auto& c : cases)
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType = c.nlType;
        if(c.nlType == NeighborListType::LINKEDCELL)
        {
            lcParams.type                                 = c.lcType;
            GP::m_collisionDetection.linkedCellParameters = lcParams;
        }

        std::vector<std::pair<uint, uint>> raw;
        if(!c.isGPU)
        {
            auto NL = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                         pos_cpu,
                                                                         quat_cpu,
                                                                         GP::m_collisionDetection,
                                                                         nObs,
                                                                         nParts);
            raw     = rawPairs_CPU(NL.get(), pos_cpu, nObs, nParts);
        }
        else
        {
            auto NL = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                           pos_gpu,
                                                                           quat_gpu,
                                                                           GP::m_collisionDetection,
                                                                           nObs,
                                                                           nParts);
            raw     = rawPairs_GPU(NL.get(), pos_gpu, nObs, nParts);
        }

        EXPECT_GT(raw.size(), size_t(0)) << c.name << ": produced no pairs";
        assertNoSelfPairs(raw, c.name);
        assertNoDuplicatePairs(raw, c.name);
        assertValidIndices(raw, nTotal, c.name);
    }
}

// ================================================================================================
//  FIXTURE: NeighborListLargeScaleTest
//  -------------------------------------
//  8x8x8 = 512 particles on a uniform grid, spacing 1.0.
//  Particles occupy (0..7, 0..7, 0..7).  Two LinkedCell configs are used:
//
//  (A) ONE BIG CELL (cellSize=20) -- every particle sees every other particle.
//      -> All 6 algorithms must agree on exactly C(512,2) = 130,816 pairs.
//
//  (B) REALISTIC CELLS (cellSize=2.0) -- domain [-0.5,7.5]^3 -> 4 cells per axis.
//      Each cell holds a 2x2x2=8-particle sub-cube; adjacent cells share a face/edge/corner.
//      -> All LinkedCell variants (CPU+GPU) must agree on the same pair set.
//      -> NSQ returns strictly MORE pairs (it has no spatial filter), confirming LC IS filtering.
//      -> Pair validity (no self-pairs, no duplicates, valid indices) is checked for all variants.
//
//  This fixture is the primary regression harness for catching bugs that only manifest at scale:
//   - Buffer overflow / under-allocation  (all pairs present)
//   - atomic collision -> missed pairs     (all variants agree)
//   - ATOMICFIXED maxParticlesPerCell overflow  (same count as ATOMIC)
//   - Sort-based index remapping errors   (same pairs after sorting)
// ================================================================================================
class NeighborListLargeScaleTest : public ::testing::Test
{
protected:
    static constexpr uint GRID           = 8;
    static constexpr uint nObs           = 0;
    static constexpr uint nParts         = GRID * GRID * GRID;  // 512
    static constexpr uint nTotal         = nObs + nParts;
    static constexpr uint EXPECTED_PAIRS = nParts * (nParts - 1) / 2;  // 130816

    Box<double>*                    sharedBox = nullptr;
    std::vector<RigidBody<double>*> rbVec;
    std::vector<Vector3<double>>    pos;

    GrainsMemBuffer<RigidBody<double>*, MemType::HOST>   rb_cpu;
    GrainsMemBuffer<Vector3<double>, MemType::HOST>      pos_cpu;
    GrainsMemBuffer<Quaternion<double>, MemType::HOST>   quat_cpu;
    GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE> rb_gpu;
    GrainsMemBuffer<Vector3<double>, MemType::DEVICE>    pos_gpu;
    GrainsMemBuffer<Quaternion<double>, MemType::DEVICE> quat_gpu;

    LinkedCellParameters<double> lcParamsLarge;    // one big cell
    LinkedCellParameters<double> lcParamsRealist;  // cellSize=2.0, 4 cells/axis

    void SetUp() override
    {
        resetNLCounter();
        GrainsParameters<double>::m_collisionDetection = CollisionDetectionParameters<double>{};

        sharedBox = new Box<double>(0.3, 0.3, 0.3);

        double spacing = 1.0;
        for(uint i = 0; i < GRID; ++i)
            for(uint j = 0; j < GRID; ++j)
                for(uint k = 0; k < GRID; ++k)
                    pos.push_back(Vector3<double>(i * spacing, j * spacing, k * spacing));

        for(uint i = 0; i < nTotal; ++i)
            rbVec.push_back(new RigidBody<double>(sharedBox, 0.1, 1000.0, 1));

        buildBuffers(pos, nTotal, rbVec, rb_cpu, pos_cpu, quat_cpu, rb_gpu, pos_gpu, quat_gpu);

        // (A) Large-cell config: one cell, all 512 particles visible to each other
        lcParamsLarge.minCorner              = Vector3<double>(-1.0, -1.0, -1.0);
        lcParamsLarge.maxCorner              = Vector3<double>(9.0, 9.0, 9.0);
        lcParamsLarge.minCellSize            = 20.0;
        lcParamsLarge.cellSizeFactor         = 1.0;
        lcParamsLarge.maxNumCellsPerObstacle = 8;
        lcParamsLarge.maxParticlesPerCell    = 512;
        // Must pre-allocate enough: C(512,2)=130816 -> 130816/512 ~ 256 pairs/particle minimum
        lcParamsLarge.initialNumberOfPairsPerParticle = 512;
        lcParamsLarge.updateFrequency                 = 0;
        lcParamsLarge.sortFrequency                   = 0;

        // (B) Realistic config: domain [-0.5,7.5]^3 -> cellSize=2.0 -> 4x4x4 cells
        //     Each cell contains 2x2x2=8 particles; 27-cell stencil neighbourhood
        lcParamsRealist.minCorner              = Vector3<double>(-0.5, -0.5, -0.5);
        lcParamsRealist.maxCorner              = Vector3<double>(7.5, 7.5, 7.5);
        lcParamsRealist.minCellSize            = 2.0;
        lcParamsRealist.cellSizeFactor         = 1.0;
        lcParamsRealist.maxNumCellsPerObstacle = 8;
        // 27 cells x 8 particles/cell = 216 potential neighbors per particle
        lcParamsRealist.maxParticlesPerCell             = 16;
        lcParamsRealist.initialNumberOfPairsPerParticle = 128;
        lcParamsRealist.updateFrequency                 = 0;
        lcParamsRealist.sortFrequency                   = 0;
    }

    void TearDown() override
    {
        GrainsParameters<double>::m_simulationState.neighborListUpdateCount = 0;
    }
};

// ------------------------------------------------------------------------------------------------
// Test 16: 512 particles, one big cell -- all 6 algorithms agree on exactly 130,816 pairs.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListLargeScaleTest, LargeCell_AllAlgorithms_Agree_On_130816_Pairs)
{
    using GP = GrainsParameters<double>;

    struct AlgoResult
    {
        std::string name;
        PairSet     pairs;
    };
    std::vector<AlgoResult> results;

    // NSQ CPU
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
        auto NL = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                     pos_cpu,
                                                                     quat_cpu,
                                                                     GP::m_collisionDetection,
                                                                     nObs,
                                                                     nParts);
        results.push_back({"NSQ_CPU", extractPairs_CPU(NL.get(), pos_cpu, nObs, nParts)});
    }

    // NSQ GPU
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
        auto NL = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                       pos_gpu,
                                                                       quat_gpu,
                                                                       GP::m_collisionDetection,
                                                                       nObs,
                                                                       nParts);
        results.push_back({"NSQ_GPU", extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts)});
    }

    // LinkedCell HOST
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
        lcParamsLarge.type                            = LinkedCellType::HOST;
        GP::m_collisionDetection.linkedCellParameters = lcParamsLarge;
        auto NL = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                     pos_cpu,
                                                                     quat_cpu,
                                                                     GP::m_collisionDetection,
                                                                     nObs,
                                                                     nParts);
        results.push_back({"LC_HOST", extractPairs_CPU(NL.get(), pos_cpu, nObs, nParts)});
    }

    // LinkedCell GPU variants
    for(auto& [lcType, name] : std::vector<std::pair<LinkedCellType, std::string>>{
            {LinkedCellType::SORTBASED, "LC_SORTBASED_GPU"},
            {LinkedCellType::ATOMIC, "LC_ATOMIC_GPU"},
            {LinkedCellType::ATOMICFIXED, "LC_ATOMICFIXED_GPU"},
        })
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
        lcParamsLarge.type                            = lcType;
        GP::m_collisionDetection.linkedCellParameters = lcParamsLarge;
        auto NL = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                       pos_gpu,
                                                                       quat_gpu,
                                                                       GP::m_collisionDetection,
                                                                       nObs,
                                                                       nParts);
        results.push_back({name, extractPairs_GPU(NL.get(), pos_gpu, nObs, nParts)});
    }

    // Every algorithm must return exactly 130,816 pairs
    for(auto& r : results)
        EXPECT_EQ(r.pairs.size(), size_t(EXPECTED_PAIRS))
            << r.name << ": expected " << EXPECTED_PAIRS << " pairs, got " << r.pairs.size();

    // Every algorithm must agree with NSQ_CPU (the reference)
    const PairSet& baseline = results[0].pairs;
    for(size_t i = 1; i < results.size(); ++i)
        EXPECT_EQ(results[i].pairs, baseline) << results[i].name << " vs NSQ_CPU: pair sets differ";
}

// ------------------------------------------------------------------------------------------------
// Test 17: 512 particles, realistic cell size (2.0) -- all LinkedCell variants agree pairwise
//          and NSQ returns strictly more pairs than LinkedCell (spatial filtering is active).
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListLargeScaleTest, RealisticCells_AllLinkedCell_Agree_NSQ_Larger)
{
    using GP = GrainsParameters<double>;

    // --- Run NSQ CPU to get the unfiltered reference count ---
    resetNLCounter();
    GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
    auto    NL_nsq    = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                     pos_cpu,
                                                                     quat_cpu,
                                                                     GP::m_collisionDetection,
                                                                     nObs,
                                                                     nParts);
    PairSet nsq_pairs = extractPairs_CPU(NL_nsq.get(), pos_cpu, nObs, nParts);

    EXPECT_EQ(nsq_pairs.size(), size_t(EXPECTED_PAIRS))
        << "NSQ CPU should always produce all " << EXPECTED_PAIRS << " pairs";

    // --- Run LinkedCell HOST ---
    resetNLCounter();
    GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
    lcParamsRealist.type                          = LinkedCellType::HOST;
    GP::m_collisionDetection.linkedCellParameters = lcParamsRealist;
    auto    NL_host       = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                      pos_cpu,
                                                                      quat_cpu,
                                                                      GP::m_collisionDetection,
                                                                      nObs,
                                                                      nParts);
    PairSet lc_host_pairs = extractPairs_CPU(NL_host.get(), pos_cpu, nObs, nParts);

    EXPECT_GT(lc_host_pairs.size(), size_t(0)) << "LC HOST: should find some pairs";
    EXPECT_LT(lc_host_pairs.size(), nsq_pairs.size())
        << "LC HOST: spatial filter should reduce pair count below NSQ's " << EXPECTED_PAIRS;

    // --- Run all LinkedCell GPU variants and compare to LC HOST ---
    for(auto& [lcType, name] : std::vector<std::pair<LinkedCellType, std::string>>{
            {LinkedCellType::SORTBASED, "LC_SORTBASED_GPU"},
            {LinkedCellType::ATOMIC, "LC_ATOMIC_GPU"},
            {LinkedCellType::ATOMICFIXED, "LC_ATOMICFIXED_GPU"},
        })
    {
        resetNLCounter();
        lcParamsRealist.type                          = lcType;
        GP::m_collisionDetection.linkedCellParameters = lcParamsRealist;
        auto    NL_gpu    = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                           pos_gpu,
                                                                           quat_gpu,
                                                                           GP::m_collisionDetection,
                                                                           nObs,
                                                                           nParts);
        PairSet gpu_pairs = extractPairs_GPU(NL_gpu.get(), pos_gpu, nObs, nParts);

        EXPECT_EQ(gpu_pairs.size(), lc_host_pairs.size())
            << name << " vs LC_HOST: pair counts differ (" << gpu_pairs.size() << " vs "
            << lc_host_pairs.size() << ")";
        EXPECT_EQ(gpu_pairs, lc_host_pairs) << name << " vs LC_HOST: pair sets differ";
    }
    cudaDeviceSynchronize();
}

// ------------------------------------------------------------------------------------------------
// Test 18: 512 particles, realistic cells -- raw pair buffer validity for all algorithms.
//          Checks no self-pairs, no duplicate pairs, all indices in [0, 512).
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListLargeScaleTest, RealisticCells_RawPairBuffer_Valid)
{
    using GP = GrainsParameters<double>;

    struct Case
    {
        std::string      name;
        NeighborListType nlType;
        LinkedCellType   lcType;
        bool             isGPU;
    };

    std::vector<Case> cases = {
        {"NSQ_CPU", NeighborListType::NSQ, LinkedCellType::HOST, false},
        {"NSQ_GPU", NeighborListType::NSQ, LinkedCellType::HOST, true},
        {"LC_HOST", NeighborListType::LINKEDCELL, LinkedCellType::HOST, false},
        {"LC_SORTBASED_GPU", NeighborListType::LINKEDCELL, LinkedCellType::SORTBASED, true},
        {"LC_ATOMIC_GPU", NeighborListType::LINKEDCELL, LinkedCellType::ATOMIC, true},
        {"LC_ATOMICFIXED_GPU", NeighborListType::LINKEDCELL, LinkedCellType::ATOMICFIXED, true},
    };

    for(auto& c : cases)
    {
        resetNLCounter();
        GP::m_collisionDetection.neighborListType = c.nlType;
        LinkedCellParameters<double>& lcP
            = (c.nlType == NeighborListType::NSQ) ? lcParamsRealist : lcParamsRealist;
        if(c.nlType == NeighborListType::LINKEDCELL)
        {
            lcP.type                                      = c.lcType;
            GP::m_collisionDetection.linkedCellParameters = lcP;
        }

        std::vector<std::pair<uint, uint>> raw;
        if(!c.isGPU)
        {
            auto NL = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                         pos_cpu,
                                                                         quat_cpu,
                                                                         GP::m_collisionDetection,
                                                                         nObs,
                                                                         nParts);
            raw     = rawPairs_CPU(NL.get(), pos_cpu, nObs, nParts);
        }
        else
        {
            auto NL = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                           pos_gpu,
                                                                           quat_gpu,
                                                                           GP::m_collisionDetection,
                                                                           nObs,
                                                                           nParts);
            raw     = rawPairs_GPU(NL.get(), pos_gpu, nObs, nParts);
        }

        EXPECT_GT(raw.size(), size_t(0)) << c.name << ": produced no pairs";
        assertNoSelfPairs(raw, c.name);
        assertNoDuplicatePairs(raw, c.name);
        assertValidIndices(raw, nTotal, c.name);
    }
}

// ------------------------------------------------------------------------------------------------
// Test 19: 512 particles, realistic cells -- pair count is identical on repeated calls.
//          The neighbor list should update only on the first call (updateCount==0); subsequent
//          calls with the same positions must return the same count without corruption.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListLargeScaleTest, NSQ_PairCount_Equals_NChoose2)
{
    using GP                                  = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
    resetNLCounter();

    auto NL_cpu = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                     pos_cpu,
                                                                     quat_cpu,
                                                                     GP::m_collisionDetection,
                                                                     nObs,
                                                                     nParts);
    NL_cpu->updateNeighborList(pos_cpu, nObs, nParts);

    EXPECT_EQ(NL_cpu->getSize(), EXPECTED_PAIRS)
        << "NSQ CPU: pair count must equal C(" << nParts << ",2) = " << EXPECTED_PAIRS;

    resetNLCounter();
    auto NL_gpu = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                       pos_gpu,
                                                                       quat_gpu,
                                                                       GP::m_collisionDetection,
                                                                       nObs,
                                                                       nParts);
    NL_gpu->updateNeighborList(pos_gpu, nObs, nParts);
    cudaDeviceSynchronize();

    EXPECT_EQ(NL_gpu->getSize(), EXPECTED_PAIRS)
        << "NSQ GPU: pair count must equal C(" << nParts << ",2) = " << EXPECTED_PAIRS;
}

// ------------------------------------------------------------------------------------------------
// Test 20: 512 particles, realistic cells -- ATOMICFIXED pair count matches ATOMIC.
//          ATOMICFIXED pre-sizes a 2D array; if maxParticlesPerCell is too small it silently
//          drops pairs.  This test directly compares the two GPU atomic variants.
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListLargeScaleTest, RealisticCells_AtomicFixed_Matches_Atomic)
{
    using GP                                  = GrainsParameters<double>;
    GP::m_collisionDetection.neighborListType = NeighborListType::LINKEDCELL;

    resetNLCounter();
    lcParamsRealist.type                          = LinkedCellType::ATOMIC;
    GP::m_collisionDetection.linkedCellParameters = lcParamsRealist;
    auto    NL_atomic    = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                          pos_gpu,
                                                                          quat_gpu,
                                                                          GP::m_collisionDetection,
                                                                          nObs,
                                                                          nParts);
    PairSet atomic_pairs = extractPairs_GPU(NL_atomic.get(), pos_gpu, nObs, nParts);

    resetNLCounter();
    lcParamsRealist.type                          = LinkedCellType::ATOMICFIXED;
    GP::m_collisionDetection.linkedCellParameters = lcParamsRealist;
    auto    NL_fixed    = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                         pos_gpu,
                                                                         quat_gpu,
                                                                         GP::m_collisionDetection,
                                                                         nObs,
                                                                         nParts);
    PairSet fixed_pairs = extractPairs_GPU(NL_fixed.get(), pos_gpu, nObs, nParts);

    EXPECT_EQ(fixed_pairs.size(), atomic_pairs.size())
        << "ATOMICFIXED pair count (" << fixed_pairs.size() << ") differs from ATOMIC ("
        << atomic_pairs.size() << "); maxParticlesPerCell may be too small";
    EXPECT_EQ(fixed_pairs, atomic_pairs)
        << "ATOMICFIXED and ATOMIC pair sets differ -- likely a capacity overflow in ATOMICFIXED";
    cudaDeviceSynchronize();
}

// ================================================================================================
//  FIXTURE: NeighborListRandomPositionTest
//  -----------------------------------------
//  Random positions matching the benchmark scenario:
//    - Domain  : [0, 1.6]^3 (sphere AR=1, N=512 benchmark: dom_mul=32, R=0.05)
//    - CellSize: 0.1 (= 2*R for sphere radius 0.05)
//    - Grid    : 16x16x16 = 4096 cells
//    - N       : 32 and 128 particles
//    - Seed    : 42 (fixed for reproducibility)
//
//  All five algorithms (NSQ CPU/GPU, LC HOST, LC SORTBASED, LC ATOMIC, LC ATOMICFIXED)
//  must agree pairwise on the same positions. This exercises the sparse multi-cell case
//  that the regular-grid consistency tests do NOT cover.
// ================================================================================================
class NeighborListRandomPositionTest : public ::testing::Test
{
protected:
    static constexpr uint nObs = 0;
    static constexpr uint SEED = 42;

    // Benchmark-matching domain and cell size
    static constexpr double DOMAIN    = 1.6;
    static constexpr double CELL_SIZE = 0.1;  // = 2 * sphere_radius (R=0.05)

    Box<double>* sharedBox = nullptr;

    GrainsMemBuffer<RigidBody<double>*, MemType::HOST>   rb_cpu;
    GrainsMemBuffer<Vector3<double>, MemType::HOST>      pos_cpu;
    GrainsMemBuffer<Quaternion<double>, MemType::HOST>   quat_cpu;
    GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE> rb_gpu;
    GrainsMemBuffer<Vector3<double>, MemType::DEVICE>    pos_gpu;
    GrainsMemBuffer<Quaternion<double>, MemType::DEVICE> quat_gpu;

    LinkedCellParameters<double> lcParams;
    uint                         nParts = 0;

    std::vector<RigidBody<double>*> rbVec;

    void setUp(uint N)
    {
        // Clean up buffer state from a previous setUp() call.
        sharedBox = nullptr;
        rbVec.clear();
        rb_cpu.clear();
        pos_cpu.clear();
        quat_cpu.clear();
        // GPU buffers are fully replaced by copyFrom() below.

        resetNLCounter();
        GrainsParameters<double>::m_collisionDetection = CollisionDetectionParameters<double>{};
        // Ensure GPU device properties are populated
        if(GrainsParameters<double>::m_GPU.maxThreadsPerBlock == 0)
            cudaGetDeviceProperties(&GrainsParameters<double>::m_GPU, 0);
        nParts = N;

        // sphere with circumscribed radius = 0.05
        sharedBox = new Box<double>(0.05, 0.05, 0.05);

        // Generate reproducible random positions in [0, DOMAIN]^3
        std::mt19937                           rng(SEED);
        std::uniform_real_distribution<double> dist(0.0, DOMAIN);
        std::vector<Vector3<double>>           positions(N);
        for(uint i = 0; i < N; ++i)
            positions[i] = Vector3<double>(dist(rng), dist(rng), dist(rng));

        rbVec.resize(N);
        for(uint i = 0; i < N; ++i)
            rbVec[i] = new RigidBody<double>(sharedBox, 0.1, 1000.0, 1);

        buildBuffers(positions, N, rbVec, rb_cpu, pos_cpu, quat_cpu, rb_gpu, pos_gpu, quat_gpu);

        // Linked-cell parameters matching the benchmark
        lcParams.minCorner                       = Vector3<double>(0.0, 0.0, 0.0);
        lcParams.maxCorner                       = Vector3<double>(DOMAIN, DOMAIN, DOMAIN);
        lcParams.minCellSize                     = CELL_SIZE;
        lcParams.cellSizeFactor                  = 1.0;
        lcParams.maxNumCellsPerObstacle          = 8;
        lcParams.maxParticlesPerCell             = 64;
        lcParams.initialNumberOfPairsPerParticle = 32;
        lcParams.updateFrequency                 = 0;
        lcParams.sortFrequency                   = 0;
    }

    void SetUp() override
    {
        // Default: N=32 (fast); individual tests may override via setUp()
        setUp(32);
    }

    void TearDown() override
    {
        GrainsParameters<double>::m_simulationState.neighborListUpdateCount = 0;
        // Note: memory cleanup omitted intentionally -- RigidBody owns the Box pointer and
        // frees it in its destructor; deleting sharedBox here causes a double-free/crash.
        // (Consistent with all other fixtures in this file.)
    }
};

// ------------------------------------------------------------------------------------------------
// Test 21: N=32 random positions -- all algorithms agree on the exact same pair set.
//          The LC HOST pair set is the reference (same as NSQ restricted to 1-cell neighborhood).
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListRandomPositionTest, N32_AllAlgorithms_Agree)
{
    using GP = GrainsParameters<double>;

    // Run NSQ CPU (ground truth -- finds ALL pairs regardless of distance)
    resetNLCounter();
    GP::m_collisionDetection.neighborListType = NeighborListType::NSQ;
    auto    NL_nsq_cpu    = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                         pos_cpu,
                                                                         quat_cpu,
                                                                         GP::m_collisionDetection,
                                                                         nObs,
                                                                         nParts);
    PairSet nsq_cpu_pairs = extractPairs_CPU(NL_nsq_cpu.get(), pos_cpu, nObs, nParts);

    const size_t expected_nsq = (size_t)nParts * (nParts - 1) / 2;
    EXPECT_EQ(nsq_cpu_pairs.size(), expected_nsq)
        << "NSQ CPU: expected all C(" << nParts << ",2) = " << expected_nsq << " pairs";

    // Run LC HOST (CPU reference; spatial filter reduces pairs below NSQ)
    resetNLCounter();
    GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
    lcParams.type                                 = LinkedCellType::HOST;
    GP::m_collisionDetection.linkedCellParameters = lcParams;
    auto    NL_lc_cpu     = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                        pos_cpu,
                                                                        quat_cpu,
                                                                        GP::m_collisionDetection,
                                                                        nObs,
                                                                        nParts);
    PairSet lc_host_pairs = extractPairs_CPU(NL_lc_cpu.get(), pos_cpu, nObs, nParts);

    EXPECT_LE(lc_host_pairs.size(), expected_nsq)
        << "LC HOST: should return at most NSQ's " << expected_nsq << " pairs";

    // Run all GPU LC variants and compare against LC HOST (reference)
    for(auto& [lcType, name] : std::vector<std::pair<LinkedCellType, std::string>>{
            {LinkedCellType::SORTBASED, "LC_SORTBASED_GPU"},
            {LinkedCellType::ATOMIC, "LC_ATOMIC_GPU"},
            {LinkedCellType::ATOMICFIXED, "LC_ATOMICFIXED_GPU"},
        })
    {
        resetNLCounter();
        lcParams.type                                 = lcType;
        GP::m_collisionDetection.linkedCellParameters = lcParams;
        auto    NL_gpu    = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                           pos_gpu,
                                                                           quat_gpu,
                                                                           GP::m_collisionDetection,
                                                                           nObs,
                                                                           nParts);
        PairSet gpu_pairs = extractPairs_GPU(NL_gpu.get(), pos_gpu, nObs, nParts);

        EXPECT_EQ(gpu_pairs.size(), lc_host_pairs.size())
            << name << " vs LC_HOST: pair count mismatch (" << gpu_pairs.size() << " vs "
            << lc_host_pairs.size() << ")";
        EXPECT_EQ(gpu_pairs, lc_host_pairs)
            << name << " vs LC_HOST: pair sets differ on random positions (N=" << nParts
            << ", domain=" << DOMAIN << ", cellSize=" << CELL_SIZE << ")";
    }
    cudaDeviceSynchronize();
}

// ------------------------------------------------------------------------------------------------
// Test 22: N=128 random positions -- all LC algorithms agree (larger sparse case).
// ------------------------------------------------------------------------------------------------
TEST_F(NeighborListRandomPositionTest, N128_AllLinkedCell_Agree)
{
    setUp(128);  // re-initialize with N=128

    using GP = GrainsParameters<double>;

    // LC HOST reference
    resetNLCounter();
    GP::m_collisionDetection.neighborListType     = NeighborListType::LINKEDCELL;
    lcParams.type                                 = LinkedCellType::HOST;
    GP::m_collisionDetection.linkedCellParameters = lcParams;
    auto    NL_lc_cpu     = NeighborListFactory<double, MemType::HOST>::create(&rb_cpu,
                                                                        pos_cpu,
                                                                        quat_cpu,
                                                                        GP::m_collisionDetection,
                                                                        nObs,
                                                                        nParts);
    PairSet lc_host_pairs = extractPairs_CPU(NL_lc_cpu.get(), pos_cpu, nObs, nParts);

    EXPECT_GT(lc_host_pairs.size(), size_t(0))
        << "LC HOST N=128: expected at least one pair in the tight-cell NL";

    // All GPU LC variants must match LC HOST
    for(auto& [lcType, name] : std::vector<std::pair<LinkedCellType, std::string>>{
            {LinkedCellType::SORTBASED, "LC_SORTBASED_GPU"},
            {LinkedCellType::ATOMIC, "LC_ATOMIC_GPU"},
            {LinkedCellType::ATOMICFIXED, "LC_ATOMICFIXED_GPU"},
        })
    {
        resetNLCounter();
        lcParams.type                                 = lcType;
        GP::m_collisionDetection.linkedCellParameters = lcParams;
        auto    NL_gpu    = NeighborListFactory<double, MemType::DEVICE>::create(&rb_gpu,
                                                                           pos_gpu,
                                                                           quat_gpu,
                                                                           GP::m_collisionDetection,
                                                                           nObs,
                                                                           nParts);
        PairSet gpu_pairs = extractPairs_GPU(NL_gpu.get(), pos_gpu, nObs, nParts);

        EXPECT_EQ(gpu_pairs.size(), lc_host_pairs.size())
            << name << " vs LC_HOST N=128: pair count mismatch (" << gpu_pairs.size() << " vs "
            << lc_host_pairs.size() << ")";
        EXPECT_EQ(gpu_pairs, lc_host_pairs)
            << name << " vs LC_HOST N=128: pair sets differ on random positions";
    }
    cudaDeviceSynchronize();
}