#include <cuda_runtime.h>
#include <gtest/gtest.h>
#include <set>

#include "ContactTable.hh"

__GLOBAL__ void testInsertionKernel(
    ContactMemoryView<float> view, uint2* pairs, uint* indices, bool* results, int numPairs)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= numPairs)
        return;
    results[idx] = view.findOrInsert(pairs[idx], indices[idx]);
}

__GLOBAL__ void testFindKernel(
    ContactMemoryView<float> view, uint2* pairs, uint* indices, bool* results, int numPairs)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= numPairs)
        return;
    results[idx] = view.find(pairs[idx], indices[idx]);
}

class ContactTableDeviceTest : public ::testing::Test
{
protected:
    ContactHashTable<float, MemType::DEVICE> table;
    ContactMemoryView<float>                 view;

    void SetUp() override
    {
        int deviceCount;
        cudaGetDeviceCount(&deviceCount);
        if(deviceCount == 0)
            GTEST_SKIP() << "No CUDA devices available";

        // Default capacity used by most tests
        table.allocate(200, 200);
        view = table.getView();
    }

    void TearDown() override
    {
        table.deallocate();
    }
};

TEST_F(ContactTableDeviceTest, DeviceBasicInsertion)
{
    const int numPairs = 100;

    uint2* h_pairs   = new uint2[numPairs];
    uint*  h_indices = new uint[numPairs];
    bool*  h_results = new bool[numPairs];
    for(int i = 0; i < numPairs; ++i)
        h_pairs[i] = make_uint2(i, i + 100);

    uint2* d_pairs;
    uint*  d_indices;
    bool*  d_results;
    cudaMalloc(&d_pairs, numPairs * sizeof(uint2));
    cudaMalloc(&d_indices, numPairs * sizeof(uint));
    cudaMalloc(&d_results, numPairs * sizeof(bool));
    cudaMemcpy(d_pairs, h_pairs, numPairs * sizeof(uint2), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (numPairs + blockSize - 1) / blockSize;
    testInsertionKernel<<<numBlocks, blockSize>>>(view, d_pairs, d_indices, d_results, numPairs);
    cudaDeviceSynchronize();

    cudaMemcpy(h_indices, d_indices, numPairs * sizeof(uint), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_results, d_results, numPairs * sizeof(bool), cudaMemcpyDeviceToHost);

    for(int i = 0; i < numPairs; ++i)
    {
        EXPECT_TRUE(h_results[i]);
        EXPECT_LT(h_indices[i], (uint)numPairs);
    }

    std::set<uint> uniqueIndices(h_indices, h_indices + numPairs);
    EXPECT_EQ(uniqueIndices.size(), numPairs);

    delete[] h_pairs;
    delete[] h_indices;
    delete[] h_results;
    cudaFree(d_pairs);
    cudaFree(d_indices);
    cudaFree(d_results);
}

TEST_F(ContactTableDeviceTest, DeviceFindAfterInsertion)
{
    const int numPairs = 50;

    uint2* h_pairs          = new uint2[numPairs];
    uint*  h_indices_insert = new uint[numPairs];
    uint*  h_indices_find   = new uint[numPairs];
    bool*  h_results_insert = new bool[numPairs];
    bool*  h_results_find   = new bool[numPairs];
    for(int i = 0; i < numPairs; ++i)
        h_pairs[i] = make_uint2(i * 2, i * 2 + 1);

    uint2* d_pairs;
    uint*  d_indices;
    bool*  d_results;
    cudaMalloc(&d_pairs, numPairs * sizeof(uint2));
    cudaMalloc(&d_indices, numPairs * sizeof(uint));
    cudaMalloc(&d_results, numPairs * sizeof(bool));
    cudaMemcpy(d_pairs, h_pairs, numPairs * sizeof(uint2), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (numPairs + blockSize - 1) / blockSize;
    testInsertionKernel<<<numBlocks, blockSize>>>(view, d_pairs, d_indices, d_results, numPairs);
    cudaDeviceSynchronize();
    cudaMemcpy(h_indices_insert, d_indices, numPairs * sizeof(uint), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_results_insert, d_results, numPairs * sizeof(bool), cudaMemcpyDeviceToHost);

    testFindKernel<<<numBlocks, blockSize>>>(view, d_pairs, d_indices, d_results, numPairs);
    cudaDeviceSynchronize();
    cudaMemcpy(h_indices_find, d_indices, numPairs * sizeof(uint), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_results_find, d_results, numPairs * sizeof(bool), cudaMemcpyDeviceToHost);

    for(int i = 0; i < numPairs; ++i)
    {
        EXPECT_TRUE(h_results_insert[i]);
        EXPECT_TRUE(h_results_find[i]);
        EXPECT_EQ(h_indices_insert[i], h_indices_find[i]);
    }

    delete[] h_pairs;
    delete[] h_indices_insert;
    delete[] h_indices_find;
    delete[] h_results_insert;
    delete[] h_results_find;
    cudaFree(d_pairs);
    cudaFree(d_indices);
    cudaFree(d_results);
}

TEST_F(ContactTableDeviceTest, ConcurrentInsertions)
{
    const int                                numPairs = 1000;
    const int                                capacity = 2000;
    ContactHashTable<float, MemType::DEVICE> table(capacity, capacity);
    ContactMemoryView<float>                 view = table.getView();

    uint2* h_pairs   = new uint2[numPairs];
    uint*  h_indices = new uint[numPairs];
    bool*  h_results = new bool[numPairs];
    for(int i = 0; i < numPairs; ++i)
        h_pairs[i] = make_uint2(i, i + 10000);

    uint2* d_pairs;
    uint*  d_indices;
    bool*  d_results;
    cudaMalloc(&d_pairs, numPairs * sizeof(uint2));
    cudaMalloc(&d_indices, numPairs * sizeof(uint));
    cudaMalloc(&d_results, numPairs * sizeof(bool));
    cudaMemcpy(d_pairs, h_pairs, numPairs * sizeof(uint2), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (numPairs + blockSize - 1) / blockSize;
    testInsertionKernel<<<numBlocks, blockSize>>>(view, d_pairs, d_indices, d_results, numPairs);
    cudaError_t err = cudaDeviceSynchronize();
    ASSERT_EQ(err, cudaSuccess) << "CUDA error: " << cudaGetErrorString(err);

    cudaMemcpy(h_indices, d_indices, numPairs * sizeof(uint), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_results, d_results, numPairs * sizeof(bool), cudaMemcpyDeviceToHost);

    int successCount = 0;
    for(int i = 0; i < numPairs; ++i)
        if(h_results[i])
            successCount++;
    EXPECT_GT(successCount, (int)(numPairs * 0.95));

    std::set<uint> uniqueIndices;
    for(int i = 0; i < numPairs; ++i)
        if(h_results[i])
            uniqueIndices.insert(h_indices[i]);
    EXPECT_EQ((int)uniqueIndices.size(), successCount);

    delete[] h_pairs;
    delete[] h_indices;
    delete[] h_results;
    cudaFree(d_pairs);
    cudaFree(d_indices);
    cudaFree(d_results);
}
