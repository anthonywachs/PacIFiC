#ifndef _GRAINSUTILS_HH_
#define _GRAINSUTILS_HH_

#include "Basic.hh"
#include "Vector3.hh"
#include <unistd.h>

// =================================================================================================
/** @brief Miscellaneous functionalities (mostly low-level) for Grains.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
/** @name Miscellaneous functions and utilities for Grains */
//@{
// Macro for outputting CUDA errors
#define cudaErrCheck(ans) cudaAssert((ans), __FILE__, __LINE__);

/** @brief Returns CUDA error.
    @param code the error code
    @param file the file name
    @param line the line number
    @param abort whether to abort the program */
__HOST__ static INLINE void
    cudaAssert(cudaError_t code, const char* file, int line, bool abort = false)
{
    if(code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if(abort)
            exit(code);
    }
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the available memory on the host in bytes */
__HOST__ static INLINE size_t getAvailableHostMemory()
{
    long pages     = sysconf(_SC_AVPHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the available memory on the device in bytes */
__HOST__ static INLINE size_t getAvailableDeviceMemory()
{
    size_t free_byte;
    size_t total_byte;
    cudaErrCheck(cudaMemGetInfo(&free_byte, &total_byte));
    return free_byte;
}

// -------------------------------------------------------------------------------------------------
/** @brief computes the optimal number of threads and blocks for a given number of elements and an
    architecture
    @param numElements the number of elements
    @param prop the device properties
    @param numBlocks the number of blocks (output)
    @param numThreads the number of threads per block (output) */
__HOST__ static INLINE void computeOptimalThreadsAndBlocks(const uint            numElements,
                                                           const cudaDeviceProp& prop,
                                                           uint&                 numBlocks,
                                                           uint&                 numThreads)
{
    if(numElements == 0)
    {
        fprintf(stderr, "[FATAL] computeOptimalThreadsAndBlocks called with numElements=0\n");
        abort();
    }
    constexpr uint maxThreads = 256;  // Avoid 1024 unless necessary
    constexpr uint minThreads = 32;
    const uint     minBlocks  = 2 * prop.multiProcessorCount;

    // Start with 128 threads and compute how many blocks we need
    numThreads = 128;
    numBlocks  = (numElements + numThreads - 1) / numThreads;

    // If we're not filling the SMs enough, reduce threads but ensure full coverage
    if(numBlocks < minBlocks && numThreads > minThreads)
    {
        numThreads = minThreads;
        numBlocks  = max(minBlocks, (numElements + numThreads - 1) / numThreads);
    }
    // If we're oversubscribing the SMs too much, increase thread count
    if(numBlocks > 4 * prop.multiProcessorCount && numThreads < maxThreads)
    {
        numThreads = maxThreads;
        numBlocks  = (numElements + numThreads - 1) / numThreads;
    }
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a uint2 object in a string
    @param os the output stream
    @param v the uint2 object */
__HOST__ static INLINE std::ostream& operator<<(std::ostream& os, const uint2& v)
{
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a real number with a prescribed number of digits in a string
    @param figure the float number
    @param size number of digits */
template <typename T>
__HOST__ static INLINE std::string realToString(const T& figure, const int size)
{
    std::ostringstream oss;
    oss.width(size);
    oss << std::left << figure;
    return (oss.str());
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a float number with a prescribed format and a prescribed number of digits after
    the decimal point in a string.
    @param format the format
    @param digits number of digits after the decimal point
    @param number the float number */
template <typename T>
__HOST__ static INLINE std::string
    realToString(std::ios_base::fmtflags format, const int digits, const T& number)
{
    std::ostringstream oss;
    if(number != T(0))
    {
        oss.setf(format, std::ios::floatfield);
        oss.precision(digits);
    }
    oss << number;
    return (oss.str());
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a vector3 object in a string
    @param vec the vector3 object */
template <typename T>
__HOST__ static INLINE std::string Vector3ToString(const Vector3<T>& vec)
{
    std::ostringstream oss;
    oss << vec;
    return (oss.str());
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a message to stdout
    @param args the output messages */
template <typename... Args>
__HOST__ static constexpr INLINE void Gout(const Args&... args)
{
    ((std::cout << args << " "), ...);
    std::cout << std::endl;
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a message to stdout with Indent (WI)
    @param numShift the number of shift characters at the beginning
    @param args the output messages */
template <typename... Args>
__HOST__ INLINE void GoutWI(const int numShift, const Args&... args)
{
    // indent
    auto shift = [](int n) { return std::string(n, ' '); };

    std::cout << shift(numShift);
    ((std::cout << args << " "), ...);
    std::cout << std::endl;
}

// -------------------------------------------------------------------------------------------------
/** @brief Writes a message to stdout with Indent (WI)
    @param numShift the number of shift characters at the beginning
    @param args the output messages */
// Helper functions for device printf with different types
__DEVICE__ INLINE void print_device_arg(const char* arg)
{
    printf("%s ", arg);
}
__DEVICE__ INLINE void print_device_arg(char* arg)
{
    printf("%s ", arg);
}
__DEVICE__ INLINE void print_device_arg(size_t arg)
{
    printf("%zu ", arg);
}
__DEVICE__ INLINE void print_device_arg(int arg)
{
    printf("%d ", arg);
}
__DEVICE__ INLINE void print_device_arg(uint arg)
{
    printf("%u ", arg);
}
__DEVICE__ INLINE void print_device_arg(long arg)
{
    printf("%ld ", arg);
}
__DEVICE__ INLINE void print_device_arg(float arg)
{
    printf("%f ", arg);
}
__DEVICE__ INLINE void print_device_arg(double arg)
{
    printf("%f ", arg);
}
__HOST__ INLINE void print_device_arg(const std::string& arg)
{
    printf("%s ", arg.c_str());
}

template <typename... Args>
__HOSTDEVICE__ INLINE void GAbort(const Args&... args)
{
#ifdef __CUDA_ARCH__
    printf("[DEVICE] ");
    (print_device_arg(args), ...);
    printf("\n");
    __trap();  // aborts the kernel
#else
    std::cerr << "[HOST] ";
    ((std::cerr << args << " "), ...);
    std::cerr << std::endl;
    std::abort();
#endif
}

// -------------------------------------------------------------------------------------------------
/** @brief Assert function that aborts the program if the condition is false
    @param condition the condition to check
    @param args the message(s) to display if the assertion fails */
template <typename... Args>
__HOSTDEVICE__ INLINE void GAssert(bool condition, const Args&... args)
{
    if(!condition)
        GAbort("GAssert failed:", args...);
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes a triangular hash for a pair of unsigned integers. Maps two IDs (x, y) to a
    unique single ID using upper-triangular indexing:
    hash = l * (l + 1) / 2 + s, where l = max(x,y) and s = min(x,y).
    This produces a symmetric hash (order-independent) for pairs.
    @param x first ID
    @param y second ID */
__HOSTDEVICE__ static INLINE uint triangularHash(uint x, uint y)
{
    uint s = min(x, y);  // smaller one
    uint l = max(x, y);  // larger one
    return (l * (l + 1) / 2 + s);
}

// -------------------------------------------------------------------------------------------------
/** @brief Computes a hash value for a pair of IDs using prime number multiplication. This provides
    good distribution of hash values to minimize collisions.
    @param k pair of IDs (i,j) where i < j */
__HOSTDEVICE__ static INLINE uint primeHash(uint2 k)
{
    // Simple, fast hash using prime numbers for good distribution
    if(k.x > k.y)
        swap(k.x, k.y);
    uint h = k.x * 73856093u ^ k.y * 19349663u;
    return h;
}
//@}

#endif
