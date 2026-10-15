#ifndef _GRAINSMEMBUFFER_HH_
#define _GRAINSMEMBUFFER_HH_

#include <type_traits>

#include "GrainsParameters.hh"
#include "GrainsUtils.hh"

enum class MemType
{
    HOST,
    DEVICE,
    MANAGED,
    PINNED,
    MAPPED,
    UNKNOWN
};

// =================================================================================================
/** @brief The class GrainsMemBuffer.

    This class provides a buffer that can be allocated in different memory spaces (host, device,
    managed, pinned) and provides methods for memory management and data transfer between these
    spaces.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
/** @name GrainsMemBuffer: External Methods */
//@{
/** @brief Fills a buffer to default.
    @param buffer the buffer to be initialized
    @param size size of the buffer
    @param value the default value to be set (default is T()) */
template <typename T>
__GLOBAL__ void fill_Kernel(T* buffer, const size_t size, const T& value = T())
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size)
        return;
    buffer[idx] = value;
}

// -------------------------------------------------------------------------------------------------
/** @brief fills a buffer with incremental values.
    @param buffer the buffer to be initialized
    @param size size of the buffer
    @param start the starting value (default is 0) */
template <typename T>
__GLOBAL__ void sequence_Kernel(T* buffer, const size_t size, const T& start = T(0))
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size)
        return;

    buffer[idx] = static_cast<T>(start + idx);
}

// -------------------------------------------------------------------------------------------------
/** @brief Kernel to free device pointers stored in a pointer array.
    @param ptrArray array of device pointers
    @param size number of pointers in the array */
template <typename T>
__GLOBAL__ void freePointedObjects_Kernel(T** ptrArray, const size_t size)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size)
        return;

    if(ptrArray[idx] != nullptr)
    {
        cudaFree(ptrArray[idx]);
        ptrArray[idx] = nullptr;
    }
}
//@}

// =================================================================================================
template <typename T, MemType M = MemType::HOST>
class GrainsMemBuffer
{
protected:
    /** @name Parameters */
    //@{
    /** \brief Pointer to the data */
    T* m_ptr = nullptr;
    /** \brief Pointer to the data on device (for zero-copy) */
    T* m_d_ptr = nullptr;
    /** \brief Number of elements in the buffer */
    size_t m_size = 0;
    /** \brief Capacity of the buffer */
    size_t m_capacity = 0;
    //@}

public:
    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    GrainsMemBuffer() = default;

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with the size */
    GrainsMemBuffer(size_t size)
    {
        initialize(size);
        fill();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with the size and default value.
        @param size size of the buffer
        @param value default value to fill the buffer */
    GrainsMemBuffer(size_t size, const T& value)
    {
        initialize(size);
        fill(value);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Destructor */
    ~GrainsMemBuffer()
    {
        free();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Copy constructor.
        @param other the other buffer to copy from */
    template <MemType srcM>
    GrainsMemBuffer(const GrainsMemBuffer<T, srcM>& other)
    {
        copyFrom(other);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Copy assignment operator.
        @param other the other buffer to copy from */
    template <MemType srcM>
    GrainsMemBuffer<T, M>& operator=(const GrainsMemBuffer<T, srcM>& other)
    {
        if(this != &other)
            copyFrom(other);
        return *this;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Move constructor */
    GrainsMemBuffer(GrainsMemBuffer<T, M>&& other) noexcept
    {
        moveFrom(other);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Move assignment operator */
    GrainsMemBuffer& operator=(GrainsMemBuffer<T, M>&& other) noexcept
    {
        if(this != &other)
        {
            // Free existing storage before taking ownership to avoid leaks
            free();
            moveFrom(other);
        }
        return (*this);
    }
    //@}

    /** @name Get methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Returns the pointer to the data */
    T* getData()
    {
        return m_ptr;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns the pointer to the data */
    const T* getData() const
    {
        return m_ptr;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns the pointer to the data on device (for zero-copy) */
    /** @note For PINNED memory, returns m_d_ptr if set, otherwise m_ptr.
              For MAPPED memory, always returns the mapped device alias m_d_ptr. */
    T* getDeviceData()
    {
        if constexpr(M == MemType::PINNED)
            return m_d_ptr ? m_d_ptr : m_ptr;
        else if constexpr(M == MemType::MAPPED)
            return m_d_ptr;
        else
            return m_ptr;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns the number of elements */
    size_t getSize() const
    {
        return m_size;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns the size in bytes */
    size_t getBytes() const
    {
        return m_size * sizeof(T);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns the capacity (maximum elements without reallocation) */
    size_t getCapacity() const
    {
        return m_capacity;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns the type of memory */
    MemType getMemType() const
    {
        if constexpr(M == MemType::HOST)
            return MemType::HOST;
        else if constexpr(M == MemType::DEVICE)
            return MemType::DEVICE;
        else if constexpr(M == MemType::MANAGED)
            return MemType::MANAGED;
        else if constexpr(M == MemType::PINNED)
            return MemType::PINNED;
        else if constexpr(M == MemType::MAPPED)
            return MemType::MAPPED;
        else
            return MemType::UNKNOWN;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns iterator to the beginning of the buffer */
    T* begin()
    {
        static_assert(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED
                          || M == MemType::MANAGED,
                      "begin() only available for HOST, PINNED, MAPPED, or MANAGED memory");
        return m_ptr;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns const iterator to the beginning of the buffer */
    const T* begin() const
    {
        static_assert(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED
                          || M == MemType::MANAGED,
                      "begin() only available for HOST, PINNED, MAPPED, or MANAGED memory");
        return m_ptr;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns iterator to the end of the buffer */
    T* end()
    {
        static_assert(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED
                          || M == MemType::MANAGED,
                      "end() only available for HOST, PINNED, MAPPED, or MANAGED memory");
        return m_ptr + m_size;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns const iterator to the end of the buffer */
    const T* end() const
    {
        static_assert(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED
                          || M == MemType::MANAGED,
                      "end() only available for HOST, PINNED, MAPPED, or MANAGED memory");
        return m_ptr + m_size;
    }
    //@}

    /** @name Set methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Sets the size of the buffer without reallocation.
        @param new_size new size of the buffer (must be <= capacity) */
    void setSize(size_t new_size)
    {
        GAssert(new_size <= m_capacity, "Size must be <= capacity in setSize()");
        m_size = new_size;
    }
    //@}

    /** @name Methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Reserves memory for the buffer.
        @param new_capacity new capacity of the buffer
        @param stream CUDA stream for asynchronous operations (default: 0) */
    void reserve(size_t new_capacity, cudaStream_t stream = 0)
    {
        GAssert(new_capacity > 0, "Capacity must be positive in reserve()");

        if(new_capacity <= m_capacity)
            return;

        // Allocate a new buffer with the requested capacity
        GrainsMemBuffer<T, M> new_buf;
        new_buf.initialize(new_capacity, new_capacity);

        // Copy existing data (only logical size worth of bytes)
        if(m_ptr && m_size)
        {
            if constexpr(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED)
            {
                std::memcpy(new_buf.getData(), m_ptr, m_size * sizeof(T));
            }
            else if constexpr(M == MemType::DEVICE || M == MemType::MANAGED)
            {
                // Device-to-device copy for device/managed memory using stream
                cudaErrCheck(cudaMemcpyAsync(new_buf.getData(),
                                             m_ptr,
                                             m_size * sizeof(T),
                                             cudaMemcpyDeviceToDevice,
                                             stream));
                // Synchronize before freeing old buffer to ensure copy completes
                cudaStreamSynchronize(stream);
            }
        }

        // Preserve logical size; capacity is already new_capacity
        new_buf.m_size = m_size;

        // Replace current storage; move assignment frees old storage now
        *this = std::move(new_buf);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Initialize/reinitialize buffer with specific size and capacity.
        This is useful for buffer initialization or complete reallocation with new size/capacity.
        @param new_size desired logical size
        @param new_capacity desired capacity (if 0, uses new_size) */
    void initialize(size_t new_size, size_t new_capacity = 0)
    {
        if(new_capacity == 0)
            new_capacity = new_size;

        GAssert(new_capacity >= new_size, "Capacity must be >= size in initialize()");

        // Free existing memory
        free();

        // Allocate with desired capacity
        m_capacity = new_capacity;
        m_size     = new_size;

        if constexpr(M == MemType::HOST)
        {
            m_ptr = static_cast<T*>(std::malloc(sizeof(T) * new_capacity));
            if(!m_ptr)
                throw std::bad_alloc();
        }
        else if constexpr(M == MemType::DEVICE)
        {
            cudaErrCheck(cudaMalloc(&m_ptr, sizeof(T) * new_capacity));
        }
        else if constexpr(M == MemType::MANAGED)
        {
            cudaErrCheck(cudaMallocManaged(&m_ptr, sizeof(T) * new_capacity));
        }
        else if constexpr(M == MemType::PINNED)
        {
            cudaErrCheck(cudaMallocHost(&m_ptr, sizeof(T) * new_capacity));
        }
        else if constexpr(M == MemType::MAPPED)
        {
            cudaErrCheck(cudaHostAlloc(&m_ptr, sizeof(T) * new_capacity, cudaHostAllocMapped));
            cudaErrCheck(cudaHostGetDevicePointer(&m_d_ptr, m_ptr, 0));
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Resizes the buffer (changes logical size, may grow capacity).
        @param new_size new size of the buffer */
    void resize(size_t new_size)
    {
        // If new size fits within current capacity, just change size
        if(new_size <= m_capacity)
        {
            m_size = new_size;
            return;
        }

        // Need to grow capacity - preserve existing data
        reserve(new_size);
        m_size = new_size;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Clears the buffer (sets size to 0, keeps capacity). Useful when you want to "empty"
        the buffer but keep memory allocated. */
    void clear()
    {
        m_size = 0;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Resets the buffer completely (frees memory, size=0, capacity=0). Useful when you want
        to completely reinitialize the buffer. */
    void reset()
    {
        free();
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Pushes a new element to the back of the buffer.
        @param value value to push */
    void push_back(const T& value)
    {
        if constexpr(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED
                     || M == MemType::MANAGED)
        {
            // Grow capacity if needed (like std::vector)
            if(m_size >= m_capacity)
            {
                size_t new_capacity = m_capacity == 0 ? 1 : m_capacity * 2;
                reserve(new_capacity);
            }
            m_ptr[m_size++] = value;
        }
        else
            std::cerr << "GrainsMemBuffer::push_back() only allowed on host, "
                         "pinned, or managed memory\n";
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Pushes a bulk of elements to the back of the buffer.
        @param values pointer to the array of values
        @param count number of elements to push */
    void push_bulk(const T* values, size_t count)
    {
        if constexpr(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED
                     || M == MemType::MANAGED)
        {
            // Grow capacity if needed
            if(m_size + count > m_capacity)
            {
                size_t new_capacity
                    = std::max(m_capacity == 0 ? 1 : m_capacity * 2, m_size + count);
                reserve(new_capacity);
            }
            T* dst = m_ptr + m_size;
            std::memcpy(dst, values, count * sizeof(T));
            m_size += count;
        }
        else
            std::cerr << "GrainsMemBuffer::push_bulk() only allowed on host, "
                         "pinned, or managed memory\n";
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Shrinks the buffer to fit the current size
        @param stream CUDA stream for asynchronous operations (default: 0) */
    void shrink_to_fit(cudaStream_t stream = 0)
    {
        if(m_size == m_capacity)
            return;

        GrainsMemBuffer<T, M> new_buf{};
        new_buf.initialize(m_size, m_size);

        if constexpr(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED)
            std::memcpy(new_buf.m_ptr, m_ptr, m_size * sizeof(T));
        else if constexpr(M == MemType::DEVICE || M == MemType::MANAGED)
        {
            cudaMemcpyKind kind = cudaMemcpyDeviceToDevice;
            cudaErrCheck(
                cudaMemcpyAsync(new_buf.getData(), m_ptr, m_size * sizeof(T), kind, stream));
        }
        else
        {
            std::cerr << "Unsupported memory type for shrink_to_fit()\n";
            return;
        }

        *this = std::move(new_buf);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns the kind of memory transfer for the given src and dst memory types. */
    template <MemType Src, MemType Dst>
    constexpr cudaMemcpyKind getMemcpyKind() const
    {
        if constexpr(Src == MemType::HOST && Dst == MemType::HOST)
            return cudaMemcpyHostToHost;
        else if constexpr(Src == MemType::HOST && Dst == MemType::DEVICE)
            return cudaMemcpyHostToDevice;
        else if constexpr(Src == MemType::DEVICE && Dst == MemType::HOST)
            return cudaMemcpyDeviceToHost;
        else if constexpr(Src == MemType::DEVICE && Dst == MemType::DEVICE)
            return cudaMemcpyDeviceToDevice;
        else if constexpr(Src == MemType::HOST && Dst == MemType::MANAGED)
            return cudaMemcpyHostToDevice;
        else if constexpr(Src == MemType::MANAGED && Dst == MemType::HOST)
            return cudaMemcpyDeviceToHost;
        else if constexpr(Src == MemType::DEVICE && Dst == MemType::MANAGED)
            return cudaMemcpyDeviceToDevice;
        else if constexpr(Src == MemType::MANAGED && Dst == MemType::DEVICE)
            return cudaMemcpyDeviceToDevice;
        else if constexpr(Src == MemType::MANAGED && Dst == MemType::MANAGED)
            return cudaMemcpyDeviceToDevice;
        else if constexpr(Src == MemType::PINNED && Dst == MemType::HOST)
            return cudaMemcpyHostToHost;
        else if constexpr(Src == MemType::HOST && Dst == MemType::PINNED)
            return cudaMemcpyHostToHost;
        else if constexpr(Src == MemType::PINNED && Dst == MemType::DEVICE)
            return cudaMemcpyHostToDevice;
        else if constexpr(Src == MemType::DEVICE && Dst == MemType::PINNED)
            return cudaMemcpyDeviceToHost;
        else if constexpr(Src == MemType::PINNED && Dst == MemType::PINNED)
            return cudaMemcpyHostToHost;
        else if constexpr(Src == MemType::MAPPED && Dst == MemType::HOST)
            return cudaMemcpyHostToHost;
        else if constexpr(Src == MemType::HOST && Dst == MemType::MAPPED)
            return cudaMemcpyHostToHost;
        else if constexpr(Src == MemType::MAPPED && Dst == MemType::DEVICE)
            return cudaMemcpyHostToDevice;
        else if constexpr(Src == MemType::DEVICE && Dst == MemType::MAPPED)
            return cudaMemcpyDeviceToHost;
        else if constexpr(Src == MemType::MAPPED && Dst == MemType::MAPPED)
            return cudaMemcpyHostToHost;
        else if constexpr(Src == MemType::PINNED && Dst == MemType::MAPPED)
            return cudaMemcpyHostToHost;
        else if constexpr(Src == MemType::MAPPED && Dst == MemType::PINNED)
            return cudaMemcpyHostToHost;
        else
            return cudaMemcpyDefault;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Copy to another buffer (host/device aware).
        @param dest destination buffer
        @param stream CUDA stream for asynchronous operations (default: 0) */
    template <MemType destM>
    void copyTo(GrainsMemBuffer<T, destM>& dest, cudaStream_t stream = 0) const
    {
        if(m_size == 0 || !m_ptr)
            return;

        GAssert(dest.getSize() >= m_size, "Destination buffer too small for copy");

        if constexpr(M == MemType::HOST && destM == MemType::HOST)
            std::memcpy(dest.getData(), m_ptr, getBytes());
        else
        {
            cudaMemcpyKind kind    = getMemcpyKind<M, destM>();
            const T*       src_ptr = m_ptr;
            T*             dst_ptr = dest.getData();
            cudaErrCheck(cudaMemcpyAsync(dst_ptr, src_ptr, getBytes(), kind, stream));
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Copy from another buffer (host/device aware).
        @param src source buffer
        @param stream CUDA stream for asynchronous operations (default: 0) */
    template <MemType srcM>
    void copyFrom(const GrainsMemBuffer<T, srcM>& src, cudaStream_t stream = 0)
    {
        if constexpr(M == srcM)
        {
            if(this == &src)
                return;
        }

        // Resize this buffer if needed
        if(m_capacity < src.getSize())
            resize(src.getSize());
        m_size = src.getSize();

        if constexpr(M == MemType::HOST && srcM == MemType::HOST)
        {
            std::memcpy(m_ptr, src.m_ptr, src.getBytes());
        }
        else
        {
            cudaMemcpyKind kind    = getMemcpyKind<srcM, M>();
            const T*       src_ptr = src.getData();
            T*             dst_ptr = m_ptr;
            cudaErrCheck(cudaMemcpyAsync(dst_ptr, src_ptr, src.getBytes(), kind, stream));
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief at method to access elements with bounds checking.
        @param index index of the element to access */
    const T& at(size_t index) const
    {
        static_assert(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED
                          || M == MemType::MANAGED,
                      "at() only available for HOST, PINNED, MAPPED, or MANAGED memory");
        GAssert(index < m_size, "Index", index, "out of bounds for size", m_size);
        return m_ptr[index];
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Frees the allocated memory */
    void free()
    {
        if(m_ptr)
        {
            if constexpr(M == MemType::HOST)
            {
                std::free(m_ptr);
            }
            else if constexpr(M == MemType::DEVICE)
            {
                cudaErrCheck(cudaFree(m_ptr));
            }
            else if constexpr(M == MemType::PINNED)
            {
                cudaErrCheck(cudaFreeHost(m_ptr));
            }
            else if constexpr(M == MemType::MAPPED)
            {
                cudaErrCheck(cudaFreeHost(m_ptr));
                // m_d_ptr is a device alias into the same physical allocation
            }
            else if constexpr(M == MemType::MANAGED)
            {
                cudaErrCheck(cudaFree(m_ptr));
            }
            else
            {
                std::cerr << "Unknown memory type for free()\n";
            }
        }

        m_ptr      = nullptr;
        m_d_ptr    = nullptr;
        m_size     = 0;
        m_capacity = 0;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Frees objects pointed to by pointers in this buffer.
        This is only valid for buffers that store pointers (e.g., GrainsMemBuffer<T*, M>).
        For HOST memory, calls delete on each pointer.
        For DEVICE memory, copies pointers to host and calls cudaFree on each.
        @param stream CUDA stream for asynchronous operations (default: 0) */
    void freePointedObjects(cudaStream_t stream = 0)
    {
        // Only allow this for pointer types
        static_assert(std::is_pointer_v<T>, "freePointedObjects only valid for pointer types");
        static_assert(M != MemType::MANAGED, "freePointedObjects not supported for MANAGED memory");

        if(m_size == 0 || !m_ptr)
            return;

        if constexpr(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED)
        {
            // Host memory: directly iterate and delete
            for(size_t i = 0; i < m_size; ++i)
            {
                if(m_ptr[i] != nullptr)
                {
                    delete m_ptr[i];
                    m_ptr[i] = nullptr;
                }
            }
        }
        else if constexpr(M == MemType::DEVICE)
        {
            // Device memory: copy pointers to host first, then free each
            std::vector<T> hostPtrs(m_size);
            cudaErrCheck(
                cudaMemcpy(hostPtrs.data(), m_ptr, m_size * sizeof(T), cudaMemcpyDeviceToHost));
            std::cout << "Freeing " << m_size << " pointed objects from device memory...\n";
            for(size_t i = 0; i < m_size; ++i)
            {
                if(hostPtrs[i] != nullptr)
                    cudaErrCheck(cudaFree(hostPtrs[i]));
            }
            std::cout << "Freed " << m_size << " pointed objects from device memory.\n";
            fill();
            cudaErrCheck(cudaDeviceSynchronize());
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Fills the buffer with default values
        @param stream CUDA stream for asynchronous operations (default: 0) */
    void fill(cudaStream_t stream = 0)
    {
        if constexpr(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED)
        {
            for(size_t i = 0; i < m_size; ++i)
            {
                m_ptr[i] = T();
            }
        }
        else if constexpr(M == MemType::DEVICE || M == MemType::MANAGED)
        {
            // Use optimized cudaMemsetAsync for zero initialization
            cudaMemsetAsync(m_ptr, 0, m_size * sizeof(T), stream);
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Fills the buffer with a user-provided value.
        @param value value to initialize with
        @param stream CUDA stream for asynchronous operations (default: 0) */
    void fill(const T& value, cudaStream_t stream = 0)
    {
        if constexpr(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED)
        {
            for(size_t i = 0; i < m_size; ++i)
                m_ptr[i] = value;
        }
        else if constexpr(M == MemType::DEVICE || M == MemType::MANAGED)
        {
            static_assert(std::is_fundamental<T>::value,
                          "T must be a primitive type for device init");

            // Use cudaMemset for common patterns: 0 and UINT_MAX
            if(value == T(0))
            {
                cudaMemsetAsync(m_ptr, 0, m_size * sizeof(T), stream);
            }
            else if constexpr(std::is_unsigned<T>::value)
            {
                if(value == static_cast<T>(~T(0)))  // UINT_MAX for this type
                {
                    cudaMemsetAsync(m_ptr, 0xFF, m_size * sizeof(T), stream);
                }
                else
                {
                    fill_Kernel<<<(m_size + 255) / 256, 256, 0, stream>>>(m_ptr, m_size, value);
                }
            }
            else
            {
                fill_Kernel<<<(m_size + 255) / 256, 256, 0, stream>>>(m_ptr, m_size, value);
            }
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Fills the buffer incrementally with values starting from a given value.
        @param start the starting value (default is 0)
        @param stream CUDA stream for asynchronous operations (default: 0) */
    void sequence(const T& start = T(0), cudaStream_t stream = 0)
    {
        if constexpr(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED)
        {
            for(size_t i = 0; i < m_size; ++i)
                m_ptr[i] = static_cast<T>(start + i);
        }
        else if constexpr(M == MemType::DEVICE || M == MemType::MANAGED)
        {
            static_assert(std::is_fundamental<T>::value,
                          "T must be a primitive type for device init");
            sequence_Kernel<<<(m_size + 255) / 256, 256, 0, stream>>>(m_ptr, m_size, start);
        }
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Prints the buffer contents (host/pinned/device/managed)
        @param label optional label to print before the buffer contents
        @param stream CUDA stream for asynchronous operations (default: 0) */
    void print(const std::string& label = "", cudaStream_t stream = 0) const
    {
        if constexpr(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED)
        {
            if(!label.empty())
                std::cout << label << ": " << "\n";
            for(size_t i = 0; i < m_size; ++i)
                std::cout << "[" << i << "]. " << m_ptr[i] << "\n";
            std::cout << std::endl;
        }
        else if constexpr(M == MemType::DEVICE || M == MemType::MANAGED)
        {
            // Copy to host and print
            std::vector<T> hostBuf(m_size);
            cudaErrCheck(
                cudaMemcpyAsync(hostBuf.data(), m_ptr, getBytes(), cudaMemcpyDeviceToHost, stream));
            cudaStreamSynchronize(stream);  // Must synchronize before printing
            if(!label.empty())
                std::cout << label << ": " << "\n";
            for(size_t i = 0; i < m_size; ++i)
                std::cout << "[" << i << "]. " << hostBuf[i] << "\n";
            std::cout << std::endl;
        }
        else
        {
            std::cerr << "print() unsupported for this memory type" << std::endl;
        }
    }
    //@}

    /** @name Operators */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Index operator
    @param i index of the element */
    T& operator[](size_t i)
    {
        static_assert(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED
                          || M == MemType::MANAGED,
                      "operator[] only available for HOST, PINNED, MAPPED, or MANAGED memory");
        GAssert(i < m_size, "Index", i, "out of bounds for size", m_size);
        return m_ptr[i];
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Index operator
        @param i index of the element */
    const T& operator[](size_t i) const
    {
        static_assert(M == MemType::HOST || M == MemType::PINNED || M == MemType::MAPPED
                          || M == MemType::MANAGED,
                      "operator[] only available for HOST, PINNED, MAPPED, or MANAGED memory");
        GAssert(i < m_size, "Index", i, "out of bounds for size", m_size);
        return m_ptr[i];
    }

private:
    // ---------------------------------------------------------------------------------------------
    /** @brief Move from another buffer
        @param other source buffer */
    void moveFrom(GrainsMemBuffer<T, M>& other)
    {
        m_ptr      = other.m_ptr;
        m_d_ptr    = other.m_d_ptr;
        m_size     = other.m_size;
        m_capacity = other.m_capacity;

        other.m_ptr      = nullptr;
        other.m_d_ptr    = nullptr;
        other.m_size     = 0;
        other.m_capacity = 0;
    }
};

#endif
