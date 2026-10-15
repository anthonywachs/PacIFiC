#ifndef _CONTACTTABLE_HH_
#define _CONTACTTABLE_HH_

#include "Basic.hh"
#include "GrainsMemBuffer.hh"
#include "GrainsUtils.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief Hash table entry for contact tracking.

    Each entry stores a key (pair of component IDs), an index pointing to the ContactForce in a
    flat array, and a flag for hash table management.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
struct ContactEntry
{
    /** @name Parameters */
    //@{
    /** \brief Contact pair key (i,j) where i < j */
    uint2 m_key;
    /** \brief Index into the ContactForce array */
    uint m_index;
    /** \brief Validity and activity flag: 0=empty, 1=stale (not touched this step), 2=active */
    uint m_valid;
    //@}

    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    __HOSTDEVICE__
    ContactEntry()
        : m_key(make_uint2(0, 0))
        , m_index(0xFFFFFFFF)  // Invalid index marker
        , m_valid(0)
    {
    }
    //@}
};

// =================================================================================================
/** @brief Contact history data for memory-enabled force models.

    Stores the cumulative tangential and rolling displacements for a contact pair.
    This data is accessed via the index returned by the hash table.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
struct ContactHistory
{
    /** @name Parameters */
    //@{
    /** \brief Cumulative tangential displacement (kt * delta) */
    Vector3<T> m_tangentialDisplacement;
    /** \brief Cumulative rolling friction spring-torque */
    Vector3<T> m_rollingDisplacement;
    /** \brief Previous contact normal (for plane rotation between timesteps) */
    Vector3<T> m_previousNormal;
    //@}

    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    __HOSTDEVICE__
    ContactHistory()
        : m_tangentialDisplacement(Vector3<T>(T(0), T(0), T(0)))
        , m_rollingDisplacement(Vector3<T>(T(0), T(0), T(0)))
        , m_previousNormal(Vector3<T>(T(0), T(0), T(0)))
    {
    }
    //@}
};

// =================================================================================================
/** @brief Complete contact memory view for passing to force models.

    This unified structure contains everything needed for contact history tracking: hash table
    pointers for looking up contact indices, and the flat array of contact history data. Contains
    only raw pointers and primitive types, making it safe to pass to kernels by value.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T>
struct ContactMemoryView
{
    /** @name Parameters */
    //@{
    /** \brief Pointer to flat array of contact history data */
    ContactHistory<T>* m_historyData;
    /** \brief Pointer to hash table entries */
    ContactEntry* m_table;
    /** \brief Total capacity of the hash table */
    uint m_capacity;
    /** \brief Pointer to the next index counter */
    uint* m_nextIndex;
    //@}

    /** @name Constructors */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Default constructor */
    __HOSTDEVICE__
    ContactMemoryView()
        : m_historyData(nullptr)
        , m_table(nullptr)
        , m_capacity(0)
        , m_nextIndex(nullptr)
    {
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Constructor with parameters
        @param historyData pointer to contact history data array
        @param table pointer to hash table entries
        @param capacity total capacity of the hash table
        @param nextIndex pointer to the next index counter */
    __HOSTDEVICE__
    ContactMemoryView(ContactHistory<T>* historyData,
                      ContactEntry*      table,
                      uint               capacity,
                      uint*              nextIndex)
        : m_historyData(historyData)
        , m_table(table)
        , m_capacity(capacity)
        , m_nextIndex(nextIndex)
    {
    }
    //@}

    /** @name Hash table operations */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Finds an existing contact index in the hash table
        @param key pair of component IDs (i,j) where i < j
        @param index output parameter for the contact state index
        @return true if found, false otherwise */
    __HOSTDEVICE__
    bool find(uint2 key, uint& index) const
    {
        if(m_capacity == 0 || m_table == nullptr)
            return false;

        uint h = primeHash(key) % m_capacity;

        // Linear probing
        for(uint i = 0; i < m_capacity; ++i)
        {
            uint                idx = (h + i) % m_capacity;
            const ContactEntry& e   = m_table[idx];

            // Empty slot found - key not in table
            if(e.m_valid == 0)
                return false;

            // Key found (either stale or active)
            if((e.m_valid == 1 || e.m_valid == 2) && e.m_key.x == key.x && e.m_key.y == key.y)
            {
#ifdef __CUDA_ARCH__
                // Ensure all writes are visible before reading (GPU only)
                __threadfence();
#endif
                index = e.m_index;
                return true;
            }
        }

        // Table full, key not found
        return false;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Finds an existing contact or inserts a new one
        @param key pair of component IDs (i,j) where i < j
        @param index output parameter for the contact state index
        @return true if found or inserted, false if table is full */
    __HOSTDEVICE__
    bool findOrInsert(uint2 key, uint& index)
    {
        if(m_capacity == 0 || m_table == nullptr || m_nextIndex == nullptr)
            return false;

        uint h = primeHash(key) % m_capacity;

        for(uint i = 0; i < m_capacity; ++i)
        {
            uint          idx = (h + i) % m_capacity;
            ContactEntry& e   = m_table[idx];

#ifdef __CUDA_ARCH__
            // GPU path
            // Empty slot - try to claim it atomically
            if(e.m_valid == 0)
            {
                // Atomic compare-and-swap: if e.m_valid is 0, set it to 2 (active)
                if(atomicCAS(&e.m_valid, 0u, 2u) == 0u)
                {
                    // Successfully claimed this slot
                    e.m_key = key;

                    // Assign new index atomically
                    uint newIndex = atomicAdd(m_nextIndex, 1u);
                    e.m_index     = newIndex;

                    // Ensure writes are visible to all threads before they can read
                    __threadfence();

                    index = newIndex;
                    return true;
                }
                // Another thread claimed it, continue probing
            }
            // Stale entry - reactivate it
            else if(e.m_valid == 1 && e.m_key.x == key.x && e.m_key.y == key.y)
            {
                // Mark as active for this timestep
                atomicExch(&e.m_valid, 2u);
                __threadfence();
                index = e.m_index;
                return true;
            }
            // Already active - just return the index
            else if(e.m_valid == 2 && e.m_key.x == key.x && e.m_key.y == key.y)
            {
                __threadfence();
                index = e.m_index;
                return true;
            }
#else
            // CPU path
            // Empty slot - claim it
            if(e.m_valid == 0)
            {
                e.m_valid = 2;  // Mark as active
                e.m_key   = key;
                e.m_index = (*m_nextIndex)++;
                index     = e.m_index;
                return true;
            }
            // Stale entry - reactivate it
            else if(e.m_valid == 1 && e.m_key.x == key.x && e.m_key.y == key.y)
            {
                e.m_valid = 2;  // Mark as active
                index     = e.m_index;
                return true;
            }
            // Already active - just return the index
            else if(e.m_valid == 2 && e.m_key.x == key.x && e.m_key.y == key.y)
            {
                index = e.m_index;
                return true;
            }
#endif
            // Different key, continue probing
        }

        // Hash table is full - cannot insert
        return false;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Removes a contact from the hash table
        @param key pair of component IDs (i,j) where i < j
        @return true if removed, false if not found */
    __HOSTDEVICE__
    bool remove(uint2 key)
    {
        if(m_capacity == 0 || m_table == nullptr)
            return false;

        uint h = primeHash(key) % m_capacity;

        for(uint i = 0; i < m_capacity; ++i)
        {
            uint          idx = (h + i) % m_capacity;
            ContactEntry& e   = m_table[idx];

            // Empty slot - key not found
            if(e.m_valid == 0)
                return false;

            // Key found - remove it
            if(e.m_key.x == key.x && e.m_key.y == key.y)
            {
#ifdef __CUDA_ARCH__
                // Mark as invalid atomically
                atomicExch(&e.m_valid, 0u);
#else
                // Mark as invalid
                e.m_valid = 0;
#endif
                return true;
            }
        }

        // Key not found
        return false;
    }
    //@}
};

// =================================================================================================
/** @brief Hash table manager for contact history tracking.

    This class manages both the hash table for fast lookups and the flat array of contact history
    data. It ensures both structures are properly sized and synchronized.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T, MemType M = MemType::HOST>
class ContactHashTable
{
private:
    /** @name Parameters */
    //@{
    /** \brief Flat array of contact history data */
    GrainsMemBuffer<ContactHistory<T>, M> m_historyData;
    /** \brief Hash table entries buffer */
    GrainsMemBuffer<ContactEntry, M> m_table;
    /** \brief Total capacity of the hash table */
    uint m_capacity;
    /** \brief Maximum number of contacts that can be stored */
    uint m_maxContacts;
    /** \brief Pointer to the next index counter (in device/host memory) */
    uint* m_nextIndex;
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    ContactHashTable();

    /** @brief Constructor with specified capacities
        @param hashCapacity number of entries in the hash table
        @param maxContacts maximum number of contacts to store */
    ContactHashTable(uint hashCapacity, uint maxContacts);

    /** @brief Destructor */
    ~ContactHashTable();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the pointer to the history data */
    const ContactHistory<T>* getHistoryData() const;

    /** @brief Gets the pointer to the history data (mutable) */
    ContactHistory<T>* getHistoryData();

    /** @brief Gets the pointer to the table data */
    const ContactEntry* getTable() const;

    /** @brief Gets the pointer to the table data (mutable) */
    ContactEntry* getTable();

    /** @brief Gets the capacity of the hash table */
    uint getCapacity() const;

    /** @brief Gets the maximum number of contacts */
    uint getMaxContacts() const;

    /** @brief Gets the current next index value (host-side read only) */
    uint getNextIndex() const;

    /** @brief Gets the pointer to the next index counter */
    uint* getNextIndexPointer();

    /** @brief Gets the pointer to the next index counter (const) */
    const uint* getNextIndexPointer() const;

    /** @brief Gets a complete view for passing to kernels */
    ContactMemoryView<T> getView();
    //@}

    /** @name Memory management methods */
    //@{
    /** @brief Allocates memory for both hash table and history data
        @param hashCapacity number of hash table entries to allocate
        @param maxContacts maximum number of contacts to store */
    void allocate(uint hashCapacity, uint maxContacts);

    /** @brief Frees memory used by hash table and history data */
    void deallocate();

    /** @brief Clears all entries in the hash table and resets history data */
    void clear();

    /** @brief Resets the index counter */
    void resetIndexCounter();

    /** @brief Performs mark-and-sweep cleanup: demotes active to stale, removes stale to empty.
        Should be called at the beginning of each timestep before contact detection */
    void markAndSweep();
    //@}
};

#endif
