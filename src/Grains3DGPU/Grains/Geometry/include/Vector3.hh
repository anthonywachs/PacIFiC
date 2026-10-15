#ifndef _VECTOR3_HH_
#define _VECTOR3_HH_

#include "Basic.hh"

// =================================================================================================
/** @brief The class Vector3.

    Vector/Point in a 3D space.

    @author A.Yazdani - 2023 - Construction */
// =================================================================================================
template <typename T>
class Vector3
{
protected:
    /**@name Parameters */
    //@{
    T m_comp[3]; /**< Array of 3 components */
    //@}

public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor
        @param def value of all 3 components */
    __HOSTDEVICE__
    Vector3(T def = T()) noexcept;

    /** @brief Constructor with the buffer
        @param buffer buffer */
    __HOSTDEVICE__
    Vector3(const T* buffer) noexcept;

    /** @brief Constructor with 3 components as inputs
        @param x 1st component
        @param y 2nd component
        @param z 3rd component*/
    __HOSTDEVICE__
    Vector3(T x, T y, T z) noexcept;

    /** @brief Copy constructor
        @param vec copied Vector3 object */
    __HOSTDEVICE__
    Vector3(const Vector3<T>& vec) noexcept;

    /** @brief Assign operator to another Vector3 object
        @param vec rhs Vector3 object */
    __HOSTDEVICE__
    Vector3<T>& operator=(const Vector3<T>& vec) noexcept;

    /** @brief Move constructor
        @param vec moved Vector3 object */
    __HOSTDEVICE__
    Vector3(Vector3<T>&& vec) noexcept;

    /** @brief Move assignment operator
        @param vec moved Vector3 object */
    __HOSTDEVICE__
    Vector3<T>& operator=(Vector3<T>&& vec) noexcept;

    /** @brief Constructor from an XML node
        @param root XML node */
    __HOST__
    Vector3(DOMNode* root) noexcept;

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Vector3() noexcept;
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the pointer to the buffer */
    __HOSTDEVICE__
    const T* getBuffer() const noexcept;
    //@}

    /** @name Set methods */
    //@{
    /** @brief Sets the vector to a 1D array of 3 values as input
        @param buffer the 1D array of values ordered as: 0=Vx, 1=Vy, 2=Vz */
    __HOSTDEVICE__
    void setValue(const T* buffer) noexcept;

    /** @brief Sets the components
        @param x the x component
        @param y the y component
        @param z the z component */
    __HOSTDEVICE__
    void setValue(const T x, const T y, const T z) noexcept;
    //@}

    /** @name Methods */
    //@{
    /** @brief Unitary nomalization operator */
    __HOSTDEVICE__
    void normalize() noexcept;

    /** @brief Returns a vector corresponding to the normalized vector */
    __HOSTDEVICE__
    Vector3<T> normalized() const noexcept;

    /** @brief set all components to zero */
    __HOSTDEVICE__
    void reset() noexcept;

    /** @brief Atomically add a vector to this vector (GPU safe)
        @param other vector to add atomically */
    __device__ void atomicAdd(const Vector3<T>& other) noexcept;
    //@}

    /** @name Operators */
    //@{
    /** @brief ith component accessor
        @param i component index */
    __HOSTDEVICE__
    T const& operator[](size_t i) const noexcept;

    /** @brief ith component accessor - modifiable lvalue
        @param i component index */
    __HOSTDEVICE__
    T& operator[](size_t i) noexcept;
    //@}
};

/** @name External Methods - I/O methods */
//@{
/** @brief Output operator
    @param fileOut output stream
    @param v vector */
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Vector3<T>& v);

/** @brief Input operator
    @param fileIn input stream
    @param v vector */
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Vector3<T>& v);
//@}

#endif
