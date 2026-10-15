#ifndef _QUATERNION_HH_
#define _QUATERNION_HH_

#include "Matrix3.hh"
#include "ReaderXML.hh"
#include "Vector3.hh"

// =================================================================================================
/** @brief The class Quaternion.

    A quaternion in a 3D space, i.e., a scalar w plus a Vector3 vector vqt as [ w, vqt ].

    @author F.PRADEL - Institut Francais du Petrole - 2000 - Modification
    @author A.WACHS  - 2019 - Modification
    @author A.Yazdani - 2024 - Major modification */
// =================================================================================================
template <typename T>
class alignas(16) Quaternion
{
protected:
    /** @name Parameters */
    //@{
    Vector3<T> m_vqt; /**< Vectorial part of the quaternion */
    T          m_w;   /**< scalar part of the quaternion */
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    __HOSTDEVICE__
    Quaternion() noexcept;

    /** @brief Constructor with 2 scalar as input parameters q and d.
        Quaternion is initialized as [ d, (q,q,q) ]
        @param q value of all 3 components of the vector
        @param w value of the scalar */
    __HOSTDEVICE__
    Quaternion(T q, T w = T(0)) noexcept;

    /** @brief Constructor with a Vector3 vec and a scalar d. Quaternion is
        initialized as [ d, vec ]
        @param vec the Vector3 vector
        @param w value of the scalar */
    __HOSTDEVICE__
    Quaternion(const Vector3<T>& vec, T w = T(0)) noexcept;

    /** @brief Constructor with a vector given by its 3 components (x,y,z) and a scalar d.
        Quaternion is initialized as [ d, (x,y,z) ].
        @param x x-component of the vector
        @param y y-component of the vector
        @param z z-component of the vector
        @param w value of the scalar */
    __HOSTDEVICE__
    Quaternion(T x, T y, T z, T w) noexcept;

    /** @brief Constructor with a buffer
        @param buffer buffer */
    __HOSTDEVICE__
    Quaternion(const T* buffer) noexcept;

    /** @brief Constructor from Euler angles (radians). Builds a quaternion from intrinsic Z-Y-X
        rotations: R = Rz(aZ) * Ry(aY) * Rx(aX).
        @param aX rotation about X (roll)
        @param aY rotation about Y (pitch)
        @param aZ rotation about Z (yaw) */
    __HOSTDEVICE__
    Quaternion(T aX, T aY, T aZ) noexcept;

    /** @brief Constructor with a rotation matrix
        @param rot rotation matrix */
    __HOSTDEVICE__
    Quaternion(const Matrix3<T>& rot) noexcept;

    /** @brief Copy constructor
        @param q copied Quaternion object */
    __HOSTDEVICE__
    Quaternion(const Quaternion<T>& q) noexcept;

    /** @brief Assign operator to another Quaternion object
        @param q rhs Quaternion object */
    __HOSTDEVICE__
    Quaternion<T>& operator=(const Quaternion<T>& q) noexcept;

    /** @brief Move constructor
        @param q moved Quaternion object */
    __HOSTDEVICE__
    Quaternion(Quaternion<T>&& q) noexcept;

    /** @brief Move assignment operator
        @param q moved Quaternion object */
    __HOSTDEVICE__
    Quaternion<T>& operator=(Quaternion<T>&& q) noexcept;

    /** @brief Constructor from an XML node
        @param root XML node */
    __HOST__
    Quaternion(DOMNode* root) noexcept;

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Quaternion() noexcept;
    //@}

    /** @name Get methods */
    //@{
    /** @brief Returns the pointer to the buffer */
    __HOSTDEVICE__
    const T* getBuffer() const noexcept;

    /** @brief Returns the vectorial part of the quaternion */
    __HOSTDEVICE__
    const Vector3<T>& getVector() const noexcept;

    /** @brief Returns the value of the scalar part of the quaternion */
    __HOSTDEVICE__
    const T& getScalar() const noexcept;
    //@}

    /**@name Set methods */
    //@{
    /** @brief Sets the vectorial part of the quaternion
        @param vec the Vector3 vector */
    __HOSTDEVICE__
    void setVector(const Vector3<T>& vec) noexcept;

    /** @brief Sets the scalar part of the quaternion
        @param w value of the scalar */
    __HOSTDEVICE__
    void setScalar(T w) noexcept;

    /** @brief Sets the quaternion with a Vector3 vector vec and a scalar d. Quaternion is set
        to [ d, vec ]
        @param vec the Vector3 vector
        @param w value of the scalar */
    __HOSTDEVICE__
    void setQuaternion(const Vector3<T>& vec, T w) noexcept;

    /** @brief Sets the quaternion with a vector given by its 3 components (x,y,z) and a scalar d.
        Quaternion is set to [ d, (x,y,z) ].
        @param x x-component of the vector
        @param y y-component of the vector
        @param z z-component of the vector
        @param w value of the scalar */
    __HOSTDEVICE__
    void setQuaternion(T x, T y, T z, T w) noexcept;

    /** @brief Sets the quaternion with a rotation matrix
        @param rot rotation matrix */
    __HOSTDEVICE__
    void setQuaternion(const Matrix3<T>& rot) noexcept;

    /** @brief Sets the quaternion from Euler angles (radians). Intrinsic Z-Y-X order:
        R = Rz(aZ) * Ry(aY) * Rx(aX).
        @param aX rotation about X (roll)
        @param aY rotation about Y (pitch)
        @param aZ rotation about Z (yaw) */
    __HOSTDEVICE__
    void setQuaternion(T aX, T aY, T aZ) noexcept;

    /** @brief Builds a unit quaternion representing the rotation, from u to v. The input vectors
        need not to be normalised.
        @param u First vector
        @param v Second vector */
    __HOSTDEVICE__
    void setRotFromTwoVectors(const Vector3<T>& u, const Vector3<T>& v) noexcept;
    //@}

    /**@name Methods */
    //@{
    /** @brief Converts the quaternion to a rotation matrix */
    __HOSTDEVICE__
    Matrix3<T> toMatrix() const noexcept;

    /** @brief Multiplies the quaternion on the right by another quaternion rhs, i.e., performs
        this x rhs, and return the vectorial part of this x rhs.
        @param q the other quaternion */
    __HOSTDEVICE__
    Vector3<T> multToVector3(const Quaternion<T>& q) const noexcept;
    //@}

    /**@name Operators */
    //@{
    /** @brief ith component accessor
        @param i component index */
    __HOSTDEVICE__
    T operator[](size_t i) const noexcept;

    /** @brief ith-component accessor: (0,1,2) for the vector components and 3 for the scalar -
        modifiable lvalue.
        @param i index */
    __HOSTDEVICE__
    T& operator[](size_t i) noexcept;
    //@}
};

/** @name External Methods - I/O methods */
//@{
/** @brief Output operator
    @param fileOut output stream
    @param q quaternion object */
template <typename T>
std::ostream& operator<<(std::ostream& fileOut, const Quaternion<T>& q);

/** @brief Input operator
    @param fileIn input stream
    @param q quaternion object */
template <typename T>
std::istream& operator>>(std::istream& fileIn, Quaternion<T>& q);
//@}

#endif
