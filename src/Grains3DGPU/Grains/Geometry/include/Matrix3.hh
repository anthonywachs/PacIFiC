#ifndef _MATRIX3_HH_
#define _MATRIX3_HH_

#include "Vector3.hh"

// =================================================================================================
/** @brief The class Matrix3.

    3x3 real matrix.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class Matrix3
{
protected:
    /**@name Parameters */
    //@{
    T m_comp[9]; /**< 3x3 array containing the matrix components */
    //@}

public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor. Matrix is initialized to the identity
        matrix */
    __HOSTDEVICE__
    Matrix3() noexcept;

    /** @brief Constructor with a 1D array of values as input
        @param buffer the 1D array of values containing the matrix components
        ordered 0=Mxx, 1=Mxy, 2=Mxz, 3=Myx, 4=Myy, 5=Myz, 6=Mzx, 7=Mzy, 8=Mzz */
    __HOSTDEVICE__
    Matrix3(T const* buffer) noexcept;

    /** @brief Constructor with 9 components as inputs
        @param xx (1,1) coefficient
        @param xy (1,2) coefficient
        @param xz (1,3) coefficient
        @param yx (2,1) coefficient
        @param yy (2,2) coefficient
        @param yz (2,3) coefficient
        @param zx (3,1) coefficient
        @param zy (3,2) coefficient
        @param zz (3,3) coefficient */
    __HOSTDEVICE__
    Matrix3(T xx, T xy, T xz, T yx, T yy, T yz, T zx, T zy, T zz) noexcept;

    /** @brief Constructor with 3 angles (in radians) as input parameters
        @param roll rotation angle around the X axis
        @param pitch rotation angle around the Y axis
        @param yaw rotation angle around the Z axis */
    __HOSTDEVICE__
    Matrix3(T roll, T pitch, T yaw) noexcept;

    /** @brief Copy constructor
        @param mat the copied matrix */
    __HOSTDEVICE__
    Matrix3(const Matrix3<T>& mat) noexcept;

    /** @brief Assign operator to another matrix
        @param mat rhs Matrix3 object */
    __HOSTDEVICE__
    Matrix3<T>& operator=(const Matrix3<T>& mat) noexcept;

    /** @brief Move constructor
        @param mat the moved matrix */
    __HOSTDEVICE__
    Matrix3(Matrix3<T>&& mat) noexcept;

    /** @brief Move assignment operator
        @param mat rhs Matrix3 object */
    __HOSTDEVICE__
    Matrix3<T>& operator=(Matrix3<T>&& mat) noexcept;

    /** @brief Constructor with an XML node
        @param root XML node */
    __HOST__
    Matrix3(DOMNode* root) noexcept;

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Matrix3() noexcept;
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the pointer to the buffer */
    __HOSTDEVICE__
    T const* getBuffer() const noexcept;
    //@}

    /** @name Set methods */
    //@{
    /** @brief Sets the matrix to a 1D array of 9 values as input
        @param buffer the 1D array of values ordered as:
        0=Mxx, 1=Mxy, 2=Mxz, 3=Myx, 4=Myy, 5=Myz, 6=Mzx, 7=Mzy, 8=Mzz */
    __HOSTDEVICE__
    void setValue(T const* buffer) noexcept;

    /** @brief Sets the matrix with all 9 components as inputs
        @param xx (1,1) coefficient
        @param xy (1,2) coefficient
        @param xz (1,3) coefficient
        @param yx (2,1) coefficient
        @param yy (2,2) coefficient
        @param yz (2,3) coefficient
        @param zx (3,1) coefficient
        @param zy (3,2) coefficient
        @param zz (3,3) coefficient */
    __HOSTDEVICE__
    void setValue(T xx, T xy, T xz, T yx, T yy, T yz, T zx, T zy, T zz) noexcept;
    //@}

    /** @name Operators */
    //@{
    /** @brief i-th row accessor
        @param i row number */
    __HOSTDEVICE__
    Vector3<T>& operator[](uint i) noexcept;

    /** @brief i-th row accessor
        @param i row number */
    __HOSTDEVICE__
    const Vector3<T>& operator[](uint i) const noexcept;

    /** @brief Element accessor
        @param i element index (0-8) */
    __HOSTDEVICE__
    T& operator()(uint i) noexcept;

    /** @brief Element accessor (const version)
        @param i element index (0-8) */
    __HOSTDEVICE__
    const T& operator()(uint i) const noexcept;

    /** @brief Element accessor
        @param i element index (0-2)
        @param j element index (0-2) */
    __HOSTDEVICE__
    T& operator()(uint i, uint j) noexcept;

    /** @brief Element accessor (const version)
        @param i element index (0-2)
        @param j element index (0-2) */
    __HOSTDEVICE__
    const T& operator()(uint i, uint j) const noexcept;
    //@}
};

/** @name External Methods - I/O methods */
//@{
/** @brief Output operator
    @param fileOut output stream
    @param m matrix object */
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Matrix3<T>& m);

/** @brief Input operator
    @param fileIn input stream
    @param m matrix object */
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Matrix3<T>& m);
//@}

#endif
