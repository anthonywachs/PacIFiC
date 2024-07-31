#ifndef _MATRIX_HH_
#define _MATRIX_HH_

#include <assert.h>
#include "Vector3.hh"
#include "Quaternion.hh"

//@{
/// 3x3 double matrix
typedef double Mat3[3][3];
//@}


class Matrix;
ostream& operator << ( ostream& fileOut, Matrix const& m );
istream& operator >> ( istream& fileIn, Matrix& m );


/** @brief The class Matrix.

    A 3x3 matrix for 3D transformations. From GJK Engine - A Fast and Robust GJK
    Implementation, Copyright (C) 1998  Gino van den Bergen.

    @author A.WACHS - 2019 - Modification */
// ============================================================================
class Matrix
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor. Matrix is initialized to the identity matrix
    */
    Matrix();

    /** @brief Constructor with an 1D array of values as inputs
    @param m the 1D array of values containing the matrix coefficients ordered
    as 0=Mxx, 1=Mxy, 2=Mxz, 3=Myx, 4=Myy, 5=Myz, 6=Mzx, 7=Mzy, 8=Mzz */
    Matrix( double const* m );

    /** @brief Constructor with a quaternion as input to initialized the matrix
    as a rotation matrix. The quaternion is not required to be unitary as it
    will be normalized by dividing by its norm.
    @param q the quaternion */
    Matrix( Quaternion const& q );

    /** @brief Constructor of a diagonal matrix with the 3 diagonal coefficients
    as inputs
    @param x (1,1) coefficient
    @param y (2,2) coefficient
    @param z (3,3) coefficient */
    Matrix( double x, double y, double z );

    /** @brief Constructor with all the 9 coefficients as inputs
    @param xx (1,1) coefficient
    @param xy (1,2) coefficient
    @param xz (1,3) coefficient
    @param yx (2,1) coefficient
    @param yy (2,2) coefficient
    @param yz (2,3) coefficient
    @param zx (3,1) coefficient
    @param zy (3,2) coefficient
    @param zz (3,3) coefficient */
    Matrix( double xx, double xy, double xz,
	          double yx, double yy, double yz,
	          double zx, double zy, double zz );

    /** @brief Copy constructor
    @param other the copied matrix */
    Matrix( Matrix const& other );

    /** @brief Destructor */
    ~Matrix();
    //@}


    /** @name Methods Get */
    //@{
    /** @brief Returns a 3x3 array describing the matrix */
    Mat3& getValue();

    /** @brief Returns a 3x3 array describing the matrix */
    Mat3 const& getValue() const;
    //@}


    /** @name Methods Set */
    //@{
    /** @brief Sets the matrix to the identity matrix */
    void setIdentity();

    /** @brief Sets the matrix to a rotation matrix with a quaternion. The
    quaternion is not required to be unitary as it will be normalized by
    dividing by its norm.
    @param q the quaternion */
    void setRotation( Quaternion const& q );

    /** @brief Sets the matrix to a diagonal matrix
    @param x (1,1) coefficient
    @param y (2,2) coefficient
    @param z (3,3) coefficient */
    void setScaling( double x, double y, double z );

    /** @brief Sets the matrix with an 1D array of 9 values as inputs.
    !!! IMPORTANT !!! the 1D array must be organized as: 0=Mxx, 1=Mxy, 2=Mxz,
    3=Myx, 4=Myy, 5=Myz, 6=Mzx, 7=Mzy, 8=Mzz
    @param m the 1D array of values containing the matrix coefficients */
    void setValue( double const* m );

    /** @brief Sets the matrix with all the 9 coefficients as inputs
    @param xx (1,1) coefficient
    @param xy (1,2) coefficient
    @param xz (1,3) coefficient
    @param yx (2,1) coefficient
    @param yy (2,2) coefficient
    @param yz (2,3) coefficient
    @param zx (3,1) coefficient
    @param zy (3,2) coefficient
    @param zz (3,3) coefficient */
    void setValue( double xx, double xy, double xz,
		double yx, double yy, double yz,
		double zx, double zy, double zz );
    //@}


    /** @name Methods */
    //@{
    /** @brief Returns a matrix with all coefficients have the absolute value of
    the coefficients of the matrix */
    Matrix absolute() const;

    /** @brief Returns the adjoint matrix. The adjoint matrix is the
    transpose of the cofactor matrix */
    Matrix adjoint() const;

    /** @brief Returns the determinant of the matrix */
    double determinant() const;

    /** @brief Returns the inverse matrix */
    Matrix inverse() const;

    /** @brief Returns the scalar product of a column of the matrix and a vector
    @param i column number
    @param v the vector */
    double tdot( int i, Vector3 const& v ) const ;

    /** @brief Returns the transposed matrix */
    Matrix transpose() const;

    /** @brief Copies the matrix in a 1D array.
    !!! IMPORTANT !!! the 1D array is organized as: Mxx, Mxy, Mxz, Myx, Myy,
    Myz, Mzx, Mzy, Mzz
    @param vit 1D array where matrix coefficients are copied
    @param i start index to copy in the 1D array */
    void copyMatrix( double *vit, int i ) const;

    /** @brief Writes the object with a high precision format given by
    POSITIONFORMAT defined in GrainsExec.hh
    @param fileOut output stream */
    void writeMatrix( ostream& fileOut ) const;

    /** @brief Writes the object with a high precision format given by
    POSITIONFORMAT defined in GrainsExec.hh and the 2014 reload format
    @param fileOut output stream */
    void writeMatrix2014( ostream& fileOut ) const;

    /** @brief Writes the object in binary format with the 2014 reload format
    @param fileOut output stream */
    void writeMatrix2014_binary( ostream& fileOut );

    /** @brief Reads the object in binary format with the 2014 reload format
    @param StreamIN input stream */
    void readMatrix2014_binary( istream& StreamIN );

    /** @brief Multiplies to the right by a scaling matrix defined by its 3
    diagonal coefficients
    @param x x-scaling (0,0) coefficient
    @param y y-scaling (1,1) coefficient
    @param z z-scaling (2,2) coefficient */
    void multiplyByScalingMatrix( double x, double y, double z );

    /** @brief Rounds components to +-tol
    @param tol Tolerance for rounding. If not given, the default is EPSILON */
    void round( double tol = EPSILON );

    /** @brief Returns whether the matrix is diagonal */
    bool isDiagonal() const;

    /** @brief Returns whether the matrix is the identity matrix
    @param tol tolerance in checking the terms */
    bool isIdentity( double tol = EPSILON2 ) const;

    /** @brief Returns whether the matrix is a rotation matrix */
    bool isRotation( double tol = EPSILON ) const;
    
    /** @brief Returns the trace of the matrix */
    double trace() const;       
    //@}


    /**@name Operators */
    //@{
    /** @brief i-th row accessor
    @param i row number */
    Vector3& operator [] ( int i );

    /** @brief i-th row accessor
    @param i row number */
    Vector3 const& operator [] ( int i ) const;

    /** @brief Operator +=
    @param m the other matrix */
    Matrix& operator += ( Matrix const& m );

    /** @brief Operator *=
    @param m the other matrix */
    Matrix& operator *= ( Matrix const& m );

    /** @brief Equal operator to another Matrix object
    @param m the other Matrix object */
    Matrix& operator = ( Matrix const& m );
    
    /** @brief Operator /= by a float number
    @param d the float number */
    Matrix& operator /= ( double const& d );        
    //@}


    /**@name Methods friend */
    //@{
    /** @brief Output operator
    @param fileOut output stream
    @param m the matrix */
    friend ostream& operator << ( ostream& fileOut, Matrix const& m );

    /** @brief Input operator
    @param fileIn input stream
    @param m the matrix */
    friend istream& operator >> ( istream& fileIn, Matrix& m );
    //@}


    /**@name Parameters */
    //@{
    static size_t m_sizeofMatrix; /** binary size of the object */
    //@}


  protected:
    /**@name Parameters */
    //@{
    Mat3 m_elem; /**< the 2D 3x3 array containing the matrix coefficients */
    //@}
};


/** @name Matrix : External methods */
//@{
/** @brief Matrces sum
@param m1 first matrix
@param m2 second matrix */
Matrix operator + ( Matrix const& m1, Matrix const& m2 );

/** @brief Matrix-vector product
@param m the matrix
@param v the vector */
Vector3 operator * ( Matrix const& m, Vector3 const& v );

/** @brief Vector-matrix product
@param m the matrix
@param v the vector */
Vector3 operator * ( Vector3 const& v, Matrix const& m );

/** @brief Matrix-matrix product
@param m1 left matrix
@param m2 right matrix */
Matrix operator * ( Matrix const& m1, Matrix const& m2 );

/** @brief Scalar-matrix product
@param c the scalar
@param m the matrix */
Matrix operator * ( double c, Matrix const& m );

/** @brief Transposed matrix-matrix product
@param m1 left matrix that the transposed is taken
@param m2 right matrix */
Matrix multTransposeLeft( Matrix const& m1, Matrix const& m2 );

/** @brief Returns the transposed matrix
@param m the matrix */
Matrix transpose( Matrix const& m );

/** @brief Returns the adjoint matrix
@param m the matrix */
Matrix adjoint( Matrix const& m );

/** @brief Returns the inverse matrix
@param m the matrix */
Matrix inverse( Matrix const& m );

/** @brief Returns a matrix with all coefficients have the absolute value of
the coefficients of the matrix
@param m the matrix */
Matrix absolute( Matrix const& m );

/** @brief Returns the matrix determinant
@param m the matrix */
double determinant( Matrix const& m );

/** @brief Returns the matrix that rotates vector src to vector dest,
i.e. dest = mat * src
@param src the source vector
@param dest the destination vector */
Matrix getRotationMatrix( Vector3 const& src, Vector3 const& dest );
//@}

#endif
