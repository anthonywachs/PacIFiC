#ifndef _QUATERNION_HH_
#define _QUATERNION_HH_

#include "Vector3.hh"
#include <iostream>
#include "WriterXML.hh"

using namespace solid;
using namespace std;


class Quaternion;
class Matrix;
ostream& operator << (ostream& fileOut, Quaternion const& q );
istream& operator >> (istream& fileIn, Quaternion& q );


/** @brief The class Quaternion.

    A quaternion in a 3D space, i.e., a scalar w plus a Vector3 vector vqt
    as [ w, vqt ].

    @author Institut Francais du Petrole - 2000 - Modification
    @author A.WACHS  - 2019 - Modification */
// ============================================================================
class Quaternion
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    Quaternion();

    /** @brief Constructor with 2 scalar as input parameters q and d.
    Quaternion is initialized as [ d, (q,q,q) ]
    @param q value of all 3 components of the vector
    @param d value of the scalar */
    Quaternion( double q, double d=0.0 );

    /** @brief Constructor with a Vector3 vector vec and a scalar d. Quaternion
    is initialized as [ d, vec ]
    @param vec the Vector3 vector
    @param d value of the scalar */
    Quaternion( Vector3 const& vec, double d=0.0 );

    /** @brief Constructor with a vector given by its 3 components (x,y,z) and
    a scalar d. Quaternion is initialized as [ d, (x,y,z) ]
    @param x x-component of the vector
    @param y y-component of the vector
    @param z z-component of the vector
    @param d value of the scalar */
    Quaternion( double x, double y, double z, double d );

    /** @brief Copy constructor
    @param q copied Quaternion object */
    Quaternion( Quaternion const& q );

    /** @brief Destructor */
    ~Quaternion();
    //@}


    /**@name Methods */
    //@{
    /** @brief Returns the conjugate of the quaternion */
    Quaternion Conjugate() const;

    /** @brief Returns the inverse of the quaternion */
    Quaternion Inverse() const;

    /** @brief Multiplies the quaternion on the left by a vector lhs, i.e.,
    performs [ 0, lhs ] x this and return the product that is a quaternion
    @param lhs the left hand side vector */
    Quaternion multLeftVec( Vector3 const& lhs ) const;

    /** @brief Multiplies the quaternion on the right by another quaternion rhs,
    i.e., performs this x rhs, and return the vectorial part of this x rhs
    @param rhs the other quaternion */
    Vector3 multToVector3( Quaternion const& rhs ) const;

    /** @brief Multiplies the quaternion on the right by the conjugate of
    another quaternion rhs, i.e., perform this x rhs^t, and return the
    vectorial part of this x rhs^t
    @param rhs the other quaternion */
    Vector3 multConjugateToVector3( Quaternion const& rhs ) const;

    /** @brief Writes the object with a high precision format given by
    FORMAT16DIGITS defined in GrainsExec.hh
    @param fileOut output stream */
    void writeQuaternion( ostream &fileOut ) const;

    /** @brief Writes the object in binary format
    @param fileOut output stream */
    void writeQuaternion_binary( ostream &fileOut );

    /** @brief Reads the object in binary format
    @param StreamIN input stream */
    void readQuaternion_binary( istream &StreamIN );
    
    /** @brief Rotates a vector using the quaternion *this
    @param v The vector to be rotated */
    Vector3 rotateVector( Vector3 const& v ) const;    
    //@}


    /**@name Methods Get */
    //@{
    /** @brief Returns the value of the scalar part of the quaternion */
    double getdouble() const;

    /** @brief Returns the pointer to the vectorial part of the quaternion */
    Vector3 const* getVector3() const;
    //@}


    /**@name Methods Set */
    //@{
    /** @brief Sets the quaternion with a Vector3 vector vec and a scalar d.
    Quaternion is set to [ d, vec ]
    @param vec the Vector3 vector
    @param d value of the scalar */
    void setQuaternion( Vector3 const& vec, double d );

    /** @brief Sets the quaternion with a vector given by its 3 components
    (x,y,z) and a scalar d. Quaternion is set to [ d, (x,y,z) ]
    @param x x-component of the vector
    @param y y-component of the vector
    @param z z-component of the vector
    @param d value of the scalar */
    void setQuaternion( double x, double y, double z, double d );

    /** @brief Sets the scalar part of the quaternion
    @param d value of the scalar */
    void setdouble( double d );

    /** @brief Sets the vectorial part of the quaternion
    @param vec the Vector3 vector */
    void setVector3( Vector3 const& vec );

    /** @brief Sets the quaternion with a rotation matrix
    @param rot rotation matrix */
    void setQuaternion( Matrix const& rot );
    //@}

    /** @brief Builds a unit quaternion representing the rotation,from u to v.
    The input vectors need not be normalised.
    @param u First vector
    @param v Second vector */
    void setRotFromTwoVectors( Vector3 const& u, Vector3 const& v );
    //@}
    

    /**@name Operators */
    //@{
    /** @brief ith-component accessor: (0,1,2) for the vector components and 3
    for the scalar
    @param i index */
    double& operator [] ( int i );

    /** @brief ith-component accessor: (0,1,2) for the vector components and 3
    for the scalar
    @param i index */
    double operator [] ( int i ) const;

    /** @brief double product this x rhs of 2 quaternions
    @param rhs the other quaternion */
    Quaternion operator * ( Quaternion const& rhs ) const;

    /** @brief double product on the right of a quaternion by a vector
    [ 0, rhs ]
    @param rhs the Vector3 vector */
    Quaternion operator , ( Vector3 const& rhs ) const;

    /** @brief Product by a scalar
    @param d scalaire */
    Quaternion operator * ( double d );

    /** @brief Unitary operator -. Return a quaternion with negative elements */
    Quaternion operator-();

    /** @brief Sum of 2 quaternions
    @param rhs the other quaternion */
    Quaternion operator + ( Quaternion const& rhs ) const;

    /** @brief Subtraction of 2 quaternions, i.e., compute this - rhs
    @param rhs the other quaternion */
    Quaternion operator - ( Quaternion const& rhs );

    /** @brief Comparison operator
    @param rhs the other quaternion */
    bool operator == ( Quaternion const& rhs );

    /** @brief Difference operator
    @param rhs the other quaternion */
    bool operator != ( Quaternion const& rhs );

    /** @brief Equal operator to another Quaternion object
    @param rhs the other Quaternion object */
    Quaternion& operator = ( Quaternion const& rhs );

    /** @brief Equal operator to a scalar d, the result is [ d, (d,d,d) ]
    @param d scalar */
    Quaternion operator = ( const double d );

    /** @brief Unitary operator *= by a scalar
    @param d multiplication factor */
    Quaternion& operator *= ( double d );

    /** @brief Operator +=
    @param rhs the other quaternion */
    Quaternion& operator += ( Quaternion const& rhs );

    /** @brief Operator -=
    @param rhs the other quaternion */
    Quaternion& operator -= ( Quaternion const& rhs );
    //@}


    /**@name Friend methods */
    //@{
    /** @brief Output operator
    @param fileOut output stream
    @param q the quaternion */
    friend ostream& operator << ( ostream& fileOut, Quaternion const& q );

    /** @brief Input operator
    @param fileIn input stream
    @param q the quaternion */
    friend istream& operator >> ( istream& fileIn, Quaternion& q );

    /** @brief Returns the norm of the quaternion
    @param q the quaternion */
    friend double Norm( Quaternion const& q );

    /** @brief Returns the norm square of the quaternion
    @param q the quaternion */
    friend double Norm2( Quaternion const& q );
    //@}


    /**@name Parameters */
    //@{
    static size_t m_sizeofQuaternion; /** binary size of the object */
    //@}
    

  protected:
    /**@name Parameters */
    //@{
    double m_w; /**< double part of the quaternion */
    Vector3 m_vqt; /**< Vectorial part of the quaternion */
    //@}
};


/**@name Quaternion : External methods */
//@{
/** @brief Product by a scalar
@param d the scalar
@param rhs the quaternion */
Quaternion operator * ( double d, Quaternion const& rhs );

/** @brief double product on the left of a vector [ 0, lhs ] by a quaternion.
Returns a quaternion
@param lhs the left hand side vector
@param q the quaternion */
Quaternion operator , ( Vector3 const& lhs, Quaternion const& q );
//@}

#endif
