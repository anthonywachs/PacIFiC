#ifndef _TRANSFORM_HH_
#define _TRANSFORM_HH_

#include "Point3.hh"
#include "Matrix.hh"
#include "Vector3.hh"
#include "Quaternion.hh"
#include "Basic.hh"
#include "ReaderXML.hh"
#include "WriterXML.hh"


class Transform;
ostream &operator << ( ostream &fileOut, Transform const& t );
istream &operator >> ( istream &fileIn, Transform &t );


/** @brief The class Transform.

    A position/origin described by a point and an orientation (i.e. a 3D angular
    position) described by a 3x3 matrix. Note that the 3D matrix can also
    contain a scaling component and is hence not necessarily unitary. 
    From GJK Engine - A Fast and Robust GJK Implementation, Copyright (C) 1998 
    Gino van den Bergen. 

    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS  - 2019 - Modification */
// ============================================================================
class Transform 
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor. Origin is (0,0,0) and matrix is identity */
    Transform();

    /** @brief Constructor with origin coordinates as inputs and matrix is 
    identity
    @param gx origin x-coordinate 
    @param gy origin y-coordinate
    @param gz origin z-coordinate */
    Transform( double gx, double gy, double gz );

    /** @brief Constructor with a 1D array of 12 values as inputs containing the
    rotation matrix coefficients following by the origin coordinates
    @param m 1D array contanining 12 values */
    Transform( double const m[12] );

    /** @brief Copy constructor
    @param other the transformation to be copied */
    Transform( Transform const& other );

    /** @brief Destructor */
    ~Transform();
    //@}
  

    /**@name Methods */
    //@{ 
    /** @brief Returns true if the transformation type is IDENTITY */
    bool isIdentity() const;
  
    /** @brief Composition with a scaling transformation: this = this o scaling
    @param x scaling factor in x
    @param y scaling factor in y
    @param z scaling factor in z */
    void composeWithScaling( double x, double y, double z );  

    /** @brief Output operator (calls <<)
    @param fileOut output stream */
    void writeShape( ostream& fileOut ) const;
  
    /** @brief Writes the object with a high precision format given by
    FORMAT16DIGITS defined in GrainsExec.hh
    @param fileOut output stream */
    void writeTransform( ostream& fileOut ) const;
  
    /** @brief Writes the object with a high precision format given by
    FORMAT16DIGITS defined in GrainsExec.hh and the 2014 reload format
    @param fileOut output stream */
    void writeTransform2014( ostream& fileOut ) const;
  
    /** @brief Writes the object in binary format with the 2014 reload format
    @param fileOut output stream */
    void writeTransform2014_binary( ostream &fileOut );  
  
    /** @brief Reads the object with the 2014 reload format
    @param StreamIN input stream */
    void readTransform2014( istream& StreamIN ); 
  
    /** @brief Reads the object in binary format with the 2014 reload format
    @param StreamIN input stream */
    void readTransform2014_binary( istream& StreamIN );        

    /** @brief Input operator (calls >>)
    @param fileIn input stream */
    void readShape( istream& fileIn );

    /** @brief Composition on the left by a rotation described by a quaternion:
    this = rot(quaternion) o this (this first followed by rot(quaternion))
    @param q quaternion descriving the rotation */
    void composeLeftByRotation( Quaternion const& q );
  
    /** @brief Composition on the left by a rotation described by a transform:
    this = rot o this (this first followed by rot). This composition leaves the
    origin unchanged but does not check that rot is indeed a rotation 
    @param rot the other transformation describing a rotation */
    void composeLeftByRotation( Transform const& rot );  

    /** @brief Composition on the left by a translation:
    this = trans(vector) o this (this first followed by trans(vector))
    @param v translation vector */
    void composeLeftByTranslation( Vector3 const& v );
  
    /** @brief Copies the transformation in a 1D array
    !!! IMPORTANT !!! the 1D array is organized as: Mxx, Mxy, Mxz, Myx, Myy,
    Myz, Mzx, Mzy, Mzz, Ox, Oy, Oz 
    @param vit 1D array where coefficients are copied
    @param i start index to copy in the 1D array */
    void copyTransform( double* vit, int i ) const; 
  
    /** @brief Copies the transformation in a 1D array composed on the left by a
    translation (useful for periodic particles in parallel)
    @param vit 1D array where coefficients are copied
    @param i start index to copy in the 1D array 
    @param vec translation vecteur */
    void copyTransform( double* vit, int i, Vector3 const& vec ) const; 
  
    /** @brief Composition on the right by another affine transformation: 
    this = this o t (t first followed by this)
    @param t the other affine transformation */
    void composeRightByTransform( Transform const& t ); 
  
    /** @brief Composition on the left by another affine transformation: 
    this = t o this (this first followed by t)
    @param t the other affine transformation */
    void composeLeftByTransform( Transform const& t );         
    //@}
 

    /**@name Methods Get */
    //@{
    /** @brief Returns the matrix part of the transformation */
    Matrix const& getBasis() const;

    /** @brief Returns the transformation type */
    unsigned int getType() const;  

    /** @brief Returns a pointer to the origin of the transformation */
    Point3 const* getOrigin() const; 
    //@}  


    /**@name Methods Set */
    //@{
    /** @brief Sets the transformation to the identity */
    void setIdentity();

    /** @brief Sets the origin of the transformation
    @param pos origin of the transformation */
    void setOrigin( double const* pos );
  
    /** @brief Sets the origin of the transformation
    @param pos origin of the transformation */
    void setOrigin( Point3 const& pos );
  
    /** @brief Sets the origin of the transformation
    @param gx x-coordinate of the origin
    @param gy y-coordinate of the origin
    @param gz z-coordinate of the origin */
    void setOrigin( double gx, double gy, double gz );
  
    /** @brief Sets the matrix part of the transformation
    @param basis_ matrix part of the transformation */
    void setBasis( Matrix const& basis_ );    

    /** @brief Sets the transformation with an 1D array of 12 values as inputs  
    !!! IMPORTANT !!! the 1D array must be organized as: 0=Mxx, 1=Mxy, 2=Mxz,
    3=Myx, 4=Myy, 5=Myz, 6=Mzx, 7=Mzy, 8=Mzz, 9=Ox, 10=Oy, 11=Oz 
    @param m the 1D array of values containing the tranformation coefficients */
    void setValue( double const m[12] );
  
    /** @brief Sets the transformation to the inverse of another transformation 
    t
    @param t the other transformation */ 
    void setToInverseTransform( Transform const& t ); 
  
    /** @brief Composition of affine transformations: this = t2 o t1 (t1 first
    followed by t2)
    @param t1 1st affine transformation  
    @param t2 2nd affine transformation */
    void setToCompositionOfTransform( Transform const& t1, 
    	Transform const& t2 );
    //@}  
  
  
    /** @name XML */
    //@{
    /** @brief Initialization of the transformation from an XML node
    @param root the XML root node */
    void load( DOMNode* root );
    //@}


    /**@name Operators */
    //@{
    /** @brief Returns a point corresponding to applying the transformation to 
    the input point
    @param p input point */
    Point3 operator () ( Point3 const& p ) const;

    /** @brief Equal operator to another Transform object
    @param transform_ the other Transform object */
    Transform& operator = ( Transform const& transform_ );
    //@}
  

    /**@name Methods friend */
    //@{
    /** @brief Output operator
    @param fileOut output stream
    @param t the transformation */
    friend ostream& operator << ( ostream& fileOut, Transform const& t );

    /** @brief Input operator
    @param fileIn input stream
    @param t the transformation */
    friend istream& operator >> ( istream& fileIn, Transform& t );
    //@}


    /**@name Parameters */
    //@{
    static size_t m_sizeofTransform; /** binary size of the object */
    //@}


  private:
    enum { 
      /**@name Transformation types */
      //@{
      /// Identity
      IDENTITY = 0x00,
      /// Translation
      TRANSLATION = 0x01,
      /// Rotation
      ROTATION = 0x02,
      /// Scaling 
      SCALING = 0x04,
      /// Linear transformation
      LINEAR = ROTATION | SCALING,
      /// Affine transformation
      AFFINE = TRANSLATION | LINEAR
      //@}
    };

    /**@name Parameters */
    //@{
    Matrix m_basis; /**< Orientation matrix */
    Point3 m_origin; /**< Center of mass position */
    unsigned int m_type; /**< Transformation type */
    //@}
};

static Transform TransformIdentity; /**< identity transformation */

#endif
