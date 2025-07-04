#ifndef _VECTOR3_HH_
#define _VECTOR3_HH_

#include "Group3.hh"
#include "WriterXML.hh"

class Quaternion;

namespace solid
{
  /** @brief The class Vector3.

  Vector in a 3D space. From GJK Engine - A Fast and Robust GJK
  Implementation, Copyright (C) 1998  Gino van den Bergen.

  @author Institut Francais du Petrole - 1999 - Creation
  @author Institut Francais du Petrole - 2000 - Modification
  @author A.WACHS  - 2019 - Modification */
  // ==========================================================================
  class Vector3 : public Group3
  {
    public:
      /**@name Constructors */
      //@{
      /** @brief Default constructor
      @param def value of all 3 components */
      Vector3( double def = 0. );

      /** @brief Constructor with 3 components as inputs
      @param x 1st component
      @param y 2nd component
      @param z 3rd component*/
      Vector3( double x, double y, double z );

      /** @brief Copy constructor
      @param g copied Group3 object */
      Vector3( Vector3 const& g );

      /** @brief Copy constructor
      @param g copied Group3 object */
      Vector3( Group3 const& g );

      /** @brief Destructor */
      ~Vector3();
      //@}


      /** @name Methods */
      //@{
      /** @brief Determines the direction of lowest absolute component */
      int closestAxis() const;

      /** @brief Unitary nomalization operator */
      void normalize();

      /** @brief Returns a vector corresponding to the normalized vector */
      Vector3 normalized() const;

      /** @brief Rotation by an unitary quaternion
      @param q unitary quaternion corresponding to the rotation */
      void Rotate( Quaternion const& q );
      //@}


      /** @name Operators */
      //@{
      /** @brief Equal operator to another Vector3 object
      @param g2 the other Vector3 object */
      Vector3& operator = ( Vector3 const& g2 );    

      /** @brief Cross product this x rhv
      @param rhv 2nd Vector3 object */
      Vector3 operator ^ ( Vector3 const& rhv ) const;
      //@}


      /** @name Friend methods */
      //@{
      /** @brief Returns whether the vector norm is less than EPSILON2
      where EPSILON2 is defined in Basic.H
      @param v Vector3 object */
      friend bool approxZero( Vector3 const& v );

      /** @brief Returns the cosine of the angle between 2 Vector3 objects
      @param v1 1st Vector3 object
      @param v2 2nd Vector3 object */
      friend double cos( Vector3 const& v1, Vector3 const& v2 );

      /** @brief Returns the norm of the vector
      @param v the Vector3 object */
      friend double Norm( Vector3 const& v );

      /** @brief Returns the norm square of the vector
      @param v the Vector3 object */
      friend double Norm2( Vector3 const& v );
    //@}
  };

  /**@name Group3 : External methods */
  //@{
  /** @brief Returns the norm of the vector
  @param v the Vector3 object */
  double Norm( Vector3 const& v );
  //@}

  static Vector3 Vector3Null; /**< Vector3 (0.,0.,0.)  */
} // namespace solid

#endif
