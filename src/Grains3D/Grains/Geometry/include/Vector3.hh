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
      inline Vector3( double def = 0. ) : Group3( def ) {}

      /** @brief Constructor with 3 components as inputs
      @param x 1st component
      @param y 2nd component
      @param z 3rd component*/
      inline Vector3( double x, double y, double z ) : Group3( x, y, z ) {}

      /** @brief Copy constructor
      @param g copied Group3 object */
      inline Vector3( Vector3 const& g ) : Group3( g ) {}

      /** @brief Copy constructor
      @param g copied Group3 object */
      inline Vector3( Group3 const& g ) : Group3( g ) {}

      /** @brief Destructor */
      inline ~Vector3() {}
      //@}


      /** @name Methods */
      //@{
      /** @brief Determines the direction of lowest absolute component */
      int closestAxis() const;

      /** @brief Unitary nomalization operator */
      inline void normalize() { *this /= Norm( *this ); }

      /** @brief Returns a vector corresponding to the normalized vector */
      inline Vector3 normalized() const { return ( *this / Norm( *this ) ); }

      /** @brief Rotation by an unitary quaternion
      @param q unitary quaternion corresponding to the rotation */
      void Rotate( Quaternion const& q );
      //@}


      /** @name Operators */
      //@{
      /** @brief Equal operator to another Vector3 object
      @param g2 the other Vector3 object */
      inline Vector3& operator = ( Vector3 const& g2 )
      {
        if ( &g2 != this )
        { m_comp[X] = g2.m_comp[X]; m_comp[Y] = g2.m_comp[Y]; m_comp[Z] = g2.m_comp[Z]; }
        return (*this);
      }

      /** @brief Cross product this x rhv
      @param rhv 2nd Vector3 object */
      inline Vector3 operator ^ ( Vector3 const& rhv ) const
      {
        return ( Vector3( m_comp[1] * rhv.m_comp[2] - m_comp[2] * rhv.m_comp[1],
          - m_comp[0] * rhv.m_comp[2] + m_comp[2] * rhv.m_comp[0],
            m_comp[0] * rhv.m_comp[1] - m_comp[1] * rhv.m_comp[0] ) );
      }
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

  /** @brief Returns the norm of the vector */
  inline double Norm( Vector3 const& v )
  { return ( sqrt( v.m_comp[X] * v.m_comp[X]
    + v.m_comp[Y] * v.m_comp[Y] + v.m_comp[Z] * v.m_comp[Z] ) ); }

  /** @brief Returns the norm square of the vector */
  inline double Norm2( Vector3 const& v )
  { return ( v.m_comp[X] * v.m_comp[X] + v.m_comp[Y] * v.m_comp[Y]
    + v.m_comp[Z] * v.m_comp[Z] ); }

  /** @brief Returns whether the vector norm is less than EPSILON2 */
  inline bool approxZero( Vector3 const& v )
  { return ( Norm2( v ) < EPSILON2 ); }

  static Vector3 Vector3Null; /**< Vector3 (0.,0.,0.)  */
} // namespace solid

#endif
