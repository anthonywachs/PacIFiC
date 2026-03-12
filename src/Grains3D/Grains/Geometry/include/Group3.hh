#ifndef _GROUP3_HH_
#define _GROUP3_HH_

#include "Basic.hh"

#include <iostream>
#include <cmath>
using namespace std;


class Group3;
ostream &operator << (ostream &fileOut, const Group3 &objet);
istream &operator >> (istream &fileIn, Group3 &objet);


namespace solid
{
  /** @brief The class Group3.

  Geometric object with 3 components. From GJK Engine - A Fast and Robust GJK
  Implementation, Copyright (C) 1998  Gino van den Bergen.

  @author Institut Francais du Petrole - 2000 - Creation
  @author A.WACHS - 2019 - Update */
  // ==========================================================================
  class Group3
  {
    public:
      /**@name Constructors */
      //@{
      /** @brief Default constructor
      @param def value of all 3 components */
      inline Group3( double def = 0. )
      { m_comp[X] = m_comp[Y] = m_comp[Z] = def; }

      /** @brief Constructor with 3 components as inputs
      @param x 1st component
      @param y 2nd component
      @param z 3rd component */
      inline Group3( double x, double y, double z )
      { m_comp[X] = x; m_comp[Y] = y; m_comp[Z] = z; }

      /** @brief Copy constructor
      @param g copied Group3 object */
      inline Group3( Group3 const& g )
      { m_comp[X] = g.m_comp[X]; m_comp[Y] = g.m_comp[Y]; m_comp[Z] = g.m_comp[Z]; }

      /** @brief Destructor */
      inline ~Group3() {}
      //@}


      /**@name Get methods */
      //@{
      /** @brief Returns the const pointer to the array */
      inline double const* getValue() const { return ( m_comp ); }

      /** @brief Returns the pointer to the array */
      inline double* getValue() { return ( m_comp ); }

      /** @brief Returns the number of components, always return 3 */
      inline int size() const { return ( 3 ); }
      //@}


      /**@name Set methods */
      //@{
      /** @brief Modifies the value of the 3 components
      @param x 1st component
      @param y 2nd component
      @param z 3rd component */
      inline void setValue( double x, double y, double z )
      { m_comp[X] = x; m_comp[Y] = y; m_comp[Z] = z; }

      /** @brief Modifies the value of the 3 components
      @param g 3-element array */
      inline void setValue( double const* g )
      { m_comp[X] = g[X]; m_comp[Y] = g[Y]; m_comp[Z] = g[Z]; }

      /** @brief Nullifies all components */
      inline void reset()
      { m_comp[X] = m_comp[Y] = m_comp[Z] = 0.; }
      //@}


      /**@name Methods */
      //@{
      /** @brief Writes the object with a high precision format given by
      FORMAT16DIGITS defined in GrainsExec.hh
      @param fileOut output stream */
      void writeGroup3( ostream& fileOut ) const;

      /** @brief Writes the object in binary format
      @param fileOut output stream */
      void writeGroup3_binary( ostream& fileOut );

      /** @brief Reads the object in binary format
      @param StreamIN input stream */
      void readGroup3_binary( istream& StreamIN );

      /** @brief Rounds components to +-tol
      @param tol Tolerance for rounding. If not given, the default is EPSILON */
      inline void round( double tol = EPSILON )
      {
        m_comp[X] = fabs( m_comp[X] ) < tol ? 0. : m_comp[X];
        m_comp[Y] = fabs( m_comp[Y] ) < tol ? 0. : m_comp[Y];
        m_comp[Z] = fabs( m_comp[Z] ) < tol ? 0. : m_comp[Z];
      }
      //@}


      /**@name Operators */
      //@{
      /** @brief ith component accessor
      @param i component index */
      inline double& operator [] ( size_t i ) { return ( m_comp[i] ); }

      /** @brief ith component accessor
      @param i component index */
      inline double const& operator [] ( size_t i ) const { return ( m_comp[i] ); }

      /** @brief Unitary operator -. Returns an object with negative
      components */
      inline Group3 operator - () const
      { return ( Group3( -m_comp[X], -m_comp[Y], -m_comp[Z] ) ); }

      /** @brief double product
      @param g 2nd Group3 object */
      inline double operator * ( Group3 const& g ) const
      { return ( m_comp[X] * g.m_comp[X] + m_comp[Y] * g.m_comp[Y]
        + m_comp[Z] * g.m_comp[Z] ); }

      /** @brief Multiplication by a scalar of the form Group3 * scalar
      @param d multiplication factor */
      inline Group3 operator * ( double d ) const
      { return ( Group3( m_comp[X] * d, m_comp[Y] * d, m_comp[Z] * d ) ); }

      /** @brief Division by a scalar
      @param d division factor */
      inline Group3 operator / ( double d ) const
      { return ( Group3( m_comp[X] / d, m_comp[Y] / d, m_comp[Z] / d ) ); }

      /** @brief Addition
      @param g2 2nd Group3 object */
      inline Group3 operator + ( Group3 const& g2 ) const
      { return ( Group3( m_comp[X] + g2.m_comp[X],
        m_comp[Y] + g2.m_comp[Y], m_comp[Z] + g2.m_comp[Z] ) ); }

      /** @brief Subtraction
      @param g2 2nd object */
      inline Group3 operator - ( Group3 const& g2 ) const
      { return ( Group3( m_comp[X] - g2.m_comp[X],
        m_comp[Y] - g2.m_comp[Y], m_comp[Z] - g2.m_comp[Z] ) ); }

      /** @brief Comparaison operator
      @param g2 2nd Group3 object */
      inline bool operator == ( Group3 const& g2) const
      { return ( m_comp[X] == g2[X] && m_comp[Y] == g2[Y] && m_comp[Z] == g2[Z] ); }

      /** @brief Difference operator
      @param g2 2nd object */
      inline bool operator != ( Group3 const& g2 )
      { return ( ! ( *this == g2 ) ); }

      /** @brief Equal operator to another Group3 object
      @param g2 the other Group3 object */
      inline Group3& operator = ( Group3 const& g2 )
      {
        if ( &g2 != this )
        { m_comp[X] = g2.m_comp[X]; m_comp[Y] = g2.m_comp[Y]; m_comp[Z] = g2.m_comp[Z]; }
        return ( *this );
      }

      /** @brief Equal operator to a scalar, all components are equal to the
      same scalar
      @param v value set to all 3 components */
      inline void operator = ( double v )
      { m_comp[X] = m_comp[Y] = m_comp[Z] = v; }

      /** @brief Unitary operator *= by a scalar
      @param d multiplication factor */
      inline Group3& operator *= ( double d )
      { m_comp[X] *= d; m_comp[Y] *= d; m_comp[Z] *= d; return ( *this ); }

      /** @brief Unitary operator /= by a scalar
      @param d multiplication factor */
      inline Group3& operator /= ( double d )
      { m_comp[X] /= d; m_comp[Y] /= d; m_comp[Z] /= d; return ( *this ); }

      /** @brief Operator +=
      @param g2 2nd Group3 object */
      inline Group3& operator += ( Group3 const& g2 )
      { m_comp[X] += g2.m_comp[X]; m_comp[Y] += g2.m_comp[Y]; m_comp[Z] += g2.m_comp[Z];
        return ( *this ); }

      /** @brief Operator -=
      @param g2 2nd Group3 object */
      inline Group3& operator -= ( Group3 const& g2 )
      { m_comp[X] -= g2.m_comp[X]; m_comp[Y] -= g2.m_comp[Y]; m_comp[Z] -= g2.m_comp[Z];
        return ( *this ); }
      //@}


      /**@name Friend methods */
      //@{
      /** @brief Mixed product of 3 Group3 objects
      @param g1 1st Group3 object
      @param g2 2nd Group3 object
      @param g3 3rd Group3 object */
      friend double triple( Group3 const& g1,
	Group3 const& g2, Group3 const& g3 );

      /** @brief Multiplication by a scalar of the form scalar * Group3
      @param d multiplication factor
      @param g Group3 object */
      friend Group3 operator * ( double d, Group3 const& g );

      /** @brief Output operator
      @param fileOut output stream
      @param g Group3 object */
      friend ostream& operator << ( ostream& fileOut, Group3 const& g );

      /** @brief Input operator
      @param fileIn input stream
      @param g Group3 object */
      friend istream& operator >> ( istream& fileIn, Group3& g );
      //@}


      /**@name Parameters */
      //@{
      static size_t m_sizeofGroup3; /** binary size of the object */
      //@}


    protected:
      /**@name Parameters */
      //@{
      double m_comp[3]; /**< array of 3 components */
      //@}
  };

  /** @brief Mixed product of 3 Group3 objects */
  inline double triple( Group3 const& g1, Group3 const& g2, Group3 const& g3 )
  {
    return (
      g1[X] * ( g2[Y] * g3[Z] - g2[Z] * g3[Y] ) +
      g1[Y] * ( g2[Z] * g3[X] - g2[X] * g3[Z] ) +
      g1[Z] * ( g2[X] * g3[Y] - g2[Y] * g3[X] ) );
  }

  /** @brief Multiplication by a scalar of the form scalar * Group3 */
  inline Group3 operator * ( double d, Group3 const& g )
  { return ( Group3( g[X] * d, g[Y] * d, g[Z] * d ) ); }

  static Group3 OrigineGroup3; /**< Origine (0.,0.,0.)  */

} // end namespace solid

#endif
