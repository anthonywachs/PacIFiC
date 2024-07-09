#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "Group3.hh"


namespace solid
{
  size_t Group3::m_sizeofGroup3 = 3 * sizeof( double );

  // --------------------------------------------------------------------------
  // Default constructor
  Group3::Group3( double def )
  {
    m_comp[X] = m_comp[Y] = m_comp[Z] = def;
  }




  // --------------------------------------------------------------------------
  // Constructor with 3 components as inputs
  Group3::Group3( double x, double y, double z )
  {
    m_comp[X] = x;
    m_comp[Y] = y;
    m_comp[Z] = z;
  }




  // --------------------------------------------------------------------------
  // Copy constructor
  Group3::Group3( Group3 const& g )
  {
    m_comp[X] = g.m_comp[X];
    m_comp[Y] = g.m_comp[Y];
    m_comp[Z] = g.m_comp[Z];
  }




  // --------------------------------------------------------------------------
  // Destructor
  Group3::~Group3()
  {
  }




  // --------------------------------------------------------------------------
  // Returns the const pointer to the array
  double const* Group3::getValue() const
  {
    return ( m_comp );
  }




  // --------------------------------------------------------------------------
  // Returns the pointer to the array
  double* Group3::getValue()
  {
    return ( m_comp );
  }




  // --------------------------------------------------------------------------
  // Modifies the value of the 3 components
  void Group3::setValue( double x, double y, double z )
  {
    m_comp[X] = x;
    m_comp[Y] = y;
    m_comp[Z] = z;
  }




  // --------------------------------------------------------------------------
  // Modifies the value of the 3 components
  void Group3::setValue( double const* g )
  {
    m_comp[X] = g[X];
    m_comp[Y] = g[Y];
    m_comp[Z] = g[Z];
  }




  // --------------------------------------------------------------------------
  // Nullifies all components
  void Group3::reset()
  {
    m_comp[X] = m_comp[Y] = m_comp[Z] = 0.;
  }




  // --------------------------------------------------------------------------
  // Returns the number of components, always return 3
  int Group3::size() const
  {
    return ( 3 );
  }




  // --------------------------------------------------------------------------
  // ith component accessor
  double& Group3::operator[]( size_t i )
  {
    return ( m_comp[i] );
  }




  // --------------------------------------------------------------------------
  // ith component accessor
  double const& Group3::operator[]( size_t i ) const
  {
    return ( m_comp[i] );
  }




  // --------------------------------------------------------------------------
  // Unitary operator -. Return an object with negative components
  Group3 Group3::operator - () const
  {
    return ( Group3( - m_comp[X], - m_comp[Y], - m_comp[Z] ) );
  }




  // --------------------------------------------------------------------------
  // double product
  double Group3::operator * ( Group3 const& g ) const
  {
    return ( m_comp[X] * g.m_comp[X] + m_comp[Y] * g.m_comp[Y]
    	+ m_comp[Z] * g.m_comp[Z] );
  }




  // --------------------------------------------------------------------------
  // Multiplication by a scalar of the form Group3 * scalar
  Group3 Group3::operator * ( double d ) const
  {
    return ( Group3( m_comp[X] * d, m_comp[Y] * d, m_comp[Z] * d ) );
  }




  // --------------------------------------------------------------------------
  // Division by a scalar
  Group3 Group3::operator / ( double d ) const
  {
    return ( Group3( m_comp[X] / d, m_comp[Y] / d, m_comp[Z] / d ) );
  }




  // --------------------------------------------------------------------------
  // Addition
  Group3 Group3::operator + ( Group3 const& g2 ) const
  {
    return ( Group3( m_comp[X] + g2.m_comp[X],
    	m_comp[Y] + g2.m_comp[Y], m_comp[Z] + g2.m_comp[Z] ) );
  }




  // --------------------------------------------------------------------------
  // Subtraction
  Group3 Group3::operator - ( Group3 const& g2 ) const
  {
    return ( Group3( m_comp[X] - g2.m_comp[X],
    	m_comp[Y] - g2.m_comp[Y], m_comp[Z] - g2.m_comp[Z] ) );
  }




  // --------------------------------------------------------------------------
  // Comparaison operator
  bool Group3::operator == ( Group3 const& g2 ) const
  {
    return ( m_comp[X] == g2[X] && m_comp[Y] == g2[Y] && m_comp[Z] == g2[Z] );
  }




  // --------------------------------------------------------------------------
  // Difference operator
  bool Group3::operator != ( Group3 const& g2 )
  {
    return ( ! ( *this == g2 ) );
  }




  // --------------------------------------------------------------------------
  // Equal operator to another Group3 object
  Group3& Group3::operator = ( Group3 const& g2 )
  {
    if ( &g2 != this )
    {
      m_comp[X] = g2.m_comp[X];
      m_comp[Y] = g2.m_comp[Y];
      m_comp[Z] = g2.m_comp[Z];
    }
    return ( *this );
  }




  // --------------------------------------------------------------------------
  // Equal operator to a scalar, all components are equal to the same scalar
  void Group3::operator = ( double v )
  {
    m_comp[X] = m_comp[Y] = m_comp[Z] = v;
  }




  // --------------------------------------------------------------------------
  // Unitary operator *= by a scalar
  Group3& Group3::operator *= ( double d )
  {
    m_comp[X] *= d;
    m_comp[Y] *= d;
    m_comp[Z] *= d;
    return ( *this );
  }




  // --------------------------------------------------------------------------
  // Unitary operator /= by a scalar
  Group3& Group3::operator /= ( double d )
  {
    m_comp[X] /= d;
    m_comp[Y] /= d;
    m_comp[Z] /= d;
    return ( *this );
  }




  // --------------------------------------------------------------------------
  // Operator +=
  Group3& Group3::operator += ( Group3 const& g2 )
  {
    m_comp[X] += g2.m_comp[X];
    m_comp[Y] += g2.m_comp[Y];
    m_comp[Z] += g2.m_comp[Z];
    return ( *this );
  }




  // --------------------------------------------------------------------------
  // Operator -=
  Group3& Group3::operator -= ( Group3 const& g2 )
  {
    m_comp[X] -= g2.m_comp[X];
    m_comp[Y] -= g2.m_comp[Y];
    m_comp[Z] -= g2.m_comp[Z];
    return ( *this );
  }




  // --------------------------------------------------------------------------
  // Mixed product of 3 Group3 objects
  double triple( Group3 const& g1, Group3 const& g2, Group3 const& g3 )
  {
    return (
    	g1.m_comp[X] * ( g2.m_comp[Y] * g3.m_comp[Z]
		- g2.m_comp[Z] * g3.m_comp[Y] ) +
      	g1.m_comp[Y] * ( g2.m_comp[Z] * g3.m_comp[X]
		- g2.m_comp[X] * g3.m_comp[Z] ) +
      	g1.m_comp[Z] * ( g2.m_comp[X] * g3.m_comp[Y]
		- g2.m_comp[Y] * g3.m_comp[X] ) );
  }




  // --------------------------------------------------------------------------
  // Multiplication by a scalar of the form scalar * Group3
  Group3 operator * ( double d, Group3 const& g )
  {
    Group3 result(g);
    result *= d;
    return ( result );
  }




  // --------------------------------------------------------------------------
  // Output operator
  ostream& operator << ( ostream& fileOut, Group3 const& g )
  {
    fileOut << g[X] << " " << g[Y] << " " << g[Z];
    return ( fileOut );
  }




  // --------------------------------------------------------------------------
  // Input operator
  istream& operator >> ( istream& fileIn, Group3& g )
  {
    fileIn >> g[X] >> g[Y] >> g[Z];
    return ( fileIn );
  }




  // --------------------------------------------------------------------------
  // Writes the object with a high precision format given by
  // POSITIONFORMAT defined in GrainsExec.hh
  void Group3::writeGroup3( ostream& fileOut ) const
  {
    fileOut << GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  	m_comp[X] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  	m_comp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  	m_comp[Z] );
  }




  // --------------------------------------------------------------------------
  // Writes the object in binary format
  void Group3::writeGroup3_binary( ostream& fileOut )
  {
    fileOut.write( reinterpret_cast<char*>( &m_comp[X] ), sizeof(double) );
    fileOut.write( reinterpret_cast<char*>( &m_comp[Y] ), sizeof(double) );
    fileOut.write( reinterpret_cast<char*>( &m_comp[Z] ), sizeof(double) );
  }




  // --------------------------------------------------------------------------
  // Reads the object in binary format
  void Group3::readGroup3_binary( istream& StreamIN )
  {
    StreamIN.read( reinterpret_cast<char*>( &m_comp[X] ), sizeof(double) );
    StreamIN.read( reinterpret_cast<char*>( &m_comp[Y] ), sizeof(double) );
    StreamIN.read( reinterpret_cast<char*>( &m_comp[Z] ), sizeof(double) );
  }




  // ---------------------------------------------------------------------------
  // Rounds components to +-tol
  void Group3::round( double tol )
  {
    m_comp[X] = fabs( m_comp[X] ) < tol ? 0. : m_comp[X];
    m_comp[Y] = fabs( m_comp[Y] ) < tol ? 0. : m_comp[Y];
    m_comp[Z] = fabs( m_comp[Z] ) < tol ? 0. : m_comp[Z];
  }
}
