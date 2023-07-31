#include "STLTriangle.hh"


// ----------------------------------------------------------------------------
// Default constructor
STLTriangle::STLTriangle( )
  : m_id( 0 )
  , m_surfacearea( 0. )
{}




// ----------------------------------------------------------------------------
// Contructor with parameters
STLTriangle::STLTriangle( tuple<STLVertex*,STLVertex*,STLVertex*> ve, 
    	Vector3 const& ne, size_t const& ide )
  : m_v( ve )
  , m_n( ne )
  , m_id( ide )
{
  computeSurfaceArea();
}




// ----------------------------------------------------------------------------
// Destructor
STLTriangle::~STLTriangle( )
{}




// ----------------------------------------------------------------------------
// Computes and returns the area of the triangle
double STLTriangle::getSurfaceArea() const
{
  return ( m_surfacearea );
}




// ----------------------------------------------------------------------------
// Computes the surface area of the triangle
void STLTriangle::computeSurfaceArea()
{
  STLVertex *v1 = std::get<0>(m_v);
  STLVertex *v2 = std::get<1>(m_v);
  STLVertex *v3 = std::get<2>(m_v);

  Vector3 u = v2->m_p - v1->m_p;
  Vector3 w = v3->m_p - v1->m_p;
  
  m_surfacearea = 0.5 * Norm( u ^ w );
}
