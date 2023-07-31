#include "STLVertex.hh"


// ----------------------------------------------------------------------------
// Default constructor
STLVertex::STLVertex( )
{}




// ----------------------------------------------------------------------------
// Contructor with parameters
STLVertex::STLVertex( double x, double y, double z, Vector3 const& ne, 
	size_t const& ide )
  : m_p( x, y, z )
  , m_n( ne )
  , m_id( ide )
{}




// ----------------------------------------------------------------------------
// Destructor
STLVertex::~STLVertex( )
{}
