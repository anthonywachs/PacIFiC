#include "VertexBase.hh"


// ----------------------------------------------------------------------------
// Default constructor
VertexBase::VertexBase() 
  : m_base(0) 
{}




// ----------------------------------------------------------------------------
// Constructor with a pointer to the array of vertices as input parameter
VertexBase::VertexBase( void const* ptr ) 
  : m_base( ptr )  
{}




// ----------------------------------------------------------------------------
// Destructor
VertexBase::~VertexBase() 
{}
