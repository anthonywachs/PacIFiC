#include "VertexBase.H"

// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
VertexBase::VertexBase() :
  base(0) 
{
}

// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
VertexBase::VertexBase(const void *ptr) : 
  base(ptr) 
{
}

// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
VertexBase::~VertexBase() 
{
}
