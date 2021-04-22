
#include "IndexArray.H"

#include <algorithm>
using namespace std;

// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
IndexArray::IndexArray() : 
  indices(0), count(0) 
{
}

// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
IndexArray::IndexArray(int n) : 
  indices(new unsigned int[n]), count(n) 
{
}

// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
IndexArray::IndexArray(int n, const unsigned int v[]) : 
  indices(new unsigned int[n]), count(n) 
{ 
  copy(&v[0], &v[n], &indices[0]); 
}  

// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
IndexArray::~IndexArray() 
{ 
  delete [] indices; 
}
