#include "IndexArray.hh"
#include <algorithm>
using namespace std;


// ----------------------------------------------------------------------------
// Default constructor
IndexArray::IndexArray() 
  : m_indices(0)
  , m_count(0) 
{}




// ----------------------------------------------------------------------------
// Constructor with the number of elements as an input parameter
IndexArray::IndexArray( int n ) 
  : m_indices( new unsigned int[n] )
  , m_count(n) 
{}




// ----------------------------------------------------------------------------
// Constructor with the number of elements and an array as input parameters
IndexArray::IndexArray( int n, unsigned int const v[] ) 
  : m_indices( new unsigned int[n] )
  , m_count(n) 
{ 
  copy( &v[0], &v[n], &m_indices[0] ); 
}  




// ----------------------------------------------------------------------------
// Destructor
IndexArray::~IndexArray() 
{ 
  delete [] m_indices; 
}
