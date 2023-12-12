#include "Polytope.hh"
#include <new>


// ----------------------------------------------------------------------------
// Constructor with an input stream, a number of vertices, a
// reference to the array of vertices and a reference to the array of indices
// as input parameters
Polytope::Polytope( istream& fileIn, int nb_point, VertexBase& ref, 
	IndexArray& ia ) 
  : m_base( ref )
  , m_index( ia ) 
{ 
  // Read the vertex coordinates from input stream
  for(int i=0; i<nb_point; i++) fileIn >> ref[i];
}




// ----------------------------------------------------------------------------
// Copy constructor
Polytope::Polytope( Polytope const& copy ) 
  : m_base( copy.m_base )
  , m_index( copy.m_index )
{
  m_fichPoly = copy.m_fichPoly ; 
}


  

// ----------------------------------------------------------------------------
// Destructor
Polytope::~Polytope()
{
  // Note: m_base and m_index are freed by the garbage collector mechanism
  // implemented in GrainsExec 
}




// ----------------------------------------------------------------------------
// ith vertex accessor
Point3 const& Polytope::operator [] ( int i ) const 
{ 
  return ( m_base[m_index[i]] ); 
}




// ----------------------------------------------------------------------------
// ith vertex accessor
Point3& Polytope::operator [] ( int i ) 
{
  return ( m_base[m_index[i]] );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices of the polytope
int Polytope::numVerts() const 
{ 
  return ( m_index.size() );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the polytope
vector<Point3> Polytope::getEnvelope() const
{
  vector<Point3> envelope;
  for (int i=0; i<numVerts(); i++) 
  {
    Point3 const& p = (*this)[i];
    envelope.push_back( p );
  }
  return ( envelope );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners
int Polytope::getNbCorners() const
{
  return ( numVerts() );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the polytope in a Paraview format
int Polytope::numberOfPoints_PARAVIEW() const 
{ 
  // Center of mass + number of vertices
  return ( 1 + m_index.size() );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the sphere in a Paraview format 
void Polytope::write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, Vector3 const* translation ) const
{
  Point3 pp, oo;
  
  // Gravity center
  pp = transform( oo );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  
  // Corners
  for (int i=0; i<numVerts(); i++) 
  {
    pp = transform( (*this)[i] );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;    
  }
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the polytope in a Paraview format
list<Point3> Polytope::get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp, oo;
  
  // Gravity center
  pp = transform( oo );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );  
  
  // Corners  
  for (int i=0;i<numVerts();++i)
  {
    pp = transform( (*this)[i] );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }
  
  return ( ParaviewPoints ); 
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the polytope
bool Polytope::isIn( Point3 const& pt ) const
{
  cout << "Warning when calling Polytope::isIn(x,y,z) "
       << "\nShould not go into this method !\n"
       << "Need for an assistance ! Stop running !\n";
  return ( false );
}
