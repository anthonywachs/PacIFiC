#include "GrainsExec.hh"
#include "Polyhedron.hh"
#include "AllComponents.hh"
#include "GrainsExec.hh"
#include <fstream>
#include <new> 
#include <sstream>
#include <vector>
using namespace std;
 
typedef vector<unsigned int> IndexBuf;


//-----------------------------------------------------------------------------
// Constructor with an input stream, a number of vertices, a
// reference to the array of vertices and a reference to the array of indices
// as input parameters
Polyhedron::Polyhedron( istream &fileIn, int nb_point, VertexBase& ref,
  	 IndexArray& ia ) 
  : Polytope( fileIn, nb_point, ref, ia )
  , m_cobound( NULL )
  , m_curr_vertex( 0 )
  , m_InertiaPoly( NULL )
  , m_VolumePoly( 0. )
{
  readFaces( fileIn );
}




// ----------------------------------------------------------------------
// Copy constructor
Polyhedron::Polyhedron( Polyhedron const& copy ) :
  Polytope( copy )
{  
  m_cobound = copy.m_cobound;
  m_curr_vertex = copy.m_curr_vertex;
  m_InertiaPoly = NULL;
  m_VolumePoly = copy.m_VolumePoly;
  m_allFaces = copy.m_allFaces;    
}




// ----------------------------------------------------------------------
// Destructor
Polyhedron::~Polyhedron() 
{
  // The inertia tensor is only constructed for the reference freely moving
  // rigid objects, i.e., the particles
  // If the particle is active, m_InertiaPoly is NULL
  if ( m_InertiaPoly ) delete m_InertiaPoly;
  
  // m_cobound and m_allFaces are freed by the garbage collector mechanism
  // implemented in GrainsExec
}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType Polyhedron::getConvexType() const 
{
  return ( POLYHEDRON );
}




// ----------------------------------------------------------------------
// Creates a polyhedron from an input stream
Polyhedron* Polyhedron::create( istream& fileIn ) 
{ 
  // Read polyhedron file name and open the corresponding file
  string fichpoly;
  fileIn >> fichpoly;
  fichpoly = GrainsExec::m_ReloadDirectory + "/" 
  	+ GrainsExec::extractFileName( fichpoly );
  ifstream *PolyIN = new ifstream( fichpoly.c_str() );
  if ( PolyIN->is_open() ) 
    GrainsExec::m_additionalDataFiles.insert( fichpoly );
  
  // Number of vertices
  int nb_point;
  *PolyIN >> nb_point >> nb_point;
  
  // Creation of the array of vertices
  Point3* point = new Point3[nb_point];
  VertexBase* vertexbase = new VertexBase( (void *)point );
  
  // Creation of the array of vertex indices in the array of vertices
  unsigned int *v = new unsigned int[nb_point];
  for(int i=0; i<nb_point; i++) v[i] = i;
  IndexArray* ia = new IndexArray( nb_point, v ); 
  delete [] v; 
    
  // Creation of the polyhedron
  Polyhedron *polyhedron = new Polyhedron( *PolyIN, nb_point, *vertexbase,
  	*ia );
  polyhedron->m_fichPoly = fichpoly;
  PolyIN->close();

  // Objects m_base, m_cobound, m_index and m_allFaces are created once per
  // particle type, other particles of same type have a pointer only to these
  // objects. They are all freed by the garbage collector mechanism implemented
  // in GrainsExec
  GrainsExec::addOnePolytopeRefPointBase( point, vertexbase );
  GrainsExec::addOnePolytopeNodeNeighbors( polyhedron->m_cobound ); 
  GrainsExec::addOnePolytopeNodeIndex( ia );
  GrainsExec::addOnePolyhedronFaceConnectivity( polyhedron->m_allFaces );
  
  delete PolyIN;
  
  return ( polyhedron );
}




// ----------------------------------------------------------------------------
// Creates a polyhedron from an XML node
Polyhedron* Polyhedron::create( DOMNode* root )
{
  // Read polyhedron file name and open the corresponding file
  string fichpoly = ReaderXML::getNodeAttr_String( root, "Name" ); 
  ifstream PolyIN( fichpoly.c_str() );
  if ( PolyIN.is_open() ) 
    GrainsExec::m_additionalDataFiles.insert( fichpoly );
  
  // Number of vertices
  int nb_point;
  PolyIN >> nb_point >> nb_point;

  // Creation of the array of vertices
  Point3* point = new Point3[nb_point];
  VertexBase* vertexbase = new VertexBase( (void *)point );

  // Creation of the array of vertex indices in the array of vertices
  unsigned int *v = new unsigned int[nb_point];
  for(int i=0; i<nb_point; i++) v[i] = i;
  IndexArray* ia = new IndexArray( nb_point, v ); 
  delete [] v; 
      
  // Creation of the polyhedron
  Polyhedron *polyhedron = new Polyhedron( PolyIN, nb_point, *vertexbase, 
  	*ia );
  polyhedron->m_fichPoly = fichpoly;

  // Objects m_base, m_cobound, m_index and m_allFaces are created once per
  // particle type, other particles of same type have a pointer only to these
  // objects. They are all freed by the garbage collector mechanism implemented
  // in GrainsExec
  GrainsExec::addOnePolytopeRefPointBase( point, vertexbase );
  GrainsExec::addOnePolytopeNodeNeighbors( polyhedron->m_cobound ); 
  GrainsExec::addOnePolytopeNodeIndex( ia );
  GrainsExec::addOnePolyhedronFaceConnectivity( polyhedron->m_allFaces );
  
  return ( polyhedron );
}




// ----------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool Polyhedron::BuildInertia( double* inertia, double* inertia_1 ) const
{
  std::copy(&m_InertiaPoly[0], &m_InertiaPoly[6], &inertia[0]);

  double determinant = inertia[0] * inertia[3] * inertia[5]
	- inertia[0] * inertia[4] * inertia[4]
	- inertia[5] * inertia[1] * inertia[1]
	- inertia[3] * inertia[2] * inertia[2]
	+ 2. * inertia[1] * inertia[2] * inertia[4];
  
  inertia_1[0] = ( inertia[3] * inertia[5] - inertia[4] * inertia[4] )
  	/ determinant;
  inertia_1[1] = ( inertia[2] * inertia[4] - inertia[1] * inertia[5] )
  	/ determinant;
  inertia_1[2] = ( inertia[1] * inertia[4] - inertia[2] * inertia[3] )
  	/ determinant;
  inertia_1[3] = ( inertia[0] * inertia[5] - inertia[2] * inertia[2] )
  	/ determinant;
  inertia_1[4] = ( inertia[1] * inertia[2] - inertia[0] * inertia[4] )
  	/ determinant;
  inertia_1[5] = ( inertia[0] * inertia[3] - inertia[1] * inertia[1] )
  	/ determinant;
 
  return ( true );
}




// ----------------------------------------------------------------------
// Returns the circumscribed radius of the reference polygon,
// i.e., without applying any transformation
double Polyhedron::computeCircumscribedRadius() const 
{
  double d , ray = ((*this)[0]) * ((*this)[0]);
  for (int i = 1; i< numVerts(); i++) 
  {
    if ((d =  (*this)[i] * (*this)[i]) > ray)
      ray = d;
  }
  return ( sqrt(ray) );
}




// -------------------------------------------------------------------
// Allocates the inertia tensor array and sets its component and the 
// volume to 0
void Polyhedron::Initialisation() 
{
  m_VolumePoly  = 0.0;
  m_InertiaPoly = new double[6];
  m_InertiaPoly[0] = m_InertiaPoly[1] = m_InertiaPoly[2] 
	= m_InertiaPoly[3] = m_InertiaPoly[4] = m_InertiaPoly[5] = 0;
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship
// between the face indices and the vertex indices
vector<vector<int> > const* Polyhedron::getFaces() const
{
  return ( m_allFaces );
}




// ----------------------------------------------------------------------
// Returns the polyhedron volume
double Polyhedron::getVolume() const 
{
  return ( m_VolumePoly );
}




// ----------------------------------------------------------------------
// Returns a clone of the polyhedron
Convex* Polyhedron::clone() const
{
  return ( new Polyhedron(*this) );
}




// ----------------------------------------------------------------------
// Polyhedron support function, returns the support point P, i.e. 
// the point on the surface of the polyhedron that satisfies max(P.v)
// Note: condition "d <= h" has been replaced by "d-h < eps" in relation to
// round error issues. For now, eps = 1.e-16 but this should be re-examined   
Point3 Polyhedron::support( Vector3 const& v ) const 
{
  double norm = Norm(v);
  if (norm > EPSILON) {

    if (m_cobound == NULL) {
      return support20(v);
    }

    int last_vertex = -1;
    double h = (*this)[m_curr_vertex] * v, d = 0.;
    for (;;) {
      IndexArray& curr_cobound = m_cobound[m_curr_vertex];
      int i = 0, n = curr_cobound.size();
      while (i != n && 
	     (curr_cobound[i] == last_vertex 
//	      || (d = (*this)[curr_cobound[i]] * v) <= h))
	      || (d = (*this)[curr_cobound[i]] * v) - h < 1.e-16 ))	      
	++i;	      
      if (i == n) break;
      last_vertex = m_curr_vertex;
      m_curr_vertex = curr_cobound[i];
      h = d;
    }
    return (*this)[m_curr_vertex];
  } else {
    return Point3();
  }
}




// ----------------------------------------------------------------------
// Output operator
void Polyhedron::writeShape( ostream& fileOut ) const 
{
  fileOut << "*Polyhedron " << GrainsExec::extractFileName( m_fichPoly ) 
  	<< " *END";
}




// ----------------------------------------------------------------------
// Input operator
void Polyhedron::readShape( istream& fileIn ) 
{
  cerr << "Program Error :\n"
       << "Polyhedron::readShape non accessible.\n";
  exit(3);
} 




// ----------------------------------------------------------------------------
// Constructs the polyhedron knowing its description, i.e.,
// vertices, faces and indices
void Polyhedron::BuildPolyhedron( int nbface, IndexArray const* face ) 
{  
  m_cobound = new IndexArray[numVerts()];
 
  IndexBuf* indexBuf = new IndexBuf[numVerts()];
  
  // Initialize volume and moment of inertia tensor to 0
  Initialisation();

  // Compute volume and moment of inertia tensor
  Point3 G_ ;
  int i, j, k;
  for(i=0; i<nbface; i++)
  {   
    for(j=0, k=face[i].size()-1; j<face[i].size(); k=j++)
      indexBuf[face[i][k]].push_back(face[i][j]);
    
    if ( face[i].size() == 3 )
    {
      GrainsExec::computeVolumeInertiaContrib( (*this)[face[i][0]], 
      	(*this)[face[i][1]], (*this)[face[i][2]], 
	m_VolumePoly, m_InertiaPoly );
    }
    else
    {
      G_.reset();
      for(j=0; j<face[i].size(); ++j) G_ += (*this)[face[i][j]];
      G_ /= face[i].size();
      for(j=0, k=face[i].size() - 1; j<face[i].size(); k=j++) 
        GrainsExec::computeVolumeInertiaContrib( G_, (*this)[face[i][k]], 
		(*this)[face[i][j]], m_VolumePoly, m_InertiaPoly);
    }
  }
  
  for (i = 0; i < numVerts(); ++i) 
    if (indexBuf[i].size()) 
      new(&m_cobound[i]) IndexArray(int(indexBuf[i].size()), &indexBuf[i][0]);
      
  m_curr_vertex = 0;
  while (indexBuf[m_curr_vertex].size() == 0) 
    ++m_curr_vertex;
  
  delete [] indexBuf;
}




// ----------------------------------------------------------------------
// Reads the face indexing and builds the polygon by calling BuildPolygon
void Polyhedron::readFaces( istream &fileIn )
{
  char buffer[long_string];
  
  int nbface;
  fileIn >> nbface;
  fileIn.getline( buffer, sizeof(buffer) ); 

  IndexArray *face = new IndexArray[nbface];
  IndexBuf facetind;
  
  // Read the face index - vertex index connectivity
  for(int i = 0;i < nbface;i++) 
  {
    fileIn.getline( buffer, sizeof(buffer) );
    istringstream lign( buffer );
    if ( !fileIn.eof() ) 
    {
      int integer;
      for(;;)
      {
	if (!lign.eof()) 
	{
	  lign >> integer;
	  facetind.push_back( integer );
	}
	else break;
      }
      new(&face[i]) IndexArray( int(facetind.size()), &facetind[0] );
      facetind.erase( facetind.begin(), facetind.end() );
    }
    else break;
  }

  BuildPolyhedron( nbface, face );
  
  // Store the face index - vertex index connectivity
  m_allFaces = new vector< vector<int> >;
  m_allFaces->reserve( nbface );
  for(int i=0;i<nbface;i++) 
  { 
    vector<int> ww( face[i].size() );
    for (int j=0;j<face[i].size();++j) ww[j] = face[i][j];
    m_allFaces->push_back( ww );
  }
  
  delete [] face;
}




// ----------------------------------------------------------------------------
// Polyhedron support function when the vertex neighbor indexing is
// not known, returns the support point P, i.e. the point on the surface of 
// the polyhedron that satisfies max(P.v). Note: Never used so far, not even 
// clear when to use it, was used in orginal SOLID-2.0 software
Point3 Polyhedron::support20( Vector3 const& v ) const
{
  int c = 0;
  double h = (*this)[0] * v, d;
  for (int i=1; i<numVerts(); ++i) {
    if ((d = (*this)[i] * v) > h) { c=i; h=d; }
  }
  return ( (*this)[c] );
}				    




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the
// polyhedron in a Paraview format
int Polyhedron::numberOfCells_PARAVIEW() const
{
  int ncorners = numVerts(), ncells = 0;

  // Tetrahedron or box
  if ( ncorners == 4 || ncorners == 8 ) ncells = 1;
  // Icosahedron
  else if ( ncorners == 12 && m_allFaces->size() == 20 ) ncells = 20;
  // Prism
  else if ( ncorners == 6 && m_allFaces->size() == 5 ) ncells = 1; 
  // Octahedron 
  else if ( ncorners == 6 && m_allFaces->size() == 8 ) ncells = 8; 
   // Dodecahedron 
  else if ( ncorners == 20 && m_allFaces->size() == 12 ) ncells = 36; 
  // Trancoctahedron 
  else if ( ncorners == 24 && m_allFaces->size() == 14 ) ncells = 44; 
  // General: not implemented yet !!
  else
  {
  
  }   
    
  return ( ncells );  
}




// ----------------------------------------------------------------------------
// Writes the polyhedron in a Paraview format
void Polyhedron::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  int ncorners = numVerts();
  
  // Tetrahedron or box or prism
  if ( ncorners == 4 || ncorners == 8 || 
  	( ncorners == 6 && m_allFaces->size() == 5 ) )
  {
    int count = firstpoint_globalnumber + 1;
    for (int i=0;i<ncorners;++i)
    {
      connectivity.push_back(count);
      ++count;
    }  
    last_offset += ncorners;
    offsets.push_back(last_offset);
    if ( ncorners == 4 ) cellstype.push_back(10);
    else if ( ncorners == 8 ) cellstype.push_back(12);
    else cellstype.push_back(13);
 	
    firstpoint_globalnumber += ncorners + 1;
  }
  // Icosahedron or octahedron
  // The icosahedron is split into 20 tetrahedrons using the center of mass and
  // 3 vertices on each triangular face
  // The octahedron is split into 8 tetrahedrons using the center of mass and
  // 3 vertices on each triangular face  
  else if ( ( ncorners == 12 && m_allFaces->size() == 20 ) ||
	  ( ncorners == 6 && m_allFaces->size() == 8 ) )
  {
    size_t nbface = m_allFaces->size();
    for (size_t i=0;i<nbface;++i)
    {
      connectivity.push_back( firstpoint_globalnumber );
      size_t nbv = (*m_allFaces)[i].size();
      for (size_t j=0;j<nbv;++j) 
        connectivity.push_back( firstpoint_globalnumber 
		+ (*m_allFaces)[i][j] + 1 );
      last_offset += 4;
      offsets.push_back(last_offset);
      cellstype.push_back(10);
    }

    firstpoint_globalnumber += ncorners + 1;       
  }
  // Dodecahedron
  // Each pentagonal face is split into 3 triangles
  // The dodecahedron is split into 12*3=36 tetrahedrons using the center of 
  // mass and each triangular sub-face of each pentagonal face
  else if ( ncorners == 20 && m_allFaces->size() == 12 )
  {
    size_t nbface = m_allFaces->size();
    for (size_t i=0;i<nbface;++i)
    {
      // First tetrahedron
      connectivity.push_back( firstpoint_globalnumber );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][0] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][1] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][2] + 1 );
      last_offset += 4;
      offsets.push_back(last_offset);
      cellstype.push_back(10);
 
      // Second tetrahedron
      connectivity.push_back( firstpoint_globalnumber );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][2] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][3] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][0] + 1 );
      last_offset += 4;
      offsets.push_back(last_offset);
      cellstype.push_back(10);
	
      // Third tetrahedron
      connectivity.push_back( firstpoint_globalnumber );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][3] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][4] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][0] + 1 );
      last_offset += 4;
      offsets.push_back(last_offset);
      cellstype.push_back(10);
    }
    
    firstpoint_globalnumber += ncorners + 1;  
  }
  // Trancoctahedron
  // Each pentagonal face is split into 3 triangles
  // The trancoctahedron is split into 8*4 + 6*2 = 44 tetrahedrons using the 
  // center of mass and each triangular sub-face of each pentagonal face  
  else if ( ncorners == 24 && m_allFaces->size() == 14 )
  {
    size_t nbface = m_allFaces->size();
    for (size_t i=0;i<8;++i)
    {
      // First tetrahedron
      connectivity.push_back( firstpoint_globalnumber );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][0] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][1] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][2] + 1 );
      last_offset += 4;
      offsets.push_back(last_offset);
      cellstype.push_back(10);
 
      // Second tetrahedron
      connectivity.push_back( firstpoint_globalnumber );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][2] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][3] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][0] + 1 );
      last_offset += 4;
      offsets.push_back(last_offset);
      cellstype.push_back(10);
	
      // Third tetrahedron
      connectivity.push_back( firstpoint_globalnumber );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][3] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][4] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][0] + 1 );
      last_offset += 4;
      offsets.push_back(last_offset);
      cellstype.push_back(10);

      // Fourth tetrahedron
      connectivity.push_back( firstpoint_globalnumber );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][4] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][5] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][0] + 1 );
      last_offset += 4;
      offsets.push_back(last_offset);
      cellstype.push_back(10);
    }
    
    for (size_t i=8;i<nbface;++i)
    {
      // First tetrahedron
      connectivity.push_back( firstpoint_globalnumber );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][0] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][1] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][2] + 1 );
      last_offset += 4;
      offsets.push_back(last_offset);
      cellstype.push_back(10);
 
      // Second tetrahedron
      connectivity.push_back( firstpoint_globalnumber );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][2] + 1 );
       connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][3] + 1 );
      connectivity.push_back( firstpoint_globalnumber
                + (*m_allFaces)[i][0] + 1 );
      last_offset += 4;
      offsets.push_back(last_offset);
      cellstype.push_back(10);
    }
    
    firstpoint_globalnumber += ncorners + 1;  
  }
  // General: not implemented yet !!
  else
  {
  
  }
}




// ----------------------------------------------------------------------------
// Writes the polyhedron in a STL format
void Polyhedron::write_convex_STL( ostream& f, Transform const& transform ) 
	const
{
  int ncorners = numVerts();

  // Box
  if ( ncorners == 8 && m_allFaces->size() == 6 )    
  {
    Point3 pp;
    Point3 GC = transform(pp);
    vector< Point3 > FC(4,pp);
    
    for (int i=0;i<int(m_allFaces->size());++i)
    {
      assert( (*m_allFaces)[i].size() == 4 ) ;      
      pp = 0.25 * ( (*this)[(*m_allFaces)[i][0]]
      	+ (*this)[(*m_allFaces)[i][1]] + (*this)[(*m_allFaces)[i][2]] 
	+ (*this)[(*m_allFaces)[i][3]] ) ;
      Point3 FaceCenterT = transform(pp);
      Vector3 outward_normal = FaceCenterT - GC;
      outward_normal.normalize(); 
      for (int j=0;j<4;++j)
        FC[j] = transform((*this)[(*m_allFaces)[i][j]]); 
	
      // Divide each faces into 2 triangles
      // Triangle 0
      f << "  facet normal " << outward_normal[X]
      	<< " " << outward_normal[Y]
      	<< " " << outward_normal[Z] << endl;
      f << "    outer loop" << endl;
      f << "      vertex " << FC[0][X] << " " << FC[0][Y] << " " 
      	<< FC[0][Z] << endl;
      f << "      vertex " << FC[1][X] << " " << FC[1][Y] << " " 
      	<< FC[1][Z] << endl;
      f << "      vertex " << FC[2][X] << " " << FC[2][Y] << " " 
      	<< FC[2][Z] << endl;
      f << "    endloop" << endl;
      f << "  endfacet" << endl;  
      // Triangle 1
      f << "  facet normal " << outward_normal[X]
      	<< " " << outward_normal[Y]
      	<< " " << outward_normal[Z] << endl;
      f << "    outer loop" << endl;
      f << "      vertex " << FC[0][X] << " " << FC[0][Y] << " " 
      	<< FC[0][Z] << endl;
      f << "      vertex " << FC[2][X] << " " << FC[2][Y] << " " 
      	<< FC[2][Z] << endl;
      f << "      vertex " << FC[3][X] << " " << FC[3][Y] << " " 
      	<< FC[3][Z] << endl;
      f << "    endloop" << endl;
      f << "  endfacet" << endl;      
    }  
  }
  else
  {
    cout << "Warning for this Convex with ncorners = " << ncorners << " the "
       << "method Polyhedron::write_convex_STL() is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
    exit(10);
  }	  
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the polyhedron
bool Polyhedron::isIn( Point3 const& pt ) const
{
  bool isIn = false;

  // Recall that by definition, center of mass is at the origin in the 
  // reference configuration
  Point3 G_, centerOfMass ;
  size_t i, j, k, nbface = m_allFaces->size();
  for (i=0; i<nbface && !isIn; i++)
  {       
    if ( (*m_allFaces)[i].size() == 3 )
    {
      isIn = GrainsExec::isPointInTetrahedron( (*this)[(*m_allFaces)[i][0]], 
      	(*this)[(*m_allFaces)[i][1]], (*this)[(*m_allFaces)[i][2]], 
	centerOfMass, pt );
    }
    else
    {
      size_t nptsface = (*m_allFaces)[i].size();
      G_.reset();
      for (j=0; j<nptsface; ++j) 
        G_ += (*this)[(*m_allFaces)[i][j]];
      G_ /= double(nptsface);
      for (j=0, k=nptsface - 1; j<nptsface && !isIn; k=j++) 
        isIn = GrainsExec::isPointInTetrahedron( G_, 
		(*this)[(*m_allFaces)[i][k]], (*this)[(*m_allFaces)[i][j]], 
		centerOfMass, pt );
    }
  }
  
  return ( isIn );
}




// ----------------------------------------------------------------------------
// Returns the bounding volume to polyhedron
BVolume* Polyhedron::computeBVolume( unsigned int type ) const
{
  // It is typically not straightforward to find the bounding volume to a 
  // general polyhedron. The method we use here is not by any means the optimal 
  // (most fitted) bounding volume to the given polyhedron.
  // Here, we simply loop over all vertices to find those with maximum abs value
  // in each direction. Persumably, for regular polyhedron it is going to be
  // enough.
  BVolume* bvol = NULL;
  Point3 pt;
  Vector3 extent( 0., 0., 0. );
  for ( int i = 0; i < numVerts(); i++ )
  {
    pt = m_base[i];
    extent[X] = fabs( pt[X] ) > extent[X] ? fabs( pt[X] ) : extent[X];
    extent[Y] = fabs( pt[Y] ) > extent[Y] ? fabs( pt[Y] ) : extent[Y];
    extent[Z] = fabs( pt[Z] ) > extent[Z] ? fabs( pt[Z] ) : extent[Z];
  }
  if ( type == 1 ) // OBB
  {
    bvol = new OBB( extent, Matrix() );
  }
  else if ( type == 2 ) // OBC
  {
    double xy = fabs( extent[X] - extent[Y] );
    double xz = fabs( extent[X] - extent[Z] );
    double yz = fabs( extent[Y] - extent[Z] );
    // pick from xy and xz, store to xy
    int zAxis = xy < xz ? Z : Y;
    xy = xy < xz ? xy : xz;
    // pick from xy and yz
    zAxis = xy < yz ? zAxis : X;

    Vector3 e( 0., 0., 0. );
    e[ zAxis ] = 1.;
    double h = 2. * extent[zAxis];
    double r = sqrt( extent[X] * extent[X] + 
                     extent[Y] * extent[Y] +
                     extent[Z] * extent[Z] -
                     extent[zAxis] * extent[zAxis] );

    bvol = new OBC( r, h, e );
  }

  return( bvol );
}



// ----------------------------------------------------------------------------
// Performs advanced comparison of the two polyhedrons and returns whether 
// they match
bool Polyhedron::equalType_level2( Convex const* other ) const
{
  // We know that other points to a Polyhedron, we dynamically cast it to 
  // actual type
  Polyhedron const* other_ = dynamic_cast<Polyhedron const*>(other);
  
  double lmin = computeCircumscribedRadius();    

  // Check the number of vertices
  int ncorners = numVerts();  
  bool same = ( ncorners == other_->numVerts() );
  
  // Check vertex coordinates
  for (int i=0;i<ncorners && same;++i) 
    same = ( (*this)[i].DistanceTo( (*other_)[i] ) < LOWEPS * lmin );     
    
  // Check face numbering
  if ( same )
  {
    size_t nbface = m_allFaces->size();
    for (size_t i=0;i<nbface;i++) 
      for (size_t j=0;j<(*m_allFaces)[i].size();++j) 
        same = ( (*m_allFaces)[i][j] == (*other_->m_allFaces)[i][j] );
  }
  
  return ( same );
} 




// ----------------------------------------------------------------------------
// Writes the polyhedron in an OBJ format
void Polyhedron::write_convex_OBJ( ostream& f, Transform  const& transform,
    	size_t& firstpoint_number ) const
{
  Point3 pp;
  int ncorners = numVerts();
  
  // Vertices
  for (int i=0; i<ncorners; i++) 
  {
    pp = transform( (*this)[i] );
    f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[Z] ) << " " << endl;	  
  }  
  
  // Faces 
  size_t nbface = m_allFaces->size();
  for (size_t i=0; i<nbface; i++) 
  {
    f << "f";
    for (size_t j=0;j<(*m_allFaces)[i].size();++j) 
      f << " " << firstpoint_number + (*m_allFaces)[i][j];
    f << endl;
  }
  
  firstpoint_number += size_t(ncorners);  			  
}
