#include "GrainsExec.hh"
#include "Polygon.hh"
#include "AllComponents.hh"
#include "Point3.hh"
#include "Vector3.hh"
#include <fstream>
#include <sstream>
#include <vector>

typedef std::vector<unsigned int> IndexBuf;



// ----------------------------------------------------------------------------
// Constructor with an input stream, a number of vertices, a
// reference to the array of vertices and a reference to the array of indices
// as input parameters
Polygon::Polygon( istream& fileIn, int nb_point, VertexBase& ref, 
	IndexArray& ia ) 
  : Polytope( fileIn, nb_point, ref, ia )
  , m_cobound( NULL )
  , m_curr_vertex( 0 ) 
  , m_InertiaPoly( NULL ) 
  , m_surfaceArea( 0. )
{
  readFaces( fileIn );
}




// ----------------------------------------------------------------------------
// Copy constructor
Polygon::Polygon( Polygon const& copy ) : 
  Polytope( copy )
{
  m_cobound = copy.m_cobound;
  m_curr_vertex = copy.m_curr_vertex;
  m_InertiaPoly = NULL;
  m_surfaceArea = copy.m_surfaceArea;
}




// ----------------------------------------------------------------------------
// Destructor
Polygon::~Polygon() 
{
  // The inertia tensor is only constructed for the reference freely moving
  // rigid objects, i.e., the particles
  // If the particle is active, m_InertiaPoly is NULL
  if ( m_InertiaPoly ) delete m_InertiaPoly;
  
  // m_cobound is freed by the garbage collector mechanism
  // implemented in GrainsExec 
}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType Polygon::getConvexType() const 
{
  return ( POLYGON );
}




// ----------------------------------------------------------------------------
// Creates a polygon from an input stream
Polygon* Polygon::create( istream& fileIn ) 
{ 
  // Read polygon file name and open the corresponding file
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
  VertexBase* vertexbase = new VertexBase((void *)point);
  
  // Creation of the array of vertex indices in the array of vertices
  unsigned int *v = new unsigned int[nb_point];
  for(int i=0; i<nb_point; i++) v[i] = i;
  IndexArray* ia = new IndexArray( nb_point, v ); 
  delete [] v; 

  // Creation of the polygon
  Polygon *polygon = new Polygon( *PolyIN, nb_point, *vertexbase, *ia );
  polygon->m_fichPoly = fichpoly;
  PolyIN->close();
  
  // Objects m_base, m_cobound and m_index are created once per
  // particle type, other particles of same type have a pointer only to these
  // objects. They are all freed by the garbage collector mechanism implemented
  // in GrainsExec
  GrainsExec::addOnePolytopeRefPointBase( point, vertexbase );
  GrainsExec::addOnePolytopeNodeNeighbors( polygon->m_cobound );
  GrainsExec::addOnePolytopeNodeIndex( ia );  
  
  delete PolyIN;
  
  return ( polygon );
}




// ----------------------------------------------------------------------------
// Creates a polygon from an XML node
Polygon* Polygon::create( DOMNode* root ) 
{
  // Read polygon file name and open the corresponding file
  string fichpoly = ReaderXML::getNodeAttr_String( root, "Name" );
  ifstream PolyIN( fichpoly.c_str() );
  if ( PolyIN.is_open() ) 
    GrainsExec::m_additionalDataFiles.insert( fichpoly );
      
  // Number of vertices
  int nb_point;
  PolyIN >> nb_point >> nb_point;

  // Creation of the array of vertices
  Point3 *point = new Point3[nb_point];
  VertexBase *vertexbase = new VertexBase((void *)point);

  // Creation of the array of vertex indices in the array of vertices
  unsigned int *v = new unsigned int[nb_point];
  for(int i=0; i<nb_point; i++) v[i] = i;
  IndexArray* ia = new IndexArray( nb_point, v ); 
  delete [] v; 

  // Creation of the polygon
  Polygon *polygon = new Polygon( PolyIN, nb_point, *vertexbase, *ia );
  polygon->m_fichPoly = fichpoly;
  
  // Objects m_base, m_cobound and m_index are created once per
  // particle type, other particles of same type have a pointer only to these
  // objects. They are all freed by the garbage collector mechanism implemented
  // in GrainsExec
  GrainsExec::addOnePolytopeRefPointBase( point, vertexbase );
  GrainsExec::addOnePolytopeNodeNeighbors( polygon->m_cobound );
  GrainsExec::addOnePolytopeNodeIndex( ia );  
 
  return ( polygon );
}




// ----------------------------------------------------------------------
// Reads the face indexing and builds the polygon by calling BuildPolygon 
void Polygon::readFaces( istream &fileIn )
{
  char buffer[long_string];
  
  int nbedge;
  fileIn >> nbedge;
  fileIn.getline(buffer,sizeof(buffer)); // pour aller a la ligne suivante
  
  IndexArray* edge = new IndexArray[nbedge];
  IndexBuf facetind;
  
  // Read the face index - vertex index connectivity
  for(int i = 0;i < nbedge;i++) 
  {
    fileIn.getline(buffer,sizeof(buffer));
    istringstream lign(buffer);
    if (!fileIn.eof()) 
    {
      int integer;
      for(;;)
      {
	if (!lign.eof()) 
	{
	  lign >> integer;
	  facetind.push_back(integer);
	}
	else break;
      }
      new(&edge[i]) IndexArray(int(facetind.size()),&facetind[0]);
      facetind.erase(facetind.begin(), facetind.end());
    }
    else break;
  }
  
  BuildPolygon( edge[0].size(), edge );
  
  delete [] edge;
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool Polygon::BuildInertia( double *inertia, double *inertia_1 ) const
{
  std::copy(&m_InertiaPoly[0], &m_InertiaPoly[6], &inertia[0]);

  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  inertia_1[5] = 1.0 / inertia[5];

  return ( true );
}
  



// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference sphere,
// i.e., without applying any transformation
double Polygon::computeCircumscribedRadius() const 
{
  double d, ray = ((*this)[0]) * ((*this)[0]);
  for (int i = 1; i< numVerts(); i++) 
  {
    if ((d =  (*this)[i] * (*this)[i]) > ray)
      ray = d;
  }
  return ( sqrt( ray ) );
}




// ----------------------------------------------------------------------------
// Returns a clone of the polygon
Convex* Polygon::clone() const
{
  return ( new Polygon(*this) );
}




// ----------------------------------------------------------------------------
// Returns the polygon surface area
double Polygon::getVolume() const 
{
  return ( m_surfaceArea );
}




// ----------------------------------------------------------------------------
// Allocates the inertia tensor array and sets its component and the surface 
// area to 0
void Polygon::Initialisation()
{
  m_surfaceArea = 0.;
  m_InertiaPoly = new double[6];
  m_InertiaPoly[0] = m_InertiaPoly[1] = m_InertiaPoly[2] = m_InertiaPoly[3] = 
	m_InertiaPoly[4] = m_InertiaPoly[5] = 0.;
}




// ----------------------------------------------------------------------------
// Polygon support function, returns the support point P, i.e. the
// point on the surface of the sphere that satisfies max(P.v)
Point3 Polygon::support( Vector3 const& v ) const 
{
  double norm = Norm(v);
  if (norm > EPSILON) {
    const Point3 ptTmp = (*this)[m_curr_vertex];
    double h = ptTmp * v;
    double d;
    int ni = int(m_curr_vertex) < numVerts()-1 ? m_curr_vertex+1 : 0;
    if ((d = (*this)[ni] * v) > h) { 
      do { 
	h = d; m_curr_vertex = ni;
	if (++ni == numVerts()) ni = 0;
      }
      while ((d = (*this)[ni] * v) > h);
    }
    else {
      ni = m_curr_vertex ? m_curr_vertex-1 : numVerts()-1;
      while ((d = (*this)[ni] * v) > h) {
	h = d; m_curr_vertex = ni;
	if (ni) --ni; else ni = numVerts()-1;
      }
    }  
    return (*this)[m_curr_vertex];
  } else {
    return Point3();
  }
}




// ----------------------------------------------------------------------------
// Output operator
void Polygon::writeShape( ostream& fileOut ) const 
{
  fileOut << "*Polygon\n";
  fileOut << GrainsExec::extractFileName( m_fichPoly ) << endl;
}




// ----------------------------------------------------------------------------
// Input operator
void Polygon::readShape( istream& fileIn ) 
{
  cerr << "Program Error :\n"
       << "Polygon::readShape non accessible.\n";
  exit(3);
}




// ----------------------------------------------------------------------
// Constructs the polygon knowing its description, i.e., vertices and indices
// Note: the array edge actually has a length of 1 and contains a single
// element with all indices as the polygon perimeter is considered as a single
// edge 
void Polygon::BuildPolygon( int nbedge, IndexArray const* edge ) 
{  
  m_cobound = new IndexArray[numVerts()];
  
  IndexBuf* indexBuf = new IndexBuf[numVerts()];
  
  // Initialize surface area and moment of inertia tensor to 0
  Initialisation();
  
  // Comppute surface area and moment of inertia tensor to 0
  Point3 G_;
  Point3 P1;
  Point3 P2;
  Vector3 u;
  
  int i;
  for (i=0; i<nbedge-1;i++) G_ += (*this)[edge[0][i]];
  G_ /=(1.0*nbedge);
  
  for(i = 0;i < nbedge-1; i++) 
  { 
    // P1 and P2 are the 2 vertices of the edge
    P1 = (*this)[edge[0][i]];
    P2 = (*this)[edge[0][i+1]];

    // Check ordering   
    u = G_;
    if ( triple( P1 - G_, P2 - G_, u ) < 0.)
    {
      Point3 aux = P1;
      P1 = P2;
      P2 = aux;
    }

    computeSurfaceInertiaContrib( P1 - G_, P2 - G_ );
  }
  
  for (i = 0; i < numVerts(); ++i) 
    if (indexBuf[i].size()) 
      new(&m_cobound[i]) IndexArray(int(indexBuf[i].size()), &indexBuf[i][0]);
  
  m_curr_vertex = 0;
  while (indexBuf[m_curr_vertex].size() == 0) ++m_curr_vertex;
  
  delete [] indexBuf;
}




// ----------------------------------------------------------------------------
// Computes the contribution to inertia and surface area of an edge
// defined by 2 consecutive vertices relative to the center of mass position
void Polygon::computeSurfaceInertiaContrib( Point3 const& P, Point3 const& Q )
{
  Vector3 x2( P ), y2( Q - P );
  Point3 gtemp( ( P + Q ) * 1. / 3. );
  Point3 Pt( P );
  Point3 Qt( Q );
  Pt = Pt - gtemp;
  Qt = Qt - gtemp;
  Vector3 z2 = x2^y2; 
    
  double S = Norm(z2)/2.0;
  m_surfaceArea += S;
   
  double delta_12 = S/12.0;
  m_InertiaPoly[1]  = m_InertiaPoly[2] = m_InertiaPoly[4] = 0.0;
  m_InertiaPoly[0]  = m_InertiaPoly[3] = 1.0;
  m_InertiaPoly[5] += delta_12 * ( Pt * Pt + Qt * Qt + gtemp * gtemp )
  	+ S * ( gtemp * gtemp );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the polygon in a 
// Paraview format
int Polygon::numberOfCells_PARAVIEW() const
{
  return ( 1 );  
}




// ----------------------------------------------------------------------------
// Writes the polygon in a Paraview format
void Polygon::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  int ncorners = numVerts();

  int count = firstpoint_globalnumber + 1;
  for (int i=0;i<ncorners;++i)
  {
    connectivity.push_back( count );
    ++count;
  }  
  last_offset += ncorners;    
  offsets.push_back( last_offset );
  cellstype.push_back( 7 );
  
  firstpoint_globalnumber += ncorners + 1;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the polygon
bool Polygon::isIn( Point3 const& pt ) const
{
  // TO DO
  
  return ( false );
}  
