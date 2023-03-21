#include "GrainsMPIWrapper.hh"
#include "STLObstacle.hh"
#include "LinkedCell.hh"
#include "Grains.hh"
#include "ContactBuilderFactory.hh"
#include "Particle.hh"
#include "PointContact.hh"
#include "Cell.hh"
#include "GrainsExec.hh"
#include "Memento.hh"
#include "ContactForceModel.hh"
#include "PointC.hh"
#include <sstream>
#include <limits>
#include <assert.h>
using namespace std;


// ----------------------------------------------------------------------------
// Constructor with name as input parameter
STLObstacle::STLObstacle( const string &s )
  : SimpleObstacle( s )
  , m_npls( 0 )
{
  m_ObstacleType = "STLObstacle";
}




// ----------------------------------------------------------------------------
// Constructor with name as input parameter
STLObstacle::STLObstacle( const string &s, string const& filename )
  : SimpleObstacle( s )
{
  m_ObstacleType = "STLObstacle";

  // The STL obstacle does not have a shape per se, its shape is made
  // of the STL triangles. Hence we defines its shape by a point
  // corresponding to its center of mass position (same as in CompositeObstacle
  // and in CompositeParticle)
  m_geoRBWC = new RigidBodyWithCrust( new PointC(), Transform() );
    
  // Read STL file and construct the triangulation
  readSTL( filename );

  // Set the bounding box
  list<STLVertex*>::iterator il;
  Point3 pmin( 1.e20, 1.e20, 1.e20 ), pmax( -1.e20, -1.e20, -1.e20 );
  for (il = m_allSTLVertices.begin(); il != m_allSTLVertices.end(); ++il)  
  {
    pmin[X] = min( (*il)->m_p[X], pmin[X] );
    pmin[Y] = min( (*il)->m_p[Y], pmin[Y] );    
    pmin[Z] = min( (*il)->m_p[Y], pmin[Y] );    
    pmax[X] = max( (*il)->m_p[X], pmax[X] );
    pmax[Y] = max( (*il)->m_p[Y], pmax[Y] );    
    pmax[Z] = max( (*il)->m_p[Y], pmax[Y] );        
  }
  m_obstacleBox.setValue( pmin, pmax );
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
STLObstacle::STLObstacle( DOMNode *root )
  : SimpleObstacle( "" )
{
  m_ObstacleType = "STLObstacle";
  assert( root != NULL );

  // The STL obstacle does not have a shape per se, its shape is made
  // of the STL triangles. Hence we defines its shape by a point
  // corresponding to its center of mass position (same as in CompositeObstacle
  // and in CompositeParticle)
  m_geoRBWC = new RigidBodyWithCrust( new PointC(), Transform() );
  
  m_name = ReaderXML::getNodeAttr_String( root, "name" );
  
  // Material
  DOMNode* materiau_ = ReaderXML::getNode( root, "Material" );
  m_materialName = ReaderXML::getNodeValue_String( materiau_ );
  ContactBuilderFactory::defineMaterial( m_materialName, true );

  // Obstacle to transfer to the fluid
  DOMNode* statut = ReaderXML::getNode( root, "Statut" );
  if ( statut )
    m_transferToFluid = ReaderXML::getNodeAttr_Int( statut, "ToFluid" );

  // m_obstacleBox = TO DO
  m_LinkUpdate_frequency = 1;
  m_LinkUpdate_counter = 0; 
  
  // File name
  DOMNode* file_ = ReaderXML::getNode( root, "File" ); 
  string filename = ReaderXML::getNodeAttr_String( file_, "Name" ); 
  
  // Read STL file and construct the triangulation
  readSTL( filename );
  
  // Translation of the STL obstacle to adjust its position
  DOMNode* translat = ReaderXML::getNode( root, "Translation" );  
  if ( translat )
  {
    Vector3 STLTranslation;
    STLTranslation[X] = ReaderXML::getNodeAttr_Double( translat, "X" );
    STLTranslation[Y] = ReaderXML::getNodeAttr_Double( translat, "Y" );
    STLTranslation[Z] = ReaderXML::getNodeAttr_Double( translat, "Z" );
    list<STLVertex*>::iterator il;
    for (il = m_allSTLVertices.begin(); il != m_allSTLVertices.end(); ++il)  
    {
      (*il)->m_p[X] += STLTranslation[X];
      (*il)->m_p[Y] += STLTranslation[Y];
      (*il)->m_p[Z] += STLTranslation[Z];     
    }     
  } 
  
  // Set the bounding box
  list<STLVertex*>::const_iterator il;
  Point3 pmin( 1.e20, 1.e20, 1.e20 ), pmax( -1.e20, -1.e20, -1.e20 );
  for (il = m_allSTLVertices.begin(); il != m_allSTLVertices.end(); ++il)  
  {
    pmin[X] = min( (*il)->m_p[X], pmin[X] );
    pmin[Y] = min( (*il)->m_p[Y], pmin[Y] );    
    pmin[Z] = min( (*il)->m_p[Y], pmin[Y] );    
    pmax[X] = max( (*il)->m_p[X], pmax[X] );
    pmax[Y] = max( (*il)->m_p[Y], pmax[Y] );    
    pmax[Z] = max( (*il)->m_p[Y], pmax[Y] );        
  }
  m_obstacleBox.setValue( pmin, pmax );       
}




// ----------------------------------------------------------------------------
// Destructor
STLObstacle::~STLObstacle()
{
  m_llvls.clear();
  m_llvns.clear();
  list<STLVertex*>::iterator il;
  for (il=m_allSTLVertices.begin();il!=m_allSTLVertices.end();il++)
    delete *il;
  m_allSTLVertices.clear();
  m_allSTLTriangles.clear();
}




// ----------------------------------------------------------------------------
// Returns whether the component is an STL obstacle ? */
bool STLObstacle::isSTLObstacle() const
{
  return ( true ); 
}




// ----------------------------------------------------------------------------
// Moves the simple obstacle and returns a list of moved obstacles (here itself)
list<SimpleObstacle*> STLObstacle::Move( double time,
	double dt, bool const& b_deplaceCine_Comp,
        bool const& b_deplaceF_Comp )
{
  // TO DO

  list<SimpleObstacle*> obstacleEnDeplacement;
  return ( obstacleEnDeplacement );
}






// ----------------------------------------------------------------------------
// Contact between a simple obstacle and a component. If contact
// exists, computes the contact force and torque and adds to each component
void STLObstacle::InterAction( Component* voisin,
	double dt, double const& time, LinkedCell* LC )
{
  cout << "STLObstacle::InterAction" << endl;
  
  try {
  list<ContactInfos*>  listContactInfos;

  // Search all contact points between the STL triangulation and the component
  // and store them in the list listContactInfos
  SearchContact( voisin, dt, time, LC, listContactInfos );

  // Loop over all contact points and compute the contact force & torque
  int nbContact = int(listContactInfos.size());
  for ( list<ContactInfos*>::iterator il=listContactInfos.begin();
      il!=listContactInfos.end(); il++ )
  {
    LC->addToContactsFeatures( time, (*il)->ContactPoint );

    if ( ContactBuilderFactory::contactForceModel(
		(*il)->p0->getMaterial(), (*il)->p1->getMaterial() )
      		->computeForces( (*il)->p0, (*il)->p1, (*il)->ContactPoint,
		LC, dt, nbContact ) )
    {
      (*il)->p0->getMasterComponent()->addToCoordinationNumber( 1 );
      (*il)->p1->getMasterComponent()->addToCoordinationNumber( 1 );
    }
    delete *il;
  }

  // Free the list
  listContactInfos.clear();
  }
  catch (const ContactError&) {
    throw ContactError();
  }
}




// ----------------------------------------------------------------------------
// Searches and stores all contact points between two components
void STLObstacle::SearchContact( Component* voisin, double dt,
      double const& time, LinkedCell *LC,
      list<ContactInfos*>& listContact )
{
  cout << "STLObstacle::SearchContact" << endl;
  bool found = false; 

  // Search for contact points between the component and the triangles
  if ( found )
  {
    ContactInfos* result = NULL ;
    PointContact closestPoint;
    result = new struct ContactInfos;
    result->ContactPoint = closestPoint;
    result->p0 = this;
    result->p1 = voisin;
    listContact.push_back( result );
  } 
  
  // Search for contact points between the component and the vertices
  if ( found )
  {
    ContactInfos* result = NULL ;
    PointContact closestPoint;
    result = new struct ContactInfos;
    result->ContactPoint = closestPoint;
    result->p0 = this;
    result->p1 = voisin;
    listContact.push_back( result );
  }   
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another component. 
bool STLObstacle::isContact( Component const* voisin ) const
{
  bool contact = false;
    
  // TO DO

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another
// component accounting for crust thickness
bool STLObstacle::isContactWithCrust( Component const* voisin ) const
{
  return ( STLObstacle::isContact( voisin ) );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric proximity with another
// component in the sense of whether their respective bounding boxes overlap 
bool STLObstacle::isClose( Component const* voisin ) const
{
  bool contact = false;
    
  // TO DO

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric proximity with another
// component in the sense of whether their respective bounding boxes minus 
// their crust thickness overlap
bool STLObstacle::isCloseWithCrust( Component const* voisin ) const
{
  bool contact = false;
    
  // TO DO

  return ( contact );
}



// ----------------------------------------------------------------------------
// Reloads the composite obstacle and links it to the higher level
// obstacle in the obstacle tree
void STLObstacle::reload( Obstacle& mother, istream& file )
{
  // TO DO
}




// ----------------------------------------------------------------------------
// Rotates the obstacle with a quaternion
void STLObstacle::Rotate( Quaternion const& rotation )
{
  // TO DO
}




// ----------------------------------------------------------------------------
// Translates the obstacle
void STLObstacle::Translate( Vector3 const& translation )
{
  // TO DO
}




// ----------------------------------------------------------------------------
// Outputs the simple obstacle for reload
void STLObstacle::write( ostream& fileSave ) const
{
  // TO DO
}




// ----------------------------------------------------------------------------
// Returns the maximum of the absolute value of the obstacle
// velocity in each direction
Vector3 STLObstacle::vitesseMaxPerDirection() const
{
  Vector3 vmax;

  // TO DO

  return ( vmax );

}




// ----------------------------------------------------------------------------
// Returns the number of points to write the STL obstacle in a
// Paraview format
int STLObstacle::numberOfPoints_PARAVIEW() const
{
  return ( int( m_allSTLVertices.size() ) );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the
// STL obstacle shape in a Paraview format
int STLObstacle::numberOfCells_PARAVIEW() const
{
  return ( int( m_allSTLTriangles.size() ) );
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the STL obstacle in a Paraview format
list<Point3> STLObstacle::get_polygonsPts_PARAVIEW(
	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  list<STLVertex*>::const_iterator il;
  Point3 p, pp;
  Transform const* transform = m_geoRBWC->getTransform();
  
  for (il=m_allSTLVertices.begin();il!=m_allSTLVertices.end();il++)
  {
    p = (*il)->m_p;    
    pp = (*transform)( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );        
  }

  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the points describing the STL obstacle in a Paraview format
void STLObstacle::write_polygonsPts_PARAVIEW( ostream &f,
  	Vector3 const* translation ) const
{
  list<STLVertex*>::const_iterator il;
  Point3 p, pp;
  Transform const* transform = m_geoRBWC->getTransform();
  
  for (il=m_allSTLVertices.begin();il!=m_allSTLVertices.end();il++)
  {
    p = (*il)->m_p;    
    pp = (*transform)( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;      
  }
}




// ----------------------------------------------------------------------------
// Writes the STL obstacle in a Paraview format
void STLObstacle::write_polygonsStr_PARAVIEW( list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  vector<STLTriangle>::const_iterator itt;
  
  for (itt=m_allSTLTriangles.begin();itt!=m_allSTLTriangles.end();itt++)
  {  
    connectivity.push_back( firstpoint_globalnumber + 
    	int(get<0>(itt->m_v)->m_id) );
    connectivity.push_back( firstpoint_globalnumber + 
    	int(get<1>(itt->m_v)->m_id) );
    connectivity.push_back( firstpoint_globalnumber + 
    	int(get<2>(itt->m_v)->m_id) );
    last_offset += 3;
    offsets.push_back( last_offset );
    cellstype.push_back( 5 );  
  }
  
  firstpoint_globalnumber += int( m_allSTLVertices.size() );  
}




// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
void STLObstacle::writePositionInFluid( ostream &fileOut )
{
  // TO DO
}




// ----------------------------------------------------------------------------
// Searches a vertex by its id number, returns a pointer to the
// vertex if found, returns NULL otherwise
STLVertex* STLObstacle::VidSearch(size_t ide)
{
  list<STLVertex*>::iterator itv;
  for (itv = m_allSTLVertices.begin(); itv != m_allSTLVertices.end(); ++itv)
    if ( (**itv).m_id == ide )
      return ( *itv );

  return ( NULL );
}




// ---------------------------------------------------------------------------
// Is a point already a vertex in the list of all vertices ?
int STLObstacle::isPInV(double x, double y, double z)
{	
  int id = 0;
  list<STLVertex*>::iterator it;

  for (it = m_allSTLVertices.begin(); it != m_allSTLVertices.end(); ++it)  
  {
    if (( fabs(x-(*it)->m_p[0]) < EPSILON ) 
    	&& ( fabs(y-(*it)->m_p[1]) < EPSILON ) 
	&& ( fabs(z-(*it)->m_p[2]) < EPSILON ))
      return ( id );
    id++;
  }

  return ( - 1 );
}




// ---------------------------------------------------------------------------
// Displays all triangles
void STLObstacle::displayAllSTLTriangles()
{
  ofstream myfile;
  myfile.open ("triangles.stl");
  
  vector<STLTriangle>::iterator itt;

  myfile << "solid dataset" << endl;
  for (itt = m_allSTLTriangles.begin(); itt != m_allSTLTriangles.end(); ++itt)  
  {
    size_t id1  = get<0>(itt->m_v)->m_id;
    size_t id2  = get<1>(itt->m_v)->m_id;
    size_t id3  = get<2>(itt->m_v)->m_id;

    STLVertex *v1 = VidSearch(id1);
    STLVertex *v2 = VidSearch(id2);
    STLVertex *v3 = VidSearch(id3);

    myfile << "   facet normal " << itt->m_n[0] << " " << itt->m_n[1] 
    	<< " " << itt->m_n[2] << endl;
    myfile << "        outer loop " << endl;
    myfile << "            vertex " << v1->m_p[0] << " " << v1->m_p[1] 
    	<< " " << v1->m_p[2] << endl;
    myfile << "            vertex " << v2->m_p[0] << " " << v2->m_p[1] 
    	<< " " << v2->m_p[2] << endl;
    myfile << "            vertex " << v3->m_p[0] << " " << v3->m_p[1] 
    	<< " " << v3->m_p[2] << endl;
    myfile << "        endloop" << endl;
    myfile << "   endfacet"	<< endl;
  }
  
  myfile << "endsolid" << endl;
  myfile.close();
}  




// ---------------------------------------------------------------------------
// Displays all vertices
void STLObstacle::displayAllSTLVertices()
{
  ofstream myfile;
  myfile.open ("vertices.obj");
  
  list<STLVertex*>::iterator itv;

  myfile << "# List of geometric vertices" << endl;
  for (itv = m_allSTLVertices.begin(); itv != m_allSTLVertices.end(); ++itv)
    myfile << "v " << (**itv).m_p << endl;
  
  myfile << "# List of vertex normals in (x,y,z) form" << endl;
  for (itv = m_allSTLVertices.begin(); itv != m_allSTLVertices.end(); ++itv)  
    myfile << "vn " << (**itv).m_n << endl;
   
  myfile.close();
}



//---------------------------------------------------------------------------
// Reads the STL from a file
void STLObstacle::readSTL( string const& filename )
{

  double xn, yn, zn;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3;
  int kk = 0;
  size_t vid = 0, tid = 0;
  string notsodummy,dummy;

  // Open the STL file
  ifstream inFilep( filename );
  cout << GrainsExec::m_shift9 << "Reading STL file: " << filename << endl;  

  // Check if the file is ASCII or binary
  int c;
  while ( ( c = inFilep.get() ) != EOF && c <= 127 );
  if ( c == EOF ) cout << GrainsExec::m_shift9 << "File is in ASCII format" 
  	<< endl;
  else cout << GrainsExec::m_shift9 << "File is in Binary format" << endl;
  inFilep.close();

  // Reading triangulation
  if ( c == EOF ) // ASCII
  {
    ifstream inFile( filename );
    string line;
    getline( inFile, line );

    while( getline( inFile, line ) )
    {
      istringstream iss1( line );
      iss1 >> notsodummy;
      if ( notsodummy.compare( "endsolid" ) == 0 ) break;      
      iss1 >> dummy;
      iss1 >> xn;
      iss1 >> yn;
      iss1 >> zn;
      
      getline( inFile, line );
      getline( inFile, line );
      
      istringstream iss2( line );
      iss2 >> dummy;
      iss2 >> x1;
      iss2 >> y1;
      iss2 >> z1;
      
      getline( inFile, line );
      istringstream iss3( line );
      iss3 >> dummy;
      iss3 >> x2;
      iss3 >> y2;
      iss3 >> z2;
	
      getline( inFile, line );
      istringstream iss4( line );
      iss4 >> dummy;
      iss4 >> x3;
      iss4 >> y3;
      iss4 >> z3;
       
      getline( inFile, line );
      getline( inFile, line );

      // Vertices
      m_llvls.push_back( make_tuple( x1, y1, z1 ) );
      m_llvls.push_back( make_tuple( x2, y2, z2 ) );
      m_llvls.push_back( make_tuple( x3, y3, z3 ) );

      // Normals
      m_llvns.push_back( make_tuple( xn, yn, zn ) );

      kk=kk+3;
    }

    inFile.close();
    m_npls=kk;
  } // end ASCII
  else // binary
  {
    ifstream inFileb( filename, ifstream::binary );

    // rdbuf returns a streambuf object associated with the
    // input fstream object ifs.
    filebuf* pbuf = inFileb.rdbuf();

    // Calculate the file's size.
    auto size = pbuf->pubseekoff( 0, inFileb.end );

    // Set the position pointer to the beginning of the file.
    pbuf->pubseekpos( 0 );

    // Allocate memory to contain file data.
    char* buffer = new char[(size_t)size];

    // Get file data. sgetn grabs all the characters from the streambuf
    // object 'pbuf'. The return value of sgetn is the number of characters
    // obtained - ordinarily, this value should be checked for equality
    // against the number of characters requested.
    pbuf->sgetn(buffer, size);

    char * bufptr = buffer;

    bufptr += 80;  // Skip past the header.
    bufptr += 4;   // Skip past the number of triangles.

    while (bufptr < buffer + size)
    {
      xn = *(float *)(bufptr);
      yn = *(float *)(bufptr + 4);
      zn = *(float *)(bufptr + 8);
      bufptr += 12;

      m_llvns.push_back(make_tuple(xn, yn, zn));

      x1 = *(float *)(bufptr);
      y1 = *(float *)(bufptr + 4);
      z1 = *(float *)(bufptr + 8);
      bufptr += 12;

      x2 = *(float *)(bufptr);
      y2 = *(float *)(bufptr + 4);
      z2 = *(float *)(bufptr + 8);
      bufptr += 12;

      x3 = *(float *)(bufptr);
      y3 = *(float *)(bufptr + 4);
      z3 = *(float *)(bufptr + 8);
      bufptr += 12;

      m_llvls.push_back(make_tuple(x1, y1, z1));
      m_llvls.push_back(make_tuple(x2, y2, z2));
      m_llvls.push_back(make_tuple(x3, y3, z3));

      // Filling m_allSTLVertices 
      Vector3 n;
      n[0] = xn; n[1] = yn; n[2] = zn;

      STLVertex *v1 = new STLVertex(x1,y1,z1,n,vid); 
      STLVertex *v2 = new STLVertex(x2,y2,z2,n,vid);  
      STLVertex *v3 = new STLVertex(x3,y3,z3,n,vid);

      if ( vid == 0 )
      {
        v1->m_id = 0;  
        v2->m_id = 1; 
        v3->m_id = 2; 
        vid += 3;

        m_allSTLVertices.push_back(v1);
        m_allSTLVertices.push_back(v2); 
        m_allSTLVertices.push_back(v3);  
      }
      else 
      {
        if (isPInV(x1, y1, z1)<0)
        {
          v1->m_id = vid; vid+=1;  
          m_allSTLVertices.push_back(v1);
        }
        else	
          v1->m_id = isPInV(x1, y1, z1);
	
        if (isPInV(x2, y2, z2)<0)
        {
          v2->m_id = vid; vid+=1; 
          m_allSTLVertices.push_back(v2);
        }
        else	
          v2->m_id = isPInV(x2, y2, z2);

        if (isPInV(x3, y3, z3)<0)
        {
          v3->m_id = vid; vid+=1;  
          m_allSTLVertices.push_back(v3);
        }
        else	
          v3->m_id = isPInV(x3, y3, z3);	   
      }
	
      // Filling m_allSTLTriangles
      tuple<STLVertex*,STLVertex*,STLVertex*> ve;
      ve = make_tuple(v1, v2, v3);
      STLTriangle tr(ve,n,tid); tid+=1;
      m_allSTLTriangles.push_back(tr);	

      kk = kk+3;
      bufptr += 2;
    }

    inFileb.close();
    m_npls=kk;
    delete[] buffer;
  }
  // end reading triangulation

  // Compute normals at vertices
  STLnormalsAtVertices();   
  
  // Outputs
  displayAllSTLVertices();
  displayAllSTLTriangles();

  cout << GrainsExec::m_shift9 << "Number of unique vertices: " << 
  	m_allSTLVertices.size() << "/" << vid << endl;
  cout << GrainsExec::m_shift9 << "Delaunay triangles (STL): " << m_npls / 3 
  	<< endl;
  cout << GrainsExec::m_shift9 << "Construction of STL object completed" 
  	<< endl;
}




//---------------------------------------------------------------------------
// Sets the STL normal vectors at all vertices
void STLObstacle::STLnormalsAtVertices()
{
  vector<STLTriangle>::iterator itt;
  double* total_area = new double[m_allSTLVertices.size()];

  for (itt = m_allSTLTriangles.begin(); itt != m_allSTLTriangles.end(); ++itt)
  {
    STLVertex *v1 = get<0>(itt->m_v);
    STLVertex *v2 = get<1>(itt->m_v);
    STLVertex *v3 = get<2>(itt->m_v);

    v1->m_n[0] = 0.0;
    v1->m_n[1] = 0.0;
    v1->m_n[2] = 0.0;

    v2->m_n[0] = 0.0;
    v2->m_n[1] = 0.0;
    v2->m_n[2] = 0.0;

    v3->m_n[0] = 0.0;
    v3->m_n[1] = 0.0;
    v3->m_n[2] = 0.0;
  }

  for (size_t i = 0; i < m_allSTLVertices.size(); i++) total_area[i] = 0.0;

  for (itt = m_allSTLTriangles.begin(); itt != m_allSTLTriangles.end(); ++itt)
  {
    size_t id1  = get<0>(itt->m_v)->m_id;
    size_t id2  = get<1>(itt->m_v)->m_id;
    size_t id3  = get<2>(itt->m_v)->m_id;

    STLVertex *v1 = VidSearch(id1);
    STLVertex *v2 = VidSearch(id2);
    STLVertex *v3 = VidSearch(id3);

    v1->m_n[0] += itt->getSurfaceArea() * itt->m_n[0];
    v1->m_n[1] += itt->getSurfaceArea() * itt->m_n[1];
    v1->m_n[2] += itt->getSurfaceArea() * itt->m_n[2];

    v2->m_n[0] += itt->getSurfaceArea() * itt->m_n[0];
    v2->m_n[1] += itt->getSurfaceArea() * itt->m_n[1];
    v2->m_n[2] += itt->getSurfaceArea() * itt->m_n[2];

    v3->m_n[0] += itt->getSurfaceArea() * itt->m_n[0];
    v3->m_n[1] += itt->getSurfaceArea() * itt->m_n[1];
    v3->m_n[2] += itt->getSurfaceArea() * itt->m_n[2];

    total_area[v1->m_id] += itt->getSurfaceArea();
    total_area[v2->m_id] += itt->getSurfaceArea();
    total_area[v3->m_id] += itt->getSurfaceArea();
  }

  list<STLVertex*>::iterator itv;

  for (itv = m_allSTLVertices.begin(); itv != m_allSTLVertices.end(); ++itv)  
  {
    (**itv).m_n[0] = (**itv).m_n[0] / total_area[(**itv).m_id];
    (**itv).m_n[1] = (**itv).m_n[1] / total_area[(**itv).m_id];
    (**itv).m_n[2] = (**itv).m_n[2] / total_area[(**itv).m_id];
  }
  
  delete [] total_area;
}





//---------------------------------------------------------------------------
// Returns the dot product of two vectors
double STLObstacle::dotProduct( double vect_A[], double vect_B[] )
{
  double product = 0;

  for (int i = 0; i < 3; i++)
     product = product + vect_A[i] * vect_B[i];

  return ( product );
}




//---------------------------------------------------------------------------
// Computes the cross product of two vectors
void STLObstacle::crossProduct( double vect_A[], double vect_B[], 
	double cross_P[] )
{
  cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
  cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
  cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}




//---------------------------------------------------------------------------
// Computes the difference between two vectors
void STLObstacle::diffProduct( double vect_A[], double vect_B[], 
	double diff_P[] )
{
  diff_P[0] = vect_A[0] - vect_B[0];
  diff_P[1] = vect_A[1] - vect_B[1];
  diff_P[2] = vect_A[2] - vect_B[2];
}




//---------------------------------------------------------------------------
// Computes the sign of the signed volume of the tetrahedron 
// (A,B,C,D): 1 -> positive or 0 -> negative
int STLObstacle::orient3d( double vect_A[], double vect_B[], 
	double vect_C[], double vect_D[] )
{
  double tmp1[3], tmp2[3], tmp3[3], tmp4[3];

  diffProduct( vect_B, vect_A, tmp1 );
  diffProduct( vect_C, vect_A, tmp2 );
  diffProduct( vect_D, vect_A, tmp3 );
  crossProduct( tmp1, tmp2, tmp4 );

  if ( dotProduct(tmp4,tmp3) > 0. ) return ( 0 );
  else return ( 1 );
}




//---------------------------------------------------------------------------
// Returns 1 if a segment [q1 q2] intersects a triangle (tri1, tri2, tri3)
int STLObstacle::intersect3d( double q1[], double q2[], double tri1[], 
	double tri2[], double tri3[] )
{
  double s1 = orient3d( q1, tri1, tri2, tri3 );
  double s2 = orient3d( q2, tri1, tri2, tri3 );
  
  // Test whether the two extermities of the segment
  // are on the same side of the supporting plane of
  // the triangle
  if ( s1 == s2 ) return ( 0 );

  // Now we know that the segment 'straddles' the supporing
  // plane. We need to test whether the three tetrahedra formed
  // by the segment and the three edges of the triangle have
  // the same orientation
  int s3 = orient3d( q1, q2, tri1, tri2 );
  int s4 = orient3d( q1, q2, tri2, tri3 );
  int s5 = orient3d( q1, q2, tri3, tri1 );

  return ( s3 == s4 && s4 == s5 );
}



//---------------------------------------------------------------------------
// ???
int STLObstacle::intersect(Point3 P, Point3 P1, Point3 P2, Point3 P3)
{
  Vector3 *u = new Vector3;
  Vector3 *v = new Vector3;
  Vector3 *w = new Vector3;
  Vector3 *n = new Vector3;
  Vector3 *tmp1 = new Vector3;
  Vector3 *tmp2 = new Vector3;

  double alpha, beta, gamma;

  (*u)[X] = P2[X] - P1[X];
  (*u)[Y] = P2[Y] - P1[Y];
  (*u)[Z] = P2[Z] - P1[Z];

  (*v)[X] = P3[X] - P1[X];
  (*v)[Y] = P3[Y] - P1[Y];
  (*v)[Z] = P3[Z] - P1[Z];

  (*w)[X] = P[X] - P1[X];
  (*w)[Y] = P[Y] - P1[Y];
  (*w)[Z] = P[Z] - P1[Z];

  *n = *u ^ *v;

  *tmp1 = *u ^ *w; 
  gamma = ( *tmp1 * *n ) / Norm2(*n);
  
  *tmp2 = *w ^ *v; 
  beta = ( *tmp2 * *n) / Norm2(*n);

  alpha = 1. - gamma - beta;

  cout << "alpha beta gamma:" << alpha << " " << beta << " " << gamma << endl;

  if ( alpha >= 0. && alpha <= 1. && beta >= 0. && beta <= 1. 
  	&& gamma >= 0. && gamma <= 1. ) return ( 1 );
  else return ( 0 );

}




/*
//---------------------------------------------------------------------------
int STLObstacle::isVinT(STLVertex &vx, STLTriangle &T)
//---------------------------------------------------------------------------
{
   STLVertex v1 = get<0>(T.v);
   STLVertex v2 = get<1>(T.v);
   STLVertex v3 = get<2>(T.v);

   if (isVinV(v1,vx) || isVinV(v2,vx) || isVinV(v3,vx))
      return 1;
   else return 0;

}

//---------------------------------------------------------------------------
int STLObstacle::isVinV(STLVertex &v1, STLVertex &v2)
//---------------------------------------------------------------------------
{
   double eps = 1.0e-12;	
   
   if (( fabs(v1.p.m_comp[0]-v2.p[0]) < eps ) && 
       ( fabs(v1.p.m_comp[1]-v2.p[1]) < eps ) && 
       ( fabs(v1.p.m_comp[2]-v2.p[2]) < eps )     )
    
       return 1;
    
   else return 0;
}*/
