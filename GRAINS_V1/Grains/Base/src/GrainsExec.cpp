#include "GrainsExec.hh"
#include "GrainsMPIWrapper.hh"
#include <sys/types.h>
#include <unistd.h>


bool GrainsExec::m_MPI = false;
string GrainsExec::m_TIScheme = "SecondOrderLeapFrog";
bool GrainsExec::m_SphereAsPolyParaview = false;
int GrainsExec::m_MPI_verbose = 0;
bool GrainsExec::m_isReloaded = false;
string GrainsExec::m_ReloadType = "new" ;
Vector3 GrainsExec::m_vgravity = Vector3Nul;
Vector3* GrainsExec::m_translationParaviewPostProcessing = NULL ;
bool GrainsExec::m_periodic = false;
bool GrainsExec::m_isGrainsCompFeatures = false;
bool GrainsExec::m_isGrainsPorosity = false;
string GrainsExec::m_ReloadDirectory = "";
string GrainsExec::m_SaveDirectory = "";
set<string> GrainsExec::m_additionalDataFiles;
bool GrainsExec::m_writingModeHybrid = false;
string GrainsExec::m_GRAINS_HOME = ".";
string GrainsExec::m_reloadFile_suffix = "B";
bool GrainsExec::m_exception_Contact = false;
bool GrainsExec::m_exception_Displacement = false;
bool GrainsExec::m_exception_Simulation = false;
string GrainsExec::m_shift0 = "";
string GrainsExec::m_shift1 = " ";
string GrainsExec::m_shift2 = "  ";
string GrainsExec::m_shift3 = "   ";
string GrainsExec::m_shift6 = "      ";
string GrainsExec::m_shift9 = "         ";
string GrainsExec::m_shift12 = "            ";
string GrainsExec::m_shift15 = "               ";
bool GrainsExec::m_output_data_at_this_time = false;
GrainsMPIWrapper* GrainsExec::m_wrapper = NULL;
list<App*> GrainsExec::m_allApp;
size_t GrainsExec::m_total_nb_physical_particles = 0;
list< pair<Point3*,VertexBase *> > GrainsExec::m_allPolytopeRefPointBase;
list<IndexArray*> GrainsExec::m_allPolytopeNodeNeighbors;
list<IndexArray*> GrainsExec::m_allPolytopeNodesIndex;
list<vector< vector<int> >*> GrainsExec::m_allPolyhedronFacesConnectivity;
string GrainsExec::m_inputFile;
int GrainsExec::m_return_syscmd = 0;
bool GrainsExec::m_preCollision_cyl = false;




// ----------------------------------------------------------------------------
// Default constructor
GrainsExec::GrainsExec()
{}




// ----------------------------------------------------------------------------
// Destructor
GrainsExec::~GrainsExec()
{}




// ----------------------------------------------------------------------------
// Frees the content of the garbage collector
void GrainsExec::GarbageCollector()
{
  if ( m_translationParaviewPostProcessing )
    delete m_translationParaviewPostProcessing;

  if ( !m_allPolytopeRefPointBase.empty() )
  {
    for (list< pair<Point3*,VertexBase *> >::iterator
    	ilPointBase=m_allPolytopeRefPointBase.begin();
  	ilPointBase!=m_allPolytopeRefPointBase.end(); ilPointBase++)
    {
      delete [] ilPointBase->first;
      delete ilPointBase->second;
    }
    m_allPolytopeRefPointBase.clear();
  }

  if ( !m_allPolytopeNodeNeighbors.empty() )
  {
    for (list<IndexArray*>::iterator ilIA=m_allPolytopeNodeNeighbors.begin();
    	ilIA!=m_allPolytopeNodeNeighbors.end(); ilIA++)
      // les pointeurs pointent sur des tableaux d'IndexArray, d'ou
      // l'utilisation de "delete []" au lieu de "delete"
      delete [] *ilIA;
    m_allPolytopeNodeNeighbors.clear();
  }

  if ( !m_allPolytopeNodesIndex.empty() )
  {
    for (list<IndexArray*>::iterator ilNI=m_allPolytopeNodesIndex.begin();
    	ilNI!=m_allPolytopeNodesIndex.end(); ilNI++)
      delete *ilNI;
    m_allPolytopeNodeNeighbors.clear();
  }

  if ( !m_allPolyhedronFacesConnectivity.empty() )
  {
    for (list<vector< vector<int> >*>::iterator
    	ilFC=m_allPolyhedronFacesConnectivity.begin();
    	ilFC!=m_allPolyhedronFacesConnectivity.end(); ilFC++)
      delete *ilFC;
    m_allPolytopeNodeNeighbors.clear();
  }
}




// ----------------------------------------------------------------------------
// Returns a pointer to the MPI wrapper
GrainsMPIWrapper* GrainsExec::getComm()
{
  return ( m_wrapper );
}




// ----------------------------------------------------------------------------
// Sets the MPI wrapper
void GrainsExec::setComm( GrainsMPIWrapper* wrapper_ )
{
  m_wrapper = wrapper_;
}




// ----------------------------------------------------------------------------
// Returns the list of applications
list<App*> GrainsExec::get_listApp()
{
  return ( m_allApp );
}




// ----------------------------------------------------------------------------
// Returns the total number of particles in the physical system 
// (i.e. on all subdomains/processes), i.e. sum of total number of active 
// particles with tag 0 or 1 and inactive particles
size_t GrainsExec::getTotalNumberPhysicalParticles()
{
  return ( m_total_nb_physical_particles );
}




// ----------------------------------------------------------------------------
// Sets the total number of particles in the physical system 
// (i.e. on all subdomains/processes), i.e. sum of total number of active 
// particles with tag 0 or 1 and inactive particles
void GrainsExec::setTotalNumberPhysicalParticles( size_t const& nb_ )
{
  m_total_nb_physical_particles = nb_;
}




// ----------------------------------------------------------------------------
// Sets the list of applications
void GrainsExec::set_listApp( list<App*> allApp_ )
{
  m_allApp = allApp_;
}




// ----------------------------------------------------------------------------
// Writes a float number with a given number of digits
string GrainsExec::doubleToString( double const& figure, int const& size )
{
  ostringstream oss;
  oss.width(size);
  oss << left << figure;

  return ( oss.str() );
}




// ----------------------------------------------------------------------------
// Writes a float number with a prescribed format and a prescribed
// number of digits after the decimal point
string GrainsExec::doubleToString( ios_base::fmtflags format, int digits,
      	double const& number )
{
  ostringstream oss;
  if ( number != 0. )
  {
    oss.setf( format, ios::floatfield );
    oss.precision( digits );
  }
  oss << number;

  return ( oss.str() );
}




// ----------------------------------------------------------------------------
// Writes an integer in a string
string GrainsExec::intToString( int const& figure )
{
  ostringstream oss;
  oss << figure;

  return ( oss.str() );
}




// ----------------------------------------------------------------------------
// Returns memory used by this process
size_t GrainsExec::used_memory( void )
{
  ostringstream os ;
  string word ;
  size_t result = 0 ;

  os << "/proc/" << getpid() << "/status" ;

  ifstream in( os.str().c_str() ) ;
  if( !in )
    cout << "GrainsExec::used_memory : unable to open " << os.str() << endl ;
  else
  {
    while( !in.eof() )
    {
      in >> word ;
      if( word == "VmSize:" )
      {
        in >> result ;
        in >> word ;
        if( !( word == "kB" ) )
          cout << "GrainsExec::used_memory : Unit is " << word << endl ;
        result *= 1000 ;
        break ;
      }
    }
    in.close() ;
  }

  return ( result ) ;
}




// ----------------------------------------------------------------------------
// Writes memory used by this process in a stream
void GrainsExec::display_memory( ostream& os, size_t memory )
{
  static size_t const mo = 1024*1024 ;
  static size_t const go = 1024*1024*1024 ;

  if( memory > go )
    os << ( (double) memory )/go << " Go" << std::flush;
  else if( memory > mo )
    os << ( (double) memory )/mo << " Mo" << std::flush ;
  else os << memory << " octets" << std::flush ;
}




// ----------------------------------------------------------------------------
// Adds a list of reference points of a polytope
void GrainsExec::addOnePolytopeRefPointBase( Point3* refPB, VertexBase* refVB )
{
  pair<Point3*,VertexBase *> pp( refPB, refVB );
  m_allPolytopeRefPointBase.push_back( pp );
}




// ----------------------------------------------------------------------------
// Adds a description of vertices of a polytope
void GrainsExec::addOnePolytopeNodeNeighbors( IndexArray* idar )
{
  m_allPolytopeNodeNeighbors.push_back( idar );
}




// ----------------------------------------------------------------------------
// Adds an array of indices of vertex of a polytope
void GrainsExec::addOnePolytopeNodeIndex( IndexArray* idar )
{
  m_allPolytopeNodesIndex.push_back( idar );
}




// ----------------------------------------------------------------------------
// Adds a face connectivity of a polyhedron
void GrainsExec::addOnePolyhedronFaceConnectivity(
	vector< vector<int> >* faceCon )
{
  m_allPolyhedronFacesConnectivity.push_back( faceCon );
}




// ----------------------------------------------------------------------------
// Checks the last output time in a file and deletes all subsequent
// lines, i.e., each line whose time is after the current time
void GrainsExec::checkTime_outputFile( string const& filename,
	const double& current_time )
{
  string tline, syscom;
  istringstream iss;
  double last_output_time, tt;

  // Get last output time
  ifstream fileIN( filename.c_str(), ios::in );
  if ( fileIN.is_open() )
  {
    getline( fileIN, tline );
    while ( !fileIN.eof() )
    {
      iss.str( tline );
      iss >> last_output_time;
      iss.clear();
      getline( fileIN, tline );
    }
    fileIN.close();

    // If last output time is greater than current time, delete all data
    // from last output time to current time
    if ( last_output_time - current_time > 1.e-12 )
    {
      cout << "Time inconsistency in output file "
    	<< filename << endl;
      ifstream fileIN_( filename.c_str(), ios::in );
      ofstream fileOUT( (filename+".tmp").c_str(), ios::out );
      getline( fileIN_, tline );
      while ( !fileIN_.eof() )
      {
        iss.str( tline );
        iss >> tt;
        iss.clear();
        if ( current_time - tt > -1.e-10 )
          fileOUT << tline << endl;
        getline( fileIN_, tline );
      }
      fileIN_.close();
      fileOUT.close();

      syscom = "mv " + filename + ".tmp" + " " + filename ;
      m_return_syscmd = system(syscom.c_str());
    }
  }
}




// ----------------------------------------------------------------------------
// Returns the root of a complete file name. Example: if input is
// Titi/tutu/toto, returns Titi/tutu
string GrainsExec::extractRoot( string const& FileName )
{
  size_t pos = FileName.rfind( "/" );
  string root ;
  if ( pos != string::npos )
  {
    root = FileName ;
    root.erase( root.begin() + pos, root.end() );
  }
  else root = "." ;

  return ( root );
}




// ----------------------------------------------------------------------------
// Returns the file name of a complete file name. Example: if input is
// Titi/tutu/toto, returns toto
string GrainsExec::extractFileName( string const& FileName )
{
  size_t pos = FileName.rfind( "/" );
  string result = FileName;
  if ( pos != string::npos )
    result.erase( result.begin(), result.begin() + pos + 1 );

  return ( result );
}




// ----------------------------------------------------------------------------
// Checks that all reload files are in the same directory (primarily
// checks that files for polyhedrons et polygons are there)
void GrainsExec::checkAllFilesForReload()
{
  if ( !m_additionalDataFiles.empty() )
  {
    set<string>::iterator is;
    string fileName ;
    string cmd = "bash " + GrainsExec::m_GRAINS_HOME
    	+ "/Tools/ExecScripts/addFiles.exec";

    for (is=m_additionalDataFiles.begin();is!=m_additionalDataFiles.end();is++)
    {
      fileName = GrainsExec::extractFileName( *is );
      cmd += " " + *is + " " + fileName + " " + m_SaveDirectory;
    }

    m_return_syscmd = system( cmd.c_str() );
  }
}




// ----------------------------------------------------------------------------
// Returns the reload file name from the restart time table
string GrainsExec::restartFileName_AorB( string const& rootName,
  	string const& RFTable_ext )
{
  string buf, rfn;
  ifstream FILE_IN( ( rootName + RFTable_ext ).c_str(), ios::in );
  if ( FILE_IN.is_open() )
  {
    while( !FILE_IN.eof() )
      FILE_IN >> buf >> rfn ;
  }
  else rfn = rootName;
  FILE_IN.close() ;

  return ( rfn );
}




// ----------------------------------------------------------------------------
// Returns a random rotation matrix
Matrix GrainsExec::RandomRotationMatrix( size_t dim )
{
  Matrix rotation;

  double angleZ = 2. * PI * (double)rand() / RAND_MAX;
  Matrix rZ( cos(angleZ), -sin(angleZ), 0.,
  	sin(angleZ), cos(angleZ), 0.,
	0., 0., 1. );

  if ( dim == 3 )
  {
    double angleX = 2. * PI * (double)rand() / RAND_MAX;
    Matrix rX( 1., 0., 0.,
   	0., cos(angleX), -sin(angleX),
	0., sin(angleX), cos(angleX) );

    double angleY = 2. * PI * (double)rand() / RAND_MAX;
    Matrix rY( cos(angleY), 0., sin(angleY),
   	0., 1., 0.,
	-sin(angleY), 0., cos(angleY) );
    Matrix tmp = rY * rZ;
    rotation = rX * tmp;
  }
  else rotation = rZ;

  return ( rotation );
}




// ----------------------------------------------------------------------------
// Returns a random unit vector
Vector3 GrainsExec::RandomUnitVector( size_t dim )
{
  Vector3 vec;

  vec[X] = 2. * (double)rand() / RAND_MAX - 1.;
  vec[Y] = 2. * (double)rand() / RAND_MAX - 1.;
  if ( dim == 3 )
    vec[Z] = 2. * (double)rand() / RAND_MAX - 1.;
  else
    vec[Z] = 0.;

  vec.normalize();

  return( vec );
}



// ----------------------------------------------------------------------------
// Returns whether a sphere is fully in, fully out or intersects
// an axis-aligned cylinder. Returned values are 0, 1 and 2 respectively
size_t GrainsExec::AACylinderSphereIntersection( Point3 const& SphereCenter,
    	double const& SphereRadius,
	Point3 const& CylBottomCentre,
	double const& CylRadius,
	double const& CylHeight,
	size_t const& CylAxisDir,
	double const& tol )
{
  size_t inter = 1;
  Vector3 CtoC = SphereCenter - CylBottomCentre;
  double radial_distance = 0.;

  switch( CylAxisDir )
  {
    case 0:
      radial_distance = sqrt( CtoC[Y] * CtoC[Y] + CtoC[Z] * CtoC[Z] );
      if ( radial_distance <= CylRadius - SphereRadius + tol
      	&& SphereCenter[X] >= CylBottomCentre[X] + SphereRadius - tol
	&& SphereCenter[X] <= CylBottomCentre[X] + CylHeight - SphereRadius
		+ tol )
	inter = 0;
      else if ( radial_distance > CylRadius + SphereRadius
      	|| SphereCenter[X] < CylBottomCentre[X] - SphereRadius
	|| SphereCenter[X] > CylBottomCentre[X] + CylHeight + SphereRadius )
	inter = 1;
      else inter = 2;
      break;

    case 1:
      radial_distance = sqrt( CtoC[X] * CtoC[X] + CtoC[Z] * CtoC[Z] );
      if ( radial_distance <= CylRadius - SphereRadius + tol
      	&& SphereCenter[Y] >= CylBottomCentre[Y] + SphereRadius - tol
	&& SphereCenter[Y] <= CylBottomCentre[Y] + CylHeight - SphereRadius
		+ tol )
	inter = 0;
      else if ( radial_distance > CylRadius + SphereRadius
      	|| SphereCenter[Y] < CylBottomCentre[Y] - SphereRadius
	|| SphereCenter[Y] > CylBottomCentre[Y] + CylHeight + SphereRadius )
	inter = 1;
      else inter = 2;
      break;

    default: // i.e. 2
      radial_distance = sqrt( CtoC[X] * CtoC[X] + CtoC[Y] * CtoC[Y] );
      if ( radial_distance <= CylRadius - SphereRadius + tol
      	&& SphereCenter[Z] >= CylBottomCentre[Z] + SphereRadius - tol
	&& SphereCenter[Z] <= CylBottomCentre[Z] + CylHeight - SphereRadius
		+ tol )
	inter = 0;
      else if ( radial_distance > CylRadius + SphereRadius
      	|| SphereCenter[Z] < CylBottomCentre[Z] - SphereRadius
	|| SphereCenter[Z] > CylBottomCentre[Z] + CylHeight + SphereRadius )
	inter = 1;
      else inter = 2;
      break;
  }

  return ( inter ) ;
}




// ----------------------------------------------------------------------------
// Returns whether a point belongs to an axis-aligned cylinder
bool GrainsExec::isPointInAACylinder( Point3 const& pt,
	Point3 const& CylBottomCentre,
	double const& CylRadius,
	double const& CylHeight,
	size_t const& CylAxisDir,
	double const& tol )
{
  bool isIn = false;
  Vector3 PtoC = pt - CylBottomCentre;
  double radial_distance = 0.;

  switch( CylAxisDir )
  {
    case 0:
      radial_distance = sqrt( PtoC[Y] * PtoC[Y] + PtoC[Z] * PtoC[Z] );
      if ( radial_distance <= CylRadius + tol
      	&& pt[X] >= CylBottomCentre[X] - tol
	&& pt[X] <= CylBottomCentre[X] + CylHeight + tol )
	isIn = true;
      break;

    case 1:
      radial_distance = sqrt( PtoC[X] * PtoC[X] + PtoC[Z] * PtoC[Z] );
      if ( radial_distance <= CylRadius + tol
      	&& pt[Y] >= CylBottomCentre[Y] - tol
	&& pt[Y] <= CylBottomCentre[Y] + CylHeight + tol )
	isIn = true;
      break;

    default: // i.e. 2
      radial_distance = sqrt( PtoC[X] * PtoC[X] + PtoC[Y] * PtoC[Y] );
      if ( radial_distance <= CylRadius + tol
      	&& pt[Z] >= CylBottomCentre[Z] - tol
	&& pt[Z] <= CylBottomCentre[Z] + CylHeight + tol )
	isIn = true;
      break;
  }

  return ( isIn ) ;
}




// ----------------------------------------------------------------------------
// Computes and returns the 4 x 4 determinant of 4 points
double GrainsExec::PointDeterm4by4( Point3 const& p1, Point3 const& p2,
    	Point3 const& p3, Point3 const& p4 )
{
  double x2 = p2[X];
  double y2 = p2[Y];
  double z2 = p2[Z];
  double x3 = p3[X];
  double y3 = p3[Y];
  double z3 = p3[Z];
  double x4 = p4[X];
  double y4 = p4[Y];
  double z4 = p4[Z];

  double det =
  	p1[X] * ( y2*z3 + y3*z4 + y4*z2 - y4*z3 - y3*z2 - y2*z4 ) -
  	p1[Y] * ( x2*z3 + x3*z4 + x4*z2 - x4*z3 - x3*z2 - x2*z4 ) +
  	p1[Z] * ( x2*y3 + x3*y4 + x4*y2 - x4*y3 - x3*y2 - x2*y4 ) -
        ( x2*y3*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2 - x3*y2*z4 - x2*y4*z3 ) ;

  return ( det );
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside a tetrahedron
bool GrainsExec::isPointInTetrahedron( Point3 const& p1, Point3 const& p2,
    	Point3 const& p3, Point3 const& p4, Point3 const& p, bool check )
{
  double detTot = PointDeterm4by4( p1, p2, p3, p4 ),
	detOne = PointDeterm4by4( p, p2, p3, p4 ),
	detTwo = PointDeterm4by4( p1, p, p3, p4 ),
	detThree = PointDeterm4by4( p1, p2, p, p4 ),
	detFour = PointDeterm4by4( p1, p2, p3, p );
  bool isIn = false;
  double sumSubElem = detOne + detTwo + detThree + detFour;

  if ( check ) 
  {
    if ( fabs( detTot - sumSubElem ) > EPSILON  )
      cout << "ERROR: summation error in determinant 3D : " << 
      	detTot - sumSubElem << endl;

    if ( detTot == 0. ) 
    {
      cout << "Degenerated tetrahedron: det == 0 " << endl;
      abort();
    }
  }

  double r1 = detOne/detTot, r2 = detTwo/detTot, r3 = detThree/detTot,
	r4 = detFour/detTot ;

  // Allowed threshold for ratio of determinant to handle points on surface
  if ( ( r1 > - EPSILON )  && ( r2 > - EPSILON ) &&
       ( r3 > - EPSILON ) && ( r4 > - EPSILON ) )
      isIn = true;

  return ( isIn );
}
