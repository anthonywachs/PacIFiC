#include "Obstacle.hh"
#include "CompositeObstacle.hh"
#include "SimpleObstacle.hh"
#include "LinkedCell.hh"
#include "ObstacleImposedVelocity.hh"
#include "Quaternion.hh"
#include "App.hh"
#include "Memento.hh"
#include <math.h>


int Obstacle::m_totalNbSingleObstacles = 0;
bool Obstacle::m_MoveObstacle = true ;
bool Obstacle::m_isConfinement = false;
int Obstacle::m_minID = 0;


//-----------------------------------------------------------------------------
// Constructor with name and autonumbering as input parameters
Obstacle::Obstacle( string const& s, bool const& autonumbering ) 
  : Component()
  , m_name( s )
  , m_ismoving( false )
  , m_indicator( 0. )
  , m_ObstacleType ( "0" )
  , m_restrict_geommotion( false )
{
  if ( autonumbering )
  {
    Obstacle::m_minID--;
    m_id = Obstacle::m_minID;
  }
}




//-----------------------------------------------------------------------------
// Destructor
Obstacle::~Obstacle()
{}




// ----------------------------------------------------------------------------
// Links imposed kinematics to the obstacle and returns true if the
// linking process is successful
bool Obstacle::LinkImposedMotion( ObstacleImposedVelocity* imposed )
{
  bool status = false;
  if ( m_name == imposed->getObstacleName() )
  {
    m_kinematics.append( imposed );
    status = true;
  }

  return ( status );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the torsor exerted on the obstacle
Torsor const* Obstacle::getTorsor()
{
  return ( &m_torsor );
}




// ----------------------------------------------------------------------------
// Links imposed kinematics to the obstacle and returns true if the
// linking process is successful
bool Obstacle::LinkImposedMotion( ObstacleImposedForce* imposed )
{
  bool status = false;
  if ( m_name == imposed->getObstacleName() )
  {
    m_confinement.append( imposed );
    Obstacle::m_isConfinement = status = true;
  }
  return ( status );
}




// ----------------------------------------------------------------------------
// Composes the obstacle kinematics with another "higher level" velocity
// kinematics
void Obstacle::Compose( ObstacleKinematicsVelocity const& other,
	Vector3 const& lever )
{
  m_kinematics.Compose( other, lever );
}




// ----------------------------------------------------------------------------
// Composes the obstacle kinematics with another "higher level" force
// kinematics
void Obstacle::Compose( ObstacleKinematicsForce const& other )
{
  m_confinement.Compose( other );
}




// ----------------------------------------------------------------------------
// Returns total number of single obstacles
int Obstacle::getTotalNbSingleObstacles()
{
  return ( Obstacle::m_totalNbSingleObstacles );
}




// ----------------------------------------------------------------------------
// Returns the velocity at a point in space based on the
// translational and angular velocity of the obstacle. This method assumes
// that the point belongs to the obstacle but this assumption is not verified.
Vector3 Obstacle::getVelocityAtPoint( Point3 const& pt ) const
{
  Vector3 lever = pt - *m_geoRBWC->getCentre();

  return ( m_translationalVelocity + ( m_angularVelocity ^ lever ) );
}




// ----------------------------------------------------------------------------
// Returns the angular velocity
Vector3 const* Obstacle::getAngularVelocity() const
{
  return ( &m_angularVelocity );
}




// ----------------------------------------------------------------------------
// Returns the translational velocity
Vector3 const* Obstacle::getTranslationalVelocity() const
{
  return ( &m_translationalVelocity );
}




// ----------------------------------------------------------------------------
// Resets kinematics to 0
void Obstacle::resetKinematics()
{
  m_kinematics.reset();
  m_confinement.reset();
  m_translationalVelocity = 0.;
  m_angularVelocity = 0.;
  m_ismoving = false;  
}




// ----------------------------------------------------------------------------
// Sets imposed velocity kinematics
void Obstacle::setKinematics( ObstacleKinematicsVelocity& kine_ )
{
  m_kinematics.set( kine_ );
}




// ----------------------------------------------------------------------------
// Sets kinematics using translational and angular velocities
void Obstacle::setVelocity( Vector3 const* vtrans, Vector3 const* vrot )
{
  m_kinematics.setVelocity( vtrans, vrot );
}




// ----------------------------------------------------------------------------
// Writes the obstacle position
void Obstacle::writePosition( ostream& position )
{
  position << "*Obstacle\n" << m_name << '\n';
  Component::writePosition( position );
  position << "*EndObstacle\n\n";
}




// ----------------------------------------------------------------------------
// Writes the obstacle's "static" data
void Obstacle::writeStatic( ostream& fileOut ) const
{
  fileOut << "*Obstacle\n" << m_name << endl;
  Component::writeStatic( fileOut );
  fileOut << "*EndObstacle\n\n";
}




// ----------------------------------------------------------------------------
// Updates indicator for Paraview post-processing
void Obstacle::updateIndicator( double time, double dt ) 
{}




// ----------------------------------------------------------------------------
// Writes the identity of a component (often id number and address)
void Obstacle::writeIdentity( ostream& file ) const
{
  file << m_name;
}




// ----------------------------------------------------------------------------
// Returns indicator for Paraview post-processing
double Obstacle::getIndicator() const
{
  return ( m_indicator );
}




// ----------------------------------------------------------------------------
// Returns whether the obstacle has moved over the last time step
bool Obstacle::hasMoved() const
{
  return ( m_ismoving );
}




// ----------------------------------------------------------------------------
// Saves obstacle state
void Obstacle::saveState()
{
  if (!m_memento)
    m_memento = new ConfigurationMemento();
  m_memento->m_position = *m_geoRBWC->getTransform();
  m_kinematics.saveState();
}




// ----------------------------------------------------------------------------
// Restores obstacle state
void Obstacle::restoreState()
{
  m_geoRBWC->setTransform( m_memento->m_position );
  m_kinematics.restoreState();
}




// ----------------------------------------------------------------------------
// Sets the boolean to actually move obstacles
void Obstacle::setMoveObstacle( bool const& depObs )
{
  Obstacle::m_MoveObstacle = depObs ;
}




// ----------------------------------------------------------------------------
// Geometric motion of obstacles
bool Obstacle::getMoveObstacle()
{
  return ( Obstacle::m_MoveObstacle ) ;
}




// ----------------------------------------------------------------------------
// Returns obstacle name
string Obstacle::getName() const
{
  return ( m_name );
}




// ----------------------------------------------------------------------------
// Sets indicator for Paraview post-processing
void Obstacle::setIndicator( double const& value )
{
  m_indicator = value;
}




// ----------------------------------------------------------------------------
// Contact between an obstacle and a component. If contact exists,
// computes the contact force and torque and adds to each component
void Obstacle::InterAction( Component* voisin,
	double dt, double const& time, LinkedCell* LC )
{
  try
  {
    cout << "Warning when calling Obstacle::InterAction() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
    exit(10);
  }
  catch ( ContactError const& ) { throw ContactError(); }
}





// ----------------------------------------------------------------------------
// Returns the number of points to write the obstacle in a Paraview format
int Obstacle::numberOfPoints_PARAVIEW() const
{
  cout << "Warning when calling Obstacle::numberOfPoints_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the obstacle
// shape in a Paraview format
int Obstacle::numberOfCells_PARAVIEW() const
{
  cout << "Warning when calling Obstacle::numberOfCells_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the obstacle in a Paraview format
list<Point3> Obstacle::get_polygonsPts_PARAVIEW( Vector3 const* translation )
const
{
  cout << "Warning when calling Obstacle::get_polygonsPts_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Writes the points describing the obstacle in a Paraview format
void Obstacle::write_polygonsPts_PARAVIEW( ostream &f,
  	Vector3 const* translation ) const
{
  cout << "Warning when calling Obstacle::write_polygonsPts_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Writes the obstacle in a Paraview format
void Obstacle::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  cout << "Warning when calling Obstacle::write_polygonsStr_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
void Obstacle::writePositionInFluid( ostream &fileOut )
{
  cout << "Warning when calling Obstacle::writePositionInFluid() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Initializes all contact map entries to false
void Obstacle::setContactMapToFalse()
{
  cout << "Warning when calling Obstacle::setContactMapToFalse() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Set contact map cumulative features to zero */
void Obstacle::setContactMapCumulativeFeaturesToZero()
{
  cout << "Warning when calling Obstacle::setContactMapCumulativeFeaturesToZero() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Updates contact map
void Obstacle::updateContactMap()
{
  cout << "Warning when calling Obstacle::updateContactMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Does the contact exist in the map? If so, return true and make
// kdelta, prev_normal and cumulSpringTorque point to the memorized info. 
// Otherwise, return false and set those pointers to NULL.
bool Obstacle::getContactMemory( std::tuple<int,int,int> const& id,
  	Vector3* &kdelta, Vector3* &prev_normal, Vector3* &cumulSpringTorque,
  	bool createContact )
{
  cout << "Warning when calling Obstacle::getContactMemory() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}	




// ----------------------------------------------------------------------------
// Adds new contact in the map
void Obstacle::addNewContactInMap( std::tuple<int,int,int> const& id,
  	Vector3 const& kdelta, Vector3 const& prev_normal,
  	Vector3 const& cumulSpringTorque )
{
  cout << "Warning when calling Obstacle::addNewContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Increases cumulative tangential motion with component id
void Obstacle::addDeplContactInMap( std::tuple<int,int,int> const& id,
  	Vector3 const& kdelta, Vector3 const& prev_normal,
  	Vector3 const& cumulSpringTorque )
{
  cout << "Warning when calling Obstacle::addDeplContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Writes the contact map information in an array of doubles
void Obstacle::copyContactMap( double* destination, int start_index )
{
  cout << "Warning when calling Obstacle::copyContactMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Adds a single contact info to the contact map
void Obstacle::copyContactInMap( std::tuple<int,int,int> const& id,
  	bool const& isActive, Vector3 const& kdelta, Vector3 const& prev_normal,
  	Vector3 const& cumulSpringTorque )
{
  cout << "Warning when calling Obstacle::copyContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}	
	
	
	

// ----------------------------------------------------------------------------
// Returns the number of contacts in the contact map */
int Obstacle::getContactMapSize()
{
  cout << "Warning when calling Obstacle::getContactMapSize() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Displays the active neighbours in the format "my_elementary_id/neighbour_id/
// neightbout_elementary_id ; ...". Useful for debugging only.
void Obstacle::printActiveNeighbors( int const& id )
{
  cout << "Warning when calling Obstacle::printActiveNeighbors() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Returns the minimum ID number of an obstacle
int Obstacle::getMinIDnumber()
{
  return ( m_minID );
}




// ----------------------------------------------------------------------------
// Adds a force exerted at a point  to the torsor (torsor adds torque
// automatically)
void Obstacle::addForce( Point3 const& point, Vector3 const& force, 
	int tagSecondComp )
{
  // In case of a contact force, if the tag of the other component is 2, this
  // component is a periodic or parallel clone particle and its contribution 
  // must not be added as its periodic/parallel master particle's contribution 
  // is already accounted for 
  if ( tagSecondComp != 2 ) m_torsor.addForce( point, force );
}




// ----------------------------------------------------------------------------
// Adds a torque to the torsor
void Obstacle::addTorque( Vector3 const& torque, int tagSecondComp )
{
  // In case of a contact force, if the tag of the other component is 2, this
  // component is a periodic or parallel clone particle and its contribution 
  // must not be added as its periodic/parallel master particle's contribution 
  // is already accounted for 
  if ( tagSecondComp != 2 ) m_torsor.addTorque( torque );
}




// ----------------------------------------------------------------------------
// Checks if there is anything special to do about periodicity and
// if there is applies periodicity
void Obstacle::periodicity( LinkedCell* LC )
{}




// ----------------------------------------------------------------------------
// Restricts the geometric directions of translational motion 
void Obstacle::setRestrictedGeomDirMotion( list<size_t> const& dir )
{
  m_restrict_geommotion = true;
  m_dir_restricted_geommotion = dir;
}




// ----------------------------------------------------------------------------
// Sets the obstacle name
void Obstacle::setName( string const& name_ )
{
  m_name = name_ ;
}




// ----------------------------------------------------------------------------
// Sets obstacle translational and angular velocities from active
// kinematics with imposed velocity and  with imposed force
void Obstacle::setVelocity()
{
  m_translationalVelocity = *(m_kinematics.getTranslationalVelocity())
  	+ *(m_confinement.getTranslationalVelocity());
  m_angularVelocity = *(m_kinematics.getAngularVelocity());
} 




// ----------------------------------------------------------------------------
// Returns whether to store the contact force for post-processing 
bool Obstacle::storePPForce( Component const* othercomp ) const
{
  // Here we know that the other component is automatically a particle
  return ( othercomp->getTag() != 2 );
}
