#include "GrainsMPIWrapper.hh"
#include "SimpleObstacle.hh"
#include "LinkedCell.hh"
#include "Grains.hh"
#include "ContactBuilderFactory.hh"
#include "Particle.hh"
#include "PointContact.hh"
#include "Cell.hh"
#include "GrainsExec.hh"
#include "Memento.hh"
#include "ContactForceModel.hh"
#include <sstream>
#include <limits>
#include <assert.h>
using namespace std;


// ----------------------------------------------------------------------------
// Constructor with name as input parameter
SimpleObstacle::SimpleObstacle( const string &s )
  : Obstacle( s, true )
  , m_transferToFluid( false )
{
  m_ObstacleType = "SimpleObstacle";
  Obstacle::m_totalNbSingleObstacles++;
  m_LinkUpdate_frequency = 1;
  m_LinkUpdate_counter = 0;
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
SimpleObstacle::SimpleObstacle( DOMNode *root )
  : Obstacle( "", true )
  , m_transferToFluid( false )
{
  m_ObstacleType = "SimpleObstacle";

  Obstacle::m_totalNbSingleObstacles++;

  m_name = ReaderXML::getNodeAttr_String( root, "name" );

  // Convex - Position & Orientation
  m_geoRBWC = new RigidBodyWithCrust( root );

  // Materiau
  DOMNode* materiau_ = ReaderXML::getNode( root, "Material" );
  m_materialName = ReaderXML::getNodeValue_String( materiau_ );
  ContactBuilderFactory::defineMaterial( m_materialName, true );

  // Obstacle to transfer to the fluid
  DOMNode* status = ReaderXML::getNode( root, "Status" );
  if ( status )
    m_transferToFluid = ReaderXML::getNodeAttr_Int( status, "ToFluid" );

  m_obstacleBox = Component::BoundingBox();
  m_LinkUpdate_frequency = 1;
  m_LinkUpdate_counter = 0;
}




// ----------------------------------------------------------------------------
// Destructor
SimpleObstacle::~SimpleObstacle()
{
  m_inCells.clear();
}




// ----------------------------------------------------------------------------
// Adds an obstacle (single or composite) to the composite obstacle
// tree. Should not be called by SimpleObstacle.
void SimpleObstacle::append( Obstacle* obstacle )
{
  cout << "Warning when calling SimpleObstacle::append(Obstacle*) "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Moves the simple obstacle and returns a list of moved obstacles (here itself)
list<SimpleObstacle*> SimpleObstacle::Move( double time,
	double dt, bool const& motherCompositeHasImposedVelocity,
        bool const& motherCompositeHasImposedForce )
{
  list<SimpleObstacle*> movingObstacles;

  // Updates the obstacle translational and angular velocity at time t and 
  // translational and angular motion from t to t+dt and returns whether the 
  // obstacle moved from t to t+dt
  m_ismoving = m_kinematics.ImposedMotion( time, dt, *m_geoRBWC->getCentre() );
  
  // Check whether the composite obstacle it belongs to has an imposed velocity  
  m_ismoving = m_ismoving || motherCompositeHasImposedVelocity;

  // Obstacle motion
  if ( m_ismoving && Obstacle::m_MoveObstacle )
  {
    Vector3 const* translation = m_kinematics.getTranslation();
    m_geoRBWC->composeLeftByTranslation( *translation );
    Quaternion const* w = m_kinematics.getQuaternionRotationOverDt();
    Rotate( *w );
  }

//   bool deplaceF = m_confinement.ImposedMotion( time, dt, this );
//   deplaceF = deplaceF || motherCompositeHasImposedForce;
// 
//   if ( deplaceF && Obstacle::m_MoveObstacle )
//   {
//     Vector3 translation = m_confinement.getTranslation( dt );
//     m_geoRBWC->composeLeftByTranslation( translation );
//     //    Quaternion w = cinematique.getRotation(dt);
//     //    Rotate(w);
//   }
//   m_ismoving = m_ismoving || deplaceF;

  // If the obstacle moved, add it to the list of moving obstacles and 
  // updates its boundind box
  if ( m_ismoving )
  {
    m_obstacleBox = Component::BoundingBox();
    movingObstacles.push_back( this );
  }

  return ( movingObstacles );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the obstacle if the name matches
const Obstacle* SimpleObstacle::getObstacleFromName( string const& nom_ ) const
{
  if ( m_name == nom_ )
    return ( (Obstacle*)this );
  else
    return ( NULL );
}




// ----------------------------------------------------------------------------
// Returns a list of simple obstacles that belong to the obstacle.
// Here returns this (i.e. itself)
list<SimpleObstacle*> SimpleObstacle::getObstacles()
{
  list<SimpleObstacle*> liste;
  liste.push_back(this);
  return ( liste );
}




// ----------------------------------------------------------------------------
// Returns a list of simple obstacles to be sent to the fluid solver
// in case of coupling with a fluid. Here returns this (i.e. itself)
list<Obstacle*> SimpleObstacle::getObstaclesToFluid()
{
  list<Obstacle*> liste;
  if ( m_transferToFluid == true ) liste.push_back(this);
  return ( liste );
}




// ----------------------------------------------------------------------------
// Contact between a simple obstacle and a component. If contact
// exists, computes the contact force and torque and adds to each component
void SimpleObstacle::InterAction( Component* voisin,
	double dt, double const& time, LinkedCell* LC )
{
  try{
  bool contactProbable = m_obstacleBox.InZone( voisin->getPosition(),
    	voisin->getCircumscribedRadius() );

  PointContact closestPoint;
  if ( contactProbable )
  {
    try {
      closestPoint = voisin->getRigidBody()->ClosestPoint( *m_geoRBWC );
    }
    catch (ContactError &erreur)
    {
      try {
	closestPoint = voisin->getRigidBody()->ClosestPoint_ErreurHandling(
		*m_geoRBWC, 10., m_id, voisin->getID() );
      }
      catch (ContactError &erreur_level2)
      {
        cout << endl << "Processor = " <<
    		(GrainsExec::m_MPI ?
			GrainsExec::getComm()->get_rank() : 0 )
		<< " has thrown a ContactError exception" <<  endl;
        erreur_level2.setMessage(
		"SimpleObstacle::InterAction : choc de croute ! a t="
		+GrainsExec::doubleToString(time,TIMEFORMAT));
        erreur_level2.setComponents( this, voisin, time );
        GrainsExec::m_exception_Contact = true;
        throw(erreur_level2);
      }
    }
  }
  else
    closestPoint = PointNoContact;

  LC->addToContactsFeatures( time, closestPoint );

  if ( closestPoint.getOverlapDistance() < 0. )
  {
    if ( ContactBuilderFactory::contactForceModel(
    	m_materialName, voisin->getMaterial() )
    	->computeForces( voisin, this, closestPoint, LC, dt ) )
      voisin->addToCoordinationNumber( 1 );
  }

  }
  catch (const ContactError&) {
    throw ContactError();
  }
}




// ----------------------------------------------------------------------------
// Reloads the simple obstacle and links it to the higher level
// obstacle in the obstacle tree
void SimpleObstacle::reload( Obstacle& mother, istream& file )
{
  string buffer;

  // Obstacle ID number
  file >> m_id;

  // Material name
  file >> buffer >> m_materialName;

  // Construction of the rigid body with crust
  // The flow is RigidBodyWithCrust( fileIn )
  // => ConvexBuilderFactory::create( cle, fileIn )
  // => xxx::create( fileIn )
  // => constructor of xxx
  m_geoRBWC = new RigidBodyWithCrust( file );

  // Whether the obstacle should be transferred to the fluid solver in
  // case of fluid-particle coupling
  file >> buffer >> m_transferToFluid;

  // Obstacle position
  file >> buffer;
  m_geoRBWC->readPosition( file );

  // Read contact map
  readContactMap2014( file );

  // Add the obstacle to the tree
  mother.append( this );
  file >> buffer;
  assert( buffer == "</Simple>" );

  // Compute the bounding box
  m_obstacleBox = Component::BoundingBox();
}




// ----------------------------------------------------------------------------
// Rotates the obstacle with a quaternion
void SimpleObstacle::Rotate( Quaternion const& rotation )
{
  m_geoRBWC->Rotate( rotation );
  m_obstacleBox = Component::BoundingBox();
}




// ----------------------------------------------------------------------------
// Deletes the obstacle if it belongs to a prescribed box
void SimpleObstacle::Suppression( BBox const& box )
{}




// ----------------------------------------------------------------------------
// Translates the obstacle
void SimpleObstacle::Translate( Vector3 const& translation )
{
  m_geoRBWC->composeLeftByTranslation( translation );
  m_obstacleBox = Component::BoundingBox();
}




// ----------------------------------------------------------------------------
// Outputs the simple obstacle for reload
void SimpleObstacle::write( ostream& fileSave ) const
{
  fileSave << "<Simple> " << m_name << " " << m_id << endl;
  fileSave << "*Material " << m_materialName << endl;
  m_geoRBWC->writeStatic( fileSave );
  fileSave << endl;
  fileSave << "*ToFluid " << m_transferToFluid << endl;
  m_geoRBWC->writePosition( fileSave );
  fileSave << endl;
  writeContactMemory2014( fileSave );
  fileSave << endl;  
  fileSave << "</Simple>";
}




// ----------------------------------------------------------------------------
// Adds a cell to the list of cells the obstacle is linked to
void SimpleObstacle::add( Cell *cel_ )
{
  m_inCells.push_back( cel_ );
}




// ----------------------------------------------------------------------------
// Empties the list of cells the obstacle is linked to and deletes
// the pointer to the obstacle is these cells
void SimpleObstacle::resetInCells()
{
  list<Cell*>::iterator il;
  for (il=m_inCells.begin();il!=m_inCells.end();il++) (*il)->remove( this );
  m_inCells.clear();
}




// ----------------------------------------------------------------------------
// Returns the list of cells the obstacle is linked to
list<Cell*> const* SimpleObstacle::getInCells() const
{
  return ( &m_inCells );
}




// ----------------------------------------------------------------------------
// Returns the bounding box of the obstacle
const BBox* SimpleObstacle::getObstacleBox() const
{
  return ( &m_obstacleBox );
}




// ----------------------------------------------------------------------------
// Returns obstacle type
string SimpleObstacle::getObstacleType()
{
  return ( m_ObstacleType );
}




// ----------------------------------------------------------------------------
// Deletes an obstacle in the obstacle tree
void SimpleObstacle::DestroyObstacle( string const& name_ )
{
  if ( m_name == name_ || name_ == "ToBeDestroyed" )
    Obstacle::m_totalNbSingleObstacles--;
}




// ----------------------------------------------------------------------------
// Deletes an obstacle in the obstacle tree and removes it from the LinkedCell
void SimpleObstacle::ClearObstacle( string const& name_, LinkedCell* LC )
{
  if ( m_name == name_ || name_ == "ToBeErased" )
    // Suppression de l'obstacle du LinkedCell
    LC->remove( this );
}




// ----------------------------------------------------------------------------
// Returns whether an update of the link between the obstacle and the
// linked-cell grid is required. If yes, returns true and sets the counter to
// 0, if no, returns false and increments the counter
bool SimpleObstacle::performLinkUpdate()
{
  bool b_doit = false;

  if ( m_LinkUpdate_counter == m_LinkUpdate_frequency )
  {
    m_LinkUpdate_counter = 1;
    b_doit = true;
  }
  else ++m_LinkUpdate_counter;

  return ( b_doit );
}




// ----------------------------------------------------------------------------
// Returns the maximum of the absolute value of the obstacle
// velocity in each direction
Vector3 SimpleObstacle::vitesseMaxPerDirection() const
{
  list<Point3> surface = m_geoRBWC->get_polygonsPts_PARAVIEW();
  list<Point3>::iterator iv;
  Vector3 vmax, v;

  for (iv=surface.begin();iv!=surface.end();iv++)
  {
    v = getVelocityAtPoint(*iv);
    vmax[X] = fabs(v[X]) > vmax[X] ? fabs(v[X]) : vmax[X];
    vmax[Y] = fabs(v[Y]) > vmax[Y] ? fabs(v[Y]) : vmax[Y];
    vmax[Z] = fabs(v[Z]) > vmax[Z] ? fabs(v[Z]) : vmax[Z];
  }

  return ( vmax );

}




// ----------------------------------------------------------------------------
// Sets the frequency at which the obstacle link to
// the cells of the linked-cell grid is updated
void SimpleObstacle::setObstacleLinkedCellUpdateFrequency(
	int const& updateFreq )
{
  m_LinkUpdate_frequency = updateFreq;
  m_LinkUpdate_counter = m_LinkUpdate_frequency;
}




// ----------------------------------------------------------------------------
// Creates obstacle state and adds state to the list of states of all obstacles
void SimpleObstacle::createState( list<struct ObstacleState*> &obsStates )
{
  struct ObstacleState* obss = new ObstacleState;
  obss->nom = m_name;
  obss->memento_config = new ConfigurationMemento();
  obss->memento_config->m_position = *m_geoRBWC->getTransform();
  obss->memento_cine = m_kinematics.createState();
  obsStates.push_back( obss );
}




// ----------------------------------------------------------------------------
// Restores obstacle state
void SimpleObstacle::restoreState( list<struct ObstacleState*>& obsStates )
{
  list<struct ObstacleState*>::iterator il;
  bool found = false ;
  for (il=obsStates.begin();il!=obsStates.end() && !found; )
    if ( (*il)->nom == m_name )
    {
      m_geoRBWC->setTransform((*il)->memento_config->m_position);
      m_kinematics.restoreState((*il)->memento_cine);
      delete (*il)->memento_config;
      delete (*il)->memento_cine;
      delete *il;
      il = obsStates.erase(il) ;
      found = true ;
    }
    else il++;
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the simple obstacle in a
// Paraview format
int SimpleObstacle::numberOfPoints_PARAVIEW() const
{
  return ( m_geoRBWC->getConvex()->numberOfPoints_PARAVIEW() );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the
// simple obstacle shape in a Paraview format
int SimpleObstacle::numberOfCells_PARAVIEW() const
{
  return ( m_geoRBWC->getConvex()->numberOfCells_PARAVIEW() );
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the simple obstacle in a Paraview format
list<Point3> SimpleObstacle::get_polygonsPts_PARAVIEW(
	Vector3 const* translation ) const
{
  return ( m_geoRBWC->get_polygonsPts_PARAVIEW( translation ) );
}




// ----------------------------------------------------------------------------
// Writes the points describing the simple obstacle in a Paraview format
void SimpleObstacle::write_polygonsPts_PARAVIEW( ostream &f,
  	Vector3 const* translation ) const
{
  m_geoRBWC->write_polygonsPts_PARAVIEW( f, translation );
}




// ----------------------------------------------------------------------------
// Writes the simple obstacle in a Paraview format
void SimpleObstacle::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  m_geoRBWC->getConvex()->write_polygonsStr_PARAVIEW( connectivity,
	offsets, cellstype, firstpoint_globalnumber, last_offset );
}




// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
void SimpleObstacle::writePositionInFluid( ostream &fileOut )
{
  m_geoRBWC->writePositionInFluid( fileOut );
}



// ----------------------------------------------------------------------------
// Initializes all contact map entries to false
void SimpleObstacle::setContactMapToFalse()
{
  Component::setContactMapToFalse();
}




// ----------------------------------------------------------------------------
// Set contact map entry features to zero
void SimpleObstacle::setContactMapFeaturesToZero()
{
  Component::setContactMapFeaturesToZero();
}




// ----------------------------------------------------------------------------
// Updates contact map
void SimpleObstacle::updateContactMap()
{
  Component::updateContactMap();
}




// ----------------------------------------------------------------------------
// Does the contact exist in the map, if yes return the pointer to the
// cumulative tangential motion
bool SimpleObstacle::getContactMemory( std::tuple<int,int,int> const& id,
  Vector3* &tangent, Vector3* &prev_normal, Vector3* &cumulSpringTorque,
  bool createContact)
{
  return ( Component::getContactMemory( id, tangent, prev_normal,
    cumulSpringTorque, createContact) );
}




// ----------------------------------------------------------------------------
// Adds new contact in the map
void SimpleObstacle::addNewContactInMap( std::tuple<int,int,int> const& id,
  Vector3 const& tangent, Vector3 const& prev_normal,
  Vector3 const& cumulSpringTorque )
{
  Component::addNewContactInMap( id, tangent, prev_normal, cumulSpringTorque );
}




// ----------------------------------------------------------------------------
// Increases cumulative tangential motion with component id
void SimpleObstacle::addDeplContactInMap( std::tuple<int,int,int> const& id,
  Vector3 const& tangent, Vector3 const& prev_normal,
  Vector3 const& cumulSpringTorque )
{
  Component::addDeplContactInMap( id, tangent, prev_normal, cumulSpringTorque );
}



// ----------------------------------------------------------------------------
// Writes the contact map information in an array of doubles
void SimpleObstacle::copyContactMap( double* destination, 
	int start_index )
{
  Component::copyContactMap( destination, start_index ) ;
}




// ----------------------------------------------------------------------------
// Adds a single contact info to the contact map
void SimpleObstacle::copyContactInMap( std::tuple<int,int,int> const& id,
  bool const& isActive, Vector3 const& tangent, Vector3 const& prev_normal,
  Vector3 const& cumulSpringTorque )
{
  Component::copyContactInMap( id, isActive, tangent, prev_normal,
    cumulSpringTorque ) ;
}



// ----------------------------------------------------------------------------
// Returns the number of contacts in the contact map
int SimpleObstacle::getContactMapSize()
{
  return ( Component::getContactMapSize() );
}




// ----------------------------------------------------------------------------
// Displays the active neighbours in the format "my_elementary_id/neighbour_id/
// neightbout_elementary_id ; ...". Useful for debugging only.
void SimpleObstacle::printActiveNeighbors( int const& id )
{
  Component::printActiveNeighbors( id );
}




// ----------------------------------------------------------------------------
// Resets the minimum ID number of an obstacle for autonumbering */
void SimpleObstacle::setMinIDnumber()
{
  m_minID = min( m_minID, m_id );
}




// ----------------------------------------------------------------------------
// Computes center of mass position
pair<Point3,double> SimpleObstacle::computeCenterOfMass()
{
  return( pair<Point3,double>( *(getPosition()), getVolume() ) );
}
