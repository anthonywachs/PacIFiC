#include "CompositeObstacle.hh"
#include "ObstacleBuilderFactory.hh"
#include "Memento.hh"
#include "LinkedCell.hh"
#include "Torsor.hh"
#include "PointC.hh"
#include "GrainsExec.hh"


// ----------------------------------------------------------------------------
// Constructor with name as input parameter
CompositeObstacle::CompositeObstacle( string const& s ) :
  Obstacle( s, false )
{
  m_id = GrainsExec::m_CompositeObstacleDefaultID;
  m_geoRBWC = new RigidBodyWithCrust( new PointC(), Transform() );
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
CompositeObstacle::CompositeObstacle( DOMNode* root ) :
  Obstacle( "obstacle", false )
{
  m_id = GrainsExec::m_CompositeObstacleDefaultID;
  assert( root != NULL );

  m_geoRBWC = new RigidBodyWithCrust( new PointC(), Transform() );

  m_name = ReaderXML::getNodeAttr_String( root, "name" );

  Obstacle *obstacle = NULL;
  DOMNodeList* allObstacles = ReaderXML::getNodes(root);
  for (XMLSize_t i=0; i<allObstacles->getLength(); i++) 
  {
    obstacle = ObstacleBuilderFactory::create( allObstacles->item( i ) );
    m_obstacles.push_back(obstacle);
  }
  computeCenterOfMass();
}




// ----------------------------------------------------------------------------
// Destructor
CompositeObstacle::~CompositeObstacle()
{
  Obstacle *obstacle;
  list<Obstacle*>::iterator iter;
  for (iter=m_obstacles.begin(); iter!=m_obstacles.end(); iter++) 
  {
    obstacle = *iter;
    delete obstacle;
  }
}




// ----------------------------------------------------------------------------
// Adds an obstacle (single or composite) to the composite obstacle tree
void CompositeObstacle::append( Obstacle* obstacle )
{
  m_obstacles.push_back( obstacle );
}




// ----------------------------------------------------------------------------
// Links imposed kinematics to the obstacle and returns true if the
// linking process is successful
bool CompositeObstacle::LinkImposedMotion( ObstacleImposedVelocity* imposed )
{
  bool status = false;
  if ( m_name == imposed->getObstacleName() ) 
  {
    m_kinematics.append( imposed );
    status = true;
  } 
  else 
  {
    list<Obstacle*>::iterator obstacle;
    for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
      status = (*obstacle)->LinkImposedMotion( imposed );
  }
  
  return ( status );  
}




// ----------------------------------------------------------------------------
// Links imposed force kinematics to the obstacle and returns true 
// if the linking process is successful
bool CompositeObstacle::LinkImposedMotion( ObstacleImposedForce* imposed )
{
  bool status = false;
  if ( m_name == imposed->getObstacleName() ) 
  {
    m_confinement.append( imposed );
    status = true;
  } 
  else 
  {
    list<Obstacle*>::iterator obstacle;
    for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
      status = (*obstacle)->LinkImposedMotion( imposed );
  }
  
  return ( status );  
}




// ----------------------------------------------------------------------------
// Moves the composite obstacle and returns a list of moved obstacles
list<SimpleObstacle*> CompositeObstacle::Move( double time, double dt, 
	bool const& motherCompositeHasImposedVelocity, 
        bool const& motherCompositeHasImposedForce )
{
  list<Obstacle*>::iterator obstacle;
  list<SimpleObstacle*> movingObstacles;
  list<SimpleObstacle*>::iterator ilo; 
  bool imposedForce = false;   

  // Updates the obstacle translational and angular velocity at time t and 
  // translational and angular motion from t to t+dt and returns whether the 
  // obstacle moved from t to t+dt
  m_ismoving = m_kinematics.ImposedMotion( time, dt, *m_geoRBWC->getCentre() );
  
  // Check whether the composite obstacle it belongs to has an imposed velocity
  m_ismoving = m_ismoving || motherCompositeHasImposedVelocity;

  // Composite center motion
  if ( m_ismoving && Obstacle::m_MoveObstacle ) 
  {
    // Translation motion
    Vector3 const* translation = m_kinematics.getTranslation();
    m_geoRBWC->composeLeftByTranslation( *translation );

    // Angular motion
    Quaternion const* w = m_kinematics.getQuaternionRotationOverDt();
    m_geoRBWC->Rotate( *w );
  }
  
  // Updated center of the composite
  Point3 centre = *getPosition();

  // Apply the composite imposed motion to its elementary obstacles
  if ( m_ismoving ) 
  {
    Vector3 lever;
    for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++) 
    {
      lever = *(*obstacle)->getPosition() - centre;
      (*obstacle)->Compose( m_kinematics, lever );
    }
  }

      
//   bool deplaceF = m_confinement.ImposedMotion( time, dt, this );
//   deplaceF = deplaceF || motherCompositeHasImposedForce;    
//   
//   if ( deplaceF ) 
//     for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
//       (*obstacle)->Compose( m_confinement, *(*obstacle)->getPosition() );
// 
//   m_ismoving = m_ismoving || deplaceF; // ??? demander à Gillos !!

  
  // Finally, move the elementary obstacles     
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++) 
  {
    list<SimpleObstacle*> lod = (*obstacle)->Move( time, dt, m_ismoving, 
    	imposedForce );
    for (ilo=lod.begin();ilo!=lod.end();ilo++) 
      movingObstacles.push_back(*ilo);
  }
  
  return ( movingObstacles );
}




// ----------------------------------------------------------------------------
// Computes center of mass position
pair<Point3,double> CompositeObstacle::computeCenterOfMass()
{
  // Notes: 
  // 1) Obstacles do not have any density, so we assume that all obstacles 
  // have the same density when computing the center of mass
  // 2) This method ignores geometric overlaps between obstacles so this
  // computation might be erroneous if obstacles overlap signigicantly

  Point3 centre;
  double mass = 0.;
  int nbre = 0;
  list<Obstacle*>::iterator obstacle;
  pair<Point3,double> pp;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); 
      nbre++, obstacle++)
  {
    pp = (*obstacle)->computeCenterOfMass();
    centre += pp.first * pp.second;
    mass += pp.second; 
  }
  centre /= mass;
  setPosition( centre );
  
  pp.first = centre;
  pp.second = mass;
  
  return ( pp );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the obstacle if the name matches. Searches all
// leaves of the composite obstacle tree
Obstacle const* CompositeObstacle::getObstacleFromName( const string& nom_ ) 
	const
{
  const Obstacle *obst = NULL;
  if ( m_name == nom_ ) obst = (Obstacle*)this;
  else 
  {
    list<Obstacle*>::const_iterator obstacle;
    for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end() && 
    	obst == NULL; obstacle++)
      obst = (*obstacle)->getObstacleFromName( nom_ );
  }
  
  return ( obst );
}




// ----------------------------------------------------------------------------
// Returns a list of simple obstacles that belong to the composite obstacle
list<SimpleObstacle*> CompositeObstacle::getObstacles()
{
  list<SimpleObstacle*> liste;
  
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++) 
  {
    list<SimpleObstacle*> tmp = (*obstacle)->getObstacles();
    liste.insert(liste.end(), tmp.begin(), tmp.end());
  }
  
  return ( liste );
}




// ----------------------------------------------------------------------------
// Returns a list of simple obstacles to be sent to the fluid solver
// in case of coupling with a fluid
list<Obstacle*> CompositeObstacle::getObstaclesToFluid()
{
  list<Obstacle*> liste;
  
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++) 
  {
    list<Obstacle*> tmp = (*obstacle)->getObstaclesToFluid();
    liste.insert(liste.end(), tmp.begin(), tmp.end());
  }
  
  return ( liste );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another component. 
bool CompositeObstacle::isContact( Component const* voisin ) const
{
  bool contact = false;
    
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); 
       obstacle!=m_obstacles.end() && !contact; obstacle++)
  {
    if ( voisin->isCompositeParticle() )
      contact = voisin->isContact( *obstacle );
    else
      contact = (*obstacle)->isContact( voisin );
  }

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another
// component accounting for crust thickness
bool CompositeObstacle::isContactWithCrust( Component const* voisin ) const
{
  bool contact = false;
    
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); 
       obstacle!=m_obstacles.end() && !contact; obstacle++)
  {
    if ( voisin->isCompositeParticle() )
      contact = voisin->isContactWithCrust( *obstacle );
    else
      contact = (*obstacle)->isContactWithCrust( voisin );
  }

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric proximity with another
// component in the sense of whether their respective bounding boxes overlap 
bool CompositeObstacle::isClose( Component const* voisin ) const
{
  bool contact = false;
    
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); 
       obstacle!=m_obstacles.end() && !contact; obstacle++)
    if ( voisin->isCompositeParticle() )
      contact = voisin->isClose( *obstacle );
    else
      contact = (*obstacle)->isClose( voisin );

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric proximity with another
// component in the sense of whether their respective bounding boxes minus 
// their crust thickness overlap
bool CompositeObstacle::isCloseWithCrust( Component const* voisin ) const
{
  bool contact = false;
    
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); 
       obstacle!=m_obstacles.end() && !contact; obstacle++)
    if ( voisin->isCompositeParticle() )
      contact = voisin->isCloseWithCrust( *obstacle );
    else
      contact = (*obstacle)->isCloseWithCrust( voisin );

  return ( contact );
}




// ----------------------------------------------------------------------------
// Reloads the composite obstacle and links it to the higher level 
// obstacle in the obstacle tree
void CompositeObstacle::reload( Obstacle& mother, istream& file )
{
  string ttag;
  file >> ttag;
  while ( ttag != "</Composite>" ) 
  {
    ObstacleBuilderFactory::reload( ttag, *this, file );
    file >> ttag;
  }
  computeCenterOfMass();
  mother.append(this);
}    




// ----------------------------------------------------------------------------
// Resets kinematics to 0
void CompositeObstacle::resetKinematics()
{
  m_kinematics.reset();
  m_confinement.reset();

  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->resetKinematics();
}




// ----------------------------------------------------------------------------
// Rotates the composite obstacle with a quaternion
void CompositeObstacle::Rotate( Quaternion const& rotation )
{
  cout << "Warning when calling CompositeObstacle::Rotate(rotation) "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Deletes the composite obstacle if it belongs to a prescribed box
void CompositeObstacle::Suppression( BBox const& box )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); ) 
  {
    (*obstacle)->Suppression( box );
    if ( (*obstacle)->isIn( box ) ) obstacle = m_obstacles.erase(obstacle);
    else obstacle++;
  }
}




// ----------------------------------------------------------------------------
// Translates the composite obstacle
void CompositeObstacle::Translate( Vector3 const& translation )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->Translate( translation );
}




// ----------------------------------------------------------------------------
// Writes the composite obstacle position
void CompositeObstacle::writePosition( ostream& position )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->writePosition( position );
}




// ----------------------------------------------------------------------------
// Writes the composite obstacle's "static" data
void CompositeObstacle::writeStatic( ostream& fileOut ) const
{
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->writeStatic( fileOut );
}




// ----------------------------------------------------------------------------
// Outputs the composite obstacle for reload
void CompositeObstacle::write( ostream& fileSave ) const
{
  fileSave << "<Composite> " << m_name << endl;
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
  {
    (*obstacle)->write( fileSave );
    fileSave << endl;
  }
  fileSave << "</Composite>";
}




// ----------------------------------------------------------------------------
// Deletes an obstacle in the obstacle tree
void CompositeObstacle::DestroyObstacle( string const& name_ )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); ) 
  {
    if ( (*obstacle)->getName() == name_ ||
    	 (*obstacle)->getName() == "ToBeDestroyed" )
    {
      (*obstacle)->DestroyObstacle( "ToBeDestroyed" );
      delete *obstacle;
      obstacle = m_obstacles.erase(obstacle);
    }
    else 
    {
      (*obstacle)->DestroyObstacle( name_ );
      obstacle++;
    }
  }       
}  




// ----------------------------------------------------------------------------
// Deletes an obstacle in the obstacle tree and removes 
// it from the LinkedCell
void CompositeObstacle::ClearObstacle( string const& name_, LinkedCell* LC )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end();obstacle++) 
  {
    if ( (*obstacle)->getName() == name_ )
    {
      list<SimpleObstacle*> allObs = (*obstacle)->getObstacles();
      list<SimpleObstacle*>::iterator il;
      for(il=allObs.begin(); il!=allObs.end(); il++)
        (*il)->ClearObstacle( "ToBeErased", LC );
    }
    else 
      (*obstacle)->ClearObstacle( name_, LC );
  }   
}  




// ----------------------------------------------------------------------------
// Updates indicator for Paraview post-processing
void CompositeObstacle::updateIndicator( double time, double dt )
{
  list<Obstacle*>::iterator obstacle;
  
  if ( m_kinematics.activAngularMotion( time, dt ) )
      getObstacles().front()->setIndicator( 1. );  
  
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)   
    (*obstacle)->updateIndicator( time, dt );
}  




// ----------------------------------------------------------------------------
// Creates obstacle state and adds state to the list of states of all obstacles 
void CompositeObstacle::createState( list<struct ObstacleState*>& obsStates )
{
  struct ObstacleState* obss = new ObstacleState;
  obss->nom = m_name;
  obss->memento_config = new ConfigurationMemento();  
  obss->memento_config->m_position = *m_geoRBWC->getTransform();
  obss->memento_cine = m_kinematics.createState();
  obsStates.push_back( obss );
  
  for (list<Obstacle*>::const_iterator obstacle=m_obstacles.begin(); 
  	obstacle!=m_obstacles.end();obstacle++)
    (*obstacle)->createState( obsStates );  
}




// ----------------------------------------------------------------------------
// Restores obstacle state
void CompositeObstacle::restoreState( list<struct ObstacleState*>& obsStates )
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
    
  for (list<Obstacle*>::iterator obstacle=m_obstacles.begin(); 
  	obstacle!=m_obstacles.end();obstacle++)
    (*obstacle)->restoreState( obsStates );       
} 




// ----------------------------------------------------------------------------
// Initializes the composite obstacle's torsor
void CompositeObstacle::InitializeForce( bool const& withWeight )
{
  for (list<Obstacle*>::iterator obstacle=m_obstacles.begin(); 
  	obstacle!=m_obstacles.end(); obstacle++)   
    (*obstacle)->InitializeForce( false );
}    




// ----------------------------------------------------------------------------
// Returns a pointer to the torsor exerted on the composite obstacle
Torsor const* CompositeObstacle::getTorsor()
{
  m_torsor.setToBodyForce( *getPosition(), Vector3Null );  

  for (list<Obstacle*>::iterator obstacle=m_obstacles.begin(); 
  	obstacle!=m_obstacles.end(); obstacle++)
    m_torsor += *(*obstacle)->getTorsor();

  return ( &m_torsor );  
} 




// ----------------------------------------------------------------------------
// Returns obstacle type
string CompositeObstacle::getObstacleType()
{
  return ( m_ObstacleType ) ;
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the composite obstacle in a
// Paraview format
int CompositeObstacle::numberOfPoints_PARAVIEW() const
{
  cout << "Warning when calling CompositeObstacle::numberOfPoints_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the 
// composite obstacle shape in a Paraview format
int CompositeObstacle::numberOfCells_PARAVIEW() const
{
  cout << "Warning when calling CompositeObstacle::numberOfCells_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the composite obstacle in a
// Paraview format
list<Point3> CompositeObstacle::get_polygonsPts_PARAVIEW( 
	Vector3 const* translation ) const
{
  cout << "Warning when calling CompositeObstacle::get_polygonsPts_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Writes the points describing the composite obstacle in a Paraview format
void CompositeObstacle::write_polygonsPts_PARAVIEW( ostream& f, 
  	Vector3 const* translation ) const
{
  cout << "Warning when calling CompositeObstacle::write_polygonsPts_PARAVIEW()"
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Writes the composite obstacle in a Paraview format
void CompositeObstacle::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  cout << "Warning when calling CompositeObstacle::write_polygonsStr_PARAVIEW()"
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
void CompositeObstacle::writePositionInFluid( ostream& fileOut )
{
  cout << "Warning when calling CompositeObstacle::writePositionInFluid() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Initialize all contact map entries to false
void CompositeObstacle::setContactMapToFalse()
{
  cout << "Warning when calling CompositeObstacle::setContactMapToFalse() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Set contact map entry features to zero */
void CompositeObstacle::setContactMapFeaturesToZero()
{
  cout << "Warning when calling CompositeObstacle::"
       << "setContactMapFeaturesToZero() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Updates contact map
void CompositeObstacle::updateContactMap()
{
  cout << "Warning when calling CompositeObstacle::updateContactMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Does the contact exist in the map? If so, return true and make
// kdelta, prev_normal and cumulSpringTorque point to the memorized info. 
// Otherwise, return false and set those pointers to NULL.
bool CompositeObstacle::getContactMemory( std::tuple<int,int,int> const& id,
  	Vector3* &kdelta, Vector3* &prev_normal, Vector3* &cumulSpringTorque,
  	bool createContact )
{
  cout << "Warning when calling CompositeObstacle::getContactMemory() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Adds new contact in the map
void CompositeObstacle::addNewContactInMap( std::tuple<int,int,int> const& id,
  	Vector3 const& kdelta, Vector3 const& prev_normal,
  	Vector3 const& cumulSpringTorque )
{
  cout << "Warning when calling CompositeObstacle::addNewContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Increases cumulative tangential motion with component id
void CompositeObstacle::addDeplContactInMap( std::tuple<int,int,int> const& id,
  	Vector3 const& kdelta, Vector3 const& prev_normal,
  	Vector3 const& cumulSpringTorque )
{
  cout << "Warning when calling CompositeObstacle::addDeplContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Writes the contact map information in an array of doubles
void CompositeObstacle::copyContactMap( double* destination, 
	int start_index )
{
  cout << "Warning when calling CompositeObstacle::copyContactMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Adds a single contact info to the contact map
void CompositeObstacle::copyContactInMap( std::tuple<int,int,int> const& id,
  	bool const& isActive, Vector3 const& kdelta, Vector3 const& prev_normal,
  	Vector3 const& cumulSpringTorque )
{
  cout << "Warning when calling CompositeObstacle::copyContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}	
	
	
	

// ----------------------------------------------------------------------------
// Returns the number of contacts in the contact map */
int CompositeObstacle::getContactMapSize()
{
  cout << "Warning when calling CompositeObstacle::getContactMapSize() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Displays the active neighbours in the format "my_elementary_id/neighbour_id/
// neightbout_elementary_id ; ...". Useful for debugging only.
void CompositeObstacle::printActiveNeighbors( int const& id )
{
  cout << "Warning when calling CompositeObstacle::printActiveNeighbors() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the composite obstacle
bool CompositeObstacle::isIn( Point3 const& pt ) const
{
  bool bisIn = false;
  
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); 
       obstacle!=m_obstacles.end() && !bisIn; obstacle++)
    bisIn = (*obstacle)->isIn( pt );

  return ( bisIn );  
}




// ----------------------------------------------------------------------------
// Returns whether the component is a composite obstacle ? */
bool CompositeObstacle::isCompositeObstacle() const
{
  return ( true ); 
}




// ----------------------------------------------------------------------------
// Resets the minimum ID number of an obstacle for autonumbering */
void CompositeObstacle::setMinIDnumber()
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++) 
    (*obstacle)->setMinIDnumber();
}
