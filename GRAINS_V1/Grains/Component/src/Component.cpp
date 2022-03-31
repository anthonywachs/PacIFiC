#include "Component.hh"
#include "Particle.hh"
#include "Obstacle.hh"
#include "Torsor.hh"
#include "Error.hh"
#include "Memento.hh"
#include <algorithm>
using namespace std;


int Component::m_nb = 0;


// ----------------------------------------------------------------------------
// Default constructor
Component::Component( bool const& autonumbering )
  : m_materialName( "" )
  , m_mass( 0. )
  , m_geoRBWC( NULL )
  , m_memento( NULL )
{
  if ( autonumbering )
  {
    m_id = Component::m_nb;
    Component::m_nb++;
  }
  else m_id = -1;
}




// ----------------------------------------------------------------------------
// Copy constructor
Component::Component( Component const& copy )
  : m_id( Component::m_nb )
  , m_materialName( copy.m_materialName )
  , m_mass( copy.m_mass )
  , m_memento( NULL )
{
  Component::m_nb++;

  m_geoRBWC = new RigidBodyWithCrust( *copy.m_geoRBWC );
}




// ----------------------------------------------------------------------------
// Destructor
Component::~Component()
{
  Component::m_nb--;
  delete m_geoRBWC;
  if ( m_memento ) delete m_memento;
}




// ----------------------------------------------------------------------------
// Adds a force exerted at a point  to the torsor (torsor adds torque
// automatically)
void Component::addForce( Point3 const& point, Vector3 const& force )
{
  m_torsor.addForce( point, force );
}




// ----------------------------------------------------------------------------
// Adds a body force exerted at the center of mass of the component
// (associated torque is 0)
void Component::addBodyForce( Vector3 const& force )
{
  m_torsor.addForce( force );
}




// ----------------------------------------------------------------------------
// Adds a torque to the torsor
void Component::addTorque( Vector3 const& torque )
{
  m_torsor.addTorque( torque );
}




// ----------------------------------------------------------------------------
// Returns the bounding box of the component
BBox Component::BoundingBox() const
{
  // We use the function BoxRigidBody from RigidBody and not from
  // RigidBodyWithCrust as the function from RigidBodyWithCrust extends the
  // bounding box by the crust thickness of the rigid body with crust.
  // m_geoRBWC is a pointer of type RigidBodyWithCrust
  return ( m_geoRBWC->RigidBody::BoxRigidBody() );
}




// ----------------------------------------------------------------------------
// Returns a point to the rigid body with crust
RigidBodyWithCrust const* Component::getRigidBody() const
{
  return ( m_geoRBWC );
}




// ----------------------------------------------------------------------------
// Returns a point to the rigid body with crust
RigidBodyWithCrust* Component::getRigidBody()
{
  return ( m_geoRBWC );
}




// ----------------------------------------------------------------------------
// Returns the component ID number
int Component::getID() const
{
  return ( m_id );
}




// ----------------------------------------------------------------------------
// Returns the mass of the component
double Component::getMass() const
{
  return m_mass;
}




// ----------------------------------------------------------------------------
// Returns a pointer to the center of mass of the component
Point3 const* Component::getPosition() const
{
  return ( m_geoRBWC->getCentre() );
}




// ----------------------------------------------------------------------------
// Copy the center of mass of the component in the array pos
void Component::getPosition( double *pos ) const
{
  Point3 const* pot = m_geoRBWC->getCentre();
  for (int i=0; i<3; i++) pos[i] = (*pot)[i];
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the rigid body
double Component::getCircumscribedRadius() const
{
  return ( m_geoRBWC->getCircumscribedRadius() );
}




// ----------------------------------------------------------------------------
// Returns the radius ofthe sphere of equivalent volume as the rigid body
double Component::getEquivalentSphereRadius() const
{
  return ( 0. );
}




// ----------------------------------------------------------------------------
// Returns the crust thickness of the rigid body
double Component::getCrustThickness() const
{
  return ( m_geoRBWC->getCrustThickness() );
}




// ----------------------------------------------------------------------------
// Returns the volume of the rigid body
double Component::getVolume() const
{
  return ( m_geoRBWC->getVolume() );
}




// ----------------------------------------------------------------------------
// Copies the center of mass of the component in a 1D array
void Component::copyPosition( double* pos, int i ) const
{
  Point3 const* pot = m_geoRBWC->getCentre();
  for (int j=0 ;j<3; j++) pos[i+j] = (*pot)[j];
}





// ----------------------------------------------------------------------------
// Copies the component transformation in a 1D array
void Component::copyTransform( double* vit, int i ) const
{
  m_geoRBWC->copyTransform( vit, i );
}




// ----------------------------------------------------------------------------
// Copies the component transformation in a 1D array with the center
// of mass translated
void Component::copyTransform( double* vit, int i, Vector3 const& vec ) const
{
  m_geoRBWC->copyTransform( vit, i, vec );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another component.
// Note: the other component must not be of the derived type CompositeObstacle
bool Component::isContact( Component const* voisin ) const
{
  if ( this != voisin )
    return m_geoRBWC->RigidBody::isContact( *voisin->m_geoRBWC );
  else return false;
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another component accounting
// for crust thickness. Note: the other component must not be of the derived
// type CompositeObstacle
bool Component::isContactWithCrust( Component const* voisin ) const
{
  if ( this != voisin )
    return m_geoRBWC->isContact( *voisin->m_geoRBWC );
  else return false;
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric proximity with another
// component in the sense of whether their respective bounding boxes overlap.
// Note: the other component must not be of the derived type CompositeObstacle
bool Component::isClose( Component const* voisin ) const
{
  if ( this != voisin )
    return m_geoRBWC->RigidBody::isClose( *voisin->m_geoRBWC );
  else return false;
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric proximity with another
// component in the sense of whether their respective bounding boxes minus
// their crust thickness overlap. Note: the other component must not be of the
// derived type CompositeObstacle
bool Component::isCloseWithCrust( Component const* voisin ) const
{
  if ( this != voisin )
    return m_geoRBWC->isClose( *voisin->m_geoRBWC );
  else return false;
}




// ----------------------------------------------------------------------------
// Returns whether the center of mass of the component is located in
// a geometric box
bool Component::isIn( BBox const& boite ) const
{
  const Point3   &origin = boite.getCenter();
  const Vector3 &extent = boite.getExtent();

  Point3 const* centre = m_geoRBWC->getCentre();
  Vector3 dist = *centre - origin;

  bool status =
    fabs( dist[X] ) <= extent[X] &&
    fabs( dist[Y] ) <= extent[Y] &&
    fabs( dist[Z] ) <= extent[Z];

  return status;
}




// ----------------------------------------------------------------------------
// Returns the component material name
string Component::getMaterial() const
{
  return ( m_materialName );
}




// ----------------------------------------------------------------------------
// Sets the material type
void Component::setMaterial( string const& mat )
{
  m_materialName = mat ;
}




// ----------------------------------------------------------------------------
// Rotates the component using a quaternion
void Component::Rotate( Quaternion const& rotation )
{
  m_geoRBWC->Rotate( rotation );
}




// ----------------------------------------------------------------------------
// Translates the component
void Component::Translate( Vector3 const& translation )
{
  m_geoRBWC->composeLeftByTranslation( translation );
}




// ----------------------------------------------------------------------------
// Sets the component's transformation with an 1D array of 12
// values (see class Transform for details)
void Component::setPosition( double const* pos )
{
  m_geoRBWC->setTransform( pos );
}




// ----------------------------------------------------------------------------
// Sets the origin of the component's transformation
void Component::setPosition( Point3 const& centre )
{
  m_geoRBWC->setOrigin( (double *) &centre );
}




// ----------------------------------------------------------------------------
// Output the data related to the position of the component
void Component::writePosition( ostream& position ) const
{
  m_geoRBWC->writePosition( position );
}




// ----------------------------------------------------------------------------
// Writes the component's "static" data
void Component::writeStatic( ostream &fileOut ) const
{
  fileOut << "*IDnumber " << m_id << endl;
  fileOut << "*Material " << m_materialName << endl;
  m_geoRBWC->writeStatic( fileOut );
  fileOut << endl;
}




// ----------------------------------------------------------------------------
// Returns a pointer to the total force exerted on the component
Vector3 const* Component::getForce() const
{
  return ( m_torsor.getForce() );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the total torque exerted on the component
Vector3 const* Component::getTorque() const
{
  return ( m_torsor.getTorque() );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the torsor exerted on the component
Torsor const* Component::getTorsor()
{
  return ( &m_torsor );
}




// ----------------------------------------------------------------------------
// Initialise le torseur des efforts sur le composant
void Component::InitializeForce( bool const& withWeight )
{
  m_torsor.setToBodyForce( *m_geoRBWC->getCentre(), Vector3Nul );
}




// ----------------------------------------------------------------------------
// Saves the component state
void Component::saveConfigState()
{
  if ( !m_memento ) m_memento = new ConfigurationMemento();
  m_memento->m_position = *m_geoRBWC->getTransform();
}




// ----------------------------------------------------------------------------
// Creates and returns the component state
ConfigurationMemento* Component::createConfigState()
{
  ConfigurationMemento* Pmemento_ = new ConfigurationMemento();
  Pmemento_->m_position = *m_geoRBWC->getTransform();

  return Pmemento_;
}




// ----------------------------------------------------------------------------
// Returns a pointer to the reference component of the component:
// this in general and the CompositeParticle for an elementary particle
Component* Component::getMasterComponent()
{
  return ( this );
}




// ----------------------------------------------------------------------------
// Returns the fluid velocity interpolated at the center of mass of
// the component
Vector3 const* Component::getTranslationalVelocity_fluide() const
{
  cout << "WARNING!!!!Component::getTranslationalVelocity_fluide()"
  "This should not be called here, refer to Particle.cpp" << endl;
  return ( NULL );
}




// ----------------------------------------------------------------------------
// Searches and stores all contact points between a composite particle and a
// component.
void Component::SearchContact( Component* voisin, double dt,
      double const& time, LinkedCell *LC, list<ContactInfos*>& listContact )
  throw ( ContactError )
{
  cout << "WARNING !!! Component::SearchContact"
  	" should never be called, this is a design problem !!!" << endl;
}




// ----------------------------------------------------------------------------
// Increments the coordination number by nc
void Component::addToCoordinationNumber( int const& nc )
{}




// ----------------------------------------------------------------------------
// Initialize all contact map entries to false
void Component::setContactMapToFalse()
{
  map<std::tuple<int,int,int>,std::tuple<bool, Vector3, Vector3,
    Vector3>,bool >::iterator it;
  for (it=m_contactMap.begin();it!=m_contactMap.end();++it)
    {
      get<0>(it->second) = false;
    }
}




// ----------------------------------------------------------------------------
// Update contact map
void Component::updateContactMap()
{
  map<std::tuple<int,int,int>,std::tuple<bool, Vector3, Vector3,
    Vector3>>::iterator it;
  list<std::tuple<int,int,int>> keywithfalse;
  list<std::tuple<int,int,int>>::const_iterator il;

  for (it=m_contactMap.begin();it!=m_contactMap.end();++it)
    if ( !(get<0>(it->second)) ) keywithfalse.push_back(it->first);

  for (il=keywithfalse.begin();il!=keywithfalse.end();il++)
  {
    m_contactMap.erase(*il);
  }
}


// ----------------------------------------------------------------------------
// Add new contact in the map
void Component::addNewContactInMap( std::tuple<int,int,int> const& id,
  Vector3 const& kdelta, Vector3 const& prev_normal,
  Vector3 const& cumulSpringTorque )
{
  m_contactMap.insert(std::make_pair( std::make_tuple(get<0>(id),get<1>(id),
    get<2>(id)), std::make_tuple(
    true, kdelta, prev_normal, cumulSpringTorque) ));
}


// ----------------------------------------------------------------------------
// Copy existing contact in the map
void Component::copyContactInMap( std::tuple<int,int,int> const& id,
  bool const& isActive, Vector3 const& kdelta, Vector3 const& prev_normal,
  Vector3 const& cumulSpringTorque )
{
  m_contactMap.insert(std::make_pair( id, std::make_tuple(
    isActive, kdelta, prev_normal, cumulSpringTorque) ));
}



// ----------------------------------------------------------------------------
// Does the contact exist in the map? If so, return true and make kdelta,
// prev_normal and cumulSpringTorque point to the memorized info. Otherwise,
// return false and set those pointers to NULL
bool Component::getContactMemory( std::tuple<int,int,int> const& id,
  Vector3* &kdelta, Vector3* &prev_normal, Vector3* &cumulSpringTorque,
  bool createContact)
{
  bool active = false;
  map<std::tuple<int,int,int>,std::tuple<bool, Vector3, Vector3, Vector3> >
    ::iterator it = m_contactMap.find(id);

  if ( it != m_contactMap.end() )
  {
    active = true;
    get<0>(it->second) = true;
    kdelta = &(get<1>(it->second));
    prev_normal = &(get<2>(it->second));
    cumulSpringTorque = &(get<3>(it->second));
  }
  else {
    if (createContact) {
      Vector3 zeroV(0.);
      std::pair< map<std::tuple<int,int,int>,std::tuple<bool, Vector3, Vector3,
        Vector3> > ::iterator, bool > ret;
      ret = m_contactMap.insert(std::make_pair( id, std::make_tuple(
        true, zeroV, zeroV, zeroV)));
      it = ret.first;
      kdelta = &(get<1>(it->second));
      prev_normal = &(get<2>(it->second));
      cumulSpringTorque = &(get<3>(it->second));
    }
    else {
      kdelta = NULL;
      prev_normal = NULL;
      cumulSpringTorque = NULL;
    }
  }
  return active;
}


// ----------------------------------------------------------------------------
// Increase cumulative tangential displacement with component id
void Component::addDeplContactInMap( std::tuple<int,int,int> const& id,
  Vector3 const& kdelta, Vector3 const& prev_normal,
  Vector3 const& cumulSpringTorque )
{
  get<0>(m_contactMap[id]) = true;
  get<1>(m_contactMap[id]) = kdelta;
  get<2>(m_contactMap[id]) = prev_normal;
  get<3>(m_contactMap[id]) = cumulSpringTorque;
}

// ---------------------------------------------------------------------------
// Print active neighbors of the particle
void Component::printActiveNeighbors(int const& id )
{
    map<std::tuple<int,int,int>,std::tuple<bool, Vector3, Vector3, Vector3> >
      ::iterator it;
    if (m_contactMap.begin() != m_contactMap.end())
    {
        cout << "Neighbors of #" << id << ": ";
        for (it=m_contactMap.begin();it!=m_contactMap.end();++it){
            if (get<0>(it->second)){
                cout << get<0>(it->first)  << "/"
                << get<1>(it->first)  << "/"
                << get<2>(it->first)  << " ; ";
            }
        }
        cout << endl;
    }
}


void Component::copyHistoryContacts( double* &destination, int start_index )
{
  int nb_contacts = (int) m_contactMap.size();
  destination[start_index] = nb_contacts;
  start_index++;
  map<std::tuple<int,int,int>,std::tuple<bool, Vector3, Vector3, Vector3> >
    ::iterator it;
  if (m_contactMap.begin() != m_contactMap.end())
  {
    for (it=m_contactMap.begin();it!=m_contactMap.end();++it)
    {
      destination[start_index] = (double)get<0>(it->first) ;
      destination[start_index + 1] = (double)get<1>(it->first) ;
      destination[start_index + 2] = (double)get<2>(it->first) ;
      destination[start_index + 3] = (double)get<0>(it->second) ;
      destination[start_index + 4] = get<1>(it->second)[0] ;
      destination[start_index + 5] = get<1>(it->second)[1] ;
      destination[start_index + 6] = get<1>(it->second)[2] ;
      destination[start_index + 7] = get<2>(it->second)[0] ;
      destination[start_index + 8] = get<2>(it->second)[1] ;
      destination[start_index + 9] = get<2>(it->second)[2] ;
      destination[start_index + 10] = get<3>(it->second)[0] ;
      destination[start_index + 11] = get<3>(it->second)[1] ;
      destination[start_index + 12] = get<3>(it->second)[2] ;
      start_index += 13 ;
    }
  }
}

int Component::getContactMapSize()
{
  return ( (int) m_contactMap.size() );
}

void Component::writeContactMemory_2014(ostream &fileOut ) const
{
  int mapSize;
  mapSize = (int) m_contactMap.size();
  fileOut << mapSize ;
  fileOut << " ";
  map<std::tuple<int,int,int>,std::tuple<bool, Vector3, Vector3, Vector3> >
    ::const_iterator it;
  if (m_contactMap.begin() != m_contactMap.end())
  {
  int c=0;
  for (it=m_contactMap.begin();it!=m_contactMap.end();++it)
    {
      c++;
      fileOut << get<0>(it->first) ;
      fileOut << " ";
      fileOut << get<1>(it->first) ;
      fileOut << " ";
      fileOut << get<2>(it->first) ;
      fileOut << " ";
      fileOut << get<0>(it->second) ;
      fileOut << " ";
      fileOut << " ";
      get<1>(it->second).writeGroup3(fileOut);
      fileOut << " ";
      get<2>(it->second).writeGroup3(fileOut);
      fileOut << " ";
      get<3>(it->second).writeGroup3(fileOut);
      if (c<mapSize) fileOut << " ";
    }
  }
}

void Component::writeContactMemory_binary( ostream &fileOut )
{
  int mapSize;
  mapSize = (int) m_contactMap.size();
  fileOut.write( reinterpret_cast<char*>( &mapSize ), sizeof(int) );
  map<std::tuple<int,int,int>,std::tuple<bool, Vector3, Vector3, Vector3> >
    ::const_iterator it;
  if (m_contactMap.begin() != m_contactMap.end())
  {
  for (it=m_contactMap.begin();it!=m_contactMap.end();++it)
    {
      int buffer_int;
      Vector3 buffer_vect;
      buffer_int = get<0>(it->first);
      fileOut.write(reinterpret_cast<char*>(&buffer_int ), sizeof(int));
      buffer_int = get<1>(it->first);
      fileOut.write(reinterpret_cast<char*>(&buffer_int ), sizeof(int));
      buffer_int = get<2>(it->first);
      fileOut.write(reinterpret_cast<char*>(&buffer_int ), sizeof(int));
      buffer_int = get<0>(it->second);
      fileOut.write(reinterpret_cast<char*>(&buffer_int ), sizeof(int));
      buffer_vect = Vector3(get<1>(it->second));
      buffer_vect.writeGroup3_binary(fileOut);
      buffer_vect = Vector3(get<2>(it->second));
      buffer_vect.writeGroup3_binary(fileOut);
      buffer_vect = Vector3(get<3>(it->second));
      buffer_vect.writeGroup3_binary(fileOut);
    }
  }
}


void Component::updateContactMapId( int prev_id, int new_id)
{
  map<std::tuple<int,int,int>,std::tuple<bool, Vector3, Vector3, Vector3> >
    ::iterator it;
  for(it=m_contactMap.begin();it!=m_contactMap.end();++it)
  {
    if (get<1>((it->first)) == prev_id)
    {
      copyContactInMap( std::make_tuple(get<0>(it->first),new_id,
            get<2>(it->first)), get<0>(it->second), get<1>(it->second), get<2>(it->second), get<3>(it->second));
      m_contactMap.erase(std::make_tuple(get<0>(it->first),get<1>(it->first),
            get<2>(it->first)));
      break;
    }
  }
}

void Component::readContactMap_2014( istream &fileSave)
{
  // Read the contact memories of the particle (if any)
  char next_char;
  int contact_map_size=-1;
  fileSave.get(next_char);
  if ( next_char != '\n' )
  {
    fileSave >> contact_map_size;
    if (contact_map_size)
    {
      for(int j=0; j<contact_map_size; j++)
      {
        // Read contact memory map here
        int id0, id1, id2;
        bool isActive ;
        double x, y, z ;
        Vector3 tangent, prev_normal, cumulSpringTorque ;
        fileSave >> id0 >> id1 >> id2 >> isActive >> x >> y >> z;
        tangent = Vector3(x,y,z);
        fileSave >> x >> y >> z ;
        prev_normal = Vector3(x,y,z);
        fileSave >> x >> y >> z ;
        cumulSpringTorque = Vector3(x,y,z);
        this->copyContactInMap(std::make_tuple(id0,id1,id2), isActive,
          tangent, prev_normal, cumulSpringTorque) ;
      }
    }
  }
}

void Component::readContactMap_binary( istream &fileSave)
{
  int mapSize;
  fileSave.read( reinterpret_cast<char*>( &mapSize ), sizeof(int) );
  for(int k=0; k<mapSize; k++)
  {
    // Read contact memory map here
    int id0, id1, id2, isActive ;
    Vector3 tangent, prev_normal, cumulSpringTorque ;
    fileSave.read( reinterpret_cast<char*>( &id0 ), sizeof(int) );
    fileSave.read( reinterpret_cast<char*>( &id1 ), sizeof(int) );
    fileSave.read( reinterpret_cast<char*>( &id2 ), sizeof(int) );
    fileSave.read( reinterpret_cast<char*>( &isActive ), sizeof(int) );
    tangent.readGroup3_binary( fileSave );
    prev_normal.readGroup3_binary( fileSave );
    cumulSpringTorque.readGroup3_binary( fileSave );
    this->copyContactInMap(std::make_tuple(id0,id1,id2), (bool)isActive,
      tangent, prev_normal, cumulSpringTorque) ;
  }
}



// ----------------------------------------------------------------------------
// Sets the component ID number
void Component::setID( int const& id_ )
{
  m_id = id_;
}




// ----------------------------------------------------------------------------
// Resets the number of created components to nb_
void Component::setNbCreatedComponents( const int &nb_ )
{
  m_nb = nb_;
}




// ----------------------------------------------------------------------------
// Returns the number of created components  */
int Component::getNbCreatedComponents()
{
  return ( m_nb );
}




// ----------------------------------------------------------------------------
// Returns whether the component is a composite particle
bool Component::isCompositeParticle() const
{
  return ( false );
}




// ----------------------------------------------------------------------------
// Returns whether the component is a composite obstacle
bool Component::isCompositeObstacle() const
{
  return ( false );
}




// ----------------------------------------------------------------------------
// Returns whether the component is an obstacle ? (use the fact that obstacle
// have a zero mass by convention)
bool Component::isObstacle() const
{
  return ( m_mass == 0. );
}




// ----------------------------------------------------------------------------
// Returns whether the component is an elementary particle
bool Component::isElementaryParticle() const
{
  return ( false );
}




// ----------------------------------------------------------------------------
// Returns the particle class
int Component::getGeometricType() const
{
  return ( -100 );
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the component
bool Component::isIn( Point3 const& pt ) const
{
  return ( m_geoRBWC->isIn( pt ) );
}
