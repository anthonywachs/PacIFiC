#include "Cell.hh"
#include "AllComponents.hh"
#include <assert.h>
#include <math.h>
#include <string>
#include <algorithm>
using namespace std;


int Cell::m_nbi = 0;
int Cell::m_nbj = 0;
int Cell::m_nbk = 0;
double Cell::m_edge_X = 0.;
double Cell::m_edge_Y = 0.;
double Cell::m_edge_Z = 0.;
double Cell::m_LC_local_xmax = 0.;
double Cell::m_LC_local_ymax = 0.;
double Cell::m_LC_local_zmax = 0.;
Point3 Cell::m_LC_local_origin;




// ----------------------------------------------------------------------------
// Default constructor
Cell::Cell() 
  : m_number( -1 ) 
  , m_tag( 0 ) 
  , m_GeoPosCell( GEOPOS_NONE )
{
  m_cel[X] = -1;
  m_cel[Y] = -1;
  m_cel[Z] = -1;
}




// ----------------------------------------------------------------------------
// Constructor with input parameters
Cell::Cell( int id1, int x, int y, int z, Point3 const& OL,
	double arete_X, double arete_Y, double arete_Z, 
	double xmax_, double ymax_, double zmax_, 
	int tag_, GeoPosition geoloc_ ) 
  : m_number( id1 ) 
  , m_tag( tag_ ) 
  , m_GeoPosCell( geoloc_ )
{
  m_cel[X] = x;
  m_cel[Y] = y;
  m_cel[Z] = z;
  Cell::m_edge_X = arete_X;
  Cell::m_edge_Y = arete_Y;  
  Cell::m_edge_Z = arete_Z; 
  Cell::m_LC_local_origin = OL;
  Cell::m_LC_local_xmax = xmax_; 
  Cell::m_LC_local_ymax = ymax_;  
  Cell::m_LC_local_zmax = zmax_; 
  setCentre();
}




// ----------------------------------------------------------------------------
// Destructor
Cell::~Cell()
{}



// ----------------------------------------------------------------------------
// Computes the cell center coordinates
void Cell::setCentre()
{
  m_centre[X] = Cell::m_LC_local_origin[X] + m_cel[X] * m_edge_X 
  	+ m_edge_X / 2.;
  m_centre[Y] = Cell::m_LC_local_origin[Y] + m_cel[Y] * m_edge_Y 
  	+ m_edge_Y / 2.;
  m_centre[Z] = Cell::m_LC_local_origin[Z] + m_cel[Z] * m_edge_Z 
  	+ m_edge_Z / 2.;
}


  
  
// ----------------------------------------------------------------------------
// Adds a particle to the list of particles in the cell. No geometric test is 
// performed, the particle pointer is simply added to the list
void Cell::add( Particle* particle_ )
{
  m_particles.push_front( particle_ );
}  




// ----------------------------------------------------------------------------
// Adds a cell to the list of neighboring cells for contact detection
void Cell::addNeighboringCellContact( Cell* neighbor )
{
  m_neighborsContact.push_front( neighbor );
}




// ----------------------------------------------------------------------------
// Ajout de la cellule voisine au voisinage complet
void Cell::addNeighboringCell( Cell* neighbor  )
{
  m_allNeighbors.push_front( neighbor );
}




// ----------------------------------------------------------------------------
// Adds an obstacle to the list of obstacles in the vicinity of the cell. 
// No geometric test is performed, the obstacle pointer is simply added to the 
// list
void Cell::addObstacle( SimpleObstacle* obstacle_ )
{
  bool alreadyInserted = false;
  list<SimpleObstacle*>::iterator il = m_obstacles.begin();
  
  while( il != m_obstacles.end() && !alreadyInserted )
  {
    if ( *il == obstacle_ ) alreadyInserted = true;
    else il++;    
  }
  
  if ( !alreadyInserted ) m_obstacles.push_back( obstacle_ );
}




// ----------------------------------------------------------------------------
// Returns the number of obstacles in the vicinity of the cell
int Cell::numberOfObstacles() const
{
  return ( int( m_obstacles.size() ) );
}




// ----------------------------------------------------------------------------
// Clears the list of particles that belong to the cell
void Cell::clearParticles()
{
  m_particles.clear();
}




// ----------------------------------------------------------------------------
// Returns whether a particle belongs to the list of particles in the cell
bool Cell::contains( Particle* particle_ )
{
  return ( find( m_particles.begin(), m_particles.end(), particle_ ) 
    != m_particles.end() );
}




// ----------------------------------------------------------------------------
// Returns the cell ijk indices that contains a point
void Cell::GetCell( Point3 const& position, int* id )
{	
  // We use floor instead of int. If x < m_LC_local_origin[X]
  // then int returns 0 while floor returns -1, and -1 is the correct value
  // as the position is out of the LinkedCell
  // Rem: floor returns a double that we recast into an int
  id[X] = int( floor( ( position[X] - Cell::m_LC_local_origin[X] ) 
  	/ Cell::m_edge_X ) );
  id[Y] = int( floor( ( position[Y] - Cell::m_LC_local_origin[Y] ) 
  	/ Cell::m_edge_Y ) );
  id[Z] = int( floor( ( position[Z] - Cell::m_LC_local_origin[Z] ) 
  	/ Cell::m_edge_Z ) );	  
}




// ----------------------------------------------------------------------------
// Returns the volume of all particles whose center of mass is in the cell
double Cell::getVolumeParticles()
{
  double volume = 0.0;
  list<Particle*>::iterator particle = m_particles.begin();
  for ( ; particle!=m_particles.end(); particle++) 
    volume += (*particle)->getVolume();

  return(volume);
}




// ----------------------------------------------------------------------------
// Returns the cell center
Point3 const* Cell::getCentre() const
{
  return ( &m_centre );
}




// ----------------------------------------------------------------------------
// Returns whether a particle is in contact with another component
// in the vicinity of the cell
bool Cell::isContact( Particle const* particle_ ) const
{
  bool contact = false;
  
  // Contact detection with other particles
  list<Particle*>::const_iterator neighbor = m_particles.begin();
  for ( ; neighbor!=m_particles.end() && !contact; neighbor++) 
    if ( *neighbor != particle_ )
    {
      if ( particle_->isCompositeParticle() )
        contact = particle_->isContact( *neighbor ); 
      else 
        contact = (*neighbor)->isContact( particle_ );
    }         

  // Contact detection with obstacles
  list<SimpleObstacle*>::const_iterator obs = m_obstacles.begin();
  for ( ; obs!=m_obstacles.end() && !contact; obs++)
    if ( particle_->isCompositeParticle() )
      contact = particle_->isContact( *obs );
    else
      contact = (*obs)->isContact( particle_ ); 
  
  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether a particle is in contact with another component
// in the vicinity of the cell. The contact detection is performed with the
// crust width
bool Cell::isContactWithCrust( Particle const* particle_ ) const
{
  bool contact = false;
  
  // Contact detection with other particles
  list<Particle*>::const_iterator neighbor = m_particles.begin();
  for ( ; neighbor!=m_particles.end() && !contact; neighbor++) 
    if ( *neighbor != particle_ ) 
    {
      if ( particle_->isCompositeParticle() )
        contact = particle_->isContactWithCrust( *neighbor ); 
      else 
        contact = (*neighbor)->isContactWithCrust( particle_ );
    }

  // Contact detection with obstacles
  list<SimpleObstacle*>::const_iterator obs = m_obstacles.begin();
  for ( ; obs!=m_obstacles.end() && !contact; obs++)
    if ( particle_->isCompositeParticle() )
      contact = particle_->isContactWithCrust( *obs );
    else
      contact = (*obs)->isContactWithCrust( particle_ ); 

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether a particle is close to another component
// in the vicinity of the cell  
bool Cell::isClose( Particle const* particle_ ) const
{
  bool contact = false;
  
  // Closeness detection with other particles
  list<Particle*>::const_iterator neighbor = m_particles.begin();
  for ( ; neighbor!=m_particles.end() && !contact; neighbor++) 
    if ( *neighbor != particle_ ) 
    {
      if ( particle_->isCompositeParticle() )
        contact = particle_->isClose( *neighbor ); 
      else 
        contact = (*neighbor)->isClose( particle_ );
    }

  // Closeness detection with obstacles
  list<SimpleObstacle*>::const_iterator obs = m_obstacles.begin();
  for ( ; obs!=m_obstacles.end() && !contact; obs++)
    if ( particle_->isCompositeParticle() )
      contact = particle_->isClose( *obs );
    else
      contact = (*obs)->isClose( particle_ );
  
  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether a particle is close to another component in the vicinity of 
// the cell. The closeness detection is performed with the crust width
bool Cell::isCloseWithCrust( Particle const* particle_ ) const
{
  bool contact = false;
  
  // Closeness detection with other particles
  list<Particle*>::const_iterator neighbor = m_particles.begin();
  for ( ; neighbor!=m_particles.end() && !contact; neighbor++) 
    if ( *neighbor != particle_ ) 
    {
      if ( particle_->isCompositeParticle() )
        contact = particle_->isCloseWithCrust( *neighbor ); 
      else 
        contact = (*neighbor)->isCloseWithCrust( particle_ );
    }

  // Closeness detection with obstacles
  list<SimpleObstacle*>::const_iterator obs = m_obstacles.begin();
  for ( ; obs!=m_obstacles.end() && !contact; obs++)
    if ( particle_->isCompositeParticle() )
      contact = particle_->isCloseWithCrust( *obs );
    else
      contact = (*obs)->isCloseWithCrust( particle_ );  
  
  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether the cell does not contain any particle
bool Cell::isEmpty() const
{
  return ( m_particles.empty() );
}




// ----------------------------------------------------------------------------
// Returns a list of particles that do not belong to the cell anymore
void Cell::LinkUpdate( list<Particle*>& particlesExit )
{
  int id[3];
  Point3 centre;
  list<Particle*>::iterator particle;
  for (particle=m_particles.begin(); particle!=m_particles.end(); ) 
  {
    centre = *(*particle)->getPosition();
    Cell::GetCell( centre, id );

    bool present = id[X] == m_cel[X] && id[Y] == m_cel[Y] && id[Z] == m_cel[Z];
    if ( present ) particle++;
    else 
    {
      particlesExit.push_back( *particle );
      particle = m_particles.erase( particle );
    } 
  }
}




// ----------------------------------------------------------------------------
// Sets the number of cells of the linked-cell in each direction
void Cell::setNbCellsPerDirection( int nbX, int nbY, int nbZ )
{
  Cell::m_nbi = nbX;
  Cell::m_nbj = nbY;
  Cell::m_nbk = nbZ;
}




// ----------------------------------------------------------------------------
// Removes a particle from the list of particles in the cell
void Cell::remove( Particle* particle_ )
{
  list<Particle*>::iterator element;
  element = find( m_particles.begin(), m_particles.end(), particle_ );

  assert( element != m_particles.end() );

  m_particles.erase( element );
}




// ----------------------------------------------------------------------------
// Removes an obstacle from the list of obstacles in the vicinity of the cell
void Cell::remove( SimpleObstacle* obs )
{
  removeObstacleFromList( m_obstacles, obs );
}




// ----------------------------------------------------------------------------
// Returns the cell index in direction dir in the linked-cell
int Cell::operator [] ( int dir ) const
{
  assert( -1 < dir && dir < 3 );
  return ( m_cel[dir] );
}


  

// ----------------------------------------------------------------------------
// Comparison operator based on addresses
bool Cell::operator == ( Cell const& cell_ ) const
{
  return ( this == &cell_ );
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& f, Cell const& C )
{
  f << "Number = " << C.m_number << endl;
  f << "Indices = (" << C.m_cel[X] << "," << C.m_cel[Y] << "," << 
  	C.m_cel[Z] << ")" << endl;
  f << "Cell size en X x Y x Z = " << C.m_edge_X << " x " 
  	<< C.m_edge_Y << " x " << C.m_edge_Z << endl;
  f << "Local origin of linked-cell grid = " << C.m_LC_local_origin;
  f << "Tag = " << C.m_tag << endl; 
  f << "Geolocalisation in the linked-cell grid = " << C.m_GeoPosCell <<
  	" " << Cell::getGeoPositionName( C.m_GeoPosCell );   
  if ( C.m_obstacles.size() )
  {
    f << endl << "Obstacles =";
    for (list<SimpleObstacle*>::const_iterator il=C.m_obstacles.begin();
    	il!=C.m_obstacles.end();il++) f << " " << (*il)->getName();
  }
       
  return ( f );
}




// ----------------------------------------------------------------------------
// Returns the complete neighborhood of the cell
list<Cell*> const* Cell::getCompleteNeighborhood() const
{
  return ( &m_allNeighbors );
}  




// ----------------------------------------------------------------------------
// Returns a pointer to the list of particles in the cell
list<Particle*>* Cell::getParticles()
{
  return ( &m_particles );
}  




// ----------------------------------------------------------------------------
// Returns the cell geographic position
GeoPosition Cell::getGeoPosition() const
{
  return ( m_GeoPosCell );
} 




// ----------------------------------------------------------------------------
// Returns the geographic position name
string Cell::getGeoPositionName( GeoPosition geoloc_ )
{
  return ( Cell::getGeoPositionName_generic( geoloc_ ) );
}




// ----------------------------------------------------------------------------
// Returns the geographic position name
string Cell::getGeoPositionName( int geoloc_ )
{
  return ( Cell::getGeoPositionName_generic( geoloc_ ) );
}




// ----------------------------------------------------------------------------
// Returns the cell tag
int Cell::getTag() const
{
  return ( m_tag );
}




// ----------------------------------------------------------------------------
// Returns the cell number
int Cell::getID() const
{
  return ( m_number );
}




// ----------------------------------------------------------------------------
// Returns the geographic position name
string Cell::getGeoPositionName_generic( int geoloc_ )
{
  string name;
  switch( geoloc_ )
  {
    case GEOPOS_NORTH:
      name = "NORTH";
      break;
    case GEOPOS_NORTH_EAST:
      name = "NORTH_EAST";
      break;      
    case GEOPOS_NORTH_WEST:
      name = "NORTH_WEST";
      break;  
    case GEOPOS_NORTH_FRONT:
      name = "NORTH_FRONT";
      break;
    case GEOPOS_NORTH_BEHIND:
      name = "NORTH_BEHIND";
      break;      
    case GEOPOS_NORTH_EAST_FRONT:
      name = "NORTH_EAST_FRONT";
      break;      
    case GEOPOS_NORTH_EAST_BEHIND:
      name = "NORTH_EAST_BEHIND";
      break;      
    case GEOPOS_NORTH_WEST_FRONT:
      name = "NORTH_WEST_FRONT";
      break;      
    case GEOPOS_NORTH_WEST_BEHIND:
      name = "NORTH_WEST_BEHIND";
      break;      
    case GEOPOS_SOUTH:
      name = "SOUTH";
      break;      
    case GEOPOS_SOUTH_EAST:
      name = "SOUTH_EAST";
      break;      
    case GEOPOS_SOUTH_WEST:
      name = "SOUTH_WEST";
      break;      
    case GEOPOS_SOUTH_FRONT:
      name = "SOUTH_FRONT";
      break;      
    case GEOPOS_SOUTH_BEHIND:
      name = "SOUTH_BEHIND";
      break;
    case GEOPOS_SOUTH_EAST_FRONT:
      name = "SOUTH_EAST_FRONT";
      break;      
    case GEOPOS_SOUTH_EAST_BEHIND:
      name = "SOUTH_EAST_BEHIND";
      break;  
    case GEOPOS_SOUTH_WEST_FRONT:
      name = "SOUTH_WEST_FRONT";
      break;
    case GEOPOS_SOUTH_WEST_BEHIND:
      name = "SOUTH_WEST_BEHIND";
      break;      
    case GEOPOS_EAST:
      name = "EAST";
      break;      
    case GEOPOS_WEST:
      name = "WEST";
      break;      
    case GEOPOS_EAST_FRONT:
      name = "EAST_FRONT";
      break;      
    case GEOPOS_EAST_BEHIND:
      name = "EAST_BEHIND";
      break;      
    case GEOPOS_WEST_FRONT:
      name = "WEST_FRONT";
      break;      
    case GEOPOS_WEST_BEHIND:
      name = "WEST_BEHIND";
      break;      
    case GEOPOS_FRONT:
      name = "FRONT";
      break;      
    case GEOPOS_BEHIND:
      name = "BEHIND";
      break;      
    case GEOPOS_NONE:
      name = "NONE";
      break; 
    case GEOPOS_EASTWEST:
      name = "EAST_WEST";
      break;      
    case GEOPOS_NORTHSOUTH:
      name = "NORTH_SOUTH";
      break;      
    case GEOPOS_EASTWESTNORTHSOUTH:
      name = "EAST_WEST_NORTH_SOUTH";
      break;              
  }      
      
  return ( name );
}




// ----------------------------------------------------------------------------
// Returns whether a point belongs to any particle that belongs to the cell
bool Cell::isInParticle( Point3 const& position ) const
{
  bool isIn = false;
  
  for ( list<Particle*>::const_iterator particle = m_particles.begin(); 
  	particle!=m_particles.end() && !isIn; particle++) 
   isIn = (*particle)->isIn( position );    
      
  return ( isIn );
}




// ----------------------------------------------------------------------------
// Returns whether a point belongs to a specific particle that 
// belongs to the cell
bool Cell::isInParticle( Point3 const& position, Particle const* particle_ ) 
	const
{
  bool isIn = false;
  
  for ( list<Particle*>::const_iterator particle = m_particles.begin(); 
  	particle!=m_particles.end() && !isIn; particle++)
    if ( *particle == particle_ ) 
      isIn = (*particle)->isIn( position );    
      
  return ( isIn );
}
