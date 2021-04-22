#include "Cellule.H"
#include "EnsComposant.H"
#include "Basic.H"
#include <assert.h>
#include <math.h>
#include <string>
#include <algorithm>
using namespace std;


// ----------------------------------------------------------------------------
// Initialisation de l'attribut static de longueur des aretes et de l'espace
int    Cellule::nbi    = 0;
int    Cellule::nbj    = 0;
int    Cellule::nbk    = 0;
double Cellule::m_arete_X = 0.;
double Cellule::m_arete_Y = 0.;
double Cellule::m_arete_Z = 0.;
double Cellule::m_xmax = 0.;
double Cellule::m_ymax = 0.;
double Cellule::m_zmax = 0.;
Point Cellule::m_LC_origine_locale;




//-----------------------------------------------------------------------------
// Constructeur par defaut
Cellule::Cellule() :
  m_number( -1 ), 
  m_tag( 0 ), 
  m_GeoLocCell( MPIGEO_NONE )
{
  m_cel[X] = -1;
  m_cel[Y] = -1;
  m_cel[Z] = -1;
}




//-----------------------------------------------------------------------------
// Constructeur avec initialisation
Cellule::Cellule( int id1, int x, int y, int z, 
	double arete_X, double arete_Y, double arete_Z,
	double xmax_, double ymax_, double zmax_ ) : 
  m_number( id1 ), 
  m_tag( 0 ), 
  m_GeoLocCell( MPIGEO_NONE )
{
  m_cel[X] = x;
  m_cel[Y] = y;
  m_cel[Z] = z;
  Cellule::m_arete_X = arete_X;
  Cellule::m_arete_Y = arete_Y;  
  Cellule::m_arete_Z = arete_Z;  
  Cellule::m_xmax = xmax_; 
  Cellule::m_ymax = ymax_;  
  Cellule::m_zmax = zmax_; 
  setCentre();  
}




//-----------------------------------------------------------------------------
// Constructeur avec initialisation
Cellule::Cellule( int id1, int x, int y, int z, double arete_X, double arete_Y,
	double arete_Z, int tag_, const Point &OL,
	double xmax_, double ymax_, double zmax_, MPIGeoLocalisation geoloc_ ) :
  m_number( id1 ), 
  m_tag( tag_ ), 
  m_GeoLocCell( geoloc_ )
{
  m_cel[X] = x;
  m_cel[Y] = y;
  m_cel[Z] = z;
  Cellule::m_arete_X = arete_X;
  Cellule::m_arete_Y = arete_Y;  
  Cellule::m_arete_Z = arete_Z; 
  Cellule::m_LC_origine_locale = OL;
  Cellule::m_xmax = xmax_; 
  Cellule::m_ymax = ymax_;  
  Cellule::m_zmax = zmax_; 
  setCentre();
}




//-----------------------------------------------------------------------------
// Destructeur
Cellule::~Cellule()
{}



// ----------------------------------------------------------------------------
// Calcule le centre de la cellule
void Cellule::setCentre()
{
  m_centre[X] = Cellule::m_LC_origine_locale[X] + m_cel[X] * m_arete_X 
  	+ m_arete_X / 2.;
  m_centre[Y] = Cellule::m_LC_origine_locale[Y] + m_cel[Y] * m_arete_Y 
  	+ m_arete_Y / 2.;
  m_centre[Z] = Cellule::m_LC_origine_locale[Z] + m_cel[Z] * m_arete_Z 
  	+ m_arete_Z / 2.;
}


  
  
// ----------------------------------------------------------------------------
// Ajout de la particule a la cellule.
void Cellule::add( Particule* particule )
{
  m_particules.push_front(particule);
}  




// ----------------------------------------------------------------------------
// Ajout de la cellule voisine au voisinage pour les contacts
void Cellule::addVoisineContact( Cellule *voisine )
{
  m_voisinesContact.push_front(voisine);
}




// ----------------------------------------------------------------------------
// Ajout de la cellule voisine au voisinage complet
void Cellule::addVoisine( Cellule *voisine )
{
  m_allVoisines.push_front( voisine );
}




// ----------------------------------------------------------------------------
// Ajout d'un obstacle dans le voisinage de la cellule 
void Cellule::addObstacle( MonObstacle* obstacle_ )
{
  bool alreadyInserted = false;
  list<MonObstacle*>::iterator il = m_obstacles.begin();
  
  while( il!=m_obstacles.end() && !alreadyInserted )
  {
    if (*il == obstacle_) alreadyInserted = true;
    else il++;    
  }
  
  if ( !alreadyInserted ) m_obstacles.push_back(obstacle_);
}




// ----------------------------------------------------------------------------
// Renvoie le nombre d'obstacles dans le voisinage de la cellule
int Cellule::nombreObstacles() const
{
  return int(m_obstacles.size());
}




// ----------------------------------------------------------------------------
// Suppression des particules dans la cellule
void Cellule::clearParticules()
{
  m_particules.clear();
}




// ----------------------------------------------------------------------------
// La cellule contient-elle le Particule indique ? 
bool Cellule::Contient( Particule* particule )
{
  return find( m_particules.begin(), m_particules.end(), particule ) 
    != m_particules.end();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Indices de la cellule associee a la position
void Cellule::GetCellule(const Point &position, int* id)
{
//   Scalar x = position[X];
//   Scalar y = position[Y];
//   Scalar z = position[Z];
// 
//   // si la particule est sur la face externe, on force son appartenance
//   // Utile mais couteux, voir plus tard si on peut faire mieux
//   if ( x == m_xmax ) x -= EPS;
//   if ( y == m_ymax ) y -= EPS;
//   if ( z == m_zmax ) z -= EPS;	
// 
// //   id[X] = int((x-Cellule::m_LC_origine_locale[X]) / Cellule::m_arete_X);
// //   id[Y] = int((y-Cellule::m_LC_origine_locale[Y]) / Cellule::m_arete_Y);
// //   id[Z] = int((z-Cellule::m_LC_origine_locale[Z]) / Cellule::m_arete_Z);
//   
//   // Utilisation de floor plutot que int car si x < m_LC_origine_locale[X]
//   // int renvoie 0 alors que floor renvoie -1, et -1 est la valeur attendue car
//   // la particule est hors du LinkedCell
//   // Rem: floor renvoie un type double, qu'on re-cast en int 
//   id[X] = int( floor( ( x - Cellule::m_LC_origine_locale[X] ) 
//   	/ Cellule::m_arete_X ) );
//   id[Y] = int( floor( ( y - Cellule::m_LC_origine_locale[Y] ) 
//   	/ Cellule::m_arete_Y ) );
//   id[Z] = int( floor( ( z - Cellule::m_LC_origine_locale[Z] ) 
//   	/ Cellule::m_arete_Z ) );
	
  // Utilisation de floor plutot que int car si x < m_LC_origine_locale[X]
  // int renvoie 0 alors que floor renvoie -1, et -1 est la valeur attendue car
  // la particule est hors du LinkedCell
  // Rem: floor renvoie un type double, qu'on re-cast en int 
  id[X] = int( floor( ( position[X] - Cellule::m_LC_origine_locale[X] ) 
  	/ Cellule::m_arete_X ) );
  id[Y] = int( floor( ( position[Y] - Cellule::m_LC_origine_locale[Y] ) 
  	/ Cellule::m_arete_Y ) );
  id[Z] = int( floor( ( position[Z] - Cellule::m_LC_origine_locale[Z] ) 
  	/ Cellule::m_arete_Z ) );	  
}




// ----------------------------------------------------------------------------
// Volume global des particules dans la cellule
double Cellule::getVolumeParticules()
{
  double volume = 0.0;
  list<Particule*>::iterator particule = m_particules.begin();
  for ( ; particule!=m_particules.end(); particule++) 
    volume += (*particule)->getVolume();

  return(volume);
}




// ----------------------------------------------------------------------------
// Position du point de gravite
Point const* Cellule::getCentre() const
{
  return &m_centre;
}




// ----------------------------------------------------------------------------
// La particule est elle en contact avec une particule de la cellule ou un
// obstacle lie à la cellule 
bool Cellule::isContact( const Particule* particule ) const
{
  bool contact = false;
  
  // Contact avec les particules voisines
  list<Particule*>::const_iterator voisine = m_particules.begin();
  for ( ; voisine!=m_particules.end() && !contact; voisine++) 
    if ( *voisine != particule ) 
      contact = particule->isContact( *voisine ); 

  // Contact avec les obstacles
  list<MonObstacle*>::const_iterator obs = m_obstacles.begin();
  for ( ; obs!=m_obstacles.end() && !contact; obs++)
    if ( (*obs)->materiau() != "periode" ) 
      contact = particule->isContact( *obs ); 
  
  return (contact);
}




// ----------------------------------------------------------------------------
// La particule est elle en contact avec une particule de la cellule ou un
// obstacle lie à la cellule, variante VdW 
bool Cellule::isContactVdW( const Particule* particule ) const
{
  bool contact = false;
  
  // Contact avec les particules voisines
  list<Particule*>::const_iterator voisine = m_particules.begin();
  for ( ; voisine!=m_particules.end() && !contact; voisine++) 
    if ( *voisine != particule ) 
      contact = particule->isContactVdW( *voisine ); 

  // Contact avec les obstacles
  list<MonObstacle*>::const_iterator obs = m_obstacles.begin();
  for ( ; obs!=m_obstacles.end() && !contact; obs++)
    if ( (*obs)->materiau() != "periode" ) 
      contact = particule->isContactVdW( *obs ); 

  return (contact);
}




// ----------------------------------------------------------------------------
// La particule est elle en contact avec un obstacle periodique
bool Cellule::isContactVdW_ObstaclePeriodique( const Particule* particule ) 
	const
{
  bool contact = false;

  // Contact avec les obstacles
  list<MonObstacle*>::const_iterator obs = m_obstacles.begin();
  for ( ; obs!=m_obstacles.end() && !contact; obs++)
    if ( (*obs)->materiau() == "periode" ) 
      contact = particule->isContactVdW( *obs ); 
  
  return (contact);
} 




// ----------------------------------------------------------------------------
// La particule est elle proche d'une particule de la cellule ou d'un
// obstacle lie à la cellule 
bool Cellule::isProche( const Particule* particule ) const
{
  bool contact = false;
  
  // Contact avec les particules voisines
  list<Particule*>::const_iterator voisine = m_particules.begin();
  for ( ; voisine!=m_particules.end() && !contact; voisine++) 
    if ( *voisine != particule ) 
      contact = particule->isProche( *voisine ); 

  // Contact avec les obstacles
  list<MonObstacle*>::const_iterator obs = m_obstacles.begin();
  for ( ; obs!=m_obstacles.end() && !contact; obs++)
    if ( (*obs)->materiau() != "periode" ) 
      contact = particule->isProche( *obs );   
  
  return (contact);
}




// ----------------------------------------------------------------------------
// La particule est elle proche d'une particule de la cellule ou d'un
// obstacle lie à la cellule, variante VdW  
bool Cellule::isProcheVdW( const Particule* particule ) const
{
  bool contact = false;
  
  // Contact avec les particules voisines
  list<Particule*>::const_iterator voisine = m_particules.begin();
  for ( ; voisine!=m_particules.end() && !contact; voisine++) 
    if ( *voisine != particule ) 
      contact = particule->isProcheVdW( *voisine ); 

  // Contact avec les obstacles
  list<MonObstacle*>::const_iterator obs = m_obstacles.begin();
  for ( ; obs!=m_obstacles.end() && !contact; obs++)
    if ( (*obs)->materiau() != "periode" ) 
      contact = particule->isProcheVdW( *obs );   
  
  return (contact);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// La cellule est elle vide (pas de Particule)
bool Cellule::isEmpty() const
{
  return m_particules.empty();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Mise a jour des particules dans la cellule.
// Les particules hors cellule sont renvoyees pour traitement.
void Cellule::LinkUpdate( list<Particule*> &particulesExit )
{
  int id[3];
  Point centre;
  list<Particule*>::iterator particule;
  for (particule=m_particules.begin(); particule!=m_particules.end(); ) 
  {
    centre = *(*particule)->getPosition();
    Cellule::GetCellule(centre, id);

    bool present = id[X] == m_cel[X] && id[Y] == m_cel[Y] && id[Z] == m_cel[Z];
    if ( present ) particule++;
    else 
    {
      particulesExit.push_back(*particule);
      particule = m_particules.erase( particule );
    } 
  }
}




// ----------------------------------------------------------------------------
// Existe t'il un contact entre le composant et les composants de la cellule.
// Le composant ne se regarde pas lui-meme.
// G.FERRER - Octo.2000 - Creation
bool Cellule::PbProximite( const Composant* composant ) const
{
  bool contact = false;
  list<Particule*>::const_iterator voisine = m_particules.begin();
  for ( ; voisine!=m_particules.end() && !contact; voisine++) 
    contact = composant->isContactVdW( *voisine );

  return contact;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definition de l'espace des cellules
void Cellule::SetEspace( int nbX, int nbY, int nbZ )
{
  Cellule::nbi = nbX;
  Cellule::nbj = nbY;
  Cellule::nbk = nbZ;
}




// ----------------------------------------------------------------------------
// Suppression d'une particule de la liste
void Cellule::remove( Particule *particule )
{
  list<Particule*>::iterator element;
  element = find( m_particules.begin(), m_particules.end(), particule );

  assert( element != m_particules.end() );

  m_particules.erase( element );
}




// ----------------------------------------------------------------------------
// Suppression d'un obstacle de la cellule
void Cellule::remove( MonObstacle *obs )
{
  removeObstacleFromList( m_obstacles, obs );
}




// ----------------------------------------------------------------------------
// Acces a l'indice de la celllule
int Cellule::operator [] ( int xyz ) const
{
  assert( -1 < xyz && xyz < 3 );
  return m_cel[xyz];
}


  

// ----------------------------------------------------------------------------
// Comparaison entre deux cellules
// G.FERRER - Fevr - Creation
bool Cellule::operator == ( const Cellule &cellule ) const
{
  return ( this == &cellule );
}




// ----------------------------------------------------------------------------
// Operateur << 
ostream& operator << ( ostream &f, const Cellule &C )
{
  f << "Numero = " << C.m_number << endl;
  f << "Indices = (" << C.m_cel[X] << "," << C.m_cel[Y] << "," << 
  	C.m_cel[Z] << ")" << endl;
  f << "Taille de la cellule en X x Y x Z = " << C.m_arete_X << " x " 
  	<< C.m_arete_Y << " x " << C.m_arete_Z << endl;
  f << "Origine locale de la grille = " << C.m_LC_origine_locale;
  f << "Tag = " << C.m_tag << endl; 
  f << "Localisation geographique dans le LinkedCell = " << C.m_GeoLocCell <<
  	" " << Cellule::getMPIGeoLocalisationName( C.m_GeoLocCell );   
  if ( C.m_obstacles.size() )
  {
    f << endl << "Obstacles =";
    for (list<MonObstacle*>::const_iterator il=C.m_obstacles.begin();
    	il!=C.m_obstacles.end();il++) f << " " << (*il)->getName();
  }
       
  return f;
}




// ----------------------------------------------------------------------------
// Renvoi la liste des cellules du voisinage complet
// A.WACHS - Aout.2009 - Creation
const list<Cellule*>* Cellule::getVoisinageComplet() const
{
  return &m_allVoisines;
}  




// ----------------------------------------------------------------------------
// Liste des particules dans la cellule
// A.WACHS - Sept.2009 - Creation
list<Particule*>* Cellule::getParticules()
{
  return &m_particules;
}  




// ----------------------------------------------------------------------------
// Localisation geographique */
MPIGeoLocalisation Cellule::getGeoLocalisation() const
{
  return m_GeoLocCell;
} 




// ----------------------------------------------------------------------------
// Renvoie la geolocalication dans une chaine de caractere
string Cellule::getMPIGeoLocalisationName( MPIGeoLocalisation geoloc_ )
{
  return( Cellule::getMPIGeoLocalisationName_generic( geoloc_ ) );
}




// ----------------------------------------------------------------------------
// Renvoie la geolocalication dans une chaine de caractere
string Cellule::getMPIGeoLocalisationName( int geoloc_ )
{
  return( Cellule::getMPIGeoLocalisationName_generic( geoloc_ ) );
}




// ----------------------------------------------------------------------------
// Renvoie la geolocalication dans une chaine de caractere
string Cellule::getMPIGeoLocalisationName_generic( int geoloc_ )
{
  string name;
  switch( geoloc_ )
  {
    case MPIGEO_NORTH:
      name = "NORTH";
      break;
    case MPIGEO_NORTH_EAST:
      name = "NORTH_EAST";
      break;      
    case MPIGEO_NORTH_WEST:
      name = "NORTH_WEST";
      break;  
    case MPIGEO_NORTH_TOP:
      name = "NORTH_TOP";
      break;
    case MPIGEO_NORTH_BOTTOM:
      name = "NORTH_BOTTOM";
      break;      
    case MPIGEO_NORTH_EAST_TOP:
      name = "NORTH_EAST_TOP";
      break;      
    case MPIGEO_NORTH_EAST_BOTTOM:
      name = "NORTH_EAST_BOTTOM";
      break;      
    case MPIGEO_NORTH_WEST_TOP:
      name = "NORTH_WEST_TOP";
      break;      
    case MPIGEO_NORTH_WEST_BOTTOM:
      name = "NORTH_WEST_BOTTOM";
      break;      
    case MPIGEO_SOUTH:
      name = "SOUTH";
      break;      
    case MPIGEO_SOUTH_EAST:
      name = "SOUTH_EAST";
      break;      
    case MPIGEO_SOUTH_WEST:
      name = "SOUTH_WEST";
      break;      
    case MPIGEO_SOUTH_TOP:
      name = "SOUTH_TOP";
      break;      
    case MPIGEO_SOUTH_BOTTOM:
      name = "SOUTH_BOTTOM";
      break;
    case MPIGEO_SOUTH_EAST_TOP:
      name = "SOUTH_EAST_TOP";
      break;      
    case MPIGEO_SOUTH_EAST_BOTTOM:
      name = "SOUTH_EAST_BOTTOM";
      break;  
    case MPIGEO_SOUTH_WEST_TOP:
      name = "SOUTH_WEST_TOP";
      break;
    case MPIGEO_SOUTH_WEST_BOTTOM:
      name = "MPIGEO_SOUTH_WEST_BOTTOM";
      break;      
    case MPIGEO_EAST:
      name = "EAST";
      break;      
    case MPIGEO_WEST:
      name = "WEST";
      break;      
    case MPIGEO_EAST_TOP:
      name = "EAST_TOP";
      break;      
    case MPIGEO_EAST_BOTTOM:
      name = "EAST_BOTTOM";
      break;      
    case MPIGEO_WEST_TOP:
      name = "WEST_TOP";
      break;      
    case MPIGEO_WEST_BOTTOM:
      name = "WEST_BOTTOM";
      break;      
    case MPIGEO_TOP:
      name = "TOP";
      break;      
    case MPIGEO_BOTTOM:
      name = "BOTTOM";
      break;      
    case MPIGEO_NONE:
      name = "NONE";
      break;  
  }      
      
  return name;
}
