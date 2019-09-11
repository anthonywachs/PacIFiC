#include "ObstacleAbsorbant.H"
#include "Particule.H"
#include "SaveTable.H"
#include <assert.h>


// ----------------------------------------------------------------------------
// Constructeur
ObstacleAbsorbant::ObstacleAbsorbant( const string &s ) :
  MonObstacle( s )
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
ObstacleAbsorbant::ObstacleAbsorbant( DOMNode *root ) :
  MonObstacle()
{
}




// ----------------------------------------------------------------------------
// Destructeur
ObstacleAbsorbant::~ObstacleAbsorbant()
{
}




// ----------------------------------------------------------------------------
// Determination du contact avec une particule
// L'obstacle absorbe la particule lorsque le centre de la particule
// penetre celui-ci.
void ObstacleAbsorbant::InterAction( Composant* particule,
	Scalar dt, double const& temps, LinkedCell *LC )
  throw (ErreurContact)
{
//   const FormeVdW* forme0 = particule->getForme();
// 
//   PointContact closestPoint;
//   try {
//     closestPoint = forme0->ClosestPoint(*geoFormeVdw);
//   } catch (ErreurContact &erreur) {
//     particule->setActivity(CLEARandWAIT);
//   }
// 
//   if (closestPoint.getDistance() < 0.) particule->setActivity(CLEARandWAIT);
// 
//   return NULL;
}




// ----------------------------------------------------------------------------
// Determination du contact avec une particule
void ObstacleAbsorbant::InterActionPostProcessing( Composant* voisin, 
	list<struct PointForcePostProcessing>* listOfContacts )
  throw (ErreurContact)
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reload de l'obstacle & publication dans son referent
void ObstacleAbsorbant::reload( Obstacle &obstacle, istream &file )
{
  string buffer, adresse;

  file >> adresse >> m_id;
  SaveTable::create( adresse, this );
  file >> buffer >> m_nomMateriau;
  m_geoFormeVdw = new FormeVdW( file );
  m_geoFormeVdw->readPosition( file );
  obstacle.append( this );
  file >> buffer;
  assert(buffer == "</Absorbant>");
}




// ----------------------------------------------------------------------------
// Sauvegarde de l'obstacle pour Reload
// G.FERRER - Juil.2003 - Creation
// D. RAKOTONIRINA - Oct. 2014 - Modification
void ObstacleAbsorbant::write( ostream &fileSave, Composant const* composant )
    const
{
  fileSave << "<Absorbant>\t" << m_nom << '\n'
	   << this << '\t' << m_id << '\n';
  fileSave << "*Materiau \n" << m_nomMateriau << '\n';  
  m_geoFormeVdw->writeStatique( fileSave );  
  m_geoFormeVdw->writePosition( fileSave );
  fileSave << "</Absorbant>\n";
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
int ObstacleAbsorbant::numberOfPoints_PARAVIEW() const
{
  cout << "Warning when calling ObstacleAbsorbant::numberOfPoints_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
int ObstacleAbsorbant::numberOfCells_PARAVIEW() const
{
  cout << "Warning when calling ObstacleAbsorbant::numberOfCells_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Nombre de polygones elementaires pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
list<Point> ObstacleAbsorbant::get_polygonsPts_PARAVIEW( Vecteur const*
    translation )
const
{
  cout << "Warning when calling ObstacleAbsorbant::get_polygonsPts_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
void ObstacleAbsorbant::write_polygonsPts_PARAVIEW( ostream &f, 
  	Vecteur const* translation ) const
{
  cout << "Warning when calling Obstacle::write_polygonsPts_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview  
// D. RAKOTONIRINA - Dec. 2014 - Creation
void ObstacleAbsorbant::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  cout << "Warning when calling ObstacleAbsorbant::write_polygonsStr_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}
