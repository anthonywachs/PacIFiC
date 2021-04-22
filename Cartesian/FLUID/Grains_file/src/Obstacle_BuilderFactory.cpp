#include "Grains_Exec.hh"
#include "Obstacle_BuilderFactory.H"
#include "CompObstacle.H"
#include "CylindricalBoxObstacle.H"
#include "MonObstacle.H"
#include "ObstacleAbsorbant.H"
#include "ObstaclePeriodique.hh"
#include <string>
using namespace std;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation de l'obstacle a partir du noeud XML
Obstacle* Obstacle_BuilderFactory::create( DOMNode *root )
{
  assert(root != NULL);

  Obstacle *obstacle = NULL;

  string type = ReaderXML::getNodeName( root );

  if ( type == "Obstacle" ) 
  {
    type = ReaderXML::getNodeAttr_String( root, "Type" );

    if ( type == "Standard" ) obstacle = new MonObstacle( root ); 
    else if ( type == "Absorbant" ) obstacle = new ObstacleAbsorbant( root );
  }
  else if ( type == "Periode" ) 
  {
    obstacle = Obstacle_BuilderFactory::createPeriode( root );
    Grains_Exec::m_periodique = true;
  }
  else if ( type == "Composite" )
    obstacle = new CompObstacle( root );

  else if ( type == "CylindricalBox" )
    obstacle = new CylindricalBoxObstacle( root ); 

  return obstacle;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation d'une periode d'obstacles.
Obstacle* Obstacle_BuilderFactory::createPeriode( DOMNode *root )
{
  Obstacle *obstacle = NULL;
  
  Vecteur direction(0.);
  DOMNode* vecteur = ReaderXML::getNode( root, "Vecteur" );
  direction[X] = ReaderXML::getNodeAttr_Double( vecteur, "X" );
  direction[Y] = ReaderXML::getNodeAttr_Double( vecteur, "Y" );
  direction[Z] = ReaderXML::getNodeAttr_Double( vecteur, "Z" );

  string name  = ReaderXML::getNodeAttr_String( root, "name" );

  obstacle = new CompObstacle(name.c_str());
  ObstaclePeriodique* obstacleA = new ObstaclePeriodique( root, name+"_A", 
	direction );
	
  ObstaclePeriodique* obstacleB = new ObstaclePeriodique( root, name+"_B",
	-direction );
  Point position = *obstacleB->getPosition();
  position += direction;
  obstacleB->setPosition( position );

  obstacleA->setAssocie( obstacleB );
  obstacleB->setAssocie( obstacleA );

  obstacle->append( obstacleA );
  obstacle->append( obstacleB );

  return obstacle;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reload de l'Obstacle sous l'Obstacle referent
void Obstacle_BuilderFactory::reload( const string &tag,
	Obstacle& obstacle, istream &file )
{ 
  if ( tag == "<Composite>" ) 
  {
    string nom;
    file >> nom;
    Obstacle *composite = new CompObstacle( nom );
    composite->reload( obstacle, file );
  } 
  else if ( tag == "<CylindricalBox>" ) 
  {
    string nom;
    double radius=0., length=0.;
    file >> nom;
    file >> radius;
    file >> length;
    Obstacle *CylindricalBox = new CylindricalBoxObstacle( nom,radius,length );
    CylindricalBox->reload( obstacle, file );
  } 
  else if ( tag == "<Simple>" ) 
  {
    string nom;
    file >> nom;
    Obstacle *simple = new MonObstacle( nom );
    simple->reload( obstacle, file );
  } 
  else if ( tag == "<Absorbant>" ) 
  {
    string nom;
    file >> nom;
    Obstacle *absorbant = new ObstacleAbsorbant( nom );
    absorbant->reload( obstacle, file );
  } 
  else if ( tag == "<Periodique>" ) 
  {
    string nom;
    file >> nom;
    ObstaclePeriodique* obstacle0 = new ObstaclePeriodique( nom );
    obstacle0->reload( obstacle, file );
    string tag1;
    file >> tag1;
    file >> nom;
    ObstaclePeriodique* obstacle1 = new ObstaclePeriodique( nom );
    obstacle1->reload( obstacle, file );
    obstacle0->setAssocie( obstacle1 );
    obstacle1->setAssocie( obstacle0 );
    Grains_Exec::m_periodique = true;
  }
}
