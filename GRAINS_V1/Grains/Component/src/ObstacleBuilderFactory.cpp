#include "GrainsExec.hh"
#include "ObstacleBuilderFactory.hh"
#include "CompositeObstacle.hh"
#include "SimpleObstacle.hh"
#include "STLObstacle.hh"
#include "CylindricalShell.hh"
#include "RoughWall.hh"
#include <string>
using namespace std;


// ----------------------------------------------------------------------------
// Creates an obstacle from an XML node
Obstacle* ObstacleBuilderFactory::create( DOMNode* root )
{
  assert(root != NULL);

  Obstacle *obstacle = NULL;

  string type = ReaderXML::getNodeName( root );

  if ( type == "Obstacle" ) 
  {
    type = ReaderXML::getNodeAttr_String( root, "Type" );

    if ( type == "Standard" ) obstacle = new SimpleObstacle( root ); 
    else if ( type == "STL" ) obstacle = new STLObstacle( root );
  }
  else if ( type == "Composite" )
    obstacle = new CompositeObstacle( root );
  else if ( type == "CylindricalShell" )
    obstacle = new CylindricalShell( root );
  else if ( type == "RoughWall" )
    obstacle = new RoughWall( root );      

  return ( obstacle );
}




// ----------------------------------------------------------------------------
// Creates and reloads an obstacle from a stream
void ObstacleBuilderFactory::reload( string const& tag, Obstacle& mother, 
	istream& file )
{ 
  string name, type;
  Obstacle *obstacle = NULL;  

  if ( tag == "<Composite>" ) 
  {
    file >> name >> type;
    if ( type == "Standard" ) obstacle = new CompositeObstacle( name );
    else if ( type == "CylindricalShell" )
      obstacle = new CylindricalShell( name );
    else if ( type == "RoughWall" )
      obstacle = new RoughWall( name );    
    obstacle->reload( mother, file );
  } 
  else if ( tag == "<Simple>" ) 
  {
    file >> name;
    obstacle = new SimpleObstacle( name );
    obstacle->reload( mother, file );
  } 
}
