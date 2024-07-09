#include "GrainsExec.hh"
#include "ObstacleBuilderFactory.hh"
#include "CompositeObstacle.hh"
#include "SimpleObstacle.hh"
#include "STLObstacle.hh"
// #include "CylindricalShell.hh"
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
  // else if ( type == "CylindricalShell" )
  //   obstacle = new CylindricalShell( root );  

  return ( obstacle );
}




// ----------------------------------------------------------------------------
// Creates and reloads an obstacle from a stream
void ObstacleBuilderFactory::reload( string const& tag, Obstacle& mother, 
	istream& file )
{ 
  string name;

  if ( tag == "<Composite>" ) 
  {
    file >> name;
    Obstacle *composite = new CompositeObstacle( name );
    composite->reload( mother, file );
  } 
  else if ( tag == "<Simple>" ) 
  {
    file >> name;
    Obstacle *simple = new SimpleObstacle( name );
    simple->reload( mother, file );
  } 
}
