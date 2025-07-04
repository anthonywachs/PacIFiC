#include "ConvexBuilderFactory.hh"
#include "Box.hh"
#include "Cylinder.hh"
#include "Disc.hh"
#include "Polygon.hh"
#include "Polyhedron.hh"
#include "Sphere.hh"
#include "Cone.hh"
#include "PointC.hh"
#include "Superquadric.hh"
#include "Rectangle.hh"
#include "TrapezoidalPrism.hh"
#include "SpheroCylinder.hh"
#include "SpheroCylindricalPrism.hh"
#include "TruncatedCone.hh"


// ----------------------------------------------------------------------------
// Construct a convex with an XML node as an input parameter
Convex* ConvexBuilderFactory::create( DOMNode* root )
{
  assert( root != NULL );

  Convex* convex = NULL;

  DOMNode* element = ReaderXML::getNodeNext( root );
  string   type    = ReaderXML::getNodeName( element );

  if ( type == "Box" ) convex = new Box( element );
  else if ( type == "Cylinder" ) convex = new Cylinder( element );
  else if ( type == "Disc" ) convex = new Disc( element );
  else if ( type == "Polygon" ) convex = Polygon::create( element );
  else if ( type == "Polyhedron" ) convex = Polyhedron::create( element );
  else if ( type == "Sphere" ) convex = new Sphere( element );
  else if ( type == "Cone" ) convex = new Cone( element );
  else if ( type == "Superquadric" ) convex = new Superquadric( element );
  else if ( type == "Rectangle" ) convex = new Rectangle( element );
  else if ( type == "TrapezoidalPrism" ) 
    convex = new TrapezoidalPrism( element );  
  else if ( type == "SpheroCylinder" || type == "SpheroCyl" ) 
    convex = new SpheroCylinder( element ); 
  else if ( type == "SpheroCylindricalPrism" ) 
    convex = new SpheroCylindricalPrism( element );
  else if ( type == "TruncatedCone" ) 
    convex = new TruncatedCone( element );         
    
  assert( convex != NULL );

  return ( convex );
}




// ----------------------------------------------------------------------------
// Construct a convex with a type and a input stream as input parameters
Convex* ConvexBuilderFactory::create( string& type, istream& fileIn )
{
  Convex *convex = NULL;

  if ( type == "*Box" ) convex = new Box( fileIn );
  else if ( type == "*Cylinder" ) convex = new Cylinder( fileIn );
  else if ( type == "*Disc" ) convex = new Disc( fileIn );
  else if ( type == "*Polygon" ) convex = Polygon::create( fileIn );
  else if ( type == "*Polyhedron" ) convex = Polyhedron::create( fileIn );
  else if ( type == "*Sphere" ) convex = new Sphere( fileIn );
  else if ( type == "*Cone" ) convex = new Cone( fileIn );
  else if ( type == "*PointC" ) convex = new PointC();
  else if ( type == "*Superquadric" ) convex = new Superquadric( fileIn );
  else if ( type == "*Rectangle" ) convex = new Rectangle( fileIn );
  else if ( type == "*TrapezoidalPrism" ) 
    convex = new TrapezoidalPrism( fileIn ); 
  else if ( type == "*SpheroCylinder" || type == "*SpheroCyl" ) 
    convex = new SpheroCylinder( fileIn );
  else if ( type == "*SpheroCylindricalPrism" ) 
    convex = new SpheroCylindricalPrism( fileIn ); 
  else if ( type == "*TruncatedCone" ) 
    convex = new TruncatedCone( fileIn );                
  else
  {
    cout << "Invalid convex type : " << type.c_str() << endl;
    exit(1);
  }

  return ( convex );
}
