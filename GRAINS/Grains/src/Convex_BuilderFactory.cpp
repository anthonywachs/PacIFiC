#include "Convex_BuilderFactory.H"
#include "Box.H"
#include "Cylinder.H"
#include "Disque.H"
#include "Polygon.H"
#include "Polyhedron.H"
#include "Sphere.H"
#include "PointC.H"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction du convexe specifie a partir d'un noeud XML
Convex* Convex_BuilderFactory::create( DOMNode* root )
{
  assert(root != NULL);

  Convex* convex = NULL;

  DOMNode* element = ReaderXML::getNodeNext(root);
  string   type    = ReaderXML::getNodeName(element);

  if ( type == "Box" ) convex = new Box( element );
  else if ( type == "Cylindre" ) convex = new Cylinder( element );
  else if ( type == "Disque" ) convex = new Disque( element );
  else if ( type == "Polygon" ) convex = Polygon::create( element );
  else if ( type == "Polyhedron" ) convex = Polyhedron::create( element );
  else if ( type == "Sphere" ) convex = new Sphere( element );

  assert(convex != NULL);

  return convex;
}




// ----------------------------------------------------------------------------
// Creation d'un convexe a partir du type indique & lecture de ses donnees.
// D. RAKOTONIRINA - Oct. 2014 - Modification
Convex* Convex_BuilderFactory::create( string &type, istream &fileIn )
{
  Convex *convex = NULL;

  if ( type == "*Box" ) convex = new Box( fileIn );    
  else if ( type == "*Cylindre" ) convex = new Cylinder( fileIn );
  else if ( type == "*Disque" ) convex = new Disque( fileIn );
  else if ( type == "*Polygon" ) convex = Polygon::create( fileIn );  
  else if ( type == "*Polyhedron" ) convex = Polyhedron::create( fileIn );
  else if ( type == "*Sphere" ) convex = new Sphere( fileIn );
  else if ( type == "*CompParticule" ) convex = new PointC();
  else
  {
    cout << "Forme a lire avec type non valide : " << type.c_str() << endl;
    exit(1);
  }

  return convex;
}
