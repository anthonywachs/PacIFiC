#include "Grains_BuilderFactory.H"
#include "Cinematique_BuilderFactory.H"
#include "Convex.H"
#include "LeapFrog_2D.H"
#include "LeapFrog_3D.hh"
#include "LeapFrog_Sphere.hh"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation de la cinematique de la particule
CineParticule* Cinematique_BuilderFactory::create( const Convex* convexe )
{
  CineParticule *cinematique = NULL;
  
  EAPPLI context = Grains_BuilderFactory::getContext();

  switch (context)
  {
    case DIM_2:
      cinematique = new LeapFrog_2D();
      break;
    case DIM_3:
      if ( convexe->getConvexType() == SPHERE ) 
        cinematique = new LeapFrog_Sphere();
      else cinematique = new LeapFrog_3D();
      break;
    case UNDEFINED:
      cinematique = NULL;
      break;  
  }

  return cinematique;
}




// ----------------------------------------------------------------------------
// Construction de la cinematique a partir de son enregistrement
CineParticule* Cinematique_BuilderFactory::read( istream &fileIn,
	const Convex* convexe )
{
  CineParticule *cinematique = NULL;

  string cle;
  fileIn >> cle;
  if ( cle == "*LeapFrog_3D" ) 
  {
    if ( convexe->getConvexType() == SPHERE ) 
      cinematique = new LeapFrog_Sphere();
    else cinematique = new LeapFrog_3D();    
  } 
  else if ( cle == "*LeapFrog_2D" ) cinematique = new LeapFrog_2D();
  
  fileIn >> *cinematique;

  return cinematique;
}
