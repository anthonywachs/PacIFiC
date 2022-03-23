#include "GrainsBuilderFactory.hh"
#include "KinematicsBuilderFactory.hh"
#include "Convex.hh"
#include "ParticleKinematics2D.hh"
#include "ParticleKinematics3D.hh"
#include "ParticleKinematicsSphere.hh"


// ----------------------------------------------------------------------------
// Creates and returns the particle kinematics depending on the 
// space dimension and the particle rigid body shape. Note about convention: 
// the pointer to the convex is NULL for a composite particle
ParticleKinematics* KinematicsBuilderFactory::create( Convex const* convex_ )
{
  ParticleKinematics *cinematique = NULL;
  
  EAPPLI context = GrainsBuilderFactory::getContext();

  switch (context)
  {
    case DIM_2:
      cinematique = new ParticleKinematics2D();
      break;
    case DIM_3:
      if ( convex_ )
      {
        if ( convex_->getConvexType() == SPHERE ) 
          cinematique = new ParticleKinematicsSphere();
	else cinematique = new ParticleKinematics3D();
      }	
      else cinematique = new ParticleKinematics3D();
      break;
    case UNDEFINED:
      cinematique = NULL;
      break;  
  }

  return ( cinematique );
}




// ----------------------------------------------------------------------------
// Construction de la cinematique a partir de son enregistrement
ParticleKinematics* KinematicsBuilderFactory::read( istream& fileIn,
	Convex const* convex_ )
{
  ParticleKinematics *cinematique = NULL;

  string cle;
  fileIn >> cle;
  if ( cle == "*ParticleKinematics3D" ) 
  {
    if ( convex_ )
    {
      if ( convex_->getConvexType() == SPHERE ) 
        cinematique = new ParticleKinematicsSphere();
      else cinematique = new ParticleKinematics3D();
    }	
    else cinematique = new ParticleKinematics3D();   
  } 
  else if ( cle == "*ParticleKinematics2D" ) 
    cinematique = new ParticleKinematics2D();
  
  fileIn >> *cinematique;

  return ( cinematique );
}
