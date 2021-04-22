#include "AddedMass.H"
#include "Torseur.H"

// ----------------------------------------------------------------------------
// Constructeur
// A.WACHS - Fev.2009 - Creation
AddedMass::AddedMass( Scalar rhoFluide_, Scalar simulTime_ ) :
  App(), rhoFluideInExplicitMass( rhoFluide_ ), simulTime( simulTime_ )
{}




// ----------------------------------------------------------------------------
// Destructeur
// A.WACHS - Fev.2009 - Creation
AddedMass::~AddedMass()
{}




// ----------------------------------------------------------------------------
// Description de la force de AddedMass exercee sur les particules
// A.WACHS - Fev.2009 - Creation
void AddedMass::CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules )
{
  Scalar masse, density;
  Vecteur addedmass_translational, addedmass_rotational,
  	RotationalVelocity_difference;
  list<Particule*>::const_iterator particule;

  for (particule=particules->begin(); particule!=particules->end();
       particule++)
  {
    density = (*particule)->getMasseVolumique();
    masse  = (*particule)->getMasse();

    // Translational part: force
    addedmass_translational =
    	( rhoFluideInExplicitMass / density ) * masse
    	* (*particule)->getTranslationalVelocityDifferencePreviousTime()
	/ simulTime;

    // Rotational part: moment
    addedmass_rotational = ( rhoFluideInExplicitMass / ( density * simulTime ) )
    	* (*particule)->getCinematique()->calculIdwExplicite(
	(*particule)->getRotationalVelocityDifferencePreviousTime(),
	(*particule)->getInertie() );
	
//    cout << simulTime << " " << addedmass_translational;
//    cout << (*particule)->getTranslationalVelocityDifferencePreviousTime();	

    (*particule)->addBodyForce( addedmass_translational );
    (*particule)->addMoment( addedmass_rotational );
  }
}




// ----------------------------------------------------------------------------
void AddedMass::setsimulTime( double simulTime_ )
{
  simulTime = simulTime_ ;
}
