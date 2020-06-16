#include "ContactLaw.hh"
#include "LinkedCell.H"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur avec un lot de parametres
ContactLaw::ContactLaw( map<string,double>& parameters )
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructeur
ContactLaw::~ContactLaw() 
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Contact entre 2 clones periodiques ? 
bool ContactLaw::is_ClonePer_ClonePer( Composant const* p0_, 
		     Composant const* p1_ )
{
  return ( p0_->getID() == -2 && p1_->getID() == -2 );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Contact entre 1 clones periodique et une particule possedant un
// clone periodique dans la meme direction (sens oppose) ? 
bool ContactLaw::is_ClonePer_ParticuleWithClonePerSameDirection(
	Composant const* p0_, 
	Composant const* p1_,
	LinkedCell const* LC )
{
  bool is_CPPWCPS = false;

  if ( p1_->getNombreClonesPeriodiques() )
    is_CPPWCPS = p1_->hasCloneInDirection( p0_->getVecteurPeriodique(), LC );

  return is_CPPWCPS;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Les directions periodiques font elles un angle de 90 ?
bool ContactLaw::are_ClonePerDirections_Perp( Composant const* p0_, 
		     Composant const* p1_ )
{
  Vecteur const* v0 = p0_->getVecteurPeriodique();
  Vecteur const* v1 = p1_->getVecteurPeriodique();
  
  return ( fabs( (*v0) * (*v1) ) < 1.e-6 ) ;
}		         
