#include "SecondOrderAdamsBashforth.hh"
#include <limits>

// ----------------------------------------------------------------------------
// Constructeur par defaut
SecondOrderAdamsBashforth::SecondOrderAdamsBashforth() :
  TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Constructeur par copie.
SecondOrderAdamsBashforth::SecondOrderAdamsBashforth(
	const SecondOrderAdamsBashforth &copie) :
  TimeIntegrator(copie)
{
  m_vitesseT_nm2 = copie.m_vitesseT_nm2;
  m_vitesseR_nm2 = copie.m_vitesseR_nm2;
  m_dUdt_nm2 = copie.m_dUdt_nm2;
  m_dOmegadt_nm2 = copie.m_dOmegadt_nm2;
  first_step = copie.first_step;
}




// ----------------------------------------------------------------------------
// Destructeur
SecondOrderAdamsBashforth::~SecondOrderAdamsBashforth()
{}




// ----------------------------------------------------------------------------
// Contruction d'un clone de l'integrateur en temps
TimeIntegrator* SecondOrderAdamsBashforth::clone() const
{
  return new SecondOrderAdamsBashforth(*this);
}




// ----------------------------------------------------------------------------
// Integration du mouvement: calcul des nouvelles vitesses & position
void SecondOrderAdamsBashforth::Deplacer(Vecteur& vitesseT,
  	Vecteur const &dUdt,
	Vecteur &deplacementTranslationnel,
	Vecteur const &dOmegadt,
	Vecteur &vitesseR,
	Vecteur &vitesseRMoyenne,
	double dt)
{
  if (first_step==false)
  {
    // Vitesse et deplacement translationnels
    deplacementTranslationnel = 0.5 * ( 3.*vitesseT - m_vitesseT_nm2 ) * dt;
    m_vitesseT_nm2 = vitesseT;
    vitesseT += 0.5 * ( 3.*dUdt - m_dUdt_nm2 ) * dt;
    m_dUdt_nm2 = dUdt;

    // Vitesse et deplacement rotationnels
    vitesseRMoyenne = 0.5 * ( 3.*vitesseR  - m_vitesseR_nm2 );
    m_vitesseR_nm2 = vitesseR;
    vitesseR += 0.5 * ( 3.*dOmegadt - m_dOmegadt_nm2 ) * dt;
    m_dOmegadt_nm2 = dOmegadt;
  }
  else
  {
    // Vitesse et deplacement translationnels
    deplacementTranslationnel = dt*vitesseT;
    m_vitesseT_nm2 = vitesseT;
    vitesseT += dt*dUdt;
    m_dUdt_nm2 = dUdt;

    // Vitesse et deplacement rotationnels
    vitesseRMoyenne = dt*vitesseR;
    m_vitesseR_nm2 = vitesseR;
    vitesseR += dt*dOmegadt;
    m_dOmegadt_nm2 = dOmegadt;
    first_step=false;
  }
}


// WARNING: NEED TO IMPLEMENT COPY OF VARIABLE first_step IN BELOW FUNCTIONS
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copie de la cinematique au temps t-2dt: vitesse de translation,
// vitese de rotation, variation de QDM translationnelle, variation de QDM
// rotationnalle
void SecondOrderAdamsBashforth::copyCinematiqueNm2(double *vit,int i) const
{
  for (int j=0 ;j<3; j++) vit[i+j] = m_vitesseT_nm2[j];
  for (int j=0 ;j<3; j++) vit[i+3+j] = m_vitesseR_nm2[j];
  for (int j=0 ;j<3; j++) vit[i+6+j] = m_dUdt_nm2[j];
  for (int j=0 ;j<3; j++) vit[i+9+j] = m_dOmegadt_nm2[j];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cinematique au temps t-2dt: vitesse de translation,
// vitese de rotation, variation de QDM translationnelle, variation de QDM
// rotationnalle
void SecondOrderAdamsBashforth::setCinematiqueNm2(double const* tab)
{
  for (int j=0 ;j<3; j++) m_vitesseT_nm2[j] = tab[j];
  for (int j=0 ;j<3; j++) m_vitesseR_nm2[j] = tab[j+3];
  for (int j=0 ;j<3; j++) m_dUdt_nm2[j] = tab[j+6];
  for (int j=0 ;j<3; j++) m_dOmegadt_nm2[j] = tab[j+9];
}
