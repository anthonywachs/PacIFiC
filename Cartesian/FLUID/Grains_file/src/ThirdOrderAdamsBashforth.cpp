#include "ThirdOrderAdamsBashforth.hh"

// ----------------------------------------------------------------------------
// Constructeur par defaut
ThirdOrderAdamsBashforth::ThirdOrderAdamsBashforth() :
  TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Constructeur par copie.
ThirdOrderAdamsBashforth::ThirdOrderAdamsBashforth(
	const ThirdOrderAdamsBashforth &copie) :
  TimeIntegrator(copie)
{
  // m_vitesseT_nm2 = copie.m_vitesseT;
  // m_vitesseR_nm2 = copie.m_vitesseR_nm2;
  m_prev_dUdt = copie.m_prev_dUdt;
  m_2prev_dUdt = copie.m_2prev_dUdt;
  m_prev_v = copie.m_prev_v;
  m_2prev_v = copie.m_2prev_v;
  // m_dOmegadt_nm2 = copie.m_dOmegadt_nm2;
}




// ----------------------------------------------------------------------------
// Destructeur
ThirdOrderAdamsBashforth::~ThirdOrderAdamsBashforth()
{}




// ----------------------------------------------------------------------------
// Contruction d'un clone de l'integrateur en temps
TimeIntegrator* ThirdOrderAdamsBashforth::clone() const
{
  return new ThirdOrderAdamsBashforth(*this);
}




// ----------------------------------------------------------------------------
// Integration du mouvement: calcul des nouvelles vitesses & position
void ThirdOrderAdamsBashforth::Deplacer(Vecteur& vitesseT,
  	Vecteur const &dUdt,
	Vecteur &deplacementTranslationnel,
	Vecteur const &dOmegadt,
	Vecteur &vitesseR,
	Vecteur &vitesseRMoyenne,
	double dt)
{
  if (remaining_use_rk2==0)
  {
    cout << "Hi, using Adams-bashforth!" << endl;
    cout << "dUdt-2=" << m_2prev_dUdt << endl;
    cout << "v-2=" << m_2prev_v << endl;
    cout << "dUdt-1=" << m_prev_dUdt << endl;
    cout << "v-1=" << m_prev_v << endl;
    cout << "dUdt=" <<dUdt <<endl;
    cout << "v=" << vitesseT << endl;
    // Vitesse et deplacement translationnels
    deplacementTranslationnel = dt/12 * ( 23*vitesseT - 16*m_prev_v + 5*m_2prev_v);
    vitesseT += dt/12 * ( 23*dUdt- 16*m_prev_dUdt + 5*m_2prev_dUdt);
    prev_translational_velocity_contribution = dt/12 * ( 23*dUdt- 16*m_prev_dUdt + 5*m_2prev_dUdt);
    m_2prev_dUdt=m_prev_dUdt;
    m_2prev_v=m_prev_v;
    m_prev_dUdt = dUdt;
    m_prev_v=vitesseT-prev_translational_velocity_contribution;
  }
  else
  {
    cout << "dUdt=" << dUdt << endl;
    cout << "v=" << vitesseT << endl;
    cout << "remaining_use_rk2=" << remaining_use_rk2 << endl;
    if (remaining_use_rk2==4 || remaining_use_rk2==2)
    {
      cout << "entering whole step rk2 process" << endl;
      if (remaining_use_rk2==4)
      {
        m_2prev_v=vitesseT;
        m_2prev_dUdt=dUdt;
      }
      else
      {
        m_prev_v=vitesseT;
        m_prev_dUdt=dUdt;
      }
      deplacementTranslationnel= dt/2. * vitesseT;
      prev_translational_displacement_contribution = dt/2. * vitesseT;
      prev_translational_velocity_contribution = dt/2. * dUdt;
      vitesseT += dt/2. * dUdt;
      remaining_use_rk2--;
    }
    else
    {
      cout << "entering intermediate step rk2 process" << endl;
      deplacementTranslationnel = -prev_translational_displacement_contribution
                                  + dt*vitesseT;
      vitesseT += -prev_translational_velocity_contribution + dt*dUdt;
      remaining_use_rk2--;
    }
  }

  // Vitesse et deplacement rotationnels
  // No rotation for now
  vitesseR += Vecteur(0.,0.,0.);
  vitesseRMoyenne = Vecteur(0.,0.,0.);

}




// // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// // Copie de la cinematique au temps t-2dt: vitesse de translation,
// // vitese de rotation, variation de QDM translationnelle, variation de QDM
// // rotationnalle
// void ThirdOrderAdamsBashforth::copyCinematiqueNm2(double *vit,int i) const
// {
//   for (int j=0 ;j<3; j++) vit[i+j] = m_vitesseT_nm2[j];
//   for (int j=0 ;j<3; j++) vit[i+3+j] = m_vitesseR_nm2[j];
//   for (int j=0 ;j<3; j++) vit[i+6+j] = m_dUdt_nm2[j];
//   for (int j=0 ;j<3; j++) vit[i+9+j] = m_dOmegadt_nm2[j];
// }
//
//
//
//
// // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// // Cinematique au temps t-2dt: vitesse de translation,
// // vitese de rotation, variation de QDM translationnelle, variation de QDM
// // rotationnalle
// void ThirdOrderAdamsBashforth::setCinematiqueNm2(double const* tab)
// {
//   for (int j=0 ;j<3; j++) m_vitesseT_nm2[j] = tab[j];
//   for (int j=0 ;j<3; j++) m_vitesseR_nm2[j] = tab[j+3];
//   for (int j=0 ;j<3; j++) m_dUdt_nm2[j] = tab[j+6];
//   for (int j=0 ;j<3; j++) m_dOmegadt_nm2[j] = tab[j+9];
// }
