#include "Grains.H"
#include "LeapFrog_3D.hh"
#include "Particule.H"
#include "Torseur.H"
#include <math.h>


// ----------------------------------------------------------------------------
// Constructeur par defaut
LeapFrog_3D::LeapFrog_3D() :
  CineParticule()
{}




// ----------------------------------------------------------------------------
// Constructeur par copie
LeapFrog_3D::LeapFrog_3D( const LeapFrog_3D &copie ) :
  CineParticule( copie )
{}




// ----------------------------------------------------------------------------
// Destructeur
LeapFrog_3D::~LeapFrog_3D()
{}




// ----------------------------------------------------------------------------
/** @brief Classname for BuilderFactory */
string LeapFrog_3D::className() const
{
  return ( "*LeapFrog_3D" );
}




// ----------------------------------------------------------------------------
// Construction d'un clone de la cinematique
CineParticule* LeapFrog_3D::clone() const
{
  return ( new LeapFrog_3D( *this ) );
}




// ----------------------------------------------------------------------------
// Evaluation du deplacement  & de la rotation de la particule
void LeapFrog_3D::CalculerVariationQDM( const Torseur &torseur,
	Scalar dt, const Particule* particule )
{
  Scalar facteurCouplage = 1.;
  if ( !Particule::getExplicitMassCorrection() )
    facteurCouplage -=
    	Particule::getFluideMasseVolumique() / particule->getMasseVolumique();

  // Quantite de mouvement translationnelle
  Scalar masseReduite = particule->getMasse() * facteurCouplage;
  m_dUdt = *torseur.getForce() / masseReduite;
  m_dUdt.round();

  // Quantite de mouvement rotationnelle
  // Calcul de dOmegadt = I^-1.(Mc + (I.omega) ^ omega)
  // Soit p le quaternion lie � la rotation de la configuration de base R
  // � la configuration courante R'. Ainsi, soit q' un quaternion dans la base
  // R', dans R, on a: q = p.q'.pConjugue, o� pConjugue est le conjugu� de p.
  // La m�me relation est valable pour les matrices.
  // Notons que:
  // # omega = 2 * m_QuaternionVitesseR * pConjugue (cf rapport D. Petit),
  // ou m_QuaternionVitesseR
  // est la d�riv�e par rapport au temps du quaternion de rotation instantan�e.
  // # l'�quation dOmegadt = I^-1.(Mc + (I.omega) ^ omega) est r�solue dans R
  // mais I et I^-1 sont calcul�s dans R', il faut donc utiliser le quaternion
  // de rotation p pour changer de rep�re.
  // Une fois d�velop�e, dOmegadt s'�crit:
  //   dOmegadt = p.I^-1.pConjugue.(Mc
  //            + 4.p.I.pConjugue.m_QuaternionVitesseR.pConjugue
  //		^ m_QuaternionVitesseR.pConjugue)
  // En pratique, les calculs type I.v ou I^-1.v s'effectuent en 3 �tapes:
  // 1. changer de rep�re pour v c.a.d. �crire v dans le rep�re de I c.a.d. R'
  // 2. effectuer le produit matrice-vecteur dans le rep�re de I c.a.d. R'
  // 3. changer de rep�re pour le vecteur r�sultat de R' � R
  // Exemple:
  // 1. v' = pConjugue.v.p
  // 2. res' = I.v'
  // 3. res = p.res'.pConjugue
  // donc au final on obtient: res = p.I.v'.pConjugue
  //                               = p.I.pConjugue.v.p.pConjugue
  //                               = p.I.pConjugue.v
  // ce qui revient � effectuer un changement de rep�re pour I car
  // p.I.pConjugue = I dans R
  // Pour calculer A=p.I.pConjugue.m_QuaternionVitesseR.pConjugue, on calcule
  // d'abord I.pConjugue.m_QuaternionVitesseR, puis on effectue un changement
  // de rep�re de R'->R
  // En fait, I.pConjugue.m_QuaternionVitesseR =
  //                        I.pConjugue.m_QuaternionVitesseR.pConjugue.p
  // et m_QuaternionVitesseR.pConjugue = 0.5.omega, donc, on a:
  // A = p.(I.(pConjugue.(0.5.omega).p)).pConjugue
  // Parenth�ses ainsi, on retrouve la m�thode de calcul g�n�rique en 3 �tapes
  // d�taillee plus haut et appliquee � 0.5.omega.
  Quaternion pConjugue = m_QuaternionRotation.Conjugate();
  const double *inertieInverse = particule->getInertieInverse();
  Vecteur vTmp,work;
  if ( Grains::isModePredictor() )
  {
    // Calcul de (I.omega) ^ omega =
    //   4.p.I.pConjugue.m_QuaternionVitesseR.pConjugue^m_QuaternionVitesseR
    //   .pConjugue

    // Calcul de vTmp = pConjugue.m_QuaternionVitesseR
    vTmp = pConjugue.multToVecteur( m_QuaternionVitesseR );
    const double *inertie = particule->getInertie();

    // Calcul de I.pConjugue.m_QuaternionVitesseR dans R'
    work[0] = inertie[0]*vTmp[0] + inertie[1]*vTmp[1] + inertie[2]*vTmp[2];
    work[1] = inertie[1]*vTmp[0] + inertie[3]*vTmp[1] + inertie[4]*vTmp[2];
    work[2] = inertie[2]*vTmp[0] + inertie[4]*vTmp[1] + inertie[5]*vTmp[2];

    // Changement de rep�re de R'->R:
    // 4.p.I.pConjugue.m_QuaternionVitesseR.pConjugue
    vTmp = 2. * m_QuaternionRotation.multToVecteur( ( work , pConjugue ) );

    // Calcul de 4.p.I.pConjugue.m_QuaternionVitesseR.pConjugue
    //    ^ m_QuaternionVitesseR.pConjugue
    //    = 2.p.I.pConjugue.m_QuaternionVitesseR.pConjugue ^ m_vitesseR
    work = vTmp ^ m_vitesseR;

    // Calcul de A=Mc+4.p.I.pConjugue.m_QuaternionVitesseR.pConjugue
    //    ^ m_QuaternionVitesseR.pConjugue
    // en ajoutant le moment des forces exterieures
    work += *torseur.getMoment();
  }
  else
    // lorque le mode est correcteur il ne faut pas calculer le terme
    // (I.omega) ^ omega = -(omega ^ (I.omega)) !!
    work = *torseur.getMoment();

  // Changement de rep�re R->R' pour A: A'=pConjugue.A.p
  vTmp = pConjugue.multToVecteur( ( work , m_QuaternionRotation  ));
  vTmp.round();

  // Calcul de I^-1.A' dans R'
  work[0] = inertieInverse[0]*vTmp[0] + inertieInverse[1]*vTmp[1]
    + inertieInverse[2]*vTmp[2];
  work[1] = inertieInverse[1]*vTmp[0] + inertieInverse[3]*vTmp[1]
    + inertieInverse[4]*vTmp[2];
  work[2] = inertieInverse[2]*vTmp[0] + inertieInverse[4]*vTmp[1]
    + inertieInverse[5]*vTmp[2];

  // Changement de rep�re de R'->R: p.I^-1.A'.pConjugue
  vTmp = m_QuaternionRotation.multToVecteur( ( work , pConjugue ) );

  // Calcul de m_dOmegadt
  m_dOmegadt = vTmp / facteurCouplage;
}




// ----------------------------------------------------------------------------
// Operateur d'ecriture
ostream &operator << ( ostream &fileOut,
	const LeapFrog_3D &cinematique )
{
  fileOut << "*LeapFrog_3D\n";
  fileOut << cinematique.m_vitesseT
	<< cinematique.m_QuaternionRotation
	<< cinematique.m_QuaternionVitesseR;

  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Operateur de lecture
istream &operator >> ( istream &fileIn,
	LeapFrog_3D &cinematique )
{
  fileIn >> cinematique.m_vitesseT
	>> cinematique.m_QuaternionRotation
	>> cinematique.m_QuaternionVitesseR;
  cinematique.m_vitesseR =
  	2.0 * cinematique.m_QuaternionVitesseR.multConjugateToVecteur(
  	cinematique.m_QuaternionRotation );

  return ( fileIn );
}




// ----------------------------------------------------------------------------
// Calcul du terme de masse ajout� explicite sur les moments
Vecteur LeapFrog_3D::calculIdwExplicite( const Vecteur &dw,
  	const double *inertie ) const
{
  Vecteur Idw;
  Quaternion pConjugue = m_QuaternionRotation.Conjugate();

  // Changement de rep�re de R->R'
  Vecteur vTmp = pConjugue.multToVecteur( ( dw , m_QuaternionRotation ) );

  // Calcul de I.pConjugue.dw.p dans R'
  Vecteur work;
  work[0] = inertie[0]*vTmp[0] + inertie[1]*vTmp[1] + inertie[2]*vTmp[2];
  work[1] = inertie[1]*vTmp[0] + inertie[3]*vTmp[1] + inertie[4]*vTmp[2];
  work[2] = inertie[2]*vTmp[0] + inertie[4]*vTmp[1] + inertie[5]*vTmp[2];

  // Changement de rep�re R'->R
  Idw = m_QuaternionRotation.multToVecteur( ( work , pConjugue ) );

  return Idw;
}
