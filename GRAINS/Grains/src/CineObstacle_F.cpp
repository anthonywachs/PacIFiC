#include "CineObstacle_F.H"
#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include <limits>
#include "Obstacle.H"
#include "Torseur.H"


// ----------------------------------------------------------------------------
// Constructeur
// G.FERRER - Aout.2003 - Creation
CineObstacle_F::CineObstacle_F() :
  m_chargement( NULL )
{
}




// ----------------------------------------------------------------------------
// Destructeur
// G.FERRER - Aout.2003 - Creation
CineObstacle_F::~CineObstacle_F()
{
}




// ----------------------------------------------------------------------------
// Ajout du chargement impose a la cinematique
// G.FERRER - Aout.2003 - Creation
void CineObstacle_F::append( ObstacleChargement_F &chargement_ )
{
  m_chargements.push_back( &chargement_ );
}




// ----------------------------------------------------------------------------
// Destruction de l'ensemble des chargements d'obstacle
// G.FERRER - Aout.2003 - Creation
void CineObstacle_F::clearAndDestroy()
{
  list<ObstacleChargement_F*>::iterator iter;

  for (iter=m_chargements.begin(); iter!=m_chargements.end(); iter++)
    delete *iter;
  m_chargements.clear();

  if ( m_chargement ) delete m_chargement;
}




// ----------------------------------------------------------------------------
// DeComposition de la cinematique d'un Obstacle avec une autre cinematique.
// G.FERRER - Aout.2003 - Creation
void CineObstacle_F::Decompose( const CineObstacle_F& voisine, 
	const Point& centre )
{
  Vecteur direction = *(voisine.m_chargement->getDirection()) - centre;
  Scalar ratio = voisine.m_vitesseD / Norm( direction );
  m_vitesseT += direction * ratio;
}




// ----------------------------------------------------------------------------
// Evaluation des vitesses Translation & Rotation dans l'intervalle de temps
// G.FERRER - Aout.2003 - Creation
// D. RAKOTONIRINA - Fev. 2017 - Modification
bool CineObstacle_F::Deplacement( Scalar temps, Scalar dt, 
	Obstacle* obstacle )
{
  // Recherche du chargement dans l'increment de temps
  Scalar fin = temps + dt;
  Scalar dtt = 0.;
  if ( m_chargement ) 
  {
    if ( m_chargement->isActif( temps, fin ) ) 
    {
      dtt = m_chargement->getTime( temps, fin );
      m_vitesseT = *m_chargement->VitesseTranslation( temps, dtt, obstacle );
    }
    else 
    {
      delete m_chargement;
      m_chargement = NULL;
    }
  }

  if ( !m_chargement ) 
  {
    list<ObstacleChargement_F*>::iterator iter;
    for (iter=m_chargements.begin(); iter!=m_chargements.end(); iter++)
      if ( (*iter)->isActif( temps, fin ) ) 
      {
	m_chargement = *iter;
	iter = m_chargements.erase( iter );
	iter = m_chargements.end();
	//dtt = m_chargement->getTime( temps, fin );
      }
  }
  
//  // Il existe un chargement dans l'increment de temps
//  if ( m_chargement ) 
//  {

    // ratio : direction du deplacement pour compenser la force
    // Cas 0 : force de reaction = force de confinement +/- epsilon
    //         deplacement - null
    // Cas 1 : force de reaction < force de confinement - epsilon
    //         deplacement + direction
    // Cas 2 : force de reaction > force de confinement + epsilon
    //         deplacement - direction
    /*
    Scalar ratio = 0.;
    if (force < chargement->getForceImpMin())
      ratio = 1.;
    else if (force > chargement->getForceImpMax())
      ratio = -1.;

    vitesseD = ratio + chargement->getDeplacement() / dt;
    Vecteur direction = chargement->getDirection();
    vitesseT += vitesseD * direction;
    */

//  }
  return ( Norm( m_vitesseT ) != 0. );
}




// ----------------------------------------------------------------------------
// G.FERRER - Aout.2003 - Creation
Vecteur CineObstacle_F::getTranslation( const Scalar dt ) const
{
  return m_vitesseT * dt;
}




// ----------------------------------------------------------------------------
// Mise a zero de la cinematique
// G.FERRER - Aout.2003 - Creation
void CineObstacle_F::reset()
{
  m_vitesseT = 0.;
  m_vitesseD = 0.;
}




// ----------------------------------------------------------------------------
// Vitesse relative de l'obstacle
// G.FERRER - Aout.2003 - Creation
Vecteur CineObstacle_F::VitesseRelative( const Vecteur &om ) const
{
  return m_vitesseT;
}
