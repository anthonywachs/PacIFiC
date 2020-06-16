#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "AppSec.H"
#include "Composant.H"
#include "Contact_BuilderFactory.hh"
#include "Memento.hh"
#include "SaveTable.H"
#include "PointContact.hh"
#include "ObstaclePeriodique.hh"
#include "MonObstacle.H"

map<const string, void*, organize> SaveTable::table;
Scalar AppSec::m_overlap_max = 0.;
Scalar AppSec::m_overlap_mean = 0.;
Scalar AppSec::m_time_overlapMax = 0.;
Scalar AppSec::m_nbIterGJK_mean = 0.;
Scalar AppSec::m_nbParticules_mean = 0.;
Scalar AppSec::m_propCoulombFriction = 0.;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur par defaut
AppSec::AppSec() :
  App(),
  m_obstacles( NULL )
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructeur
AppSec::~AppSec()
{}




// ----------------------------------------------------------------------------
// Association des obstacles a l'algorithme
// G.FERRER - Mai 2003 - Creation
void AppSec::Link( Obstacle *obstacle )
{
  if ( !m_obstacles )
  {
    m_obstacles = obstacle;
    m_allObstacles = obstacle->getObstacles();
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Suppression de l'obstacle de l'APPlication
void AppSec::remove( MonObstacle* obs )
{
  list<MonObstacle*>::iterator il;
  for (il=m_allObstacles.begin();il!=m_allObstacles.end();)
    if ( *il == obs ) il = m_allObstacles.erase( il );
    else il++;
}




// ----------------------------------------------------------------------------
// Caracteristiques des contacts
// A.WACHS - Sep .2010 - Creation
void AppSec::computeMeanNbParticules( size_t const& nbPart )
{
  static Scalar global_counter_nbParticules_mean = 0.;
  m_nbParticules_mean = ( m_nbParticules_mean * global_counter_nbParticules_mean
  	+ double(nbPart) ) / ( global_counter_nbParticules_mean + 1. );
  global_counter_nbParticules_mean += 1.;
}




// ----------------------------------------------------------------------------
// Caracteristiques des contacts
void AppSec::addToContactsFeatures( Scalar time,
	PointContact & contactPoint )
{
  static Scalar global_counter_overlap_mean = 0.;
  Scalar overlap = contactPoint.getDistance();
  if ( overlap < 0. )
  {
    overlap = fabs( overlap ) ;

    // Overlap moyen
    m_overlap_mean = ( m_overlap_mean * global_counter_overlap_mean + overlap )
    	/ ( global_counter_overlap_mean + 1. );
    global_counter_overlap_mean += 1.;

    // Overlap max
    if ( overlap > m_overlap_max )
    {
      m_overlap_max = overlap;
      m_time_overlapMax = time;
    }
  }

  static Scalar global_counter_nbIterGJK_mean = 0.;
  int nbGJK = contactPoint.getNbIterGJK();
  if ( nbGJK )
  {
    m_nbIterGJK_mean = ( m_nbIterGJK_mean * global_counter_nbIterGJK_mean
    	+ Scalar(nbGJK) ) / ( global_counter_nbIterGJK_mean + 1. );
    global_counter_nbIterGJK_mean += 1.;
  }

  // static Scalar global_counter_contacts = 0.;
  // float nbCoulombRegime = contactPoint.getNbCoulombRegimes();
  // // if (nbCoulombRegime) cout << "nbCoulombRegime=" << nbCoulombRegime <<endl;
  // // cout << "global_counter_prop_coulomb_friction=" << global_counter_prop_coulomb_friction << endl;
  // if (nbGJK)
  // {
  //       m_propCoulombFriction = (m_propCoulombFriction *
  //           global_counter_contacts + nbCoulombRegime ) /
  //           (global_counter_contacts + 1.);
  //       global_counter_contacts += 1.;
  //   }
}




// ----------------------------------------------------------------------------
// Renvoie la penetration moyenne
Scalar AppSec::getOverlapMean()
{
  return m_overlap_mean;
}




// ----------------------------------------------------------------------------
// Renvoie la penetration maximale
Scalar AppSec::getOverlapMax()
{
  return m_overlap_max;
}




// ----------------------------------------------------------------------------
// Renvoie l'instant de l'overlap max
Scalar AppSec::getTimeOverlapMax()
{
  return m_time_overlapMax;
}




// ----------------------------------------------------------------------------
// Renvoie le nombre d'iterations moyen de GJK
Scalar AppSec::getNbIterGJKMean()
{
  return m_nbIterGJK_mean;
}




// ----------------------------------------------------------------------------
// Renvoie le nombre moyen de particules par proc
Scalar AppSec::getNbParticulesPerProcMean()
{
  return m_nbParticules_mean;
}



// ----------------------------------------------------------------------------
// Renvoie la proportion de régime de coulomb
Scalar AppSec::getNbCoulombRegimes()
{
  return m_propCoulombFriction;
}




// ----------------------------------------------------------------------------
// Affecte les caracteristiques globales des contacts
void AppSec::setContactsFeatures( Scalar const& overlap_max_,
	  Scalar const& overlap_mean_,
	  Scalar const& time_overlapMax_,
	  Scalar const& nbIterGJK_ )
{
  m_overlap_max = overlap_max_;
  m_overlap_mean = overlap_mean_;
  m_time_overlapMax = time_overlapMax_;
  m_nbIterGJK_mean = nbIterGJK_;
}




// ----------------------------------------------------------------------------
// Renvoie un obstacle en fonction de son numero
// A.WACHS - Aout .2010 - Creation
MonObstacle* AppSec::getMonObstacle( int const& num ) const
{
  list<MonObstacle*>::const_iterator il;
  MonObstacle* pp = NULL;
  bool found = false;
  for (il=m_allObstacles.begin();il!=m_allObstacles.end() && !found;il++)
    if ( (*il)->getID() == num )
    {
      pp = *il;
      found = true;
    }

  return pp;
}




// ----------------------------------------------------------------------------
// Renvoie un obstacle periodique en fonction de son numero
// A.WACHS - Aout .2010 - Creation
ObstaclePeriodique const* AppSec::getObstaclePeriodique( int const& num ) const
{
  list<MonObstacle*>::const_iterator il;
  ObstaclePeriodique const* pp = NULL;
  bool found = false;
  for (il=m_allObstacles.begin();il!=m_allObstacles.end() && !found;il++)
    if ( (*il)->getID() == num )
    {
      pp = (*il)->getObstaclePeriodic();
      found = true;
    }

  return pp;
}
