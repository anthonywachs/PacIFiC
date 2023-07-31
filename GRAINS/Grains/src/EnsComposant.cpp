#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "EnsComposant.H"
#include "App.H"
#include "AppSec.H"
#include "FormeVdW.H"
#include "CompObstacle.H"
#include "Obstacle_BuilderFactory.H"
#include "PostProcessingWriter.hh"
#include "Paraview_PostProcessingWriter.hh"
#include "Matlab_PostProcessingWriter.hh"
#include "Grains_BuilderFactory.H"
#include <math.h>
#include <stdlib.h>


double EnsComposant::SauterMeanDiameter = 0.;
double EnsComposant::sumXY = 0.;

// ----------------------------------------------------------------------------
// Constructeur
EnsComposant::EnsComposant() :
  m_wait( NULL ),
  m_nb_total_particules( 0 ),
  m_obstacle( NULL ),
  m_hasSerialPostProcessors( false ),
  m_outputTorseurObstacles_counter( 1 ),
  m_outputTorseurObstacles_frequency( 0 ),
  m_initial_time( 0. )
{
  m_obstacle = new CompObstacle( "__AllObstacles___" );
  Composant::setNbComposantsCrees( 0 );
}




// ----------------------------------------------------------------------------
// Destructeur
EnsComposant::~EnsComposant()
{
  list<Particule*>::iterator particule;

  // Remarque: les particules periodiques (id = -2) sont detruites
  // par ce destructeur.
  // Dans le destructeur de Particule, la liste de pointeurs
  // sur les particules periodiques est simplement vid�e mais la destruction
  // des objets point�s est r�alis�e ici
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    if ( (*particule)->getID() != -2 ) delete *particule;
  for (particule=m_pwait.begin(); particule!=m_pwait.end();
       	particule++)  delete *particule;
  for (particule=m_particulesClonesPeriodiques.begin();
  	particule!=m_particulesClonesPeriodiques.end();
       	particule++)  delete *particule;
  m_particulesActives.clear();
  m_pwait.clear();
  m_particulesClonesPeriodiques.clear();

  m_particulesHalozone.clear();
  m_particulesClones.clear();

  vector<Particule*>::iterator ivp;
  for (ivp=m_ParticuleClassesReference.begin();
  	ivp!=m_ParticuleClassesReference.end(); ivp++)
    delete *ivp;
  m_ParticuleClassesReference.clear();

  delete m_obstacle;

  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    delete *pp;
  m_postProcessors.clear();

  m_ChargementsCinematiques.clear();
  m_ChargementsForces.clear();
}




// ----------------------------------------------------------------------------
// Actualisation des relations entre composants
// Si particule inactive, bascule vers liste attente
void EnsComposant::Actualiser()
{
  list<Particule*>::iterator particule;
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); )
  {
    switch ( (*particule)->getActivity() )
    {
      case COMPUTE :
        particule++;
        break;

      case CLEARandWAIT:
        (*particule)->reset();
        m_pwait.push_back(*particule);
        particule = m_particulesActives.erase( particule );
        break;

      default:
        break;
    }
  }
}




// ----------------------------------------------------------------------------
// Ajout d'une particule
// G.FERRER - Aout.2001 - Creation
void EnsComposant::Ajouter( Particule* particule )
{
  switch ( particule->getActivity() )
  {
    case WAIT:
      m_pwait.push_back(particule);
      break;

    case COMPUTE:
      m_particulesActives.push_back(particule);
      break;

    default:
      break;
  }
}




// ----------------------------------------------------------------------------
// Ajout d'une classe de particules
void EnsComposant::AjouterClasseParticules( Particule *particule )
{
  m_ParticuleClassesReference.reserve( m_ParticuleClassesReference.size() + 1 );
  m_ParticuleClassesReference.push_back( particule );
}




// ----------------------------------------------------------------------------
// Ajout d'un obstacle
// G.FERRER - Nove.2000 - Creation
void EnsComposant::Ajouter( Obstacle* obstacle_ )
{
  if ( !m_obstacle ) m_obstacle = obstacle_;
  else m_obstacle->append(obstacle_);
}




// ----------------------------------------------------------------------------
// Association du chargement a son obstacle
// G.FERRER - Nove.2000 - Creation
void EnsComposant::Associer( ObstacleChargement &chargement )
{
  m_obstacle->Associer( chargement );
  m_ChargementsCinematiques.push_back( &chargement );
}




// ----------------------------------------------------------------------------
// Association du chargement a son obstacle
// G.FERRER - Nove.2000 - Creation
void EnsComposant::Associer( ObstacleChargement_F &chargement )
{
  m_obstacle->Associer( chargement );
  m_ChargementsForces.push_back( &chargement );
}




// ----------------------------------------------------------------------------
// Nombre de particules actives de tag 0 ou 1
size_t EnsComposant::nbreParticulesActivesOnProc() const
{
  size_t nb_part = m_particulesActives.size();
  for (list<Particule*>::const_iterator il=m_particulesActives.begin();
  	il!=m_particulesActives.end();il++)
    if ( (*il)->getTag() == 2 || (*il)->getID() == -2 ) nb_part--;

  return nb_part;
}




// ----------------------------------------------------------------------------
// Deplacement des particules & obstacles
// G.FERRER - Janv.2000 - Creation
list<MonObstacle*> EnsComposant::Deplacer( Scalar temps, Scalar dt )
  throw(ErreurDeplacement)
{
  // Deplacement des particules
  list<Particule*>::iterator particule;
  for (particule=m_particulesActives.begin();
      particule!=m_particulesActives.end(); particule++)
    if ( (*particule)->getTag() != 2 && (*particule)->getMobilite() )
      (*particule)->Deplacer( temps, dt );

  // Deplacement des obstacles
  list<MonObstacle*> obstaclesDeplaces;
  if ( !m_ChargementsCinematiques.empty() || !m_ChargementsForces.empty() )
  {
    m_obstacle->resetCinematique();
    obstaclesDeplaces = m_obstacle->Deplacer( temps, dt, false, false );

    list<ObstacleChargement*>::iterator chargement;
    for (chargement=m_ChargementsCinematiques.begin();
  	chargement!=m_ChargementsCinematiques.end(); )
      if ( (*chargement)->isCompleted( temps, dt ) )
        chargement = m_ChargementsCinematiques.erase( chargement );
      else chargement++;

    list<ObstacleChargement_F*>::iterator chargement_F;
    for (chargement_F=m_ChargementsForces.begin();
  	chargement_F!=m_ChargementsForces.end(); )
    {
      if ( (*chargement_F)->isCompleted( temps, dt ) )
        chargement_F = m_ChargementsForces.erase( chargement_F );
      else chargement_F++;
    }
  }

  return obstaclesDeplaces;
}




// ----------------------------------------------------------------------------
// Post-processing des contraintes
// D. RAKOTONIRINA - Mars 2017 - Creation
void EnsComposant::computeStressTensor( vector<Scalar> &stressTensor,
  MPIWrapperGrains const* wrapper )
{
  vector<Scalar> const* internalMoment = NULL;
  // Particules
  list<Particule*>::const_iterator particule;
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
  {
    if ( (*particule)->getTag() != 2 )
    {
      internalMoment = (*particule)->getInternalMoment();
      stressTensor[0] += (*internalMoment)[0];
      stressTensor[1] += (*internalMoment)[1];
      stressTensor[2] += (*internalMoment)[2];
      stressTensor[3] += (*internalMoment)[3];
      stressTensor[4] += (*internalMoment)[4];
      stressTensor[5] += (*internalMoment)[5];
      stressTensor[6] += (*internalMoment)[6];
      stressTensor[7] += (*internalMoment)[7];
      stressTensor[8] += (*internalMoment)[8];
    }
  }

  if ( wrapper )
  {
    stressTensor[0] = wrapper->sum_DOUBLE_master( stressTensor[0] );
    stressTensor[1] = wrapper->sum_DOUBLE_master( stressTensor[1] );
    stressTensor[2] = wrapper->sum_DOUBLE_master( stressTensor[2] );
    stressTensor[3] = wrapper->sum_DOUBLE_master( stressTensor[3] );
    stressTensor[4] = wrapper->sum_DOUBLE_master( stressTensor[4] );
    stressTensor[5] = wrapper->sum_DOUBLE_master( stressTensor[5] );
    stressTensor[6] = wrapper->sum_DOUBLE_master( stressTensor[6] );
    stressTensor[7] = wrapper->sum_DOUBLE_master( stressTensor[7] );
    stressTensor[8] = wrapper->sum_DOUBLE_master( stressTensor[8] );
  }
}




// ----------------------------------------------------------------------------
// Post-processing des contraintes
// D. RAKOTONIRINA - Fev 2017 - Creation
void EnsComposant::setStressTensorDomain( vector<Fenetre> fenetres )
{
  // Domaine des particules en entier sans fenetre
  if ( fenetres.size() == 0 )
  {
    double ox, oy, oz, lx, ly, lz;
    App::getOrigineGlobale( ox, oy, oz );
    App::getDimesionsGlobales( lx, ly, lz );
    Grains_Exec::m_stressTensorDomain[0] = ox;
    Grains_Exec::m_stressTensorDomain[1] = lx + ox;
    Grains_Exec::m_stressTensorDomain[2] = oy;
    Grains_Exec::m_stressTensorDomain[3] = ly + oy;

    // Take the current position of MonObstacle to compute the volume
    list<MonObstacle*> obstacles = m_obstacle->getObstacles();
    //list<MonObstacle*>::const_iterator myObs;
    //for (myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++)
    //{
    //  tmp = (*(*myObs)->getPosition())[Z] -
    //      (*myObs)->getObstacleBox()->getExtent()[Z];
    //  posZ = posZ > tmp ? posZ : tmp;
    //}
    Scalar posZ1, posZ2;
    posZ1 = (*(*obstacles.begin())->getPosition())[Z];
    posZ2 = (*(obstacles.back())->getPosition())[Z];
    Scalar extent1 = (*obstacles.begin())->getObstacleBox()->getLower(2);
    Scalar extent2 = (obstacles.back())->getObstacleBox()->getUpper(2);

    if ( posZ1 > posZ2 )
    {
      // Position de l'obstacle superieur
      Grains_Exec::m_stressTensorDomain[5] = posZ1 - extent1;
      // Position de l'obstacle inferieur
      Grains_Exec::m_stressTensorDomain[4] = posZ2 + extent2;
    }
    else
    {
      // Position de l'obstacle inferieur
      Grains_Exec::m_stressTensorDomain[4] = posZ1 + extent1;
      // Position de l'obstacle superieur
      Grains_Exec::m_stressTensorDomain[5] = posZ2 - extent2;
    }
  }
  else
  {
    // Pour le moment une seule fenetre
    Point fA, fB;
    fA = fenetres[0].ptA;
    fB = fenetres[0].ptB;
    Grains_Exec::m_stressTensorDomain[0] = fA[X];
    Grains_Exec::m_stressTensorDomain[1] = fB[X];
    Grains_Exec::m_stressTensorDomain[2] = fA[Y];
    Grains_Exec::m_stressTensorDomain[3] = fB[Y];
    Grains_Exec::m_stressTensorDomain[4] = fA[Z];
    Grains_Exec::m_stressTensorDomain[5] = fB[Z];
  }

//  // Pour le moment une seule fenetre
//  Point fA, fB;
//  Scalar tmp = 0., posZ = 0.;
//  fA = fenetres[0].ptA;
//  fB = fenetres[0].ptB;
//  Grains_Exec::m_stressTensorDomain[0] = fA[X];
//  Grains_Exec::m_stressTensorDomain[1] = fA[Y];
//  Grains_Exec::m_stressTensorDomain[2] = fB[X];
//  Grains_Exec::m_stressTensorDomain[3] = fB[Y];
//
//  // Take the position of MonObstacle to compute the volume
//  list<MonObstacle*> obstacles = m_obstacle->getObstacles();
//  list<MonObstacle*>::const_iterator myObs;
//  for (myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++)
//  {
//    tmp = (*(*myObs)->getPosition())[Z] -
//	(*myObs)->getObstacleBox()->getExtent()[Z];
//    posZ = posZ > tmp ? posZ : tmp;
//  }
//  Grains_Exec::m_stressTensorDomain[4] = posZ;
}




// ----------------------------------------------------------------------------
// Compute the void ratio : Volume of void / Volume of particles
// D. RAKOTONIRINA - Mai 2017 - Creation
void EnsComposant::ComputeVoidRatio( Scalar &volume, vector<Point> &obsPos,
    MPIWrapperGrains const* wrapper ) const
{
  volume = getVolumeIn();

  list<MonObstacle*> obstacles = m_obstacle->getObstacles();
  obsPos = vector<Point>( obstacles.size(), Point(0.,0.,0.) );
  list<MonObstacle*>::const_iterator myObs;
  size_t ii = 0;
  for (myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++)
  {
    obsPos[ii] = *(*myObs)->getPosition();
    ++ii;
  }

//  LX = fabs( Grains_Exec::m_stressTensorDomain[1]
//           - Grains_Exec::m_stressTensorDomain[0] );
//  LY = fabs( Grains_Exec::m_stressTensorDomain[3]
//           - Grains_Exec::m_stressTensorDomain[2] );
//  LZ = fabs( Grains_Exec::m_stressTensorDomain[5]
//           - Grains_Exec::m_stressTensorDomain[4] );
//  volTotal = LX*LY*LZ;
//
//  volVoid = volTotal - volParticle;
//  voidRatio = volVoid / volTotal;
  if ( wrapper )
  {
    volume = wrapper->sum_DOUBLE_master( volume );
    volume = wrapper->Broadcast_DOUBLE( volume );
//    volVoid = volTotal - volParticle;
//    voidRatio = volVoid / volTotal;
//    voidRatio = wrapper->Broadcast_DOUBLE( voidRatio );
  }
}




// ----------------------------------------------------------------------------
// Solid-Bodies temperature evolution
// Manuel - July 2015 - Creation
void EnsComposant::ComputeTemperature( Scalar temps, Scalar dt )
{
  list<Particule*>::iterator il;
  for (il=m_particulesActives.begin(); il!=m_particulesActives.end(); il++)
    if ( (*il)->getTag() != 2 )
    {
     if ((*il)->getTempEvolution())
     {
       (*il)->ComputeTemperature( temps, dt );
     }
    }

}




// ----------------------------------------------------------------------------
// Initialise la cinematique des obstacles
void EnsComposant::setCinematiqueObstacleSansDeplacement(
	Scalar temps, Scalar dt )
{
  bool bbb = Obstacle::getDeplaceObstacle() ;

  Obstacle::setDeplaceObstacle( false ) ;

  // Deplacement des obstacles
  if ( !m_ChargementsCinematiques.empty() )
  {
    m_obstacle->resetCinematique();
    m_obstacle->Deplacer( temps, dt, false, false );

    list<ObstacleChargement*>::iterator chargement;
    for (chargement=m_ChargementsCinematiques.begin();
  	chargement!=m_ChargementsCinematiques.end(); )
      if ( (*chargement)->isCompleted( temps, dt ) )
        chargement = m_ChargementsCinematiques.erase( chargement );
      else chargement++;
  }
  Obstacle::setDeplaceObstacle( bbb ) ;
}




// ----------------------------------------------------------------------------
// Deplacement des particules & obstacles
// A.WACHS - Aout.2010 - Creation
void EnsComposant::AddForcesFromPeriodicClonesToParticules(
	Scalar temps, Scalar dt )
{
  list<Particule*>::iterator particule;
  for (particule=m_particulesClonesPeriodiques.begin();
  	particule!=m_particulesClonesPeriodiques.end(); particule++)
    (*particule)->AddForcesFromPeriodicCloneToParticule( temps, dt,
    	Grains_Exec::m_ContactforceOutput );
}




// ----------------------------------------------------------------------------
// Initialise le torseur de force des particules
// A.WACHS - Aout.2011 - Creation
void EnsComposant::InitializeForces( Scalar temps, Scalar dt,
	bool const& withWeight )
{
  list<Particule*>::iterator particule;

  // Particules actives
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    (*particule)->InitializeForce( withWeight );

  // Clones periodiques
  for (particule=m_particulesClonesPeriodiques.begin();
  	particule!=m_particulesClonesPeriodiques.end(); particule++)
    (*particule)->InitializeForce( false );

  // Obstacles
  m_obstacle->InitializeForce( false );

  // Initialisation of forces at the contact point for post-processing of
  // the stress tensor
//  if ( Grains_Exec::m_stressTensor || Grains_Exec::m_particleStressTensor )
//    InitializeForcesAtContactPoint();
}




// ----------------------------------------------------------------------------
// Initialise force at the contact point for the post-processing
// of the stress tensor
// D. RAKOTONIRINA - Fev. 2017 - Creation
void EnsComposant::InitializeForcesAtContactPoint()
{
  list<Particule*>::iterator particule;

  // Particules actives
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    (*particule)->InitializeForceAtContactPoint( );

  // Clones periodiques
  for (particule=m_particulesClonesPeriodiques.begin();
  	particule!=m_particulesClonesPeriodiques.end(); particule++)
    (*particule)->InitializeForceAtContactPoint( );
}




// ----------------------------------------------------------------------------
// Initialise le torseur de force PostProcessing
// A.ESTEGHAMATIAN - Aout.2015 - Creation
void EnsComposant::InitializePostProcessingForces()
{
  list<Particule*>::iterator particule;

  // Particules actives
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    (*particule)->InitializePostProcessingForce( );

  // Clones periodiques
  for (particule=m_particulesClonesPeriodiques.begin();
  	particule!=m_particulesClonesPeriodiques.end(); particule++)
    (*particule)->InitializePostProcessingForce( );

  // Obstacles
  m_obstacle->InitializePostProcessingForce( );
}




// ----------------------------------------------------------------------------
// Initialise le torseur de force PostProcessing
// A.ESTEGHAMATIAN - Aout.2015 - Creation
void EnsComposant::IntegrateContactPostProcessingForces(Scalar nb)
{
  list<Particule*>::iterator particule;

  // Particules actives
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    (*particule)->IntegrateContactPostProcessingForce( nb );

  // Clones periodiques
  for (particule=m_particulesClonesPeriodiques.begin();
  	particule!=m_particulesClonesPeriodiques.end(); particule++)
    (*particule)->IntegrateContactPostProcessingForce( nb );

  // Obstacles
  m_obstacle->IntegrateContactPostProcessingForce( nb );
}




// ----------------------------------------------------------------------------
// Initialise l'indicateur de calcul de la transformation
// avec epaiseur de croute a faux pour tous les composants
// A.WACHS - Fev.2012 - Creation
void EnsComposant::InitializeVdWState( Scalar temps, Scalar dt )
{
  list<Particule*>::iterator particule;

  // Particules actives
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    (*particule)->initializeVdWtransform_to_notComputed();
//    (*particule)->getForme()->initializeVdWtransform_to_notComputed();

  // Clones periodiques
  for (particule=m_particulesClonesPeriodiques.begin();
  	particule!=m_particulesClonesPeriodiques.end(); particule++)
    (*particule)->initializeVdWtransform_to_notComputed();
//    (*particule)->getForme()->initializeVdWtransform_to_notComputed();

  // Obstacles
  list<MonObstacle*> list_obstacles = m_obstacle->getObstacles();
  list<MonObstacle*>::iterator myObs;
  for (myObs=list_obstacles.begin();myObs!=list_obstacles.end();myObs++)
    (*myObs)->getForme()->initializeVdWtransform_to_notComputed();
}




// ----------------------------------------------------------------------------
// Calcul le poids de toutes les particules
// A.WACHS - Fev.2012 - Creation
void EnsComposant::computeWeight( Scalar temps, Scalar dt )
{
  list<Particule*>::iterator particule;
  vector<Particule*>::iterator ivp;

  // Classes de reference
  for (ivp=m_ParticuleClassesReference.begin();
  	ivp!=m_ParticuleClassesReference.end(); ivp++)
    (*ivp)->computeWeight();

  // Particules en attente
  for (particule=m_pwait.begin(); particule!=m_pwait.end();
       	particule++)
    (*particule)->computeWeight();

  // Particules actives
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    (*particule)->computeWeight();
}




// ----------------------------------------------------------------------------
// Recherche de l'obstacle "nom"
// G.FERRER - Juil.2003 - Creation
const Obstacle* EnsComposant::getObstacle( const string &nom ) const
{
  return m_obstacle->getNom( nom );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Acces aux obstacles
// G.FERRER - Fevr.2004 - Creation
Obstacle* EnsComposant::getObstacles()
{
  return m_obstacle;
}




// ----------------------------------------------------------------------------
// Particule demandee a l'indice id
// G.FERRER - Janv.2000 - Creation
// A.WACHS - Dec.2009 - Modification
Particule* EnsComposant::getParticule( int id )
{
  Particule *particule = NULL;
  list<Particule*>::iterator iter;
  bool found = false;
  for (iter=m_particulesActives.begin();iter!=m_particulesActives.end()
  	&& !found;iter++)
    if ( (*iter)->getID() == id )
    {
      particule = *iter;
      found = true;
    }

  return ( particule );
}




// ----------------------------------------------------------------------------
// Composant demandee a l'indice id
// A.WACHS - Dec.2009 - Creation
Composant* EnsComposant::getComposant( int id )
{
  Composant *composant = getParticule( id );
  if ( composant == NULL )
  {
    bool found = false;
    list<MonObstacle*> obstacles = m_obstacle->getObstacles();
    list<MonObstacle*>::iterator myObs;
    for (myObs=obstacles.begin(); myObs!=obstacles.end() && !found; myObs++)
      if ( (*myObs)->getID() == id )
      {
        composant = *myObs;
        found = true;
      }
  }

  return ( composant );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Particule cliente en fonction du mode de tirage.
Particule* EnsComposant::getParticule( PullMode mode,
	MPIWrapperGrains const* wrapper )
{
  if ( !m_pwait.empty() )
  {
    switch ( mode )
    {
      case ORDER:
        m_wait = m_pwait.front();
        break;

      case RANDOM:
        if ( m_wait == NULL )
	{
	  Scalar v = double(random()) / double(INT_MAX);
	  int id = int( double(m_pwait.size()) * v );

	  // Parall�le: afin que le tirage al�atoire soit le m�me sur tous les
	  // procs, seul le master envoie la position tiree au hasard dans la
	  // liste aux autres procs
	  if ( wrapper ) id = wrapper->Broadcast_INT( id );

	  list<Particule*>::iterator p = m_pwait.begin();
	  for (int i=0; i<id && p!=m_pwait.end(); i++, p++) {}
	  m_wait = *p;
        }
        break;

      case NONE:
        m_wait = NULL;
        break;
      }
  }
  else m_wait = NULL;

  return m_wait;
}




// ----------------------------------------------------------------------------
// Liste des particules
// G.FERRER - Juil.2000 - Creation
list<Particule*>* EnsComposant::getParticulesActives()
{
  return ( &m_particulesActives );
}




// ----------------------------------------------------------------------------
// Liste des particules
// A.WACHS - Nov.2009 - Creation
list<Particule*> const* EnsComposant::getParticulesActives() const
{
  return ( &m_particulesActives );
}




// ----------------------------------------------------------------------------
// Liste des particules en attente
// G.FERRER - Juil.2003 - Creation
list<Particule*>* EnsComposant::getParticulesWait()
{
  return ( &m_pwait );
}




// ----------------------------------------------------------------------------
// Liste des particules en attente
// A.WACHS - Nov.2009 - Creation
list<Particule*> const* EnsComposant::getParticulesWait() const
{
  return ( &m_pwait );
}




// ----------------------------------------------------------------------------
// Liste des particules dans la zone de recouvrement d'un autre processeur.
list<Particule*>* EnsComposant::getParticulesHalozone()
{
  return ( &m_particulesHalozone );
}




// ----------------------------------------------------------------------------
// Liste des particules clones.
list<Particule*>* EnsComposant::getParticulesClones()
{
  return ( &m_particulesClones );
}




// ----------------------------------------------------------------------------
// Liste des particules clones periodiques
list<Particule*>* EnsComposant::getParticulesClonesPeriodiques()
{
  return ( &m_particulesClonesPeriodiques );
}




// ----------------------------------------------------------------------------
// Liste des particules possedant des clones periodiques
set<Particule*>* EnsComposant::getParticulesReferencesPeriodiques()
{
  return ( &m_particulesReferencesPeriodiques );
}




// ----------------------------------------------------------------------------
// Vecteur des particules de r�f�rence des classes.
vector<Particule*>* EnsComposant::getParticuleClassesReference()
{
  return ( &m_ParticuleClassesReference );
}




// ----------------------------------------------------------------------------
// Vecteur des particules de r�f�rence des classes.
vector<Particule*> const* EnsComposant::getParticuleClassesReference() const
{
  return ( &m_ParticuleClassesReference );
}



// ----------------------------------------------------------------------------
// Modif Manu 04/2015 list<Obstacle*> instead of list<MonObstacle*>
// for CylindricalBox
list<Obstacle*> EnsComposant::getObstaclesToFluid() const
{
  return ( m_obstacle->getObstaclesToFluid() );
}




// ----------------------------------------------------------------------------
// Rayon maximal
// G.FERRER - Aout.2001 - Creation
// A.WACHS - D�c.2009 - Modification
Scalar EnsComposant::getRayonMax()
{
  Scalar rayonMax = 0.0;
  Scalar rayon;

  vector<Particule*>::iterator particule;
  for (particule=m_ParticuleClassesReference.begin();
  	particule!=m_ParticuleClassesReference.end(); particule++)
  {
    rayon = (*particule)->getRayon();
    rayonMax = rayonMax > rayon ? rayonMax : rayon;
  }

  return ( rayonMax );
}




// ----------------------------------------------------------------------------
// Rayon minimal
// G.FERRER - Aout.2001 - Creation
// A.WACHS - D�c.2009 - Modification
Scalar EnsComposant::getRayonMin()
{
  Scalar rayonMin = 1.e10;
  Scalar rayon;

  vector<Particule*>::iterator particule;
  for (particule=m_ParticuleClassesReference.begin();
  	particule!=m_ParticuleClassesReference.end(); particule++)
  {
    rayon = (*particule)->getRayon();
    rayonMin = rayonMin > rayon ? rayon : rayonMin;
  }

  return ( rayonMin );
}




// ----------------------------------------------------------------------------
// Call shrinking functions if shrinking
// choice is 1
// M.SULAIMAN - Nov.2015 - Creation
void EnsComposant::ShrinkingRate( Scalar CurrentTime )
{
  //Scalar CurrentRadius ,CurrentWeight;
  //CurrentRadius = update_shrinking_radius (exp(-CurrentTime));
  //CurrentWeight = compute_shrinking_weight ();
}




// ----------------------------------------------------------------------------
// retourne l'information si la particule est retrecissante
// M.SULAIMAN - Nov.2015 - Creation
// I'm assuming here that  we have only one class
//or many classes of the same initial radius
bool EnsComposant::IsShrinking()
{
  int Choice=3;
  int tempvar=0;
  vector<Particule*>::iterator particule;

  for (particule=m_ParticuleClassesReference.begin();
        particule!=m_ParticuleClassesReference.end(); particule++)
  {
    tempvar++;
    Choice = (*particule)->getShrinkingMode();
  }

  if( Choice == 1 ) return true; else return false;
}




// ----------------------------------------------------------------------------
// recupere le rayon initial
// M.SULAIMAN - Nov.2015 - Creation
Scalar EnsComposant::get_initial_radius()
{
  Scalar initial_radius=0.0;
  vector<Particule*>::iterator particule;

  for (particule=m_ParticuleClassesReference.begin();
        particule!=m_ParticuleClassesReference.end(); particule++)
  {
    initial_radius = (*particule)->getRayon();
  }

  return ( initial_radius );
}




// ----------------------------------------------------------------------------
// recupere la masse initiale
// M.SULAIMAN - Nov.2015 - Creation
Scalar EnsComposant::get_initial_mass()
{

  Scalar initial_mass=0.0;

  vector<Particule*>::iterator particule;
  for (particule=m_ParticuleClassesReference.begin();
        particule!=m_ParticuleClassesReference.end(); particule++)
    initial_mass = (*particule)->getMasse();

  return ( initial_mass );
}




// ----------------------------------------------------------------------------
// update le rayon
// M.SULAIMAN - Nov.2015 - Creation
Scalar EnsComposant::update_shrinking_radius( Scalar rate)
{
  Scalar initial_radius=get_initial_radius();
  Scalar rayon=0.0;
  list<Particule*>::iterator particule;

  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
  {
    (*particule)->set_shrinking_radius(initial_radius*rate);
    rayon = (*particule)->getRayon();
  }

  return rayon;
}




// ----------------------------------------------------------------------------
// update la masse
// M.SULAIMAN - Nov.2015 - Creation
Scalar EnsComposant::compute_shrinking_weight( )
{
  Scalar CurrentMass = 0.0;
  list<Particule*>::iterator particule;

  for (particule=m_particulesActives.begin();
                 particule!=m_particulesActives.end(); particule++)
    CurrentMass = (*particule)->set_shrinking_mass();

  return CurrentMass ;
}




// ----------------------------------------------------------------------------
// Rayon d'interaction maximal
// A.WACHS - D�c.2009 - Modification
Scalar EnsComposant::getRayonInteractionMax()
{
  Scalar rayonInteractionMax = 0.0;
  Scalar rayonInteraction;

  vector<Particule*>::iterator particule;
  for (particule=m_ParticuleClassesReference.begin();
  	particule!=m_ParticuleClassesReference.end(); particule++)
  {
    rayonInteraction = (*particule)->getRayonInteraction();
    rayonInteractionMax = rayonInteractionMax > rayonInteraction ?
    	rayonInteractionMax : rayonInteraction;
  }

  return ( rayonInteractionMax );
}




// ----------------------------------------------------------------------------
// Rayon d'interaction minimal
// A.WACHS - D�c.2009 - Modification
Scalar EnsComposant::getRayonInteractionMin()
{
  Scalar rayonInteractionMin = 1.e10;
  Scalar rayonInteraction;

  vector<Particule*>::iterator particule;
  for (particule=m_ParticuleClassesReference.begin();
  	particule!=m_ParticuleClassesReference.end(); particule++)
  {
    rayonInteraction = (*particule)->getRayonInteraction();
    rayonInteractionMin = rayonInteractionMin > rayonInteraction ?
    	rayonInteraction : rayonInteractionMin;
  }

  return ( rayonInteractionMin );
}




// ----------------------------------------------------------------------------
// Volume total de particules
// G.FERRER - Octo.2003 - Creation
Scalar EnsComposant::getVolume() const
{
  Scalar volume = 0.;
  list<Particule*>::const_iterator particule;
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    volume += (*particule)->getVolume();

  for (particule=m_pwait.begin(); particule!=m_pwait.end(); particule++)
    volume += (*particule)->getVolume();

  return volume;
}




// ----------------------------------------------------------------------------
// Volume total de particules inserees
// G.FERRER - Octo.2003 - Creation
Scalar EnsComposant::getVolumeIn() const
{
  Scalar volume = 0.;
  list<Particule*>::const_iterator particule;
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    volume += (*particule)->getVolume();

  return volume;
}




// ----------------------------------------------------------------------------
// Volume total de particules non-inserees
// G.FERRER - Octo.2003 - Creation
Scalar EnsComposant::getVolumeOut() const
{
  Scalar volume = 0.;
  list<Particule*>::const_iterator particule;
  for (particule=m_pwait.begin(); particule!=m_pwait.end(); particule++)
    volume += (*particule)->getVolume();

  return volume;
}




// ----------------------------------------------------------------------------
// Association des composants avec l'algorithme
// G.FERRER - Aout.2002 - Creation
// G.FERRER - Octo.2003 - Pour tout algorithme
void EnsComposant::Link( AppSec &app )
{
  list<Particule*>::iterator particule;
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    app.Link(*particule);
}




// ----------------------------------------------------------------------------
// Mise a zero de la cinematique si necessaire.
// G.FERRER - Mars.2000 - Creation
void EnsComposant::ResetCinematique(string &reset)
{
  if ( reset == "Reset" )
  {
    list<Particule*>::iterator particule;
    for (particule=m_particulesActives.begin();
    	particule!=m_particulesActives.end(); particule++)
      (*particule)->ResetCinematique();
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Transfert de la particule inactive en attente vers les actives.
// La particule inactive est definie par la fonction getParticule(PullMode)
// G.FERRER - Octo.2003 - Creation
void EnsComposant::ShiftParticuleOutIn()
{
  m_wait->setActivity( COMPUTE );
  removeParticuleFromList( m_pwait, m_wait );
  m_particulesActives.push_back(m_wait);
  if ( m_wait->getTag() == 1 ) m_particulesHalozone.push_back(m_wait);
  else if ( m_wait->getTag() == 2 ) m_particulesClones.push_back(m_wait);
  m_wait = NULL;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cas o� la particule � ins�rer n'est pas dans le sous-domaine:
// on efface le pointeur wait de la liste pwait et on d�truit
// l'objet point� par wait
void EnsComposant::DeleteAndDestroyWait()
{
  removeParticuleFromList( m_pwait, m_wait );
  delete m_wait;
  m_wait = NULL;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Supprime la 1ere instance d'une valeur, si elle existe, dans une
// liste de pointeurs de particule
bool removeParticuleFromList( list<Particule*> &pointerslist, Particule* value )
{
  list<Particule*>::iterator particule;
  bool found =false;
  for (particule=pointerslist.begin();particule!=pointerslist.end() && !found; )
    if ( *particule == value )
    {
      particule = pointerslist.erase( particule );
      found = true;
    }
    else particule++;

  return found;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Supprime l'instance d'une valeur, si elle existe, dans un
// set de pointeurs de particule
bool removeParticuleFromSet( set<Particule*> &pointersSet, Particule* value )
{
  set<Particule*>::iterator particule;
  bool found =false;
  for (particule=pointersSet.begin();particule!=pointersSet.end() && !found; )
    if ( *particule == value )
    {
      pointersSet.erase( particule );
      found = true;
    }
    else particule++;

  return found;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Supprime la 1ere instance d'une valeur, si elle existe, dans une
// liste de pointeurs d'obstacle
bool removeObstacleFromList( list<MonObstacle*> &pointerslist,
	MonObstacle* value )
{
  list<MonObstacle*>::iterator obs;
  bool found =false;
  for (obs=pointerslist.begin();obs!=pointerslist.end() && !found; )
    if (*obs == value)
    {
      obs = pointerslist.erase( obs );
      found = true;
    }
    else obs++;

  return found;
}




// ----------------------------------------------------------------------------
// Mise � jour des particules clones periodiques
// A.WACHS - Sept.2009 - Creation
void EnsComposant::updateClonesPeriodiques( LinkedCell* LC )
{
  if ( !m_particulesClonesPeriodiques.empty() )
  {
    list<Particule*>::iterator particule;
    for (particule=m_particulesClonesPeriodiques.begin();
  	particule!=m_particulesClonesPeriodiques.end();particule++)
    {
      // Suppression du lien avec la la cellule car les clones multi-periodiques
      // ne sont pas mis a jour par LinkUpdate
      if ( LC && (*particule)->getNbPeriodes() > 1 )
        LC->remove( *particule );

      (*particule)->updateVitessePositionPeriodique();
    }
  }

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture des composants pour Reload
// G.FERRER - Janv.2004 - Creation
// A. WACHS - Aout 2014 - Modification
// D. RAKOTONIRINA - Sept 2014 - Modification
void EnsComposant::read( istream &fileSave, string const& filename,
	bool const& new_reload_format )
{
  string buffer, particuleCle, adresse, particuleRefID, materiauCle,
      matType, particuleType, readingMode ;
  int nbreParticules_, nbreParticulesClasses_, ParticuleMPITag;
  Particule *particule;

  // Lecture du nombre de classe de particules
  fileSave >> buffer >> nbreParticulesClasses_;

  // Lecture de tous les elements jusqu'au type de particule puis remise a
  // l'etat initial du flux fileSave. Ceci est du au cas particulier des
  // particules composites
  streampos length = fileSave.tellg();
  fileSave >> particuleCle >> adresse
    >> particuleRefID >> materiauCle >> matType >> particuleType;
  fileSave.seekg( length );

  // Lecture des particule de reference des classes
  // Pour celles-ci, au niveau des formes de particules, on passe par la methode
  // ::create( istream &fileIn ) qui appellent le constructeur standard
  // d'ou NULL au 2nd argument de "read"
  m_ParticuleClassesReference.reserve( nbreParticulesClasses_ );

  for (int i=0; i<nbreParticulesClasses_; i++)
  {
    if ( particuleType == "*CompParticule" )
      particule = new CompParticule( false );
    else particule = new Particule( false );

    particule->read( fileSave, NULL );
    m_ParticuleClassesReference.push_back( particule );
  }

  // Lecture du nb de particule
  fileSave >> readingMode >> nbreParticules_;
  if ( readingMode == "Hybride" ) Grains_Exec::m_writingModeHybrid = true ;

  // Reload des particules
  // Pour celles-ci, au niveau des formes de particules, on passe par le
  // constructeur de recopie (qui pour les polytopes en particule copie les
  // pointeurs plutot que recreer inutilement des tableaux qui sont les
  // memes pour toutes les particules d'une classe)
  // Pour cela, on a besoin de la particule de reference de classe
  // d'ou m_ParticuleClassesReference au 2nd argument de "read"
  ifstream FILEbin;
  if ( Grains_Exec::m_writingModeHybrid )
  {
    string binary_filename = filename + ".bin";
    FILEbin.open( binary_filename.c_str(), ios::in | ios::binary );
  }

  for (int i=0; i<nbreParticules_; i++)
  {
    if ( particule->isCompParticule() ) particule = new CompParticule( false );
    else particule = new Particule( false );

    if ( new_reload_format )
    {
      if ( Grains_Exec::m_writingModeHybrid )
        particule->read2014_binary( FILEbin, &m_ParticuleClassesReference );
      else
        particule->read2014( fileSave, &m_ParticuleClassesReference );
    }
    else particule->read( fileSave, &m_ParticuleClassesReference );
    switch ( particule->getActivity() )
    {
      case COMPUTE:
        m_particulesActives.push_back( particule );
        ParticuleMPITag = particule->getTag();
        switch ( ParticuleMPITag )
        {
          case 1:
            m_particulesHalozone.push_back( particule );
            break;
          case 2:
            m_particulesClones.push_back( particule );
            break;
          default:
            break;
        }
        break;
      default:
        m_pwait.push_back( particule );
        break;
    }
  }
  if ( Grains_Exec::m_writingModeHybrid ) FILEbin.close();

  // Reload de l'arbre des Obstacles
  string stag, nom;
  fileSave >> stag;

  fileSave >> stag >> nom;
  m_obstacle = new CompObstacle( nom );
  fileSave >> stag;
  while ( stag != "</Composite>" )
  {
    Obstacle_BuilderFactory::reload( stag, *m_obstacle, fileSave );
    fileSave >> stag;
  }
  fileSave >> stag;

  assert( stag == "</Obstacle>" );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde des composants pour Reload
// G.FERRER - Janv.2004 - Creation
void EnsComposant::write( ostream &fileSave, string const& filename ) const
{
  fileSave << endl << "nbreParticulesClasses\t" <<
  	m_ParticuleClassesReference.size() << endl;

  vector<Particule*>::const_iterator iv;
  for (iv=m_ParticuleClassesReference.begin();
  	iv!=m_ParticuleClassesReference.end();iv++) (*iv)->write( fileSave );

  // Remarque: ne pas sauvegarder les clones periodiques car ils sont crees
  // au reload (sinon ils existent 2 fois)
  size_t nbActivesNonPer = m_particulesActives.size();
  list<Particule*>::const_iterator particule;

  fileSave << endl << endl << ( Grains_Exec::m_writingModeHybrid ? "Hybride" :
  	"Texte" ) << "\t" << nbActivesNonPer + m_pwait.size() << endl;

  if ( Grains_Exec::m_writingModeHybrid )
  {
    string binary_filename = filename + ".bin";
    ofstream FILEbin( binary_filename.c_str(), ios::out | ios::binary );

    for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
      if ( (*particule)->getID() != -2 )
        (*particule)->write2014_binary( FILEbin );

    for (particule=m_pwait.begin(); particule!=m_pwait.end();
       particule++)
      (*particule)->write2014_binary( FILEbin );

    FILEbin.close();
  }
  else
  {
    for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
      if ( (*particule)->getID() != -2 ) (*particule)->write2014( fileSave );

    for (particule=m_pwait.begin(); particule!=m_pwait.end();
       particule++)
      (*particule)->write2014( fileSave );
  }

  fileSave << "\n<Obstacle>\n";
  m_obstacle->write( fileSave );
  fileSave << "</Obstacle>\n";
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde de l'etat
void EnsComposant::saveState()
{
  list<Particule*>::iterator particule;
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    (*particule)->saveState();

  for (particule=m_pwait.begin(); particule!=m_pwait.end(); particule++)
    (*particule)->saveState();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void EnsComposant::restaureState()
{
  list<Particule*>::iterator particule;
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    (*particule)->restaureState();

  for (particule=m_pwait.begin(); particule!=m_pwait.end(); particule++)
    (*particule)->restaureState();
}




// ----------------------------------------------------------------------------
void EnsComposant::debug( char* s )
{
  cout << s << '\n'
       << "Particules " << m_particulesActives.size() + m_pwait.size()  << '\n'
       << "   Actives " << m_particulesActives.size() << '\t'
       << "   Wait    " << m_pwait.size()      << endl;
}




// ----------------------------------------------------------------------------
// Operateur <<
ostream& operator << ( ostream &f, const EnsComposant &EC )
{
  f << "Nombre total de particules sur tous les processeurs = "
  	<< EC.m_nb_total_particules << endl;
  f << "Nombre de particules sur le processeur = " <<
  	EC.m_particulesActives.size() + EC.m_pwait.size() << endl;
  f << "Nombre de particules actives = " << EC.m_particulesActives.size()
  	<< endl;
  f << "Nombre de particules en attente = " << EC.m_pwait.size() << endl;
  f << "Nombre de particules dans la zone de recouvrement = " <<
  	EC.m_particulesHalozone.size() << endl;
  f << "Nombre de clones = " << EC.m_particulesClones.size() << endl;
  list<Particule*>::const_iterator il;
  for (il=EC.m_particulesActives.begin();il!=EC.m_particulesActives.end();il++)
    f << *(*il) << endl;

  return f;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture des composants pour Post-Processing (Ecriture de reference)
// G.FERRER - Janv.2004 - Creation
// A.WACHS - Fev.2010 - Modification
// D. RAKOTONIRINA - Juil.2014 - Modification
void EnsComposant::PostProcessing_start( Scalar temps, Scalar dt,
	LinkedCell const* LC, vector<Fenetre> const& insert_windows,
	int rang, int nprocs,
	MPIWrapperGrains const* wrapper,
	size_t indent_width )
{
  list<Particule*>* postProcessingParticules = NULL;
  list<Particule*>* postProcessingWait = NULL;
  list<Particule*>* postProcessingPeriodiques = NULL;
  vector<Particule*>* particulespost = NULL;
  list<PostProcessingWriter*>::iterator pp;
  bool written = false ;
  string const siw( indent_width, ' ' ) ;

  if ( rang == 0 )
    for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
      if ( !(*pp)->isCompFeaturesWriter() && !written )
      {
	cout << siw << "Sortie resultats: START" << endl;
	written = true;
      }

  // Dans le cas d'une simulation periodique ou la periodicite est geree
  // par le pattern MPI, il faut determiner les clones periodiques
  // Ce sont ceux qui sont hors du domaine de calcul tel que defini dans App
  // soit sans les cellules supplementaires aux extremites pour gerer la
  // periodicite
  if ( Grains_Exec::m_MPIperiodique )
    for (list<Particule*>::iterator  particule=m_particulesClones.begin();
  	particule!=m_particulesClones.end(); particule++)
      if ( !App::isInDomain( (*particule)->getPosition() ) )
        m_particulesClonesPeriodiques.push_back( *particule );

  // Est ce que un des post-processors est sequentiel
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end()
    && !m_hasSerialPostProcessors;pp++)
  {
    if ( !(*pp)->isParallelWriter() ) m_hasSerialPostProcessors = true;
  }

  if ( nprocs > 1 && m_hasSerialPostProcessors )
  {
    if ( rang == 0 )
      cout << siw << "Copie des particules sur le master pour Post-processing"
      	<< endl;

    // Collecte les particules de tous les processeurs dans un vecteur
    // sur le master
    particulespost=wrapper->GatherParticules_PostProcessing(
  	m_particulesActives, m_pwait, m_ParticuleClassesReference,
	m_nb_total_particules );

    // Creation de la liste
    if ( particulespost )
    {
      postProcessingParticules = new list<Particule*>;
      vector<Particule*>::iterator iv;
      for (iv=particulespost->begin();iv!=particulespost->end();iv++)
        postProcessingParticules->push_back(*iv);
    }

    if ( Grains_Exec::m_periodique || Grains_Exec::m_MPIperiodique )
      // Collecte les clones periodiques de tous les processeurs dans une liste
      // sur le master
      postProcessingPeriodiques =
      	wrapper->GatherClonesPeriodiques_PostProcessing(
  	m_particulesClonesPeriodiques, m_ParticuleClassesReference );
    else
      postProcessingPeriodiques = new list<Particule*>;
  }

  // Initialise force at the contact point for the post-processing
  // of the stress tensor
  if ( Grains_Exec::m_stressTensor || Grains_Exec::m_particleStressTensor )
    InitializeForcesAtContactPoint();

  // Post processing writers
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    if (nprocs > 1 && m_hasSerialPostProcessors && !(*pp)->isParallelWriter())
    {
      (*pp)->PostProcessing_start( temps,dt,
  	postProcessingParticules,
	postProcessingWait,
	postProcessingPeriodiques,
	&m_ParticuleClassesReference,
	m_obstacle,
	LC,
	insert_windows );
    }
    else
    {
      (*pp)->PostProcessing_start( temps,dt,
  	&m_particulesActives,
	&m_pwait,
	&m_particulesClonesPeriodiques,
	&m_ParticuleClassesReference,
	m_obstacle,
	LC,
	insert_windows );
    }

  if ( nprocs > 1 && m_hasSerialPostProcessors )
  {
    // Destruction du vecteur
    if ( particulespost )
    {
      vector<Particule*>::iterator iv;
      for (iv=particulespost->begin();iv!=particulespost->end();iv++)
        delete *iv;
      particulespost->clear();
      delete particulespost;
      postProcessingParticules->clear();
      delete postProcessingParticules;
    }

    // Destruction de la liste
    if ( postProcessingPeriodiques )
    {
      list<Particule*>::iterator il;
      for (il=postProcessingPeriodiques->begin();
      	il!=postProcessingPeriodiques->end();il++)
        delete *il;
      postProcessingPeriodiques->clear();
      delete postProcessingPeriodiques;
    }
  }

  // Dans le cas d'une simulation periodique ou la periodicite est geree
  // par le pattern MPI, on vide la liste
  if ( Grains_Exec::m_MPIperiodique ) m_particulesClonesPeriodiques.clear();

  written = false ;
  if ( rang == 0 )
    for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
      if ( !(*pp)->isCompFeaturesWriter() && !written )
      {
	cout << siw << "Sortie resultats: COMPLETED" << endl;
	written = true;
      }
  if ( rang == 0 && Grains_Exec::m_ContactforceOutput )
    m_initial_time = temps;

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture des composants pour Post-Processing (Ecriture d'evolution)
// G.FERRER - Janv.2004 - Creation
// A.WACHS - Fev.2010 - Modification
void EnsComposant::PostProcessing( Scalar temps, Scalar dt,
	LinkedCell const* LC, int rang,
	int nprocs, MPIWrapperGrains const* wrapper,
	size_t indent_width )
{
  list<Particule*>* postProcessingParticules = NULL;
  list<Particule*>* postProcessingWait = NULL;
  list<Particule*>* postProcessingPeriodiques = NULL;
  vector<Particule*>* particulespost = NULL;
  string const siw( indent_width, ' ' ) ;

  if ( rang == 0 )
    cout << siw << "Sortie resultats: START" << endl;

  // Dans le cas d'une simulation periodique ou la periodicite est geree
  // par le pattern MPI, il faut determiner les clones periodiques
  // Ce sont ceux qui sont hors du domaine de calcul tel que defini dans App
  // soit sans les cellules supplementaires aux extremites pour gerer la
  // periodicite
  if ( Grains_Exec::m_MPIperiodique )
    for (list<Particule*>::iterator  particule=m_particulesClones.begin();
  	particule!=m_particulesClones.end(); particule++)
      if ( !App::isInDomain( (*particule)->getPosition() ) )
        m_particulesClonesPeriodiques.push_back( *particule );

  if ( nprocs > 1 && m_hasSerialPostProcessors )
  {
    if ( rang == 0 )
      cout << siw << "Copie des particules sur le master pour Post-processing"
      	<< endl;

    // Collecte les particules de tous les processeurs dans un vecteur
    // sur le master
    particulespost = wrapper->GatherParticules_PostProcessing(
  	m_particulesActives, m_pwait, m_ParticuleClassesReference,
	m_nb_total_particules );

    // Creation de la liste
    if ( particulespost )
    {
      postProcessingParticules = new list<Particule*>;
      vector<Particule*>::iterator iv;
      for (iv=particulespost->begin();iv!=particulespost->end();iv++)
        postProcessingParticules->push_back(*iv);
    }

    if ( Grains_Exec::m_periodique || Grains_Exec::m_MPIperiodique )
      // Collecte les clones periodiques de tous les processeurs dans une liste
      // sur le master
      postProcessingPeriodiques =
      	wrapper->GatherClonesPeriodiques_PostProcessing(
  	m_particulesClonesPeriodiques, m_ParticuleClassesReference );
    else
      postProcessingPeriodiques = new list<Particule*>;
  }

  if ( Grains_Exec::m_ContactforceOutput )
  {
    IntegrateContactPostProcessingForces((temps - m_initial_time)/dt);
    m_initial_time = temps;
  }

  // Post processing writers
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    if ( nprocs > 1 && m_hasSerialPostProcessors && !(*pp)->isParallelWriter() )
    {
      (*pp)->PostProcessing( temps,dt,
  	postProcessingParticules,
	postProcessingWait,
	postProcessingPeriodiques,
	&m_ParticuleClassesReference,
	m_obstacle,
	LC );
    }
    else
    {
      (*pp)->PostProcessing( temps,dt,
  	&m_particulesActives,
	&m_pwait,
	&m_particulesClonesPeriodiques,
	&m_ParticuleClassesReference,
	m_obstacle,
	LC );
    }

  if ( nprocs > 1 && m_hasSerialPostProcessors )
  {
    // Destruction du vecteur
    if ( particulespost )
    {
      vector<Particule*>::iterator iv;
      for (iv=particulespost->begin();iv!=particulespost->end();iv++)
        delete *iv;
      particulespost->clear();
      delete particulespost;
      postProcessingParticules->clear();
      delete postProcessingParticules;
    }

    // Destruction de la liste
    if ( postProcessingPeriodiques )
    {
      list<Particule*>::iterator il;
      for (il=postProcessingPeriodiques->begin();
      	il!=postProcessingPeriodiques->end();il++)
        delete *il;
      postProcessingPeriodiques->clear();
      delete postProcessingPeriodiques;
    }
  }

  // if Contact force is written as output, the PP force should be initialized
  // at the end of each output saving
  if ( Grains_Exec::m_ContactforceOutput || Grains_Exec::m_withlubrication ||
            Grains_Exec::m_ContactforceOutput_instantaneous  )
     InitializePostProcessingForces();

  // Initialise force at the contact point for the post-processing
  // of the stress tensor
  if ( Grains_Exec::m_stressTensor || Grains_Exec::m_particleStressTensor )
    InitializeForcesAtContactPoint();

  // Dans le cas d'une simulation periodique ou la periodicite est geree
  // par le pattern MPI, on vide la liste
  if ( Grains_Exec::m_MPIperiodique ) m_particulesClonesPeriodiques.clear();

  if ( rang == 0 )
    cout << siw << "Sortie resultats: COMPLETED" << endl;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Fermeture du Post-Processing
// G.FERRER - Janv.2004 - Creation
// A.WACHS - Fev.2010 - Modification
void EnsComposant::PostProcessing_end()
{
  // Post processing writers
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    (*pp)->PostProcessing_end();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture des composants pour Post-Processing d'une erreur de contact
void EnsComposant::PostProcessingErreurComposants( string const& filename,
	list<Composant*> const& errcomposants )
{
  // Post processing writers
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    (*pp)->writeErreurComposantsPostProcessing( filename, errcomposants );
}




// ----------------------------------------------------------------------------
// Renvoie les vitesses max et moyenne de l'ensemble des composants
// A.WACHS - Sept 2009 - Creation
void EnsComposant::ComputeMaxMeanVelocity( Scalar &vmax, Scalar &vmean,
  	MPIWrapperGrains const* wrapper ) const
{
  Scalar vit;
  vmax = vmean = 0. ;
  size_t ncomp = 0 ;

  // Particules
  list<Particule*>::const_iterator particule;
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
  {
    vit = Norm( *(*particule)->getVitesseTranslation() );
    vmax = vit > vmax ? vit : vmax;
    vmean += vit;
    ++ncomp;
  }

  // Obstacles
  list<MonObstacle*> obstacles = m_obstacle->getObstacles();
  list<MonObstacle*>::const_iterator myObs;
  for (myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++)
  {
    if ( (*myObs)->hasMoved() )
    {
      vit = Norm( *(*myObs)->getVitesseTranslation() );
      vmax = vit > vmax ? vit : vmax;
      vmean += vit;
      ++ncomp;
    }
  }

  if ( wrapper )
  {
    vmax = wrapper->max_DOUBLE_master( vmax );
    ncomp = wrapper->sum_UNSIGNED_INT_master( ncomp );
    vmean = wrapper->sum_DOUBLE_master( vmean );
  }

  if ( ncomp ) vmean /= double(ncomp) ;
}


// ----------------------------------------------------------------------------
// Ecrit dans un fichier les vitesses min, max et moy des particles en
// fonction du temps
// A.WACHS - Dec 2014 - Creation
void EnsComposant::monitorParticlesVelocity( Scalar temps, ofstream& fileOut,
  	int rang, MPIWrapperGrains const* wrapper ) const
{
  Vecteur vmin( 1.20 ), vmax( -1.e20 ), vmean;
  size_t ncomp = m_particulesActives.size() ;
  Vecteur const* vtrans = NULL ;
  Scalar vit = 0. ;

  // Particules
  list<Particule*>::const_iterator particule;
  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
  {
    vtrans = (*particule)->getVitesseTranslation();
    for (int i=0;i<3;++i)
    {
      vit = (*vtrans)[i];
      vmin[i] = vmin[i] < vit ? vmin[i] : vit;
      vmax[i] = vmax[i] > vit ? vmax[i] : vit;
      vmean[i] += vit;
    }
  }

  if ( wrapper )
  {
    ncomp = wrapper->sum_UNSIGNED_INT_master( ncomp );
    for (int i=0;i<3;++i)
    {
      vmin[i] = wrapper->min_DOUBLE_master( vmin[i] );
      vmax[i] = wrapper->max_DOUBLE_master( vmax[i] );
      vmean[i] = wrapper->sum_DOUBLE_master( vmean[i] );
    }
  }

  if ( rang == 0 )
  {
    fileOut << temps;
    for (int i=0;i<3;++i)
    {
      vmean[i] /= double(ncomp) ;
      fileOut << " " << vmin[i] << " " << vmax[i] << " " << vmean[i];
    }
    fileOut << endl;
  }
}




// ----------------------------------------------------------------------------
// Mise � jour de la localisation geographique des particules de
// la zone de recouvrement
// A.WACHS - Oct 2009 - Creation
void EnsComposant::updateGeoLocalisationParticulesHalozone()
{
  list<Particule*>::iterator particule;
  for (particule=m_particulesHalozone.begin();
  	particule!=m_particulesHalozone.end();particule++)
    (*particule)->updateGeoLocalisation();
}




// ----------------------------------------------------------------------------
// Ajoute un post-processing writer � la liste
void EnsComposant::addPostProcessingWriter( PostProcessingWriter* ppw )
{
  m_postProcessors.push_back(ppw);
}




// ----------------------------------------------------------------------------
// Numero de cycle initial
void EnsComposant::setInitialCycleNumber( const int& cycle0 )
{
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    (*pp)->setInitialCycleNumber( cycle0 );
}




// ----------------------------------------------------------------------------
// Verifie que le post processing Paraview est actif, sinon le cree
void EnsComposant::checkParaviewPostProcessing(
	int const& rang,
  	int const& nprocs,
	const string &name_,
	const string &root_,
  	const bool &isBinary )
{
  PostProcessingWriter *ParaviewPP = NULL;
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    if ( (*pp)->getPostProcessingWriterType() == "Paraview" )
      ParaviewPP = *pp;

  // Si aucun post processing writer n'existe, initialiser le vecteur
  // des fenetres de post processing a true sur tous les procs
  if ( m_postProcessors.empty() )
    PostProcessingWriter::allocate_PostProcessingWindow( nprocs );

  // Si un post processing writer Paraview existe deja, on le detruit
  // et on le recreer avec les bons parametres
  if ( ParaviewPP )
  {
    for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
      if ( (*pp) == ParaviewPP )
      {
        delete *pp;
	pp = m_postProcessors.erase( pp );
      }
      else pp++;
    ParaviewPP = NULL;
  }

  // Creation du post processing writer Paraview avec les bons parametres
  ParaviewPP = new Paraview_PostProcessingWriter( rang, nprocs, name_, root_,
  	isBinary );
  m_postProcessors.push_back( ParaviewPP );
}




// ----------------------------------------------------------------------------
// Verifie que le post processing Matlab est actif, sinon le cree
void EnsComposant::checkMatlabPostProcessing(
	int const& rang,
  	int const& nprocs,
	const string &name_,
	const string &root_,
  	const bool &isBinary )
{
  PostProcessingWriter *MatlabPP = NULL;
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    if ( (*pp)->getPostProcessingWriterType() == "Matlab" )
      MatlabPP = *pp;

  // Si aucun post processing writer n'existe, initialiser le vecteur
  // des fenetres de post processing a true sur tous les procs
  if ( m_postProcessors.empty() )
    PostProcessingWriter::allocate_PostProcessingWindow( nprocs );

  // Si un post processing writer Matlab existe deja, on le detruit
  // et on le recreer avec les bons parametres
  if ( MatlabPP )
  {
    for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
      if ( (*pp) == MatlabPP )
      {
        delete *pp;
	pp = m_postProcessors.erase( pp );
      }
      else pp++;
    MatlabPP = NULL;
  }

  // Creation du post processing writer Matlab avec les bons parametres
  MatlabPP = new Matlab_PostProcessingWriter( rang, nprocs, name_, root_,
  	isBinary );
  m_postProcessors.push_back( MatlabPP );
}




// ----------------------------------------------------------------------------
// Numero maximum de particules
// A.WACHS - Mars.2010 - Creation
int EnsComposant::numeroMaxParticules() const
{
  int numeroMax = 0;
  list<Particule*>::const_iterator particule;

  for (particule=m_particulesActives.begin();
  	particule!=m_particulesActives.end(); particule++)
    numeroMax = numeroMax < (*particule)->getID() ?
    	(*particule)->getID() : numeroMax;

  for (particule=m_pwait.begin(); particule!=m_pwait.end(); particule++)
    numeroMax = numeroMax < (*particule)->getID() ?
    	(*particule)->getID() : numeroMax;

  return numeroMax;
}




// ----------------------------------------------------------------------------
// Fr�quence de mise � jour du lien entre les obstacles et le LinkedCell
// A.WACHS - Nov.2010 - Creation
void EnsComposant::setObstaclesLinkedCellUpdateFrequency(
	int const &updateFreq )
{
  list<MonObstacle*> obstacles = m_obstacle->getObstacles();
  list<MonObstacle*>::iterator myObs;
  for (myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++)
    (*myObs)->setObstacleLinkedCellUpdateFrequency( updateFreq );
}




// ----------------------------------------------------------------------------
// Mise � jour du statut des particules clones periodiques en parallele
// A.WACHS - Janv.2011 - Creation
void EnsComposant::statutClonesPeriodiques_MPI_Step3( Scalar temps,
  	list<int>& ClonestoDestroy,
  	list<int>& ClonestoParticules,
	LinkedCell* LC )
{
  list<Particule*>::iterator il;
  list<int>::iterator ilD,ilS;
  int refID = 0;
  Particule* pdestroy = NULL;
  bool b_EraseClonesFromList = false;

  for (il=m_particulesClonesPeriodiques.begin();
  	il!=m_particulesClonesPeriodiques.end(); )
  {
    refID = (*il)->getPeriodicReferenceID();
    b_EraseClonesFromList = false;

    // Clones � detruire
    for (ilD=ClonestoDestroy.begin();ilD!=ClonestoDestroy.end()
    	&& !b_EraseClonesFromList; )
    {
      if ( refID == *ilD )
      {
	if ( Grains_Exec::m_MPI_verbose )
	{
	  ostringstream oss;
	  Point const* gc = (*il)->getPosition();
          oss << "   t=" << Grains_Exec::doubleToString( temps, TIMEFORMAT )
      		<< " Periodic clone out of domain COM            Id = " <<
      		(*il)->getPeriodicReferenceID() << "  " <<
		(*gc)[X] << " " << (*gc)[Y] << " " << (*gc)[Z] << endl;
      	  MPIWrapperGrains::addToMPIString( oss.str() );
	}
	pdestroy = *il;

        // Suppression de la particule dans le LinkedCell
        LC->remove( pdestroy );

        // Destruction de l'objet point�
        delete pdestroy;

        // Pour suppression de la liste des clones
        b_EraseClonesFromList = true;
	ilD = ClonestoDestroy.erase( ilD );
      }
      else ilD++;
    }

    // Clones dont le statut change en particule classique
    if (!b_EraseClonesFromList)
    {
      for (ilS=ClonestoParticules.begin();ilS!=ClonestoParticules.end()
    	&& !b_EraseClonesFromList; )
      {
        if ( refID == *ilS )
        {
	  if ( Grains_Exec::m_MPI_verbose )
	  {
	    ostringstream oss;
	    Point const* gc = (*il)->getPosition();
            oss << "   t=" << Grains_Exec::doubleToString( temps, TIMEFORMAT )
      		<< " Periodic clone -> particule COM             Id = " <<
      		(*il)->getPeriodicReferenceID() << "  " <<
		(*gc)[X] << " " << (*gc)[Y] << " " << (*gc)[Z] << endl;
      	    MPIWrapperGrains::addToMPIString( oss.str() );
	  }
	  pdestroy = *il;

	  // Creer une nouvelle particule active
	  Particule* pa = new Particule( pdestroy->getPeriodicReferenceID(),
	    	m_ParticuleClassesReference[pdestroy->getParticuleClasse()],
		*(pdestroy->getVitesseTranslation()),
		*(pdestroy->getCinematique()->getRotation()),
		*(pdestroy->getVitesseRotation()),
		*(pdestroy->getForme()->getTransform()),
		COMPUTE,
		0 );

          // Suppression de la particule dans le LinkedCell
          LC->remove( pdestroy );

          // Ajout de la nouvelle particule active dans le LinkedCell
          LC->Link( pa );

          // Ajout de la particule active dans les differentes listes
	  m_particulesActives.push_back(pa);
	  if (pa->getTag() == 1) m_particulesHalozone.push_back(pa);
	  else if (pdestroy->getTag() == 2) m_particulesClones.push_back(pa);

          // Destruction du clone
          delete pdestroy;

          // Pour suppression de la liste des clones
          b_EraseClonesFromList = true;
	  ilS = ClonestoParticules.erase( ilS );
        }
        else ilS++;
      }
    }

    // Bilan: suppresion ou non
    if ( b_EraseClonesFromList ) il = m_particulesClonesPeriodiques.erase( il );
    else il++;
  }
}




// ----------------------------------------------------------------------------
// Mise � jour du statut des particules clones paralleles qui sont
// egalement des references periodiques
// A.WACHS - Janv.2011 - Creation
void EnsComposant::statutClonesReferencesPeriodiques_MPI_Step2( Scalar temps,
  	list<int>& PartRefPerHalozone,
  	list<int>& PartRefPerOutDomainHalozone,
  	list<int>& InNotRefPerHalozone,
	LinkedCell* LC )
{
  list<Particule*>::iterator il;
  list<int>::iterator ilP,ilO,ilN;
  int refID = 0;
  Particule* pdestroy = NULL;
  bool found = false, erase = false;

  for (il=m_particulesClones.begin();il!=m_particulesClones.end(); )
  {
    refID = (*il)->getID();
    erase = false;
    found = false;

    // Nouvelle reference periodique
    for (ilP=PartRefPerHalozone.begin();ilP!=PartRefPerHalozone.end() &&
    	!found; )
    {
      if ( refID == *ilP )
      {
        if ( !(*il)->getNombreClonesPeriodiques() )
        {
          // Ajout dans la liste des particules de reference periodiques
	  m_particulesReferencesPeriodiques.insert(*il);
	  ilP = PartRefPerHalozone.erase( ilP );

	  // Ajout dans la liste des clones de la particule de reference
	  (*il)->addPeriodicObstacleID( *ilP );
	  ilP = PartRefPerHalozone.erase( ilP );

          if ( Grains_Exec::m_MPI_verbose )
	  {
	    ostringstream oss;
	    Point const* gc = (*il)->getPosition();
            oss << "   t=" << Grains_Exec::doubleToString( temps, TIMEFORMAT )
      		<< " New periodic reference COM                  Id = " <<
      		refID << "  " <<
		(*gc)[X] << " " << (*gc)[Y] << " " << (*gc)[Z] << endl;
      	    MPIWrapperGrains::addToMPIString( oss.str() );
	  }
        }
	found = true;
      }
      else {ilP++;ilP++;}
    }

    // Particule hors du domaine
    if ( !found )
      for (ilO=PartRefPerOutDomainHalozone.begin();
    	ilO!=PartRefPerOutDomainHalozone.end() && !found; )
      {
        if ( refID == *ilO )
        {
          if ( Grains_Exec::m_MPI_verbose )
	  {
	    ostringstream oss;
            Point const* gc = (*il)->getPosition();
            oss << "   t=" << Grains_Exec::doubleToString( temps, TIMEFORMAT )
      		<< " Periodic reference out of domain COM        Id = " <<
      		refID << "  " <<
		(*gc)[X] << " " << (*gc)[Y] << " " << (*gc)[Z] << endl;
      	    MPIWrapperGrains::addToMPIString( oss.str() );
	  }
	  pdestroy = *il;

          // Suppression de la particule du LinkedCell
          LC->remove(pdestroy);

          // Suppression des differentes listes
          removeParticuleFromSet( m_particulesReferencesPeriodiques, pdestroy );
	  removeParticuleFromList( m_particulesActives, pdestroy );
	  erase=true;

          // Destruction de l'objet point�
          delete pdestroy;

	  ilO = PartRefPerOutDomainHalozone.erase( ilO );
	  found = true;
        }
        else ilO++;
      }

    // Particules ayant perdu leur statut de reference periodique
    if ( !found )
      for (ilN=InNotRefPerHalozone.begin();
    	ilN!=InNotRefPerHalozone.end() && !found; )
      {
        if ( refID == *ilN )
        {
	  if ( Grains_Exec::m_MPI_verbose )
	  {
	    ostringstream oss;
	    Point const* gc = (*il)->getPosition();
            oss << "   t="
	      	<< Grains_Exec::doubleToString( temps, TIMEFORMAT )
      		<< " Particule is not periodic reference anymore COM Id = " <<
      		refID << "  " <<
		(*gc)[X] << " " << (*gc)[Y] << " " << (*gc)[Z] << endl;
      	    MPIWrapperGrains::addToMPIString( oss.str() );
          }
          // supprimer le numero de l'obstacle periodique
	  ilN = InNotRefPerHalozone.erase( ilN );
          (*il)->erasePeriodicObstacleID( *ilN );

          // Supprimer la particule de la liste des particules de reference
          removeParticuleFromSet( m_particulesReferencesPeriodiques, *il );

	  ilN = InNotRefPerHalozone.erase( ilN );
	  found = true;
        }
        else {ilN++;ilN++;}
      }

    if ( erase ) il = m_particulesClones.erase( il );
    else il++;
  }
}




// ----------------------------------------------------------------------------
// Cree une sauvegarde de l'etat des obstacles
// A.WACHS - Fev.2012 - Creation
void EnsComposant::createStateObstacles( list<struct ObstacleState*>
	&obsStates ) const
{
  m_obstacle->createState( obsStates );
}




// ----------------------------------------------------------------------------
// Restaure l'etat des obstacles
// A.WACHS - Fev.2012 - Creation
void EnsComposant::restaureStateObstacles( list<struct ObstacleState*>
	&obsStates )
{
  m_obstacle->restaureState( obsStates );
}




// ----------------------------------------------------------------------------
// Parametres de post-processing des efforts sur les obstacles
// A.WACHS - Mai.2012 - Creation
void EnsComposant::setOutputObstaclesLoadParameters( string const& root_,
  	int const& freq_,
	list<string> const& ObsNames )
{
  m_outputTorseurObstacles_dir = root_;
  m_outputTorseurObstacles_frequency = freq_;

  for (list<string>::const_iterator il=ObsNames.begin();il!=ObsNames.end();il++)
  {
    Obstacle* pobs = const_cast<Obstacle*>(m_obstacle->getNom( *il ));
    if ( pobs )
      m_outputTorseurObstacles.push_back( pobs );
  }
}




// ----------------------------------------------------------------------------
// Post-processing des efforts sur les obstacles
// A.WACHS - Mai.2012 - Creation
void EnsComposant::outputObstaclesLoad( Scalar temps, Scalar dt,
	bool enforceOutput, bool increaseCounterOnly, int rang, int nprocs,
	MPIWrapperGrains const* wrapper )
{
  if ( !increaseCounterOnly )
  {
    if ( m_outputTorseurObstacles_counter == 0 || enforceOutput )
    {
      if ( nprocs > 1 )
        wrapper->sumObstaclesLoad( m_obstacle->getObstacles() );

      if ( rang == 0 )
      {
        Torseur const* torseur = NULL;
        Vecteur const* force = NULL;
        Vecteur const* torque = NULL;
        for (list<Obstacle*>::iterator
    		obstacle=m_outputTorseurObstacles.begin();
		obstacle!=m_outputTorseurObstacles.end();obstacle++)
        {
          ofstream OUT( ( m_outputTorseurObstacles_dir
      	+ "/Loading_" + (*obstacle)->getName() + ".res" ).c_str(), ios::app );
          torseur = (*obstacle)->getTorseur();
	  force = torseur->getForce();
          torque = torseur->getMoment();
	  OUT << temps << " " <<
		Grains_Exec::doubleToString( ios::scientific, 6, (*force)[X] )
		<< " " <<
		Grains_Exec::doubleToString( ios::scientific, 6, (*force)[Y] )
		<< " " <<
		Grains_Exec::doubleToString( ios::scientific, 6, (*force)[Z] )
		<< " " <<
		Grains_Exec::doubleToString( ios::scientific, 6, (*torque)[X] )
		<< " " <<
		Grains_Exec::doubleToString( ios::scientific, 6, (*torque)[Y] )
		<< " " <<
		Grains_Exec::doubleToString( ios::scientific, 6, (*torque)[Z] )
		<< " " << endl;
          OUT.close();
        }
      }
    }
  }

  if ( !enforceOutput )
  {
    ++m_outputTorseurObstacles_counter;
    if ( m_outputTorseurObstacles_counter ==
    	m_outputTorseurObstacles_frequency )
    m_outputTorseurObstacles_counter = 0 ;
  }
}




// ----------------------------------------------------------------------------
// Initialise les fichiers de sortie des efforts sur les obstacles
void EnsComposant::initialiseOutputObstaclesLoadFiles( int rang,
	bool coupledFluid, Scalar temps )
{
  m_outputTorseurObstacles_counter = coupledFluid ;

  if ( rang == 0 )
  {
    if ( Grains_Exec::m_ReloadType == "new" )
    {
      string cmd = "bash " + Grains_Exec::m_GRAINS_HOME
     	+ "/ExecScripts/ObstaclesLoadFiles_clear.exec "
	+ m_outputTorseurObstacles_dir;
      system( cmd.c_str() );
    }
    else
      for (list<Obstacle*>::iterator
    	obstacle=m_outputTorseurObstacles.begin();
	obstacle!=m_outputTorseurObstacles.end();obstacle++)
        Grains_Exec::checkTime_outputFile( m_outputTorseurObstacles_dir
      		+ "/Loading_" + (*obstacle)->getName() + ".res",
    		temps ) ;
  }

}




// ----------------------------------------------------------------------------
// Nombre total de particules sur tous les procs
void EnsComposant::setNbreParticulesOnAllProc( const size_t &nb_ )
{
  m_nb_total_particules = int(nb_);
  Grains_Exec::setNbreParticulesOnAllProc( m_nb_total_particules );
}




// ----------------------------------------------------------------------------
// Initialisation de la vitesse translationnelle et rotationnelle des
// particules au pas de temps precedent en cas de restart (lecture dans les
// fichiers resultats generes par PeliGRIFF, solution temporaire)
void EnsComposant::setVelocityAndVelocityDifferencePreviousTimeRestart(
	string const& PelDirRes )
{
  string linet, linetnm1, sbuffer ;
  istringstream iss;
  double time_ ;
  vector<double> vxnm1( m_nb_total_particules, 0. ),
  	vynm1( m_nb_total_particules, 0. ),
	vznm1( m_nb_total_particules, 0. ),
	omxnm1( m_nb_total_particules, 0. ),
  	omynm1( m_nb_total_particules, 0. ),
	omznm1( m_nb_total_particules, 0. ) ;
  int i = 0 ;
  list<Particule*>::iterator particule;

  // Lecture des vitesse au temps t^n-1 dans les fichiers de resultat PeliGRIFF
  // !!! Rem: cela suppose que l'ecriture a ete realisee a chaque pas de temps
  // fluide !!!
  // Vitesse translationnelle
  // Composante x
  ifstream vtx( ( PelDirRes + "/x-velocity_time.res" ).c_str(), ios::in );
  while ( !vtx.eof() )
  {
    getline( vtx, sbuffer, '\n' );
    if ( !vtx.eof() )
    {
      linetnm1 = linet ;
      linet = sbuffer ;
    }
  }
  vtx.close() ;

  iss.str( linetnm1 );
  iss >> time_ ;
  for (size_t j=0;j<m_nb_total_particules;++j) iss >> vxnm1[j];
  iss.clear();

  // Composante y
  ifstream vty( ( PelDirRes + "/y-velocity_time.res" ).c_str(), ios::in );
  while ( !vty.eof() )
  {
    getline( vty, sbuffer, '\n' );
    if ( !vty.eof() )
    {
      linetnm1 = linet ;
      linet = sbuffer ;
    }
  }
  vty.close() ;

  iss.str( linetnm1 );
  iss >> time_ ;
  for (size_t j=0;j<m_nb_total_particules;++j) iss >> vynm1[j];
  iss.clear();

  // Vitesse rotationnelle
  // Composante z
  ifstream omz( ( PelDirRes + "/z-omega_time.res" ).c_str(), ios::in );
  while ( !omz.eof() )
  {
    getline( omz, sbuffer, '\n' );
    if ( !omz.eof() )
    {
      linetnm1 = linet ;
      linet = sbuffer ;
    }
  }
  omz.close() ;

  iss.str( linetnm1 );
  iss >> time_ ;
  for (size_t j=0;j<m_nb_total_particules;++j) iss >> omznm1[j];
  iss.clear();

  if ( Grains_BuilderFactory::getContext() == DIM_3 )
  {
    // Vitesse translationnelle
    // Composante z
    ifstream vtz( ( PelDirRes + "/z-velocity_time.res" ).c_str(), ios::in );
    while ( !vtz.eof() )
    {
      getline( vtz, sbuffer, '\n' );
      if ( !vtz.eof() )
      {
        linetnm1 = linet ;
        linet = sbuffer ;
      }
    }
    vtz.close() ;

    iss.str( linetnm1 );
    iss >> time_ ;
    for (size_t j=0;j<m_nb_total_particules;++j) iss >> vznm1[j];
    iss.clear();

    // Vitesse rotationnelle
    // Composante x
    ifstream omx( ( PelDirRes + "/x-omega_time.res" ).c_str(), ios::in );
    while ( !omx.eof() )
    {
      getline( omx, sbuffer, '\n' );
      if ( !omx.eof() )
      {
        linetnm1 = linet ;
        linet = sbuffer ;
      }
    }
    omx.close() ;

    iss.str( linetnm1 );
    iss >> time_ ;
    for (size_t j=0;j<m_nb_total_particules;++j) iss >> omxnm1[j];
    iss.clear();

    // Composante y
    ifstream omy( ( PelDirRes + "/y-omega_time.res" ).c_str(), ios::in );
    while ( !omy.eof() )
    {
      getline( omy, sbuffer, '\n' );
      if ( !omy.eof() )
      {
        linetnm1 = linet ;
        linet = sbuffer ;
      }
    }
    omy.close() ;

    iss.str( linetnm1 );
    iss >> time_ ;
    for (size_t j=0;j<m_nb_total_particules;++j) iss >> omynm1[j];
    iss.clear();
  }

  // Affectation des vitesses a t^n-1 a chaque particule
  // et initialisation des differences de vitesse
  // On boucle sur les particules sans tenir compte de leur numero
  // comme dans InterfaceFluide car lorsque l'insertion des particules
  // s'est faite en "Aleatoire", les numeros cote Grains3D et PeliGRIFF
  // ne correspondent pas
  for (particule=m_particulesActives.begin(),i=0;
  	particule!=m_particulesActives.end(); particule++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
    {
      (*particule)->setVelocityPreviousTimeRestart( vxnm1[i],
  	vynm1[i], vznm1[i], omxnm1[i], omynm1[i], omznm1[i] ) ;
      (*particule)->setVelocityAndVelocityDifferencePreviousTime() ;
      ++i;
    }
  }
}

double EnsComposant::
	setVelocityAndVelocityDifferencePreviousTimeRestart_Basilisk(
	string const& rootfilename )
{
  string linet, linetnm1, sbuffer ;
  istringstream iss;
  double time_n, time_nm1, previousdtfluid = 1, vxnm1, vynm1, vznm1, omxnm1,
  	omynm1, omznm1;
  int i = 0 ;
  list<Particule*>::iterator particule;

  for (particule=m_particulesActives.begin(),i=0;
  	particule!=m_particulesActives.end(); particule++)
  {
    if ( (*particule)->getActivity() == COMPUTE
    	&& (*particule)->getID() >= 0 )
    {
      ifstream pdata( ( rootfilename + "_" + Grains_Exec::intToString(i)
      	+ ".dat" ).c_str(), ios::in );

      while ( !pdata.eof() )
      {
        getline( pdata, sbuffer, '\n' );
        if ( !pdata.eof() )
        {
          linetnm1 = linet ;
          linet = sbuffer ;
        }
      }

      iss.str( linetnm1 );
      iss >> time_nm1 >> sbuffer >> sbuffer >> sbuffer >> vxnm1 >> vynm1 >>
	vznm1 >> omxnm1 >> omynm1 >> omznm1;
      iss.clear();

      iss.str( linet );
      iss >> time_n;
      iss.clear();

      previousdtfluid = time_n - time_nm1;

//      cout << previousdtfluid << " " << vznm1 << endl;

      (*particule)->setVelocityPreviousTimeRestart( vxnm1,
  	vynm1, vznm1, omxnm1, omynm1, omznm1 ) ;
      (*particule)->setVelocityAndVelocityDifferencePreviousTime() ;

      pdata.close();
    }
  }

  return ( previousdtfluid );
}



// ----------------------------------------------------------------------------
// Adjust the size of the voctor nbParticlesPerClass to the number of class
void EnsComposant::prepare_Polydisperse_vectors( int NbClass )
{
  m_nbParticlesPerClass.resize( NbClass, 0 );
  ParticleClassesConcentration.resize( NbClass, 0. );
  diameterRatio.resize( NbClass, 0. );
}




// ----------------------------------------------------------------------------
// set the number of particle per class
void EnsComposant::set_NbParticlesPerClass( int classe, int nb )
{
//  m_nbParticlesPerClass.push_back(nb) ;
  m_nbParticlesPerClass[classe] = nb ;
}




// ----------------------------------------------------------------------------
// Set the number of particle per class for restarts or coupled simulations.
// Used for dry granular simulation, when we know initialy the number of class

void EnsComposant::set_NbParticlesPerClass( void )
{
  list<Particule*>::const_iterator il;
  int i=0, classe=0;
  for (il=m_particulesActives.begin();
  	il!=m_particulesActives.end(); il++, i++)
  {
    if( (*il)->getTag() !=2 )
    {
      classe = (*il)->getParticuleClasse();
      m_nbParticlesPerClass[classe]+=1;
    }
  }
}




// ----------------------------------------------------------------------------
// Compute Sauter mean diameter
void EnsComposant::compute_SauterMeanDiameter( void )
{
  double diameter=0., num=0., den=0.;
  vector<Particule*>::iterator iv;
  int i=0;
  MPIWrapperGrains const* wrapper = Grains_Exec::getComm() ;

  for (iv=m_ParticuleClassesReference.begin();
  	iv!=m_ParticuleClassesReference.end(); iv++, i++)
  {
    diameter = 2*((*iv)->getRayon());
    num += m_nbParticlesPerClass[i] * pow(diameter,3) ;
    den += m_nbParticlesPerClass[i] * pow(diameter,2) ;
  }
  if( wrapper )
    SauterMeanDiameter = wrapper->sum_DOUBLE( num )/wrapper->sum_DOUBLE( den );
  else
    SauterMeanDiameter = num / den;
}




// ----------------------------------------------------------------------------
// Compute Sauter mean diameter
void EnsComposant::compute_ParticleClassesConcentration( void )
{
  // When we pass through this method from GRAINS-DRY, den and num are correct
  // before gathering. (because we use ParticleRef and we read the number of
  // particles directly in insert.xml). Thus we do Nproc*num / Nproc*den.
  // Whereas for grainsCoupledWithFluid(MPI)
  // we count them on each proc, then gather them.

  double diameter=0., num=0., den=0.;
  vector<Particule*>::iterator iv;
  int i=0;
  MPIWrapperGrains const* wrapper = Grains_Exec::getComm() ;

  for (iv=m_ParticuleClassesReference.begin();
  	iv!=m_ParticuleClassesReference.end(); iv++, i++)
  {
    diameter = 2*((*iv)->getRayon());
    den += m_nbParticlesPerClass[i] * pow(diameter,3) ;
  }

  i=0;
  for (iv=m_ParticuleClassesReference.begin();
  	iv!=m_ParticuleClassesReference.end(); iv++, i++)
  {
    diameter = 2*((*iv)->getRayon());

    num = m_nbParticlesPerClass[i] * pow(diameter,3) ;
    if( wrapper )
      ParticleClassesConcentration[i] = wrapper->sum_DOUBLE( num )
      	/ wrapper->sum_DOUBLE( den );
    else
      ParticleClassesConcentration[i] = num / den;

    diameterRatio[i] = diameter / SauterMeanDiameter;

    sumXY += ParticleClassesConcentration[i] * diameterRatio[i];
  }
}




// ----------------------------------------------------------------------------
// Return Sauter mean diameter
double EnsComposant::get_SauterMeanDiameter( )
{
  return( SauterMeanDiameter );
}




// ----------------------------------------------------------------------------
// Return Sauter mean diameter
double EnsComposant::get_sumXY( )
{
  return( sumXY );
}




// ----------------------------------------------------------------------------
// Mouvement aleatoire sur les particules actives
void EnsComposant::setRandomMotion( double const& coefTrans,
	double const& coefRot )
{
  for (list<Particule*>::iterator particule=m_particulesActives.begin();
	particule!=m_particulesActives.end();particule++)
    (*particule)->setRandomMotion( coefTrans, coefRot );
}




// ----------------------------------------------------------------------------
// Initialize all contact map entries to false in all particles
// and all elementary obstacles
void EnsComposant::setAllContactMapToFalse()
{
  for (list<Particule*>::iterator particule=m_particulesActives.begin();
	particule!=m_particulesActives.end();particule++)
    (*particule)->setContactMapToFalse();

  list<MonObstacle*> obstacles = m_obstacle->getObstacles();
  list<MonObstacle*>::iterator myObs;
  for( myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++ )
    (*myObs)->setContactMapToFalse();
}




// ----------------------------------------------------------------------------
// Update all contact map entries in all particles
// and all elementary obstacles
void EnsComposant::updateAllContactMaps()
{
  for (list<Particule*>::iterator particule=m_particulesActives.begin();
	particule!=m_particulesActives.end();particule++){
        (*particule)->updateContactMap();
    }


  list<MonObstacle*> obstacles = m_obstacle->getObstacles();
  list<MonObstacle*>::iterator myObs;
  for( myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++ )
    (*myObs)->updateContactMap();
}
