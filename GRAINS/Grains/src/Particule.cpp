#include "MPIWrapperGrains.hh"
#include "Particule.H"
#include "Cinematique_BuilderFactory.H"
#include "Contact_BuilderFactory.hh"
#include "GrainsCoupledWithFluid.hh"
#include "Grains_BuilderFactory.H"
#include "AppFluide_LubricationCorrection.H"
#include "Obstacle.H"
#include "MonObstacle.H"
#include "SaveTable.H"
#include "Memento.hh"
#include "Sphere.H"
#include "Disque.H"
#include "LinkedCell.H"
#include "ParticulePeriodique.hh"
#include "CompParticulePeriodique.hh"
#include "ObstaclePeriodique.hh"
#include "Grains_Exec.hh"
#include "ContactLaw.hh"
#include <algorithm>
#include <sstream>
#include <string>
// #include "AppFluide_Temperature.H"
using namespace std;


// Initialisation des attributs static
Scalar Particule::m_fluideMasseVolumique = 0.;
bool Particule::m_MassCorrection = false ;
bool Particule::m_explicitAddedMass = false ;
Scalar Particule::m_fluidViscosity = 0.;
double Particule::m_fluidInitialTemperature = 0.;


// ----------------------------------------------------------------------------
// Constructeur par defaut
Particule::Particule( const bool &autonumbering ) :
  Composant( autonumbering ),
  m_cinematique( NULL ),
  m_masseVolumique( 0.0 ),
  is_mobile( true ),
  m_temp_evolution( true ),
  m_energie( 0.0 ),
  m_activity( WAIT ),
  m_addedMassInfos( NULL ),
  m_tag( 0 ),
  m_GeoLoc( MPIGEO_NONE ),
  m_cellule_nm1( NULL ),
  m_ParticuleClasse( 0 ),
  m_coordination_number( 0 ),
  m_fluidInfos( NULL ),
  m_solidTemperature( 0. ),
  m_solidNusselt( 0. ),
  m_g_rnd(0., 0., 0.),
  m_g_rnd_Nu( 1., 1., 0. )
{
  // Masse & Inertie calculees
  for (int i=0; i<6; i++)
  {
    m_inertie[i] = 0.0;
    m_inertie_1[i] = 0.0;
  }

  // Type "P" (standard) ou "PP" (periodique)
  setType();
}




// ----------------------------------------------------------------------------
// Constructeur par copie
Particule::Particule( const Particule &copie ) :
  Composant( copie ),
  m_masseVolumique( copie.m_masseVolumique ),
  is_mobile(copie.is_mobile),
  m_temp_evolution(copie.m_temp_evolution),
  m_energie( copie.m_energie ),
  m_activity( WAIT ),
  m_addedMassInfos( NULL ),
  m_tag( copie.m_tag ),
  m_GeoLoc( copie.m_GeoLoc ),
  m_cellule_nm1( copie.m_cellule_nm1 ),
  m_ParticuleClasse( copie.m_ParticuleClasse ),
  m_coordination_number( copie.m_coordination_number ),
  m_weight( copie.m_weight ),
  m_fluidInfos( NULL ),
  m_solidTemperature( copie.m_solidTemperature ),
  m_solidNusselt( copie.m_solidNusselt ),
  m_g_rnd(copie.m_g_rnd),
  m_g_rnd_Nu(copie.m_g_rnd_Nu)
{
  m_cinematique = copie.m_cinematique->clone();
  copy( &copie.m_inertie[0], &copie.m_inertie[6], &m_inertie[0] );
  copy( &copie.m_inertie_1[0], &copie.m_inertie_1[6], &m_inertie_1[0] );

  setType( copie.m_type );

  if ( copie.m_addedMassInfos )
  {
    m_addedMassInfos = new struct AddedMassInfos;
    m_addedMassInfos->TranslationalVelocity_nm1 = copie.m_addedMassInfos->
    	TranslationalVelocity_nm1;
    m_addedMassInfos->TranslationalVelocity_difference =
    	copie.m_addedMassInfos->TranslationalVelocity_difference;
    m_addedMassInfos->RotationalVelocity_nm1 = copie.m_addedMassInfos->
    	RotationalVelocity_nm1;
    m_addedMassInfos->RotationalVelocity_difference = copie.m_addedMassInfos->
    	RotationalVelocity_difference;
    m_addedMassInfos->TranslationalSlipVelocity_nm1 = copie.m_addedMassInfos->
    	TranslationalSlipVelocity_nm1;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
Particule::Particule( DOMNode* root, const bool &autonumbering, const int &pc ):
  Composant( autonumbering ),
  m_cinematique( NULL ),
  m_masseVolumique( 2500. ),
  is_mobile( true ),
  m_temp_evolution( true ),
  m_energie( 0.0 ),
  m_activity( WAIT ),
  m_addedMassInfos( NULL ),
  m_tag( 0 ),
  m_GeoLoc( MPIGEO_NONE ),
  m_cellule_nm1( NULL ),
  m_ParticuleClasse( pc ),
  m_coordination_number( 0 ),
  m_fluidInfos( NULL ),
  m_solidTemperature( 0. ),
  m_solidNusselt( 0. ),
  m_g_rnd(0., 0., 0.),
  m_g_rnd_Nu( 1., 1., 0. )
{
  for (int i=0; i<6; i++)
  {
    m_inertie[i] = 0.0;
    m_inertie_1[i] = 0.0;
  }

  m_geoFormeVdw = new FormeVdW(root);
  m_cinematique = Cinematique_BuilderFactory::create(
  	m_geoFormeVdw->getConvex() );

  // Masse Volumique
  if ( ReaderXML::hasNodeAttr_Double( root, "MasseVolumique" ) )
  {
    m_masseVolumique = ReaderXML::getNodeAttr_Double( root, "MasseVolumique" );
  }

  // Materiau
  DOMNode* materiau_ = ReaderXML::getNode( root, "Materiau" );
  if ( materiau_ )
  {
    m_nomMateriau = ReaderXML::getNodeValue_String( materiau_ );
    Contact_BuilderFactory::defineMaterial( m_nomMateriau, false );
  }

  // Masse & Inertie calculees
  m_masse = m_masseVolumique * m_geoFormeVdw->getVolume();
  m_geoFormeVdw->BuildInertie( m_inertie, m_inertie_1 );
  for (int i=0; i<6; i++)
  {
    m_inertie[i] *= m_masseVolumique;
    m_inertie_1[i] /= m_masseVolumique;
  }

  //Does it move ?
  if ( ReaderXML::hasNodeAttr_String( root, "Mobilite" ) )
  {
    if (ReaderXML::getNodeAttr_String(root, "Mobilite" ) == "fix" )
      is_mobile = false;
    else if (ReaderXML::getNodeAttr_String(root, "Mobilite" ) == "mobile" )
      is_mobile = true;
    else
    {
      cout << "Mobility not recognized, particle is considered as mobile"<<endl;
      is_mobile = true;
    }
  }
  else is_mobile = true;

  // In case of temperature, is the temperature constant with time ?
  if ( ReaderXML::hasNodeAttr_String( root, "TemperatureEvolution" ) )
  {
    if (ReaderXML::getNodeAttr_String(root, "TemperatureEvolution" ) == "fixe" )
      m_temp_evolution = false;
    else if (ReaderXML::getNodeAttr_String(root, "TemperatureEvolution" ) == "evolutive" )
      m_temp_evolution = true;
    else
    {
      cout << "Temperature evolution not recognized, particle is considered as time evolving"<<endl;
      m_temp_evolution = true;
    }
  }
  else m_temp_evolution = true;

  // Type "P" (standard) ou "PP" (periodique)
  setType();

  // Calcul du poids
  computeWeight();

  // Masse ajoutee explicite
  if ( Particule::m_explicitAddedMass || Grains_Exec::m_addedmass_demcfd )
  	 createAddedMassInfos();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
Particule::Particule( const int &id_, Particule const* ParticuleRef,
	const double &vx, const double &vy, const double &vz,
	const double &qrotationx, const double &qrotationy,
	const double &qrotationz, const double &qrotations,
	const double &rx, const double &ry, const double &rz,
	const Scalar m[16],
	const ParticuleActivity &activ,
	const int &tag_,
	const int &coordination_number_ ) :
  Composant( false ),
  m_cinematique( NULL ),
  m_energie( 0.0 ),
  m_activity( activ ),
  m_addedMassInfos( NULL ),
  m_tag( tag_ ),
  m_GeoLoc( MPIGEO_NONE ),
  m_cellule_nm1( NULL ),
  m_coordination_number( coordination_number_ ),
  m_fluidInfos( NULL ),
  m_solidTemperature( 0. ),
  m_solidNusselt( 0. ),
  m_g_rnd(0., 0., 0.),
  m_g_rnd_Nu( 1., 1., 0. )
{
  // Initialisation de l'inertie
  for (int i=0; i<6; i++)
  {
    m_inertie[i] = 0.0;
    m_inertie_1[i] = 0.0;
 }

  // Numero de la particule
  m_id = id_;

  // Classe
  m_ParticuleClasse = ParticuleRef->m_ParticuleClasse;

  // Construction de la forme
  m_geoFormeVdw = new FormeVdW( *ParticuleRef->m_geoFormeVdw );

  // Transform
  m_geoFormeVdw->setPosition( m );

  // Cinematique
  m_cinematique = ParticuleRef->m_cinematique->clone();
  m_cinematique->setVitesseTranslation( Vecteur( vx, vy, vz ) );
  m_cinematique->setQuaternionRotation( qrotationx, qrotationy,
	qrotationz, qrotations );
  m_cinematique->setVitesseRotation( Vecteur( rx, ry, rz ) );

  // Materiau
  m_nomMateriau = ParticuleRef->m_nomMateriau;

  // Masse & Inertie calculees
  m_masseVolumique = ParticuleRef->m_masseVolumique;
  m_masse = ParticuleRef->m_masse;
  is_mobile = ParticuleRef->is_mobile;
  m_temp_evolution = ParticuleRef->m_temp_evolution;
  for (int i=0; i<6; i++)
  {
    m_inertie[i] = ParticuleRef->m_inertie[i];
    m_inertie_1[i] = ParticuleRef->m_inertie_1[i];
  }

  // Type "P" (standard) ou "PP" (periodique)
  setType();

  // Calcul du poids
  computeWeight();

  // Masse ajoutee explicite
  if ( Particule::m_explicitAddedMass || Grains_Exec::m_addedmass_demcfd )
  	 createAddedMassInfos();

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
Particule::Particule( const int &id_, Particule const* ParticuleRef,
	const Vecteur &vtrans,
	const Quaternion &qrot,
	const Vecteur &vrot,
	const Transform &config,
	const ParticuleActivity &activ,
	const int &tag_ ) :
  Composant( false ),
  m_cinematique( NULL ),
  m_energie( 0.0 ),
  m_activity( activ ),
  m_addedMassInfos( NULL ),
  m_tag( tag_ ),
  m_GeoLoc( MPIGEO_NONE ),
  m_cellule_nm1( NULL ),
  m_coordination_number( 0 ),
  m_fluidInfos( NULL ),
  m_solidTemperature( 0. ),
  m_solidNusselt( 0. ),
  m_g_rnd(0., 0., 0.),
  m_g_rnd_Nu( 1., 1., 0. )
{
  // Initialisation de l'inertie
  for (int i=0; i<6; i++)
  {
    m_inertie[i] = 0.0;
    m_inertie_1[i] = 0.0;
  }

  // Numero de la particule
  m_id = id_;

  // Classe
  m_ParticuleClasse = ParticuleRef->m_ParticuleClasse;

  // Construction de la forme
  m_geoFormeVdw = new FormeVdW( *ParticuleRef->m_geoFormeVdw );

  // Transform
  m_geoFormeVdw->setTransform( config );

  // Cinematique
  m_cinematique = ParticuleRef->m_cinematique->clone();
  m_cinematique->setVitesseTranslation( vtrans );
  m_cinematique->setQuaternionRotation( qrot );
  m_cinematique->setVitesseRotation( vrot );

  // Materiau
  m_nomMateriau = ParticuleRef->m_nomMateriau;

  // Masse & Inertie calculees
  m_masseVolumique = ParticuleRef->m_masseVolumique;
  m_masse = ParticuleRef->m_masse;
  is_mobile = ParticuleRef->is_mobile;
  m_temp_evolution = ParticuleRef->m_temp_evolution;
  for (int i=0; i<6; i++)
  {
    m_inertie[i] = ParticuleRef->m_inertie[i];
    m_inertie_1[i] = ParticuleRef->m_inertie_1[i];
  }

  // Type "P" (standard) ou "PP" (periodique)
  setType();

  // Calcul du poids
  computeWeight();

  // Masse ajoutee explicite
  if ( Particule::m_explicitAddedMass || Grains_Exec::m_addedmass_demcfd )
  	 createAddedMassInfos();
}




// ----------------------------------------------------------------------------
// Constructeur d'un clone par copie
Particule* Particule::createCloneCopy() const
{
  Particule* particule = new Particule( *this );
  return particule;
}




// ----------------------------------------------------------------------------
// Calcul du poids de la particule
void Particule::computeWeight()
{
  m_weight = m_masse
  	* ( 1. - Particule::m_fluideMasseVolumique / m_masseVolumique )
  	* Grains_Exec::m_vgravite ;

}




// ----------------------------------------------------------------------------
// Destructeur
Particule::~Particule()
{
  delete m_cinematique;

  // Ici, la liste de pointeurs sur les particules periodiques
  // est simplement videe; la destruction des objets pointes est realisee
  // par le destructeur de EnsComposant car les particules periodiques
  // font parties de la liste "particulesClonesPeriodiques" de EnsComposant
  m_periodicClones.clear();
  m_periodicObstaclesID.clear();

  if ( m_addedMassInfos ) delete m_addedMassInfos;

  if ( m_fluidInfos ) delete m_fluidInfos;
}




// ----------------------------------------------------------------------------
// Resolution des equations de la dynamique et double integration pour
// obtenir la nouvelle vitesse et position.
void Particule::Deplacer( Scalar temps, double dt ) throw(ErreurDeplacement)
{
  // Double integration du principe fondamental de la dynamique.
  m_cinematique->CalculerVariationQDM( m_somme, dt, this );
  double depl = m_cinematique->Deplacer( this, dt );

  // On regarde si le deplacement est superieur au rayon d'interaction
  double rayon = m_geoFormeVdw->getRayonInterAction();
  if ( depl > rayon )
  {
    cout << endl << "Processor = " <<
    	(Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
	<< " has thrown an ErreurDeplacement exception" <<  endl;
    Grains_Exec::m_exception_Deplacement = true;
    ErreurDeplacement erreur( this, depl, rayon, temps );
    throw erreur;
  }
}




// ----------------------------------------------------------------------------
// Solve temperature problem
void Particule::ComputeTemperature( Scalar temps, double dt )
{
  m_solidTemperature = m_solidTemperature
    + dt/(m_masse*AppFluide_Temperature::m_heatCapacityS) * m_sum_HeatFlux;
  // cout << " TEMPORARY : m_solidTemperature "
  //      << dt/(m_masse*AppFluide_Temperature::m_heatCapacityS)* m_sum_HeatFlux << endl;

}




// ----------------------------------------------------------------------------
// Contact entre la particule et une particule.
void Particule::InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC )
  throw ( ErreurContact )
{
  PointContact closestPoint;
  double delta=0.;
  //double R=0., A=0.;

  try {
    closestPoint = m_geoFormeVdw->ClosestPoint( *(voisin->getForme()) );
  }
  catch ( ErreurContact &erreur )
  {
    try {
      closestPoint = m_geoFormeVdw->ClosestPoint_ErreurHandling(
      	*(voisin->getForme()), 10., m_id, voisin->getID() );
    }
    catch (ErreurContact &erreur_level2)
    {
      cout << endl << "Processor = "
           << ( Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
	<< " has thrown an ErreurContact exception" <<  endl;
      erreur_level2.setMessage( "Particule::InterAction : choc de croute a t="
      	+ Grains_Exec::doubleToString( temps, TIMEFORMAT ) );
      erreur_level2.setComposants( this, voisin, temps );
      Grains_Exec::m_exception_Contact = true;
      throw(erreur_level2);
    }
  }

  LC->addToContactsFeatures( temps, closestPoint );
  delta = closestPoint.getDistance();
  if( delta < 0. )
  {
    if ( Contact_BuilderFactory::contactForceModel(
		m_nomMateriau, voisin->materiau() )
    		->computeForces( this, voisin, closestPoint, LC, dt ) )
    {
      this->addToCoordinationNumber( 1 );
      voisin->addToCoordinationNumber( 1 );
    }
  }
  else if ( Grains_Exec::m_withlubrication )
      (GrainsCoupledWithFluid::LubricationCorrection())->computeforce( this,
		     voisin, LC, dt );

  if( Grains_Exec::m_withSolidTemperature )
  {
    cout << "TEMP : Call 'app Solid Temperature compute flux' here" << endl;
    //R = getRayonSphereEquivalente();
    //A = PI * pow(R,2.) * pow(sin(acos((R-delta)/R)),2.);
//    app_Temperature->computeSolidHeatTransfer( this, voisin, closestPoint, LC, dt );
  }
}




// ----------------------------------------------------------------------------
// Contact entre la particule et une particule.
void Particule::InterActionCohesiveInit( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC )
  throw ( ErreurContact )
{
  PointContact closestPoint;
  try {
    closestPoint = m_geoFormeVdw->ClosestPoint( *(voisin->getForme()) );
  }
  catch ( ErreurContact &erreur )
  {
    try {
      closestPoint = m_geoFormeVdw->ClosestPoint_ErreurHandling(
          *(voisin->getForme()), 10., m_id, voisin->getID() );
    }
    catch( ErreurContact &erreur_level2 )
    {
      cout << endl << "Processor = "
           << ( Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
           << " has thrown an ErreurContact exception" <<  endl;
      erreur_level2.setMessage( "Particule::InterAction : choc de croute a t="
          + Grains_Exec::doubleToString( temps, TIMEFORMAT ) );
      erreur_level2.setComposants( this, voisin, temps );
      Grains_Exec::m_exception_Contact = true;
      throw(erreur_level2);
    }
  }

  LC->addToContactsFeatures( temps, closestPoint );

  if( closestPoint.getDistance() < 0. )
  {
    if( Contact_BuilderFactory::contactForceModel(
        m_nomMateriau, voisin->materiau() )
        ->InitializeCohesiveObjects( this, voisin, closestPoint, LC, dt ) )
    {
      this->addToCoordinationNumber( 1 );
      voisin->addToCoordinationNumber( 1 );
    }
  }
}




// ----------------------------------------------------------------------------
// Contact entre particule composite et une particule et/ou un obstacle
// D. RAKOTONIRINA - Dec. 2014 - Creation
void Particule::SearchContact( Composant* voisin, double dt,
      double const& temps, LinkedCell *LC,
      list<ContactInfos*> &listContact )
{
  PointContact closestPoint;
  try {
    closestPoint = m_geoFormeVdw->ClosestPoint( *(voisin->getForme()) );
  }
  catch ( ErreurContact &erreur )
  {
    try {
      closestPoint = m_geoFormeVdw->ClosestPoint_ErreurHandling(
      	*( voisin->getForme() ), 10., m_id, voisin->getID() );
    }
    catch ( ErreurContact &erreur_level2 )
    {
      cout << endl << "Processor = " <<
    	( Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
	<< " has thrown an ErreurContact exception" <<  endl;
      erreur_level2.setMessage( "Particule::InterAction : choc de croute a t="
      	+Grains_Exec::doubleToString( temps, TIMEFORMAT ) );
      erreur_level2.setComposants( this, voisin, temps );
      Grains_Exec::m_exception_Contact = true;
      throw( erreur_level2 );
    }
  }

  ContactInfos* result = NULL ;
  if ( closestPoint.getDistance() < 0. )
  {
    result = new struct ContactInfos;
    result->ContactPoint = closestPoint;
    result->p0 = this;
    result->p1 = voisin;
    listContact.push_back( result );
  }

//  return ( result );
}




// ----------------------------------------------------------------------------
// Contact entre la particule et une particule pour le post processing
void Particule::InterActionPostProcessing( Composant* voisin, Scalar dt,
	list<struct PointForcePostProcessing>* listOfContacts )
  throw ( ErreurContact )
{
  PointContact closestPoint;
  try {
    closestPoint = m_geoFormeVdw->ClosestPoint( *(voisin->getForme()) );
  }
  catch ( ErreurContact &erreur )
  {
    closestPoint = m_geoFormeVdw->ClosestPoint_ErreurHandling(
      	*( voisin->getForme() ), 10., m_id, voisin->getID() );
  }

  if ( closestPoint.getDistance() < 0. )
    Contact_BuilderFactory::contactForceModel(
    	m_nomMateriau, voisin->materiau() )
    	->computeForcesPostProcessing( this, voisin, dt, closestPoint,
	listOfContacts );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Acces a la cinematique de la particule
const CineParticule* Particule::getCinematique() const
{
  return m_cinematique;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Masse volumique du fluide en interaction avec les particules
double Particule::getFluideMasseVolumique()
{
  return Particule::m_fluideMasseVolumique;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Viscosite du fluide en interaction avec les particules
double Particule::getFluidViscosity()
{
  return Particule::m_fluidViscosity ;
}




// ----------------------------------------------------------------------------
// Inertie de la particule.
const double* Particule::getInertie() const
{
  return m_inertie;
}




// ----------------------------------------------------------------------------
// Inertie inverse de la particule.
const double* Particule::getInertieInverse() const
{
  return m_inertie_1;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Masse volumique de la particule
Scalar Particule::getMasseVolumique() const
{
  return m_masseVolumique;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Energie de la particule
// D. RAKOTONIRINA - Oct. 2014 - Creation
double Particule::getEnergie() const
{
  return m_energie;
}




// ----------------------------------------------------------------------------
// Rayon de la sphere de meme volume que le composant
// A. WACHS- Sept.2014 - Creation
Scalar Particule::getRayonSphereEquivalente() const
{
  return ( pow( ( 0.75 / PI ) * m_masse / m_masseVolumique, 1. / 3. ) );
}




// ----------------------------------------------------------------------------
// Vitesse en un point de la particule
// G.FERRER - Janv.2002 - Creation
// A.WACHS - Octo.2010 - Modif
Vecteur Particule::getVitesse( const Point &pt ) const
{
  Vecteur levier = pt - *m_geoFormeVdw->getCentre();
  return m_cinematique->Vitesse( levier );
}




// ----------------------------------------------------------------------------
// Vitesse de rotation
// G.FERRER - Janv.2002 - Creation
Vecteur const* Particule::getVitesseRotation() const
{
  return m_cinematique->getVitesseRotation();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Vitesse de Translation
// G.FERRER - Aout.2004 - Creation
Vecteur const* Particule::getVitesseTranslation() const
{
  return m_cinematique->getVitesseTranslation();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Force torseur
// A.HAMMOUTI - Aout.2013 - Creation
Vecteur const* Particule::getForce() const
{
  return m_somme.getForce();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Force torseur
// A. Esteghamatian - Aout.2015 - Creation
Vecteur const* Particule::getForceContactPP() const
{
  return &m_ForceContactPP;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Force torseur
// A. Esteghamatian - Aout.2015 - Creation
Vecteur const* Particule::getForceContactPP_instantaneous() const
{
  return &m_ForceContactPP_instantaneous;
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Force torseur
// A. Esteghamatian - Aout.2015 - Creation
Vecteur const* Particule::getForceLubriPP() const
{
  return &m_ForceLubriPP;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Quaternion
// A.HAMMOUTI - Aout.2013 - Creation
Quaternion const* Particule::getRotation() const
{
  return m_cinematique->getRotation();
}




// ----------------------------------------------------------------------------
// Get body temperature
double const* Particule::get_solidTemperature() const
{
  return &m_solidTemperature;
}




// ----------------------------------------------------------------------------
// Quaternion de rotation dans le vecteur vit en debutant a la position i
void Particule::copyQuaternionRotation( double *vit, int i ) const
{
  Quaternion const* qr = m_cinematique->getRotation();
  Vecteur const* vqr = qr->getVecteur();
  for (int j=0 ;j<3; j++) vit[i+j] = (*vqr)[j];
  vit[i+3] = qr->getScalaire();
}




// ----------------------------------------------------------------------------
// Copie force & moment exerces sur la particule dans le vecteur fm
// en debutant a la position i
void Particule::copyForceMoment( double *fm, int i ) const
{
  m_somme.copyForceMoment( fm, i );
}




// ----------------------------------------------------------------------------
// Vitesse de rotation dans le vecteur vit en debutant a la position i
void Particule::copyVitesseRotation( double *vit, int i ) const
{
  Vecteur const* vr = m_cinematique->getVitesseRotation();
  for (int j=0 ;j<3; j++) vit[i+j] = (*vr)[j];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Vitesse de translation dans le vecteur vit en debutant a la position i
void Particule::copyVitesseTranslation( double *vit, int i ) const
{
  Vecteur const* vt = m_cinematique->getVitesseTranslation();
  for (int j=0 ;j<3; j++) vit[i+j] = (*vt)[j];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copie de la cinematique au temps t-2dt: vitesse de translation,
// vitese de rotation, variation de QDM translationnelle, variation de QDM
// rotationnalle
void Particule::copyCinematiqueNm2( double *vit, int i ) const
{
  m_cinematique->copyCinematiqueNm2( vit, i );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Vitesse de translation + Gradient de pression + Fraction volumique du
// fluide dans le vecteur vit en debutant a la position i
void Particule::copyFluidInformations( double *vit, int i ) const
{
//  bool b_isLiftForce = Grains_Exec::m_withLiftForce ;

  vit[i]   = m_fluidInfos->m_demcfd_epsilon;
  for (int j=0; j<3; j++)
  {
    vit[i+1+j]   = m_fluidInfos->m_vitesseTr_fluide[j];
//    if( b_isLiftForce )
//    {
//      vit[i+4+j] = m_fluidInfos->m_vorticity_fluide[j];
//      vit[i+7+j] = m_fluidInfos->m_gradientPression_fluide[j];
//    }
//    else
    vit[i+4+j] = m_fluidInfos->m_gradientPression_fluide[j];
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy Fluid Vorticity interpolated at particle center into buffer
// Used in case of lift force in DEMCFD
void Particule::copyFluidVorticity( double *vit, int i ) const
{
  for( int j=0; j<3; j++ )
    vit[i+j] = m_fluidInfos->m_vorticity_fluide[j];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy Fluid Temperature into buffer
void Particule::copy_fluidTemperature( double *vit, int i ) const
{
  vit[i] = m_fluidInfos->m_demcfd_fluidTemperature;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy Solid Temperature into buffer
void Particule::copy_solidTemperature( double *vit, int i ) const
{
  vit[i] = m_solidTemperature;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy Solid Nusselt into buffer
void Particule::copy_solidNusselt( double *vit, int i ) const
{
  vit[i] = m_solidNusselt;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy rnd number for stochastic drag into buffer
void Particule::copy_rnd( double *vit, int i ) const
{
  for( int j=0; j<3; j++ )
    vit[i+j] = m_g_rnd[j];
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy rnd number for stochastic Nusselt into buffer
void Particule::copy_rnd_Nu( double *vit, int i ) const
{
  for( int j=0; j<3; j++ )
    vit[i+j] = m_g_rnd_Nu[j];
}

// ----------------------------------------------------------------------------
// Ajout d'une force au torseur des efforts
void Particule::addForce( const Point &point, const Vecteur &force )
{
  m_somme.addForce( point, force );
}




// ----------------------------------------------------------------------------
// Force at the contact point for the stress tensor sans fenetre
// D. RAKOTONIRINA - Fev. 2017 - Creation
void Particule::computeInternalMoments( const Point &point, const Vecteur &force )
{
  if ( getTag() != 2 )
  {
    Point gc = *getPosition();
    Scalar px = point[X] - gc[X];
    Scalar py = point[Y] - gc[Y];
    Scalar pz = point[Z] - gc[Z];
    Scalar fx = force[X];
    Scalar fy = force[Y];
    Scalar fz = force[Z];

    // Compute the internal moments of the particle
    m_mxx += fx*px;  m_mxy += fx*py; m_mxz += fx*pz;
    m_myx += fy*px;  m_myy += fy*py; m_myz += fy*pz;
    m_mzx += fz*px;  m_mzy += fz*py; m_mzz += fz*pz;

  }
}




// ----------------------------------------------------------------------------
// Get the stress tensor of the particle divided its volume
// J. F. Peters et al. Characterization of force chains in granular material.
// Phys. Rev. 2005
// D. RAKOTONIRINA - Fev. 2017 - Creation
vector<Scalar> const* Particule::getStressTensor()
{
  vector<Scalar>* stressTensor = new vector<Scalar>(9,0.);
  Scalar volume = m_geoFormeVdw->getVolume();
  (*stressTensor)[0] = m_mxx / volume;
  (*stressTensor)[1] = m_mxy / volume;
  (*stressTensor)[2] = m_mxz / volume;
  (*stressTensor)[3] = m_myx / volume;
  (*stressTensor)[4] = m_myy / volume;
  (*stressTensor)[5] = m_myz / volume;
  (*stressTensor)[6] = m_mzx / volume;
  (*stressTensor)[7] = m_mzy / volume;
  (*stressTensor)[8] = m_mzz / volume;

  return stressTensor;
}




// ----------------------------------------------------------------------------
// Get the internal moments of the particle divided by the volume of the sample
// Note: the volume varies according to the applied stress
// D. RAKOTONIRINA - Fev. 2017 - Creation
vector<Scalar> const* Particule::getInternalMoment()
{
  m_InternalMoment[0] = m_mxx;
  m_InternalMoment[1] = m_mxy;
  m_InternalMoment[2] = m_mxz;
  m_InternalMoment[3] = m_myx;
  m_InternalMoment[4] = m_myy;
  m_InternalMoment[5] = m_myz;
  m_InternalMoment[6] = m_mzx;
  m_InternalMoment[7] = m_mzy;
  m_InternalMoment[8] = m_mzz;

  return &m_InternalMoment;
}




// ----------------------------------------------------------------------------
// Initialise force at the contact point for the post-processing
// of the stress tensor
// D. RAKOTONIRINA - Fev. 2017 - Creation
void Particule::InitializeForceAtContactPoint()
{
  m_InternalMoment = vector<Scalar>(9, 0.);
  m_mxx = 0., m_mxy = 0., m_mxz = 0.;
  m_myx = 0., m_myy = 0., m_myz = 0.;
  m_mzx = 0., m_mzy = 0., m_mzz = 0.;
}




// ----------------------------------------------------------------------------
// Add contact force on each particle for postprocessing purposes
void Particule::addContactForcePP( const Vecteur &force )
{
  Vecteur absforce;
  absforce[X] = fabs(force[X]);
  absforce[Y] = fabs(force[Y]);
  absforce[Z] = fabs(force[Z]);
  m_ForceContactPP += absforce;
}




// ----------------------------------------------------------------------------
// Add contact force on each particle for postprocessing purposes
void Particule::addContactForcePP_instantaneous( const Vecteur &force )
{
  m_ForceContactPP_instantaneous += force;
}




// ----------------------------------------------------------------------------
// Add contact force on each particle for postprocessing purposes
void Particule::addLubriForcePP( const Vecteur &force )
{
  Vecteur absforce;
  absforce[X] = fabs(force[X]);
  absforce[Y] = fabs(force[Y]);
  absforce[Z] = fabs(force[Z]);
  m_ForceLubriPP += absforce;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajoute une force et un moment au torseur des efforts exerces sur la
// particule (utile en simulation periodique)
void Particule::addForceMoment( const double &fx, const double &fy,
	const double &fz, const double &mx, const double &my, const double &mz )
{
  m_somme.addForce( fx, fy, fz );
  m_somme.addMoment( mx, my, mz );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Remise a zero de la particule
// Usage : suppression de la particule pour reinsertion
void Particule::reset()
{
  Point zero;
  m_geoFormeVdw->setOrigin( (Scalar *) &zero );
  m_cinematique->reset();
}




// ----------------------------------------------------------------------------
// Remise a zero arbitraire de la cinematique
void Particule::ResetCinematique()
{
  m_cinematique->reset();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde de l'etat
void Particule::saveState()
{
  if (!m_memento)
    m_memento = new ConfigurationMemento();
  m_memento->m_position = *m_geoFormeVdw->getTransform();
  m_cinematique->saveState();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cree et renvoie l'etat
pair<ConfigurationMemento*,CineParticuleMemento*> Particule::createState()
{
  ConfigurationMemento* Pmemento_ = new ConfigurationMemento();
  Pmemento_->m_position = *m_geoFormeVdw->getTransform();
  CineParticuleMemento* Cmemento_ = m_cinematique->createState();
  pair<ConfigurationMemento*,CineParticuleMemento*> ppp( Pmemento_, Cmemento_ );

  return ppp;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void Particule::restaureState()
{
  m_geoFormeVdw->setTransform( m_memento->m_position );
  m_cinematique->restaureState();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void Particule::restaureState( ConfigurationMemento const* Pmemento_,
  	CineParticuleMemento const* Cmemento_,
	ObstaclePeriodique const* m_contactReferenceObstacle,
	Particule* reference_ )
{
  m_geoFormeVdw->setTransform( Pmemento_->m_position );
  m_cinematique->restaureState( Cmemento_ );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation de la masse volumique du fluide
void Particule::setFluideMasseVolumique( double rho )
{
  m_fluideMasseVolumique = rho;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Viscosite du fluide en interaction avec les particules
void Particule::setFluidViscosity(double mu)
{
  m_fluidViscosity = mu;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Temperature du fluide en interaction avec les particules
void Particule::setFluidInitialTemperature( double tempF )
{
  m_fluidInitialTemperature = tempF;
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// GET vecteur vitesse de translation du fluide interpole
// au centre de gravite de la particule
// M. Bernard - Janvier 2012 - Creation
Vecteur const* Particule::getVitesseTranslation_fluide() const
{
  return &(m_fluidInfos->m_vitesseTr_fluide);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// GET vecteur vorticite du fluide interpole
// au centre de gravite de la particule
// M. Bernard - Janvier 2014 - Creation
Vecteur const* Particule::getVorticity_fluide() const
{
  return &(m_fluidInfos->m_vorticity_fluide);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Gradient de pression du fluide
// M. Bernard - Janvier 2012 - Creation
Vecteur const* Particule::getGradientPression_fluide() const
{
  return &(m_fluidInfos->m_gradientPression_fluide);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// GET Fraction volumique
// M. Bernard - Janvier 2012 - Creation
double Particule::get_DEMCFD_volumeFraction() const
{
  return m_fluidInfos->m_demcfd_epsilon;
}




// ----------------------------------------------------------------------------
// Get DEMCFD surounding fluid temperature
double const* Particule::get_DEMCFD_fluidTemperature() const
{
  return &m_fluidInfos->m_demcfd_fluidTemperature;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// GET Force hydrodynamique ==> to apply backward force
// M. Bernard - Fevrier 2012 - Creation
void Particule::getParticleHydroForce(double& Fx, double& Fy, double& Fz)
{
  Fx = m_fluidInfos->m_hydroForce[X];
  Fy = m_fluidInfos->m_hydroForce[Y];
  Fz = m_fluidInfos->m_hydroForce[Z];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// GET Force hydrodynamique ==> for post-processing purposes
// M. Bernard - Septembre 2014 - Creation
Vecteur const* Particule::getParticleHydroForce() const
{
  return &(m_hydroForce_PP);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// GET Slip Velocity
// A. ESTEGHAMATIAN - March 2015 - Creation
Vecteur const* Particule::getParticleSlipVel() const
{
  return &(m_fluidInfos->m_slipvel);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation du vecteur vitesse de translation du fluide interpole
// au centre de gravite de la particule
// M. Bernard - Janvier 2012 - Creation
void Particule::setVitesseTr_fluide(
    double const& vitesseTr_fluide_x,
    double const& vitesseTr_fluide_y,
    double const& vitesseTr_fluide_z )
{
  m_fluidInfos->m_vitesseTr_fluide[X] = vitesseTr_fluide_x;
  m_fluidInfos->m_vitesseTr_fluide[Y] = vitesseTr_fluide_y;
  m_fluidInfos->m_vitesseTr_fluide[Z] = vitesseTr_fluide_z;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation de la vorticite du fluide interpolee
// au centre de gravite de la particule
// M. Bernard - Janvier 2014 - Creation
void Particule::setVorticity_fluide(
    double const& vorticity_fluide_x,
    double const& vorticity_fluide_y,
    double const& vorticity_fluide_z )
{
  m_fluidInfos->m_vorticity_fluide[X] = vorticity_fluide_x;
  m_fluidInfos->m_vorticity_fluide[Y] = vorticity_fluide_y;
  m_fluidInfos->m_vorticity_fluide[Z] = vorticity_fluide_z;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation du gradient de pression du fluide interpole sur la particule
// M. Bernard - Janvier 2012 - Creation
void Particule::setGradientPression_fluide(
	double const& gradP_fluid_x,
  	double const& gradP_fluid_y,
	double const& gradP_fluid_z )
{
  m_fluidInfos->m_gradientPression_fluide[X] = gradP_fluid_x;
  m_fluidInfos->m_gradientPression_fluide[Y] = gradP_fluid_y;
  m_fluidInfos->m_gradientPression_fluide[Z] = gradP_fluid_z;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation de la fraction volumique interpolee sur la particule
// M. Bernard - Fevrier 2012 - Creation
void Particule::set_DEMCFD_volumeFraction( double const& fluidVolumeFraction_ )
{
  m_fluidInfos->m_demcfd_epsilon = fluidVolumeFraction_ ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation de la force hydrodynamique
// M. Bernard - Fevrier 2012 - Creation
void Particule::setHydroForce( Vecteur const* hydroForce_ )
{
  m_fluidInfos->m_hydroForce = *hydroForce_ ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Keep hydroforce for post-processing
// A. Esteghamatian - Creation
void Particule::setHydroForce_PP( Vecteur const* hydroForce_ )
{
  m_hydroForce_PP = *hydroForce_ ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set slip velocity
// A. Esteghamatian - March2015 - Creation
void Particule::setSlipVel( Vecteur const* slipVel_ )
{
  m_fluidInfos->m_slipvel = *slipVel_ ;
}




// ----------------------------------------------------------------------------
// set DEMCFD surounding fluid temperature
void Particule::set_DEMCFD_fluidTemperature( const Scalar fluidTemperature_ )
{
  m_fluidInfos->m_demcfd_fluidTemperature = fluidTemperature_ ;
}




// ----------------------------------------------------------------------------
// set body temperature
void Particule::set_solidTemperature( const Scalar solidTemperature_ )
{
  m_solidTemperature = solidTemperature_ ;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set the random number vector for fluctuating drag
// A. ESTEGHAMATIAN 2016
void Particule::set_rnd(
	double const& g_rnd_x,
  	double const& g_rnd_y,
	double const& g_rnd_z )
{
  m_g_rnd[0] = g_rnd_x;
  m_g_rnd[1] = g_rnd_y;
  m_g_rnd[2] = g_rnd_z;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set the random number vector for fluctuating drag
// A. ESTEGHAMATIAN 2016
Vecteur const * Particule::get_rnd() const
{
  return &m_g_rnd;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set the random number vector for fluctuating drag
// F. EUZENAT 2017
void Particule::set_rnd_Nu(
	double const& g_rnd_Nux,
        double const& g_rnd_Nuy,
        double const& g_rnd_Nuz )
{
  m_g_rnd_Nu[0] = g_rnd_Nux;
  m_g_rnd_Nu[1] = g_rnd_Nuy;
  m_g_rnd_Nu[2] = g_rnd_Nuz;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set the random number vector for fluctuating drag
// F. EUZENAT 2017
Vecteur const * Particule::get_rnd_Nu() const
{
  return &m_g_rnd_Nu;
}
// ----------------------------------------------------------------------------
// Set Fluid-Solid heat flux (to compute backward flux on fluid)
void Particule::set_fluidSolidHeatFlux( const double heatFlux_ )
{
  m_fluidInfos->fluidSolid_heatFlux = heatFlux_;
}



// ----------------------------------------------------------------------------
// Get Fluid-Solid heat flux (to compute backward flux on fluid)
double const* Particule::get_fluidSolidHeatFlux( void ) const
{
  return &m_fluidInfos->fluidSolid_heatFlux;
}


// ----------------------------------------------------------------------------
// Set Solid Nusselt number
void Particule::set_solidNusselt( const Scalar Nusselt_ )
{
  m_solidNusselt = Nusselt_;
}




// ----------------------------------------------------------------------------
// Get Solid Nusselt number
double const* Particule::get_solidNusselt( void ) const
{
  return &m_solidNusselt;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation d'un champ de vitesse de rotation
// G.FERRER - Aout.2004 - Creation
void Particule::setVitesseRotation( const Vecteur &rotation )
{
  m_cinematique->setVitesseRotation( rotation );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation d'un champ de vitesse de translation
// G.FERRER - Aout.2004 - Creation
void Particule::setVitesseTranslation( const Vecteur &translation )
{
  m_cinematique->setVitesseTranslation( translation );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cinematique au temps t-2dt: vitesse de translation,
// vitese de rotation, variation de QDM translationnelle, variation de QDM
// rotationnalle
void Particule::setCinematiqueNm2( double const* tab )
{
  m_cinematique->setCinematiqueNm2( tab );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture de la particule pour Reload. Utilise pour les particules de
// classe de reference dans l'entete du fichier de reload et pour l'ancien
// format de reload
// G.FERRER - Janv.2004 - Creation
void Particule::read( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference )
{
  string buffer, adresse;

  fileSave >> buffer
	   >> adresse >> m_id;
  SaveTable::create( adresse, this );

  fileSave >> buffer
	   >> m_nomMateriau;

  // Si c'est une particule de reference de classe, on utilise le constructeur
  // standard de la forme afin de creer tous les attributs de ce type
  // FormeVdW( fileSave ) => Convex_BuilderFactory::create( cle, fileSave )
  // => xxx::create( fileSave ) => constructeur de xxx
  // Ceci est du au cas particulier des polytopes
  if ( !ParticuleClassesReference ) m_geoFormeVdw = new FormeVdW( fileSave );
  else
  {
    string cle;
    fileSave >> cle;
    while ( cle != "*END" ) fileSave >> cle;
  }

  fileSave >> buffer >> m_ParticuleClasse;

  // Si c'est une particule active, on utilise le constructeur de recopie
  if ( ParticuleClassesReference )
    m_geoFormeVdw = new FormeVdW(
    *(*ParticuleClassesReference)[m_ParticuleClasse]->getForme() );

  fileSave >> buffer
	   >> m_tag;
  fileSave >> buffer
	   >> m_masse >> m_energie;
   // We check if the mobility is precised in the text reload file then we read
   // it if it exist, otherwhise we put back the input stream to the initial
   // position ( to be conformant with both old and new reload formats (2016))
  streampos length = fileSave.tellg();
  string is_mobilite_precised;
  fileSave >> is_mobilite_precised;
  if ( is_mobilite_precised=="*Mobilite" )
    fileSave >> is_mobile;
  else
   fileSave.seekg( length );

   // We check if the mobility is precised in the text reload file then we read
   // it if it exist, otherwhise we put back the input stream to the initial
   // position ( to be conformant with both old and new reload formats (2016))
  length = fileSave.tellg();
  string is_temperature_evolution_precised;
  fileSave >> is_temperature_evolution_precised;
  if ( is_temperature_evolution_precised=="*TempEvolution" )
    fileSave >> m_temp_evolution;
  else
   fileSave.seekg( length );

  fileSave >> buffer
	   >> m_inertie[0] >> m_inertie[1] >> m_inertie[2]
	   >> m_inertie[3] >> m_inertie[4] >> m_inertie[5];
  fileSave >> buffer
	   >> m_inertie_1[0] >> m_inertie_1[1] >> m_inertie_1[2]
	   >> m_inertie_1[3] >> m_inertie_1[4] >> m_inertie_1[5];

  int buf = 0;
  fileSave >> buffer >> buf >> buf >> buf;
  m_geoFormeVdw->readPosition( fileSave );

  bool actif;
  fileSave >> buffer >> actif;
  m_activity = ( actif == true ) ? COMPUTE : WAIT;

  if ( m_cinematique ) delete m_cinematique;
  m_cinematique = Cinematique_BuilderFactory::read( fileSave,
  	m_geoFormeVdw->getConvex() );

  m_masseVolumique = m_masse / m_geoFormeVdw->getVolume();
  computeWeight();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture de la particule pour Reload. Utilise pour les particules
// dans la simulation (actives ou en attentes) pour le format de reload 2014
// A.WACHS - Aout 2014 - Creation
void Particule::read2014( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference )
{
  // Lecture du n\B0 et de la classe de reference
  fileSave >> m_id >> m_ParticuleClasse;

  // Creation de la forme (convex + transform) en utilisant le constructeur
  // de recopie
  m_geoFormeVdw = new FormeVdW(
    *(*ParticuleClassesReference)[m_ParticuleClasse]->getForme() );

  // Materiau, masse, energie et inertie de la particule de class de reference
  m_nomMateriau =
  	(*ParticuleClassesReference)[m_ParticuleClasse]->m_nomMateriau ;
  m_masse = (*ParticuleClassesReference)[m_ParticuleClasse]->m_masse ;
  m_energie = (*ParticuleClassesReference)[m_ParticuleClasse]->m_energie ;
  is_mobile = (*ParticuleClassesReference)[m_ParticuleClasse]->is_mobile ;
  m_temp_evolution = (*ParticuleClassesReference)[m_ParticuleClasse]->m_temp_evolution;
  for (size_t i=0;i<6;++i)
  {
    m_inertie[i] =
    	(*ParticuleClassesReference)[m_ParticuleClasse]->m_inertie[i] ;
    m_inertie_1[i] =
    	(*ParticuleClassesReference)[m_ParticuleClasse]->m_inertie_1[i] ;
  }

  // Lecture du tag
  fileSave >> m_tag;

  // Lecture de la transformation
  m_geoFormeVdw->readPosition2014( fileSave );

  // Lecture de l'activite
  bool actif;
  fileSave >> actif;
  m_activity = ( actif == true ) ? COMPUTE : WAIT;

  // Lecture de la cinematique
  if ( m_cinematique ) delete m_cinematique;
  m_cinematique = Cinematique_BuilderFactory::create(
  	m_geoFormeVdw->getConvex() );
  fileSave >> *m_cinematique;

  // Calcul de la masse volumique et du poids
  m_masseVolumique = m_masse / m_geoFormeVdw->getVolume();
  computeWeight();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture de la particule pour Reload. Utilise pour les particules
// dans la simulation (actives ou en attentes) pour le format de reload 2014
// en binaire
// A.WACHS - Dec 2014 - Creation
void Particule::read2014_binary( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference )
{
  // Lecture du n\B0 et de la classe de reference
  fileSave.read( reinterpret_cast<char*>( &m_id ), sizeof(int) );
  fileSave.read( reinterpret_cast<char*>( &m_ParticuleClasse ), sizeof(int) );

  // Creation de la forme (convex + transform) en utilisant le constructeur
  // de recopie
  m_geoFormeVdw = new FormeVdW(
    *(*ParticuleClassesReference)[m_ParticuleClasse]->getForme() );

  // Materiau, masse, energie et inertie de la particule de class de reference
  m_nomMateriau =
  	(*ParticuleClassesReference)[m_ParticuleClasse]->m_nomMateriau ;
  m_masse = (*ParticuleClassesReference)[m_ParticuleClasse]->m_masse ;
  m_energie = (*ParticuleClassesReference)[m_ParticuleClasse]->m_energie ;
  is_mobile = (*ParticuleClassesReference)[m_ParticuleClasse]->is_mobile ;
  m_temp_evolution = (*ParticuleClassesReference)[m_ParticuleClasse]->m_temp_evolution;

  for (size_t i=0;i<6;++i)
  {
    m_inertie[i] =
    	(*ParticuleClassesReference)[m_ParticuleClasse]->m_inertie[i] ;
    m_inertie_1[i] =
    	(*ParticuleClassesReference)[m_ParticuleClasse]->m_inertie_1[i] ;
  }

  // Lecture du tag
  fileSave.read( reinterpret_cast<char*>( &m_tag ), sizeof(int) );

  // Lecture de la transformation
  m_geoFormeVdw->readPosition2014_binary( fileSave );

  // Lecture de l'activite
  unsigned int iact = 0;
  fileSave.read( reinterpret_cast<char*>( &iact ), sizeof(unsigned int) );
  if ( iact == 1 ) m_activity = COMPUTE;
  else m_activity = WAIT;

  // Lecture de la cinematique
  if ( m_cinematique ) delete m_cinematique;
  m_cinematique = Cinematique_BuilderFactory::create(
  	m_geoFormeVdw->getConvex() );
  m_cinematique->readCineParticule2014_binary( fileSave );

  // Calcul de la masse volumique et du poids
  m_masseVolumique = m_masse / m_geoFormeVdw->getVolume();
  computeWeight();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde du Particule pour Reload
// G.FERRER - Janv.2004 - Creation
// D. RAKOTONIRINA - Sept. 2014 - Modification
void Particule::write( ostream &fileSave, Composant const* composant ) const
{
  fileSave << "\n*Particule\n";
  Composant::writeStatique( fileSave, composant );
  fileSave << "*ParticuleClasse\n"
  	   << m_ParticuleClasse << '\n';
  fileSave << "*Tag\n"
  	   << m_tag << endl;
  fileSave << "*Masse&Energie\n"
	   << Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_masse ) << " "
	   << Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_energie ) << endl;
  fileSave << "*Mobilite\n"
  	   << is_mobile << endl;
  fileSave << "*TempEvolution\n"
  	   << m_temp_evolution << endl;
  fileSave << "*Inertie" << endl
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie[0] ) << " "
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie[1] ) << " "
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie[2] ) << " "
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie[3] ) << " "
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie[4] ) << " "
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie[5] ) << endl;
  fileSave << "*Inertie_Inverse" << endl
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie_1[0] ) << " "
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie_1[1] ) << " "
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie_1[2] ) << " "
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie_1[3] ) << " "
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie_1[4] ) << " "
	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertie_1[5] ) << endl;

  Composant::writePosition( fileSave );
  fileSave << endl;
  bool b_activ = ( m_activity == COMPUTE ) ? true : false;
  fileSave << "*Activite\t" << b_activ << endl;
  m_cinematique->writeCineParticule( fileSave );
  fileSave << endl;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde de la particule pour Reload
// A.WACHS - Aout 2014 - Creation
void Particule::write2014( ostream &fileSave ) const
{
  fileSave << m_id << " " << m_ParticuleClasse << " " << m_tag << " " ;
  m_geoFormeVdw->getTransform()->writeTransform2014( fileSave );
  fileSave << " " << ( m_activity == COMPUTE ? true : false ) << " ";
  m_cinematique->writeCineParticule2014( fileSave );
  fileSave << endl;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde en binaire du Particule pour Reload
// A.WACHS - Aout 2014 - Creation
void Particule::write2014_binary( ostream &fileSave )
{
  fileSave.write( reinterpret_cast<char*>( &m_id ), sizeof(int) );
  fileSave.write( reinterpret_cast<char*>( &m_ParticuleClasse ), sizeof(int) );
  fileSave.write( reinterpret_cast<char*>( &m_tag ), sizeof(int) );
  m_geoFormeVdw->getTransform()->writeTransform2014_binary( fileSave );
  unsigned int iact = m_activity == COMPUTE ? true : false ;
  fileSave.write( reinterpret_cast<char*>( &iact ), sizeof(unsigned int) );
  m_cinematique->writeCineParticule2014_binary( fileSave );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de l'identite du monolithe
// G.FERRER - Janv.2004 - Creation
void Particule::writeIdentity( ostream &file ) const
{
  file << m_id << '\t' << this << '\n';
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Difference de vitesse de translation au temps precedent
// A.WACHS - Janvier.2009 - Creation
Vecteur Particule::getTranslationalVelocityDifferencePreviousTime() const
{
  return m_addedMassInfos->TranslationalVelocity_difference;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Difference de vitesse de rotation au temps precedent
// A.WACHS - Janvier 2009 - Creation
Vecteur Particule::getRotationalVelocityDifferencePreviousTime() const
{
  return m_addedMassInfos->RotationalVelocity_difference;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Mise a jour de la vitesse au pas de temps precedent et de la
// difference de vitesse au pas de temps precedent
// A.WACHS - Janvier 2009 - Creation
void Particule::setVelocityAndVelocityDifferencePreviousTime()
{
  m_addedMassInfos->TranslationalVelocity_difference =
  	*m_cinematique->getVitesseTranslation()
  	- m_addedMassInfos->TranslationalVelocity_nm1;
  m_addedMassInfos->TranslationalVelocity_nm1 =
  	*m_cinematique->getVitesseTranslation();
  m_addedMassInfos->RotationalVelocity_difference =
  	*m_cinematique->getVitesseRotation()
  	- m_addedMassInfos->RotationalVelocity_nm1;
  m_addedMassInfos->RotationalVelocity_nm1 =
  	*m_cinematique->getVitesseRotation();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Update the particle slip velocity at the precedent time step
//(applied in DEMCFD addedmass term)
// A.ESTEGHAMATIAN - dec 2015 - Creation
void Particule::setVelocitySlipPreviousTime(Vecteur const Vr )
{
  m_addedMassInfos->TranslationalSlipVelocity_nm1 = Vr;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Get particle slip velocity at the previous time step
// A.ESTEGHAMATIAN - dec 2015 - Creation
Vecteur Particule::getTranslationalSlipVelocityPreviousTime() const
{
  Vecteur SlipVel(0.,0.,0.);
  if ( m_addedMassInfos )
    SlipVel = m_addedMassInfos->TranslationalSlipVelocity_nm1;
  else
  {
    cout << "Error in getTranslationalSlipVelocityPreviousTime"
    " for particle tag = " << getTag() << endl;
    SlipVel = m_addedMassInfos->TranslationalSlipVelocity_nm1;
  }
  return( SlipVel );

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie this pour une particule convex et m_masterComposite pour une
// particule elemenatire de CompParticule
// D. RAKOTONIRINA - Sept 2014 - Creation
//Composant const* Particule::ReferenceComposant() const
//{
//  return( this );
//}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Mise a jour de la vitesse au pas de temps precedent
void Particule::setVelocityPreviousTimeRestart(
  	double const& vx, double const& vy, double const& vz,
  	double const& omx, double const& omy, double const& omz )
{
  m_addedMassInfos->TranslationalVelocity_nm1[X] = vx ;
  m_addedMassInfos->TranslationalVelocity_nm1[Y] = vy ;
  m_addedMassInfos->TranslationalVelocity_nm1[Z] = vz ;
  m_addedMassInfos->RotationalVelocity_nm1[X] = omx ;
  m_addedMassInfos->RotationalVelocity_nm1[Y] = omy ;
  m_addedMassInfos->RotationalVelocity_nm1[Z] = omz ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Correction ou non de la masse devant dv/dt
// M. BERNARD - Octobre 2012 - Creation
void Particule::setMassCorrection( bool is_MassCorrection )
{
  Particule::m_MassCorrection = is_MassCorrection;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Correction ou non de la masse devant dv/dt
// M. BERNARD - Octobre 2012 - Creation
bool Particule::getMassCorrection()
{
  return Particule::m_MassCorrection;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Traitement explicite ou non de la masse ajoutee
// A.WACHS - Janvier 2009 - Creation
void Particule::setExplicitMassCorrection( bool is_explicit )
{
  Particule::m_explicitAddedMass = is_explicit;
//  Particule::m_explicitMassCorrection = is_explicit;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Traitement explicite ou non de la masse ajoutee
// A.WACHS - Janvier 2009 - Creation
bool Particule::getExplicitMassCorrection()
{
  return Particule::m_explicitAddedMass;
}




// ----------------------------------------------------------------------------
// Operateur <<
ostream& operator << ( ostream &f, const Particule &P )
{
  f << "Numero = " << P.m_id << endl;
  f << "Position = " << *P.getPosition();
  f << "Tag = " << P.m_tag;
  return f;
}




// ----------------------------------------------------------------------------
// Modification du quaternion de rotation.
void Particule::setQuaternionRotation(
    const Scalar &vecteur0,
    const Scalar &vecteur1,
    const Scalar &vecteur2,
    const Scalar &scalaire )
{
  m_cinematique->setQuaternionRotation( vecteur0, vecteur1, vecteur2,
    scalaire );
}




// ----------------------------------------------------------------------------
// Modification du quaternion de rotation.
void Particule::setQuaternionRotation( const Quaternion &qrot )
{
  m_cinematique->setQuaternionRotation( qrot );
}




// ----------------------------------------------------------------------------
// Renvoie si un clone periodique existe ou non
bool Particule::periodicCloneExistence( const Vecteur &trans ) const
{
  bool exist = false;
  list<Particule*>::const_iterator ipart;
  Point translatedParticule( *getPosition() + trans );
  Point const* cloneGC;

  for (ipart=m_periodicClones.begin();ipart!=m_periodicClones.end() && !exist;
  	ipart++)
  {
    cloneGC = (*ipart)->getPosition();
    if ( fabs( (*cloneGC)[X] - translatedParticule[X] ) < EPSILON
    	&& fabs( (*cloneGC)[Y] - translatedParticule[Y] ) < EPSILON
    	&& fabs( (*cloneGC)[Z] - translatedParticule[Z] ) < EPSILON )
      exist = true;
  }

  return exist;
}




// ----------------------------------------------------------------------------
// Renvoie si un clone periodique existe ou non
bool Particule::periodicCloneExistence( const int& ObsID ) const
{
  bool exist = false;
  list<int>::const_iterator il;

  for (il=m_periodicObstaclesID.begin();il!=m_periodicObstaclesID.end()
  	&& !exist; il++)
    if ( *il == ObsID ) exist = true;

  return exist;
}




// ----------------------------------------------------------------------------
// Renvoie si un clone periodique existe ou non
bool Particule::periodicCloneExistence( Particule const* pClone ) const
{
  bool exist = false;
  list<Particule*>::const_iterator ipart;

  for (ipart=m_periodicClones.begin();ipart!=m_periodicClones.end() && !exist;
  	ipart++)
    if ( *ipart == pClone ) exist = true;

  return exist;
}




// ----------------------------------------------------------------------------
// Ajoute un clone periodique
void Particule::addPeriodicClone( Particule *clone )
{
  m_periodicClones.push_back( clone );
}




// ----------------------------------------------------------------------------
// Ajoute un vecteur de periodicite de clone periodique
void Particule::addPeriodicObstacleID( int const& ObsID )
{
  m_periodicObstaclesID.push_back( ObsID );
}




// ----------------------------------------------------------------------------
// Renvoie le clone dont la reference est en contact avec l'obstacle
// et l'efface de la liste
Particule* Particule::getPeriodicCloneAndErase( const ObstaclePeriodique* obs )
{
  Particule* part=NULL;
  list<Particule*>::iterator ipart;
  for (ipart=m_periodicClones.begin();ipart!=m_periodicClones.end();)
    if ( (*ipart)->getObstacle() == obs )
    {
      part=*ipart;
      ipart = m_periodicClones.erase(ipart);
    }
    else ipart++;

  return part;
}




// ----------------------------------------------------------------------------
// Gestion des clones multi-periodiques
// D. RAKOTONIRINA - Juil 2014 - Modification
void Particule::createMultiPeriodicClones(
	list<Particule*>* particulesClonesPeriodiques,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC )
{
  // Rem: Qd cette methode est appelee, il est suppose que periodicClones est
  // soit vide soit contient uniquement des clones uni-periodiques
  size_t nbClones = m_periodicClones.size();
  if ( nbClones == 2 )
  {
    list<Particule*>::iterator clone;
    Vecteur periode;

    // Vecteur de bi-periodicite
    for (clone=m_periodicClones.begin();clone!=m_periodicClones.end();clone++)
      periode += *((*clone)->getVecteurPeriodique());

    // Creation du clone bi-periodique
      if ( isCompParticule() )
      {
	Particule* biperclone = new CompParticulePeriodique( this,
	    (*ParticuleClassesReference)[m_ParticuleClasse],
	    NULL, periode, 2 );

	// Ajout a la liste des clones de la particule
    	// et a la liste de tous les clones
    	m_periodicClones.push_back(biperclone);
    	particulesClonesPeriodiques->push_back(biperclone);

    	// Link avec le LinkedCell
    	LC->Link( biperclone );
      }
      else
      {
	Particule* biperclone = new ParticulePeriodique( this,
	    (*ParticuleClassesReference)[m_ParticuleClasse],
	    NULL, periode, 2 );
	// Ajout a la liste des clones de la particule
    	// et a la liste de tous les clones
    	m_periodicClones.push_back(biperclone);
    	particulesClonesPeriodiques->push_back(biperclone);

    	// Link avec le LinkedCell
    	LC->Link( biperclone );
      }
  }
  else if ( nbClones == 3 )
  {
    list<Particule*>::iterator clone;
    Vecteur periode;
    Particule* pNullPart = NULL;
    vector<Particule*> uniclones( 3, pNullPart );
    int i=0;
    for (clone=m_periodicClones.begin();clone!=m_periodicClones.end();
    	clone++,i++) uniclones[i] = *clone;

    // Creation des 3 clones bi-periodiques
    pair<Particule*,Particule*> pairclones;
    list< pair<Particule*,Particule*> > allpairs;
    pairclones.first = uniclones[0];
    pairclones.second = uniclones[1];
    allpairs.push_back(pairclones);
    pairclones.first = uniclones[0];
    pairclones.second = uniclones[2];
    allpairs.push_back(pairclones);
    pairclones.first = uniclones[1];
    pairclones.second = uniclones[2];
    allpairs.push_back(pairclones);
    list< pair<Particule*,Particule*> >::iterator ilp;

    for (ilp=allpairs.begin();ilp!=allpairs.end();ilp++)
    {
      // Vecteur de bi-periodicite
      periode = *(ilp->first->getVecteurPeriodique())
      	+ *(ilp->second->getVecteurPeriodique());

      if ( isCompParticule() )
      {
	// Creation du clone bi-periodique
      	Particule* biperclone = new CompParticulePeriodique( this,
      	  (*ParticuleClassesReference)[m_ParticuleClasse],
      	  NULL, periode, 2 );

      	// Initialisation de son torseur de force a zero
      	biperclone->InitializeForce( false );

      	// Ajout a la liste des clones de la particule
      	// et a la liste de tous les clones
      	m_periodicClones.push_back(biperclone);
      	particulesClonesPeriodiques->push_back(biperclone);

      	// Link avec le LinkedCell
      	LC->Link( biperclone );
      	}
      else
      {
	// Creation du clone bi-periodique
      	Particule* biperclone = new ParticulePeriodique( this,
      	  (*ParticuleClassesReference)[m_ParticuleClasse],
      	  NULL, periode, 2 );

      	// Initialisation de son torseur de force a zero
      	biperclone->InitializeForce( false );

      	// Ajout a la liste des clones de la particule
      	// et a la liste de tous les clones
      	m_periodicClones.push_back(biperclone);
      	particulesClonesPeriodiques->push_back(biperclone);

      	// Link avec le LinkedCell
      	LC->Link( biperclone );
      	}
    }

    // Creation du clone tri-periodique
    // Vecteur de tri-periodicite
    periode = *(uniclones[0]->getVecteurPeriodique())
      	+ *(uniclones[1]->getVecteurPeriodique())
      	+ *(uniclones[2]->getVecteurPeriodique());

    if ( isCompParticule() )
    {
      // Creation du clone tri-periodique
      Particule* triperclone = new CompParticulePeriodique( this,
          (*ParticuleClassesReference)[m_ParticuleClasse],
          NULL, periode, 3 );

      // Initialisation de son torseur de force a zero
      triperclone->InitializeForce( false );

      // Ajout a la liste des clones de la particule
      // et a la liste de tous les clones
      m_periodicClones.push_back(triperclone);
      particulesClonesPeriodiques->push_back(triperclone);

      // Link avec le LinkedCell
      LC->Link( triperclone );
    }
    else
    {
      // Creation du clone tri-periodique
      Particule* triperclone = new ParticulePeriodique( this,
          (*ParticuleClassesReference)[m_ParticuleClasse],
          NULL, periode, 3 );

      // Initialisation de son torseur de force a zero
      triperclone->InitializeForce( false );

      // Ajout a la liste des clones de la particule
      // et a la liste de tous les clones
      m_periodicClones.push_back(triperclone);
      particulesClonesPeriodiques->push_back(triperclone);

      // Link avec le LinkedCell
      LC->Link( triperclone );
    }
  }
}




// ----------------------------------------------------------------------------
// Supprime le clone de la liste
void Particule::erasePeriodicClone( Particule *clone )
{
  list<Particule*>::iterator il;
  for (il=m_periodicClones.begin();il!=m_periodicClones.end(); )
    if ( *il == clone ) il = m_periodicClones.erase(il);
    else il++;
}




// ----------------------------------------------------------------------------
// Supprime l'obstacle periodique de la liste
void Particule::erasePeriodicObstacleID( int const& ObsID )
{
  list<int>::iterator il;
  for (il=m_periodicObstaclesID.begin();il!=m_periodicObstaclesID.end(); )
    if ( *il == ObsID ) il = m_periodicObstaclesID.erase(il);
    else il++;
}




// ----------------------------------------------------------------------------
// Translate les clones et met a jour leur lien dans le LinkedCell
void Particule::translateAndUpdatePeriodicClones( Vecteur const& translation,
	LinkedCell* LC )
{
  for (list<Particule*>::iterator clone=m_periodicClones.begin();
  	clone!=m_periodicClones.end();clone++)
  {
    (*clone)->Translate( translation );
    LC->LinkUpdateActiveParticule( *clone );
  }
}




// ----------------------------------------------------------------------------
// Mise a jour de la localisation geographique de la particule
// dans le LinkedCell en utilisant cellule_nm1. Si la methode est appelee apres
// le LinkUpdate du LinkedCell, cellule_nm1 = cellule courante a laquelle
// appartient la particule
void Particule::updateGeoLocalisation()
{
  m_GeoLoc = m_cellule_nm1->getGeoLocalisation();
}




// ----------------------------------------------------------------------------
// Localisation geographique */
MPIGeoLocalisation Particule::getGeoLocalisation() const
{
  return m_GeoLoc;
}




// ----------------------------------------------------------------------------
// La particule possede un clone dans une direction donnee ?
bool Particule::hasCloneInDirection( Vecteur const* direction,
	LinkedCell const* LC ) const
{
  bool found = false;

  // Meme direction pour 2 vecteurs => produit vectoriel de norme 0

  // Version sequentielle
  if ( !m_periodicClones.empty() )
  {
    Vecteur const* vp = NULL;
    for (list<Particule*>::const_iterator il=m_periodicClones.begin();
    	il!=m_periodicClones.end() && !found;il++)
    {
      vp = (*il)->getVecteurPeriodique();
      if ( Norm((*direction) ^ (*vp)) < EPS ) found = true;
    }
  }
  // Version parallele
  else if ( !m_periodicObstaclesID.empty() )
  {
    for (list<int>::const_iterator ip=m_periodicObstaclesID.begin();
    	ip!=m_periodicObstaclesID.end() && !found;ip++)
    {
      if ( Norm((*direction) ^ *(LC->getObstaclePeriodique(*ip)->getPeriode()) )
       < EPS ) found = true;
    }
  }

  return found;
}




// ----------------------------------------------------------------------------
// Creer la structure AddedMassInfos
void Particule::createAddedMassInfos()
{
  if ( !m_addedMassInfos )
  {
    m_addedMassInfos = new struct AddedMassInfos;
    m_addedMassInfos->TranslationalSlipVelocity_nm1.setValue(0.,0.,0.);
  }
}




// ----------------------------------------------------------------------------
// Renvoie un vecteur orientation de la particule
Vecteur Particule::vecteurOrientation() const
{
  return( m_geoFormeVdw->getConvex()->vecteurOrientation(
  	m_geoFormeVdw->getTransform() ) );
}




// ----------------------------------------------------------------------------
// Renvoie le nombre de contacts de la particule
int Particule::getCoordinationNumber() const
{
  return m_coordination_number;
}




// ----------------------------------------------------------------------------
// Ajoute un nombre de contacts au nombre de contacts de la particule;
// Utilisation: ajoute a la particule de reference periodique les contacts
// de son clone periodique
void Particule::addToCoordinationNumber( int const& nc )
{
  m_coordination_number += nc;
}




// ----------------------------------------------------------------------------
// Initialise le torseur des efforts sur le composant
void Particule::InitializeForce( bool const& withWeight )
{
  if ( withWeight )
    m_somme.setToBodyForce( *m_geoFormeVdw->getCentre(), m_weight );
  else m_somme.setToBodyForce( *m_geoFormeVdw->getCentre(), VecteurNul );
  m_coordination_number = 0 ;

  if( Grains_Exec::m_withFluidTemperature ||
      Grains_Exec::m_withSolidTemperature )
    m_sum_HeatFlux = 0.;

}




// ----------------------------------------------------------------------------
// Cree la structure DEM-CFD_FluidInfos
void Particule::allocateDEMCFD_FluidInfos()
{
  if( !m_fluidInfos )
  {
    m_fluidInfos = new struct DEMCFD_FluidInfos;
    setVitesseTr_fluide( 0., 0., 0. );
    set_DEMCFD_volumeFraction( 1. );

    if( Grains_Exec::m_withFluidTemperature )
      set_DEMCFD_fluidTemperature( m_fluidInitialTemperature );
  }
}




// ----------------------------------------------------------------------------
// Affectation de la massse volumique
void Particule::setMasseVolumique( double const& density_ )
{
  m_masse *= density_ / m_masseVolumique ;
  for (int i=0; i<6; i++)
  {
    m_inertie[i] *= density_ / m_masseVolumique;
    m_inertie_1[i] /= density_ / m_masseVolumique;
  }
  m_masseVolumique = density_ ;
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
void Particule::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  m_geoFormeVdw->getConvex()->write_polygonsStr_PARAVIEW(connectivity,
	offsets, cellstype, firstpoint_globalnumber, last_offset);
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
void Particule::write_polygonsPts_PARAVIEW( ostream &f,
	const Transform &transform, Vecteur const* translation ) const
{
  m_geoFormeVdw->getConvex()->write_polygonsPts_PARAVIEW(f,
	transform, translation);
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
int Particule::numberOfPoints_PARAVIEW() const
{
  return ( m_geoFormeVdw->getConvex()->numberOfPoints_PARAVIEW() );
}




// ----------------------------------------------------------------------------
// Nombre de polygones elementaires pour post-processing avec Paraview
int Particule::numberOfCells_PARAVIEW() const
{
  return ( m_geoFormeVdw->getConvex()->numberOfCells_PARAVIEW() );
}




// ----------------------------------------------------------------------------
// Nombre de polygones elementaires pour post-processing avec Paraview
list<Point> Particule::get_polygonsPts_PARAVIEW( Vecteur const* translation )
const
{
  return ( getForme()->get_polygonsPts_PARAVIEW( translation ) );
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
void Particule::write_polygonsPts_PARAVIEW( ostream &f,
  	Vecteur const* translation ) const
{
  getForme()->write_polygonsPts_PARAVIEW( f, translation );
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de la position de la particule pour Basilisk
void Particule::writePositionInFluid (BasiliskDataStructure * b)
{
  m_geoFormeVdw->writePositionInFluid (b);
}



// ----------------------------------------------------------------------------
// Initialise a faux le boolean correspondant au calcul de la
// transformation avec scaling par l'epaisseur de croute
// D. RAKOTONIRINA - Fev. 2014 - Creation
void Particule::initializeVdWtransform_to_notComputed()
{
  getForme()->initializeVdWtransform_to_notComputed();
}




// ----------------------------------------------------------------------------
// Utile pour le typage dynamique (acces a la methode dans CompParticule)
// Renvoie les limites du domaine pour la discretisation de la
// particule composite afin d'approximer les elements d'inertie
// D. RAKOTONIRINA - Mars 2014 - Creation
vector<Scalar> Particule::getCompFeaturesBox()
{
  vector<Scalar> Vdim(6, 0.);

  return( Vdim );
}




// ----------------------------------------------------------------------------
// Mouvement aleatoire sur les particules actives
// D. RAKOTONIRINA - Nov. 2014 - Creation
void Particule::setRandomMotion( double const& coefTrans,
	double const& coefRot )
{
  if ( coefTrans )
    if ( m_tag != 2 )
    {
      Vecteur rvel;
      rvel[X] = coefTrans * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      rvel[Y] = coefTrans * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      if ( Grains_BuilderFactory::getContext() == DIM_3 )
        rvel[Z] = coefTrans * (
	      	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      setVitesseTranslation( rvel );
    }

  if ( coefTrans )
    if ( m_tag != 2 )
    {
      Vecteur rvel;
      rvel[X] = coefTrans * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      rvel[Y] = coefTrans * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      if ( Grains_BuilderFactory::getContext() == DIM_3 )
        rvel[Z] = coefTrans * (
	      	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      setVitesseRotation( rvel );
    }
}




// ----------------------------------------------------------------------------
// Utile pour le typage dynamique (acces a la methode dans CompParticule)
// !!! Seulement utilisee dans le cas d'une CompParticule !!!
// D. RAKOTONIRINA - Avril 2014 - Creation
vector<ElementParticule*> Particule::getElementParticules() const
{
  cout << "WARNING!!! Particule::getElementParticules()"
       << "is not implemented!\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
}




// ----------------------------------------------------------------------------
// Utile pour le typage dynamique (acces a la methode dans CompParticule)
// !!! Seulement utilisee dans le cas d'une CompParticule !!!
// D. RAKOTONIRINA - Avril 2014 - Creation
vector<Vecteur> Particule::getInitialRelativePositions() const
{
  cout << "WARNING!!! Particule::getInitialRelativePositions()"
       << "is not implemented!\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
}




// ----------------------------------------------------------------------------
// Utile pour le typage dynamique (acces a la methode dans CompParticule)
// !!! Seulement utilisee dans le cas d'une CompParticule !!!
// D. RAKOTONIRINA - Avril 2014 - Creation
vector<Vecteur> Particule::getRelativePositions() const
{
  cout << "WARNING!!! Particule::getRelativePositions()"
       << "is not implemented!\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les matrices de rotations intiales des particules elementaires
// Utile pour le typage dynamique (acces a la methode dans CompParticule)
// !!! Seulement utilisee dans le cas d'une CompParticule !!!
// D. RAKOTONIRINA - Avril 2014 - Creation
vector<Matrix> Particule::getInitialMatrix() const
{
  cout << "WARNING!!! Particule::getInitialMatrix()"
       << "is not implemented!\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre total de particules elementaires
// Utile pour le typage dynamique (acces a la methode dans CompParticule)
// !!! Seulement utilisee dans le cas d'une CompParticule !!!
// D. RAKOTONIRINA - Juin 2014 - Creation
size_t Particule::getNbreElemPart() const
{
  cout << "WARNING!!! Particule::getNbreElemPart()"
       << "is not implemented!\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
}




// ----------------------------------------------------------------------------
// Positionne les particules elementaires dans l'espace
// Utile pour le typage dynamique (acces a la methode dans CompParticule)
// !!! Seulement utilisee dans le cas d'une CompParticule !!!
// D. RAKOTONIRINA - Juin 2014 - Creation
void Particule::setElementPosition()
{
  cout << "WARNING!!! Particule::setElementPosition()"
       << "is not implemented!\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
}




// ----------------------------------------------------------------------------
// Renvoie le nom associe a la particule composite
// Utile pour le typage dynamique (acces a la methode dans CompParticule)
// !!! Seulement utilisee dans le cas d'une CompParticule !!!
// D. RAKOTONIRINA - Sept. 2014 - Creation
string Particule::getPartName() const
{
  cout << "WARNING!!! Particule::getPartName()"
       << "is not implemented!\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de sommets ou un code equivalent
// D. RAKOTONIRINA - Avril. 2015 - Creation
int Particule::getNbCorners() const
{
  return getForme()->getConvex()->getNbCorners();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de la position de la particule pour le Fluide
// D. RAKOTONIRINA - Avril. 2015 - Creation
void Particule::writePositionInFluid( ostream &fileOut )
{
  m_geoFormeVdw->writePositionInFluid( fileOut );
}




// ----------------------------------------------------------------------------
// Renvoie le nombre de sommets ou un code equivalent
// M. SULAIMAN - Nov.2015 - Creation
Scalar Particule::set_shrinking_mass()
{
  m_masse = m_masseVolumique * m_geoFormeVdw->getVolume();
  return( m_masse );
}
