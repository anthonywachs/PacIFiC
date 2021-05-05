#ifndef _Particule
#define _Particule

#include "Composant.H"
#include "CineParticule.H"
#include "PointContact.hh"
#include "Vecteur.H"
#include "Cellule.H"
#include "ReaderXML.hh"
#include "WriterXML.hh"
using namespace solid;


class Cellule;
class Obstacle;
class MonObstacle;
class ObstaclePeriodique;
class Quaternion;
class LinkedCell;
class ElementParticule;
struct PointForcePostProcessing;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Description des activites de Particule
    @author GRAINS Project - IFP - 2007 */
enum ParticuleActivity
  {
    /** @brief En Attente */
    WAIT,
    /** @brief En Simulation */
    COMPUTE,
    /** @brief Bascule de Simulation vers Attente */
    CLEARandWAIT
  };



/** @brief Stockage des vitesses et des differences de vitesse au temps fluide
    precedent pour calcul explicite de la masse ajoutee

@author A.WACHS - IFP Energies nouvelles - 2011 - Creation */
// ============================================================================
struct AddedMassInfos
{
  Vecteur TranslationalVelocity_nm1; /**< Vitesse de translation au pas de
  	temps fluide pr�c�dent */
  Vecteur TranslationalVelocity_difference; /**< Difference de vitesse de
  	translation au pas de temps fluide pr�c�dent */
  Vecteur RotationalVelocity_nm1; /**< Vitesse de rotation au pas de temps
  	fluide pr�c�dent */
  Vecteur RotationalVelocity_difference; /**< Difference de vitesse de rotation
  	au pas de temps fluide pr�c�dent */
  Vecteur TranslationalSlipVelocity_nm1; /**< Vitesse de translation relative
        au pas de temps fluide pr�c�dent (utilis� dans le method DEMCFD)*/
};



/** @brief Fluid informations
@author M.BERNARD - Creation */
// ============================================================================
struct DEMCFD_FluidInfos
{
  Vecteur m_vitesseTr_fluide; /** Vitesse de translation du fluide
         (default 0.) */
  Vecteur m_vorticity_fluide; /** fluid vorticity
         (default 0.) */
  Vecteur m_gradientPression_fluide; /** Gradient de pression fluide
         (default 0.) */
  double m_demcfd_epsilon; /** fraction volumique (default 1.) */
  Vecteur m_hydroForce; /** Force hydro (default 0.) */
  Vecteur m_slipvel; /** Slip Velocity (default 0.) */
  double m_demcfd_fluidTemperature; /** fluid temperature (default 0.) */
  double fluidSolid_heatFlux; /** fluid temperature (default 0.) */
};


/** @brief Une Particule est un composant possedant une cinematique dependant
    des forces calculees par la simulation au sein des differents objets APP
    auxquels est associee la Particule.

    La Particule correspond a l'element principal du logiciel.
    Elle contient sa position physique dans le probleme ainsi que la
    cinematique associee. Elle sait comment interagir avec d'autres
    particules ou des composants (obstacles).
    Une particule est le seul composant inserable en cours de simulation.
    Pour cela est possede un etat : actif ou inactif (par defaut : inactif).

    @author F.PRADEL - Institut Francais du Petrole - 1999 - Creation */
// ============================================================================
class Particule : public Composant
{
public:
  /**@name Constructors & Destructor */
  //@{
  /** @brief Constructeur par defaut
  @param autonumbering numerotation automatique ou non */
  Particule( const bool &autonumbering = true );

  /** @brief Constructeur
  @param root Noeud "<Particule>"
  @param autonumbering numerotation automatique ou non
  @param pc classe de particules */
  Particule( DOMNode* root, const bool &autonumbering = true,
  	const int &pc = 0 );

  /** @brief Constructeur avec arguments
  @param id_ numero
  @param ParticuleRef particule de r�f�rence
  @param vx composante x de la vitesse de translation
  @param vy composante y de la vitesse de translation
  @param vz composante z de la vitesse de translation
  @param rx composante x de la vitesse de rotation
  @param ry composante y de la vitesse de rotation
  @param rz composante z de la vitesse de rotation
  @param qrotationx composante x du quaternion de rotation
  @param qrotationy composante y du quaternion de rotation
  @param qrotationz composante z du quaternion de rotation
  @param qrotations scalaire du quaternion de rotation
  @param m position & configuration de la particule
  @param activ activite
  @param tag_ tag
  @param coordination_number_ nombre de contacts de la particule */
  Particule( const int &id_, Particule const* ParticuleRef,
	const double &vx, const double &vy, const double &vz,
	const double &qrotationx, const double &qrotationy,
	const double &qrotationz, const double &qrotations,
	const double &rx, const double &ry, const double &rz,
	const Scalar m[16],
	const ParticuleActivity &activ,
	const int &tag_,
	const int &coordination_number_ = 0 );

  /** @brief Constructeur avec arguments
  @param id_ numero
  @param ParticuleRef particule de r�f�rence
  @param vtrans vitesse de translation
  @param vrot vitesse de rotation
  @param qrot quaternion de rotation
  @param config transformation li�e a la particule
  @param activ activite
  @param tag_ tag */
  Particule( const int &id_, Particule const* ParticuleRef,
	const Vecteur &vtrans,
	const Quaternion &qrot,
	const Vecteur &vrot,
	const Transform &config,
	const ParticuleActivity &activ,
	const int &tag_ );

  /** @brief Constructeur par Copie.
  Le torseur des forces appliquees est vide.
  Les types de Forme et de Cinematique sont instancies.
  @param copie Particule de reference. */
  Particule( const Particule &copie );

  /** @brief Destructeur */
  virtual ~Particule();
  //@}


  /** @name Methods Static */
  //@{
  /** @brief Masse volumique du fluide en interaction avec les particules */
  static double getFluideMasseVolumique();

  /** @brief Masse volumique du fluide en interaction avec les particules
  @param rho masse volumique du fluide */
  static void setFluideMasseVolumique( double rho );

  /** @brief Correction ou non de la masse devant dv/dt
  @param is_MassCorrection vrai si Correction de la masse, faux sinon */
  static void setMassCorrection( bool is_MassCorrection );

  /** @brief Correction ou non de la masse devant dv/dt */
  static bool getMassCorrection();

  /** @brief Traitement explicite ou non de la masse ajout�e
  @param is_explicit vrai si explicite, faux sinon */
  static void setExplicitMassCorrection( bool is_explicit );

  /** @brief Traitement explicite ou non de la masse ajout�e */
  static bool getExplicitMassCorrection();

  /** @brief Viscosite du fluide en interaction avec les particules */
  static double getFluidViscosity();

  /** @brief Viscosite du fluide en interaction avec les particules
  @param mu viscosite du fluide */
  static void setFluidViscosity( double mu );

  /** @brief Set fluid initial temperature
  @param tempF Fluid temperature */
  static void setFluidInitialTemperature( double tempF );
  //@}


  /** @name Methods Set */
  //@{
  /** @brief Remise a zero de la particule.
  Usage : suppression du particule pour reinsertion */
  void reset();

  /** @brief Remise a zero de la cinematique. */
  void ResetCinematique();

  /** @brief Affectation d'un champ de vitesse de rotation.
  @param rotation vecteur de rotation */
  virtual void setVitesseRotation( const Vecteur &rotation );

  /** @brief Affectation d'un champ de vitesse de translation.
  @param translation vecteur de translation */
  virtual void setVitesseTranslation( const Vecteur &translation );

  /** @brief Indication de la nouvelle activite de la particule.
  @param activity Le nouvel etat */
  void setActivity( ParticuleActivity activity ) { m_activity = activity; }

  /** @brief Affectation du type a la Particule.
  Dans certaine configuration (Periodique) celui-ci peut etre modifie.
  @param type Le nouveau type (defaut : "P") */
  void setType( const string &type = "P" ) { m_type = type; }

  /** @brief Mise a jour de la vitesse au pas de temps pr�c�dent et de la
  difference de vitesse au pas de temps pr�c�dent */
  void setVelocityAndVelocityDifferencePreviousTime();

  /** @brief Update the particle slip velocity at the precedent time step
  (applied in DEMCFD addedmass term) */
  void setVelocitySlipPreviousTime(Vecteur Vr);


  /** @brief Mise a jour de la vitesse au pas de temps pr�c�dent
  @param vx composante x de la vitesse translationnelle
  @param vy composante y de la vitesse translationnelle
  @param vz composante z de la vitesse translationnelle
  @param omx composante x de la vitesse rotationnelle
  @param omy composante y de la vitesse rotationnelle
  @param omz composante z de la vitesse rotationnelle */
  void setVelocityPreviousTimeRestart(
      double const& vx, double const& vy, double const& vz,
      double const& omx, double const& omy, double const& omz );

  /** @brief Affectation du tag de la particule
  @param tag_ tag */
  int setTag( int tag_ ) { m_tag = tag_; return m_tag; }

  /** @brief Affectation de la localisation geographique de la particule
  dans le LinkedCell
  @param geoloc_ localisation geographique de la particule dans le LinkedCell */
  void setGeoLocalisation( MPIGeoLocalisation const& geoloc_ )
  	{ m_GeoLoc = geoloc_; }

  /** @brief Cellule dans laquelle la particule se situait au
  pas de temps pr�c�dent
  @param cel la celulle au pas dz temps pr�c�dent */
  void setCelluleNm1( Cellule* cel ) { m_cellule_nm1 = cel; }

  /** @brief Modification du quaternion de rotation.
  @param vecteur0 Composante 0 du vecteur
  @param vecteur1 Composante 1 du vecteur
  @param vecteur2 Composante 2 du vecteur
  @param scalaire Composante du scalaire */
  virtual void setQuaternionRotation( const Scalar &vecteur0,
	const Scalar &vecteur1,
	const Scalar &vecteur2,
	const Scalar &scalaire );

  /** @brief Modification du quaternion de rotation.
  @param qrot quaternion de rotation */
  virtual void setQuaternionRotation( const Quaternion &qrot );

  /** @brief Affectation de la classe de la particule
  @param pc classe de particules */
  void setParticuleClasse( const int &pc ) { m_ParticuleClasse = pc; }

  /** @brief Cinematique au temps t-2dt: vitesse de translation,
  vitese de rotation, variation de QDM translationnelle, variation de QDM
  rotationnalle
  @param tab tableau de 12 doubles contenant les valeurs des 4 vecteurs vitesse
  de translation, vitese de rotation, variation de QDM translationnelle,
  variation de QDM rotationnalle */
  void setCinematiqueNm2( double const* tab );

  /**
    @brief Affectation d'un champ de vitesse de translation du fluide sur
      la particule
    @details Methode par passage d'adresse
    @param vitesseTr_fluide_x composante x de la vitesse fluide
    @param vitesseTr_fluide_y composante y de la vitesse fluide
    @param vitesseTr_fluide_z composante z de la vitesse fluide
  */
  void setVitesseTr_fluide(
  	double const& vitesseTr_fluide_x,
  	double const& vitesseTr_fluide_y,
	double const& vitesseTr_fluide_z);

  /** @brief Affectation d'un champ de vorticite du fluide sur
  la particule
  @param vorticity_fluide_x first component of the fluid vorticity
  @param vorticity_fluide_y second component of the fluid vorticity
  @param vorticity_fluide_z third component of the fluid vorticity */
  void setVorticity_fluide(
    double const& vorticity_fluide_x,
    double const& vorticity_fluide_y,
    double const& vorticity_fluide_z);

  /** @brief Affectation d'un gradient de pression sur
  la particule
  @param gradP_fluid_x composante x du gradient de pression
  @param gradP_fluid_y composante y du gradient de pression
  @param gradP_fluid_z composante z du gradient de pression */
  void setGradientPression_fluide(
  	double const& gradP_fluid_x,
  	double const& gradP_fluid_y,
	double const& gradP_fluid_z);

  /** @brief Affectation de la fraction volumique interpolee sur
  la particule
  @param fluidVolumeFraction fraction volumique */
  void set_DEMCFD_volumeFraction( double const& fluidVolumeFraction );

  /**
    @brief Set DEMCFD surounding fluid temperature
  */
  void set_DEMCFD_fluidTemperature( const Scalar fluidTemperature_ );

  /** @brief Affectation de la force hydrodynamique
  @param hydroForce_ force hydro sur la particule */
  void setHydroForce( Vecteur const* hydroForce_ );

//
  /** @brief Keep hydroforce for post-processing
  @param hydroForce_ force hydro sur la particule */
  void setHydroForce_PP( Vecteur const* hydroForce_ );

  /** @brief Affectation de la slip velocity
  @param slipVel_ particle slip velocity */
  void setSlipVel( Vecteur const* slipVel_ );

  /** @brief Affectation de la massse volumique
  @param density_  masse volumique */
  void setMasseVolumique( double const& density_ );

  /** @brief Mouvement aleatoire sur les particules actives
  @param coefTrans amplitude max du mouvement translationnel
  @param coefRot amplitude max du mouvement rotationnel	*/
  void setRandomMotion( double const& coefTrans,
	double const& coefRot );

  /**
    @brief Set solid-body temperature
    @param fluidTemperature_ Solid-body temperature
  */
  void set_solidTemperature( const Scalar solidTemperature_ );

  /**
    @brief Set fluid-solid heat flux to m_infos structure
    @param heatFlux_ Fluid-Solid heat flux computed in AppFluide_Temperature
    @remarks Called in AppFluide_Temperature::compute_fluidHeatFlux()
    @remarks This heat flux will be used to compute backard volumic
        heat flux on the fluid phase
  */
  void set_fluidSolidHeatFlux( const double heatFlux_ );

   /**
    @brief Set solid-body Nusselt number
    @param fluidTemperature_ Solid-body temperature
  */
  void set_solidNusselt( const Scalar Nusselt_ );




  /**
    @brief Set rnd variables required for stochastic drag law to buffer
    or at AppFluid_Drag.cc
    @param first comp
    @param 2nd comp
    @param 3rd comp
   */
  void set_rnd( double const& g_rnd_x,
                double const& g_rnd_y,
		double const& g_rnd_z );
  /**
    @brief Set rnd variable required for stochastic Nusselt number law to buffer
    or at AppFluid_Temperature.cc
    @param first comp
   */
  void set_rnd_Nu( double const& g_rnd_Nux,
                   double const& g_rnd_Nuy,
                   double const& g_rnd_Nuz );

  //@}

  /** @name Methods */
  //@{
  /** @brief Resolution des equations de la dynamique et double integration
  pour obtenir la nouvelle vitesse et position.
  <ul>
    <li> v^{n+0.5} & =  &v_i^{n-0.5} + Dt * f_{total}(i)^{n}
    <li> x^{n+1}_i & = & x^{n}_i + Dt * v^{n+0.5}_i
  </ul>
  De plus, cette fonction ne vide pas la liste des torseurs ;
  cette action doit etre realisee par la methode "clearForces()".
  @exception ErreurDeplacement si le deplacement est trop grand.
  @param temps Temps de simulation
  @param dt Pas de temps de resolution */
  virtual void Deplacer( Scalar temps, double dt ) throw(ErreurDeplacement);

  /**
    @brief Solve temperature problem
    @param temps Temps de simulation.
    @param dt Pas de temps.
  */
  void ComputeTemperature( Scalar temps, Scalar dt ) ;

  /** @brief Contact entre la particule et une particule.
  Si le contact est realise, on ajoute le contact au deux composants.
  @exception ErreurContact Si un choc entre croute est detecte
  @param voisin Le composant (monolithe) a etudier
  @param dt Increment de temps
  @param temps temps physique
  @param LC grille de cellules */
  virtual void InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC ) throw (ErreurContact);

  /**
    @brief Contact entre la particule et une particule pour l'initialisation
        des conditions de cohesion. (similaire a InterAction)
    @exception ErreurContact Si un choc entre croute est detecte
    @param voisin Le composant (monolithe) a etudier
    @param dt Increment de temps
    @param temps temps physique
    @param LC grille de cellules
  */
  virtual void InterActionCohesiveInit( Composant* voisin,
      double dt, double const& temps, LinkedCell *LC ) throw (ErreurContact);

  /** @brief Contact entre la particule et une particule pour le post processing
  Si le contact est realise, on ajoute le contact au deux composants.
  @param voisin Le composant (monolithe) a etudier
  @param listOfContacts liste de contacts */
  virtual void InterActionPostProcessing( Composant* voisin, Scalar dt,
	list<struct PointForcePostProcessing>* listOfContacts )
	throw (ErreurContact);

  /** @brief Copie du vitesse de rotation dans le vecteur vit en d�butant � la
  position i
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  virtual void copyVitesseRotation( double *vit, int i ) const;

  /** @brief Copie du vitesse de translation dans le vecteur vit en d�butant �
  la position i
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  virtual void copyVitesseTranslation( double *vit, int i ) const;

  /** @brief Copie du quaternion de rotation dans le vecteur vit en d�butant �
  la position i
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  virtual void copyQuaternionRotation( double *vit, int i ) const;

  /** @brief Copie force & moment exerc�s sur la particule dans le vecteur fm
  en d�butant � la position i
  @param fm vecteur de copie
  @param i position dans le vecteur fm */
  void copyForceMoment( double *fm, int i ) const;

  /** @brief Copie de la cinematique au temps t-2dt: vitesse de translation,
  vitese de rotation, variation de QDM translationnelle, variation de QDM
  rotationnalle
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  virtual void copyCinematiqueNm2(double *vit,int i) const;

  /** @brief Copie de la vitesse de translation, du gradient de pression et
  de la fraction volumique du fluide dans le vecteur vit en d�butant �
  la position i
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  void copyFluidInformations( double *vit, int i ) const;

  /**
    @brief Copie de la vorticite du fluide
    @param vit vecteur de copie
    @param i position dans le vecteur vit
  */
  void copyFluidVorticity( double *vit, int i ) const;

  /**
    @brief Copy Fluid Temperature
    @param vit vecteur de copie
    @param i position dans le vecteur vit
  */
  void copy_fluidTemperature( double *vit, int i ) const;

  /**
    @brief Copy Solid Temperature
    @param vit vecteur de copie
    @param i position dans le vecteur vit
  */
  void copy_solidTemperature( double *vit, int i ) const;

  /**
    @brief Copy Solid Nusselt number
    @param vit vecteur de copie
    @param i position dans le vecteur vit
  */
  void copy_solidNusselt( double *vit, int i ) const;

  /**
   Copy rnd number for stochastic drag into buffer
   @param vit vecteur de copie
   @param i position dans le vecteur vit
   */
  void copy_rnd(double *vit, int i ) const;

     /**
   Copy rnd number for stochastic Nusselt into buffer
   @param vit vecteur de copie
   @param i position dans le vecteur vit
   */
  void copy_rnd_Nu(double *vit, int i ) const;

  /** @brief Ajout d'une force au torseur des efforts
  @param force la nouvelle force
  @param point point d'action de la force */
  virtual void addForce( const Point &point, const Vecteur &force );

  //@{
  /** @brief Force at the contact point for the stress tensor
  @param force la nouvelle force
  @param point point d'action de la force */
  virtual void computeInternalMoments( const Point &point, const Vecteur &force );

  /** @brief Add Contact force for postprocessing purposes
  @param force the new force */
  virtual void addContactForcePP( const Vecteur &force );

  /** @brief Add Contact force for postprocessing purposes
  @param force the new force */
  virtual void addContactForcePP_instantaneous( const Vecteur &force );

  /** @brief Add Contact force for postprocessing purposes
  @param force the new force */
  virtual void addLubriForcePP( const Vecteur &force );

  /** @brief Ajoute une force et un moment au torseur des efforts exerces sur la
  particule (utile en simulation periodique)
  @param fx composante x de la force
  @param fy composante y de la force
  @param fz composante z de la force
  @param mx composante x du moment
  @param my composante y du moment
  @param mz composante z du moment */
  virtual void addForceMoment( const double &fx, const double &fy,
  	const double &fz, const double &mx, const double &my,
	const double &mz ) ;

  /** @brief Renvoie si un clone periodique existe ou non
  @param trans vecteur de translation */
  bool periodicCloneExistence( const Vecteur &trans ) const;

  /** @brief Renvoie si un clone periodique existe ou non
  @param ObsID numero d'obstacle periodique */
  bool periodicCloneExistence( const int& ObsID ) const;

  /** @brief Renvoie si un clone periodique existe ou non
  @param pClone addresse du clone */
  bool periodicCloneExistence( Particule const* pClone ) const;

  /** @brief Ajoute un clone periodique
  @param clone Le clone */
  void addPeriodicClone( Particule *clone );

  /** @brief Ajoute un obstacle periodique en contact avec la particule
  @param ObsID numero d'obstacle */
  void addPeriodicObstacleID( int const& ObsID );

  /** @brief Gestion des clones multi-periodiques
  @param particulesClonesPeriodiques clones periodiques
  @param ParticuleClassesReference classes de r�f�rence de particules
  @param LC grille de cellules */
  void createMultiPeriodicClones(
  	list<Particule*>* particulesClonesPeriodiques,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC );

  /** @brief Translate les clones et met � jour leur lien dans le LinkedCell
  @param translation vecteur de translation
  @param LC grille de cellules */
  void translateAndUpdatePeriodicClones( Vecteur const& translation,
  	LinkedCell* LC );

  /** @brief Mise a jour de la localisation geographique de la particule
  dans le LinkedCell en utilisant cellule_nm1. Si la methode est appelee apres
  le LinkUpdate du LinkedCell, cellule_nm1 = cellule courante a laquelle
  appartient la particule */
  void updateGeoLocalisation();

  /** @brief Supprime le clone de la liste
  @param clone particule */
  void erasePeriodicClone( Particule *clone );

  /** @brief Supprime l'obstacle periodique de la liste
  @param ObsID numero d'obstacle */
  void erasePeriodicObstacleID( int const& ObsID );

  /** @brief Creer la structure AddedMassInfos */
  void createAddedMassInfos();

  /** @brief Renvoie un vecteur orientation de la particule */
  Vecteur vecteurOrientation() const;

  /** Ajoute un nombre de contacts au nombre de contacts de la particule;
  Utilisation: ajoute a la particule de reference periodique les contacts
  de son clone periodique */
  virtual void addToCoordinationNumber( int const& nc );

  /** @brief Initialise le torseur des efforts sur la particule
  @return La somme des forces de contact
  @param withWeight avec ou sans le poids de la particule */
  void InitializeForce( bool const& withWeight );

  /** @brief Initialise force at the contact point for the post-processing
  of the stress tensor */
  virtual void InitializeForceAtContactPoint( );

  /** @brief Calcul du poids de la particule */
  virtual void computeWeight();

  /** @brief Alloue la structure pour le couplage DEM-CFD */
  void allocateDEMCFD_FluidInfos();

  /** @brief Contact entre particule composite et une particule ou un obstacle
  Si le contact est realise, on ajoute le contact au deux composants.
  @exception ErreurContact Si un choc entre croute est detecte
  @param voisin Le composant (monolithe) a etudier
  @param dt Increment de temps
  @param temps temps physique
  @param LC grille de cellules
  @param listContact list containing contacts information */
  virtual void SearchContact( Composant* voisin, double dt,
      double const& temps, LinkedCell *LC,
      list<ContactInfos*> &listContact );
  //@}


  /** @name Accessors */
  //@{
  /** @brief Get vitesse de translation du fluide */
  Vecteur const* get_VitesseTranslation_fluide() const;

  /** @brief Get vorticite du fluide */
  Vecteur const* get_Vorticity_fluide() const;

  /** @brief Get gradient de pression du fluide */
  Vecteur const* get_GradientPression_fluide() const;

  /** @brief Acces a la cinematique de la particule
  @return La cinematique */
  const CineParticule* getCinematique() const;

  /** @brief Qu'elle est l'activite actuelle de la particule ?
  @return L'etat d'activite */
  ParticuleActivity getActivity() const { return m_activity; }

  /** @brief Type de Particule
  @return "P" pour les Particule, "PP" pour les Periodique */
  string getType() const { return m_type; }

  /** @brief Tag de la particule */
  int getTag() const { return m_tag; }

  /** @brief Mobilite de la particule */
  bool getMobilite() const { return is_mobile; }

  /** @brief Temperature evolution de la particule */
  bool getTempEvolution() const { return m_temp_evolution; }

  /** @brief Difference de vitesse de translation au temps pr�c�dent.*/
  Vecteur getTranslationalVelocityDifferencePreviousTime() const;

  /** @brief Difference de vitesse de rotation au temps pr�c�dent. */
  Vecteur getRotationalVelocityDifferencePreviousTime() const;

  /** @brief Get particle velocity at the previous time step. */
  Vecteur getTranslationalSlipVelocityPreviousTime() const;

  /** @brief Inertie de la particule.
  @return Le vecteur inertie */
  virtual const double* getInertie() const;

  /** @brief Inertie inverse de la particule.
  @return Le vecteur inertie inverse */
  const double* getInertieInverse() const;

  /** @brief Masse volumique de la particule
  @return La masse volumique */
  Scalar getMasseVolumique() const;

  /** @brief Rayon de la sphere de meme volume que le composant.
  @return Le rayon. */
  Scalar getRayonSphereEquivalente() const;

  /** @brief Energie de la particule
  @return L'energie de la particule */
  double getEnergie() const;

  /** @brief Vitesse du composant en un point (appartenant au composant a priori
  mais cette condition n'est pas v�rifi�e)
  @return La vitesse
  @param pt le point de definir la partie rotation avec la bras de levier par
  rapport au centre de gravite du composant */
  Vecteur getVitesse( const Point &pt ) const;

  //@{
  /** @brief Particle stress tensor */
  virtual vector<Scalar> const* getStressTensor( );

  /** @brief Internal moment divided by the volume of the system */
  virtual vector<Scalar> const* getInternalMoment( );

   /** @brief Rotation. */
  Quaternion const* getRotation() const;

  /** @brief Vitesse de rotation.
  @return La composante demandee. */
  virtual Vecteur const* getVitesseRotation() const;

  /** @brief Vitesse de translation. */
  virtual Vecteur const* getVitesseTranslation() const;

  /** @brief Force. */
  Vecteur const* getForce() const;

  /** @brief Get sum of contact force for Post Processing purposes. */
  Vecteur const* getForceContactPP() const;

  /** @brief Get sum of contact force for Post Processing purposes. */
  Vecteur const* getForceContactPP_instantaneous() const;

  /** @brief Get sum of lubri force for Post Processing purposes. */
  Vecteur const* getForceLubriPP() const;

  /** @brief Cellule dans laquelle la particule se situait au
  pas de temps pr�c�dent */
  Cellule* getCelluleNm1() const { return m_cellule_nm1; }

  /** @brief Classe de la particule */
  int getParticuleClasse() const { return m_ParticuleClasse; }

  /** @brief Renvoie le clone periodique dont la reference est en contact avec
  l'obstacle et le supprime de la liste
  @param obs obstacle en contact avec la reference */
  Particule* getPeriodicCloneAndErase( const ObstaclePeriodique* obs );

  /** @brief Renvoie le nombre de clones periodiques */
  size_t getNombreClonesPeriodiques() const
  	{ return max( m_periodicClones.size(), m_periodicObstaclesID.size() ); }

  /** @brief Localisation geographique */
  MPIGeoLocalisation getGeoLocalisation() const;

  /** @brief Renvoie la liste des numeros d'obstacle periodique
  avec lesquels la particule est en contact */
  list<int> const* getPeriodicObstaclesID() const
  	{ return &m_periodicObstaclesID; }

  /** @brief Renvoie la liste des clones periodiques */
  list<Particule*> const* getPeriodicClones() const
  	{ return &m_periodicClones; }

  /** @brief La particule poss�de un clone dans une direction donnee ?
  @param direction la direction de periodicite
  @param LC grille de cellules */
  bool hasCloneInDirection( Vecteur const* direction,
	LinkedCell const* LC ) const;

  /** @brief Renvoie le nombre de contacts de la particule */
  virtual int getCoordinationNumber() const;

  /** @brief Fluid volume fraction */
  double get_DEMCFD_volumeFraction() const;

  /**
    @brief Get DEMCFD surounding fluid temperature
    @return Returns DEMCFD surounding fluid temperature
  */
  double const* get_DEMCFD_fluidTemperature() const;

  /**
    @brief Get body temperature
    @return Returns solid body temperature
  */
  double const* get_solidTemperature() const;

  /**
    @brief Get fluid-solid heat flux store in m_infos
    @return Returns fluid-solid heat flux
  */
  double const* get_fluidSolidHeatFlux() const;

  /**
    @brief Get fluid-solid Nusselt number store in m_infos
    @return Returns fluid-solid heat flux
  */
  double const* get_solidNusselt() const;

  /** @brief Vitesse de translation du fluide */
  Vecteur const* getVitesseTranslation_fluide() const;

  /** @brief Vorticite du fluide */
  Vecteur const* getVorticity_fluide() const;

  /** @brief Gradient de pression du fluide */
  Vecteur const* getGradientPression_fluide() const;

  /**
    @brief Get particle hydroForce
    @param Fx x component of hydro force
    @param Fy y component of hydro force
    @param Fz z component of hydro force
  */
  void getParticleHydroForce( double& Fx, double& Fy, double& Fz );

  /** @brief Get particle hydroForce BIS */
  virtual Vecteur const* getParticleHydroForce() const;
  //@}

/** @brief Get particle SlipVel BIS */
  virtual Vecteur const* getParticleSlipVel() const;

  /** @brief Renvoie le nombre de sommets ou un code �quivalent */
  virtual int getNbCorners() const;
  //@}

  /** @brief get rnd variables for the stochastic drag law */
  Vecteur const* get_rnd() const;
  //@}

  /** @brief get rnd variable for the stochastic Nusselt law */
  Vecteur const* get_rnd_Nu() const;
  //@}

  /** @name Processus de restauration */
  //@{
  /** @brief Sauvegarde de l'etat */
  void saveState();

  /** @brief Cree et renvoie l'etat */
  pair<ConfigurationMemento*,CineParticuleMemento*> createState();

  /** @brief Restauration de l'etat */
  void restaureState();

  /** @brief Restauration de l'etat
  @param Pmemento_ configuration
  @param Cmemento_ cinematique
  @param contactReferenceObstacle_ pointeur sur l'obstacle periodique de la
  	particule de reference
  @param reference_ particule de reference periodique */
  virtual void restaureState( ConfigurationMemento const* Pmemento_,
  	CineParticuleMemento const* Cmemento_,
	ObstaclePeriodique const* contactReferenceObstacle_ = NULL,
	Particule* reference_ = NULL );
  //@}


  /**@name Methods I/O */
  //@{
  /** @brief Lecture de la particule pour Reload. Utilise pour les particules de
  classe de reference dans l'entete du fichier de reload et pour l'ancien format
  de reload
  @param fileSave Flux de lecture
  @param ParticuleClassesReference particules de r�f�rence
  	pour chaque classe de particules */
  virtual void read( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference = NULL );

  /** @brief Lecture de la particule pour Reload. Utilise pour les particules
  dans la simulation (actives ou en attentes) pour le format de reload 2014
  @param fileSave Flux de lecture
  @param ParticuleClassesReference particules de r�f�rence
  	pour chaque classe de particules */
  virtual void read2014( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference );

  /** @brief Lecture de la particule pour Reload. Utilise pour les particules
  dans la simulation (actives ou en attentes) pour le format de reload 2014
  en binaire
  @param fileSave Flux de lecture
  @param ParticuleClassesReference particules de r�f�rence
  	pour chaque classe de particules */
  virtual void read2014_binary( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference );

  /** @brief Sauvegarde de la particule pour Reload
  @param fileSave Flux de sauvegarde
  @param composant Composant de ref : utile dans le cas de particule composite
  et NULL dans le cans d'une particule convexe */
  virtual void write( ostream &fileSave, Composant const* composant = NULL )
      const;

  /** @brief Sauvegarde de la particule pour Reload
  @param fileSave Flux de sauvegarde */
  void write2014( ostream &fileSave ) const;

  /** @brief Sauvegarde en binaire de la particule pour Reload
  @param fileSave Flux de sauvegarde */
  void write2014_binary( ostream &fileSave );

  /** @brief Ecriture de l'identite du monolithe */
  void writeIdentity( ostream &file ) const;

  /** @brief Ecriture de l'information de position pour le Fluide */
  virtual void writePositionInFluid( ostream &fluid );

  /** @brief Operateur <<
  @param f flux de sortie
  @param P particule */
  friend ostream& operator << ( ostream &f, const Particule &P );
  //@}


  /** @name Methods Virtual pour la Periodicite */
  //@{
  /** @brief Mise � jour de la vitesse & position de la particule periodique */
  virtual void updateVitessePositionPeriodique() {}

  /** @brief Renvoie l'obstacle en contact avec la reference */
  virtual ObstaclePeriodique const* getObstacle() const { return NULL; }

  /** @brief Renvoie le nombre de periodes de la particule periodique: 1 =
  uni-periodicite, 2 = bi-periodicite, 3 = tri-periodicite */
  virtual int getNbPeriodes() const { return 0; }

  /** @brief Renvoie le vecteur de periodicite */
  virtual Vecteur const* getVecteurPeriodique() const { return NULL; }

  /** @brief Renvoie la particule de reference */
  virtual Particule* getPeriodicReference() const { return NULL; }

  /** @brief Ajout des forces & moments de contact � la particule maitre du
  clone
  @param temps Temps de simulation.
  @param dt Pas de temps */
  virtual void AddForcesFromPeriodicCloneToParticule( Scalar temps,
  	Scalar dt, bool ContactforceOutput ) const {}

  /** @brief Definit l'etat de restauration du clone
  @param done_ restaure ou non */
  virtual void setRestaureState( bool const& done_ ) {}

  /** @brief Recupere l'etat de restauration du clone */
  virtual bool getRestaureState() { return true; }
  //@}

  /** @brief Ecrit le convexe pour post-processing avec Paraview
  @param connectivity connectivite des polygones Paraview
  @param offsets decalage dans la connectivite
  @param cellstype type de polygones Paraview
  @param firstpoint_globalnumber numero global du 1er point
  @param last_offset dernier offset utilise pour le convexe precedent */
  virtual void write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const ;

  /** @brief Ecrit les points du convexe pour post-processing avec Paraview
  @param f Flux de sortie
  @param transform transformation courante
  @param translation translation du centre de gravite */
  virtual void write_polygonsPts_PARAVIEW( ostream &f,
  	const Transform &transform,
  	Vecteur const* translation = NULL ) const;

  /** @brief Nombre de points pour post-processing avec Paraview */
  virtual int numberOfPoints_PARAVIEW() const;

  /** @brief Nombre de polygones elementaires pour post-processing avec
  Paraview */
  virtual int numberOfCells_PARAVIEW() const;

  /** @brief Constructeur d'un clone par copie */
  virtual Particule* createCloneCopy() const ;

  /** @brief Renvoie les points du convexe pour post-processing avec Paraview
  @param translation translation du centre de gravite */
  virtual list<Point> get_polygonsPts_PARAVIEW( Vecteur const*
	translation = NULL ) const;

  /** @brief Ecrit les points du convexe pour post-processing avec Paraview
  @param f Flux de sortie
  @param translation translation du centre de gravite */
  virtual void write_polygonsPts_PARAVIEW( ostream &f,
  	Vecteur const* translation = NULL ) const;

  /** @brief Initialize a faux le boolean correspondant au calcul de la
  transformation avec scaling par l'epaisseur de croute */
  virtual void initializeVdWtransform_to_notComputed();

  //@}

  /** @name Methods CompParticule */
  //@{

  /** @brief Renvoie les limites du domaine pour la discretisation
  de la particule composite afin d'approximer les elements d'inertie
  v[0] = x_max, v[1] = x_min
  v[2] = y_max, v[3] = y_min
  v[4] = z_max, v[5] = z_min */
  virtual vector<Scalar> getCompFeaturesBox() ;

  /** @brief Acces aux particules elementaires d'une particule composite */
  virtual vector<ElementParticule*> getElementParticules() const;

  /** @brief Renvoie les positions initiales des particules elementaires a
  l'insertion */
  virtual vector<Vecteur> getInitialRelativePositions() const;

  /** @brief Renvoie les positions des particules elementaires
  apres mise a jour */
  virtual vector<Vecteur> getRelativePositions() const;

  /** @brief Renvoie les matrices de rotations intiales des particules
  elementaires */
  virtual vector<Matrix> getInitialMatrix() const;

  /** @brief Renvoie le nombre total de particules elementaires pour
  la particule composite */
  virtual size_t getNbreElemPart() const;

  /** @brief Renvoie le nom associe a la particule composite */
  virtual string getPartName() const;

  /** @brief Positionne les particules elementaires dans l'espace */
  virtual void setElementPosition();

  Scalar set_shrinking_mass(); /** set the current mass */
    //@}

protected:
  /**@name Parameters */
  //@{
  CineParticule *m_cinematique; /** Algorithme de resolution du deplacement et
  	de la rotation */
  Scalar m_masseVolumique; /**< Masse volumique */
  bool is_mobile; /**< If particle moves or not */
  bool m_temp_evolution; /**< Temperature evolution with time */
  static Scalar m_fluideMasseVolumique; /**< Masse volumique du fluide
  	(default 0.) */
  static Scalar m_fluidViscosity; /**< Viscosite du fluide
  static Scalar m_fluidViscosity; /** Viscosite du fluide
  	(default 0.) */
  static double m_fluidInitialTemperature; /** Fluid initial temperature
  	(default 0.) */
  static bool m_MassCorrection; /** Correction de la masse devant dv/dt ? */
  static bool m_explicitAddedMass; /** Traitement explicite ou non de la
  	masse ajout�e */
  Scalar m_inertie[6]; /** Matrice d'inertie I={I(1,1), I(1,2), I(1,3),
  	I(2,2), I(2,3), I(3,3)} */
  Scalar m_inertie_1[6]; /** Matrice d'inertie inverse */
  double m_energie; /** Energie de la Particule */
  ParticuleActivity m_activity; /** Activite de la Particule */
  string m_type; /** Type de la Particule */
  struct AddedMassInfos *m_addedMassInfos; /** vitesses et differences de
  	vitesse au temps precedent pour calcul explicite de la masse ajoutee,
	utilise en couplage avec un fluide quand les particules sont de meme
	densite ou plus legeres que le fluide */
  int m_tag; /** Tag de la cellule: 0=interieur, 1=dans la zone de recouvrement
  	d'un autre processeur et geometriquement sur ce processeur,
  	2=dans la zone de recouvrement de ce processeur et geometriquement sur
  	un autre processeur */
  MPIGeoLocalisation m_GeoLoc; /** localisation geographique de la particule
  	dans le LinkedCell, utile pour les particules taggees 1 uniquement */
  Cellule* m_cellule_nm1; /** cellule dans laquelle la particule se situait au
  	pas de temps pr�c�dent */
  int m_ParticuleClasse; /** classe de particules � laquelle la particule
  	appartient */
  list<Particule*> m_periodicClones; /** clones periodiques de la
  	particule dans le cas d'une simulation periodique (utilises
	en sequentiel) */
  list<int> m_periodicObstaclesID; /** liste des numeros d'obstacle periodique
  	avec lesquels la particule est en contact
	dans le cas d'une simulation periodique */
  int m_coordination_number; /** nombre de voisins avec qui la particule est
  	en contact */
  Vecteur m_weight; /** poids de la particule */
  struct DEMCFD_FluidInfos *m_fluidInfos; /** infos provenant du fluide pour le
  	couplage DEM-CFD */
  Scalar m_solidTemperature; /** Solid body temperature */
  double m_solidNusselt; /** Solid Nusselt number */

  Vecteur m_g_rnd;/** random variable vector for stochastic drag laws*/
  Vecteur m_g_rnd_Nu;/** random variable vector for stochastic Nusselt laws*/

  Vecteur m_hydroForce_PP;/** hydroforce for post-processing*/
  //@}

  /** Elements of the internal moment of the particle*/
  Scalar m_mxx, m_mxy, m_mxz;
  Scalar m_myx, m_myy, m_myz;
  Scalar m_mzx, m_mzy, m_mzz;
  vector<Scalar> m_InternalMoment; /** Internal moments for the stress tensor */

  //@}

private:

  /**@name Constructors & Destructor */
  //@{
  /** @brief Constructeur par defaut */
  Particule();
  //@}
};

#endif
