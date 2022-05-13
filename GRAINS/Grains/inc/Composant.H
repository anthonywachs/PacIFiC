#ifndef _Composant
#define _Composant

#include "Basic.H"
#include "FormeVdW.H"
#include "Point.H"
#include "PointContact.hh"
#include "Torseur.H"
#include <list>
#include <string>
#include <map>
using namespace std;

#define NB_STEPS_REMEMBER_ENDED_CONTACT 1


class BBox;
class ErreurContact;
class ConfigurationMemento;
class LinkedCell;

/** @brief Gestion des contacts dans le cas de particule composite
@author D. RAKOTONIRINA - Dec. 2014 - Creation */
// ============================================================================
struct ContactInfos
{
  PointContact ContactPoint; /** Point de contact */
  Composant* p0; /** composant0 a etudier */
  Composant* p1; /** composant1 a etudier */
};

/** @brief Info au point de contact pour le post-processing de la contrainte
@author D. RAKOTONIRINA - Fev. 2017 - Creation */
// ============================================================================

/** @brief Gestion des composants physiques des applications
    en mecanique granulaire.

    Ces composants sont divises en fonction du type de cinematique et
    d'interaction qui leurs sont appliquees.

    @author G.FERRER - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Composant
{
public:
  /** @name Constructors */
  //@{
  /** @brief Destructeur */
  virtual ~Composant();
  //@}


  /** @name Methods Virtual */
  //@{
  /** @brief Vitesse du composant en un point (appartenant au composant a priori
  mais cette condition n'est pas v�rifi�e)
  @return La vitesse
  @param pt le point de definir la partie rotation avec la bras de levier par
  rapport au centre de gravite du composant */
  virtual Vecteur getVitesse( const Point &pt ) const = 0;

  /** @brief Vitesse de rotation.
  @return La composante demandee. */
  virtual Vecteur const* getVitesseRotation() const = 0;

  /** @brief Vitesse de translation. */
  virtual Vecteur const* getVitesseTranslation() const = 0;

  /** @brief Renvoie le nom associe a la particule composite */
  virtual string getPartName() const;

  /** @brief Les deux composants sont ils en contact ?
  Si les deux composants sont identiques (==), le contact n'existe pas.
  Limitation: le composant voisin ne doit PAS etre de type deriv� CompObstacle.
  @return TRUE Si un contact existe
  @param voisin Le deuxieme composant etudie. */
  virtual bool isContact( const Composant* voisin ) const;

  /** @brief Les deux composants sont ils en contact ? Variante avec les rayons
  de VdW. Si les deux composants sont identiques (==), le contact n'existe pas.
  Limitation: le composant voisin ne doit PAS etre de type deriv� CompObstacle.
  @return TRUE Si un contact existe
  @param voisin Le deuxieme composant etudie. */
  virtual bool isContactVdW( const Composant* voisin ) const;

  /** @brief Les deux composants sont ils proches ? La proximite est liee a
  l'intersection des BBox des deux formes. Si les deux composants sont
  identiques (==), le contact n'existe pas.
  Limitation: le composant voisin ne doit PAS etre de type deriv� CompObstacle.
  @return TRUE Si un contact existe
  @param voisin Le deuxieme composant etudie. */
  virtual bool isProche( const Composant* voisin ) const;

  /** @brief Les deux composants sont ils en contact ? Variante avec les rayons
  de VdW. La proximite est liee a l'intersection des BBox des deux formes.
  Si les deux composants sont identiques (==), le contact n'existe pas.
  Limitation: le composant voisin ne doit PAS etre de type deriv� CompObstacle.
  @return TRUE Si un contact existe
  @param voisin Le deuxieme composant etudie. */
  virtual bool isProcheVdW( const Composant* voisin ) const;

  /** @brief Rotation de la particule.
  @param rotation Quaternion decrivant la rotation a effectuer. */
  virtual void Rotate( const Quaternion &rotation );

  /** Translation de la particule
  @param translation Le vecteur de translation. */
  virtual void Translate( const Vecteur &translation );

  /** @brief Initialise le torseur des efforts sur le composant
  @param withWeight avec ou sans le poids du composant */
  virtual void InitializeForce( bool const& withWeight );

  /** @brief Initialise force at the contact point for the post-processing
  of the stress tensor */
  virtual void InitializeForceAtContactPoint();

  /** @brief Initialize PP forces */
  void InitializePostProcessingForce( );

  /** @brief Initialize PP forces
  @param nb number of analyzed solid time steps between two saving output */
  void IntegrateContactPostProcessingForce(Scalar nb);

  /** @brief Determine s'il y a intersection entre une Forme et
  un Composant
  @return Vrai si Forme et Composant s'intersectent
  @param b Une Forme B */
  virtual bool intersectWithForme( const Forme &b ) ;

  /** @brief Update contact map */
  virtual void updateContactMap();

  /** @brief Does the contact exist in the map? If so, return true and make
  kdelta, prev_normal and cumulSpringTorque point to the memorized info. Otherwise, return false and set those pointers to NULL.
  @param id key in the map
  @param kdelta pointer to the memory of the vector kt * delta_t
  @param prev_normal pointer to the previous vector normal to the contact plane
  @param cumulSpringTorque pointer to the memory of the spring-like component of the friction torque */
  virtual bool getContactMemory( std::tuple<int,int,int> const& id,
    Vecteur* &kdelta, Vecteur* &prev_normal, Vecteur* &cumulSpringTorque );

  /** @brief Add new contact in the map
  @param id key in the map
  @param kdelta pointer to the memory of the vector kt * delta_t
  @param prev_normal pointer to the previous vector normal to the contact plane
  @param cumulSpringTorque pointer to the memory of the spring-like component of the friction torque */
  virtual void addNewContactInMap( std::tuple<int,int,int> const& id,
    Vecteur const& kdelta, Vecteur const& prev_normal,
    Vecteur const& cumulSpringTorque );

  /** @brief Store memory of the contact with component id: increase cumulative tangential displacement and cumulative spring torque, remember contact normal.
  @param id key in the map
  @param kdelta pointer to the memory of the vector kt * delta_t
  @param prev_normal pointer to the previous vector normal to the contact plane
  @param cumulSpringTorque pointer to the memory of the spring-like component of the friction torque */
  virtual void addDeplContactInMap( std::tuple<int,int,int> const& id,
    Vecteur const& kdelta, Vecteur const& prev_normal,
    Vecteur const& cumulSpringTorque );

  /** @brief copyHistoryContacts write the contact map information in an array of doubles
  @param destination the array of double where the contact map should be stored
  @param start index the index of destination where the copy should start */
  virtual void copyHistoryContacts( double* &destination, int start_index );

  /** @brief copyContactInMap add a single contact info to the contact map
  @param id key in the map
  @param isActive boolean: true if the contact is active, false otherwise
  @param kdelta pointer to the memory of the vector kt * delta_t
  @param prev_normal pointer to the previous vector normal to the contact plane
  @param cumulSpringTorque pointer to the memory of the spring-like component of the friction torque */
  virtual void copyContactInMap( std::tuple<int,int,int> const& id,
    bool const& isActive, Vecteur const& kdelta, Vecteur const& prev_normal,
    Vecteur const& cumulSpringTorque );

  /** @brief getContactMapSize Return the number of contacts in the contact map
  */
  virtual int getContactMapSize();

  /** @brief printActiveNeighbors displays the active neighbours in the format "my_elementary_id/neighbour_id/neightbout_elementary_id ; ...". Useful for debugging only.
  @param id id of this component
  */
  virtual void printActiveNeighbors(int const& id );
  //@}


  /** @name Methods Set */
  //@{

  /**
    @brief Set "composant" identifier (@ref m_id).
    @param id_ l'ID du composant.
  */
  void setID( const int &id_ ) { m_id = id_;}

  /** @brief Positionne un composant dans l'espace.
  @param pos Description de la matrice de base
  (12 coefficients : 3+1+3+1+3+1) puis la position du point,
  soit un vecteur de 16 doubles */
  virtual void setPosition( const Scalar *pos );

  /** @brief Positionne une particule a partir d'un point.
  @param centre position du centre de gravit� de la particule */
  virtual void setPosition( const Point &centre );

  /** @brief Definit le type de materiau
  @param mat nom du materiau */
  void setMateriau( const string &mat );

  //@{
  /** @brief Force at the contact point for the stress tensor
  @param force la nouvelle force
  @param point point d'action de la force */
  virtual void computeInternalMoments( const Point &point, const Vecteur &force );

/*
    @brief Set solid-body heat capacity
    @param Solid-body heat capacity
  static void set_heatCapacity( const Scalar &heatCapacity_ );
  */

  /**
    @brief Definit la distance ou la force d'interaction est maximale
    @param id_particle Colliding particle identifier
    @param distance Distance for which cohesive force is maximum
    @param kn_elast Normal stiffness coefficient
    @param kt_elast Tangencial stiffness coefficient
    @param overlap Distance between particles surfaces
  */
  void initializeCohesiveProperties( int id_particle, Scalar distance,
      Scalar kn_elast, Scalar kt_elast , Scalar overlap );

  /**
    @brief Definit la distance ou la force d'interaction est maximale
    @param id_particle Colliding particle identifier
    @param distance Distance for which cohesive force is maximum
    @param kn_elast Normal stiffness coefficient
    @param kt_elast Tangencial stiffness coefficient
  */
  void set_CohesiveProperties( int id_particle, Scalar distance,
      Scalar kn_elast, Scalar kt_elast );

  /** @brief Initialize all contact map entries to false */
  virtual void setContactMapToFalse();
  //@}


  // ===========================================================================
  /** @name Accessors */
  //@{
  /**
    @brief Determination de la "Bounding Box" du composant.
    @return La boite contenant le composant.
  */
  BBox Boite() const;

  /**
    @brief get "composant" identifier.
    @return l'ID du composant (@ref m_id).
  */
  int getID() const;

  /** @brief Acces a la forme du composant.
  Nota : cf. Particule::Interaction(Obstacle& ...)
  @return Le pointeur d'acces a la forme du composant. */
  const FormeVdW* getForme() const;

  /** @brief Acces a la forme du composant.
  Nota : cf. Particule::Interaction(Obstacle& ...)
  @return Le pointeur d'acces a la forme du composant. */
  FormeVdW* getForme();

  /** @brief Masse du Composant.
  @return Valeur de masse */
  double getMasse() const;

  /** Point du centre de gravite du composant
      @return Le point demande */
  Point const* getPosition() const;

  /** @brief Position de la particule dans le pointeur pos.*/
  void getPosition( Scalar *pos ) const;

  /** @brief Rayon externe du composant.
  @return Le rayon. */
  virtual Scalar getRayon() const;

  /** @brief return shrinking choice of composant.
  @return shrinking. */
  virtual int getShrinkingMode();

  /** @brief Rayon de la sphere de meme volume que le composant.
  @return Le rayon. */
  virtual Scalar getRayonSphereEquivalente() const;

  /** @brief Rayon d'interaction du composant.
  @return Le rayon. */
  Scalar getRayonInteraction() const;


  /** @brief shrinking radius of composant.
  @return shrinking_radius. */
  void set_shrinking_radius(Scalar CurrentRadius);


  /** @brief Volume du composant.
  @return Le volume */
  virtual Scalar getVolume() const;

  /** @brief Nom du materiau
  @return Le materiau defini */
  const string& materiau() const;

  /** @brief Renvoie la force appliquee */
  Vecteur const* getForce() const;

  /** @brief Renvoie le moment applique */
  Vecteur const* getMoment() const;

  /** @brief Renvoie le torseur applique */
  virtual Torseur const* getTorseur();

  /** @brief Vitesse de translation du fluide */
  virtual Vecteur const* getVitesseTranslation_fluide() const;

  //@{
  /** @brief Particle stress tensor */
  virtual vector<Scalar> const* getStressTensor( );

  /** @brief Internal moment divided by the volume of the system */
  virtual vector<Scalar> const* getInternalMoment( );

  /**
    @brief Get distance for which cohesive force is maximum
    @param id_particle Particle identifier
  */
  Scalar get_FmaxDist( int id_particle ) const;

  /**
    @brief Get normal stiffness coefficient
    @param id_particle Particle identifier
  */
  Scalar get_KnElast( int id_particle ) const;

  /**
    @brief Get tangential stiffness coefficient
    @param id_particle Particle identifier
  */
  Scalar get_KtElast( int id_particle ) const;

  /**
    @brief Get initial overlap
    @param id_particle Particle identifier
  */
  Scalar get_InitialOverlap( int id_particle ) const;

  /** @brief Gestion des contacts en simmulation PERIODIQUE SEQ:
  Reference du composant lors du calcul des forces de contact
  @return this si contact avec composant autre que CompParticule
  sinon m_masterComposite si contact avec CompParticule */
  virtual Composant* ReferenceComposant();

  /** @brief Contact entre particule composite et une particule ou un obstacle
  Si le contact est realise, on ajoute le contact au deux composants.
  @exception ErreurContact Si un choc entre croute est detecte
  @param voisin Le composant (monolithe) a etudier
  @param dt Increment de temps
  @param temps temps physique
  @param LC grille de cellules
  @param listContact list of contact information */
  virtual void SearchContact( Composant* voisin, double dt,
      double const& temps, LinkedCell *LC,
      list<ContactInfos*> &listContact ){ };

  /** @brief Nombre de points pour post-processing avec Paraview */
  virtual int numberOfPoints_PARAVIEW() const = 0;

  /** @brief Nombre de polygones elementaires pour post-processing avec
  Paraview */
  virtual int numberOfCells_PARAVIEW() const = 0;

  /** @brief Renvoie les points du convexe pour post-processing avec Paraview
  @param transform transformation courante
  @param translation translation du centre de gravite */
  virtual list<Point> get_polygonsPts_PARAVIEW(
  	Vecteur const* translation = NULL ) const = 0;

  /** @brief Ecrit les points du convexe pour post-processing avec Paraview
  @param f Flux de sortie
  @param translation translation du centre de gravite */
  virtual void write_polygonsPts_PARAVIEW( ostream &f,
  	Vecteur const* translation = NULL ) const = 0;

  /** @brief Ecrit le convexe pour post-processing avec Paraview
  @param connectivity connectivite des polygones Paraview
  @param offsets decalage dans la connectivite
  @param cellstype type de polygones Paraview
  @param firstpoint_globalnumber numero global du 1er point
  @param last_offset dernier offset utilise pour le convexe precedent */
  virtual void write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const = 0;
  //@}


  /** @name Methods */
  //@{
  /** @brief Ajout d'une force au torseur des efforts
  @param force la nouvelle force
  @param point point d'action de la force */
  virtual void addForce( const Point &point, const Vecteur &force );

  /** @brief Add Contact force for postprocessing purposes
  @param force the new force */
  virtual void addContactForcePP( const Vecteur &force );

  /** @brief Add Contact force for postprocessing purposes
  @param force the new force */
  virtual void addContactForcePP_instantaneous( const Vecteur &force );

  /** @brief Add lubrication force for postprocessing purposes
  @param force the new force */
  virtual void addLubriForcePP( const Vecteur &force );

  /** @brief Ajout d'une force s'exercant au centre de gravite au torseur des
  efforts
  @param force la nouvelle force */
  virtual void addBodyForce( const Vecteur &force );

  /** @brief Ajouter un moment
  @param moment La nouvelle force */
  virtual void addMoment( const Vecteur &moment );

  /**
    @brief Add heat flux to solid-body
    @param heatFlux_ heat flux adress
  */
  virtual void add_heatFlux( const Scalar &heatFlux_ );

  /**
    @brief Position de la particule dans le pointeur pos en d�butant � la
  position i
  @param pos vecteur de valeurs
    @param i position de d�part
  */
  void copyPosition( double *pos, int i ) const;

  /** @brief Copie de la transformation dans le vecteur vit en d�butant �
  la position i
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  virtual void copyTransform( double *vit, int i ) const;

  /**
    @brief A FAIRE
  */
  virtual void copy_VectIdParticle( double* id_part, int i ) const;

  /**
    @brief A FAIRE
  */
  virtual void copy_VectFmaxDist( double* max_dist, int i ) const;

  /**
    @brief A FAIRE
  */
  virtual void copy_VectKnElast( double* kn_elast, int i ) const;

  /**
    @brief A FAIRE
  */
  virtual void copy_VectKtElast( double* kt_elast, int i ) const;

  /**
    @brief A FAIRE
  */
  virtual void copy_VectInitialOverlap ( double* initialOverlap, int i ) const;

  /**
    @brief A FAIRE
  */
  virtual void set_VectIdParticle( double id_part, int i );

  /**
    @brief A FAIRE
  */
  virtual void set_VectFmaxDist( double max_dist, int i );

  /**
    @brief A FAIRE
  */
  virtual void set_VectKnElast( double kn_elast, int i );

  /**
    @brief A FAIRE
  */
  virtual void set_VectKtElast( double kt_elast, int i );

  /**
    @brief A FAIRE
  */
  virtual void set_VectInitialOverlap( double initialOverlap, int i );

  /** @brief Copie de la transformation dans le vecteur vit en d�butant �
  la position i, avec une translation du centre de gravite (utile pour les
  particules periodiques en parallele)
  @param vit vecteur de copie
  @param i position dans le vecteur vit
  @param vec vecteur de translation */
  virtual void copyTransform( double *vit, int i, Vecteur const& vec ) const;

  /** Le composant est il dans la Boite ?
  @return Vrai si le centre est dans la boite */
  bool isIn( const BBox &boite ) const;

  /** @brief Applique la transformation trot au composant
  @param trot transformation � appliquer */
  virtual void composePosition( const Transform &trot );

  /** @brief Contact entre composant et composant.
  Si le contact est realise, on ajoute le contact au deux composants.
  @exception ErreurContact Si un choc entre croute est detecte
  @param voisin Le composant (monolithe) a etudier
  @param dt Increment de temps
  @param temps temps physique
  @param LC grille de cellules */
  virtual void InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC ) throw (ErreurContact) = 0;

  /** @brief Contact entre l'obstacle et une particule pour le post processing
  Si le contact est realise, on ajoute le contact au deux composants.
  @param voisin Le composant (monolithe) a etudier
  @param listOfContacts liste de contacts */
  virtual void InterActionPostProcessing( Composant* voisin, Scalar dt,
	list<struct PointForcePostProcessing>* listOfContacts )
	throw (ErreurContact) = 0;

  /** Ajoute un nombre de contacts au nombre de contacts de la particule;
  Utilisation: ajoute a la particule de reference periodique les contacts
  de son clone periodique */
  virtual void addToCoordinationNumber( int const& nc );

  /** @brief Le composant est il une particule composite ?
  @return Vrai si la particule est un composite */
  virtual bool isCompParticule() const { return ( false ); }

  /** @brief Est ce que le composant est un obstacle periodique ? */
  virtual bool isObsPeriodique() const { return ( false ); }

  /** Le composant est il un compObstacle ?
  @return Vrai si le composant est un compObstacle */
  virtual bool isCompObstacle() const {return ( false ); }

  /** @brief Est ce que le composant est un obstacle ? (utilise le fait que les
  obstacles ont une masse nulle par convention) */
  bool isObstacle() const {return ( m_masse == 0. ); }

  /** @brief Est ce que le composant est une particule elementaire ?
  ( Voir CompParticule ) */
  virtual bool isElementParticule() const { return ( false ); }

  /** @brief Update the ids of the contact map: in the case of a reload with insersion, the obstacle's ids are reset. This function keeps track of that change.
  @param prev_id previous id that should be updated
  @param new_id updated id
  */
  void updateContactMapId( int prev_id, int new_id);

  //@}


  /**@name Operators */
  //@{
  /** @brief Comparaison entre deux composants.
  @return Vrai si identite des deux composants : meme adresse.
  @param composant Le composant a comparer. */
  bool operator == ( const Composant &composant ) const;

  /** @brief Comparaison entre deux composants.
  @return Vrai si non identite des deux composants : adresse differente.
  @param composant Le composant a comparer. */
  bool operator != ( const Composant &composant ) const;
  //@}


  /**@name Methods I/O */
  //@{
  /** @brief Sauvegarde du composant pour Reload
  @param fileSave Flux de sauvegarde
  @param composant solid object */
  virtual void write( ostream &fileSave,
	Composant const* composant = NULL ) const = 0;

  /**
    @brief Ecriture de l'information de position pour le Fluide
    @param fluid Stream containing information to be transfered to fluid
    @remark Called by Composant::writePositionInFluid
    @par Content
      Calls Forme::writePositionInFluid
  */
  virtual void writePositionInFluid( ostream &fluid ) = 0;

  /** @brief Ecriture de l'identite du composant */
  virtual void writeIdentity( ostream &file ) const = 0;

  /** @brief Ecriture de l'objet dans sa configuration courante pour
  post-processing avec GMV
  @param fileOut Flux de sortie. */
  virtual void GMVoutput( ostream &fileOut ) const;
  //@}


  /** @name Methods Virtual pour les particules */
  //@{
  /** @brief Renvoie le numero de la particule de reference */
  virtual int getPeriodicReferenceID() const { return -100; }

  /** @brief Classe de la particule */
  virtual int getParticuleClasse() const { return -100; }

  /** @brief Renvoie le nombre de periodes de la particule periodique: 1 =
  uni-periodicite, 2 = bi-periodicite, 3 = tri-periodicite */
  virtual int getNbPeriodes() const { return 0; }

  /** @brief La particule poss�de un clone dans une direction donnee ?
  @param direction la direction de periodicite
  @param LC grille de cellules */
  virtual bool hasCloneInDirection( Vecteur const* direction,
	LinkedCell const* LC ) const { return false; }

  /** @brief Renvoie le nombre de clones periodiques */
  virtual size_t getNombreClonesPeriodiques() const { return 0; }

  /** @brief Renvoie le vecteur de periodicite */
  virtual Vecteur const* getVecteurPeriodique() const { return NULL; }
  //@}


  /** @name Methods Static */
  //@{
  /** @brief Definit le nombre de composants crees (generalement pour
  re-numerotation avant creation)
  @param nb_ nombre de composants crees */
  static void setNbComposantsCrees( const int &nb_ ) { m_nb = nb_; }

  /** @brief Renvoie le nombre de composants crees */
  static int getNbComposantsCrees() { return m_nb; }
  //@}


  /** @name Processus de restauration */
  //@{
  /** @brief Sauvegarde de l'etat de configuration */
  void saveConfigState();

  /** @brief Cree et renvoie l'etat de configuration */
  ConfigurationMemento* createConfigState();
  //@}

  /** @name Parameters */
  //@{
  //@}


protected:
  /** @name Constructors */
  //@{
  /** Constructeur par defaut
  @param autonumbering numerotation automatique ou non */
  Composant( const bool &autonumbering = true );

  /** Constructeur par copie.
  La liste des contacts de l'objet cree reste vide.
  @param copie Le composant a dupliquer. */
  Composant( const Composant &copie );
  //@}


  /** @name Methods I/O */
  //@{
  /** @brief Ecriture de l'information de position.
  @param position Flux recepteur
  @param composant solid object */
  void writePosition( ostream &position, Composant* composant = NULL ) const;

  /** @brief Ecriture de l'information statique.
  @param statique Flux recepteur
  @param composant solid object */
  virtual void writeStatique( ostream &statique,
	Composant const* composant = NULL ) const;

  /** @brief Write the contact map to file in plain 2014 format
  @param fileSave output file stream */
  void writeContactMemory_2014(ostream &fileSave ) const;

  /** @brief Write the contact map to file in binary format
  @param fileSave output file stream */
  void writeContactMemory_binary( ostream &fileOut );

  /** @brief Read the contact map to file in plain 2014 format
  @param fileSave input file stream */
  void readContactMap_2014( istream &fileSave);

  /** @brief Read the contact map to file in binary format
  @param fileSave input file stream */
  void readContactMap_binary( istream &fileSave );
  //@}


  /** @name Parameters */
  //@{
  int m_id; /** Identificateur */
  string m_nomMateriau; /** Nom du materiau associe au composant */

  Scalar m_masse; /** Masse du composant */
  FormeVdW *m_geoFormeVdw; /** forme geometrique avec epaisseur de croute type
    Van der Walls */
  vector<int> m_vectIdParticle; /** Particle Identifier vector */
  vector<Scalar> m_vectFmaxDist; /** Distance where the cohesive Force is
        maximum */
  vector<Scalar> m_vectKnElast; /** Normal elastic coef for cohesive force */
  vector<Scalar> m_vectKtElast; /** Tangential elastic coef for cohesive force */
  vector<Scalar> m_vectInitialOverlap; /** Initial overlap for cohesive force */
  Torseur m_somme; /** Somme des forces appliquees sur le composant */
  ConfigurationMemento *m_memento; /** Memento (stockage d'un etat de
  	restauration de la Particule) */
  Scalar m_sum_HeatFlux; /** Sum of all heat fluxes acting on a solid-body */
  Vecteur m_ForceContactPP; /** Sum of contact force acting on each particle
        for postprocessing purposes */
  Vecteur m_ForceContactPP_instantaneous; /** Sum of contact force acting on
  	each particle for postprocessing purposes */
  Vecteur m_ForceLubriPP; /** Sum of lubrication force acting on each particle
        for postprocessing purposes */
  map < std::tuple<int,int,int>, std::tuple<bool, Vecteur, Vecteur, Vecteur> > m_contactMap;
  /** List of active contacts with other components. It reads as follows:
  map<tuple<own elementary particle id, neighbour id, neighbour elementary
  particle id>, tuple<isContactActive, kt * cumulative tangential dispacement,
  previous normal vector, kr * cumulative rotational displacement> > */
  //@}

  /** @name Parameters Static */
  //@{
  static int m_nb; /** Nombre de composants crees */
//  static double m_heatCapacity; /** Solid-bodies heat capacity */
  //@}
};

#endif
