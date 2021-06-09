#ifndef _MonObstacle
#define _MonObstacle

#include "Obstacle.H"
#include <list>
#include <set>
using namespace std;

#include "ReaderXML.hh"

class Particule;
class App;
class Cellule;
class LinkedCell;
class ObstaclePeriodique;


/** @brief Gestion de l'obstacle de base : MON OBSTACLE.

    Un obstacle de base est un obstacle compose d'une Forme unique.
    Il est le composant de reference a partir duquel on peut construire des
    obstacles complexes. C'est a partir de lui qu'on evalue les contacts
    avec les particules.

    @author G.FERRER - Institut Francais du Petrole - 2002 - Creation */
// ============================================================================
class MonObstacle : public Obstacle
{
public:
  /** @name Constructeurs */
  //@{
  /** @brief Constructeur par defaut.
  @param s Nom eventuel de l'obstacle */
  MonObstacle( const string &s = "" );

  /** @brief Constructeur par copie d'un Composant. */
  MonObstacle( Composant &copie, char const* s = "obstacle" );

  /** @brief Constructeur par decodage d'un noeud XML
  @param root Le noeud "Obstacle" */
  MonObstacle( DOMNode *root );

  /**
    @brief Constructeur a partir d'une forme.
    @param geoFormeVdw Elementary form used to build the obstacle.
    @param name Obstacle name.
    @param materialName material name used.
    @param transferToFluid_ Boolean which says if the obstacle is seen by fluid.
    @remark Method called in CylindricalBoxObstacle constructor
  */
  MonObstacle( FormeVdW *geoFormeVdw, const string &name="",
      const string &materialName="", const bool &transferToFluid_=false );

  /** @brief Destructeur.*/
  ~MonObstacle();
  //@}


  /** @name Methods Get */
  //@{
  /** @brief L'obstacle a t'il le nom indique.
  @return this si le nom correspond a l'obstacle.
  @param nom_ nom de l'obstacle */
  const Obstacle* getNom( const string &nom_ ) const;

  /** @brief Liste des obstacles.
  @return On renvoie lui-meme. */
  list<MonObstacle*> getObstacles();

  /** @brief Liste des obstacles a transmettre au fluide
  @return Pour un obstacle composite on renvoie les obstacles de reference,
  pour un obstacle simple on renvoie lui-meme. */
//  list<MonObstacle*> getObstaclesToFluid();
  list<Obstacle*> getObstaclesToFluid();

  /** @brief Renvoie la liste des cellules auxquelles l'obstacle est lie */
  const list<Cellule*>* getInCells() const;

  /** @brief Renvoie la BBox de l'obstacle */
  const BBox* getObstacleBox() const;

  /** @brief Renvoie le pointeur sur l'obstacle si il est periodique, renvoie
  NULL sinon */
  virtual ObstaclePeriodique const* getObstaclePeriodic() const { return NULL; }

  /** @brief Renvoie la frequence de mise a jour du lien entre
  l'obstacle et le LinkedCell */
  int getObstacleLinkedCellUpdateFrequency() const
  	{ return m_LinkUpdate_frequency; }

  /** @brief Renvoie le type de l'obstacle */
  string getObstacleType();
  //@}


  /** @name Methods */
  //@{
  /** @brief Ajout d'un obstacle (unique ou composite) a l'arborescence
  Obstacle. En fait, sans action dans le cas d'un obstacle unique.
  Renvoie un message de warning.
  @param obstacle obstacle � ajouter */
  void append( Obstacle* obstacle );

  /** @brief Deplacement de l'obstacle.
  @return liste des obstacles deplaces
  @param temps Temps de l'increment precedent
  @param dt Increment de temps
  @param b_deplaceCine_Comp deplacement du composite dont l'obstacle fait
  partie (cinematique)
  @param b_deplaceF_Comp deplacement du composite dont l'obstacle fait
  partie (force) */
  list<MonObstacle*> Deplacer( Scalar temps,
	Scalar dt, const bool &b_deplaceCine_Comp,
	const bool &b_deplaceF_Comp );

  /** @brief Contact entre l'obstacle et une particule.
  Si le contact est realise, on ajoute le contact au deux composants.
  @param voisin Le composant (monolithe) a etudier
  @param dt Increment de temps
  @param temps temps physique
  @param LC grille de cellules
  @exception ErreurContact Si un choc entre croute est detecte */
  virtual void InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC ) throw (ErreurContact);

  /** @brief Contact entre l'obstacle et une particule pour le post processing
  Si le contact est realise, on ajoute le contact au deux composants.
  @param voisin Le composant (monolithe) a etudier
  @param listOfContacts liste de contacts */
  virtual void InterActionPostProcessing( Composant* voisin, Scalar dt,
	list<struct PointForcePostProcessing>* listOfContacts )
	throw (ErreurContact);

  /** @brief Construction d'un objet de type obstacle par lecture sur fichier.
  Cette methode est utilisee pour reconstruire un objet persistant.
  @param file Fichier de description */
  void reload( istream &file );

  /** @brief Rotation de la particule.
  @param rotation Quaternion decrivant la rotation a effectuer. */
  void Rotate( const Quaternion &rotation );

  /** @brief Suppression de l'obstacle dans la zone specifie. */
  void Suppression( const BBox &box );

  /** @brief Translation de la particule
  @param translation Le vecteur de translation. */
  void Translate( const Vecteur &translation );

  /** @brief Destruction l'obstacle dans l'arborescence des obstacles
  @param name_ nom de l'obstacle */
  void DestroyObstacle( const string &name_ );

  /** @brief Supprime l'obstacle de nom name_
  @param name_ nom de l'obstacle
  @param LC grille de cellules */
  void SupprimeObstacle( const string &name_, LinkedCell *LC );

  /** @brief Ajout d'une cellule � la liste de cellules auxquelles l'obstacle
  est lie
  @param cel_ la cellule */
  void add( Cellule *cel_ );

  /** @brief Vide la liste des cellules auxquelles l'obstacle est lie et
  supprime le pointeur sur l'obstacle dans les cellules de la liste */
  void resetInCells();

  /** @brief Traitement des particules periodiques: creation et destrcution des
  clones periodiques de ces particules
  @param particulesClonesPeriodiques clones periodiques
  @param ParticuleClassesReference classes de r�f�rence de particules
  @param LC grille de cellules */
  virtual void LinkUpdateParticulesPeriodiques(
  	list<Particule*>*particulesClonesPeriodiques,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC ) {}

  /** @brief Traitement des particules periodiques en parallele: seules les
  particules de tag 0 ou 1 sont parcourues, ensuite les infos sont envoyees aux
  autres processeurs pour traitement soit des references periodiques soit des
  clones periodiques
  @param time temps
  @param ClonestoDestroy liste des IDs des clones a detruire
  @param ClonestoParticules liste des IDs des clones a transformer en particules
  	standards
  @param PartRefPerHalozone liste (ID,obstacle ID) de nouvelles particules de
  	reference periodique de la zone de recouvrement
  @param PartRefPerOutDomainHalozone list ID de particules de
  	reference periodique hors du domaine de la zone de recouvrement
  @param InNotRefPerHalozone list (ID,obstacle ID)  de particules ayant perdu
  	leur statut de reference periodique mais toujours dans le domaine
	et la zone de recouvrement
  @param particulesReferencesPeriodiques liste des particules possedant
  	des clones periodiques
  @param ENSparticules particules actives
  @param particulesHalozone particules dans la zone de recouvrement
  @param LC grille de cellules */
  virtual void LinkUpdateParticulesPeriodiques_MPI(
  	Scalar time,
  	list<int>& ClonestoDestroy,
  	list<int>& ClonestoParticules,
  	list<int>& PartRefPerHalozone,
  	list<int>& PartRefPerOutDomainHalozone,
  	list<int>& InNotRefPerHalozone,
  	set<Particule*>* particulesReferencesPeriodiques,
	list<Particule*>* ENSparticules,
  	list<Particule*>* particulesHalozone,
	LinkedCell* LC ) {}

  /** @brief Ajout des nouvelles references periodiques
  @param time temps
  @param particulesReferencesPeriodiques liste des particules possedant
  	des clones periodiques */
  virtual void addNewPeriodicReference_MPI( Scalar time,
  	set<Particule*>* particulesReferencesPeriodiques ) {}

  /** @brief Indique si un obstacle est periodique et intersecte une boite
  @param boite la boite */
  virtual bool obstaclePeridiqueIntersect( BBox const& boite ) const
  	{ return false; }

  /** @brief Indique si il y a lieu de mettre a jour le lien entre l'obstacle et
  le LinkedCell: si oui, la methode renvoie true et met le compteur a zero, si
  non la methode renvoie false et met le compteur est incremente */
  bool performLinkUpdate();

  /** @brief Fr�quence de mise � jour du lien entre l'obstacle
  et le LinkedCell
  @param updateFreq la frequence */
  void setObstacleLinkedCellUpdateFrequency( int const &updateFreq );

  /** @brief Maximum de la valeur absolue de vitesse de l'obstacle
  (en parcourant les points de sa surface) dans chaque direction */
  Vecteur vitesseMaxPerDirection() const;

  /** @brief Update contact map */
  virtual void updateContactMap();

  virtual void setContactMapToFalse();

  /** @brief Does the contact exist in the map, if yes return the pointer to the
  cumulative tangential displacement
  @param tangentialDepl pointer to the cumulative tangential displacement
  @param id id number of the other component */
  virtual bool getContactMemory( std::tuple<int,int,int> const& id,
    Vecteur* &tangent, Vecteur* &prev_normal, Vecteur* &cumulSpringTorque );

  /** @brief Add new contact in the map
  @param tangentialDepl initial tangential displacement
  @param id id number of the other component */
  virtual void addNewContactInMap( std::tuple<int,int,int> const& id,
    Vecteur const& tangent, Vecteur const& prev_normal,
    Vecteur const& cumulSpringTorque );

  /** @brief Increase cumulative tangential displacement with component id
  @param tangentialDepl additional tangential displacement
  @param id id number of the other component */
  virtual void addDeplContactInMap( std::tuple<int,int,int> const& id,
    Vecteur const& tangent, Vecteur const& prev_normal,
    Vecteur const& cumulSpringTorque );

  //@}

  void updateContactMapId( int prev_id, int new_id);


  virtual void copyHistoryContacts( double* &destination, int start_index );


  // Copy existing contact in the map
  virtual void copyContactInMap( std::tuple<int,int,int> const& id,
    bool const& isActive, Vecteur const& tangent, Vecteur const& prev_normal,
    Vecteur const& cumulSpringTorque );

  virtual int getContactMapSize();


  /** @name Processus de restauration */
  //@{
  /** @brief Cree et ajoute l'etat
  @param obsStates liste des etats de tous les obstacles */
  void createState( list<struct ObstacleState*> &obsStates );

  /** @brief Restauration de l'etat
  @param obsStates liste des etats de tous les obstacles */
  void restaureState( list<struct ObstacleState*>& obsStates );
  //@}


  /** @name Methods I/O */
  //@{
  /** @brief Reload de l'obstacle & publication dans son referent
  @param obstacle L'obstacle referent du courant
  @param file Le fichier de persistance */
  virtual void reload( Obstacle &obstacle, istream &file );

  /** @brief Sauvegarde de l'obstacle pour Reload
  @param fileSave Flux de sauvegarde
  @param composant solid object*/
  virtual void write( ostream &fileSave, Composant const* composant = NULL )
      const;

  /** @brief Nombre de points pour post-processing avec Paraview */
  virtual int numberOfPoints_PARAVIEW() const;

  /** @brief Nombre de polygones elementaires pour post-processing avec
  Paraview */
  virtual int numberOfCells_PARAVIEW() const;

  /** @brief Renvoie les points du convexe pour post-processing avec Paraview
  @param transform transformation courante
  @param translation translation du centre de gravite */
  virtual list<Point> get_polygonsPts_PARAVIEW(
  	Vecteur const* translation = NULL ) const;

  /** @brief Ecrit les points du convexe pour post-processing avec Paraview
  @param f Flux de sortie
  @param translation translation du centre de gravite */
  virtual void write_polygonsPts_PARAVIEW( ostream &f,
  	Vecteur const* translation = NULL ) const;

  /** @brief Ecrit le convexe pour post-processing avec Paraview
  @param connectivity connectivite des polygones Paraview
  @param offsets decalage dans la connectivite
  @param cellstype type de polygones Paraview
  @param firstpoint_globalnumber numero global du 1er point
  @param last_offset dernier offset utilise pour le convexe precedent */
  virtual void write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const ;

  /** @brief Ecriture de l'information de position pour le Fluide */
  virtual void writePositionInFluid( ostream &fluid );
  //@}


protected:
  /**@name Parameters */
  //@{
  BBox m_obstacleBox; /**< BBox contenant l'obstacle */
  list<Cellule*> m_inCells; /**< cellules auxquelles l'obstacle est lie */
  int m_LinkUpdate_frequency; /**< fr�quence de mise � jour du lien entre
  	l'obstacle et le LinkedCell */
  int m_LinkUpdate_counter; /**< compteur de mise � jour du lien entre
  	l'obstacle et le LinkedCell */
  bool m_transferToFluid; /**< Obstacle � transf�rer au fluide dans le cas d'un
  	couplage */
  //@}

};

#endif
