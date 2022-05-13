#ifndef _CompObstacle
#define _CompObstacle

#include "Obstacle.H"
#include "MonObstacle.H"

#include <list>
using namespace std;

#include "ReaderXML.hh"


/** @brief Gestion des obstacles composites : COMPosed OBSTACLE.

    Un obstacle composite est un obstacle pour lequel chaque composant est
    associe a une meme cinematique. Il represente par exemple un cylindre
    decrit par ces n facettes auquel on desire appliquer une rotation.
    A partir d'obstacles unitaires ou composites, on decrit un obstacle plus
    complexe auquel on peut associer un mouvement de rotation et de
    translation. Toutefois, chaque composant conserve en plus de la
    cinematique globale sa cinematique propre.

    @author G.FERRER - Institut Francais du Petrole - 2002 - Creation
    @author D. RAKOTONIRINA - IFP Energies nouvelles - 2014 - Modification */
// ============================================================================
class CompObstacle : public Obstacle
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut
  @param s Nom eventuel de l'obstacle */
  CompObstacle( const string &s = "" );

  /**
    @brief Constructeur du groupe d'obstacles.
    @param root Le noeud "MultiObstacles" de groupes d'obstacles.
    @details Dynamically allocates an instance of FormeVdW
        (@ref m_geoFormeVdw).<BR>
    Reads the %CompObstacle name.<BR>
    For each composing-obstacle (can be %Obstacle, Composite or Periode type)
        it calls Obstacle_BuilderFactory::create, then add it to @ref
        m_obstacles<BR>
    Finally calls @ref EvalPosition to compute Composant geometrical-center
        and set it (by using @ref m_geoFormeVdw object)
  */
  CompObstacle( DOMNode* root );

  /** @brief Destructeur */
  ~CompObstacle();
  //@}


  /** @name Methods Get */
  //@{
  /**
    @brief L'obstacle a t'il le nom indique ?
    @detail Si le nom n'est pas valide, on poursuit la recherche dans les feuilles.
    @return this si le nom correspond a l'obstacle.
    @param nom_ nom de l'obstacle
  */
  const Obstacle* getNom( const string &nom_ ) const;

  /**
    @brief Liste des obstacles.
    @return On renvoie les obstacles de reference.
  */
  list<MonObstacle*> getObstacles();

  /**
    @brief Liste des obstacles � transmettre au fluide
    @return Pour un obstacle composite on renvoie les obstacles de reference,
      pour un obstacle simple on renvoie lui-meme.
  */
//  list<MonObstacle*> getObstaclesToFluid();
  list<Obstacle*> getObstaclesToFluid();

  /**
    @brief Renvoie le type de l'obstacle
  */
  string getObstacleType();
  //@}

  /** @name Methods Set */
  //@{
  /** @brief Initialize all contact map entries to false */
  virtual void setContactMapToFalse();
  //@}

  /** @name Methods */
  //@{
  /** @brief Association du chargement a l'obstacle
  ou un de ses sous-ensembles.
  L'association est definie par le nom de l'obstacle dans le chargement.
  @return Vrai si l'obstacle et le chargement sont lies.
  @param chargement Le chargement candidat. */
  bool Associer( ObstacleChargement &chargement );

  /** @brief Association du chargement a l'obstacle
  ou un de ses sous-ensembles.
  L'association est definie par le nom de l'obstacle dans le chargement.
  @return Vrai si l'obstacle et le chargement sont lies.
  @param chargement Le chargement candidat. */
  bool Associer( ObstacleChargement_F &chargement );

 /** @brief Ajout d'un obstacle (unique ou composite)
  a l'arborescence Obstacle.
  @param obstacle L'arborescence Obstacle construit au prealable. */
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

  /** @brief Le composant est-il un obstacle composite?
  @return Vrai si c'est un obstacle composite. */
  virtual bool isCompObstacle() const {return ( true ); }

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

  /** @brief Remise a 0 de la cinematique de l'obstacle. */
  void resetCinematique();

  /** @brief Construction d'un objet de type obstacle par lecture sur fichier.
  Cette methode est utilisee pour reconstruire un objet persistant. */
  // void reload(istream &file);

  /** @brief Rotation de la particule.
  @param rotation Quaternion decrivant la rotation a effectuer. */
  void Rotate( const Quaternion &rotation );

  /** @brief Suppression de l'obstacle dans la zone specifie.
  Supprime les obstacles contenus se trouvant dans la zone. */
  void Suppression( const BBox &box );

  /** @brief Translation de l'obstacle
  @param translation Le vecteur de translation. */
  void Translate( const Vecteur &translation );

  /** @brief Initialise le torseur des efforts sur le composant
  @param withWeight avec ou sans le poids du composant */
  void InitializeForce( bool const& withWeight );

  /** @brief Destruction l'obstacle dans l'arborescence des obstacles
  @param name_ nom de l'obstacle */
  void DestroyObstacle( const string &name_ );

  /** @brief Supprime l'obstacle de nom name_
  @param name_ nom de l'obstacle
  @param LC grille de cellules */
  void SupprimeObstacle( const string &name_, LinkedCell *LC );

  /** @brief Renvoie le torseur applique */
  Torseur const* getTorseur();

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

  /** @brief Update contact map */
  virtual void updateContactMap();

  /** @brief Does the contact exist in the map, if yes return the pointer to the
  cumulative tangential displacement
  @param tangentialDepl pointer to the cumulative tangential displacement
  @param id id number of the other component */
  virtual bool ContactInMapIsActive( double* &tangentialDepl, int const& id );

  /** @brief Add new contact in the map
  @param tangentialDepl initial tangential displacement
  @param id id number of the other component */
  virtual void addNewContactInMap( double const& tangentialDepl,
  	int const& id );

  /** @brief Increase cumulative tangential displacement with component id
  @param tangentialDepl additional tangential displacement
  @param id id number of the other component */
  virtual void addDeplContactInMap( double const& tangentialDepl,
  	int const& id );
  //@}

  virtual void copyHistoryContacts( double* &destination, int start_index );

  // Copy existing contact in the map
  virtual void copyContactInMap( std::tuple<int,int,int> const& id,
    bool const& isActive, int const& nbDtRemember, int const& nbCumulTangent,
    Vecteur const& tangent, Vecteur const& prev_normal,
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
  /** @brief Reload du groupe d'Obstacle & publication dans son referent
  @param obstacle L'obstacle referent du courant
  @param file Le fichier de persistance */
  void reload( Obstacle &obstacle, istream &file );

  /** @brief Sauvegarde des obstacles pour Reload
  @param fileSave Flux de sauvegarde
  @param composant solid object */
  virtual void write( ostream &fileSave, Composant const* composant = NULL )
      const;

  /** @brief Ecriture de l'information de position.
  @param position Flux recepteur. */
  void writePosition( ostream &position );

  /** @brief Ecriture de l'information statique.
  @param statique Flux recepteur
  @param composant solid object */
  void writeStatique( ostream &statique, Composant* composant = NULL );

  /** @brief Ecriture de l'objet dans sa configuration courante pour
  post-processing avec GMV
  @param fileOut Flux de sortie. */
  void GMVoutput( ostream &fileOut ) const;

  /** @brief Mise a jour de l'indicateur pour post-traitement Paraview
  @param temps temps physique
  @param dt pas de temps */
  void updateIndicator( Scalar temps,Scalar dt );

  /** @brief Ecriture de l'information de position pour le Fluide */
  virtual void writePositionInFluid( ostream &fluid );
  //@}


protected:
  /** @name Methods */
  //@{
  /** @brief Definition du centre de gravite */
  void EvalPosition();
  //@}

  /** @name Parameters */
  //@{
  list<Obstacle*> m_obstacles; /**< Liste des obstacles composant
  	le CompObstacle */
  //@}


private:
  /** @name Methods */
  //@{

  //@}
};

#endif
