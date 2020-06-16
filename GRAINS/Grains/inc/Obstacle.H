#ifndef _Obstacle
#define _Obstacle

#include "Composant.H"
#include "CineObstacle.H"
#include "CineObstacle_F.H"
#include "Basic.H"
#include "Vecteur.H"
using namespace solid;

class MonObstacle;
class ObstacleChargement;
class Quaternion;
class LinkedCell;
class LinkedBox;
class App;


/** @brief Sauvegarde de l'etat de l'obstacle
    @author A.WACHS - IFP Energies nouvelles - 2012 - Creation */
struct ObstacleState
{
  string nom; /**< nom de l'obstacle */
  ConfigurationMemento* memento_config; /**< configuration de l'obstacle */
  CineObstacleMemento* memento_cine; /**< cinematique de l'obstacle */
};


/** @brief Description d'un obstacle.

    Ce type d'objet ne permet pas son franchissement par une particule
    qui rebondit a son contact. Un obstacle peut etre fixe ou mobile.
    Dans ce dernier cas, on connait par donnees sa cinematique.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation */
// ============================================================================
class Obstacle : public Composant
{
public:
  /** @name Constructeurs */
  //@{
  /** @brief Constructeur par defaut.
  @param s Nom eventuel de l'obstacle
  @param autonumbering obstacle autonumbering */
  Obstacle( const string &s = "obstacle", const bool &autonumbering = true );

  /** @brief Constructeur par copie d'un Composant
  @param copie objet � copier
  @param s type de composant, ici "obstacle" */
  Obstacle( Composant &copie, char const* s = "obstacle" );

  /** @brief Destructeur.*/
  virtual ~Obstacle();
  //@}


  /** @name Methods Pure Virtual */
  //@{
  /** @brief Ajout d'un obstacle (unique ou composite)
  a l'arborescence Obstacle.
  @param obstacle L'arborescence Obstacle construit au prealable. */
  virtual void append( Obstacle* obstacle ) = 0;

  /** @brief Deplacement de l'obstacle.
  @return liste des obstacles deplaces
  @param temps temps physique
  @param dt pas de temps
  @param b_deplaceCine_Comp deplacement du composite dont l'obstacle fait
  partie (cinematique)
  @param b_deplaceF_Comp deplacement du composite dont l'obstacle fait
  partie (force) */
  virtual list<MonObstacle*> Deplacer( Scalar temps,
	Scalar dt, const bool &b_deplaceCine_Comp,
	const bool &b_deplaceF_Comp ) = 0;

  /** @brief L'obstacle a t'il le nom indique.
  @return this si le nom correspond a l'obstacle. */
  virtual const Obstacle* getNom( const string &nom_ ) const = 0;

  /** @brief Liste des obstacles.
  @return Pour un obstacle composite on renvoie les obstacles de reference,
  pour un obstacle simple on renvoie lui-meme. */
  virtual list<MonObstacle*> getObstacles() = 0;

  /** @brief Liste des obstacles � transmettre au fluide
  @return Pour un obstacle composite on renvoie les obstacles de reference,
  pour un obstacle simple on renvoie lui-meme. */
//  virtual list<MonObstacle*> getObstaclesToFluid() = 0;
  virtual list<Obstacle*> getObstaclesToFluid() = 0;

  /** @brief Rotation de l'obstacle.
  @param rotation Quaternion decrivant la rotation a effectuer. */
  virtual void Rotate( const Quaternion &rotation ) = 0;

  /** @brief Suppression de l'obstacle dans la zone specifie.
  @param box boite contenant l'obstacle */
  virtual void Suppression( const BBox &box ) = 0;

  /** @brief Translation de l'obstacle
  @param translation Le vecteur de translation. */
  virtual void Translate( const Vecteur &translation ) = 0;

  /** @brief Destruction l'obstacle dans l'arborescence des obstacles
  @param name_ nom de l'obstacle */
  virtual void DestroyObstacle( const string &name_ ) = 0;

  /** @brief Supprime l'obstacle de nom name_
  @param name_ nom de l'obstacle
  @param LC grille de cellules */
  virtual void SupprimeObstacle( const string &name_, LinkedCell *LC ) = 0;
  //@}


  /** @name Methods Set */
  //@{
  /** @brief Affectation de la cinematique
  @param cinematique La cinematique a affecter a l'obstacle (vitesse de
  translation et de rotation) */
  void setCinematique( CineObstacle &cinematique );

  /** @brief Affectation de la vitesse de l'obstacle
  @param vitesseTranslation vitesse de translation
  @param vitesseRotation vitesse de rotation */
  void setVitesse( Vecteur const* vitesseTranslation,
  	Vecteur const* vitesseRotation );

  /** @brief Affectation de l'indicateur pour post-traitement Paraview */
  void setIndicator( double const& value ) { m_indicator = value; }

  /** @brief Initialize all contact map entries to false */
  virtual void setContactMapToFalse();
  //@}


  /** @name Accessors */
  //@{
  /** @brief Nom de l'Obstacle
  @return Le nom */
  string getName() const { return m_nom; }

  /** @brief Vitesse du composant en un point (appartenant au composant a priori
  mais cette condition n'est pas v�rifi�e)
  @return La vitesse
  @param pt le point de definir la partie rotation avec la bras de levier par
  rapport au centre de gravite du composant */
  Vecteur getVitesse( const Point &pt ) const;

  /** @brief Vitesse de rotation. */
  Vecteur const* getVitesseRotation() const;

  /** @brief Vitesse de translation. */
  Vecteur const* getVitesseTranslation() const;

  /** @brief Renvoie l'indicateur pour post-traitement Paraview */
  double getIndicator() const;

  /** @brief Renvoie le type de l'obstacle */
  virtual string getObstacleType() = 0;

  /** @brief Renvoie le torseur applique */
  virtual Torseur const* getTorseur();
  //@}


  /** @name Processus de restauration */
  //@{
  /** @brief Sauvegarde de l'etat */
  virtual void saveState();

  /** @brief Cree et ajoute l'etat
  @param obsStates liste des etats de tous les obstacles */
  virtual void createState( list<struct ObstacleState*> &obsStates ) = 0;

  /** @brief Restauration de l'etat */
  virtual void restaureState();

  /** @brief Restauration de l'etat
  @param obsStates liste des etats de tous les obstacles */
  virtual void restaureState( list<struct ObstacleState*>& obsStates ) = 0;
  //@}


  /** @name Methods Virtual */
  //@{
  /** @brief Association du chargement a l'obstacle.
  @return Vrai si l'obstacle et le chargement sont lies.
  @param chargement Le chargement candidat. */
  virtual bool Associer( ObstacleChargement &chargement );

  /** @brief Association du chargement a l'obstacle.
  @return Vrai si l'obstacle et le chargement sont lies.
  @param chargement Le chargement candidat. */
  virtual bool Associer( ObstacleChargement_F &chargement );

  /** @brief Remise a 0 de la cinematique de l'obstacle. */
  virtual void resetCinematique();

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

  /** @name Methods */
  //@{
  /** @brief Decomposition de la cinematique sur la Cinematique de l'obstacle.
  @param voisine Cinematique a appliquer sur l'obstacle
  @param levier Bras de levier de la cinematique a appliquer */
  void Decompose( const CineObstacle &voisine, const Vecteur &levier );

  /** @brief Decomposition de la cinematique sur la cinematique de l'obstacle.
  @param voisine Cinematique a appliquer sur l'obstacle
  @param centre point */
  void Decompose( const CineObstacle_F &voisine, const Point &centre );

  /** @brief L'obstacle s'est il deplace ?
  @return TRUE si deplacement */
  bool hasMoved() const;
  //@}


  /**@name Methods Static */
  //@{
  /** @brief Nombre d'obstacles unitaires */
  static int getNombre();

  /** @brief Deplacement geometrique des obstacles
  param depObs doit on deplacer geometriquement les obstacles */
  static void setDeplaceObstacle( bool const& depObs );

  /** @brief Deplacement geometrique des obstacles */
  static bool getDeplaceObstacle();
  //@}


  /**@name Methods I/O */
  //@{
  /** @brief Reload du groupe d'Obstacle & publication dans son referent
  @param obstacle L'obstacle referent du courant
  @param file Le fichier de persistance */
  virtual void reload( Obstacle &obstacle, istream &file ) = 0;

  /** @brief Ecriture de l'identite du composant */
  void writeIdentity( ostream &file ) const;

  /** @brief Ecriture de l'information de position.
  @param position Flux recepteur. */
  virtual void writePosition( ostream &position );

  /** @brief Ecriture de l'information statique.
  @param statique Flux recepteur. */
  virtual void writeStatique( ostream &statique );

  /** @brief Mise a jour de l'indicateur pour post-traitement Paraview
  @param temps temps physique
  @param dt pas de temps */
  virtual void updateIndicator( Scalar temps, Scalar dt ) {}

  /** @brief Contact entre composant et composant.
  Si le contact est realise, on ajoute le contact au deux composants.
  @exception ErreurContact Si un choc entre croute est detecte
  @param voisin Le composant (monolithe) a etudier
  @param dt Increment de temps
  @param temps temps physique
  @param LC grille de cellules */
  virtual void InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC ) throw (ErreurContact);

  /** @brief Contact entre l'obstacle et une particule pour le post processing
  Si le contact est realise, on ajoute le contact au deux composants.
  @param voisin Le composant (monolithe) a etudier
  @param listOfContacts liste de contacts */
  virtual void InterActionPostProcessing( Composant* voisin, Scalar dt,
	list<struct PointForcePostProcessing>* listOfContacts )
	throw (ErreurContact);

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
  /** @name Parameters */
  //@{
  CineObstacle m_cinematique; /**< Cinematique en deplacement */
  CineObstacle_F m_confinement; /**< Cinematique en force */
  string m_nom; /**< Nom de l'obstacle pour utiliser en reference */
  bool m_deplace; /**< Indicateur d'activation de deplacement */
  double m_indicator; /**< indicateur pour le post-traitement; jusqu'ici,
  	utilis� pour le mouvement de rotation de composite dans Paraview */
  string m_ObstacleType;
  //@}

  /** @name Parameters Static */
  //@{
  static int m_nbreObstacles; /**< Nombre d'obstacles unitaires */
  static bool m_DeplaceObstacle; /**< si false, permet d'imposer une cinematique
  	sans deplacer les obstacles, utile dans les cas periodiques */
  static bool m_isConfinement; /**< si true, chargement en force */
  //@}

};

#endif
