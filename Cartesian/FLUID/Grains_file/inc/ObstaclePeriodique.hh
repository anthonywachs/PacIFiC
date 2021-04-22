#ifndef _ObstaclePeriodique
#define _ObstaclePeriodique

#include "MonObstacle.H"
#include "Vecteur.H"
#include "ReaderXML.hh"

class App;


/** @brief Obstacle Periodique.

    Ce type d'obstacle ne genere pas de Contact reel sur les Particules.
    Il entraine la creation d'une Particule "fantome" reliee a la Particule
    de reference.

    @author GRAINS Project - IFP - 2007 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ObstaclePeriodique : public MonObstacle
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut
  @param s Nom eventuel de l'obstacle */
  ObstaclePeriodique( const string &s = "" );

  /** @brief Constructeur
  @param root Le noeud "<Periodique>"
  @param name Nom de l'obstacle periodique
  @param direction Vecteur de periode utilise par l'obstacle */
  ObstaclePeriodique( DOMNode* root,
	const string  &name,
	const Vecteur &direction );

  /** @brief Destructeur */
  ~ObstaclePeriodique();
  //@}


  /** @name Methods */
  //@{
  /** @brief Deplacement de l'obstacle
  (pas de deplacement pour obstacle periodique)
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

  /** @brief Acces a l'obstacle associe
  @return L'obstacle associe */
  ObstaclePeriodique* getAssocie() const;

  /** @brief Acces au vecteur de translation periodique
  @return le Vecteur de periode */
  Vecteur const* getPeriode() const;

  /** @brief Renvoie le pointeur sur l'obstacle si il est periodique, renvoie
  NULL sinon */
  ObstaclePeriodique const* getObstaclePeriodic() const { return this; }

  /** @brief Determination du contact avec la Particule
  @param voisin La particule etudiee
  @param dt Increment de temps
  @param temps temps physique
  @param LC grille de cellules
  @exception ErreurContact jamais declanche pour un obstacle periodique */
  virtual void InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC ) throw (ErreurContact);

  /** @brief Contact entre l'obstacle et une particule pour le post processing
  Si le contact est realise, on ajoute le contact au deux composants.
  @param voisin Le composant (monolithe) a etudier
  @param listOfContacts liste de contacts */
  virtual void InterActionPostProcessing( Composant* voisin,
	list<struct PointForcePostProcessing>* listOfContacts )
	throw (ErreurContact);

  /** @brief Ajout de la reference a l'ObstaclePeriodique associe
  @param obstacle L'Obstacle periodique asssocie */
  void setAssocie( ObstaclePeriodique *obstacle );

  /** @brief Traitement des particules periodiques: creation et destruction des
  clones periodiques de ces particules
  @param particulesClonesPeriodiques clones periodiques
  @param ParticuleClassesReference classes de r�f�rence de particules
  @param LC grille de cellules */
  void LinkUpdateParticulesPeriodiques(
  	list<Particule*>* particulesClonesPeriodiques,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC );

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
  void LinkUpdateParticulesPeriodiques_MPI(
  	Scalar time,
  	list<int>& ClonestoDestroy,
  	list<int>& ClonestoParticules,
  	list<int>& PartRefPerHalozone,
  	list<int>& PartRefPerOutDomainHalozone,
  	list<int>& InNotRefPerHalozone,
  	set<Particule*>* particulesReferencesPeriodiques,
	list<Particule*>* ENSparticules,
  	list<Particule*>* particulesHalozone,
	LinkedCell* LC );

  /** @brief Ajout des nouvelles references periodiques
  @param time temps
  @param particulesReferencesPeriodiques liste des particules possedant
  	des clones periodiques */
  void addNewPeriodicReference_MPI( Scalar time,
  	set<Particule*>* particulesReferencesPeriodiques );

  /** @brief Indique si un obstacle periodique intersecte une boite
  @param boite la boite */
  virtual bool obstaclePeridiqueIntersect( BBox const& boite ) const;
  //@}


  /** @name Methods I/O */
  //@{
  /** @brief Reload de l'obstacle & publication dans son referent
  @param obstacle L'obstacle referent du courant
  @param file Le fichier de persistance */
  void reload( Obstacle &obstacle, istream &file );

  /** @brief Sauvegarde de l'obstacle pour Reload
  @param fileSave Flux de sauvegarde
  @param composant Reference Composant : utile en cas de particule composite
  et NULL pour les particules convexes */
  void write( ostream &fileSave, Composant const* composant = NULL ) const;

  /** @brief Est ce que le composant est un obstacle periodique ? */
  bool isObsPeriodique() const {return ( true ); }

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
  //@}


private:
  /** @name Methods */
  //@{
  /** @brief Renvoie si la particule est hors ou dans le domaine periodique
  @param time temps
  @param reference la particule */
  bool isIn( Particule* reference, Scalar time );
  //@}


protected:
  /** @name Parameters */
  //@{
  Vecteur m_vecteur; /**< Vecteur definissant la periode de l'obstacle */
  ObstaclePeriodique *m_associe; /**< Obstacle periodique associe */
  //@}
};

#endif
