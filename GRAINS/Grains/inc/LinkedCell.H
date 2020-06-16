#ifndef _LinkedCell
#define _LinkedCell

#include "AppSec.H"
#include "BBox.H"
#include "FormeVdW.H"

class Cellule;
class Voisins;
class MPIWrapperGrains;

#include <vector>
using namespace std;

// D. RAKOTONIRINA - June 2017 - Modification
struct PointForcePostProcessing
{
  Point geometricPointOfContact;
  Vecteur contactForce;
  Scalar contactDissip;
  Composant* comp0;
  Composant* comp1;
};


/** @brief Gestion de la recherche des contacts par algorithme de Linked-Cell.

    Package APP.

    @author G.FERRER - Institut Francais du Petrole - 2002 - Creation */
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class LinkedCell : public AppSec
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut. */
  LinkedCell();

  /** @brief Destructeur */
  ~LinkedCell();
  //@}


  /** @name Methods */
  //@{
  /** @brief Calcul des forces d'interaction
  Recherche des particules voisines
  <ul>
    <li> Dans la cellule locale: On recherche les voisines situees au dela
    dans la liste. Cela permet de ne considerer un contact entre 2 particules.
    Exemple: soit 4 particules dans la cellule aux positions 1,2,3,4;
    Si on s'interesse a la particule 3, le seul contact possible est 3-4
    Puis quand on considere 4 plus tard dans la boucle, plus aucun contact
    possible n'est recherche et donc le couple 3-4 = 4-3 n'est pris en compte
    qu'une fois
    <li> Dans les cellules adjacentes: Grace au voisinage de contact
    (9 derriere, 3 au-dessus, 1 droite), les
    couples de contact possibles ne sont egalement comptes qu'une seule fois
  </ul>
  @param time temps physique
  @param dt pas de temps granulaire
  @param particules liste de particules */
  virtual void CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules );

  /**
    @brief Recherche des particules voisines pour l'initialisation des
        conditions de cohesion
    @param time temps physique
    @param dt pas de temps granulaire
    @param particules liste de particules
  */
  virtual void InitializeCohesiveForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules );

  /** @brief Find the spheres subjected to lubrication correction
  @param particules liste de particules */

  /** @brief Find the periodic spheres subjected to lubrication correction
  @param particulesClonesPeriodiques liste des clones periodiques */

  /** @brief Calcul des forces d'interaction
  @param particules liste de particules */
  list<struct PointForcePostProcessing>* CalculerForcesPostProcessing(
  	list<Particule*> const* particules, Scalar dt ) const;

  /** @brief Calcul des forces d'interaction pour les clones periodiques
  @param time temps physique
  @param dt pas de temps granulaire
  @param particulesClonesPeriodiques liste des clones periodiques */
  void CalculerForcesClonesPeriodiques( Scalar time, Scalar dt,
  	list<Particule*>* particulesClonesPeriodiques );

  /** @brief La particule est elle en contact avec un composant de l'App ?
  @return TRUE si contact
  @param particule particule */
  bool isContact( const Particule *particule ) const;

  /** @brief La particule est elle en contact avec un composant de l'App ?
  variante VdW
  @return TRUE si contact
  @param particule particule */
  bool isContactVdW( const Particule *particule ) const;

  /** @brief La particule est elle proche d'un composant de l'App ?
  @return TRUE si contact
  @param particule particule */
  bool isProche( const Particule *particule ) const;

  /** @brief La particule est elle proche d'un composant de l'App ?
  variante VdW
  @return TRUE si contact
  @param particule particule */
  bool isProcheVdW( const Particule *particule ) const;

  /** @brief La particule est elle en contact avec un obstacle periodique
  @return Un contact eventuel.
  @param particule La particule etudiee. */
  bool isContactVdW_ObstaclePeriodique( const Particule* particule ) const;

  /** @brief Association de la particule avec l'algorithme
  sans controle de contact.
  @param particule particule */
  void Link( Particule* particule );

  /** @brief Association de l'obstacle avec l'algorithme.
  @param obstacle obstacle */
  void Link( Obstacle *obstacle );

  /** @brief Actualisation des associations des particules avec l'algorithme
  @param time temps physique
  @param dt pas de temps granulaire
  @param particules liste de particules */
  void LinkUpdate( Scalar time, Scalar dt,
  	list<Particule*>* particules ) throw (ErreurSimulation);

  /** @brief Actualisation d'une particule active avec l'algorithme.
  @param particule particule */
  void LinkUpdateActiveParticule( Particule* particule )
  	throw ( ErreurSimulation );

  /** @brief Transfert de la particule de la cellule A vers B
  @param particule La particule a transferer
  @param vecteur   Le vecteur de translation entrainant le transfert */
  void shift( Particule& particule, const Vecteur& vecteur );

  /** @brief Operateur <<
  @param f flux de sortie
  @param LC objet Linked Cell */
  friend ostream& operator << ( ostream &f, const LinkedCell &LC );

  /** @brief Appartenance du centre de gravite d'une particule au Linked Cell
  @param position position du centre de gravite */
  bool isInLinkedCell( const Point &position ) const;

  /** @brief Appartenance du centre de gravite d'une particule au Linked Cell
  @param gx coordonn�e x du centre de gravite
  @param gy coordonn�e y du centre de gravite
  @param gz coordonn�e z du centre de gravite */
  bool isInLinkedCell( const double &gx, const double &gy, const double &gz )
  	const;

  /** @brief Suppression de la particule de l'APPlication
  @param particule La particule a supprimer */
  void remove( Particule* particule );

  /** @brief Suppression de l'obstacle de l'APPlication
  @param obs obstacle a supprimer */
  void remove( MonObstacle* obs );

  /** @brief Mise a jour des tags des particules interieures
  @param time temps physique
  @param ENSparticules particules actives
  @param particulesHalozone particules actives dans la zone de recouvrement
  @param wrapper wrapper MPI pour les simulations parall�les */
  void updateInteriorTag( Scalar time,
  	list<Particule*>* ENSparticules,
  	list<Particule*>* particulesHalozone,
	MPIWrapperGrains const* wrapper=NULL );

  /** @brief Mise a jour des tags des particules de la zone de recouvrement et
  des clones
  @param time temps physique
  @param particulesHalozone particules actives dans la zone de recouvrement
  @param particulesClones clones
  @param wrapper wrapper MPI pour les simulations parall�les */
  void updateHalozoneCloneTag( Scalar time,
  	list<Particule*>* particulesHalozone,
	list<Particule*>* particulesClones,
	MPIWrapperGrains const* wrapper=NULL );

  /** @brief Suppression des clones sortis de la grille
  @param time temps physique
  @param particulesClones clones
  @param particulesReferencesPeriodiques liste des particules possedant
  	des clones periodiques
  @param ENSparticules particules actives
  @param wrapper wrapper MPI pour les simulations parall�les */
  void DestroyOutOfDomainClones( Scalar time,
  	list<Particule*>* particulesClones,
  	set<Particule*>* particulesReferencesPeriodiques,
	list<Particule*>* ENSparticules,
	MPIWrapperGrains const* wrapper=NULL );

  /** @brief Traitement des particules periodiques: creation et destruction des
  clones periodiques de ces particules
  @param time temps physique
  @param ENSparticules particules actives
  @param particulesClonesPeriodiques clones periodiques
  @param ParticuleClassesReference particules de r�f�rence */
  void LinkUpdateParticulesPeriodiques( Scalar time,
  	list<Particule*>* ENSparticules,
	list<Particule*>* particulesClonesPeriodiques,
	vector<Particule*> const* ParticuleClassesReference );

  /** @brief Traitement des particules periodiques en parallele
  @param time temps physique
  @param ClonestoDestroy liste des IDs des clones a detruire
  @param ClonestoParticules liste des IDs des clones a transformer en particules
  	standards
  @param PartRefPerHalozone liste (ID,obstacle ID) de particules de
  	reference periodique de la zone de recouvrement
  @param PartRefPerOutDomainHalozone list ID de particules de
  	reference periodique hors du domaine de la zone de recouvrement
  @param InNotRefPerHalozone list (ID,obstacle ID)  de particules ayant perdu
  	leur statut de reference periodique mais toujours dans le domaine
	et la zone de recouvrement
  @param particulesClonesPeriodiques clones periodiques
  @param particulesReferencesPeriodiques liste des particules possedant
  	des clones periodiques
  @param ENSparticules particules actives
  @param particulesHalozone particules dans la zone de recouvrement */
  void LinkUpdateParticulesPeriodiques_MPI_Step1( Scalar time,
	list<int>& ClonestoDestroy,
  	list<int>& ClonestoParticules,
  	list<int>& PartRefPerHalozone,
  	list<int>& PartRefPerOutDomainHalozone,
  	list<int>& InNotRefPerHalozone,
	list<Particule*>* particulesClonesPeriodiques,
  	set<Particule*>* particulesReferencesPeriodiques,
	list<Particule*>* ENSparticules,
  	list<Particule*>* particulesHalozone );

  /** @brief D�finition le voisinage complet des cellules */
  void affecteVoisinageComplet();

  /** @brief Ajout des nouvelles references periodiques
  @param time temps physique
  @param particulesReferencesPeriodiques liste des particules possedant
  	des clones periodiques */
  void addNewPeriodicReference_MPI( Scalar time,
  	set<Particule*>* particulesReferencesPeriodiques );

  /** @brief Contient une paroi periodique */
  bool intersectObstaclePeriodique() const;
  //@}


  /**@name Methods Get */
  //@{
  /** @brief Cellule associee a la position designee
  @param position La position demandee
  @return La cellule contenant la particule */
  Cellule* getCellule( const Point &position ) const;

  /** @brief Renvoie la taille d'une cellule dans une direction
  @param dir direction */
  Scalar getCelluleSize( int const& dir ) const;

  /** @brief Renvoie un pointeur sur le vecteur de cellules */
  vector<Cellule*> const* getAllCellules() const;
  //@}


  /**@name Methods Set */
  //@{
  /** @brief Definition de l'espace du Linked_Cell.
  @param arete_ Longueur de l'arete de reference des cellules. */
  void define( double arete_ );

  /** @brief Definition de l'espace du Linked_Cell.
  @param arete_ Longueur de l'arete de reference des cellules.
  @param nprocsdir nb de procs dans chaque direction
  @param MPIcoords coordonnes dans la topologie cartesienne MPI
  @param voisins voisinage dans la topologie cartesienne MPI
  @param MPIperiod periodicite de la topologie MPI */
  void define( double arete_, int const* nprocsdir, int const* MPIcoords,
  	Voisins const* voisins, int const* MPIperiod );
  //@}


private:
  /** @name Methods */
  //@{
  /** @brief Reference de la cellule liee a la position indiquee.
  @param i indice dans la direction X
  @param j indice dans la direction Y
  @param k indice dans la direction Z */
  Cellule* getCellule( int i, int j, int k ) const;

  /** @brief Numero de la cellule liee a la position indiquee.
  @param i indice dans la direction X
  @param j indice dans la direction Y
  @param k indice dans la direction Z */
  int getCelluleNumber( int i, int j, int k ) const;

  /** @brief D�finition le voisinage des cellules */
  void affecteVoisinage();

  /** @brief Mise a jour de l'association de l'obstacle avec l'algorithme.
  @param time temps physique
  @param dt pas de temps
  @param myObs obstacle */
  void LinkUpdate( Scalar time, Scalar dt, MonObstacle *myObs );
  //@}


  /** @name Parameters */
  //@{
  vector<Cellule*> m_cellules; /**< Cellules definissant le recipient */
  int m_nb; /**< Nombre total de cellules */
  int m_nbi; /**< Nombre de cellules direction x */
  int m_nbj; /**< Nombre de cellules direction y */
  int m_nbk; /**< Nombre de cellules direction z */
  double m_arete_X; /**< Longueur de reference des cellules
  	dans la direction X */
  double m_arete_Y; /**< Longueur de reference des cellules
  	dans la direction Y */
  double m_arete_Z; /**< Longueur de reference des cellules
  	dans la direction Z */
  Point m_LC_origine_globale; /**< Origine globale de la grille */
  Point m_LC_origine_locale; /**< Origine locale de la grille */
  double m_xmin; /**< coordonnee minimale dans la direction X */
  double m_ymin; /**< coordonnee minimale dans la direction Y */
  double m_zmin; /**< coordonnee minimale dans la direction Z */
  double m_xmax; /**< coordonnee maximale dans la direction X */
  double m_ymax; /**< coordonnee maximale dans la direction Y */
  double m_zmax; /**< coordonnee maximale dans la direction Z */
  BBox *m_paveEtendu; /**< la boite union de toutes les cellules, agrandie d'une
  	demi-cellule dans toutes les directions */
  //@}
};

#endif
