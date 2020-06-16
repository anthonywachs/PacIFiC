#ifndef _ObstacleAbsorbant
#define _ObstacleAbsorbant

#include "MonObstacle.H"

#include "ReaderXML.hh"



/** @brief Gestion d'un obstacle ayant la propriete d'absorber les particules.

    Une particule en contact avec ce type d'obstacle ne recoit pas de force
    de repulsion. Lorsque le critere de penetration est atteint, l'etat 
    d'activite de la particule devient inactif. Le logiciel est alors en 
    charge de l'action a accomplir.

    Pour simplifier la mise en donnees, on continue a definir l'ensemble des
    proprietes des obstacles.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation */
// ============================================================================
class ObstacleAbsorbant : public MonObstacle
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut
  @param s Nom eventuel de l'obstacle  */
  ObstacleAbsorbant( const string &s = "" );

  /** @brief Constructeur
  @param root Le noeud "<Obstacle Type="Absorbant">" */
  ObstacleAbsorbant( DOMNode *root );

  /** @brief Destructeur */
  ~ObstacleAbsorbant();
  //@}


  /**@name Methods */
  //@{
  /** @brief Contact entre l'obstacle et une particule.
  Si le contact est actif, pas d'ajout de contact au deux composants.
  Si un choc de croute est detecte, le monolithe devient inactif.     
  @param voisin Le composant (monolithe) a etudier
  @param dt Increment de temps
  @param temps temps physique
  @param LC pointer to the Linked-Cell
  @exception ErreurContact Si un choc entre croute est detecte
  @return Le contact determine (NULL pour ce composant) */
  virtual void InterAction( Composant* voisin, 
	double dt, double const& temps, LinkedCell *LC ) throw (ErreurContact);

  /** @brief Contact entre l'obstacle et une particule pour le post processing
  Si le contact est realise, on ajoute le contact au deux composants. 
  @param voisin Le composant (monolithe) a etudier
  @param listOfContacts liste de contacts */
  virtual void InterActionPostProcessing( Composant* voisin, 
	list<struct PointForcePostProcessing>* listOfContacts ) 
	throw (ErreurContact);
  //@}


  /**@name Methods I/O */
  //@{
  /** @brief Reload de l'obstacle & publication dans son referent
  @param obstacle L'obstacle referent du courant
  @param file Le fichier de persistance */
  void reload( Obstacle &obstacle, istream &file );

  /** @brief Sauvegarde de l'obstacle pour Reload
  @param fileSave Flux de sauvegarde
  @param composant Reference Composant utile en cas de particule composite
  et NULL pour les particules convexes */
  void write( ostream &fileSave, Composant const* composant = NULL ) const;

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
};

#endif
