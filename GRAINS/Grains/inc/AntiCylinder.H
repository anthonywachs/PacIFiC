#ifndef _ANTICYLINDER_H_
#define _ANTICYLINDER_H_

#include "Convex.H"

#include "ReaderXML.hh"


/** @brief Gestion des convexes de forme cylindrique.
    
    @author M.BERNARD - April 2015 - Creation */
// ============================================================================
class AntiCylinder : public Convex 
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur avec initialisation
  @param r le rayon du cylindre
  @param h la hauteur du cylindre */
  AntiCylinder( Scalar r = 0, Scalar h = 0 );

  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant */
  AntiCylinder( istream& fileIn );

  /** @brief Constructeur
  @param root Le noeud "Cylindre" */
  AntiCylinder( DOMNode* root );

  /** Destructeur */
  ~AntiCylinder();
  //@}
  

  /** @name Methods */
  //@{
  /** @brief Determine l'inertie d'un cylindre
  @param inertie L'inertie du convexe
  @param inertie_1 L'inertie inverse du convexe*/
  bool BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const;
  
  /** @brief Renvoie le choix de retrecissement ou pas*/
  int getShrinkingChoice() const;
  
  /** @brief fixe le rayon courant 
  @param CurrentRadius rayon actuel */
  void setShrinkingRadius(Scalar CurrentRadius); 

  /** @brief Clone d'un cylindre
  @return le clone d'un cylindre */
  Convex* clone() const;

  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return ANTICYLINDER;};

  /** @brief Points decrivant l'enveloppe
  @return La liste des points sommets */
  vector<Point> getEnveloppe() const;

  /** @brief Description des faces
  @return La description des faces via les indices des points sommets */
  vector<vector<int> > const* getFaces() const;

  /** @brief Renvoie le nombre de sommets ou un code équivalent */
  int getNbCorners() const; 
     
  /** @brief Fonction support : Determine le point du support convexe de la 
  forme se situant dans une direction et un sens donne. 
  Pour le cylindre : On choisit le point appartenant a l'un des deux 
  cercles de base le plus 'proche' du vecteur.
  @return le point sur l'enveloppe convexe obtenu
  @param v la direction et le sens */
  Point support( const Vecteur& v ) const;

  /** @brief Determine le volume d'un cylindre
  @return le volume*/
  Scalar getVolume() const;

  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  void printClass( ostream &fileOut ) const;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  void readClass( istream &fileIn );
  
  /** @brief Nombre de points pour post-processing avec Paraview */
  int numberOfPoints_PARAVIEW() const;
  
  /** @brief Nombre de polygones elementaires pour post-processing avec 
  Paraview */
  int numberOfCells_PARAVIEW() const;

  /** @brief Ecrit les points du convexe pour post-processing avec Paraview  
  @param f Flux de sortie 
  @param transform transformation courante 
  @param translation translation du centre de gravite */
  void write_polygonsPts_PARAVIEW( ostream &f, 
  	const Transform &transform, 
  	Vecteur const* translation = NULL ) const;
	
  /** @brief Renvoie les points du convexe pour post-processing avec Paraview  
  @param transform transformation courante 
  @param translation translation du centre de gravite */
  list<Point> get_polygonsPts_PARAVIEW( const Transform &transform,
  	Vecteur const* translation = NULL ) const;
  
  /** @brief Ecrit le convexe pour post-processing avec Paraview 
  @param connectivity connectivite des polygones Paraview
  @param offsets decalage dans la connectivite
  @param cellstype type de polygones Paraview
  @param firstpoint_globalnumber numero global du 1er point
  @param last_offset dernier offset utilise pour le convexe precedent */
  void write_polygonsStr_PARAVIEW( list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const; 
	
  /** @brief La forme est elle un cylindre ? */
  virtual bool isCylindre() const {return true;};	  
   //@}
  

protected:
  /** @name Parameters */
  //@{  
  Scalar m_radius; /**< le rayon du cylindre */  
  Scalar m_halfHeight; /**< la demi hauteur du cylindre */
  static int visuNodeNbOnPer; /**< nb de points sur le perimetre du cercle pour 
  	post-processing Paraview */  
  //@}


  /**@name Methods */
  //@{
  /** @brief Rayon circonscrit de la sphere
  @return Le rayon */
  Scalar BuildRayonRef() const;
  //@}   
};

#endif
