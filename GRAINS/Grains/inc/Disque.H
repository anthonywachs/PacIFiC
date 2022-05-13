#ifndef _DISQUE_H_
#define _DISQUE_H_

#include "Convex.H"

#include "ReaderXML.hh"


/** @brief Gestion des convexes de forme disque (forme 2D).
    
    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation */
// ============================================================================
class Disque : public Convex 
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur avec initialisation
  @param r le rayon du disque */
  Disque( Scalar r = 0 );

  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant */
  Disque( istream &fileIn );

  /** @brief Constructeur
  @param root Le noeud "<Disque>" */
  Disque( DOMNode* root );

  /** @brief Destructeur */
  ~Disque();
  //@}
  

  /** @name Methods */
  //@{
  /** @brief Determine l'inertie d'une sphere 
  @param inertie L'inertie de la sphere
  @param inertie_1 L'inertie inverse de la sphere */
  bool BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const;

  /** @brief Clone d'une sphere
  @return Le clone d'une sphere */
  Convex* clone() const;

  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return DISQUE2D;};

  /** @brief Points decrivant l'enveloppe du Disque.
  Pour un disque pas de prise en compte de points sommets -> 0
  @return Un vecteur de points sur l'enveloppe du disque. */
  vector<Point> getEnveloppe() const;

  /** @brief Renvoie le nombre de sommets ou un code équivalent */
  int getNbCorners() const;  

  /** @brief Determine le volume d'un disque.
  @return Le volume (en fait pour un disque, la surface) */
  Scalar getVolume() const;

  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  void printClass( ostream &fileOut ) const;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  void readClass( istream &fileIn );

  /** @brief Determine le point du support convexe de la forme se
  situant dans une direction et un sens donne.
  Pour la sphere : on prend l'intersection de la sphere 
  et de la demi droite dirigee par le vecteur et d'origine 
  l'origne de la sphere
  @return le point sur l'enveloppe convexe obtenu.
  @param v la direction et le sens. */
  Point support( const Vecteur& v ) const;
  
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

  /** @brief Nombre de points sur le perimetre pour post-processing
  @param nbpts nombre de points sur le perimetre */
  static void SetvisuNodeNb( int nbpts ) { visuNodeNb = nbpts; };
  //@}
  

protected:
  /**@name Methods */
  //@{
  /** @brief Rayon circonscrit de la sphere
  @return Le rayon */
  Scalar BuildRayonRef() const;
  //@} 

/** @brief get shrinking choice
  @ return Shrinking */
  int getShrinkingChoice()const;

/** @brief set current radius
  */
  void setShrinkingRadius(Scalar CurrentRadius=0);
  
  static int visuNodeNb; /**< nb de points sur le perimetre du cercle pour 
  	post-processing Paraview */
  //@}  



  Scalar radius ;/**< radius of the disk */
  int Shrinking;/**< the shrinking choice to be read from xml*/
  Scalar initial_radius;/**< the initial radius for shrinking particles */

};

#endif
