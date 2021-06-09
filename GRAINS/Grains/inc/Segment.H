/* D.PETIT - Creation
   
   Complement de l'algorithme GJK - Engine
   Nouvelle classe de convexe : Le segment .
*/
#ifndef _SEGMENT_H_
#define _SEGMENT_H_

#include "Convex.H"


/** Gestion des convexes de la forme d'un segment.
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Segment : public Convex 
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur avec initialisation
  @param x la longueur du segment */
  Segment( Scalar x = 0 );

  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant */
  Segment( istream &fileIn );

  /** @brief Destructeur */
  ~Segment() {};
  //@}


  /** @name Methods */
  //@{
  /** @brief Determine l'inertie d'un segment
  @param inertie L'inertie du convexe
  @param inertie_1 L'inertie inverse du convexe*/
  bool BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const;

  /** @brief gets the shrinking choice 
  @param shrink */
  int getShrinkingChoice()const;
   

  /** @brief gets the Shrinking radius 
  @param Shrinking_radius */
  void setShrinkingRadius(Scalar CurrentRadius); 
  
  /** @brief Clone d'un Segment
  @return le clone d'une sphere */
  Convex* clone() const;

  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return SEGMENT;};

  /** @brief Renvoi la longueur du segment a partir de son centre de gravite
  @return la longueur */
  Scalar getLength() const;

  /** @brief Determine le volume d'un segment 
  @return le volume*/
  Scalar getVolume() const;

  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  void printClass( ostream &fileOut ) const;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  void readClass( istream &fileIn );

  /** @brief Fonction support : Determine le point du support convexe de la 
  forme se situant dans une direction et un sens donne. 
  Pour le segment : on prend le sommet le plus 'proche' de la 'pointe' 
  du vecteur
  @return le point sur l'enveloppe convexe obtenu
  @param v la direction et le sens */
  Point support( const Vecteur& v ) const;
  
  /** @brief Renvoi la tranformation pour un segment associé au vecteur v et au
  point gc
  @param v vecteur direction
  @param gc centre de gravité du segment */
  static Transform computeTransform( const Vecteur& v, const Point &gc );

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
  	Vecteur const* translation = NULL )
  	const;
  
  /** @brief Ecrit le convexe pour post-processing avec Paraview 
  @param connectivity connectivite des polygones Paraview
  @param offsets decalage dans la connectivite
  @param cellstype type de polygones Paraview
  @param firstpoint_globalnumber numero global du 1er point
  @param last_offset dernier offset utilise pour le convexe precedent */
  void write_polygonsStr_PARAVIEW( list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;
  //@}


protected:
  /** @name Parameters */
  //@{
  Scalar halflength; /**< la demi longueur du segment */
  //@}


  /**@name Methods */
  //@{
  /** @brief Rayon circonscrit du segment
  @return Le rayon */
  Scalar BuildRayonRef() const;
  //@}   
};

#endif
