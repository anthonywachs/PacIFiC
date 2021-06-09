/* D.PETIT - Creation
   
   Complement de l'algorithme GJK - Engine
   Nouvelle classe de convexe : Le point.
*/
#ifndef _POINTC_H_
#define _POINTC_H_

#include "Convex.H"


/** Gestion des convexes de la forme d'un point.
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class PointC : public Convex 
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  PointC();

  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant */
  PointC(istream &fileIn);

  /** @brief Destructeur */
  ~PointC();
  //@}

  
  /** @name Methods */
  //@{
  /** @brief Determine l'inertie d'un point 
  @param inertie L'inertie du convexe
  @param inertie_1 L'inertie inverse du convexe*/
  bool BuildInertie(Scalar *inertie,Scalar *inertie_1) const;
  
   /** @brief gets the shrinking choice*/ 
  int getShrinkingChoice()const;

   /** @brief sets the Shrinking radius */
  void setShrinkingRadius(Scalar CurrentRadius); 

  /** @brief Clone d'un point
  @return le clone d'un point */
  Convex* clone() const;

  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return POINT;};

  /** @brief Determine le volume d'un point
  @return le volume*/
  Scalar getVolume() const;

  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  void printClass(ostream &fileOut) const;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  void readClass(istream &fileIn);

  /** @brief Fonction support : Determine le point du support convexe de la 
  forme se situant dans une direction et un sens donne. 
  Pour le point : le point lui-meme
  @return le point sur l'enveloppe convexe obtenu
  @param v la direction et le sens */
  Point support(const Vecteur& v) const;
  //@}


private:
  /**@name Methods */
  //@{
  /** @brief Rayon circonscrit du point
  @return Le rayon */
  Scalar BuildRayonRef() const;
  //@}   
};

#endif
