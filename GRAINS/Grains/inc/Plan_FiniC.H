/* D.PETIT - Creation
   
   Complement de l'algorithme GJK - Engine
   Nouvelle classe de convexe : Le Plan fini .
*/

#ifndef _PLAN_FINIC_
#define _PLAN_FINIC_

#include "Convex.H"


/** Gestion des convexes de la forme plane (et finie).
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Plan_FiniC : public Convex 
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur avec initialisation
  @param x la longueur du plan suivant l'axe x
  @param y la longueur du plan suivant l'axe y */
  Plan_FiniC(const Scalar x = 0.0, 
	     const Scalar y = 0.0);

  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant */
  Plan_FiniC(istream &fileIn);

  /** @brief Destructeur */
  ~Plan_FiniC();
  //@}

  
  /** @name Methods */
  //@{
  /** @brief Determine l'inertie d'un plan fini
  @param inertie L'inertie du convexe
  @param inertie_1 L'inertie inverse du convexe*/
  bool BuildInertie(Scalar *inertie,
	Scalar *inertie_1) const;
  
  /** @brief gets the shrinking choice*/
  int getShrinkingChoice()const;  
  
  /** @brief set the Shrinking radius*/
  void setShrinkingRadius(Scalar CurrentRadius);
 
  /** @brief Clone d'un plan fini
  @return le clone d'un plan fini */
  Convex* clone() const;

  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return PLAN;};

  /** @brief Renvoi la longueur du plan suivant l'axe x
  @return la longueur */
  Scalar getProjx() const;

  /** @brief Renvoi la longueur du plan suivant l'axe y
  @return la longueur */
  Scalar getProjy() const;

  /** @brief Determine le volume d'un plan fini 
  @return le volume */
  Scalar getVolume() const;

  /** @brief Fonction support : Determine le point du support convexe de la 
  forme se situant dans une direction et un sens donne. 
  Pour le Plan fini : on prend le sommet le plus proche du vecteur
  @return le point sur l'enveloppe convexe obtenu
  @param v la direction et le sens */
  Point support(const Vecteur& v) const;
  //@}

  
  /**@name Methods virtual */
  //@{
  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  virtual void printClass(ostream &fileOut) const;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  virtual void readClass(istream &fileIn);
  //@}


protected:
  /** @name Parameters */
  //@{
  Scalar projx; /**< la longueur du plan suivant l'axe x */
  Scalar projy; /**< la longueur du plan suivant l'axe y */
  //@}


  /**@name Methods */
  //@{
  /** @brief Rayon circonscrit du plan
  @return Le rayon */
  Scalar BuildRayonRef() const;
  //@}   
};

#endif
