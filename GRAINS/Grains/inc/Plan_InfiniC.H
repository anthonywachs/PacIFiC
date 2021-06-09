/* D.PETIT - Creation
   
   Complement de l'algorithme GJK - Engine
   Nouvelle classe de convexe : Le Plan infini .
*/

#ifndef _PLAN_INFINIC_
#define _PLAN_INFINIC_

#include "Plan_FiniC.H"


/** Gestion des convexes de la forme plane (et infinie).
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Plan_InfiniC : public Plan_FiniC 
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  Plan_InfiniC();

  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant */
  Plan_InfiniC(istream &fileIn);

  /** @brief Destructeur */
  ~Plan_InfiniC();
  //@}
   /** @brief get shrinking choice */
  int getShrinkingChoice()const;


  /** @name Methods */
  //@{
  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return PLAN;};

  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  void printClass(ostream &fileOut) const;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  void readClass(istream &fileIn);
  //@}
  

private:
  /**@name Methods */
  //@{
  /** @brief Rayon circonscrit du plan
  @return Le rayon */
  Scalar BuildRayonRef() const;
  //@}   
};

#endif
