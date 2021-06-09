#ifndef _Cinematique_BuilderFactory
#define _Cinematique_BuilderFactory

#include <string>
#include <iostream>
using namespace std;

class CineParticule;
class Convex;


/** @brief Fabrique des Cinematiques.

    Permet la construction d'une cinematique a partir des informations de
    fichier de mise en donnees.
    En version 2D -> LeapFrog_2D, 3D -> TR_LeapFrog.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation */
// ============================================================================
class Cinematique_BuilderFactory
{
public:
  /**@name Methods Static */
  //@{
  /** @brief Creation de la cinematique de la particule
  !! IMPORTANT : si le pointeur convexe est nul, cela signifie que la particule
  correspondante est une particule composite !! (convention
  @return La cinematique associee a la particule 
  @param convexe type de convexe */
  static CineParticule* create( const Convex* convexe );

  /** @brief Construction de la cinematique a partir de son enregistrement.
  @return La cinematique construite.
  @param fileIn Flux de lecture. 
  @param convexe type de convexe */
  static CineParticule* read( istream &fileIn, const Convex* convexe );
  //@}


private:
  /**@name Constructors */
  //@{
  /** @brief Classe Statique : pas de contructeur disponible. */
  Cinematique_BuilderFactory() {};

  /** @brief Classe Statique : pas de destructeur disponible. */
  ~Cinematique_BuilderFactory() {};
  //@}
};

#endif
    
