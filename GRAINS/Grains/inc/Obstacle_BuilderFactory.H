#ifndef _Obstacle_BuilderFactory
#define _Obstacle_builderFactory

#include "ReaderXML.hh"
#include <iostream>
using namespace std;

class Obstacle;


/** @brief Fabrique des obstacles.

    Construction des obstacles en fonction des options definies.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation */
// ============================================================================
class Obstacle_BuilderFactory
{
public:
  /**@name Methods Static */
  //@{
  /**
    @brief Creation d'un obstacle par decodage d'un noeud XML
    @param root Le noeud XML.
    @return L'obstacle cree.
    @details %Obstacle type can be :
      - Obstacle, then type Standard (MonObstacle) or Absorbant (ObstacleAbsorbant)
      - Periode (Obstacle_BuilderFactory::createPeriode)
      - Composite (CompObstacle)
  */
  static Obstacle* create( DOMNode *root );

  /** @brief Construction d'un obstacle par lecture sur fichier.
  Cette methode est utilisee pour reconstruire un objet persistant. */
  //static Obstacle* read(istream &fileIn);

  /** @brief Reload de l'Obstacle sous l'Obstacle referent
  @param tag Le type d'Obstacle à reloader
  @param obstacle L'Obstacle referent
  @param file Le fichier de persistance */
  static void reload( const string &tag, Obstacle& obstacle, 
		istream &file );
  //@}


private:
  /** @name Methods */
  //@{
  /** @brief Creation d'une periode d'obstacles.
  Une periode correspond a deux obstacles periodiques 
  regroupes dans un CompObstacle.
  @return La periode d'obstacle generee */
  static Obstacle* createPeriode( DOMNode *root );  
  //@}


  /**@name Constructors */
  //@{
  /** @brief Classe Statique : pas de contructeur disponible. */
  Obstacle_BuilderFactory() {};

  /** @brief Classe Statique : pas de destructeur disponible. */
  ~Obstacle_BuilderFactory() {};
  //@}
};

#endif
