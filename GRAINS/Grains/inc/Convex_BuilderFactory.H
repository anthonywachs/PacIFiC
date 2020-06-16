#ifndef _Convex_BuilderFactory
#define _Convex_BuilderFactory

class Convex;

#include <iostream>
#include <string>
using namespace std;

#include "ReaderXML.hh"


/** @brief Fabrique pour des Convexes.

    Permet de contruire un convexe quelconque a partir des informations 
    de type et lecture de ses parametres

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation
    @author D. RAKOTONIRINA - IFP Energies nouvelles - Oct. 2014
    - Modificaion */
// ============================================================================
class Convex_BuilderFactory
{
public:
  /** @name Methods Static */
  //@{
  /** @brief Construction du convexe specifie a partir du noeud
  @param root Le noeud "Convex"
  @return Le nouveau Convex */
  static Convex* create(DOMNode *root);

  /** @brief Creation d'un convexe a partir du type indique 
  & lecture de ses donnees.
  @param type Nom du convexe a cree
  @param fileIn Flux de lecture
  @return Le convexe specifie */
  static Convex* create(string &type, istream &fileIn);
  //@}


private:
  /** @name Constructors */
  //@{
  /** @brief Classe Statique : pas de contructeur disponible. */
  Convex_BuilderFactory() {};

  /** @brief Classe Statique : pas de destructeur disponible. */
  ~Convex_BuilderFactory() {};
  //@}
};

#endif
