/*
  SOLID - Software Library for Interference Detection
  Copyright (C) 1997-1998  Gino van den Bergen

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Library General Public
  License as published by the Free Software Foundation; either
  version 2 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Library General Public License for more details.

  You should have received a copy of the GNU Library General Public
  License along with this library; if not, write to the Free
  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  Please send remarks, questions and bug reports to gino@win.tue.nl,
  or write to:
  Gino van den Bergen
  Department of Mathematics and Computing Science
  Eindhoven University of Technology
  P.O. Box 513, 5600 MB Eindhoven, The Netherlands
*/
// ============================================================================
#ifndef _BBOX_H_
#define _BBOX_H_

#include "Vecteur.H"
#include "Point.H"
using namespace solid;


/** Boite englobant des formes geometriques.

Une BBox est un pave ou une "boite" permettant de contenir 
des formes a geometrie complexe. Cet objet peut donc s'averer 
necessaire pour optimiser certains algorithmes.
    
@author 
Gino van den Bergen - Eindhoven University of Technology - Creation */
// ==========================================================================
class BBox 
{
public:
  /**@name Constructeurs*/
  //@{
  /** @brief Constructeur par defaut */
  BBox(); 
  
  /** @brief Constructeur avec initialisation par deux points 
  @param min le point dont les composantes sont minimum
  @param max le point dont les composantes sont maximum */
  BBox(const Point& min, const Point& max);

  /** @brief Constructeur par copie. 
  @param bbox_ La boite de reference */
  BBox(const BBox &bbox_);
  
  /** @brief Destructeur */
  ~BBox();
  //@}

  /**@name Methods */
  //@{
  /** @brief Engendre la boite intersection des deux autres boites : a et b.
  Si les deux boites ne sont pas en intersection, 
  la boite d'origine n'est pas modifiee. 
  @param a boite 1 
  @param b boite 2 */
  void closest(const BBox &a, const BBox &b);

  /** @brief Engendre la boite contenant deux autres boites
  @param a boite 1
  @param b boite 2 */
  void enclose(const BBox& a, const BBox& b); 

  /** @brief Renvoi le centre de la boite 
  @return le centre de la boite*/
  const Point& getCenter() const;

  /** @brief Renvoi l'etendue de la boite a partir de son centre 
  @return le vecteur representant l'etendue de la boite */
  const Vecteur& getExtent() const;

  /** @brief Renvoi la coordonnee numero i minimum de la boite
  @return la coordonnee minimum
  @param i indice de la composante */
  Scalar getLower(int i) const;

  /** @brief Renvoi la coordonnee numero i maximum de la boite
  @return la coordonnee maximum
  @param i indice de la composante */
  Scalar getUpper(int i) const;

  /** @brief Construit une boite contenant la boite courante et un point
  @param p le point a inclure */
  void include (const Point& p);
  
  /** @brief Construit une boite contenant la boite courante et une autre boite
  @param b une boite */
  void include (const BBox& b);

  /** @brief Determine si une cellule est dans la zone d'influence de la boite
  @return Vrai si la cellule est dans la zone d'influence
  @param p Le centre de la cellule
  @param arete la demi-taille de la cellule */
  bool InZone(Point const* p, const Scalar &arete) const;

  /** @brief Determine si une cellule est dans la zone d'influence de la boite
  @return Vrai si la cellule est dans la zone d'influence
  @param p Le centre de la cellule
  @param arete_X la demi-taille de la cellule dans la direction X 
  @param arete_Y la demi-taille de la cellule dans la direction Y   
  @param arete_Z la demi-taille de la cellule dans la direction Z */
  bool InZone(Point const* p, const Scalar &arete_X, const Scalar &arete_Y,
  	const Scalar &arete_Z) const;

  /** @brief Renvoi l'axe des coordonnees ou la boite est le plus etendue
  @return le numero de l'axe des coordonnees */
  int longestAxis() const;

  /** @brief Positionne le centre de la boite 
  @param p le point ou va se positionner la boite */
  void setCenter(const Point& p);

  /** @brief Construit une boite infini */
  void setEmpty(); 

  /** @brief Change l'etendue de la boite
  @param v le vecteur representant la nouvelle etendue */
  void setExtent(const Vecteur& v);

  /** @brief Modifie la boite a l'aide de deux points
  @param min le point dont les composantes sont minimum
  @param max le point dont les composantes sont maximum */
  void setValue(const Point& min, const Point& max);

  /** @brief Renvoi la taille de la boite. 
  En fait, la demi plus grande arete de la boite.
  @return Taille de la boite */
  Scalar size() const ;
  
  /** @brief Debug
  @param s message de debug */
  void debug(char *s) const;
  //@}
  

  /** @name Methods friend */
  //@{
  /** @brief Determine si deux boites sont en contact ou non
  @return Vrai si les boites sont en contact
  @param a boite 1
  @param b boite 2 */
  friend bool intersect(const BBox& a, const BBox& b);
  
  /** @brief Operateur << 
  @param f flux de sortie 
  @param B la boite */
  friend ostream& operator << (ostream &f, const BBox &B);   
  //@}


  /**@name Operators */
  //@{
  /** @brief Operateur d'affectation 
  @param rhs boite à affecter */
  BBox& operator= (const BBox& rhs);
  //@}


private:
  /** @name Parameters */
  //@{  
  Point center; /**< Centre de la boite */
  Vecteur extent; /**< Etendue de la boite a partir de son centre */
  //@}
};

#endif


