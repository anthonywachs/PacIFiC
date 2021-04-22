/*
  GJK Engine - A Fast and Robust GJK Implementation 
  Copyright (C) 1998  Gino van den Bergen

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

#ifndef _CONE_H_
#define _CONE_H_

#include "Convex.H"


/** Gestion des convexes de la forme coniques.
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Cone : public Convex {
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur avec initialisation
  @param r le rayon de la base du cone
  @param h la hauteur du cone */
  Cone(Scalar r = 0, Scalar h = 0); 
  
  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant */
  Cone(istream &fileIn);
  
  /** @brief Destructeur */
  ~Cone();
  //@}
  

  /** @name Methods */
  //@{
  /** @brief Determine l'inertie d'un cone
  @param inertie L'inertie du convexe
  @param inertie_1 L'inertie inverse du convexe*/
  bool BuildInertie(Scalar *inertie,Scalar *inertie_1) const;
  
 /** @name Methods */
  //@{
  /** @ get shrinking choice*/
  int getShrinkingChoice()const;

 /** @name Methods */
  //@{
  /** @ set the Shrinking radius*/
  void setShrinkingRadius(Scalar CurrentRadius); 

  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return CONE;};

  /** @brief Fonction support : Determine le point du support convexe de la 
  forme se situant dans une direction et un sens donne
  Pour le cone : ...
  @return le point sur l'enveloppe convexe obtenu
  @param v la direction et le sens */
  Point support(const Vecteur& v) const;

  /** @brief Clone d'un cone
  @return le clone d'un cone */
  Convex* clone() const;

  /** @brief Determine le volume d'un cone
  @return le volume*/
  Scalar getVolume() const;

  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  void printClass(ostream &fileOut) const;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  void readClass(istream &fileIn);
   //@}
  

protected:
  /**@name Methods */
  //@{
  /** @brief Rayon circonscrit de la sphere
  @return Le rayon */
  Scalar BuildRayonRef() const;
  //@} 


  /** @name Parameters */
  //@{
  Scalar bottomRadius; /**< le rayon de la base du cone */
  Scalar quaterHeight; /**< le quart de la hauteur du cone */
  Scalar sinAngle; /**< le sinus du demi angle au centre */
  //@}  
};

#endif
