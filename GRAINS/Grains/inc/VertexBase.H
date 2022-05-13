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

#ifndef _VERTEXBASE_H_
#define _VERTEXBASE_H_

#include "Point.H"
using namespace solid;

/** @brief Gestion d'un ensemble de points.

    Cette classe permet de pointer sur un ensemble de points
 @author Gino van den Bergen -  Eindhoven University of Technology - Creation*/
// =========================================================================== 
class VertexBase 
{
public:
  /** @name Constructeurs */
  //@{
  /** @brief Constructeur par defaut. */
  VertexBase();
  
  /** @brief Constructeur avec initialisation du pointeur 
  @param ptr un pointeur sur un tableau de point */
  VertexBase(const void *ptr);

  /** @brief Destructeur */
  ~VertexBase();
  //@}

  
  /** @name Methods */
  //@{
  /** @brief Renvoi l'adresse du pointeur 
  @return l'adresse du pointeur*/
  const void *getPointer() const;

  /** @brief Renvoi le point d'indice i dans le tableau de point
  @return le point a l'indice i
  @param i l'indice */
  Point& operator[](int i) const ;
  //@}
  
private:
  /** @name Parameters */
  //@{
  const void *base; /**< le pointeur sur la zone memoire contenant les points */
  //@}
};


// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
inline const void *VertexBase::getPointer() const { return base; }

// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
inline Point& VertexBase::operator[](int i) const 
{
  return ((Point *)base)[i]; 
}

#endif
