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

#ifndef _INDEXARRAY_H_
#define _INDEXARRAY_H_

#include "Basic.H"

/** Gestion d'un tableau d'indices
    @author 
    Gino van den Bergen -  Eindhoven University of Technology - Creation */
// ============================================================================
class IndexArray 
{
public:
  /** @name Constructeurs */
  //@{
  /** @brief Constructeur par defaut */
  IndexArray();

  /** @brief Constructeur avec initialisation du nombre d'element
  @param n le nombre d'elements */
  IndexArray(int n);

  /** @brief Constructeur avec initialisation du nombre d'elements 
  et du contenu du tableau
  @param n le nombre de d'elements
  @param v tableau a copie */
  IndexArray(int n, const unsigned int v[]);

  /** @brief Destructeur */
  ~IndexArray();
  //@}
  

  /**@name Methods */
  //@{
  /** @brief Renvoi l'adresse du tableau 
  @return l'adresse du tableau d'entier */
  unsigned int* getAdress() const;

  /** @brief Renvoi la taille du tableau 
  @return la taille du tableau*/
  int size() const;
  //@}


  /**@name Operators */
  //@{
  /** @brief Determine l'indice contenu dans le tableau a l'indice i
  @return l'indice contenu dans le tableau
  @param i recherche dans le tableau a l'indice i*/
  int operator[](int i) const ;
  //@}
  

private:
  /**@name Constructors */
  //@{
  /** @brief Constructeur de copie interdit 
  @param Ia objet a copier */
  IndexArray(const IndexArray &Ia);
  //@}


  /**@name Operators */
  //@{
  /** Operateur d'egalite interdit  
  @param Ia objet a copier */
  IndexArray& operator=(const IndexArray &Ia);
  //@}

  
  /** @name Parameters */
  //@{
  unsigned int *indices; /**< Pointeur sur le tableau d'indices */
  int count; /**< nombre d'elements du tableau */
  //@}
};
  
// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
inline int IndexArray::operator[](int i) const { return indices[i]; }

// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
inline unsigned int* IndexArray::getAdress() const { return indices; }

// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
inline int IndexArray::size() const { return count; }

#endif

