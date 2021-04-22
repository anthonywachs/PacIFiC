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
#ifndef _POLYTOPE_H_
#define _POLYTOPE_H_

#include "Convex.H"
#include "IndexArray.H"
#include "VertexBase.H"
#include <string>
using namespace std;


/** Gestion des convexes de la forme polyedrique (Simplex, Polygone, Polyedre).
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Polytope : public Convex 
{
public:
  /** @name Methods Virtual */
  //@{
  /** @brief Determine l'inertie d'un polytope
  @param inertie L'inertie du convexe
  @param inertie_1 L'inertie inverse du convexe*/
  virtual bool BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const = 0;

  /** @brief Clone d'un polytope
  @return le clone d'un polytope */
  virtual Convex* clone() const = 0;

  /** @brief Determine le volume d'un polytope
  @return le volume*/
  virtual Scalar getVolume() const = 0;

  /** @brief Fonction support : Determine le point du support convexe de la 
  forme se situant dans une direction et un sens donne. 
  Pour le polytope : le sommet ayant la plus grande projection sur le 
  vecteur v (ie. le sommet le plus 'proche' de v)
  @return le point sur l'enveloppe convexe obtenu
  @param v la direction et le sens */
  virtual Point support( const Vecteur& v ) const = 0;
  //@}
  

  /**@name Methods */
  //@{
  /** @brief Renvoi le nombre de sommets du polytope
  @return le nombre de sommets */
  int numVerts() const;

  /** @brief Points decrivant l'enveloppe
  @return Un vecteur de points contenant les neouds sommets. */
  vector<Point> getEnveloppe() const;
  
  /** @brief Renvoie le nombre de sommets ou un code équivalent */
  int getNbCorners() const;  
  
  /** @brief Nombre de points pour post-processing avec Paraview */
  virtual int numberOfPoints_PARAVIEW() const;
  
  /** @brief Ecrit les points du convexe pour post-processing avec Paraview  
  @param f Flux de sortie 
  @param transform transformation courante 
  @param translation translation du centre de gravite */
  void write_polygonsPts_PARAVIEW( ostream &f, 
  	const Transform &transform, 
  	Vecteur const* translation = NULL ) const;
	
  /** @brief Renvoie les points du convexe pour post-processing avec Paraview  
  @param transform transformation courante 
  @param translation translation du centre de gravite */
  list<Point> get_polygonsPts_PARAVIEW( const Transform &transform,
  	Vecteur const* translation = NULL ) const; 
  //@}

  
  /**@name Operators */
  //@{
  /** @brief Renvoi le i sommet du polytope
  @return le point etant le i sommet du polytope
  @param i l'indices du sommet a retourner */
  const Point& operator[] (int i) const;

  /** @brief Renvoi le i sommet du polytope
  @return le point etant le i sommet du polytope
  @param i l'indices du sommet a retourner */
  Point& operator[] (int i);
  //@}


protected:
  /**@name Parameters */
  //@{
  const VertexBase& m_base; /**< Alias sur le tableau de points ou se trouvent 
  	les sommets du polytope */
  const IndexArray& m_index; /**< Tableau d'indices permettant de situer les sommets 
      dans le tableau de points. */
  string m_fichPoly; /**< Fichier contenant les points et les faces du 
  	simplex */
  //@}


  /**@name Constructors */
  //@{
  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant
  @param nb_point le nombre points constituant le polytope
  @param ref la reference du tableau de points 
  @param ia la reference du tableau d'indices */
  Polytope( istream &fileIn, int nb_point, VertexBase& ref,
  	 IndexArray& ia );

  /** @brief Constructeur par copie 
  @param copie le polytope a copier */
  Polytope( const Polytope& copie );
 
  /** @brief Destructeur */
  virtual ~Polytope();
  //@}

  
private:
  /**@name Constructor */
  //@{
  /** @brief Constructeur par defaut interdit */
  Polytope();
  //@}
  
};

#endif
