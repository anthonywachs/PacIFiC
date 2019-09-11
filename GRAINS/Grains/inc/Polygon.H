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
#ifndef _POLYGON_H_
#define _POLYGON_H_

#include "Polytope.H"
#include "Basic.H"
#include "IndexArray.H"
#include "ReaderXML.hh"
#include <string>
using namespace std;


/** @brief Gestion des convexes de la forme d'un polytope : Le polygone (objet 2D).
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Polygon : public Polytope 
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant
  @param nb_point le nombre points constituant le polytope
  @param ref la reference du tableau de points 
  @param ia la reference du tableau d'indices */
  Polygon( istream &fileIn, int nb_point, VertexBase& ref,
  	 IndexArray& ia );

  /** @brief Constructeur par copie 
  @param copie le polygone a copier */
  Polygon( const Polygon &copie );

  /** @brief sets the Shrinking radius 
  @param Shrinking radius */
  void setShrinkingRadius(Scalar CurrentRadius);   
  
  /** @brief gets the shrinking choice 
  @param shrinking*/
  int getShrinkingChoice()const;
   
  /** @brief Destructeur */
  ~Polygon();
  //@}
  

  /** @name Methods */
  //@{
  /** @brief Determine l'inertie d'un polygone
  @param inertie L'inertie du convexe
  @param inertie_1 L'inertie inverse du convexe*/
  bool BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const;

  /** @brief Clone d'un polygone
  @return le clone d'un polygone */
  Convex* clone() const;

  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return POLYGON;};  

  /** @brief Determine le volume d'un polygone
  @return le volume*/
  Scalar getVolume() const;

  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  void printClass( ostream &fileOut ) const;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  void readClass( istream &fileIn );

  /** @brief Fonction support : Determine le point du support convexe de la 
  forme se situant dans une direction et un sens donne. 
  Pour le polygone : le sommet ayant la plus grande projection sur le 
  vecteur v (ie. le sommet le plus 'proche' de v)
  @return le point sur l'enveloppe convexe obtenu
  @param v la direction et le sens */
  Point support( const Vecteur& v ) const;
  
  /** @brief Nombre de polygones elementaires pour post-processing avec 
  Paraview */
  int numberOfCells_PARAVIEW() const;

  /** @brief Ecrit le convexe pour post-processing avec Paraview 
  @param connectivity connectivite des polygones Paraview
  @param offsets decalage dans la connectivite
  @param cellstype type de polygones Paraview
  @param firstpoint_globalnumber numero global du 1er point
  @param last_offset dernier offset utilise pour le convexe precedent */
  void write_polygonsStr_PARAVIEW( list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;
  //@}


  /** @name Methods static */
  //@{
  /** @brief Creation d'un polygone a partir d'un flux entrant
  @return le polygone cree
  @param fileIn le flux entrant */
  static Polygon* create( istream &fileIn );

  /** @brief Creation d'un polygone
  @param root Le noeud "<Polygon>"
  @return Le nouveau Polygon */
  static Polygon* create( DOMNode* root );
  //@}
  

private:
  /** @name Parameters */
  //@{
  IndexArray *m_cobound; /**< Description des points sommets. Noeuds relies au 
  	sommet courant */
  mutable unsigned int m_curr_vertex; /**< indice du dernier point obtenu par 
  	la fonction support */
  Scalar *m_InertiePoly; /**< Inertie du polygone */
  Scalar m_surface; /**< volume du polygone */
  //@}


  /**@name Constructor */
  //@{
  /** @brief Constructeur par defaut interdit */
  Polygon();
  //@}


  /**@name Methods */
  //@{
  /** @brief Constructions du polygone , de son volume (en réalité sa surface) 
  et de son inertie a partir de ses sommets (Si qhull n'est pas interface)
  @param nbedge le nombre d'arête du polygone
  @param edge represente un pointeur sur des tableaux representant les 
  indices des sommets de chaque arête */
  void BuildPolygon( const int nbedge, const IndexArray *edge );
  
  /** @brief Calcul la surface et l'inertie d'un polygon. 
  @param P1 point 1 (signification ?)
  @param P2 point 2 (signification ?) */
  void CalculSurfaceInertie( const Point &P1, const Point &P2 );
    
  /** @brief Lecture des sommets de chaque face a partir d'un flux entrant
  (Si qhull n'est pas interface).
  @param fileIn le flux entrant*/ 
  void readface( istream &fileIn );

  /** @brief Rayon circonscrit du polygone
  @return Le rayon */
  Scalar BuildRayonRef() const;
  
  /** @brief Remise a zero du volume et de l'inertie*/
  void Initialisation();
  //@}   
};

#endif
