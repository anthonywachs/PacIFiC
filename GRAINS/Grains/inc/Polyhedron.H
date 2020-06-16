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
#ifndef _POLYHEDRON_H_
#define _POLYHEDRON_H_

#include "Polytope.H"
#include "IndexArray.H"
#include "Basic.H"
#include "ReaderXML.hh"


/** @brief Gestion des convexes de la forme d'un polytope : Le Polyedre.

    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Polyhedron : public Polytope 
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant
  @param nb_point le nombre points constituant le polytope
  @param ref la reference du tableau de points 
  @param ia la reference du tableau d'indices */
  Polyhedron( istream &fileIn, int nb_point, VertexBase& ref,
  	 IndexArray& ia );

  /** @brief Constructeur par copie 
  @param polyedre le Polyedre a copier */
  Polyhedron( const Polyhedron& polyedre );

  /** @brief Destructeur */
  ~Polyhedron();
  //@}
  
  /** @brief gets the shrinking choice 
  @param shrink*/ 
  int getShrinkingChoice()const;
   
  /** @brief sets the Shrinking radius 
  @param Shrinking_radius */
  void setShrinkingRadius(Scalar CurrentRadius); 
  
 /** @name Methods */
  //@{
  /** @brief Determine l'inertie d'un polyedre
  @param inertie L'inertie du convexe
  @param inertie_1 L'inertie inverse du convexe*/
  bool BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const;

  /** @brief Constructions du polyedre , de son volume et de son inertie a 
  partir de ses sommets (Si qhull n'est pas interface)
  @param nbface le nombre de face du polyedre
  @param face represente un pointeur sur des tableaux representant les i
  ndices des sommets de chaque face */
  void BuildPolyhedron( const int nbface, const IndexArray *face );
  
  /** @brief Lecture des sommets de chaque face a partir d'un flux entrant
  (Si qhull n'est pas interface).
  @param fileIn le flux entrant*/ 
  void readface(istream &fileIn);
  
  /** @brief Fonction support : Determine le point du support convexe de la 
  forme se situant dans une direction et un sens donne. 
  Accessible si la description des faces n'est pas indiquee (original
  SOLID 2.0 Software).
  Pour le polyedre : le sommet ayant la plus grande projection sur le 
  vecteur v (ie le sommet le plus 'proche' de v)
  @return le point sur l'enveloppe convexe obtenu
  @param v la direction et le sens */
  Point support20( const Vecteur& v ) const;

  /** @brief Clone d'un polyedre
  @return le clone d'un polyedre */
  Convex* clone() const;

  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return POLYHEDRON;};

  /** @brief Description des faces
  @return La description des faces via les indices des points sommets */
  vector<vector<int> > const* getFaces() const;

  /** @brief Determine le volume d'un polyedre
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
  Pour le polyedre : le sommet ayant la plus grande projection sur le 
  vecteur v (ie le sommet le plus 'proche' de v)
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
	
  /** @brief Ecrit le convexe au format STL pour lien avec openFoam 
  @param f Flux de sortie 
  @param transform transformation courante */
  void write_convex_STL( ostream &f, const Transform &transform ) const;	
  //@}
  

  /** @name Methods statics*/
  //@{
  /** @brief Creation d'un polyedre a partir d'un flux entrant
  @return le polyedre cree
  @param fileIn le flux entrant */  
  static Polyhedron* create( istream &fileIn );

  /** @brief Construction d'un Polyhedron
  @param root Le noeud "<Polyhedron>" */
  static Polyhedron* create( DOMNode* root );
  //@}


private:
  /** @name Parameters */
  //@{
  IndexArray *m_cobound; /**< Description des points sommets. Noeuds relies au 
  	sommet courant */
  mutable unsigned int m_curr_vertex; /**< indice du dernier point obtenu par la
  	fonction support */
  Scalar *m_InertiePoly; /**< Inertie du polyedre */
  Scalar m_VolumePoly; /**< volume du polyedre */
  vector< vector<int> >* m_allFaces; /**< description des faces par les numéros 
  	des sommets */ 
  //@}


  /** @name Methods */
  //@{
  /** @brief Rayon circonscrit du polyèdre
  @return Le rayon */
  Scalar BuildRayonRef() const;
  
  /** @brief Determine le volume et de l'inertie d'un polyedre
  @param H la projection de l'origine du polyedre sur une face
  @param P1 un point de la face ou se trouve H
  @param P2 le point voisin de P1 (en parcourant la face suivant son 
      orientation) sur la face ou se trouve H. */
  void CalculVolumeInertie( const Point &H, const Point &P1, const Point &P2 );

  /** @brief Remise a zero du volume et de l'inertie*/
  void Initialisation();
  //@}


  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut interdit */
  Polyhedron();
  //@}
};

#endif
