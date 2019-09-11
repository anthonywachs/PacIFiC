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
// ============================================================================
#ifndef _SPHERE_H_
#define _SPHERE_H_

#include "Convex.H"

#include "Point.H"
using namespace solid;


#include "ReaderXML.hh"


/** @brief Gestion des convexes de forme spherique.

    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Sphere : public Convex 
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur avec initialisation
  @param r le rayon de la sphere */
  Sphere( Scalar r = 0 );

  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant */
  Sphere( istream &fileIn );

  /** @brief Constructeur
  @param root Le noeud "<Sphere>" */
  Sphere( DOMNode* root );

  /** @brief Destructeur */
  ~Sphere();
  //@}
  

  /** @name Methods */
  //@{
  /** @brief Determine l'inertie d'une sphere 
  @param inertie L'inertie de la sphere
  @param inertie_1 L'inertie inverse de la sphere */
  bool BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const;


  /** @brief gets the shrinking choice 
  @param shrink */
  int getShrinkingChoice()const;

  /** @brief gets the Shrinking radius 
  @param Shrinking radius*/
  void setShrinkingRadius(Scalar CurrentRadius=0); 

  /** @brief Clone d'une sphere
  @return Le clone d'une sphere */
  Convex* clone() const;

  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return SPHERE;};

  /** @brief Points decrivant l'enveloppe
  @return La liste des points sommets */
  vector<Point> getEnveloppe() const;

  /** @brief Description des faces
  @return La description des faces via les indices des points sommets */
  vector<vector<int> > const* getFaces() const;

  /** @brief Renvoie le nombre de sommets ou un code équivalent */
  int getNbCorners() const;  

  /** @brief Determine le volume d'une sphere.
  @return le volume */
  Scalar getVolume() const;

  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  void printClass( ostream &fileOut ) const;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  void readClass( istream &fileIn );

  /** @brief Determine le point du support convexe de la forme se
  situant dans une direction et un sens donne.
  Pour la sphere : on prend l'intersection de la sphere 
  et de la demi droite dirigee par le vecteur et d'origine 
  l'origne de la sphere
  @return le point sur l'enveloppe convexe obtenu.
  @param v la direction et le sens. */
  Point support( const Vecteur& v ) const;
  
  /** @brief Ecriture de l'objet dans sa configuration courante pour
  post-processing avec GMV
  @param fileOut flux de sortie
  @param transform transformation courante */
  void GMVoutput( ostream &fileOut, const Transform &transform ) const;
  
  /** @brief Nombre de points pour post-processing avec Paraview */
  int numberOfPoints_PARAVIEW() const;
  
  /** @brief Nombre de polygones elementaires pour post-processing avec 
  Paraview */
  int numberOfCells_PARAVIEW() const;

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
	
  /** @brief Ecrit le convexe au format STL pour lien avec openFoam 
  @param f Flux de sortie 
  @param transform transformation courante */
  void write_convex_STL(ostream &f,const Transform &transform) const;		
  
  /** @brief Ecrit le convexe pour post-processing avec Paraview 
  @param connectivity connectivite des polygones Paraview
  @param offsets decalage dans la connectivite
  @param cellstype type de polygones Paraview
  @param firstpoint_globalnumber numero global du 1er point
  @param last_offset dernier offset utilise pour le convexe precedent */
  void write_polygonsStr_PARAVIEW( list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const; 
	
  /** @brief Nombre de points par quart de cercle pour post-processing 
  @param nbpts nombre de points par quart de cercle */
  static void SetvisuNodeNbPerQar( int nbpts ) { visuNodeNbPerQar = nbpts; };
  
  /** @brief La forme est elle une sphere ? */
  bool isSphere() const {return true;};  
  
  /** @brief Renvoie un vecteur orientation du convex 
  @param transform transformation courante */
  Vecteur vecteurOrientation( Transform const* transform ) const;     
  //@}
  

protected:
  /**@name Parameters */
  //@{
  Scalar radius; /**< le rayon de la sphere */
  int Shrinking;  /**< the shrinking choice */
  Scalar initial_radius; /**< the current radius */
  static int visuNodeNbPerQar; /**< nb de points par quart de cercle pour 
  	post-processing Paraview */
  //@}


  /**@name Methods */
  //@{
  /** @brief Rayon circonscrit de la sphere
  @return Le rayon */
  Scalar BuildRayonRef() const;
  //@} 
  
private:
  /**@name Methods */
  //@{
  /** @brief Ecriture d'une facette triangulaire au format STL à partir 
  des 3 points qui la constituent et du centre de gravite de la sphere
  @param f Flux de sortie 
  @param GC centre de gravite
  @param pp1 point 1
  @param pp2 point 2  
  @param pp3 point 3 */
  void write_STLfacet_sphere( ostream &f, Point const& GC,
  	Point const& pp1,
  	Point const& pp2,	
  	Point const& pp3 ) const;
  //@}     
};

#endif
