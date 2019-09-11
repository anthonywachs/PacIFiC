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
#ifndef _BOX_H_
#define _BOX_H_

#include "Convex.H"
#include "ReaderXML.hh"
#include "Vecteur.H"
#include "Erreur.H"
using namespace solid;


/** @brief Gestion des convexes de la forme d'un pave.
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Box : public Convex 
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur avec initialisation
  @param x Longueur de la Box axe x
  @param y Longueur de la Box axe y
  @param z Longueur de la Box axe z */
  Box(Scalar x = 0., Scalar y = 0., Scalar z = 0.);

  /** @brief sets current radius
  @param ShrinkingRadius */
  void setShrinkingRadius(Scalar CurrentRadius); 

  /** @brief  get the shrinking choice */
  int getShrinkingChoice()const;
  
  /**@brief  Constructeur avec initialisation.
  @param extent_ Le vecteur de description de la Box. */
  Box(const Vecteur extent_);

  /** @brief Constructeur a partir d'un flux entrant
  @param fileIn le flux entrant */
  Box(istream &fileIn);

  /** @brief Constructeur
  @param root Le noeud "Box" */
  Box(DOMNode* root);

  /** @brief Destructeur */
  ~Box();
  //@}


  /** @name Methods */
  //@{
  /** @brief Determine l'inertie d'un pave
  @param inertie L'inertie du pave
  @param inertie_1 L'inertie inverse du pave */
  bool BuildInertie(Scalar *inertie,Scalar *inertie_1) const;

  /** @brief Clone d'un pave
  @return le clone d'un pave */
  Convex* clone() const;

  /** @brief Convex Type
  @return The convex type */
  ConvexType getConvexType() const {return BOX;};

  /** @brief Renvoi les longueurs du pave dans chaque direction 
  a partir son origine
  @return le vecteur representant les longueurs */
  const Vecteur& getExtent() const { return extent; }

  /** @brief Points decrivant l'enveloppe */
  vector<Point> getEnveloppe() const;
  
  /** @brief Renvoi le point[0] de la face i
  @param i numero de face */
  Point getFirstPointFace(int i) const;   

  /** @brief Description des faces
  @return La description des faces via les indices des points sommets */
  vector< vector<int> > const* getFaces() const;

  /** @brief Renvoie le nombre de sommets ou un code équivalent */
  int getNbCorners() const;  

  /** @brief Determine le volume d'un pave. */
  Scalar getVolume() const;

   /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  void printClass(ostream &fileOut) const;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  void readClass(istream &fileIn);

  /** @brief Determine le point du support convexe de la forme se
  situant dans une direction et un sens donne
  Pour le pave : on choisit le sommet le plus 'proche' du vecteur 
  ou de projection sur le vecteur maximal
  @return le point sur l'enveloppe convexe obtenu
  @param v la direction et le sens */
  Point support(const Vecteur& v) const;
  
  /** @brief Ecriture de l'objet dans sa configuration courante pour
  post-processing avec GMV
  @param fileOut flux de sortie
  @param transform transformation courante */
  void GMVoutput(ostream &fileOut,const Transform &transform) const;  

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
	
  /** @brief Ecrit le convexe pour post-processing avec Paraview 
  @param connectivity connectivite des polygones Paraview
  @param offsets decalage dans la connectivite
  @param cellstype type de polygones Paraview
  @param firstpoint_globalnumber numero global du 1er point
  @param last_offset dernier offset utilise pour le convexe precedent */
  void write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const;  
  
  /** @brief Renvoie le point de contact dans l'espace du convexe et la distance
  de recouvrement entre la boite et une sphere. Si le contact n'existe pas,
  renvoie le point origine (0,0,0)
  @param SphereCenter centre de la sphere dans l'espace du convexe Box
  @param SphereRadius rayon de la sphere
  @param overlap recouvrement (negatif si contact) */
  Point IntersectionPointSPHERE(const Point& SphereCenter,
  	const Scalar& SphereRadius,
	double &overlap) const     
    throw(ErreurContact); 
    
  /** @brief Renvoie le point d'intersection avec la coque du polyèdre d'un
  segment défini par un point interne et un point externe 
  @param PtIN point interne
  @param PtOut point externe */
  Point intersectionToShell(const Point &PtIN, const Point &PtOut) const;
  //@}
  
  /** @brief Send the projected sphere center on the obstacle and the gap
  between sphere and obstacle. If it is contact, it sends back origin (0,0,0)
  @param SphereCenter center of sphere in Box space
  @param SphereRadius Radius of the sphere  
  @param gap normal distance between sphere surface and the wall (obstacle face 
  only !) */
  Point ProjectedPointSPHERE( const Point& SphereCenter,
  	const Scalar& SphereRadius,
	double &gap ) const;

  
protected:
  /**@name Methods */
  //@{
  /** @brief Rayon circonscrit de la sphere
  @return Le rayon */
  Scalar BuildRayonRef() const;
  //@}   


  /** @name Parameters*/
  //@{  
  Vecteur extent; /**< Represente la demi-longueur du pave 
  	dans chaque direction */
  vector<Point> corners; /**< les 8 sommets du pave */
  vector<Point>* corners2D_XY; /**< les 4 sommets du pave dans le plan XY */  
  static vector< vector<int> > allFaces; /**< description des faces 
  	formant le pavé par les numéros des sommets */
  //@}
  
  
private:
  /**@name Methods */
  //@{
  /** @brief Construit les sommets du pavé et les faces */
  void setCornersFaces();
  
  /** @brief Renvoie le point de contact dans l'espace du convexe et la distance
  de recouvrement entre un corner de la boite et une sphere. Si le contact 
  n'existe pas, renvoie le point origine (0,0,0)
  @param SphereCenter centre de la sphere dans l'espace du convexe Box
  @param SphereRadius rayon de la sphere
  @param cornerNumber numero du corner
  @param overlap recouvrement (negatif si contact) */
  Point ContactCornerSPHERE(const Point& SphereCenter,
  	const Scalar& SphereRadius,
	const int &cornerNumber,
	double &overlap) const; 
	
  /** @brief Renvoie le point de contact dans l'espace du convexe et la distance
  de recouvrement entre une arete de la boite et une sphere. Si le contact 
  n'existe pas, renvoie le point origine (0,0,0)
  @param SphereCenter centre de la sphere dans l'espace du convexe Box 
  @param SphereRadius rayon de la sphere
  @param cornerNumber numero d'un des 2 corners de l'arete
  @param projectionDirection plan sur lequel on projete le centre de la sphere
  @param overlap recouvrement (negatif si contact) */
  Point ContactEdgeSPHERE(const Point& SphereCenter,
  	const Scalar& SphereRadius,
	const int &cornerNumber,
	const int &projectionDirection,
	double &overlap) const;	       
  //@} 
};

#endif
