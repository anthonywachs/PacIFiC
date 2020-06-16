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
#ifndef _CONVEX_H_
#define _CONVEX_H_

#include "Shape.H"

#include "BBox.H"
#include "Transform.H"

#include "Point.H"
using namespace solid;

#include <vector>
#include <list>
using namespace std;

#include "ReaderXML.hh"
#include "WriterXML.hh"


class Convex;
ostream &operator << (ostream &fileOut, const Convex &convex);
istream &operator >> (istream &fileIn, Convex &convex);


enum ConvexType {
  SPHERE,
  DISQUE2D,
  POLYHEDRON,
  POLYGON,
  BOX,
  CONE,
  CYLINDER,
  ANTICYLINDER,
  PLAN,
  POINT,
  SEGMENT,
  SIMPLEX
};


/** @brief Gestion des objets de forme convexe.

    @author 
    Gino van den Bergen - Eindhoven University of Technology - Creation.
    @author 
    D.PETIT - Institut Francais du Petrole - 2000 - Modification. 
    D. RAKOTONIRINA - IFP Energies Nouvelles - 2014 - Modification. */
// ============================================================================
class Convex : public Shape
{
  // ~~~~~~~~~~ Original SOLID Software ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
  /**@name Constructors */
  //@{
  /** @brief Destructeur */
  virtual ~Convex();
  //@}


  /**@name Methods Original SOLID Software */
  //@{
  /** @brief Type d'objet */
  ShapeType getType() const { return CONVEX; } 
  //@}
  

  /**@name Methods Virtual Original SOLID Software */
  //@{
  /** @brief Determine la BBox d'une forme convexe suivant sa position 
  et son orientation dans l'espace
  @return La BBox englobant le convexe
  @param t La position et l'orientation dans l'espace */
  virtual BBox bbox( const Transform &t ) const;

  /** @brief Determine le point du support convexe de la forme se
  situant dans une direction et un sens donne
  @return le point sur l'enveloppe convexe obtenu
  @param v la direction et le sens */
  virtual Point support( const Vecteur& v ) const = 0;
  //@}


  // ~~~~~~~~~~ New GRAINS Software ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
  /**@name Methods */
  //@{
  /** @brief Determine le rayon circonscrit du convexe
  avec application de la matrice de transformation sur celui-ci.
  IMPORTANT: En fait, l'utilisation de la transformation n'est necessaire que 
  si on deforme la particule (allongement, compression ou grossissement) par 
  rapport a la forme de base definie dans le fichier de mise en donnee. 
  Ex: on definit un cube puis on l'etire pour en faire un parallelepipede
  rectangle ou on le grossit pour en faire un cube plus gros.
  En pratique, cette option n'est QUASIMENT JAMAIS utilisee !!
  Pour le cas general, la tranformation est une (ROTATION o TRANSLATION), soit
  la matrice de tranformation associee est de type ROTATION !! 
  et le rayon circonscrit est INDEPENDANT d'une transformation type 
  (ROTATION o TRANSLATION).
  @return Le rayon
  @param t position La transformation a appliquer. */
  Scalar BuildRayon( const Transform &t ) const;
  //@}
  /** @ get shrinking choice.
  @param shrink. */
  int getShrinkingMode(); 

   /** @ set the shrinking radius.
  @param CurrentRadius. */
  void set_shrinking_radius(Scalar CurrentRadius=0); 
 
  /**@name Methods Virtual */
  //@{
  /** @brief Convex Type
  @return The convex type */
  virtual ConvexType getConvexType() const = 0;

  /** @brief Determine l'inertie d'un convexe 
  @return Faux si la fonction est appelee au moins deux fois 
  dans le cas du polyedre.
  @param inertie L'inertie du convexe
  @param inertie_1 L'inertie inverse du convexe */
  virtual bool BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const = 0;

  /** @brief Clone un convexe 
  @return Une copie du convexe */
  virtual Convex* clone() const = 0;

  /** @brief Points decrivant l'enveloppe du Convex.
  @return Un vecteur de points decrivant l'enveloppe du convexe. */
  virtual vector<Point> getEnveloppe() const;

  /** @brief Description des faces
  @return La description des faces via les indices des points sommets */
  virtual vector<vector<int> > const* getFaces() const;
  
  /** @brief Renvoie le nombre de sommets ou un code équivalent */
  virtual int getNbCorners() const;  

  /** @brief Volume du convexe 
  @return Le volume */
  virtual Scalar getVolume() const = 0;
  
  /** @brief Ecriture de l'objet dans sa configuration courante pour
  post-processing avec GMV
  @param fileOut Flux de sortie. 
  @param transform transformation courante */
  virtual void GMVoutput( ostream &fileOut, const Transform &transform ) const; 
	
  /** @brief Nombre de points pour post-processing avec Paraview */
  virtual int numberOfPoints_PARAVIEW() const;
  
  /** @brief Nombre de polygones elementaires pour post-processing avec 
  Paraview */
  virtual int numberOfCells_PARAVIEW() const;

  /** @brief Ecrit les points du convexe pour post-processing avec Paraview  
  @param f Flux de sortie 
  @param transform transformation courante 
  @param translation translation du centre de gravite */
  virtual void write_polygonsPts_PARAVIEW( ostream &f, 
  	const Transform &transform, 
  	Vecteur const* translation = NULL ) const;
	
  /** @brief Ecrit le convexe au format STL pour lien avec openFoam 
  @param f Flux de sortie 
  @param transform transformation courante */
  virtual void write_convex_STL( ostream &f, const Transform &transform )
  	const;	

  /** @brief Renvoie les points du convexe pour post-processing avec Paraview  
  @param transform transformation courante 
  @param translation translation du centre de gravite */
  virtual list<Point> get_polygonsPts_PARAVIEW( const Transform &transform,
  	Vecteur const* translation = NULL )
  	const;
  
  /** @brief Ecrit le convexe pour post-processing avec Paraview 
  @param connectivity connectivite des polygones Paraview
  @param offsets decalage dans la connectivite
  @param cellstype type de polygones Paraview
  @param firstpoint_globalnumber numero global du 1er point
  @param last_offset dernier offset utilise pour le convexe precedent */
  virtual void write_polygonsStr_PARAVIEW( list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;

  /** @brief Renvoie le point d'intersection avec la coque du polyèdre d'un
  segment défini par un point interne et un point externe 
  @param PtIN point interne
  @param PtOut point externe */	
  virtual Point intersectionToShell( const Point &PtIN, const Point &PtOut ) 
  	const {return Point();};     	   

  /** @brief La forme est elle une sphere ? */
  virtual bool isSphere() const { return false; };
  
  /** @brief La forme est elle un cylindre ? */
  virtual bool isCylindre() const { return false; };

  /** @brief Determine si il y a intersection entre deux convexes 
  @return Vrai si les deux convexes sont en contact
  @param a Un convexe A
  @param b Un convexe B
  @param a2w Le repere (position + base) du convexe A dans le repere absolu
  @param b2w Le repere du convexe B dans le repere absolu
  @param v vecteur arbitraire (permet l'initialisation de GJK) */
  virtual bool isIn(const Convex& a, const Convex& b, const Transform& a2w,
	      const Transform& b2w, Vecteur& v); 

  /** @brief Determine si il y a intersection entre deux convexes 
  @return Vrai si les deux convexes sont en contact
  @param a un convexe A
  @param b un convexe B
  @param b2a le repere (position + base) du convexe B 
      dans le repere du convexe a
  @param v vecteur arbitraire (permet l'initialisation de GJK) */
  virtual bool isIn(const Convex& a, const Convex& b, 
	      const Transform& b2a, Vecteur& v);


  /** @brief Renvoie un vecteur orientation du convex 
  @param transform transformation courante */
  virtual Vecteur vecteurOrientation( Transform const* transform ) const;
  //@}


  /** @name Methods Friend */
  //@{
  /** @brief Operateur d'ecriture.
  @return Flux recepteur
  @param fileOut Flux recepteur
  @param convex Objet courant. */
  friend ostream &operator << ( ostream &fileOut, const Convex &convex );
  
  /** @brief Operateur de lecture.
  @return Flux emetteur
  @param fileIn Flux emetteur
  @param convex Objet courant */
  friend istream &operator >> ( istream &fileIn, Convex &convex );
  //@}


protected:
  /**@name Contructors */
  //@{
  /** @brief Constructeur par defaut interdit. */
  Convex();
  //@}


  /**@name Methods virtual */
  //@{
  /** @brief Determine le rayon circonscrit du convexe de reference, 
  sans application de la matrice de transformation.
  @return Le rayon */
  virtual Scalar BuildRayonRef() const = 0;
  
  /**@name Methods virtual */
  //@{
  /** @get the shrinking choice form the forme.
  @return shrink */
  virtual int getShrinkingChoice()const=0;  

  /**@name Methods virtual */
  //@{
  /** @set the Shrinking radius in the active particle's class.
   */
  virtual void setShrinkingRadius(Scalar CurrentRadius=0)=0;
  
  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  virtual void printClass( ostream &fileOut ) const = 0;

  /** @brief Lecture de l'objet sur le flux d'entree.
  @param fileIn Flux d'entree. */
  virtual void readClass( istream &fileIn ) = 0;  
  //@}
};


// ~~~~~~~~~~ Original SOLID Software ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**@name Convex : Methods externals Original SOLID Software */
//@{
/** @brief Determine si il y a intersection entre deux convexes 
@return Vrai si les deux convexes sont en contact
@param a Un convexe A
@param b Un convexe B
@param a2w Le repere (position + base) du convexe A dans le repere absolu
@param b2w Le repere du convexe B dans le repere absolu
@param v vecteur arbitraire (permet l'initialisation de GJK) */
bool intersect(const Convex& a, const Convex& b,
	       const Transform& a2w, const Transform& b2w,
	       Vecteur& v);

/** @brief Determine si il y a intersection entre deux convexes 
@return Vrai si les deux convexes sont en contact
@param a un convexe A
@param b un convexe B
@param b2a le repere (position + base) du convexe B 
    dans le repere du convexe a
@param v vecteur arbitraire (permet l'initialisation de GJK) */
bool intersect(const Convex& a, const Convex& b, 
	       const Transform& b2a, 
	       Vecteur& v);

/** @brief Determine si il y a intersection entre deux convexes 
et renvoi un point d'intersection dans le repere de chaque convexe 
(si il y a intersection)
@return Vrai si les deux convexes sont en contact
@param a un convexe A
@param b un convexe B
@param a2w le repere (position + base) du convexe A dans le repere absolu
@param b2w le repere du convexe B dans le repere absolu
@param v  vecteur arbitraire (permet l'initialisation de GJK)
@param pa le point d'intersection dans le repere du convexe A
@param pb le point d'intersection dans le repere du convexe B */
bool common_point(const Convex& a, const Convex& b, 
		  const Transform& a2w, const Transform& b2w,
		  Vecteur& v, 
		  Point& pa, Point& pb);

/** @brief Determine si il y a intersection entre deux convexes 
et renvoi un point d'intersection dans le repere de chaque convexe 
(si il y a intersection)
@return Vrai si les convexes sont en contact
@param a un convexe A
@param b un convexe B
@param b2a le repere (position + base) du convexe B
@param v  vecteur arbitraire (permet l'initialisation de GJK)
@param pa le point d'intersection dans le repere du convexe A
@param pb le point d'intersection dans le repere du convexe B */
bool common_point(const Convex& a, const Convex& b, 
		  const Transform& b2a,
		  Vecteur& v, 
		  Point& pa, Point& pb);

/** @brief Determine la distance entre deux convexes et renvoi un point 
de chaque convexe donnant le distance minimum
@return La distance entre les deux convexes
@param a Un convexe A
@param b Un convexe B
@param a2w Le repere (position + base) du convexe A dans le repere absolu
@param b2w Le repere du convexe B dans le repere absolu
@param pa Un point de a formant la plus courte distance 
	dans le repere du convexe A
@param pb Un point de b formant la plus courte distance 
dans le repere du convexe B 
@param nbIter nombre d'iterations de l'algorithme GJK */
Scalar closest_points(const Convex& a, const Convex& b, 
		      const Transform& a2w, const Transform& b2w,
		      Point& pa, Point& pb,int& nbIter);
//@}

#endif
