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
#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include "Point.H"
#include "Matrix.H"
#include "Vecteur.H"
#include "Quaternion.H"
#include "Basic.H"

#include "ReaderXML.hh"
#include "WriterXML.hh"


class Transform;
ostream &operator << ( ostream &fileOut, const Transform &t );
istream &operator >> ( istream &fileIn, Transform &t );


/** @brief Gestion de la position et de l'orientation d'une forme : 
    modelise par une transformation affine.
    
    La position de la forme est donnee par un point et l'orientation par une 
    matrice ou les colonnes representent les vecteurs de base de la forme 
    dans le repere absolu. 

    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class Transform 
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut.
  le type est IDENTITY, l'origine est positionnee a {0.,0.,0.}, 
  la matrice de base est initialisee a la matrice identite. */
  Transform();

  /** @brief Constructeur avec initialisation de l'origine
  le type est IDENTITY et
  la matrice de base est initialisee a la matrice identite. 
  @param gx coordonnée x de l'origine 
  @param gy coordonnée y de l'origine 
  @param gz coordonnée z de l'origine */
  Transform( const double &gx, const double &gy, const double &gz );

  /** @brief Constructeur avec initialisation 
  @param m Le repere et la position de la forme */
  Transform( const Scalar m[16] );

  /** @brief Constructeur par copie
  @param other la transformation rhs */
  Transform( const Transform &other );

  /** @brief Destructeur */
  ~Transform() {}
  //Removed ~Transform() {} and placed the above ~Transform(); instead
  //@}
  

  /**@name Methods */
  //@{ 
  /** @brief Inversion de la transformation
      @param t La transformation a inverser */ 
  void invert( const Transform& t );

  /** @brief La matrice est-elle matrice identite.
  @return Vrai si identite. */
  bool isIdentity() const;

  /** @brief Composition de transformation affine ( this = t1 o t2 )
  @param t1 une transformation affine 
  @param t2 une transformation affine */
  void mult( const Transform& t1, const Transform& t2 );
  
  /** @brief Composition avec une transformation de scaling 
  this = this o scaling
  @param x scaling en x
  @param y scaling en y
  @param z scaling en z */
  void composeWithScaling( Scalar x, Scalar y, Scalar z );  

  /** @brief Composition de transformation affine ( this = ( t1 )-1 o t2 )
  @param t1 une transformation affine ( celle qui sera inversee )
  @param t2 une transformation affine */
  void multInverseLeft( const Transform& t1, const Transform& t2 ); 
  
  /** @brief Scaling de l'objet ( homothetie dans les trois directions )
  @param x homothetie suivant l'axe x dans le repere locale
  @param y homothetie suivant l'axe y dans le repere locale
  @param z homothetie suivant l'axe z dans le repere locale*/
  void scale( Scalar x, Scalar y, Scalar z );

  /** @brief Ecriture de l'objet sur le flux de sortie.
  @param fileOut Flux de sortie. */
  void printClass( ostream &fileOut ) const;
  
  /** @brief Ecriture de l'objet sur le flux de sortie avec un format de
  precision elevee
  @param fileOut Flux de sortie. */
  void writeTransform( ostream &fileOut ) const;
  
  /** @brief Ecriture de l'objet sur le flux de sortie avec un format de
  precision elevee avec le format de reload 2014
  @param fileOut Flux de sortie. */
  void writeTransform2014( ostream &fileOut ) const;
  
  /** @brief Ecriture en binaire de l'objet sur le flux de sortie 
  avec le format de reload 2014
  @param fileOut Flux de sortie. */
  void writeTransform2014_binary( ostream &fileOut );  
  
  /** @brief Lecture de l'objet sur le flux d'entre avec le format 
  de reload 2014
  @param StreamIN Flux d'entree. */
  void readTransform2014( istream &StreamIN ); 
  
  /** @brief Lecture de l'objet sur le flux d'entree avec le format 
  de reload 2014 en binaire
  @param StreamIN Flux d'entree. */
  void readTransform2014_binary( istream &StreamIN );        

  /** @brief Lecture de la transformation sur le flux d'entree.
  @param fileIn Flux d'entree. */
  void readClass( istream &fileIn );

  /** @brief Composition par une rotation 
  @param q le quaternion representant la rotation */
  void rotate( const Quaternion& q );

  /** @brief Rotation de l'objet
  @param q le quaternion representant la rotation */
  void rotateOrientation( const Quaternion& q );

  /** @brief Translation de la position de l'objet 
  @param v le vecteur de translation */
  void translate( const Vecteur& v );
  
  /** @brief Copie la Transformation du Convex decrivant la Forme
  !!! IMPORTANT !!! cette methode impose de copier la transformation
  dans le tableau de copie comme un tableau de 16 valeurs telle une suite de 
  3 *(1 vecteur colonne + 1 double) (le double est inutile, 0 par defaut) pour
  la matrice de transformation, puis
  le centre de gravite + 1 double (le double est inutile, 0 par defaut);
  ainsi l'ordre des elements dans le tableau est le suivant: 
  0=Mat(1,1), 1=Mat(2,1), 2=Mat(3,1), 3=0, 4=Mat(1,2), 
  5=Mat(2,2), 6=Mat(3,2), 7=0, 8=Mat(1,3), 9=Mat(2,3), 10=Mat(3,3), 11=0, 
  12=GC(X), 13=GC(Y), 14=GC(Z), 15=0.
  @param vit tableau de copie 
  @param i position de depart dans le tableau de copie */
  void copyTransform( double *vit, int i ) const; 
  
  /** @brief Copie la Transformation du Convex decrivant la Forme translatée
  d'un vecteur (utile pour les particules periodiques en parallele)
  @param vit tableau de copie 
  @param i position de depart dans le tableau de copie 
  @param vec vecteur de translation */
  void copyTransform( double *vit, int i, Vecteur const& vec ) const; 
  
  /** @brief Composition a gauche par une autre transformation affine
  !!! Resultat: d'abord t puis this !!! 
  @param t une transformation affine */
  void composeTransformLeft( const Transform& t ); 
  
  /** @brief Composition a droite par une autre transformation affine
  !!! Resultat: d'abord this puis t !!! 
  @param t une transformation affine */
  void composeTransformRight( const Transform& t ); 
  
  /** @brief Composition a droite par une rotation rot par rapport a l'origine 
  de this
  !!! Resultat: d'abord this puis r !!! 
  @param rot une transformation affine */
  void composeRotationRight( const Transform& rot );          
  //@}


  /**@name Methods Get */
  //@{
  /** @brief Renvoi les vecteurs de base sous forme de matrice
  @return Les vecteurs de base en colonne */
  const Matrix& getBasis() const { return basis; }


  /** @brief Renvoi le type de l'objet */
  unsigned int getType() const { return type; }  

  /** @brief Renvoi la position de l'objet
  @return la position */
  const Point* getOrigin() const { return &origin; }  
  //@}  


  /**@name Methods Set */
  //@{
  /** @brief Initialise la position a l'origine du repere absolue.
  et modifie la base en base canonique */
  void setIdentity();

  /** @brief Modifie la position de l'objet.
  @param pos la nouvelle position de l'objet */
  void setOrigin( const Scalar* pos );
  
  /** @brief Modifie la position de l'objet.
  @param pos Nouvelle position de l'objet */
  void setOrigin( const Point &pos );
  
  /** @brief Modifie la position de l'objet.
  @param gx coordonnée x de la position
  @param gy coordonnée y de la position  
  @param gz coordonnée z de la position */
  void setOrigin( const double &gx, const double &gy, const double &gz );
  
  /** @brief Modifie la matrice de transformation
  @param basis_ Nouvelle position de l'objet */
  void setBasis( const Matrix &basis_ );    

  /** @brief Modifie la position et la base de l'objet.
  !!! IMPORTANT !!! cette methode impose de passer la transformation
  comme un tableau de 16 valeurs telle une suite de 
  3 *(1 vecteur colonne + 1 double) (le double est inutile, 0 par defaut) pour
  la matrice de transformation, puis
  le centre de gravite + 1 double (le double est inutile, 0 par defaut);
  ainsi l'ordre des elements dans le tableau est le suivant: 
  0=Mat(1,1), 1=Mat(2,1), 2=Mat(3,1), 3=0, 4=Mat(1,2), 
  5=Mat(2,2), 6=Mat(3,2), 7=0, 8=Mat(1,3), 9=Mat(2,3), 10=Mat(3,3), 11=0, 
  12=GC(X), 13=GC(Y), 14=GC(Z), 15=0.
  @param m Contient le repere et la position de l'objet */
  void setValue( const Scalar m[16] );
  //@}  
  
  
  /** @name Persistance */
  //@{
  /** @brief Persistance du Transform
  @param root Le noeud racine "Orientation" */
  void load( DOMNode* root );
  //@}


  /**@name Operators */
  //@{
  /** @brief Passage d'un point du repere local au repere global
  @return le point dans le repere global
  @param p un point dans le repere local */
  Point operator() ( const Point& p ) const;

  /** @brief Operateur d'egalite : affecte la transformation donnee a *this 
  @param transform transformation à affecter */
  Transform& operator= ( const Transform& transform );
  //@}
  

  /**@name Methods friend */
  //@{
  /** @brief Operateur d'ecriture.
  @param fileOut Flux recepteur
  @param t Objet courant. 
  @return Flux recepteur */
  friend ostream &operator << ( ostream &fileOut, const Transform &t );

  /** @brief Operateur de lecture.
  @param fileIn Flux emetteur
  @param t Objet courant 
  @return Flux emetteur */
  friend istream &operator >> ( istream &fileIn, Transform &t );
  //@}


private:
  enum { 
    /**@name Enumeration */
    //@{
    /** Identite */
    IDENTITY = 0x00,
    /** Translation */
    TRANSLATION = 0x01,
    /// Rotation
    ROTATION = 0x02,
    /// Scaling (homothetie dans les trois directions)
    SCALING = 0x04,
    /// transformation lineaire
    LINEAR = ROTATION | SCALING,
    /// transformation affine
    AFFINE = TRANSLATION | LINEAR
    //@}
  }; /**< Type de transformation possible */

  /**@name Parameters */
  //@{
  Matrix basis; /**< Matrice representant les vecteurs de base de l'objet */
  Point origin; /**< Position de l'objet */
  unsigned int type; /**< Type de transformation deja effectuee */
  //@}
};

static Transform TransformIdentity; /* transformation a l'identique */

#endif
