/*
  3D - C++ Class Library for 3D Transformations
  Copyright (C) 1996-1998  Gino van den Bergen

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
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <assert.h>

#include "Vecteur.H"
#include "Quaternion.H"

//@{
/// Matrice reels 3x3
typedef Scalar Mat3[3][3];
//@}


class Matrix;
ostream &operator << ( ostream &fileOut, const  Matrix &m );
istream &operator >> ( istream &fileIn, Matrix &m );


/** Classe de gestion des matrices 3x3.

    @author 
    Gino van den Bergen -  Eindhoven University of Technology - Creation */
// ============================================================================
class Matrix 
{
public:
  /**@name Constructeurs */
  //@{
  /** @brief Constructeurs par defaut.
  La matrice est initialisee a la matrice identite. */
  Matrix();
  
  /** @brief Constructeur avec initialisation par un tableau de valeurs
  @param m le tableau contenant les coefficients de la matrice */
  Matrix( const Scalar *m );
  
  /** @brief Constructeur avec initialisation par un quaternion (norme)
  @param q le quaternion */
  Matrix( const Quaternion& q );
  
  /** @brief Constructeur d'une matrice diagonale
  @param x le coefficient a(1,1) de la matrice a
  @param y le coefficient a(2,2) de la matrice a
  @param z le coefficient a(3,3) de la matrice a */
  Matrix( Scalar x, Scalar y, Scalar z );
  
  /** @brief Constructeur d'une matrice avec tous les coefficients en parametres
  @param xx le coefficient a(1,1) de la matrice a
  @param xy le coefficient a(1,2) de la matrice a
  @param xz le coefficient a(1,3) de la matrice a
  @param yx le coefficient a(2,1) de la matrice a
  @param yy le coefficient a(2,2) de la matrice a
  @param yz le coefficient a(2,3) de la matrice a
  @param zx le coefficient a(3,1) de la matrice a
  @param zy le coefficient a(3,2) de la matrice a
  @param zz le coefficient a(3,3) de la matrice a */
  Matrix( Scalar xx, Scalar xy, Scalar xz,
	 Scalar yx, Scalar yy, Scalar yz,
	 Scalar zx, Scalar zy, Scalar zz ); 

  /** @brief Constructeur par copie
  @param other la matrice rhs */
  Matrix( const Matrix &other );

  /** @brief Destructeur */  
  ~Matrix();
  //@}

  
  /** @name Methods Get */
  //@{
  /** @brief Renvoi la matrice courante.
  @return Alias de la matrice */
  Mat3& getValue();

  /** @brief Renvoi la matrice courante.
  @return Alias de la matrice */
  const Mat3& getValue() const;  
  //@}


  /** @name Methods Set */
  //@{
  /** @brief Modifie la matrice courante par la matrice identite. */
  void setIdentity();  
  
  /** @brief Modifie la matrice courante par une matrice rotation
  a partir d'un quaternion.
  @param q Le quaternion de rotation */
  void setRotation( const Quaternion& q );

  /** @brief Modifie la matrice courante par une matrice diagonale.
  @param x
  @param y
  @param z */
  void setScaling( Scalar x, Scalar y, Scalar z ); 

  /** @brief Modifie la matrice en courante.
  !!! IMPORTANT !!! cette methode impose de passer la matrice comme un tableau
  de 12 valeurs telle une suite de 3 *(1 vecteur colonne + 1 double) 
  (le double est inutile, 0 par defaut); ainsi l'ordre des
  elements dans le tableau est le suivant: 0=(1,1), 1=(2,1), 2=(3,1), 3=0, 
  4=(1,2), 5=(2,2), 6=(3,2), 7=0, 8=(1,3), 9=(2,3), 10=(3,3), 11=0.
  @param m le tableau de valeur */
  void setValue( const Scalar *m );
  
  /** @brief Modifie la matrice en courante
  @param xx le coefficient a(1,1) de la matrice
  @param xy le coefficient a(1,2) de la matrice
  @param xz le coefficient a(1,3) de la matrice
  @param yx le coefficient a(2,1) de la matrice
  @param yy le coefficient a(2,2) de la matrice
  @param yz le coefficient a(2,3) de la matrice
  @param zx le coefficient a(3,1) de la matrice
  @param zy le coefficient a(3,2) de la matrice
  @param zz le coefficient a(3,3) de la matrice */
  void setValue( Scalar xx, Scalar xy, Scalar xz, 
		Scalar yx, Scalar yy, Scalar yz, 
		Scalar zx, Scalar zy, Scalar zz );
  //@}


  /** @name Methods */
  //@{
  /** @brief Renvoie la matrice construite a partir de la matrice courante et 
  ayant tous ses coefficients positifs.
  @return La matrice obtenue */
  Matrix absolute() const;

  /** @brief Renvoi la matrice adjointe de la matrice courante 
  @return La matrice adjointe */
  Matrix adjoint() const;

  /** @brief Determinant de la matrice courante 
  return le determinant de la matrice*/    
  Scalar determinant() const;

  /** @brief Renvoi la matrice inverse de la matrice courante.
  @return la matrice inverse */
  Matrix inverse() const; 

  /** @brief Produit scalaire entre une colone de la matrice et un vecteur
  @return le produit scalaire entre la ligne et le vecteur 
  @param i
  @param v Le vecteur*/
  Scalar tdot( int i, const Vecteur& v ) const ;

  /** @brief Renvoi la matrice transpose de la matrice courante 
  @return la matrice transpose */
  Matrix transpose() const;  
  
  /** @brief Copie de la matrice de transformation 
  !!! IMPORTANT !!! cette methode impose de copier la matrice dans le tableau
  de copie comme un tableau
  de 12 valeurs telle une suite de 3 *(1 vecteur colonne + 1 double) 
  (le double est inutile, 0 par defaut); ainsi l'ordre des
  elements dans le tableau est le suivant: 0=(1,1), 1=(2,1), 2=(3,1), 3=0, 
  4=(1,2), 5=(2,2), 6=(3,2), 7=0, 8=(1,3), 9=(2,3), 10=(3,3), 11=0.
  @param vit tableau de copie 
  @param i position de depart dans le tableau de copie */
  void copyMatrix( double *vit, int i ) const;
  
  /** @brief Ecriture de l'objet sur le flux de sortie avec un format de
  precision elevee
  @param fileOut Flux de sortie. */
  void writeMatrix( ostream &fileOut ) const;
  
  /** @brief Ecriture de l'objet sur le flux de sortie avec un format de
  precision elevee avec le format de reload 2014
  @param fileOut Flux de sortie. */
  void writeMatrix2014( ostream &fileOut ) const; 
  
  /** @brief Ecriture en binaire de l'objet avec le format de reload 2014
  @param fileOut Flux de sortie. */
  void writeMatrix2014_binary( ostream &fileOut ); 
  
  /** @brief Lecture en binaire de l'objet avec le format de reload 2014 en
  binaire
  @param StreamIN Flux d'entree. */
  void readMatrix2014_binary( istream &StreamIN );        
  
  /** @brief Multiplie à droite par une matrice de scaling
  @param x scaling en x
  @param y scaling en y
  @param z scaling en z */
  void multiplyByScalingMatrix( Scalar x, Scalar y, Scalar z );    
  //@}


  /**@name Operators */
  //@{
  /** @brief Operateur d'indexation (pour affecter). 
  @return La premiere ligne de la matrice.
  @param i Indice de la ligne */
  Vecteur& operator[]( int i );
  
  /** @brief Operateur d'indexation (pour lire). 
  @return La premiere ligne de la matrice.
  @param i Indice de la ligne */
  const Vecteur& operator[]( int i ) const;

  /** @brief Addition matricielle
  @return la matrice somme
  @param m une matrice */
  Matrix& operator+=( const Matrix& m );

  /** @brief Multiplication matricielle
  @return la matrice produit
  @param m une matrice */
  Matrix& operator*=( const Matrix& m );
  
  /** @brief Operation =
  @return la matrice affectee
  @param other la matrice rhs */
  Matrix& operator=( const Matrix& other );  
  //@}
  

  /**@name Methods friend */
  //@{
  /** @brief Operateur d'ecriture.
  @return Flux recepteur
  @param fileOut Flux recepteur
  @param m Objet courant. */
  friend ostream &operator << ( ostream &fileOut, const  Matrix &m );
			       
  /** @brief Operateur de lecture.
  @return Flux emetteur
  @param fileIn Flux emetteur
  @param m Objet courant */
  friend istream &operator >> ( istream &fileIn, Matrix &m );
  //@}

protected:
  /**@name Parameters */
  //@{
  Mat3 elem; /**< Element de type matrice 3x3 */
  //@}
};


/** @name Matrix : Methods external*/
//@{
/** @brief Multiplication entre une matrice et un vecteur
@return le vecteur obtenu par le produit
@param m la matrice
@param v le vecteur */
Vecteur operator*( const Matrix& m, const Vecteur& v );

/** @brief Multiplication entre un vecteur et une matrice
@return le vecteur obtenu par le produit
@param v le vecteur
@param m la matrice */
Vecteur operator*( const Vecteur& v, const Matrix& m );

/** @brief Multiplication entre deux matrice
@return la matrice obtenu par le produit
@param m1 une matrice
@param m2 une matrice */
Matrix operator*( const Matrix& m1, const Matrix& m2 );

/** @brief Multiplication a gauche par une matrice transpose
@return la matrice obtenu par le produit
@param m1 la matrice a transpose dans le produit
@param m2 une matrice */
Matrix multTransposeLeft( const Matrix& m1, const Matrix& m2 );

/** @brief Transposition de matrice
@return la matrice transpose
@param m la matrice a transposer */
Matrix transpose( const Matrix& m );

/** @brief Obtention d'un matrice adjointe
@return la matrice adjointe
@param m une matrice */
Matrix adjoint( const Matrix& m );

/** @brief Inversion de matrice
@return la matrice inverse
@param m une matrice */
Matrix inverse( const Matrix& m );

/** @brief Obtention d'une matrice dont tous les coefficients sont positifs
@return la matrice a coefficients positifs
@param m une matrice */
Matrix absolute( const Matrix& m );
//@}
 
#endif
