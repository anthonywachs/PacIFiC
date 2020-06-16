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
#ifndef _BASIC_H_
#define _BASIC_H_

#include <math.h>
#include <stdlib.h>


//@{
/// type de reel utiliser dans le code
typedef double Scalar;
//@}

/**@name Constants */
//@{
/// Dimension de l'espace : 3
//const int DIM = 3;
/// taille d'une ligne de caracteres : 255
const int long_chaine = 255;
/// constante de scaling (homothetie) suivant l'axe des X
const Scalar X_Scale = 3.0/4.0;
/// constante de scaling (homothetie) suivant l'axe des Y
const Scalar Y_Scale = 3.0/4.0;
/// constante de scaling (homothetie) suivant l'axe des X
const Scalar Z_Scale = 3.0/4.0;
/// constante de scaling (homothetie)
const Scalar SCALE = 3.0/4.0;
/// Constante de gravitation
const Scalar G = 9.81;
/// constante de conversion des degres en radians
const Scalar DEGS_PER_RAD = 57.29577951308232286465;
/// constante de conversion des radians en degres
const Scalar RADS_PER_DEG =  0.01745329251994329547;
/// Pi
const Scalar PI =            3.14159265358979323846;
/// 2.Pi
const Scalar TWO_PI =        6.28318530717958623200;
/// erreur approximation
const Scalar EPS = 1.0e-10;
/// erreur commise (faible precision)
const Scalar EPSILON =       1.0e-10;
/// erreur commise (haute precision)
const Scalar EPSILON2 =      1.0e-15;
/// erreur commise (haute precision)
const Scalar EPSILON3 =      1.0e-20;
/// constante simulant l'infini pour les reels
//const Scalar INFINITY =      1.0e50;
/// nombre de caractere pour l'ecriture du temps
const int TIMEFORMAT = 10;
/// nombre de caractere pour l'ecriture de la position
/// dans les fichiers de reload
const int POSITIONFORMAT = 16;
//@}

/** @name Methods de base */
//@{
/** Valeur absolue avec type Scalar */
//inline Scalar abs(Scalar x) { return x < 0 ? -x : x; }
/** valeur aleatoire reel */
inline Scalar rnd() { return (Scalar(rand()) + 0.5) / (Scalar(RAND_MAX) + 1); }
/** signe d'un reel */ 
inline int    sgn(Scalar x) { return x < 0 ? -1 : x > 0 ? 1 : 0; }
/** test si un reel est approximativement nul */
inline bool   eqz(Scalar x) { return fabs(x) <= EPSILON; }
/** minimum entre deux reels */
inline Scalar min(Scalar x, Scalar y) { return x > y ? y : x; }
/** maximun entre deux reels */
inline Scalar max(Scalar x, Scalar y) { return x < y ? y : x; }
/** la valeur minimun entre deux reels est affectee a chacun de ces reels*/
inline void set_min(Scalar& x, Scalar y) { if (x > y) x = y; }
/** la valeur maximum entre deux reels est affectee a chacun de ces reels*/
inline void set_max(Scalar& x, Scalar y) { if (x < y) x = y; }
/** converti les degres en radians */
inline Scalar rads(Scalar x) { return x * RADS_PER_DEG; }
/** converti les radians en degres */
inline Scalar degs(Scalar x) { return x * DEGS_PER_RAD; }
//@}


/** Direction de l'espace utilise.

    Permet par exemple de remplacer une ecriture v[0] par v[X].
    
    @author G.FERRER - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
enum Direction 
{
  /**@name Enumeration */
  //@{
  /// Direction suivant l'axe des X (longueur)
  X,
  /// Direction suivant l'axe des Y (profondeur)
  Y,
  /// Direction suivant l'axe des Z (hauteur)
  Z,
  /// premiere composante des quaternions
  W
  //@}
};

#endif
