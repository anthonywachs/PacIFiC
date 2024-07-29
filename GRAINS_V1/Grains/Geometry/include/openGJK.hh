//                           _____      _ _  __                                   //
//                          / ____|    | | |/ /                                   //
//    ___  _ __   ___ _ __ | |  __     | | ' /                                    //
//   / _ \| '_ \ / _ \ '_ \| | |_ |_   | |  <                                     //
//  | (_) | |_) |  __/ | | | |__| | |__| | . \                                    //
//   \___/| .__/ \___|_| |_|\_____|\____/|_|\_\                                   //
//        | |                                                                     //
//        |_|                                                                     //
//                                                                                //
// Copyright 2022 Mattia Montanari, University of Oxford                          //
//                                                                               //
// This program is free software: you can redistribute it and/or modify it under  //
// the terms of the GNU General Public License as published by the Free Software  //
// Foundation, either version 3 of the License. You should have received a copy   //
// of the GNU General Public License along with this program. If not, visit       //
//                                                                                //
//     https://www.gnu.org/licenses/                                              //
//                                                                                //
// This program is distributed in the hope that it will be useful, but WITHOUT    //
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  //
// FOR A PARTICULAR PURPOSE. See GNU General Public License for details.          //

/**
 * @file openGJK.h
 * @author Mattia Montanari
 * @date 1 Jan 2023
 * @brief Main interface of OpenGJK containing quick reference and API documentation.
 *
 * @see https://www.mattiamontanari.com/opengjk/
 */

#ifndef OPENGJK_H__
#define OPENGJK_H__

#include <float.h>
#include "Transform.hh"
#include "Convex.hh"

/*! @brief Precision of floating-point numbers.
 *
 * Default is set to 64-bit (Double). Change this to quickly play around with 16- and 32-bit. */
// #ifdef USE_32BITS
#define gkFloat   double
#define gkEpsilon DBL_EPSILON
// #else
// #define gkFloat   double
// #define gkEpsilon DBL_EPSILON
// #endif

/*! @brief Data structure for simplex.
   *
   * The simplex is updated at each GJK-iteration. For the first iteration this value is a guess - and this guess not irrelevant. */
typedef struct gkSimplex_ {
  int nvrtx;          /*!< Number of points defining the simplex. */
  gkFloat vrtx[4][3]; /*!< Coordinates of the points of the simplex. */
} gkSimplex;

/*! @brief Invoke the GJK algorithm to compute the minimum distance between two polytopes.
   *
   * The simplex has to be initialised prior the call to this function. */
double closest_points_GJK_SV2( Convex const& a, 
                       Convex const& b, 
                       Transform const& a2w,
                       Transform const& b2w, 
                       Point3& pa,
                       Point3& pb,
                       int& nbIter );

#endif // OPENGJK_H__
