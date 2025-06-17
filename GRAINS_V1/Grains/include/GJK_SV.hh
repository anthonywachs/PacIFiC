#ifndef _GJK_SV_HH_
#define _GJK_SV_HH_


#include "Transform.hh"
#include "Convex.hh"


// =============================================================================
//                          / ____|    | | |/ /                               //
//    ___  _ __   ___ _ __ | |  __     | | ' /                                //
//   / _ \| '_ \ / _ \ '_ \| | |_ |_   | |  <                                 //
//  | (_) | |_) |  __/ | | | |__| | |__| | . \                                //
//   \___/| .__/ \___|_| |_|\_____|\____/|_|\_\                               //
//        | |                                                                 //
//        |_|                                                                 //
//                                                                            //
// Copyright 2022 Mattia Montanari, University of Oxford                      //
//                                                                            //
// This program is free software: you can redistribute it and/or modify under //
// the terms of the GNU General Public License as published by Free Software  //
// Foundation, either version 3 of the License. You should have received copy //
// of the GNU General Public License along with this program. If not, visit   //
//                                                                            //
//     https://www.gnu.org/licenses/                                          //
//                                                                            //
// This program is distributed in the hope that it will be useful, WITHOUT    //
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY         //
// FOR A PARTICULAR PURPOSE. See GNU General Public License for details.      //
// =============================================================================
/** @brief The header for the signed volume GJK distance query algorithm.

    The signed volume GJK distance query algorithm using the Johnson subalgorithm 
    with the backup procedure. It supports both single and double floating point
    operations.  see the following for more info on the algorithm:
    https://www.mattiamontanari.com/opengjk/

    @author Mattia Montanari - Jan 2023 - Construction
    @author A.Yazdani - 2024 - Modifying for Grains, adding witness points, and
    implementing device code along with templating */
// =============================================================================
/** @name GJK_SV : External methods */
//@{
/** @brief Returns the minimal distance between 2 convex shapes and a point per
convex shape that represents the tips of the minimal distance segment
@param a convex shape A
@param b convex shape B
@param a2w geometric tramsformation describing convex A in the world reference
frame
@param b2w geometric tramsformation describing convex B in the world reference
frame
@param pa point representing one tip of the minimal distance segment on A
@param pb point representing the other tip of the minimal distance segment on B
@param nbIter number of iterations of GJK for convergence */
double closest_points_GJK_SV( Convex const& a, 
                              Convex const& b, 
                              Transform const& a2w,
                              Transform const& b2w, 
                              Point3& pa,
                              Point3& pb,
                              int& nbIter );

/** @brief Returns the minimal distance between 2 convex shapes and a point per
convex shape that represents the tips of the minimal distance segment
@param a convex shape A
@param b convex shape B
@param a2w geometric tramsformation describing convex A in the world reference
frame
@param b2w geometric tramsformation describing convex B in the world reference
frame
@param v initial search direction
@param pa point representing one tip of the minimal distance segment on A
@param pb point representing the other tip of the minimal distance segment on B
@param nbIter number of iterations of GJK for convergence */
double closest_points_GJK_SV( Convex const& a, 
                              Convex const& b, 
                              Transform const& a2w,
                              Transform const& b2w, 
                              Vector3& v,
                              Point3& pa,
                              Point3& pb,
                              int& nbIter );
//@}


#endif
