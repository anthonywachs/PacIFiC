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

#include "Cone.H"

// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Constructeur 
Cone::Cone(Scalar r, Scalar h) : 
    bottomRadius(r), 
    quaterHeight(h / 4), 
    sinAngle(r / sqrt(r * r + h * h))  
{
} 


// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Constructeur avec un fichier
Cone::Cone(istream &fileIn)
{
  readClass(fileIn);
}


// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Destructeur
Cone::~Cone() 
{
}


// ----------------------------------------------------------------------
Point Cone::support(const Vecteur& v) const 
{
  Scalar norm = Norm(v);
  if (norm > EPSILON) {
    if (v[Y] > norm * sinAngle) {
      return Point(0, 3.0*quaterHeight, 0);
    } else {
      Scalar s = sqrt(v[X] * v[X] + v[Z] * v[Z]);
      if (s > EPSILON) {
	Scalar d = bottomRadius / s;  
	return Point(v[X] * d, -quaterHeight, v[Z] * d);
      } else  {
	return Point(0, -quaterHeight, 0);
      }
    }
  } else {
    return Point();
  }
}


// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction d'un clone d'un cone
Convex* Cone::clone() const 
{
  return new Cone(bottomRadius,4.0*quaterHeight);
}

// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le Volume d'un cone
Scalar Cone::getVolume() const
{
  return 4.0*quaterHeight*PI*bottomRadius*bottomRadius/3.0;
}

// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine l'inertie et l'inertie inverse d'un cone
bool Cone::BuildInertie(Scalar *inertie,Scalar *inertie_1) const 
{
  inertie[1]=inertie[2]=inertie[4]= 0.0;
  const Scalar constante = 0.2*quaterHeight*bottomRadius
    *bottomRadius*PI;
  inertie[0]=inertie[5]= constante*
    (4*quaterHeight*quaterHeight + bottomRadius*bottomRadius);
  inertie[3] = 2. * constante * bottomRadius*bottomRadius;

  inertie_1[1]=inertie_1[2]=inertie_1[4]= 0.0;
  inertie_1[5]=inertie_1[0]= 1.0/inertie[0];
  inertie_1[3]=1.0/inertie[3];
  return true;
}


// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
//  Determine le rayon circonscrit d'un cone
Scalar Cone::BuildRayonRef() const 
{
  return max(sqrt(bottomRadius*bottomRadius + quaterHeight*quaterHeight)
	     , 3. * quaterHeight);
}


// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov. 2015 - Creation
//Determine le rayon circonscrit d'un cone retrecissant
int Cone::getShrinkingChoice() const
{
  return(0);
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov. 2015 - Creation
// Determine le rayon circonscrit d'un cone retrecissant
void Cone:: setShrinkingRadius(Scalar CurrentRadius)
{
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture d'un Cone
void Cone::printClass(ostream &fileOut) const 
{
  fileOut << "*Cone\n";
  fileOut << bottomRadius     << '\t' 
	  << 4.0*quaterHeight << '\n' ;
}


// ----------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// lecture d'un Cone
void Cone::readClass(istream &fileIn) 
{
  fileIn >> bottomRadius
	 >> quaterHeight;
  quaterHeight /= 4.0;
  sinAngle = bottomRadius / sqrt(bottomRadius * bottomRadius + 
				 quaterHeight * quaterHeight);
}
