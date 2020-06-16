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

#include "Cylinder.H"

int Cylinder::visuNodeNbOnPer = 36;

// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Cylinder::Cylinder(Scalar r , Scalar h ) : 
    radius(r), halfHeight(h / 2) 
{
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction avec un fichier
Cylinder::Cylinder(istream& fileIn) 
{
  readClass(fileIn);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructor
Cylinder::Cylinder(DOMNode* root)
{
  radius     = ReaderXML::getNodeAttr_Double(root, "Radius" );
  halfHeight = ReaderXML::getNodeAttr_Double(root, "Hauteur") / 2.;
}


  

// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Cylinder::~Cylinder() 
{
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine l'inertie et l'inertie inverse d'un cylindre
bool Cylinder::BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const 
{
  inertie[1] = inertie[2] = inertie[4] = 0.0;
  const Scalar constante = 0.5 * halfHeight * radius * radius * PI;
  inertie[0] = inertie[5] = 
    constante * ( 4.0 * halfHeight * halfHeight / 3.0 + radius * radius );
  inertie[3] = 2.0 * constante * radius * radius;

  inertie_1[1] = inertie_1[2] = inertie_1[4] = 0.0;
  inertie_1[5] = inertie_1[0] = 1.0 / inertie[0];
  inertie_1[3] = 1.0 / inertie[3];
  return true;
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// Renvoie le choix de retrecissement
int Cylinder:: getShrinkingChoice()const
 {

return (0);
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// set the shrinking radius
void Cylinder::setShrinkingRadius(Scalar CurrentRadius)
{

}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
//  Determine le rayon circonscrit d'un cylindre
Scalar Cylinder::BuildRayonRef() const 
{
  return sqrt(radius*radius+halfHeight*halfHeight);
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction d'un clone du cylindre
Convex* Cylinder::clone() const
{
  return new Cylinder(radius,2*halfHeight);
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le Volume d'un cylindre
Scalar Cylinder::getVolume() const
{
  return 2*halfHeight*PI*radius*radius;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Point Cylinder::support(const Vecteur& v) const 
{
  Scalar norm = Norm(v);
  if (norm > EPSILON) {
    Scalar s = sqrt(v[X] * v[X] + v[Z] * v[Z]);
    if (s > EPSILON) {
      Scalar d = radius / s;  
      return Point(v[X] * d, v[Y] < 0. ? -halfHeight : halfHeight, v[Z] * d);
    } else {
      return Point(0., v[Y] < 0. ? -halfHeight : halfHeight, 0.);
    }
  } else {
    return Point();
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Points decrivant le Convex pour visu dans une appli externe (Fluide)
// Pour un disque pas de prise en compte de points sommets -> 0
vector<Point> Cylinder::getEnveloppe() const
{
  Point point(0.,0.,0.);
  vector<Point> enveloppe(3,point);
  enveloppe[0][Y] = - halfHeight;
  enveloppe[1][Y] = - halfHeight;
  enveloppe[1][X] = radius;  
  enveloppe[2][Y] = halfHeight;    
  return enveloppe;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de sommets ou un code équivalent
int Cylinder::getNbCorners() const
{
  return 777;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Description des faces
vector< vector<int> > const* Cylinder::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return allFaces;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture d'un Cylindre
void Cylinder::printClass(ostream &fileOut) const 
{
  fileOut << "*Cylindre\n";
  fileOut << radius << '\t' 
	  << 2.0*halfHeight << '\n';
}



// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// lecture d'un Cylindre
void Cylinder::readClass(istream &fileIn) 
{
  fileIn >> radius 
	 >> halfHeight;
  halfHeight /= 2.0;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de points pour post-processing avec Paraview
int Cylinder::numberOfPoints_PARAVIEW() const
{
  return 2*visuNodeNbOnPer+2;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de polygones elementaires pour post-processing avec Paraview
int Cylinder::numberOfCells_PARAVIEW() const
{
  return visuNodeNbOnPer;  
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit les points du convexe pour post-processing avec Paraview  
void Cylinder::write_polygonsPts_PARAVIEW( ostream &f,
	const Transform &transform, Vecteur const* translation ) const
{
  Point pp,p;
  Scalar dtheta = 2.* PI / visuNodeNbOnPer;
  
  // Couronne inferieure
  p[Y] = -halfHeight;
  for (int i=0;i<visuNodeNbOnPer;++i)
  {
    p[X] = radius * cos ( i*dtheta );
    p[Z] = radius * sin ( i*dtheta ); 
    pp = transform(p);
    if ( translation ) pp += *translation;    
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }
  
  // Couronne superieure
  p[Y] = halfHeight;
  for (int i=0;i<visuNodeNbOnPer;++i)
  {
    p[X] = radius * cos ( i*dtheta );
    p[Z] = radius * sin ( i*dtheta ); 
    pp = transform(p);
    if ( translation ) pp += *translation;    
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }
  
  // Centre inferieur
  p[X] = 0.;
  p[Y] = -halfHeight;
  p[Z] = 0.;
  pp = transform(p);

  if ( translation ) pp += *translation;  
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  
  // Centre superieur
  p[Y] = halfHeight;
  pp = transform(p);
  if ( translation ) pp += *translation;    
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les points du convexe pour post-processing avec Paraview  
list<Point> Cylinder::get_polygonsPts_PARAVIEW( const Transform &transform,
	Vecteur const* translation ) const
{
  list<Point> ParaviewPoints;
  Point pp,p;
  Scalar dtheta = 2.* PI / visuNodeNbOnPer;
  
  // Couronne inferieure
  p[Y] = -halfHeight;  
  for (int i=0;i<visuNodeNbOnPer;++i)
  {
    p[X] = radius * cos ( i*dtheta );
    p[Z] = radius * sin ( i*dtheta ); 
    pp = transform(p);
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back(pp);
  }
  
  // Couronne superieure
  p[Y] = halfHeight;
  for (int i=0;i<visuNodeNbOnPer;++i)
  {
    p[X] = radius * cos ( i*dtheta );
    p[Z] = radius * sin ( i*dtheta ); 
    pp = transform(p);
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back(pp);
  }
  
  // Centre inferieur
  p[X] = 0.;
  p[Y] = -halfHeight;
  p[Z] = 0.;
  pp = transform(p);
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back(pp);
  
  // Centre superieur
  p[Y] = halfHeight;
  pp = transform(p);
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back(pp);
  
  return ParaviewPoints; 
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe pour post-processing avec Paraview 
void Cylinder::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  for (int i=0;i<visuNodeNbOnPer-1;++i)
  {
    connectivity.push_back(firstpoint_globalnumber+i);
    connectivity.push_back(firstpoint_globalnumber+i+1);    
    connectivity.push_back(firstpoint_globalnumber+2*visuNodeNbOnPer);    
    connectivity.push_back(firstpoint_globalnumber+i+visuNodeNbOnPer);    
    connectivity.push_back(firstpoint_globalnumber+i+visuNodeNbOnPer+1);
    connectivity.push_back(firstpoint_globalnumber+2*visuNodeNbOnPer+1);        
    last_offset+=6;    
    offsets.push_back(last_offset);
    cellstype.push_back(13);
  }
  connectivity.push_back(firstpoint_globalnumber+visuNodeNbOnPer-1);
  connectivity.push_back(firstpoint_globalnumber);    
  connectivity.push_back(firstpoint_globalnumber+2*visuNodeNbOnPer);    
  connectivity.push_back(firstpoint_globalnumber+2*visuNodeNbOnPer-1);    
  connectivity.push_back(firstpoint_globalnumber+visuNodeNbOnPer);
  connectivity.push_back(firstpoint_globalnumber+2*visuNodeNbOnPer+1);        
  last_offset+=6;    
  offsets.push_back(last_offset);
  cellstype.push_back(13);  
  
  firstpoint_globalnumber+=2*visuNodeNbOnPer+2;
}
