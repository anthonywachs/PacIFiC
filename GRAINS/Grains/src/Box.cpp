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
#include "Grains_Exec.hh"
#include "Grains_BuilderFactory.H"
#include "Box.H"
#include <sstream>

// Attribut statique
vector< vector<int> > Box::allFaces;


// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Box::Box(const Vecteur extent_) 
  : extent(extent_),
  corners2D_XY( NULL )
{  
  setCornersFaces();
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Box::Box(Scalar x, Scalar y, Scalar z) 
  : extent(fabs(x)/2., fabs(y)/2., fabs(z)/2.),
  corners2D_XY( NULL )
{
  setCornersFaces();
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction avec un fichier
Box::Box(istream &fileIn)
  : corners2D_XY( NULL )
{
  readClass(fileIn);
  setCornersFaces();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
Box::Box(DOMNode* root)
  : corners2D_XY( NULL )
{
  extent[X] = ReaderXML::getNodeAttr_Double(root, "LX") / 2.;
  extent[Y] = ReaderXML::getNodeAttr_Double(root, "LY") / 2.;
  extent[Z] = ReaderXML::getNodeAttr_Double(root, "LZ") / 2.;
  setCornersFaces();
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Box::~Box() 
{
  corners.clear();
  if ( corners2D_XY )
  {
    corners2D_XY->clear();
    delete corners2D_XY;    
  } 
}




// ----------------------------------------------------------------------------
// Construit les sommets du pavé et les faces
void Box::setCornersFaces()
{
  corners.reserve(8);
  Point sommet;
  sommet.setValue(-extent[X],-extent[Y],-extent[Z]);
  corners.push_back(sommet);
  sommet.setValue(extent[X],-extent[Y],-extent[Z]);
  corners.push_back(sommet);  
  sommet.setValue(extent[X],-extent[Y],extent[Z]);
  corners.push_back(sommet);  
  sommet.setValue(-extent[X],-extent[Y],extent[Z]);
  corners.push_back(sommet);  
  sommet.setValue(-extent[X],extent[Y],-extent[Z]);
  corners.push_back(sommet);  
  sommet.setValue(extent[X],extent[Y],-extent[Z]);
  corners.push_back(sommet);  
  sommet.setValue(extent[X],extent[Y],extent[Z]);
  corners.push_back(sommet);  
  sommet.setValue(-extent[X],extent[Y],extent[Z]);
  corners.push_back(sommet);
  
  if (allFaces.empty())
  {
    vector<int> oneFace(4,0);
    allFaces.reserve(6);
    for (int i=0;i<6;++i) allFaces.push_back(oneFace);

    allFaces[0][0]=4;
    allFaces[0][1]=5;    
    allFaces[0][2]=6;    
    allFaces[0][3]=7;    
    
    allFaces[1][0]=5;
    allFaces[1][1]=1;    
    allFaces[1][2]=2;    
    allFaces[1][3]=6;        
    
    allFaces[2][0]=1;
    allFaces[2][1]=0;    
    allFaces[2][2]=3;    
    allFaces[2][3]=2;        
        
    allFaces[3][0]=0;
    allFaces[3][1]=4;    
    allFaces[3][2]=7;    
    allFaces[3][3]=3;        
    
    allFaces[4][0]=4;
    allFaces[4][1]=0;    
    allFaces[4][2]=1;    
    allFaces[4][3]=5;        
    
    allFaces[5][0]=3;
    allFaces[5][1]=7;    
    allFaces[5][2]=6;    
    allFaces[5][3]=2;        
  }
  
  // Dans le cas où la box est utlisée dans une simu 2D
  if ( Grains_BuilderFactory::getContext() == DIM_2 )
  {
    corners2D_XY = new vector<Point>(4,sommet);
    (*corners2D_XY)[0] = corners[0] ;
    (*corners2D_XY)[0][Z] = 0.;
    (*corners2D_XY)[1] = corners[1] ;
    (*corners2D_XY)[1][Z] = 0.;
    (*corners2D_XY)[2] = corners[5] ;
    (*corners2D_XY)[2][Z] = 0.;
    (*corners2D_XY)[3] = corners[7] ;
    (*corners2D_XY)[3][Z] = 0.;   
  }
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// Renvoie le choix de retrecissement
int Box:: getShrinkingChoice()const
{
  return (0);
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// set the shrinking radius
void Box:: setShrinkingRadius(Scalar CurrentRadius)
{

}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction d'un clone d'un box
Convex* Box::clone() const
{
  return new Box(2.0*extent[X],
		 2.0*extent[Y],
		 2.0*extent[Z]);
}




// ----------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine l'inertie et l'inertie inverse d'un box
bool Box::BuildInertie(Scalar *inertie,Scalar *inertie_1) const
{
  inertie[1] = inertie[2] = inertie[4] = 0.0;
  inertie[0] = 8.0*extent[X]*extent[Y]*extent[Z]
    *(extent[Y]*extent[Y]+extent[Z]*extent[Z])/3.0;
  inertie[3] = 8.0*extent[X]*extent[Y]*extent[Z]
    *(extent[X]*extent[X]+extent[Z]*extent[Z])/3.0;
  inertie[5] = 8.0*extent[X]*extent[Y]*extent[Z]
    *(extent[Y]*extent[Y]+extent[X]*extent[X])/3.0;

  inertie_1[1] = inertie_1[2] = inertie_1[4] = 0.0;
  inertie_1[0] = 1.0/inertie[0];
  inertie_1[3] = 1.0/inertie[3];
  inertie_1[5] = 1.0/inertie[5];
  return true;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le rayon circonscrit d'une boite
Scalar Box::BuildRayonRef() const 
{
  return Norm(extent);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Points decrivant l'enveloppe
vector<Point> Box::getEnveloppe() const
{
  if ( Grains_BuilderFactory::getContext() == DIM_2 ) return *corners2D_XY;
  else return corners;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Description des faces
vector<vector<int> > const* Box::getFaces() const
{
  return &allFaces;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de sommets ou un code équivalent
int Box::getNbCorners() const
{
  if ( Grains_BuilderFactory::getContext() == DIM_2 ) return 444;  
  else return 666;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le Volume d'un box
Scalar Box::getVolume() const
{
  if ( Grains_BuilderFactory::getContext() == DIM_2 ) 
    return 4.0*extent[X]*extent[Y];  
  else return 8.0*extent[X]*extent[Y]*extent[Z];
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture d'un Box
void Box::printClass(ostream &fileOut) const 
{
  fileOut << "*Box\n"
	  << extent;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// lecture d'un Box
void Box::readClass(istream &fileIn) 
{
  fileIn >> extent;
}




// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
Point Box::support(const Vecteur& v) const 
{
  Scalar norm = Norm(v);
  if (norm < EPSILON) {
    return Point();
  } else {
    return Point(v[X] < 0 ? -extent[X] : extent[X],
		 v[Y] < 0 ? -extent[Y] : extent[Y],
		 v[Z] < 0 ? -extent[Z] : extent[Z]); 
  }
}




// ----------------------------------------------------------------------------
// Ecriture de l'objet dans sa configuration courante pour
// post-processing avec GMV
void Box::GMVoutput(ostream &fileOut,const Transform &transform) const
{
  int i,j,coor;
  
  for (i=0;i<6;++i)
  {
    vector<Point> PointFace;
    PointFace.reserve(5);
    for (j=0;j<4;++j) 
      PointFace.push_back(transform(corners[allFaces[i][j]]));  
    PointFace.push_back(transform(corners[allFaces[i][0]]));
    fileOut << "2 5";
    for (coor=0;coor<3;++coor)
      for (j=0;j<5;++j)
        fileOut << " " << PointFace[j][coor];
    fileOut << endl;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de points pour post-processing avec Paraview
int Box::numberOfPoints_PARAVIEW() const
{
  return 8;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de polygones elementaires pour post-processing avec Paraview
int Box::numberOfCells_PARAVIEW() const
{
  return 1;  
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit les points du convexe pour post-processing avec Paraview  
void Box::write_polygonsPts_PARAVIEW( ostream &f,
	const Transform &transform, Vecteur const* translation ) const
{
  Point pp;
  for (int i=0;i<8;++i)
  {
    pp = transform(corners[i]);
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les points du convexe pour post-processing avec Paraview  
list<Point> Box::get_polygonsPts_PARAVIEW( const Transform &transform,
	Vecteur const* translation ) const
{
  list<Point> ParaviewPoints;
  Point pp;
  for (int i=0;i<8;++i)
  {
    pp = transform(corners[i]);
    if ( translation ) pp += *translation;    
    ParaviewPoints.push_back(pp);
  }
  
  return ParaviewPoints; 
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe pour post-processing avec Paraview 
void Box::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  int count=firstpoint_globalnumber;
  for (int i=0;i<8;++i)
  {
    connectivity.push_back(count);
    ++count;
  }
  last_offset+=8;    
  offsets.push_back(last_offset);
  cellstype.push_back(12);
  
  firstpoint_globalnumber+=8;
}




// ----------------------------------------------------------------------------
// Renvoi le point[0] de la face i
Point Box::getFirstPointFace(int i) const
{
  return corners[allFaces[i][0]];
}




// ----------------------------------------------------------------------------
// Renvoie le point de contact dans l'espace du convexe et la distance
// de recouvrement entre la boite et une sphere. Si le contact n'existe pas,
// renvoie le point origine (0,0,0)
Point Box::IntersectionPointSPHERE(const Point& SphereCenter,
  	const Scalar& SphereRadius,
	double &overlap) const
  throw(ErreurContact)
{
  Point contactPoint;
  Scalar gx = SphereCenter[X], gy = SphereCenter[Y], gz = SphereCenter[Z];
  Vecteur distance;
  Scalar normDistance = 0.;
  overlap = 1.;
  
  if ( gx > extent[X] )
  {
    if ( gy > extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance au corner 6
	contactPoint = ContactCornerSPHERE(SphereCenter,SphereRadius,6,overlap);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance au corner 5
        contactPoint = ContactCornerSPHERE(SphereCenter,SphereRadius,5,overlap);
      }
      else
      {
        // Distance a l'arete 5-6
        contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,5,Z,overlap);
      }
    }
    else if ( gy < -extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance au corner 2
	contactPoint = ContactCornerSPHERE(SphereCenter,SphereRadius,2,overlap);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance au corner 1
	contactPoint = ContactCornerSPHERE(SphereCenter,SphereRadius,1,overlap);
      }
      else
      {
        // Distance a l'arete 1-2
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,1,Z,overlap);
      }    
    }
    else
    {
      if ( gz > extent[Z] )
      {
        // Distance a l'arete 2-6
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,6,Y,overlap);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance a l'arete 1-5
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,1,Y,overlap);
      }
      else
      {
        // Distance a la face 1
	normDistance = gx - extent[X];
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue(extent[X]+0.5*overlap,SphereCenter[Y],
	  	SphereCenter[Z]);
        }   
      }    
    }
  }
  else if ( gx < -extent[X] )
  {
    if ( gy > extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance au corner 7
	contactPoint = ContactCornerSPHERE(SphereCenter,SphereRadius,7,overlap);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance au corner 4
	contactPoint = ContactCornerSPHERE(SphereCenter,SphereRadius,4,overlap);
      }
      else
      {
        // Distance a l'arete 4-7
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,4,Z,overlap);
      }
    }
    else if ( gy < -extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance au corner 3
	contactPoint = ContactCornerSPHERE(SphereCenter,SphereRadius,3,overlap);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance au corner 0
	contactPoint = ContactCornerSPHERE(SphereCenter,SphereRadius,0,overlap);
      }
      else
      {
        // Distance a l'arete 0-3
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,0,Z,overlap);
      }    
    }
    else
    {
      if ( gz > extent[Z] )
      {
        // Distance a l'arete 3-7
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,3,Y,overlap);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance a l'arete 0-4
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,0,Y,overlap);
      }
      else
      {
        // Distance a la face 3
	normDistance = - extent[X] - gx;
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue(-extent[X]-0.5*overlap,SphereCenter[Y],
	  	SphereCenter[Z]);
        }   	       
      }    
    }
  }
  else
  {
    if ( gy > extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance a l'arete 6-7
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,6,X,overlap);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance a l'arete 4-5
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,4,X,overlap);
      }
      else
      {
        // Distance a la face 0
	normDistance = gy - extent[Y];
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue(SphereCenter[X],extent[Y]+0.5*overlap,
	  	SphereCenter[Z]);
        } 	       
      }
    }
    else if ( gy < -extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance a l'arete 3-2
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,3,X,overlap);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance a l'arete 0-1
	contactPoint = ContactEdgeSPHERE(SphereCenter,SphereRadius,0,X,overlap);
      }
      else
      {
        // Distance a la face 2
	normDistance = - extent[Y] - gy;
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue(SphereCenter[X],-extent[Y]-0.5*overlap,
	  	SphereCenter[Z]);
        } 	         
      }    
    }
    else
    {
      if ( gz > extent[Z] )
      {
        // Distance a la face 5
	normDistance = gz - extent[Z];
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue(SphereCenter[X],SphereCenter[Y],
	  	extent[Z]+0.5*overlap);
        } 	         
      }
      else if ( gz < -extent[Z] )
      {
        // Distance a la face 4
	normDistance = - extent[Z] - gz;
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue(SphereCenter[X],SphereCenter[Y],
	  	-extent[Z]-0.5*overlap);
        } 		          
      }
      else
      {
        cout << "Warning: sphere center in box in Box::IntersectionPointSPHERE"
		<< endl;
	Grains_Exec::m_exception_Contact = true;
        throw ErreurContact();
      }    
    }
  }
  
  return contactPoint;

}  




// ----------------------------------------------------------------------------
// Renvoie le point de contact dans l'espace du convexe et la distance
// de recouvrement entre un corner de la boite et une sphere. Si le contact 
// n'existe pas, renvoie le point origine (0,0,0)
Point Box::ContactCornerSPHERE(const Point& SphereCenter,
  	const Scalar& SphereRadius,
	const int &cornerNumber,
	double &overlap) const
{
  Point contactPoint;  
  Vecteur distance = SphereCenter - corners[cornerNumber];
  Scalar normDistance = Norm(distance);
  if ( normDistance < SphereRadius )
  {
    overlap = normDistance - SphereRadius;
    contactPoint = SphereCenter - (1.+0.5*overlap/normDistance) * distance;
  }
	
  return contactPoint;	
} 




// ----------------------------------------------------------------------------
// Renvoie le point de contact dans l'espace du convexe et la distance
// de recouvrement entre une arete de la boite et une sphere. Si le contact 
// n'existe pas, renvoie le point origine (0,0,0)
Point Box::ContactEdgeSPHERE(const Point& SphereCenter,
  	const Scalar& SphereRadius,
	const int &cornerNumber,
	const int &projectionDirection,
	double &overlap) const
{
  Point contactPoint,SphereCenterProjected(SphereCenter); 
  SphereCenterProjected[projectionDirection] = 
  	corners[cornerNumber][projectionDirection];
  Vecteur distance = SphereCenterProjected - corners[cornerNumber];
  Scalar normDistance = Norm(distance);
  if ( normDistance < SphereRadius )
  {
    overlap = normDistance - SphereRadius;
    contactPoint = SphereCenter - (1.+0.5*overlap/normDistance) * distance;
  }
	
  return contactPoint;	
} 




// ----------------------------------------------------------------------------
// Renvoie le point d'intersection avec la coque du polyèdre d'un
// segment défini par un point interne et un point externe
Point Box::intersectionToShell(const Point &PtIN, const Point &PtOut) const
{
  Scalar parameter=0.;
  bool b_found = false;
  Vecteur InOut = PtOut - PtIN;
  Point OnFace;

  // Face 0 : y=Ly/2
  if ( !b_found )
  {
    if ( fabs(InOut[Y]) > EPSILON )
    {
      parameter = (extent[Y] - PtIN[Y])/InOut[Y];
      OnFace[X] = PtIN[X]+parameter*InOut[X];
      OnFace[Y] = extent[Y];  
      OnFace[Z] = PtIN[Z]+parameter*InOut[Z]; 
      if ( OnFace[X] <= extent[X] && OnFace[X] >= -extent[X] 
  	&& OnFace[Z] <= extent[Z] && OnFace[Z] >= -extent[Z] 
	&& PtOut[Y] >= extent[Y]) b_found = true;
    }
  }

  // Face 1 : x=Lx/2
  if ( !b_found )
  {
    if ( fabs(InOut[X]) > EPSILON )
    {
      parameter = (extent[X] - PtIN[X])/InOut[X];
      OnFace[X] = extent[X];
      OnFace[Y] = PtIN[Y]+parameter*InOut[Y];  
      OnFace[Z] = PtIN[Z]+parameter*InOut[Z]; 
      if ( OnFace[Y] <= extent[Y] && OnFace[Y] >= -extent[Y]
  	&& OnFace[Z] <= extent[Z] && OnFace[Z] >= -extent[Z]
	&& PtOut[X] >= extent[X]) b_found = true;
    }	
  }

  // Face 2 : y=-Ly/2
  if ( !b_found )
  {
    if ( fabs(InOut[Y]) > EPSILON )
    {
      parameter = (-extent[Y] - PtIN[Y])/InOut[Y];
      OnFace[X] = PtIN[X]+parameter*InOut[X];
      OnFace[Y] = -extent[Y];  
      OnFace[Z] = PtIN[Z]+parameter*InOut[Z]; 
      if ( OnFace[X] <= extent[X] && OnFace[X] >= -extent[X] 
  	&& OnFace[Z] <= extent[Z] && OnFace[Z] >= -extent[Z]
	&& PtOut[Y] <= -extent[Y]) b_found = true;
    }	
  }

  // Face 3 : x=-Lx/2
  if ( !b_found )
  {
    if ( fabs(InOut[X]) > EPSILON )
    {
      parameter = (-extent[X] - PtIN[X])/InOut[X];
      OnFace[X] = -extent[X];
      OnFace[Y] = PtIN[Y]+parameter*InOut[Y];  
      OnFace[Z] = PtIN[Z]+parameter*InOut[Z]; 
      if ( OnFace[Y] <= extent[Y] && OnFace[Y] >= -extent[Y]
  	&& OnFace[Z] <= extent[Z] && OnFace[Z] >= -extent[Z]
	&& PtOut[X] <= -extent[X]) b_found = true;
    }	
  }

  // Face 4 : z=-Lz/2
  if ( !b_found )
  {
    if ( fabs(InOut[Z]) > EPSILON )
    {
      parameter = (-extent[Z] - PtIN[Z])/InOut[Z];
      OnFace[X] = PtIN[X]+parameter*InOut[X];
      OnFace[Y] = PtIN[Y]+parameter*InOut[Y];  
      OnFace[Z] = -extent[Z]; 
      if ( OnFace[X] <= extent[X] && OnFace[X] >= -extent[X] 
  	&& OnFace[Y] <= extent[Y] && OnFace[Y] >= -extent[Y]
	&& PtOut[Z] <= -extent[Z]) b_found = true;
    }		
  }

  // Face 5 : z=Lz/2
  if ( !b_found )
  {
    if ( fabs(InOut[Z]) > EPSILON )
    {
      parameter = (extent[Z] - PtIN[Z])/InOut[Z];
      OnFace[X] = PtIN[X]+parameter*InOut[X];
      OnFace[Y] = PtIN[Y]+parameter*InOut[Y];  
      OnFace[Z] = extent[Z]; 
      if ( OnFace[X] <= extent[X] && OnFace[X] >= -extent[X] 
  	&& OnFace[Y] <= extent[Y] && OnFace[Y] >= -extent[Y]
	&& PtOut[Z] >= extent[Z]) b_found = true;
    }	
  }  		  

  if ( !b_found ) 
    cout << "!!! Box::intersectionToShell not found !!!" << endl;
  
  return OnFace;
}




// ----------------------------------------------------------------------------
// This method sends back the projection of sphere center on obstacle's faces, 
// edges or corners (based on local position of the sphere with respect to the 
// obstacle). 
// Pay attention that only in the case of projection on faces, the gap (normal 
// distance between sphere and the wall) is  modified. In case of projection on 
// edges and corners, the gap is set to zero ( can be modified later) */ 
Point Box::ProjectedPointSPHERE( const Point& SphereCenter,
  	const Scalar& SphereRadius,
	double &gap ) const
{
  Point ProjectedPoint;
  Scalar gx = SphereCenter[X], gy = SphereCenter[Y], gz = SphereCenter[Z];
  Scalar normDistance = 0.;
  gap = 0.;  
  
  if ( gx > extent[X] )
  {
    if ( gy > extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance au corner 6
	ProjectedPoint = corners[6];
      }
      else if ( gz < -extent[Z] )
      {
        // Distance au corner 5
        ProjectedPoint = (corners[5]);
      }
      else
      {
        // Distance a l'arete 5-6
        ProjectedPoint = SphereCenter; 
        ProjectedPoint[Z] = corners[5][Z];
      }
    }
    else if ( gy < -extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance au corner 2
	ProjectedPoint = (corners[2]);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance au corner 1
	ProjectedPoint = (corners[1]);
      }
      else
      {
        // Distance a l'arete 1-2
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[Z] = corners[1][Z];
      }    
    }
    else
    {
      if ( gz > extent[Z] )
      {
        // Distance a l'arete 2-6
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[Y] = corners[6][Y];
      }
      else if ( gz < -extent[Z] )
      {
        // Distance a l'arete 1-5
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[Y] = corners[1][Y];
      }
      else
      {
        // Distance a la face 1
	normDistance = gx - extent[X];
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue(extent[X],SphereCenter[Y],
	  	SphereCenter[Z]);
        }   
      }    
    }
  }
  else if ( gx < -extent[X] )
  {
    if ( gy > extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance au corner 7
	ProjectedPoint = (corners[7]);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance au corner 4
	ProjectedPoint = (corners[4]);
      }
      else
      {
        // Distance a l'arete 4-7
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[Z] = corners[4][Z];
      }
    }
    else if ( gy < -extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance au corner 3
	ProjectedPoint = (corners[3]);
      }
      else if ( gz < -extent[Z] )
      {
        // Distance au corner 0
	ProjectedPoint = (corners[0]);
      }
      else
      {
        // Distance a l'arete 0-3
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[Z] = corners[0][Z];
      }    
    }
    else
    {
      if ( gz > extent[Z] )
      {
        // Distance a l'arete 3-7
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[Y] = corners[3][Y];
      }
      else if ( gz < -extent[Z] )
      {
        // Distance a l'arete 0-4
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[Y] = corners[0][Y];
      }
      else
      {
        // Distance a la face 3
	normDistance = - extent[X] - gx;
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue(-extent[X],SphereCenter[Y],
	  	SphereCenter[Z]);
        }   	       
      }    
    }
  }
  else
  {
    if ( gy > extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance a l'arete 6-7
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[X] = corners[6][X];
      }
      else if ( gz < -extent[Z] )
      {
        // Distance a l'arete 4-5
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[X] = corners[4][X];
      }
      else
      {
        // Distance a la face 0
	normDistance = gy - extent[Y];
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue(SphereCenter[X],extent[Y],
	  	SphereCenter[Z]);
        } 	       
      }
    }
    else if ( gy < -extent[Y] )
    {
      if ( gz > extent[Z] )
      {
        // Distance a l'arete 3-2
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[X] = corners[3][X];
      }
      else if ( gz < -extent[Z] )
      {
        // Distance a l'arete 0-1
	ProjectedPoint = SphereCenter; 
        ProjectedPoint[X] = corners[0][X];
      }
      else
      {
        // Distance a la face 2
	normDistance = - extent[Y] - gy;
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue(SphereCenter[X],-extent[Y],
	  	SphereCenter[Z]);
        } 	         
      }    
    }
    else
    {
      if ( gz > extent[Z] )
      {
        // Distance a la face 5
	normDistance = gz - extent[Z];
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue(SphereCenter[X],SphereCenter[Y],
	  	extent[Z]);
        } 	         
      }
      else if ( gz < -extent[Z] )
      {
        // Distance a la face 4
	normDistance = - extent[Z] - gz;
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue(SphereCenter[X],SphereCenter[Y],
	  	-extent[Z]);
        } 		          
      }
      else
      {
        cout << "Warning: sphere center in box in Box::ProjectedPointSPHERE"
		<< endl;
      }    
    }
  }
  
  return ProjectedPoint;

}
