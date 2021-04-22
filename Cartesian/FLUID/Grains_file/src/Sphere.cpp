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
#include "Sphere.H"
#include "sstream"


int Sphere::visuNodeNbPerQar = 8;

// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// M. SULAIMAN - Nov.2015 - MODIF
Sphere::Sphere( Scalar r ) :
  radius(r),
  Shrinking(0),
  initial_radius(r)
{
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Constructeur avec fichier
Sphere::Sphere( istream &fileIn ) :
  Shrinking(0),
  initial_radius(0)
{
  readClass(fileIn);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
Sphere::Sphere( DOMNode* root ) :
  Shrinking(0)
{
  initial_radius = ReaderXML::getNodeAttr_Double(root, "Radius");
  if ( ReaderXML::hasNodeAttr_String( root, "Shrink" ) )
  Shrinking = ReaderXML::getNodeAttr_Int(root,"Shrink");
  radius = initial_radius;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Sphere::~Sphere() 
{
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine l'inertie et l'inertie inverse d'une sphere 
bool Sphere::BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const
{
  inertie[1]=inertie[2]=inertie[4]= 0.0;
  inertie[5]=inertie[3]
    =inertie[0]= 8.0*PI*radius*radius*radius*radius*radius/15.0;
  
  inertie_1[1]=inertie_1[2]=inertie_1[4]= 0.0;
  inertie_1[5]=inertie_1[3]=inertie_1[0]= 
    1.0/inertie[0];
  return true;
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov. 2015 - Creation
// Determine le rayon circonscrit a la sphere retrecissante
void Sphere::setShrinkingRadius(Scalar CurrentRadius)
{
  radius = CurrentRadius;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le rayon circonscrit a la sphere
Scalar Sphere::BuildRayonRef() const 
{
  return radius;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction d'un clone de la sphere
Convex* Sphere::clone() const
{
  return new Sphere(radius);
}




// ----------------------------------------------------------------------------
Point Sphere::support( const Vecteur& v ) const 
{
  Scalar s = Norm(v);
  if (s > EPSILON) {
    Scalar r = radius / s;
    return Point(v[X] * r, v[Y] * r, v[Z] * r);
  } else 
    return Point();    
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Points decrivant le Convex pour visu dans une appli externe (Fluide)
// Pour un disque pas de prise en compte de points sommets -> 0
vector<Point> Sphere::getEnveloppe() const
{
  vector<Point> enveloppe;
  Point point(0.,0.,0.);
  enveloppe.push_back(point);
  return enveloppe;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de sommets ou un code équivalent
int Sphere::getNbCorners() const
{
  return 1;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Description des faces
vector< vector<int> > const* Sphere::getFaces() const
{
  vector< vector<int> > const* allFaces = NULL;
  return allFaces;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le Volume d'une sphere
Scalar Sphere::getVolume() const
{
  return (4.0 * PI *radius*radius*radius / 3.0);
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov. 2015 - Creation
// Renvoie le choix si la particule est retrecissante
int Sphere:: getShrinkingChoice() const
{
  return Shrinking;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// M.SULAIMAN - Nov. 2015 - Modification
// ecriture de la sphere
void Sphere::printClass( ostream &fileOut ) const 
{
  fileOut << "*Sphere\n";
  
  if(Shrinking!=1)
    fileOut << radius <<endl; 
  else
  {
    fileOut << radius;
    fileOut << " " << Shrinking <<'\n';
  }
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// M.SULAIMAN - Nov. 2015 - Modification
// lecture de la sphere
void Sphere::readClass( istream &fileIn ) 
{
  fileIn >> radius;
  if ( Shrinking == 1 )
    fileIn >> Shrinking; 
}




// ----------------------------------------------------------------------------
// Ecriture de l'objet dans sa configuration courante pour
// post-processing avec GMV
void Sphere::GMVoutput( ostream &fileOut, const Transform &transform ) const
{
  const int nq = visuNodeNbPerQar;
  Point ptcourant;
  vector<Point> vecp(2*nq+1,ptcourant);
  vector< vector<Point> > pts(4*nq+1,vecp);
  double angle = PI / (2. * nq),theta,beta,z,rxy,x,y;
  int i,j;
  
  // Calcul des coordonnées de pts 
  // La sphère est vue comme un tableau 2D de points de type Point 
  // (projection d'une sphère ouverte sur un plan, type Mercator)
  // La 1ère et la dernière ligne correspondent aux 2 poles de la sphère
  // Pour ces 2 lignes du tableau, les coordonnées x,y,z sont les mêmes
  // mais cela simplifie la définition des facettes
  for (j=0;j<2*nq+1;++j)
  {
    theta = PI - j * angle;
    z = radius * cos(theta);
    rxy = radius * sin(theta);
    for (i=0;i<4*nq+1;++i)
    {
      beta = i * angle; 
      x = rxy * cos(beta);
      y = rxy * sin(beta);
      ptcourant.setValue(x,y,z);
      pts[i][j] = transform(ptcourant);
    }
  }

  // Ecriture de la sphère décomposée en facettes au format GMV
  // Chaque facette est un polygone à 4 côtés
  for (i=0;i<4*nq;++i)
  {
    for (j=0;j<2*nq;++j)
    {
      // 1 est le type de matériau et 5 le nombre de points de la facette
      // le 1er et le dernier point étant le même
      fileOut << "1 5 ";
      
      // Coordonnées X des pts de la facette
      fileOut << pts[i][j][X] << " " << pts[i+1][j][X] << " " 
      	<< pts[i+1][j+1][X] << " " << pts[i][j+1][X] << " " 
	<< pts[i][j][X] << " ";
	
      // Coordonnées Y des pts de la facette
      fileOut << pts[i][j][Y] << " " << pts[i+1][j][Y] << " " 
      	<< pts[i+1][j+1][Y] << " " << pts[i][j+1][Y] << " " 
	<< pts[i][j][Y] << " ";
		
      // Coordonnées X des pts de la facette
      fileOut << pts[i][j][Z] << " " << pts[i+1][j][Z] << " " 
      	<< pts[i+1][j+1][Z] << " " << pts[i][j+1][Z] << " " 
	<< pts[i][j][Z] << endl;		
    }
  }
   
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview 
int Sphere::numberOfPoints_PARAVIEW() const 
{ 
  return (4*visuNodeNbPerQar*(2*visuNodeNbPerQar-1)+3);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit les points du convexe pour post-processing avec Paraview  
void Sphere::write_polygonsPts_PARAVIEW( ostream &f,
	const Transform &transform, Vecteur const* translation ) const
{
  double pi=acos(-1.);
  double angle = pi/(2.*visuNodeNbPerQar) ;
  double angleZ = 0., local_radius = 0.;
  int k,i,ptsPerlevel = 4*visuNodeNbPerQar;
  Point pp,pptrans;
  
  // Regular points on the surface
  for ( k = 0; k < 2*visuNodeNbPerQar-1 ; ++k ) 
  {  
    angleZ = -pi/2.+(k+1)*angle;
    local_radius = radius*cos(angleZ);
    pp[Z] = radius*sin(angleZ);
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius*cos(i*angle);
      pp[Y] = local_radius*sin(i*angle);
      pptrans = transform(pp);
      if ( translation ) pptrans += *translation;
      f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
    }
  }
   
  pp[X] = 0.;
  pp[Y] = 0.;
  // Bottom point
  pp[Z] = -radius;
  pptrans = transform(pp);
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
	
  // Top point
  pp[Z] = radius;
  pptrans = transform(pp);
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
	
  // Gravity center
  pp[Z] = 0.;
  pptrans = transform(pp);
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe au format STL pour lien avec openFoam 
void Sphere::write_convex_STL( ostream &f,const Transform &transform ) 
	const
{
  int visuNodeNbPerQar_STL = 16;
  double pi=acos(-1.);
  double angle = pi/(2.*visuNodeNbPerQar_STL) ;
  double angleZ = 0., local_radius = 0.;
  int k,i,ptsPerlevel = 4*visuNodeNbPerQar_STL, 
  	nbLevels = 2*visuNodeNbPerQar_STL-1;
  Point pp,pptrans,ppTop,ppBottom,GC;
  vector<Point> work(ptsPerlevel,pp);
  vector< vector<Point> > STLPoints(nbLevels,work);  
  
  // Regular points on the surface
  for ( k = 0; k < nbLevels ; ++k ) 
  {  
    angleZ = -pi/2.+(k+1)*angle;
    local_radius = radius*cos(angleZ);
    pp[Z] = radius*sin(angleZ);
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius*cos(i*angle);
      pp[Y] = local_radius*sin(i*angle);
      pptrans = transform(pp);
      STLPoints[k][i] = pptrans;
    }
  }
   
  pp[X] = 0.;
  pp[Y] = 0.;
  // Bottom point
  pp[Z] = -radius;
  ppBottom = transform(pp);
	
  // Top point
  pp[Z] = radius;
  ppTop = transform(pp);
	
  // Gravity center
  pp[Z] = 0.;
  GC = transform(pp);  
  
  
  // Writing facets
  // Regular facets
  for ( k = 0; k < nbLevels-1 ; ++k ) 
  {
    for ( i = 0; i < ptsPerlevel-1 ; ++i )
    {
      // Bottom triangle
      write_STLfacet_sphere(f,GC,STLPoints[k][i],STLPoints[k][i+1],
      	STLPoints[k+1][i]);
            
      // Top triangle
      write_STLfacet_sphere(f,GC,STLPoints[k][i+1],STLPoints[k+1][i+1],
      	STLPoints[k+1][i]);
    }
    
    // Last bottom triangle
    write_STLfacet_sphere(f,GC,STLPoints[k][ptsPerlevel-1],STLPoints[k][0],
      	STLPoints[k+1][ptsPerlevel-1]); 
	
    // Last top triangle
    write_STLfacet_sphere(f,GC,STLPoints[k][0],STLPoints[k+1][0],
      	STLPoints[k+1][ptsPerlevel-1]);	 
  }

  // Bottom facets
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
    write_STLfacet_sphere(f,GC,ppBottom,STLPoints[0][i],STLPoints[0][i+1]);
  write_STLfacet_sphere(f,GC,ppBottom,STLPoints[0][ptsPerlevel-1],
  	STLPoints[0][0]);  

  // Top facets
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
    write_STLfacet_sphere(f,GC,ppTop,STLPoints[nbLevels-1][i],
    	STLPoints[nbLevels-1][i+1]);
  write_STLfacet_sphere(f,GC,ppTop,STLPoints[nbLevels-1][ptsPerlevel-1],
    	STLPoints[nbLevels-1][0]);	  
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture d'une facette triangulaire au format STL à partir 
// des 3 points qui la constituent et du centre de gravite de la sphere
void Sphere::write_STLfacet_sphere( ostream &f, Point const& GC,
  	Point const& pp1,
  	Point const& pp2,	
  	Point const& pp3 ) const
{  
  Point triangleGC = pp1 + (2./3.) * ( 0.5*(pp2+pp3) - pp1 );
  Vecteur outward_normal = triangleGC - GC;
  outward_normal.normalize();
  f << "  facet normal " << outward_normal[X]
      	<< " " << outward_normal[Y]
      	<< " " << outward_normal[Z] << endl;
  f << "    outer loop" << endl;
  f << "      vertex " << pp1[X] << " " << pp1[Y] << " " << pp1[Z] << endl;
  f << "      vertex " << pp2[X] << " " << pp2[Y] << " " << pp2[Z] << endl;
  f << "      vertex " << pp3[X] << " " << pp3[Y] << " " << pp3[Z] << endl;
  f << "    endloop" << endl;
  f << "  endfacet" << endl;     
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les points du convexe pour post-processing avec Paraview  
list<Point> Sphere::get_polygonsPts_PARAVIEW( const Transform &transform,
	Vecteur const* translation ) const
{
  list<Point> ParaviewPoints;
  double pi=acos(-1.);
  double angle = pi/(2.*visuNodeNbPerQar) ;
  double angleZ = 0., local_radius = 0.;
  int k,i,ptsPerlevel = 4*visuNodeNbPerQar;
  Point pp,pptrans;
  
  // Regular points on the surface
  for ( k = 0; k < 2*visuNodeNbPerQar-1 ; ++k ) 
  {  
    angleZ = -pi/2.+(k+1)*angle;
    local_radius = radius*cos(angleZ);
    pp[Z] = radius*sin(angleZ);
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius*cos(i*angle);
      pp[Y] = local_radius*sin(i*angle);
      pptrans = transform(pp);
      if ( translation ) pptrans += *translation;
      ParaviewPoints.push_back(pptrans);
    }
  }
   
  pp[X] = 0.;
  pp[Y] = 0.;
  // Bottom point
  pp[Z] = -radius;
  pptrans = transform(pp);
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back(pptrans);
	
  // Top point
  pp[Z] = radius;
  pptrans = transform(pp);
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back(pptrans);
	
  // Gravity center
  pp[Z] = 0.;
  pptrans = transform(pp);  
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back(pptrans);
  
  return ParaviewPoints; 
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de polygones elementaires pour post-processing avec Paraview
int Sphere::numberOfCells_PARAVIEW() const
{
  return (8*visuNodeNbPerQar*visuNodeNbPerQar); 
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe pour post-processing avec Paraview 
void Sphere::write_polygonsStr_PARAVIEW( list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  int i,k,ptsPerlevel = 4*visuNodeNbPerQar,
	Bottom_number = ptsPerlevel*(2*visuNodeNbPerQar-1),
	Top_number = ptsPerlevel*(2*visuNodeNbPerQar-1)+1,  
  	GC_number = ptsPerlevel*(2*visuNodeNbPerQar-1)+2;
  
  // Regular cells: Pyramid
  for ( k = 0; k < 2*visuNodeNbPerQar-2 ; ++k ) 
  {  
    for ( i = 0; i < ptsPerlevel-1 ; ++i )
    {
      connectivity.push_back(firstpoint_globalnumber+k*ptsPerlevel+i); 
      connectivity.push_back(firstpoint_globalnumber+k*ptsPerlevel+i+1);
      connectivity.push_back(firstpoint_globalnumber+(k+1)*ptsPerlevel+i+1);	
      connectivity.push_back(firstpoint_globalnumber+(k+1)*ptsPerlevel+i);
      connectivity.push_back(firstpoint_globalnumber+GC_number);
      last_offset+=5;
      offsets.push_back(last_offset);
      cellstype.push_back(14);		
    }
    connectivity.push_back(firstpoint_globalnumber+k*ptsPerlevel+ptsPerlevel-1);
    connectivity.push_back(firstpoint_globalnumber+k*ptsPerlevel);
    connectivity.push_back(firstpoint_globalnumber+(k+1)*ptsPerlevel);	
    connectivity.push_back(firstpoint_globalnumber
    	+(k+1)*ptsPerlevel+ptsPerlevel-1);
    connectivity.push_back(firstpoint_globalnumber+GC_number);
    last_offset+=5;
    offsets.push_back(last_offset);
    cellstype.push_back(14);    
  }  

  // Bottom cells: tetraedron
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back(firstpoint_globalnumber+i); 
    connectivity.push_back(firstpoint_globalnumber+i+1);
    connectivity.push_back(firstpoint_globalnumber+Bottom_number);	
    connectivity.push_back(firstpoint_globalnumber+GC_number);
    last_offset+=4;
    offsets.push_back(last_offset);
    cellstype.push_back(10);   
  }
  connectivity.push_back(firstpoint_globalnumber+ptsPerlevel-1);
  connectivity.push_back(firstpoint_globalnumber);
  connectivity.push_back(firstpoint_globalnumber+Bottom_number);	
  connectivity.push_back(firstpoint_globalnumber+GC_number);
  last_offset+=4;
  offsets.push_back(last_offset);
  cellstype.push_back(10);  
  
  // Top cells: tetraedron  
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back(firstpoint_globalnumber
    	+(2*visuNodeNbPerQar-2)*ptsPerlevel+i);
    connectivity.push_back(firstpoint_globalnumber
    	+(2*visuNodeNbPerQar-2)*ptsPerlevel+i+1);
    connectivity.push_back(firstpoint_globalnumber+Top_number);	
    connectivity.push_back(firstpoint_globalnumber+GC_number);
    last_offset+=4;
    offsets.push_back(last_offset);
    cellstype.push_back(10); 
  }
  connectivity.push_back(firstpoint_globalnumber
  	+(2*visuNodeNbPerQar-1)*ptsPerlevel-1);
  connectivity.push_back(firstpoint_globalnumber
  	+(2*visuNodeNbPerQar-2)*ptsPerlevel);
  connectivity.push_back(firstpoint_globalnumber+Top_number);	
  connectivity.push_back(firstpoint_globalnumber+GC_number);
  last_offset+=4;
  offsets.push_back(last_offset);
  cellstype.push_back(10); 

  firstpoint_globalnumber += 4*visuNodeNbPerQar*(2*visuNodeNbPerQar-1)+3;
}




// ----------------------------------------------------------------------------
// Renvoie un vecteur orientation de la particule
Vecteur Sphere::vecteurOrientation( Transform const* transform ) const
{
  Point pp(0.,radius,0.);
  Point pptrans = (*transform)(pp); 

  return ( pptrans - *transform->getOrigin() );
} 
