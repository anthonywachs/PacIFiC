/* 
   D.PETIT - Creation
   
   Complement de l'algorithme GJK - Engine
   Nouvelle classe de convexe : Le segment.
*/

#include "Segment.H"
#include "Transform.H"
#include "Basic.H"


// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction avec un initialisation
Segment::Segment(Scalar x) :
  Convex(), halflength(x / 2.)
{
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction avec un fichier
Segment::Segment(istream &fileIn) :
  Convex()
{
  readClass(fileIn);
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine l'inertie et l'inertie inverse d'un Segment
bool Segment::BuildInertie(Scalar *inertie,Scalar *inertie_1) const
{
  inertie[0]=inertie[1]=inertie[2]=inertie[3]=inertie[4]
    =inertie[5]=0.0;
  inertie_1[0]=inertie_1[1]=inertie_1[2]=inertie_1[3]=inertie_1[4]
    =inertie_1[5]=0.0;
  return true;
}




// ----------------------------------------------------------------------------
// M. SULAIMAN - Nov.2015 - Creation
// get the shrinking choice
 int Segment::getShrinkingChoice()const
{
  return(0);
}




// ----------------------------------------------------------------------------
// M. SULAIMAN - Nov.2015 - Creation
// Set the shrinking radius value
void Segment::setShrinkingRadius(Scalar CurrentRadius)
{
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
//  Determine le rayon circonscrit d'un segment
Scalar Segment::BuildRayonRef() const 
{
  return halflength;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction d'un clone d'un Segment
Convex* Segment::clone() const 
{
  return new Segment(2.0*halflength);
}




// ----------------------------------------------------------------------------
Scalar Segment::getLength() const 
{
  return 2.0*halflength;
}




// ----------------------------------------------------------------------------
Point Segment::support(const Vecteur& v) const 
{
  Scalar norm = Norm(v);
  if (norm > EPSILON) {
    return (v[X] > 0.0 ? 
	    Point(halflength,0.0,0.0) : Point(-halflength,0.0,0.0));	    
  } else {
    return Point();
  }
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le Volume d'un Segment
Scalar Segment::getVolume() const 
{
  return 0;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture d'un Segment
void Segment::printClass(ostream &fileOut) const 
{
  fileOut << "*Segment\n";
  fileOut << 2.0*halflength;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// lecture d'un Segment
void Segment::readClass(istream &fileIn) 
{
  fileIn >> halflength;
  halflength /= 2.0;
}




// ----------------------------------------------------------------------------
// A.WACHS - Juil. 2010 - Creation
// Renvoi la tranformation pour un segment associé au vecteur v et au point gc
Transform Segment::computeTransform( const Vecteur& v, const Point &gc )
{
  Transform tr;

  double normxy = sqrt( v[X]*v[X] + v[Y]*v[Y] ),
  	normxz = sqrt( v[X]*v[X] + v[Z]*v[Z] );
  double bx = fabs(v[X]);

  // Angle pr la rotation par rapport à Z
  double angleZ = 0.;
  if ( normxy > 1.e-12 )
  { 
    angleZ = acos( bx / normxy );
    if ( v[X] > 0. )
    {
     if ( v[Y] < 0. ) angleZ *= -1.;
    }
    else
    {
      if ( v[Y] > 0. ) angleZ = PI - angleZ;
      else angleZ += PI;
    }
  }
   
  // Angle pr la rotation par rapport à Y 
  double angleY = 0.;
  if ( normxz > 1.e-12 )
  {
    angleY = acos( bx / normxz );
    if ( v[X] > 0. )
    {
     if ( v[Z] < 0. ) angleY *= -1.;
    }
    else
    {
      if ( v[Z] > 0. ) angleY = PI - angleY;
      else angleY += PI; 
    }
  } 

//   cout << "AngleZ = " << angleZ*180/PI 
//   	<< "   AngleY = " << angleY*180/PI << endl;

  // Matrice de rotation
  Matrix rZ(cos(angleZ),-sin(angleZ),0.,sin(angleZ),cos(angleZ),0.,0.,0.,1.);
  Matrix rY(cos(angleY),0.,sin(angleY),0.,1.,0.,-sin(angleY),0.,cos(angleY));
  Matrix rotation = rY * rZ;
  
  // Set transform
  tr.setBasis(rotation);
  tr.setOrigin(gc);
   
  return tr;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de points pour post-processing avec Paraview
int Segment::numberOfPoints_PARAVIEW() const
{
  return 2;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de polygones elementaires pour post-processing avec Paraview
int Segment::numberOfCells_PARAVIEW() const
{
  return 1;  
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit les points du convexe pour post-processing avec Paraview  
void Segment::write_polygonsPts_PARAVIEW( ostream &f,
	const Transform &transform, Vecteur const* translation ) const
{
  Point pp, p;
  
  // Point -
  p[X] = - halflength;
  pp = transform(p);
  if ( translation ) pp += *translation;    
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  
  // Point +
  p[X] = halflength;
  pp = transform(p);
  if ( translation ) pp += *translation;    
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les points du convexe pour post-processing avec Paraview  
list<Point> Segment::get_polygonsPts_PARAVIEW( const Transform &transform,
	Vecteur const* translation ) const
{
  list<Point> ParaviewPoints;
  Point pp, p;
  
  // Point -
  p[X] = - halflength;
  pp = transform(p);
  if ( translation ) pp += *translation;    
  ParaviewPoints.push_back(pp);
  
  // Point +
  p[X] = halflength;
  pp = transform(p);
  if ( translation ) pp += *translation;    
  ParaviewPoints.push_back(pp);
  
  return ParaviewPoints; 
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe pour post-processing avec Paraview 
void Segment::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  connectivity.push_back(firstpoint_globalnumber);
  connectivity.push_back(firstpoint_globalnumber+1); 
  last_offset+=2;  
  offsets.push_back(last_offset);
  cellstype.push_back(3);
  
  firstpoint_globalnumber+=2;
}
