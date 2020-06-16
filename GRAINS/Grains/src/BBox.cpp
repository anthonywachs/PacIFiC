#include "BBox.H"


// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
BBox::BBox() 
{
} 




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
BBox::BBox(const Point& min, const Point& max) 
{ 
  setValue(min, max);
}




// --------------------------------------------------------------------
// Constructeur par copie
// A.WACHS - Aout.2009 - Creation
BBox::BBox(const BBox &bbox_)
{
  center = bbox_.center;
  extent = bbox_.extent;
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
BBox::~BBox()
{
}




// --------------------------------------------------------------------
// Operateur d'affectation 
// A.WACHS - Aout.2009 - Creation
BBox& BBox::operator= (const BBox& rhs)
{
  if ( &rhs != this )
  {      
    center = rhs.center;
    extent = rhs.extent;   
  }
  return (*this);
}




// ----------------------------------------------------------------------------
// BBox intersection des deux BBox
// G.FERRER - Juin.2003 - Creation
void BBox::closest(const BBox &a, const BBox &b)
{
  if ( intersect(a,b) ) {
    Point lower(max(a.getLower(X), b.getLower(X)),
		max(a.getLower(Y), b.getLower(Y)),
		max(a.getLower(Z), b.getLower(Z)));
    Point upper(min(a.getUpper(X), b.getUpper(X)),
		min(a.getUpper(Y), b.getUpper(Y)),
		min(a.getUpper(Z), b.getUpper(Z)));
    setValue(lower, upper);
  }
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
void BBox::enclose(const BBox& a, const BBox& b) 
{
  Point lower(min(a.getLower(X), b.getLower(X)),
	      min(a.getLower(Y), b.getLower(Y)),
	      min(a.getLower(Z), b.getLower(Z)));
  Point upper(max(a.getUpper(X), b.getUpper(X)),
	      max(a.getUpper(Y), b.getUpper(Y)),
	      max(a.getUpper(Z), b.getUpper(Z)));
  setValue(lower, upper);
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
const Point& BBox::getCenter() const 
{
  return center; 
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
const Vecteur& BBox::getExtent() const 
{
  return extent; 
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
Scalar BBox::getLower(int i) const 
{ 
  return center[i] - extent[i]; 
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
Scalar BBox::getUpper(int i) const 
{
  return center[i] + extent[i]; 
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
void BBox::include(const Point& p) 
{
  Point lower(min(getLower(X), p[X]),
	      min(getLower(Y), p[Y]),
	      min(getLower(Z), p[Z]));
  Point upper(max(getUpper(X), p[X]),
	      max(getUpper(Y), p[Y]),
	      max(getUpper(Z), p[Z]));
  setValue(lower, upper);
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
void BBox::include(const BBox& b) 
{
  enclose(*this, b); 
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
bool BBox::InZone(Point const* p, const Scalar &arete) const 
{
  return (fabs(center[X] - (*p)[X]) <= extent[X] + arete &&
    fabs(center[Y] - (*p)[Y]) <= extent[Y] + arete &&
    fabs(center[Z] - (*p)[Z]) <= extent[Z] + arete);
}




// --------------------------------------------------------------------
// A.WACHS - Aout.2009 - Creation
bool BBox::InZone(Point const* p, const Scalar &arete_X, const Scalar &arete_Y,
  	const Scalar &arete_Z) const 
{
  return (fabs(center[X] - (*p)[X]) <= extent[X] + arete_X &&
    fabs(center[Y] - (*p)[Y]) <= extent[Y] + arete_Y &&
    fabs(center[Z] - (*p)[Z]) <= extent[Z] + arete_Z);
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
int BBox::longestAxis() const 
{
  return extent.closestAxis(); 
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
void BBox::setCenter(const Point& p)  
{
  center = p; 
}




// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
void BBox::setEmpty() 
{ 
  center.setValue(0, 0, 0); 
  extent.setValue(-INFINITY, -INFINITY, -INFINITY);
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
void BBox::setExtent(const Vecteur& v) 
{
  extent = v; 
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation  
void BBox::setValue(const Point& min, const Point& max) 
{ 
  extent = (max - min) / 2.;
  center = min + extent; 
}




// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
Scalar BBox::size() const 
{
  return max(max(extent[X], extent[Y]), extent[Z]); 
}




// ----------------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
bool intersect(const BBox& a, const BBox& b) 
{
  return fabs(a.center[X] - b.center[X]) <= a.extent[X] + b.extent[X] &&
    fabs(a.center[Y] - b.center[Y]) <= a.extent[Y] + b.extent[Y] &&
    fabs(a.center[Z] - b.center[Z]) <= a.extent[Z] + b.extent[Z];
}




#include <iostream>
using namespace std;
// ----------------------------------------------------------------------------
// Debug
void BBox::debug(char *s) const
{
  cerr << s;
  cerr << "BBox Center " << center
       << "     Extent " << extent;
}




// ----------------------------------------------------------------------------
// Operateur << 
// A.WACHS - Aout.2009 - Creation
ostream& operator << (ostream &f, const BBox &B)
{
  f << "BBox: Center = " << B.center;
  f << "      Extent = " << B.extent;
  return f;  
}
