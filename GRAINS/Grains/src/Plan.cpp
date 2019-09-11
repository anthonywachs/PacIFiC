// G.FERRER - Dece.1999 - Creation
// ============================================================================
#include "Plan.H"

#include <math.h>
#include <stdlib.h>


// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
// Constructeur par defaut
Plan::Plan() :
  nbrePts(0), a(0.0), b(0.0), c(0.0), d(0.0)
{
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin 2000 - Creation
// Constructeur avec initialisation
Plan::Plan(const Point& p1, const Point& p2, const Point& p3)
{
  // Calcul de la normale
  a =  (p2[Y]-p1[Y])*(p3[Z]-p1[Z]) - (p2[Z]-p1[Z])*(p3[Y]-p1[Y]);
  b = -(p2[X]-p1[X])*(p3[Z]-p1[Z]) + (p2[Z]-p1[Z])*(p3[X]-p1[X]);
  c =  (p2[X]-p1[X])*(p3[Y]-p1[Y]) - (p2[Y]-p1[Y])*(p3[X]-p1[X]);
  d = sqrt(a*a + b*b + c*c);
  a = a / d;
  b = b / d;
  c = c / d;
  d = a*p1[X] + b*p1[Y] + c*p1[Z];
}




// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
// Constructeur avec initialisation
Plan::Plan(int nbre, const Scalar *xPt, const Scalar *yPt, const Scalar *zPt)
{
  if (nbre < 3) {
    cerr << "Erreur dans construction du plan, 3 points minimum !\n";
    exit(1);
  }

  nbrePts = nbre;

  // Coordonnees des points
  x = new Scalar[nbrePts];
  y = new Scalar[nbrePts];
  z = new Scalar[nbrePts];
  for (int i=0; i<nbrePts; i++) {
    x[i] = xPt[i];
    y[i] = yPt[i];
    z[i] = zPt[i];
  }

  // Calcul de la normale
  a =  (y[Y]-y[X])*(z[Z]-z[X]) - (z[Y]-z[X])*(y[Z]-y[X]);
  b = -(x[Y]-x[X])*(z[Z]-z[X]) + (z[Y]-z[X])*(x[Z]-x[X]);
  c =  (x[Y]-x[X])*(y[Z]-y[X]) - (y[Y]-y[X])*(x[Z]-x[X]);
  d = sqrt(a*a + b*b + c*c);
  a = a / d;
  b = b / d;
  c = c / d;
  d = a*x[X] + b*y[X] + c*z[X];
}




// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
// Constructeur avec initialisation
Plan::Plan(const Scalar *normale)
{
  nbrePts = 0;
  a = normale[X];
  b = normale[Y];
  c = normale[Z];
  d = normale[3];
}




// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
// Destructeur
Plan::~Plan()
{
  if (nbrePts > 0) {
    delete x;
    delete y;
    delete z;
  }
}




// -----------------------------------------------------------
void Plan::setPlan(const Point& p1, const Point& p2, const Point& p3){
  // Calcul de la normale
  a =  (p2[Y]-p1[Y])*(p3[Z]-p1[Z]) - (p2[Z]-p1[Z])*(p3[Y]-p1[Y]);
  b = -(p2[X]-p1[X])*(p3[Z]-p1[Z]) + (p2[Z]-p1[Z])*(p3[X]-p1[X]);
  c =  (p2[X]-p1[X])*(p3[Y]-p1[Y]) - (p2[Y]-p1[Y])*(p3[X]-p1[X]);
  d = sqrt(a*a + b*b + c*c);
  a = a / d;
  b = b / d;
  c = c / d;
  d = a*p1[X] + b*p1[Y] + c*p1[Z];
}




// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
// Deplacement de dist normale au plan
void Plan::Deplacer(Scalar dist)
{
  if (nbrePts > 0) {
    Scalar n = sqrt(a*a + b*b + c*c); 
    for (int i=0; i<nbrePts; i++) {
      x[i] += dist*a / n;
      y[i] += dist*b / n;
      z[i] += dist*c / n;
    }
    d = a*x[X] + b*y[X] + c*z[X]; 
  }
}




// ----------------------------------------------------------------------------
void Plan::Deplacer(const Scalar *dist)
{
  if (nbrePts > 0) {
    for (int i=0; i<nbrePts; i++) {
      x[i] += dist[X];
      y[i] += dist[Y];
      z[i] += dist[Z];
    }
    d = a*x[X] + b*y[X] + c*z[X]; 
  }
}




// ----------------------------------------------------------------------------
void Plan::Deplacer(const Vecteur &dist)
{
  if (nbrePts > 0) {
    for (int i=0; i<nbrePts; i++) {
      x[i] += dist[X];
      y[i] += dist[Y];
      z[i] += dist[Z];
    }
    d = a*x[X] + b*y[X] + c*z[X]; 
  }
}




// ----------------------------------------------------------------------------
void Plan::Deplacer(Scalar distX, Scalar distY, Scalar distZ)
{
  if (nbrePts > 0) {
    for (int i=0; i<nbrePts; i++) {
      x[i] += distX;
      y[i] += distY;
      z[i] += distZ;
    }
    d = a*x[X] + b*y[X] + c*z[X]; 
  }
}




// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
Scalar Plan::DistanteDe(const Point &point) const
{
  Scalar dist;
  dist = (a*point[X] + b*point[Y] + c*point[Z] - d) / sqrt(a*a + b*b + c*c);
  return (dist);
}




// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
Scalar Plan::DistanteDe(const Scalar *pt) const
{
  Scalar dist;
  dist = (a*pt[X] + b*pt[Y] + c*pt[Z] - d) / sqrt(a*a + b*b + c*c);
  return (dist);
}




// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
Scalar Plan::DistanteDe(Scalar x_, Scalar y_, Scalar z_) const
{
  Scalar dist;
  dist = (a*x_ + b*y_ + c*z_ - d) / sqrt(a*a + b*b + c*c);
  return (dist);
}




// ----------------------------------------------------------------------------
// G.FERRER - Fevr.2000 - Creation
// Coordonnees d'un point en projection sur le plan
Point Plan::Projection(const Point &point) const
{  
  Scalar proj = -(a*point[X]+b*point[Y]+c*point[Z]-d)/sqrt(a*a + b*b + c*c);
  return Point(point[X]+proj*a,point[Y]+proj*b,point[Z]+proj*c);
}




// ----------------------------------------------------------------------------
Point Plan::Projection(const Scalar *point) const
{
  Scalar proj = -(a*point[X]+b*point[Y]+c*point[Z]-d)/sqrt(a*a + b*b + c*c);
  return Point(point[X]+proj*a,point[Y]+proj*b,point[Z]+proj*c);
}




// ----------------------------------------------------------------------------
Point Plan::Projection(const Scalar x_,const Scalar y_,const Scalar z_) const{
  Scalar proj = -(a*x_+b*y_+c*z_-d)/sqrt(a*a + b*b + c*c);
  return Point(x_+proj*a,y_+proj*b,z_+proj*c);
}




// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
// Validation du plan par rapport a un point
void Plan::Valider(Scalar *point)
{
  Scalar dist;
  dist = DistanteDe(point[X], point[Y], point[Z]);

  if (dist < 0.0) {
    Scalar xi, yi, zi;
    xi=x[Y];   yi=y[Y];   zi=z[Y];
    x[Y]=x[Z]; y[Y]=y[Z]; z[Y]=z[Z];
    x[Z]=xi;   y[Z]=yi;   z[Z]=zi;

    // Calcul de la normale
    a =  (y[Y]-y[X])*(z[Z]-z[X]) - (z[Y]-z[X])*(y[Z]-y[X]);
    b = -(x[Y]-x[X])*(z[Z]-z[X]) + (z[Y]-z[X])*(x[Z]-x[X]);
    c =  (x[Y]-x[X])*(y[Z]-y[X]) - (y[Y]-y[X])*(x[Z]-x[X]);
    d = sqrt(a*a + b*b + c*c);
    a = a / d;
    b = b / d;
    c = c / d;
    d = a*x[X] + b*y[X] + c*z[X];
  }
}




//-----------------------------------------------------------------------------
// G.FERRER - Nove.1999 - Creation
// Position du point par rapport au plan
PositionPlan Plan::WhereIs(const Point &point) const
{
  PositionPlan position;
  Scalar dist;
  dist = (a*point[X] + b*point[Y] + c*point[Z] - d) / sqrt(a*a + b*b + c*c);
  if (dist < 0.0)
    position = down;
  else if (dist == 0.0)
    position = in;
  else 
    position = up;
  return position;
}




//-----------------------------------------------------------------------------
PositionPlan Plan::WhereIs(const Scalar *pt) const
{
  PositionPlan position;
  Scalar dist;
  dist = (a*pt[X] + b*pt[Y] + c*pt[Z] - d) / sqrt(a*a + b*b + c*c);
  if (dist< 0.0)
    position = down;
  else if (dist == 0.0)
    position = in;
  else 
    position = up;
  return position;
}




//-----------------------------------------------------------------------------
PositionPlan Plan::WhereIs(Scalar x_, Scalar y_, Scalar z_) const
{
  PositionPlan position;
  Scalar dist;
  dist = (a*x_ + b*y_ + c*z_ - d) / sqrt(a*a + b*b + c*c);
  if (dist< 0.0)
    position = down;
  else if (dist == 0.0)
    position = in;
  else 
    position = up;
  return position;
}
