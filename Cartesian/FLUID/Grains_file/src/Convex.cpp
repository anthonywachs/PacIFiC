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
#include "Convex.H"

#include "Basic.H"
#include "BBox.H"
#include "Sphere.H"
#include "Transform.H"

Scalar rel_error = EPSILON;   // relative error in the computed distance
Scalar abs_error = EPSILON2;  // absolute error if the distance is almost zero
int num_iterations = 0;

static Point p[4];         // support points of object A in local coordinates
static Point q[4];         // support points of object B in local coordinates
static Vecteur y[4];       // support points of A - B in world coordinates

static int bits;           // identifies current simplex
static int last;           // identifies last found support point
static int last_bit;       // last_bit = 1<<last
static int all_bits;       // all_bits = bits|last_bit

static Scalar det[16][4];  // cached sub-determinants


// ----------------------------------------------------------------------------
// Constructeur par defaut
// G.FERRER - Octo.2000 - Creation
Convex::Convex() : Shape()
{}




// ----------------------------------------------------------------------------
Convex::~Convex()
{}




// ----------------------------------------------------------------------------
BBox Convex::bbox(const Transform& t) const
{
  Point const* ori = t.getOrigin() ;

  Point min( (*ori)[X] + t.getBasis()[X] * support(-t.getBasis()[X]),
	    (*ori)[Y] + t.getBasis()[Y] * support(-t.getBasis()[Y]),
	    (*ori)[Z] + t.getBasis()[Z] * support(-t.getBasis()[Z]) );
  Point max( (*ori)[X] + t.getBasis()[X] * support(t.getBasis()[X]),
	    (*ori)[Y] + t.getBasis()[Y] * support(t.getBasis()[Y]),
	    (*ori)[Z] + t.getBasis()[Z] * support(t.getBasis()[Z]) );

  return BBox(min, max);
}




// ----------------------------------------------------------------------------
// Determination du rayon du convexe via bricolage
// G.FERRER - Avri.2001 - Creation
Scalar Convex::BuildRayon(const Transform &t) const
{
  Scalar rayon = BuildRayonRef();
  // L'utilisation de la transformation n'est n�cessaire que si on d�forme la
  // particule (allongement, compression ou grossissement) par rapport � la
  // forme de base d�finie dans le fichier de mise en donn�e.
  // Ex: on d�finit un cube puis on l'�tire pour en faire un parall�l�pip�de
  // rectangle ou on le grossit pour en faire un cube plus gros.
  // En pratique, cette option n'est QUASIMENT JAMAIS utilis�e !!
  // Donc la partie ci-dessous est a priori inutile dans le cas general ou
  // la matrice de la tranformation est une ROTATION, la transformation
  // complete etant de type (ROTATION o TRANSLATION), car
  // le rayon circonscrit est INDEPENDANT d'une transformation de type
  // (ROTATION o TRANSLATION).
  // N�anmoins, cette operation ne change pas la valeur du rayon circonscrit
  // dans le cas general ou la matrice de la tranformation est une ROTATION.
  if (!t.isIdentity()) {
    Sphere sphere(rayon);
    BBox   box = sphere.bbox(t);
    rayon = box.size();
  }
  return (rayon);
}



// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// Renvoie le choix de retrecissement
int Convex::getShrinkingMode()
{
  int Shrinking = getShrinkingChoice() ;

  return (Shrinking);
}




// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// set the shrinking radius
void Convex::set_shrinking_radius(Scalar CurrentRadius)
{
  setShrinkingRadius(CurrentRadius);
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Points decrivant l'enveloppe du Convex
// G.FERRER - Aout.2004 - Creation
vector<Point> Convex::getEnveloppe() const
{
  cout << "Warning for this Convex the method Convex::getEnveloppe() "
       << "is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);

  vector<Point> enveloppe;
  return enveloppe;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Description des faces
vector<vector<int> > const* Convex::getFaces() const
{
  cout << "Warning for this Convex the method Convex::getFaces() "
       << "is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);

  vector<vector<int> >* allFaces = NULL;
  return allFaces;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de sommets ou un code �quivalent
int Convex::getNbCorners() const
{
  cout << "Warning for this Convex the method Convex::getNbCorners() "
       << "is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);

  return 0;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de l'objet dans sa configuration courante pour
// post-processing avec GMV
void Convex::GMVoutput(ostream &fileOut,const Transform &transform) const
{
  cout << "Warning for this Convex the method Convex::GMVoutput() "
       << "is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de points pour post-processing avec Paraview
int Convex::numberOfPoints_PARAVIEW() const
{
  cout << "Warning for this Convex the method Convex::numberOfPoints_PARAVIEW()"
       << " is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
  return 0;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de polygones elementaires pour post-processing avec Paraview
int Convex::numberOfCells_PARAVIEW() const
{
  cout << "Warning for this Convex the method Convex::numberOfCells_PARAVIEW()"
       << " is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
  return 0;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit les points du convexe pour post-processing avec Paraview
// D. RAKOTONIRINA - Janv. 2014 - Creation
void Convex::write_polygonsPts_PARAVIEW( ostream &f,
	const Transform &transform, Vecteur const* translation ) const
{
  cout << "Warning for this Convex the method "
       << "Convex::write_polygonsPts_PARAVIEW() is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe au format STL pour lien avec openFoam
// D. RAKOTONIRINA - Janv. 2014 - Creation
void Convex::write_convex_STL( ostream &f, const Transform &transform ) const
{
  cout << "Warning for this Convex the method "
       << "Convex::write_convex_STL() is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les points du convexe pour post-processing avec Paraview
// D. RAKOTONIRINA - Janv. 2014 - Creation
list<Point> Convex::get_polygonsPts_PARAVIEW( const Transform &transform,
	Vecteur const* translation )
	const
{
  cout << "Warning for this Convex the method "
       << "Convex::get_polygonsPts_PARAVIEW() is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe pour post-processing avec Paraview
// D. RAKOTONIRINA - Janv. 2014 - Creation
void Convex::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  cout << "Warning for this Convex the method "
       << "Convex::write_polygonsStr_PARAVIEW() is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie un vecteur orientation du convex
Vecteur Convex::vecteurOrientation(Transform const* transform) const
{
  Vecteur vecNul(0.);
  return vecNul;
}




// ----------------------------------------------------------------------------
void compute_det()
{
  static Scalar dp[4][4];

  for (int i = 0, bit = 1; i < 4; ++i, bit <<=1)
    if (bits & bit) dp[i][last] = dp[last][i] = y[i] * y[last];
  dp[last][last] = y[last] * y[last];

  det[last_bit][last] = 1;
  for (int j = 0, sj = 1; j < 4; ++j, sj <<= 1) {
    if (bits & sj) {
      int s2 = sj|last_bit;
      det[s2][j] = dp[last][last] - dp[last][j];
      det[s2][last] = dp[j][j] - dp[j][last];
      for (int k = 0, sk = 1; k < j; ++k, sk <<= 1) {
	if (bits & sk) {
	  int s3 = sk|s2;
	  det[s3][k] = det[s2][j] * (dp[j][j] - dp[j][k]) +
	               det[s2][last] * (dp[last][j] - dp[last][k]);
	  det[s3][j] = det[sk|last_bit][k] * (dp[k][k] - dp[k][j]) +
	               det[sk|last_bit][last] * (dp[last][k] - dp[last][j]);
	  det[s3][last] = det[sk|sj][k] * (dp[k][k] - dp[k][last]) +
	                  det[sk|sj][j] * (dp[j][k] - dp[j][last]);
	}
      }
    }
  }
  if (all_bits == 15) {
    det[15][0] = det[14][1] * (dp[1][1] - dp[1][0]) +
                 det[14][2] * (dp[2][1] - dp[2][0]) +
                 det[14][3] * (dp[3][1] - dp[3][0]);
    det[15][1] = det[13][0] * (dp[0][0] - dp[0][1]) +
                 det[13][2] * (dp[2][0] - dp[2][1]) +
                 det[13][3] * (dp[3][0] - dp[3][1]);
    det[15][2] = det[11][0] * (dp[0][0] - dp[0][2]) +
                 det[11][1] * (dp[1][0] - dp[1][2]) +
                 det[11][3] * (dp[3][0] - dp[3][2]);
    det[15][3] = det[7][0] * (dp[0][0] - dp[0][3]) +
                 det[7][1] * (dp[1][0] - dp[1][3]) +
                 det[7][2] * (dp[2][0] - dp[2][3]);
  }
}




// ----------------------------------------------------------------------------
inline bool valid(int s)
{
  // checks if the simplex fulfills the 2 conditions for being the smallest one to include v
  for (int i = 0, bit = 1; i < 4; ++i, bit <<= 1) {
    if (all_bits & bit) {
      if (s & bit) {
	if (det[s][i] <= EPSILON3)
	  return false;
      } else if (det[s|bit][i] > 0.)
	return false;
    }
  }
  return true;
}




// ----------------------------------------------------------------------------
inline void compute_vector(int bits_, Vecteur& v)
{
  Scalar sum = 0.;
  v.setValue(0., 0., 0.);
  for (int i = 0, bit = 1; i < 4; ++i, bit <<= 1) {
    if (bits_ & bit) {
      sum += det[bits_][i];
      v += y[i] * det[bits_][i];
    }
  }
  v *= 1 / sum;
}




// ----------------------------------------------------------------------------
inline void compute_points(int bits_, Point& p1, Point& p2)
{
  Scalar sum = 0;
  p1.setValue(0, 0, 0);
  p2.setValue(0, 0, 0);
  for (int i = 0, bit = 1; i < 4; ++i, bit <<= 1) {
    if (bits_ & bit) {
      sum += det[bits_][i];
      p1 += p[i] * det[bits_][i];
      p2 += q[i] * det[bits_][i];
    }
  }
  Scalar s = 1 / sum;
  p1 *= s;
  p2 *= s;
}




// ----------------------------------------------------------------------------
inline bool proper(int s)
{
  for (int i = 0, bit = 1; i < 4; ++i, bit <<= 1)
    if ((s & bit) && det[s][i] <= EPSILON3) return false;
  return true;
}




// ----------------------------------------------------------------------------
inline bool closest(Vecteur& v)
{
  int s;
  compute_det();    // compute once and for all the 16 determinants needed in the following
  for (s = bits; s; --s) {  // go through the different sub-simplices (from many to low vertices)
    if ((s & bits) == s) {  // (make sure the sub-simplex considered is included in the main simplex)
      if (valid(s|last_bit)) {    // does the considered simplex + the last point found fulfills the
                                                // criteria for being the smallest sub-simplex to include last point?
	bits = s|last_bit;             // if yes, validate this simplex by adding the last point to it
 	compute_vector(bits, v);  // and compute this shortest vector from the simplex to the origin
	return true;
      }
    }
  }
  if (valid(last_bit)) {    // last simplex we need to consider is the last found point alone
    bits = last_bit;
    v = y[last];
    return true;
  }
  // Original GJK calls the backup procedure at this point.
  Scalar min_dist2 = INFINITY;
  for (s = all_bits; s; --s) {
    if ((s & all_bits) == s) {
      if (proper(s)) {
	Vecteur u;
 	compute_vector(s, u);
	Scalar dist2 = Norm2(u);
	if (dist2 < min_dist2) {
	  min_dist2 = dist2;
	  bits = s;
	  v = u;
	}
      }
    }
  }

  return false;
}




// ----------------------------------------------------------------------------
// The next function is used for detecting degenerate cases that cause
// termination problems due to rounding errors.
inline bool degenerate(const Vecteur& w)
{
  for (int i = 0, bit = 1; i < 4; ++i, bit <<= 1)
    if ((all_bits & bit) && y[i] == w)  return true;
  return false;
}




// ----------------------------------------------------------------------------
bool intersect(const Convex& a, const Convex& b,
	       const Transform& a2w, const Transform& b2w,
	       Vecteur& v)
{
  Vecteur w;

  bits = 0;
  all_bits = 0;

  Scalar prod;

  do {
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }
    w = a2w(a.support((-v) * a2w.getBasis())) -
      b2w(b.support(v * b2w.getBasis()));
    prod = v * w;
    if (prod > 0. || fabs(prod) < EPSILON2) {
      return false;
    }
    if (degenerate(w)) {
      return false;
    }
    y[last] = w;
    all_bits = bits|last_bit;
    if (!closest(v)) {
      return false;
    }
  }
  while (bits < 15 && !approxZero(v));
  return true;
}




// ----------------------------------------------------------------------------
bool intersect(const Convex& a, const Convex& b, const Transform& b2a,
	       Vecteur& v)
{
  Vecteur w;

  bits = 0;
  all_bits = 0;
  Scalar prod;

  do {
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }
    w = a.support(-v) - b2a(b.support(v * b2a.getBasis()));
    prod = v * w;
    if (prod > 0. || fabs(prod) < EPSILON2) {
      return false;
    }
    if (degenerate(w)) {
      return false;
    }
    y[last] = w;
    all_bits = bits|last_bit;
    if (!closest(v)) {
      return false;
    }
  }
  while (bits < 15 && !approxZero(v));
  return true;
}




// ----------------------------------------------------------------------------
bool common_point(const Convex& a, const Convex& b,
		  const Transform& a2w, const Transform& b2w,
		  Vecteur& v, Point& pa, Point& pb)
{
  Vecteur w;

  bits = 0;
  all_bits = 0;
  Scalar prod;

  do {
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }
    p[last] = a.support((-v) * a2w.getBasis());
    q[last] = b.support(v * b2w.getBasis());
    w = a2w(p[last]) - b2w(q[last]);
    prod = v * w;
    if (prod > 0. || fabs(prod) < EPSILON2) {
      return false;
    }
    if (degenerate(w)) {
      return false;
    }
    y[last] = w;
    all_bits = bits|last_bit;
    if (!closest(v)) {
      return false;
    }
  }
  while (bits < 15 && !approxZero(v) ) ;
  compute_points(bits, pa, pb);
  return true;
}




// ----------------------------------------------------------------------------
bool common_point(const Convex& a, const Convex& b, const Transform& b2a,
		  Vecteur& v, Point& pa, Point& pb)
{
  Vecteur w;

  bits = 0;
  all_bits = 0;
  Scalar prod;

  do {
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }
    p[last] = a.support(-v);
    q[last] = b.support(v * b2a.getBasis());
    w = p[last] - b2a(q[last]);
    prod = v * w;
    if (prod > 0. || fabs(prod) < EPSILON2) {
      return false;
    }
    if (degenerate(w)) {
      return false;
    }
    y[last] = w;
    all_bits = bits|last_bit;
    if (!closest(v)) {
      return false;
    }
  }
  while (bits < 15 && !approxZero(v) );
  compute_points(bits, pa, pb);
  return true;
}




// ----------------------------------------------------------------------------
// For num_iterations > 1000
void catch_me()
{
  cerr << "closest_points : Out on iteration > 1000\n";
}




// ----------------------------------------------------------------------------
// G.FERRER - Juil.2003 - Ajout de num_iterations
Scalar closest_points(const Convex& a, const Convex& b,
		      const Transform& a2w, const Transform& b2w,
		      Point& pa, Point& pb,int& nbIter)
{
//  static Vecteur zero(-EPSILON, EPSILON, EPSILON);
  static Vecteur zero(0., 0., 0.);

  Vecteur v = a2w(a.support(zero)) - b2w(b.support(zero));
  Scalar dist = Norm(v);

  Vecteur w;

  bits = 0;
  all_bits = 0;
  Scalar mu = 0;

  num_iterations = 0;

  while (bits < 15 && dist > abs_error && num_iterations < 1000) {
    // -> bits<15 because a simplex should contain 4 vertices maximum at all times (so no more
    // than 4 bits needed)
    // -> dist > abs_error is the standard break condition: when no significant change is made to
    // the vector we consider this is a close enough solution
    // -> num_iterations is just ensuring we don't spend too much time on a difficult case
    last = 0;
    last_bit = 1;
    // the line below finds where the new support point should be stored
    while (bits & last_bit) { ++last; last_bit <<= 1; }
    p[last] = a.support((-v) * a2w.getBasis()); // necessary for compute_points routine

    q[last] = b.support(v * b2w.getBasis());    // necessary for compute_points routine
    w = a2w(p[last]) - b2w(q[last]);    // compute the new vertex of the simplex

    set_max(mu, v*w / dist);
    if (dist - mu <= dist * rel_error) break;   // If the algorithm is close enough the the solution
    if (degenerate(w)) {
      break;
    }
    y[last] = w;    // adds the newly found vertex to the simplex y and (below) to the bit array
    all_bits = bits|last_bit;

    ++num_iterations;

    if (!closest(v)) {  // compute the smallest simplex to contain v (the smallest Norm2)
      break;    // if this fails, backup procedure has occured and exit main loop to compute points
    }
    dist = Norm(v);
  }
  compute_points(bits, pa, pb); // compute the closest points of convexes A and B
  if (num_iterations > 1000) catch_me();    // if no solution has been found, return an error
  else nbIter=num_iterations;

  return dist;
}




// ---------------------------------------------------------------------
// Operateur d'ecriture
// G.FERRER - Juin.2000 - Creation
ostream &operator << (ostream &fileOut,
                      const Convex &convex)
{
  convex.printClass(fileOut);
  return (fileOut);
}




// ----------------------------------------------------------------------------
// Operateur de lecture
// G.FERRER - Juin.2000 - Creation
istream &operator >> (istream &fileIn,
                      Convex &convex)
{
  convex.readClass(fileIn);
  return (fileIn);
}




// ----------------------------------------------------------------------------
// Determine si il y a intersection entre deux convexes
// D. RAKOTONIRINA - Mars 2014 - Creation
bool Convex::isIn(const Convex& a, const Convex& b, const Transform& a2w,
	const Transform& b2w, Vecteur& v)
{
  return( intersect(a, b, a2w, b2w, v) );
}




// ----------------------------------------------------------------------------
// Determine si il y a intersection entre deux convexes
// D. RAKOTONIRINA - Mars 2014 - Creation
bool Convex::isIn(const Convex& a, const Convex& b, const Transform& b2a,
	Vecteur& v)
{
  return( intersect(a, b, b2a, v) );
}
