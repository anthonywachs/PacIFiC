#include "Vecteur.H"

#include "Quaternion.H"

#include <math.h>


namespace solid
{
  // ---------------------------------------------------------------------------
  // Constructeur par defaut
  Vecteur::Vecteur(Scalar def) :
    Group3(def)
  {
  }




  // ---------------------------------------------------------------------------
  // Constructeur avec initialisation
  Vecteur::Vecteur(Scalar x, Scalar y, Scalar z) :
    Group3(x, y, z)
  {
  }




  // ---------------------------------------------------------------------------
  // Constructeur par copie
  Vecteur::Vecteur(const Group3 &g) :
    Group3(g)
  {
  }




  // ---------------------------------------------------------------------------
  // Destructeur
  Vecteur::~Vecteur()
  {
  }




  // ---------------------------------------------------------------------------
  // Determine l'axe la plus proche du Vecteur. 
  int Vecteur::closestAxis() const
  {
    Scalar a[2];
    int axis = (a[X] = fabs(comp[X])) < (a[Y] = fabs(comp[Y])) ? Y : X;
    return a[axis] < fabs(comp[Z]) ? Z : axis;
  }




  // ---------------------------------------------------------------------------
  // Normalise un Vecteur. 
  void Vecteur::normalize()
  {
    *this /= Norm(*this);
  }




  // ---------------------------------------------------------------------------
  // Normalise un Vecteur. 
  Vecteur Vecteur::normalized() const
  {
    return *this/Norm(*this);
  } 




  // ---------------------------------------------------------------------------
  // arrondi des valeurs proches de zeros
  void Vecteur::round() 
  {
    comp[X] = fabs(comp[X])<EPS ? 0. : comp[X];
    comp[Y] = fabs(comp[Y])<EPS ? 0. : comp[Y];
    comp[Z] = fabs(comp[Z])<EPS ? 0. : comp[Z];
  }




  // ---------------------------------------------------------------------------
  // creation de la rotation d'un vecteur 
  void Vecteur::Tourne(const Quaternion& q)
  {
    Quaternion tmp(*this);
    tmp = (q*tmp) * q.Conjugate();  
    // il serait plus approprie d'utiliser tmp=(q*tmp)*q.Conjug(); 
    // car q est unitaire dans le cas d'une rotation .
    *this = *tmp.getVecteur();
  }




  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Vecteur& Vecteur::operator = (const Vecteur &g2)
  {
    if ( &g2 != this )
    {     
      comp[X] = g2.comp[X];
      comp[Y] = g2.comp[Y];
      comp[Z] = g2.comp[Z];
    }
    return (*this);
  }




  // ---------------------------------------------------------------------------
  Vecteur Vecteur::operator^(const Vecteur& rhs) const
  {
    return Vecteur( comp[1]*rhs.comp[2] - comp[2]*rhs.comp[1],
		    -comp[0]*rhs.comp[2] + comp[2]*rhs.comp[0], 
		    comp[0]*rhs.comp[1] - comp[1]*rhs.comp[0]);
  }
  



  // ---------------------------------------------------------------------------
  // test si la longueur du vecteur est nul. 
  bool approxZero(const Vecteur &v)
  {
    return Norm2(v) < EPSILON2;
  }




  // --------------------------------------------------------------------
  // Determine le cos de l'angle formee entre x et y
  Scalar cos(const Vecteur& x, const Vecteur& y)
  {
    return x*y/(Norm(x) * Norm(y));
  }




  // ---------------------------------------------------------------------------
  // Determine la norme d'un vecteur
  Scalar Norm(const Vecteur& v)
  {
    return sqrt(v.comp[X]*v.comp[X]+v.comp[Y]*v.comp[Y]+v.comp[Z]*v.comp[Z]);
  }




  // ---------------------------------------------------------------------------
  // Determine la norme au carre d'un vecteur
  Scalar Norm2(const Vecteur& v)
  {
    return (v.comp[X]*v.comp[X]+v.comp[Y]*v.comp[Y]+v.comp[Z]*v.comp[Z]);
  }  
}



