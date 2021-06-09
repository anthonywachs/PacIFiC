// G.FERRER - Dece.1999 - Creation
// ============================================================================
#include "Point.H"

#include "Vecteur.H"

#include <math.h>


namespace solid
{
  // ---------------------------------------------------------------------------
  // G.FERRER - Dece.1999 - Creation
  // Constructeur par defaut
  Point::Point(Scalar def) :
    Group3(def)
  {
  }




  // ---------------------------------------------------------------------------
  // G.FERRER - Dece.1999 - Creation
  // Constructeur avec initialisation
  Point::Point(Scalar x, Scalar y, Scalar z) :
    Group3(x, y, z)
  {
  }



  // ---------------------------------------------------------------------------
  // M.SULAIMAN - Nov.2015 - Creation
  // return choice
  int Point::getShrinkingChoice()const
  {
    return(0);
  }
  


  
  // ---------------------------------------------------------------------------
  // G.FERRER - Dece.1999 - Creation
  // Constructeur par copie
  Point::Point(const Group3 &point) :
    Group3(point)
  {
  }




  // ---------------------------------------------------------------------------
  // G.FERRER - Dece.1999 - Creation
  // Destructeur
  Point::~Point()
  {
  }




  // ---------------------------------------------------------------------------
  // G.FERRER - Dece.1999 - Creation
  // Deplacement de l'objet
  void Point::Deplacer(Scalar dist)
  {
    for (int i=0; i<3; i++) comp[i] += dist;
  }




  // ---------------------------------------------------------------------------
  void Point::Deplacer(const Scalar *dist)
  {
    for (int i=0; i<3; i++) comp[i] += dist[i];
  }




  // ---------------------------------------------------------------------------
  void Point::Deplacer(Scalar distX, Scalar distY, Scalar distZ)
  {
    comp[0] += distX;
    comp[1] += distY;
    comp[2] += distZ;
  }




  // ---------------------------------------------------------------------------
  // G.FERRER - Dece.1999 - Creation
  // Distance entre deux points
  Scalar Point::DistanteDe(const Point &point) const
  {
    Scalar a = comp[0] - point[0];
    Scalar b = comp[1] - point[1];
    Scalar c = comp[2] - point[2];
    return (sqrt(a*a + b*b + c*c));
  }




  // ---------------------------------------------------------------------------
  Scalar Point::DistanteDe(const Scalar *point) const
  {
    Scalar a = comp[0] - point[0];
    Scalar b = comp[1] - point[1];
    Scalar c = comp[2] - point[2];
    return (sqrt(a*a + b*b + c*c));
  }




  // ---------------------------------------------------------------------------
  Scalar Point::DistanteDe(Scalar x, Scalar y, Scalar z) const
  {
    Scalar a = comp[0] - x;
    Scalar b = comp[1] - y;
    Scalar c = comp[2] - z;
    return (sqrt(a*a + b*b + c*c));
  }
}
