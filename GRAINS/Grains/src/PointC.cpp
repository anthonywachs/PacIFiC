/* D.PETIT - Creation
   
   Complement de l'algorithme GJK - Engine
   Nouvelle classe de convexe : Le Point .
*/
#include "PointC.H"

// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
PointC::PointC() 
{
}


// ----------------------------------------------------------------------------
// Constructeur avec lecture des donnees
// G.FERRER - Juil.2001 - Creation
PointC::PointC(istream &fileIn)
{
  readClass(fileIn);
}

// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
PointC::~PointC() 
{
}


// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// Renvoie le choix de retrecissement
int PointC:: getShrinkingChoice()const
{
  return 0;
}




// ----------------------------------------------------------------------------
// M. SULAIMAN - Nov.2015 - Creation
// Set the shrinking radius value
void PointC::setShrinkingRadius(Scalar CurrentRadius)
{
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine l'inertie et l'inertie inverse d'un Point
bool PointC::BuildInertie(Scalar *inertie,Scalar *inertie_1) const
{
  inertie[0]=inertie[1]=inertie[2]=inertie[3]=inertie[4]
    =inertie[5]=0.0;
  inertie_1[0]=inertie_1[1]=inertie_1[2]=inertie_1[3]=inertie_1[4]
    =inertie_1[5]=0.0;
  return true;
}


// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
//  Determine le rayon circonscrit d'un point
Scalar PointC::BuildRayonRef() const 
{
  return 0.;
}


// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction d'un clone d'un Point
Convex* PointC::clone() const 
{
  return new PointC();
}


// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le Volume d'un Point
Scalar PointC::getVolume() const
{
  return 0;
}


// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
Point PointC::support(const Vecteur& v) const 
{
  return Point(0.0, 0.0, 0.0);
}


// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture d'un PointC
void PointC::printClass(ostream &fileOut) const 
{
  fileOut << "*PointC\n";
}

  
// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// lecture d'un PointC
void PointC::readClass(istream &fileIn) 
{
}
