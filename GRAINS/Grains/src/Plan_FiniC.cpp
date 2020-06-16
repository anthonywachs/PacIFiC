/* D.PETIT - Creation
   
   Complement de l'algorithme GJK - Engine
   Nouvelle classe de conevexe : Le segment .
*/
#include "Plan_FiniC.H"

// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Plan_FiniC::Plan_FiniC(const Scalar x, const Scalar y) : 
  projx(x),projy(y) 
{
}

// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction avec un fichier
Plan_FiniC::Plan_FiniC(istream &fileIn) 
{
  readClass(fileIn);
}


// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Plan_FiniC::~Plan_FiniC() 
{
}


// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine l'inertie et l'inertie inverse d'un Segment
bool Plan_FiniC::BuildInertie(Scalar *inertie,Scalar *inertie_1) const
{
  inertie[0] = inertie[1] = inertie[2] = 
    inertie[3] = inertie[4] = inertie[5] = 0.0;
  inertie_1[0] = inertie_1[1] = inertie_1[2] = 
    inertie_1[3] = inertie_1[4] = inertie_1[5] = 0.0;
  return true;
}


// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
//  Determine le rayon circonscrit d'un Plan fini
Scalar Plan_FiniC::BuildRayonRef() const 
{
  return sqrt(projx*projx + projy*projy);
}


// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Construction d'un clone d'un Plan fini
Convex* Plan_FiniC::clone() const 
{
  return new Plan_FiniC(projx,projy);
}


// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// Renvoie le choix de retrecissement
int Plan_FiniC::getShrinkingChoice()const
{

return(0);

}



// ----------------------------------------------------------------------------
// M. SULAIMAN - Nov.2015 - Creation
// Set the shrinking radius value
void Plan_FiniC:: setShrinkingRadius(Scalar CurrentRadius)
{

}

// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Scalar Plan_FiniC::getProjx() const 
{
  return projx;
}

// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Scalar Plan_FiniC::getProjy() const 
{
  return projy;
}


// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le Volume d'un Plan
Scalar Plan_FiniC::getVolume() const 
{
  return 0;
}


// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Point Plan_FiniC::support(const Vecteur& v) const 
{
  Scalar norm = Norm(v);
  if (norm > EPSILON) {
    return Point(v[X]>0.0 ? projx:-projx, v[Y]>0.0 ? projy:-projy,0.0);
  } else {
    return Point();
  }
}


// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture d'un Plan fini
void Plan_FiniC::printClass(ostream &fileOut) const 
{
  fileOut << "*Plan_FiniC\n";
  fileOut << projx << '\t' 
	  << projy << '\n';
}


// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture d'un plan fini
void Plan_FiniC::readClass(istream &fileIn) 
{
  fileIn >> projx >> projy;
}
