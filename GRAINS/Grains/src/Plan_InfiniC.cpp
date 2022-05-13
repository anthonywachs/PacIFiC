#include "Plan_InfiniC.H"

// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
Plan_InfiniC::Plan_InfiniC() : 
  Plan_FiniC(INFINITY, INFINITY) 
{
}


// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
Plan_InfiniC::Plan_InfiniC(istream &fileIn) : 
  Plan_FiniC(INFINITY, INFINITY) 
{
}


// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
Plan_InfiniC::~Plan_InfiniC() 
{
}


// ----------------------------------------------------------------------------
// M.SULAIMAN - Nov.2015 - Creation
// Renvoie le choix de retrecissement
int Plan_InfiniC::getShrinkingChoice()const
{

return(0);

}



// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
//  Determine le rayon circonscrit d'un plan infini
Scalar Plan_InfiniC::BuildRayonRef() const 
{
  return INFINITY;
}

// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// ecriture d'un plan infini
void Plan_InfiniC::printClass(ostream &fileOut) const 
{
  fileOut << "*Plan_InfiniC\n";
}

// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// lecture d'un plan infini
void Plan_InfiniC::readClass(istream &fileIn) 
{
}
