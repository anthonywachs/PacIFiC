#include "AppFluide_Puits.H"
#include "Torseur.H"


// ----------------------------------------------------------------------------
// Constructeur
// G.FERRER - Juil.2003 - Creation
AppFluide_Puits::AppFluide_Puits() :
  App()
{
}




// ----------------------------------------------------------------------------
// Constructeur
// G.FERRER - Juil.2003 - Creation
AppFluide_Puits::AppFluide_Puits( istream &fileIn ) :
  App()
{
  read( fileIn );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
AppFluide_Puits::AppFluide_Puits( DOMNode* root ) :
  App()
{
  // Point
  DOMNode* nPoint = ReaderXML::getNode(root, "Puits");
  puits[X] = ReaderXML::getNodeAttr_Double(nPoint, "X");
  puits[Y] = ReaderXML::getNodeAttr_Double(nPoint, "Y");
  puits[Z] = ReaderXML::getNodeAttr_Double(nPoint, "Z");

  // Force au puits
  DOMNode* nForce = ReaderXML::getNode(root, "ForceAuPuits");
  f0 = ReaderXML::getNodeAttr_Double(nForce, "f");

  // Distance critique
  DOMNode* nDist = ReaderXML::getNode(root, "DistanceCritique");
  distance = ReaderXML::getNodeAttr_Double(nDist, "D");

  // Gradient de pression
  DOMNode* nGradP = ReaderXML::getNode(root, "Equilibre");
  deltaP = ReaderXML::getNodeAttr_Double(nGradP, "dP");
}




// ----------------------------------------------------------------------------
// Destructeur
// G.FERRER - Juil.2003 - Creation
AppFluide_Puits::~AppFluide_Puits()
{
}




// ----------------------------------------------------------------------------
// Calcul des forces hydrodynamique exercees sur les particules
// G.FERRER - Juil.2003 - Creation
void AppFluide_Puits::CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules )
{
  Vecteur direction, om;
  Scalar dist;
  list<Particule*>::const_iterator particule;

  for (particule=particules->begin(); particule!=particules->end();
       particule++)
  {
    direction = puits - *(*particule)->getPosition();
    dist = Norm(direction);
    direction /= dist;
    if ( f0 - dist/distance * deltaP > 0. )
      om = (f0 - dist/distance * deltaP) * (*particule)->getVolume()
      	* direction;
    (*particule)->addBodyForce(om);
  }
}




// ----------------------------------------------------------------------------
// Lecture des donnees fluide
// G.FERRER - Mai .2003 Creation
void AppFluide_Puits::read( istream &fileIn )
{
  string cle;
  fileIn >> cle >> puits >> f0;
  fileIn >> cle >> distance >> deltaP;
}
