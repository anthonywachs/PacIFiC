#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "Erreur.H"
#include "Composant.H"


// ----------------------------------------------------------------------------
// Contructeur par defaut
ErreurContact::ErreurContact() : 
  id0( NULL ), 
  id1( NULL ) 
{}




// ----------------------------------------------------------------------------
// Destructeur
ErreurContact::~ErreurContact()
{}




// ----------------------------------------------------------------------------
// Message de l'exception : identificateur des composants en contact 
void ErreurContact::Message( ostream &fileOut ) const
{
  fileOut << "ERR Contact : " << message << endl; 

  fileOut << "Composant 0 : " << endl;
  fileOut << "  Numero = " << id0->ReferenceComposant()->getID() << endl;
  //fileOut << "  Numero = " << id0->getID() << endl;
  if ( id0->getID() == -2 ) 
    fileOut << "  Numero de la particule de reference = " 
   	<< id0->getPeriodicReferenceID() << endl;
  fileOut << "  Classe = ";
  if ( id0->getParticuleClasse() != -100 ) 
    fileOut << id0->getParticuleClasse() << endl;
  else fileOut << "Obstacle" << endl;
  fileOut << "  Rayon d'interaction = " << id0->getRayonInteraction() << endl;  
  fileOut << "  Position = " << *(id0->getPosition());
  fileOut << "  Vitesse de translation = " << *(id0->getVitesseTranslation());
  fileOut << "    Norme = " << Norm( *(id0->getVitesseTranslation()) ) << endl;  
  fileOut << "  Vitesse de rotation = " << *(id0->getVitesseRotation());    
  fileOut << "    Norme = " << Norm( *(id0->getVitesseRotation()) ) << endl;    

  fileOut << "Composant 1 : " << endl;
  fileOut << "  Numero = " << id1->ReferenceComposant()->getID() << endl; 
  //fileOut << "  Numero = " << id1->getID() << endl; 
  if ( id1->getID() == -2 ) 
    fileOut << "  Numero de la particule de reference = " 
   	<< id1->getPeriodicReferenceID() << endl
	<< "  Nb de periodes = " << id1->getNbPeriodes() << endl;   
  fileOut << "  Classe = ";
  if ( id1->getParticuleClasse() != -100 ) 
    fileOut << id1->getParticuleClasse() << endl;
  else fileOut << "Obstacle" << endl;
  fileOut << "  Rayon d'interaction = " << id1->getRayonInteraction() << endl;
  fileOut << "  Position = " << *(id1->getPosition());
  fileOut << "  Vitesse de translation = " << *(id1->getVitesseTranslation());
  fileOut << "    Norme = " << Norm( *(id1->getVitesseTranslation()) ) << endl;  
  fileOut << "  Vitesse de rotation = " << *(id1->getVitesseRotation());    
  fileOut << "    Norme = " << Norm( *(id1->getVitesseRotation()) ) << endl; 
  
  fileOut << "Distance entre GC = " << 
  	Norm( *id0->getPosition() - *id1->getPosition() ) << endl;
  fileOut << "Penetration max autorise = " << 
  	id0->getRayonInteraction() + id1->getRayonInteraction() << endl;
}




// ----------------------------------------------------------------------------
// Affectation des identificateurs des composants en contact
void ErreurContact::setComposants( Composant* id0_, Composant* id1_,
	double temps_ ) 
{
  id0 = id0_;
  id1 = id1_;
  m_temps = temps_;
}




// ----------------------------------------------------------------------------
// Affectation du message decrivant l'exception
void ErreurContact::setMessage( const string &mes )
{
  message = mes;
}




// ----------------------------------------------------------------------------
// Renvoie les 2 composants dans une liste pour post-processing
list<Composant*> ErreurContact::getComposants()
{
  list<Composant*> lc;
  lc.push_back(id0);
  lc.push_back(id1);
  
  return lc;
}







// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Constructeur
ErreurDeplacement::ErreurDeplacement( Composant* id0_, double depl, 
	double deplMax, double temps_ ) : 
  id0( id0_ ), 
  m_depl( depl ), 
  m_deplMax( deplMax ), 
  m_temps( temps_ ) 
{}



  
// ----------------------------------------------------------------------------
ErreurDeplacement::~ErreurDeplacement() 
{}




// ----------------------------------------------------------------------------
// Message de l'exception : identificateur des composants en contact 
void ErreurDeplacement::Message( ostream &fileOut ) const
{
  fileOut << "ERR Deplacement : " << m_depl 
	<< " pour " << m_deplMax << " autorise a t=" 
	<< Grains_Exec::doubleToString(m_temps,TIMEFORMAT) << endl;
  fileOut << "Composant : " << endl;
  fileOut << "  Numero = " << id0->getID() << endl;
  if ( id0->getID() == -2 ) 
    fileOut << "  Numero de la particule de reference = " 
   	<< id0->getPeriodicReferenceID() << endl;
  fileOut << "  Classe = ";
  if ( id0->getParticuleClasse() != -100 ) 
    fileOut << id0->getParticuleClasse() << endl;
  else fileOut << "Obstacle" << endl;
  fileOut << "  Position = " << *id0->getPosition(); 
  fileOut << "  Vitesse de translation = " << *(id0->getVitesseTranslation());
  fileOut << "    Norme = " << Norm( *(id0->getVitesseTranslation()) ) << endl;
  fileOut << "  Vitesse de rotation = " << *(id0->getVitesseRotation());    
  fileOut << "    Norme = " << Norm( *(id0->getVitesseRotation()) ) << endl;
}




// ----------------------------------------------------------------------------
// Renvoie le 2 composant dans une liste pour post-processing
list<Composant*> ErreurDeplacement::getComposant()
{
  list<Composant*> lc;
  lc.push_back(id0);
  
  return lc;
}








// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Contructeur par defaut
ErreurSimulation::ErreurSimulation() 
{}




// ----------------------------------------------------------------------------
// Constructeur avec indication de la methode de levee d'erreur 
ErreurSimulation::ErreurSimulation( const string &str ) : 
  methode( str ) 
{}




// ----------------------------------------------------------------------------
// Contructeur par defaut
ErreurSimulation::~ErreurSimulation() 
{}




// ----------------------------------------------------------------------------
//  Message d'erreur pour arret simulation 
void ErreurSimulation::Message( ostream &fileOut ) const
{
  fileOut << "ERR simulation : in function " << methode << '\n'
	<< "Stop execution\n";
}
