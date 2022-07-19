#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "Error.hh"
#include "Component.hh"


// ContactError
// ----------------------------------------------------------------------------
// Default contructor
ContactError::ContactError() 
  : m_id0( NULL ) 
  , m_id1( NULL ) 
{}




// ----------------------------------------------------------------------------
// Destructor
ContactError::~ContactError()
{}




// ----------------------------------------------------------------------------
// Outputs message when exception is caught
void ContactError::Message( ostream& fileOut ) const
{
  fileOut << "ERR Contact : " << m_message << endl; 

  fileOut << "Component 0 : " << endl;
  fileOut << "  Numero = " << m_id0->getMasterComponent()->getID() << endl;
  //fileOut << "  Numero = " << id0->getID() << endl;
//   if ( m_id0->getID() == -2 ) 
//     fileOut << "  Numero de la particle de reference = " 
//    	<< m_id0->getPeriodicReferenceID() << endl;
  fileOut << "  Classe = ";
  if ( m_id0->getGeometricType() != -100 ) 
    fileOut << m_id0->getGeometricType() << endl;
  else fileOut << "Obstacle" << endl;
  fileOut << "  Rayon d'interaction = " << m_id0->getCrustThickness() << endl;  
  fileOut << "  Position = " << *(m_id0->getPosition());
  fileOut << "  Velocity de translation = " << 
  	*(m_id0->getTranslationalVelocity());
  fileOut << "    Norme = " << Norm( *(m_id0->getTranslationalVelocity()) ) 
  	<< endl;  
  fileOut << "  Velocity de rotation = " << *(m_id0->getAngularVelocity());    
  fileOut << "    Norme = " << Norm( *(m_id0->getAngularVelocity()) ) << endl;    

  fileOut << "Component 1 : " << endl;
  fileOut << "  Numero = " << m_id1->getMasterComponent()->getID() << endl; 
  //fileOut << "  Numero = " << id1->getID() << endl; 
//   if ( m_id1->getID() == -2 ) 
//     fileOut << "  Numero de la particle de reference = " 
//    	<< m_id1->getPeriodicReferenceID() << endl
// 	<< "  Nb de periodes = " << m_id1->getNbPeriodes() << endl;   
  fileOut << "  Classe = ";
  if ( m_id1->getGeometricType() != -100 ) 
    fileOut << m_id1->getGeometricType() << endl;
  else fileOut << "Obstacle" << endl;
  fileOut << "  Rayon d'interaction = " << m_id1->getCrustThickness() << endl;
  fileOut << "  Position = " << *(m_id1->getPosition());
  fileOut << "  Velocity de translation = " << 
  	*(m_id1->getTranslationalVelocity());
  fileOut << "    Norme = " << Norm( *(m_id1->getTranslationalVelocity()) ) 
  	<< endl;  
  fileOut << "  Velocity de rotation = " << *(m_id1->getAngularVelocity());    
  fileOut << "    Norme = " << Norm( *(m_id1->getAngularVelocity()) ) << endl; 
  
  fileOut << "Distance entre GC = " << 
  	Norm( *m_id0->getPosition() - *m_id1->getPosition() ) << endl;
  fileOut << "Penetration max autorise = " << 
  	m_id0->getCrustThickness() + m_id1->getCrustThickness() << endl;
}




// ----------------------------------------------------------------------------
// Sets the pointers to the components invovled in the contact error
// and the time of the contact error 
void ContactError::setComponents( Component* id0_, Component* id1_,
	double time_ ) 
{
  m_id0 = id0_;
  m_id1 = id1_;
  m_time = time_;
}




// ----------------------------------------------------------------------------
// Sets header message
void ContactError::setMessage( string const& mes )
{
  m_message = mes;
}




// ----------------------------------------------------------------------------
// Returns the pointers to the 2 components involved in the contact 
// error as a list for further post-processing
list<Component*> ContactError::getComponents()
{
  list<Component*> lc;
  lc.push_back(m_id0);
  lc.push_back(m_id1);
  
  return ( lc );
}








// DisplacementError
// ----------------------------------------------------------------------------
// Constructor with input parameters
DisplacementError::DisplacementError( Component* id0_, double depl, 
	double deplMax, double time_ ) 
  : m_id0( id0_ ) 
  , m_depl( depl ) 
  , m_deplMax( deplMax ) 
  , m_time( time_ ) 
{}




// ----------------------------------------------------------------------------
DisplacementError::DisplacementError() 
{}



  
// ----------------------------------------------------------------------------
DisplacementError::~DisplacementError() 
{}




// ----------------------------------------------------------------------------
// Outputs message when exception is caught
void DisplacementError::Message( ostream& fileOut ) const
{
  fileOut << "ERR Deplacement : " << m_depl 
	<< " pour " << m_deplMax << " autorise a t=" 
	<< GrainsExec::doubleToString(m_time,TIMEFORMAT) << endl;
  fileOut << "Component : " << endl;
  fileOut << "  Numero = " << m_id0->getID() << endl;
//   if ( m_id0->getID() == -2 ) 
//     fileOut << "  Numero de la particle de reference = " 
//    	<< m_id0->getPeriodicReferenceID() << endl;
  fileOut << "  Classe = ";
  if ( m_id0->getGeometricType() != -100 ) 
    fileOut << m_id0->getGeometricType() << endl;
  else fileOut << "Obstacle" << endl;
  fileOut << "  Position = " << *m_id0->getPosition(); 
  fileOut << "  Velocity de translation = " << 
  	*(m_id0->getTranslationalVelocity());
  fileOut << "    Norme = " << Norm( *(m_id0->getTranslationalVelocity()) ) 
  	<< endl;
  fileOut << "  Velocity de rotation = " << *(m_id0->getAngularVelocity());    
  fileOut << "    Norme = " << Norm( *(m_id0->getAngularVelocity()) ) << endl;
}




// ----------------------------------------------------------------------------
// Returns the pointer to the component involved in the displacement
// error in a list for further post-processing
list<Component*> DisplacementError::getComponent()
{
  list<Component*> lc;
  lc.push_back(m_id0);
  
  return ( lc );
}








// SimulationError
// ----------------------------------------------------------------------------
// Default constructor
SimulationError::SimulationError() 
{}




// ----------------------------------------------------------------------------
// Constructor with the name of the method involved in the error 
// as an input parameter
SimulationError::SimulationError( string const& str ) 
  : m_method( str ) 
{}




// ----------------------------------------------------------------------------
// Destructor
SimulationError::~SimulationError() 
{}




// ----------------------------------------------------------------------------
// Outputs message when exception is caught
void SimulationError::Message( ostream& fileOut ) const
{
  fileOut << "ERR simulation : in function " << m_method << '\n'
	<< "Stop execution\n";
}
