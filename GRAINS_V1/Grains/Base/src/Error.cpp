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
  fileOut << "  Number = " << m_id0->getMasterComponent()->getID() << endl;
  fileOut << "  Class = ";
  if ( m_id0->getGeometricType() != -100 ) 
    fileOut << m_id0->getGeometricType() << endl;
  else fileOut << "Obstacle" << endl;
  fileOut << "  Crust thickness = " << m_id0->getCrustThickness() << endl;  
  fileOut << "  Position = " << *(m_id0->getPosition()) << endl;
  fileOut << "  Translational velocity = " << 
  	*(m_id0->getTranslationalVelocity());
  fileOut << "    Norm = " << Norm( *(m_id0->getTranslationalVelocity()) ) 
  	<< endl;  
  fileOut << "  Angular velocity = " << *(m_id0->getAngularVelocity());    
  fileOut << "    Norm = " << Norm( *(m_id0->getAngularVelocity()) ) << endl;

  fileOut << "Component 1 : " << endl;
  fileOut << "  Numero = " << m_id1->getMasterComponent()->getID() << endl;    
  fileOut << "  Class = ";
  if ( m_id1->getGeometricType() != -100 ) 
    fileOut << m_id1->getGeometricType() << endl;
  else fileOut << "Obstacle" << endl;
  fileOut << "  Crust thickness = " << m_id1->getCrustThickness() << endl;
  fileOut << "  Position = " << *(m_id1->getPosition()) << endl;
  fileOut << "  Translational velocity = " << 
  	*(m_id1->getTranslationalVelocity());
  fileOut << "    Norm = " << Norm( *(m_id1->getTranslationalVelocity()) ) 
  	<< endl;  
  fileOut << "  Angular velocity = " << *(m_id1->getAngularVelocity());    
  fileOut << "    Norm = " << Norm( *(m_id1->getAngularVelocity()) ) << endl; 
  
  fileOut << "Distance between GC = " << 
  	Norm( *m_id0->getPosition() - *m_id1->getPosition() ) << endl;
  fileOut << "Max overlap allowed = " << 
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








// MotionError
// ----------------------------------------------------------------------------
// Constructor with input parameters
MotionError::MotionError( Component* id0_, double depl, 
	double deplMax, double time_ ) 
  : m_id0( id0_ ) 
  , m_depl( depl ) 
  , m_deplMax( deplMax ) 
  , m_time( time_ ) 
{}




// ----------------------------------------------------------------------------
MotionError::MotionError() 
{}



  
// ----------------------------------------------------------------------------
MotionError::~MotionError() 
{}




// ----------------------------------------------------------------------------
// Outputs message when exception is caught
void MotionError::Message( ostream& fileOut ) const
{
  fileOut << "ERR Motion : " << m_depl 
	<< " compared to " << m_deplMax << " allowed at t=" 
	<< GrainsExec::doubleToString(m_time,TIMEFORMAT) << endl;
  fileOut << "Component : " << endl;
  fileOut << "  Number = " << m_id0->getID() << endl;
  fileOut << "  Class = ";
  if ( m_id0->getGeometricType() != -100 ) 
    fileOut << m_id0->getGeometricType() << endl;
  else fileOut << "Obstacle" << endl;
  fileOut << "  Position = " << *m_id0->getPosition() << endl;; 
  fileOut << "  Translational velocity = " << 
  	*(m_id0->getTranslationalVelocity());
  fileOut << "    Norm = " << Norm( *(m_id0->getTranslationalVelocity()) ) 
  	<< endl;
  fileOut << "  Angular velocity = " << *(m_id0->getAngularVelocity());    
  fileOut << "    Norm = " << Norm( *(m_id0->getAngularVelocity()) ) << endl;
}




// ----------------------------------------------------------------------------
// Returns the pointer to the component involved in the motion
// error in a list for further post-processing
list<Component*> MotionError::getComponent()
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
