#include "Brownian.hh"
#include "GrainsExec.hh"
#include "GrainsBuilderFactory.hh"


// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
Brownian::Brownian( DOMNode* root, int rank, size_t& error )
  : App()
  , m_mag( 0 )
{
  if ( ReaderXML::hasNodeAttr( root, "Magnitude" ) )
  {
    m_mag = ReaderXML::getNodeAttr_Double( root, "Magnitude" );
    if ( rank == 0 ) cout << GrainsExec::m_shift6 
    	<< "Brownian Force magnitude = " << m_mag << endl;
    error = 0;
  }
  else
  {
    if ( rank == 0 ) cout << GrainsExec::m_shift6 
    	<< "Brownian Force: magnitude is mandatory!" << endl;
    error = 1;  
  }  
}




// ----------------------------------------------------------------------------
// Destructor
Brownian::~Brownian()
{}




// ----------------------------------------------------------------------------
// Computes forces and torques exerted on rigid bodies
void Brownian::ComputeForces( double time, double dt,
  	list<Particle*> const* particles )
{
  list<Particle*>::const_iterator particle;
  Vector3 force;
  double alpha = 0, beta = 0;

  if ( GrainsBuilderFactory::getContext() == DIM_2 )
    for (particle=particles->cbegin(); particle!=particles->cend(); particle++)
    {
      alpha = ( double(random()) / double(INT_MAX) ) * 2. * PI;
      force[X] = m_mag * cos( alpha );
      force[Y] = m_mag * sin( alpha ); 
      (*particle)->addBodyForce( force );     
    }
  else
    for (particle=particles->cbegin(); particle!=particles->cend(); particle++)
    {
      alpha = ( double(random()) / double(INT_MAX) ) * 2. * PI;
      beta = ( double(random()) / double(INT_MAX) ) * 2. * PI;      
      force[X] = m_mag * cos( beta ) * cos( alpha );
      force[Y] = m_mag * cos( beta ) * sin( alpha );
      force[Y] = m_mag * sin( beta );       
      (*particle)->addBodyForce( force );     
    }  
}
