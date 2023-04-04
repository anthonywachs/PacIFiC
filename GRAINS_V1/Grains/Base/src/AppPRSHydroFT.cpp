#include "AppPRSHydroFT.hh"
#include "AllComponents.hh"
#include "GrainsExec.hh"


// ----------------------------------------------------------------------------
// Default constructor
AppPRSHydroFT::AppPRSHydroFT()
{}




// ----------------------------------------------------------------------------
// Destructeur
AppPRSHydroFT::~AppPRSHydroFT()
{
  m_PRSHydroForce.clear();
  m_PRSHydroTorque.clear();    
}




// ----------------------------------------------------------------------------
// Computes forces and torques exerted on rigid bodies
void AppPRSHydroFT::ComputeForces( double time, double dt,
    	list<Particle*> const* particles )
{
  list<Particle*>::const_iterator particle;
  size_t i = 0;
  
  for( particle=particles->begin(); particle!=particles->end();particle++)
  {
    if ( (*particle)->getTag() < 2 )
    {
      (*particle)->addBodyForce( m_PRSHydroForce[i] );
      (*particle)->addTorque( m_PRSHydroTorque[i] );
      i++;
    }
  }        
}  




// ----------------------------------------------------------------------------
// Allocates the arrays of hydro force and torque 
void AppPRSHydroFT::allocateHydroFT( size_t const& nbPart )
{
  Vector3 vnul;
  m_PRSHydroForce.reserve( nbPart );
  for (size_t i=0;i<nbPart;++i) m_PRSHydroForce.push_back( vnul ); 
  m_PRSHydroTorque.reserve( nbPart );
  for (size_t i=0;i<nbPart;++i) m_PRSHydroTorque.push_back( vnul );      
}




// ----------------------------------------------------------------------------
// Sets the arrays of hydro force and torque 
void AppPRSHydroFT::setHydroFT( vector< vector<double> > const* hydrovec )
{
  size_t i, j, nbPart = hydrovec->size();
  
  
  for (i=0;i<nbPart;++i)
    for (j=0;j<3;++j) 
    {
      m_PRSHydroForce[i][j] = (*hydrovec)[i][j];      
      m_PRSHydroTorque[i][j] = (*hydrovec)[i][j+3];  
    }
}
