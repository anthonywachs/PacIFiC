#include "SecondOrderAdamsBashforth.hh"


// ----------------------------------------------------------------------------
// Default constructor
SecondOrderAdamsBashforth::SecondOrderAdamsBashforth() 
  : TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Copy constructor
SecondOrderAdamsBashforth::SecondOrderAdamsBashforth(
	 SecondOrderAdamsBashforth const& copy ) 
  : TimeIntegrator( copy )
{
  m_translationalVelocity_nm2 = copy.m_translationalVelocity_nm2;
  m_angularVelocity_nm2 = copy.m_angularVelocity_nm2;
  m_dUdt_nm2 = copy.m_dUdt_nm2;
  m_dOmegadt_nm2 = copy.m_dOmegadt_nm2;
}




// ----------------------------------------------------------------------------
// Destructor
SecondOrderAdamsBashforth::~SecondOrderAdamsBashforth()
{}




// ----------------------------------------------------------------------------
// Creates and returns a clone of the time integrator
TimeIntegrator* SecondOrderAdamsBashforth::clone() const
{
  return ( new SecondOrderAdamsBashforth(*this) );
}




// ----------------------------------------------------------------------------
// Computes the new velocity and position at time t+dt
void SecondOrderAdamsBashforth::Move( Vector3& vtrans, Vector3 const& dUdt,
	Vector3& transDisplacement, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double dt )
{      
  // Velocity et deplacement translationnels
  transDisplacement = 0.5 * ( 3.*vtrans 
  	- m_translationalVelocity_nm2 ) * dt;
  m_translationalVelocity_nm2 = vtrans;
  vtrans += 0.5 * ( 3.*dUdt - m_dUdt_nm2 ) * dt;
  m_dUdt_nm2 = dUdt;  

  // Velocity et deplacement rotationnels
  meanVRot = 0.5 * ( 3.*vrot  - m_angularVelocity_nm2 );
  m_angularVelocity_nm2 = vrot;  
  vrot += 0.5 * ( 3.*dOmegadt - m_dOmegadt_nm2 ) * dt;
  m_dOmegadt_nm2 = dOmegadt;  
}




// ----------------------------------------------------------------------------
// Copies kinematics at time t-2dt (translational velocity, angular 
// velocity, variation of translational velocity, variation of angular 
// velocity) in a 1D array 
void SecondOrderAdamsBashforth::copyKinematicsNm2(double *vit,int i) const
{
  for (int j=0 ;j<3; j++) vit[i+j] = m_translationalVelocity_nm2[j];
  for (int j=0 ;j<3; j++) vit[i+3+j] = m_angularVelocity_nm2[j];  
  for (int j=0 ;j<3; j++) vit[i+6+j] = m_dUdt_nm2[j];  
  for (int j=0 ;j<3; j++) vit[i+9+j] = m_dOmegadt_nm2[j];  
}




// ----------------------------------------------------------------------------
// Sets kinematics at time t-2dt from a 1D array of 12 scalars
// (translational velocity, angular velocity, variation of translational
// velocity, variation of angular velocity)
void SecondOrderAdamsBashforth::setKinematicsNm2(double const* tab)
{
  for (int j=0 ;j<3; j++) m_translationalVelocity_nm2[j] = tab[j];
  for (int j=0 ;j<3; j++) m_angularVelocity_nm2[j] = tab[j+3];  
  for (int j=0 ;j<3; j++) m_dUdt_nm2[j] = tab[j+6];  
  for (int j=0 ;j<3; j++) m_dOmegadt_nm2[j] = tab[j+9];   
} 
