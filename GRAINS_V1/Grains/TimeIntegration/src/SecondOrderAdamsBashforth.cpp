#include "SecondOrderAdamsBashforth.hh"
#include "ParticleKinematics.hh"
#include "Particle.hh"
#include "GrainsExec.hh"


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
  m_angularVelocity_bf_nm2 = copy.m_angularVelocity_bf_nm2;
  m_dUdt_nm2 = copy.m_dUdt_nm2;
  m_dOmegadt_bf_nm2 = copy.m_dOmegadt_bf_nm2;
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
void SecondOrderAdamsBashforth::Move( Particle* particle, 
	ParticleKinematics* kine, double const& coupling_factor, 
	Vector3 const& torque_bf, Vector3& vtrans, Vector3& transMotion, 
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{      
  // Note: at first discrete time when the simulation is not reloaded, 
  // values at t^n-2 are unknown and the scheme is First Order Euler only

  // Angular velocity
  if ( GrainsExec::m_time_counter == 1 && !GrainsExec::m_isReloaded )
  {
    // Translational velocity and motion
    Vector3 dUdt_nm1 = *(particle->getForce()) / 
  	( particle->getMass() * coupling_factor );
    transMotion = vtrans * dt_particle_disp;
    m_translationalVelocity_nm2 = vtrans;
    vtrans += dUdt_nm1 * dt_particle_vel;
    m_dUdt_nm2 = dUdt_nm1;
    
    // Angular velocity
    Vector3 dOmegadt_bf_nm1;
    kine->computeAngularAccelerationBodyFixed( particle, torque_bf, vrot, 
  	dOmegadt_bf_nm1 );
    meanVRot = vrot;
    m_angularVelocity_bf_nm2 = vrot; 
    vrot += dOmegadt_bf_nm1 * dt_particle_vel;
    m_dOmegadt_bf_nm2 = dOmegadt_bf_nm1 ;  
  }
  else
  {
    // Translational velocity and motion
    Vector3 dUdt_nm1 = *(particle->getForce()) / 
  	( particle->getMass() * coupling_factor );
    transMotion = 0.5 * ( 3. * vtrans 
  	- m_translationalVelocity_nm2 ) * dt_particle_disp;
    m_translationalVelocity_nm2 = vtrans;
    vtrans += 0.5 * ( 3. * dUdt_nm1 - m_dUdt_nm2 ) * dt_particle_vel;
    m_dUdt_nm2 = dUdt_nm1;

    // Angular velocity
    Vector3 dOmegadt_bf_nm1;
    kine->computeAngularAccelerationBodyFixed( particle, torque_bf, vrot, 
  	dOmegadt_bf_nm1 );
    meanVRot = 0.5 * ( 3. * vrot  - m_angularVelocity_bf_nm2 );
    m_angularVelocity_bf_nm2 = vrot; 
    vrot += 0.5 * ( 3. * dOmegadt_bf_nm1 - m_dOmegadt_bf_nm2 ) 
    	* dt_particle_vel;
    m_dOmegadt_bf_nm2 = dOmegadt_bf_nm1 ;
  }  
}




// ----------------------------------------------------------------------------
// Copies kinematics at time t-2dt (translational velocity, angular 
// velocity, variation of translational velocity, variation of angular 
// velocity) in a 1D array 
void SecondOrderAdamsBashforth::copyKinematicsNm2(double *vit,int i) const
{
  for (int j=0 ;j<3; j++) vit[i+j] = m_translationalVelocity_nm2[j];
  for (int j=0 ;j<3; j++) vit[i+3+j] = m_angularVelocity_bf_nm2[j];  
  for (int j=0 ;j<3; j++) vit[i+6+j] = m_dUdt_nm2[j];  
  for (int j=0 ;j<3; j++) vit[i+9+j] = m_dOmegadt_bf_nm2[j];  
}




// ----------------------------------------------------------------------------
// Sets kinematics at time t-2dt from a 1D array of 12 scalars
// (translational velocity, angular velocity, variation of translational
// velocity, variation of angular velocity)
void SecondOrderAdamsBashforth::setKinematicsNm2(double const* tab)
{
  for (int j=0 ;j<3; j++) m_translationalVelocity_nm2[j] = tab[j];
  for (int j=0 ;j<3; j++) m_angularVelocity_bf_nm2[j] = tab[j+3];  
  for (int j=0 ;j<3; j++) m_dUdt_nm2[j] = tab[j+6];  
  for (int j=0 ;j<3; j++) m_dOmegadt_bf_nm2[j] = tab[j+9];   
} 




// ----------------------------------------------------------------------------
// Writes time integrator data in an output stream with a high
// precision and 2014 format
void SecondOrderAdamsBashforth::writeParticleKinematics2014( ostream& fileOut, 
    	Particle const* particle ) const
{
  fileOut << " "; 
  m_translationalVelocity_nm2.writeGroup3( fileOut ); 
  fileOut << " "; 
  m_angularVelocity_bf_nm2.writeGroup3( fileOut ); 
  fileOut << " ";  
  m_dUdt_nm2.writeGroup3( fileOut ); 
  fileOut << " ";  
  m_dOmegadt_bf_nm2.writeGroup3( fileOut );    
} 



  
// ----------------------------------------------------------------------------
// Writes time integrator data in an output stream with a binary and 2014 format
void SecondOrderAdamsBashforth::writeParticleKinematics2014_binary( 
	ostream& fileOut, Particle const* particle )
{
  m_translationalVelocity_nm2.writeGroup3_binary( fileOut ); 
  m_angularVelocity_bf_nm2.writeGroup3_binary( fileOut );   
  m_dUdt_nm2.writeGroup3_binary( fileOut );  
  m_dOmegadt_bf_nm2.writeGroup3_binary( fileOut );
}




// ----------------------------------------------------------------------------
// Reads time integrator data from a stream in the 2014 format 
void SecondOrderAdamsBashforth::readParticleKinematics2014( istream& StreamIN,
	Particle* particle )
{
  StreamIN >> m_translationalVelocity_nm2 >> m_angularVelocity_bf_nm2
  	>> m_dUdt_nm2 >> m_dOmegadt_bf_nm2;
} 



  
// ----------------------------------------------------------------------------
// Reads time integrator data from a stream in a binary form in the 2014 format 
void SecondOrderAdamsBashforth::readParticleKinematics2014_binary( 
	istream& StreamIN, Particle* particle )
{
  m_translationalVelocity_nm2.readGroup3_binary( StreamIN ); 
  m_angularVelocity_bf_nm2.readGroup3_binary( StreamIN );   
  m_dUdt_nm2.readGroup3_binary( StreamIN );  
  m_dOmegadt_bf_nm2.readGroup3_binary( StreamIN );
}




// ----------------------------------------------------------------------------
// Returns the number of bytes of the time integrator data when written in a 
// binary format to an output stream
size_t SecondOrderAdamsBashforth::get_numberOfBytes() const
{
  return ( 4 * solid::Group3::m_sizeofGroup3 );
}
