#include "SecondOrderLeapFrog.hh"
#include "ParticleKinematics.hh"
#include "Particle.hh"


// ----------------------------------------------------------------------------
// Default constructeur
SecondOrderLeapFrog::SecondOrderLeapFrog() 
  : TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Copy constructor
SecondOrderLeapFrog::SecondOrderLeapFrog( SecondOrderLeapFrog const& copy ) 
  : TimeIntegrator( copy )
{}




// ----------------------------------------------------------------------------
// Destructor
SecondOrderLeapFrog::~SecondOrderLeapFrog()
{}




// ----------------------------------------------------------------------------
// Creates and returns a clone of the time integrator
TimeIntegrator* SecondOrderLeapFrog::clone() const
{
  return ( new SecondOrderLeapFrog(*this) );
}




// ----------------------------------------------------------------------------
// Computes the new velocity and position at time t+dt
void SecondOrderLeapFrog::Move( Particle* particle, ParticleKinematics* kine,
	double const& coupling_factor, Vector3 const& torque_bf,
	Vector3& vtrans, Vector3& transMotion, 
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{
  // Translational velocity and motion 
  vtrans += ( *(particle->getForce()) / 
  	( particle->getMass() * coupling_factor ) ) * dt_particle_vel;
  transMotion = vtrans * dt_particle_disp; 
  
  // Angular velocity integrated by a Strong Stability Preserving RK3  
  Vector3 K1, K2, K3;
  kine->computeAngularAccelerationBodyFixed( particle, torque_bf, vrot, K1 );
  K1 *= dt_particle_vel;
  kine->computeAngularAccelerationBodyFixed( particle, torque_bf, vrot + K1, 
  	K2 );
  K2 *= dt_particle_vel;  
  kine->computeAngularAccelerationBodyFixed( particle, torque_bf, 
  	vrot + 0.25 * ( K1 + K2 ), K3 ); 
  K3 *= dt_particle_vel;
  vrot += ( K1 + K2 + 4. * K3 ) / 6.;
  meanVRot = vrot;	  
}




// ----------------------------------------------------------------------------
// Advances velocity over dt_particle_vel
void SecondOrderLeapFrog::advanceVelocity( Particle* particle, 
	ParticleKinematics* kine, double const& coupling_factor, 
	Vector3 const& torque_bf, Vector3& vtrans, Vector3& vrot, 
	double const& dt_particle_vel )
{
  // Translational velocity
  vtrans += ( *(particle->getForce()) / 
  	( particle->getMass() * coupling_factor ) ) * dt_particle_vel;
	
  // Angular velocity integrated by a Strong Stability Preserving RK3
  Vector3 K1, K2, K3;
  kine->computeAngularAccelerationBodyFixed( particle, torque_bf, vrot, K1 );
  K1 *= dt_particle_vel;
  kine->computeAngularAccelerationBodyFixed( particle, torque_bf, vrot + K1, 
  	K2 );
  K2 *= dt_particle_vel;  
  kine->computeAngularAccelerationBodyFixed( particle, torque_bf, 
  	vrot + 0.25 * ( K1 + K2 ), K3 ); 
  K3 *= dt_particle_vel;
  vrot += ( K1 + K2 + 4. * K3 ) / 6.;	  	  
}




// ----------------------------------------------------------------------------
// Writes time integrator data in an output stream with a high
// precision and 2014 format
void SecondOrderLeapFrog::writeParticleKinematics2014( ostream& fileOut,
    	Particle const* particle ) const
{
  fileOut << " "; 
  particle->getForce()->writeGroup3( fileOut ); 
  fileOut << " "; 
  particle->getTorque()->writeGroup3( fileOut );   
} 



  
// ----------------------------------------------------------------------------
// Writes time integrator data in an output stream with a binary and 2014 format
void SecondOrderLeapFrog::writeParticleKinematics2014_binary( 
	ostream& fileOut, Particle const* particle )
{
  Vector3 force = *(particle->getForce()), torque = *(particle->getTorque());
  force.writeGroup3_binary( fileOut ); 
  torque.writeGroup3_binary( fileOut );   
}




// ----------------------------------------------------------------------------
// Reads time integrator data from a stream in the 2014 format 
void SecondOrderLeapFrog::readParticleKinematics2014( istream& StreamIN,
    	Particle* particle )
{
  Vector3 force, torque;
  StreamIN >> force >> torque;
  particle->setForce( force );
  particle->setTorque( torque );
} 



  
// ----------------------------------------------------------------------------
// Reads time integrator data from a stream in a binary form in the 2014 format 
void SecondOrderLeapFrog::readParticleKinematics2014_binary( 
	istream& StreamIN, Particle* particle )
{
  Vector3 force, torque;
  force.readGroup3_binary( StreamIN ); 
  torque.readGroup3_binary( StreamIN ); 
  particle->setForce( force );
  particle->setTorque( torque );    
}




// ----------------------------------------------------------------------------
// Returns the number of bytes of the time integrator data when written in a 
// binary format to an output stream
size_t SecondOrderLeapFrog::get_numberOfBytes() const
{
  return ( 2 * solid::Group3::m_sizeofGroup3 );
}
