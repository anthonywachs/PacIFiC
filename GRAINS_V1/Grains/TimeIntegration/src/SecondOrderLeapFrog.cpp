#include "SecondOrderLeapFrog.hh"


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
void SecondOrderLeapFrog::Move( Vector3 const& dUdt, Vector3& vtrans, 
	Vector3& transMotion, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{
  // Translational velocity and motion
  vtrans += dUdt * dt_particle_vel;
  transMotion = vtrans * dt_particle_disp;

  // Angular velocity and motion
  vrot += dOmegadt * dt_particle_vel;
  meanVRot = vrot;    
}




// ----------------------------------------------------------------------------
// Advances velocity over dt_particle_vel
void SecondOrderLeapFrog::advanceVelocity( Vector3 const& dUdt, Vector3& vtrans,
	Vector3 const& dOmegadt, Vector3& vrot, double const& dt_particle_vel )
{
  // Translational velocity and motion
  vtrans += dUdt * dt_particle_vel;

  // Angular velocity and motion
  vrot += dOmegadt * dt_particle_vel;
}




// ----------------------------------------------------------------------------
// Writes time integrator data in an output stream with a high
// precision and 2014 format
void SecondOrderLeapFrog::writeParticleKinematics2014( ostream& fileOut,
    	Vector3 const& dUdt, Vector3 const& dOmegadt ) const
{
  fileOut << " "; 
  dUdt.writeGroup3( fileOut ); 
  fileOut << " "; 
  dOmegadt.writeGroup3( fileOut );   
} 



  
// ----------------------------------------------------------------------------
// Writes time integrator data in an output stream with a binary and 2014 format
void SecondOrderLeapFrog::writeParticleKinematics2014_binary( 
	ostream& fileOut, Vector3& dUdt, Vector3& dOmegadt )
{
  dUdt.writeGroup3_binary( fileOut ); 
  dOmegadt.writeGroup3_binary( fileOut );   
}




// ----------------------------------------------------------------------------
// Reads time integrator data from a stream in the 2014 format 
void SecondOrderLeapFrog::readParticleKinematics2014( istream& StreamIN,
    	Vector3& dUdt, Vector3& dOmegadt )
{
  StreamIN >> dUdt >> dOmegadt;
} 



  
// ----------------------------------------------------------------------------
// Reads time integrator data from a stream in a binary form in the 2014 format 
void SecondOrderLeapFrog::readParticleKinematics2014_binary( 
	istream& StreamIN, Vector3& dUdt, Vector3& dOmegadt )
{
  dUdt.readGroup3_binary( StreamIN ); 
  dOmegadt.readGroup3_binary( StreamIN );   
}




// ----------------------------------------------------------------------------
// Returns the number of bytes of the time integrator data when written in a 
// binary format to an output stream
size_t SecondOrderLeapFrog::get_numberOfBytes() const
{
  return ( 2 * solid::Group3::m_sizeofGroup3 );
}
