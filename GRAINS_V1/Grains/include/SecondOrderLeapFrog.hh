#ifndef _SECONDORDERLEAPFROG_HH_
#define _SECONDORDERLEAPFROG_HH_

#include "TimeIntegrator.hh"


/** @brief The class SecondOrderExplicit.

    Second order Leap Frog integration scheme: v(t+dt/2)=v(t-dt/2)+dt*a(t) and 
    x(t+dt)=x(t)+dt*v(t+dt/2).

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class SecondOrderLeapFrog : public TimeIntegrator
{
  public:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Destructor */
    ~SecondOrderLeapFrog();
  
    /** @brief Default constructor */
    SecondOrderLeapFrog();  
    //@}


    /** @name Methods */
    //@{
    /** @brief Creates and returns a clone of the time integrator */
    TimeIntegrator* clone() const ;

    /** @brief Computes the new velocity and position at time t
    @param particle particle
    @param kine particle kinematics
    @param coupling_factor coupling factor 
    @param torque_bf torque exerted on the particle in body-fixed coordinates 
    system
    @param vtrans translational velocity 
    @param transMotion translation motion
    @param vrot angular velocity in body-fixed coordinates system 
    @param meanVRot average angular velocity in body-fixed coordinates system 
    in interval [t,t-dt]
    @param dt_particle_vel velocity time step magnitude 
    @param dt_particle_disp motion time step magnitude */        
    void Move( Particle* particle, ParticleKinematics* kine,
	double const& coupling_factor, Vector3 const& torque_bf,
	Vector3& vtrans, Vector3& transMotion, 
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp );

    /** @brief Advances velocity over dt_particle_vel
    @param particle particle
    @param kine particle kinematics
    @param coupling_factor coupling factor
    @param torque_bf torque exerted on the particle in body-fixed coordinates 
    system     
    @param vtrans translational velocity 
    @param vrot angular velocity in body-fixed coordinates system 
    @param dt_particle_vel velocity time step magnitude */  
    void advanceVelocity( Particle* particle, ParticleKinematics* kine,
	double const& coupling_factor, Vector3 const& torque_bf, 
	Vector3& vtrans, Vector3& vrot, 
	double const& dt_particle_vel );

    /** @brief Writes time integrator data in an output stream with a high
    precision and 2014 format
    @param fileOut output stream 
    @param particle particle related to the time integrator */
    void writeParticleKinematics2014( ostream& fileOut, 
    	Particle const* particle ) const; 
  
    /** @brief Writes time integrator data in an output stream with a binary 
    and 2014 format
    @param fileOut output stream 
    @param particle particle related to the time integrator */
    void writeParticleKinematics2014_binary( ostream& fileOut, 
    	Particle const* particle );

    /** @brief Reads time integrator data from a stream in the 2014 format 
    @param StreamIN input stream 
    @param particle particle related to the time integrator */
    void readParticleKinematics2014( istream& StreamIN,
    	Particle* particle ); 
  
    /** @brief Reads time integrator data from a stream in a binary form in the
    2014 format 
    @param StreamIN input stream
    @param particle particle related to the time integrator */
    void readParticleKinematics2014_binary( istream& StreamIN,
    	Particle* particle );
	
    /** @brief Returns the number of bytes of the time integrator data when 
    written in a binary format to an output stream */
    size_t get_numberOfBytes() const ;		
    //@}


  protected:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Copy constructor
    @param copy copied SecondOrderLeapFrog object */
    SecondOrderLeapFrog( SecondOrderLeapFrog const& copy );
    //@}
};

#endif
