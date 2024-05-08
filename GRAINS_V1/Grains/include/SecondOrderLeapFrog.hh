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

    /** @brief Computes the new velocity and position at time t+dt
    @param dUdt Translational acceleration dU/dt
    @param vtrans translational velocity 
    @param transDisplacement translation displacement
    @param dOmegadt Angular ecceleration dom/dt
    @param vrot angular velocity 
    @param meanVRot average angular velocity in interval [t,t+dt]
    @param dt_particle_vel velocity time step magnitude 
    @param dt_particle_disp displacement time step magnitude */        
    void Move( Vector3 const& dUdt, Vector3& vtrans, 
	Vector3& transDisplacement, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp );

    /** @brief Advances velocity over dt_particle_vel
    @param dUdt Translational acceleration dU/dt
    @param vtrans translational velocity 
    @param dOmegadt Angular ecceleration dom/dt
    @param vrot angular velocity 
    @param dt_particle_vel velocity time step magnitude */
    void advanceVelocity( Vector3 const& dUdt, Vector3& vtrans, 
	Vector3 const& dOmegadt, Vector3& vrot, double const& dt_particle_vel );

    /** @brief Writes time integrator data in an output stream with a high
    precision and 2014 format
    @param fileOut output stream 
    @param dUdt particle translational acceleration 
    @param dOmegadt particle angular acceleration */
    void writeParticleKinematics2014( ostream& fileOut,
    	Vector3 const& dUdt, Vector3 const& dOmegadt ) const; 
  
    /** @brief Writes time integrator data in an output stream with a binary 
    and 2014 format
    @param fileOut output stream 
    @param dUdt particle translational acceleration 
    @param dOmegadt particle angular acceleration */
    virtual void writeParticleKinematics2014_binary( ostream& fileOut,
    	Vector3& dUdt, Vector3& dOmegadt );

    /** @brief Reads time integrator data from a stream in the 2014 format 
    @param StreamIN input stream 
    @param dUdt particle translational acceleration 
    @param dOmegadt particle angular acceleration */
    void readParticleKinematics2014( istream& StreamIN,
    	Vector3& dUdt, Vector3& dOmegadt ); 
  
    /** @brief Reads time integrator data from a stream in a binary form in the
    2014 format 
    @param StreamIN input stream 
    @param dUdt particle translational acceleration 
    @param dOmegadt particle angular acceleration */
    void readParticleKinematics2014_binary( istream& StreamIN,
    	Vector3& dUdt, Vector3& dOmegadt );
	
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
