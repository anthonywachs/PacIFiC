#ifndef _SECONDORDERADAMSBASHFORTH_HH_
#define _SECONDORDERADAMSBASHFORTH_HH_

#include "TimeIntegrator.hh"


/** @brief The class SecondOrderAdamsBashforth.

    Second order Adams-Bashforth integration scheme:
    x(t+dt)=x(t)+0.5*(3*v(t)-v(t-dt)) and v(t+dt)=v(t)+0.5*(3*a(t)-a(t-dt)).

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class SecondOrderAdamsBashforth : public TimeIntegrator
{
  public:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Destructor */
    ~SecondOrderAdamsBashforth();
  
    /** @brief Default constructor */
    SecondOrderAdamsBashforth();  
    //@}


    /** @name Methods */
    //@{
    /** @brief Creates and returns a clone of the time integrator */
    TimeIntegrator* clone() const ;

    /** @brief Computes the new velocity and position at time t+dt
    @param dUdt Translational acceleration dU/dt
    @param vtrans translational velocity 
    @param transMotion translation motion
    @param dOmegadt Angular ecceleration dom/dt
    @param vrot angular velocity 
    @param meanVRot average angular velocity in interval [t,t+dt]
    @param dt_particle_vel velocity time step magnitude 
    @param dt_particle_disp motion time step magnitude */        
    void Move( Vector3 const& dUdt, Vector3& vtrans, 
	Vector3& transMotion, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp );

    /** @brief Copies kinematics at time t-2dt (translational velocity, angular 
    velocity, variation of translational velocity, variation of angular 
    velocity) in a 1D array 
    @param vit 1D array where kinematics at time t-2dt is copied
    @param i start index to copy in the 1D array */
    void copyKinematicsNm2( double* vit, int i ) const;
    
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
    void writeParticleKinematics2014_binary( ostream& fileOut,
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

    
    /** @name Methods Set */
    //@{  
    /** @brief Sets kinematics at time t-2dt from a 1D array of 12 scalars
    (translational velocity, angular velocity, variation of translational
    velocity, variation of angular velocity)
    @param tab 1D array of 4 vectors containing translational velocity, angular
    velocity, variation of translational velocity and variation of angular 
    velocity */
    void setKinematicsNm2( double const* tab );   
    //@}    


  protected:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Copy constructor
    @param copy copied SecondOrderAdamsBashforth object */
    SecondOrderAdamsBashforth( SecondOrderAdamsBashforth const& copy );
    //@}

    
  private:
    /** @name Parameters */
    //@{
    Vector3 m_translationalVelocity_nm2; /**< Translational velocity at t-dt */
    Vector3 m_angularVelocity_nm2; /**< Angular velocity at t-dt */
    Vector3 m_dUdt_nm2; /**< Translational acceleration at t-dt */
    Vector3 m_dOmegadt_nm2; /**< Angular acceleration at t-dt */
    //@}      
};

#endif
