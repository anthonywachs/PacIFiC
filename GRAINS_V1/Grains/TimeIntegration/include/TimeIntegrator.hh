#ifndef _TIMEINTEGRATOR_HH_
#define _TIMEINTEGRATOR_HH_

#include "Vector3.hh"
#include "Quaternion.hh"


/** @brief The class TimeIntegrator.

    Numerical scheme for the time integration of the Newton's law and the
    kinematic equations. 

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class TimeIntegrator
{
  public:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Destructor */
    virtual ~TimeIntegrator();
    //@}


    /** @name Methods */
    //@{
    /** @brief Creates and returns a clone of the time integrator */
    virtual TimeIntegrator* clone() const = 0;

    /** @brief Computes the new velocity and position at time t+dt
    @param dUdt Translational acceleration dU/dt
    @param vtrans translational velocity 
    @param transMotion translation motion
    @param dOmegadt Angular ecceleration dom/dt
    @param vrot angular velocity 
    @param meanVRot average angular velocity in interval [t,t+dt]
    @param dt_particle_vel velocity time step magnitude 
    @param dt_particle_disp motion time step magnitude */        
    virtual void Move( Vector3 const& dUdt, Vector3& vtrans, 
	Vector3& transMotion, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp ) = 0;

    /** @brief Advances velocity over dt_particle_vel
    @param dUdt Translational acceleration dU/dt
    @param vtrans translational velocity 
    @param dOmegadt Angular ecceleration dom/dt
    @param vrot angular velocity 
    @param dt_particle_vel velocity time step magnitude */
    virtual void advanceVelocity( Vector3 const& dUdt, Vector3& vtrans, 
	Vector3 const& dOmegadt, Vector3& vrot, double const& dt_particle_vel );

    /** @brief Copies kinematics at time t-2dt (translational velocity, angular 
    velocity, variation of translational velocity, variation of angular 
    velocity) in a 1D array 
    @param vit 1D array where kinematics at time t-2dt is copied
    @param i start index to copy in the 1D array */
    virtual void copyKinematicsNm2( double* vit, int i ) const;
    
    /** @brief Writes time integrator data in an output stream with a high
    precision and 2014 format
    @param fileOut output stream 
    @param dUdt particle translational acceleration 
    @param dOmegadt particle angular acceleration */
    virtual void writeParticleKinematics2014( ostream& fileOut,
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
    virtual void readParticleKinematics2014( istream& StreamIN,
    	Vector3& dUdt, Vector3& dOmegadt ); 
  
    /** @brief Reads time integrator data from a stream in a binary form in the
    2014 format 
    @param StreamIN input stream 
    @param dUdt particle translational acceleration 
    @param dOmegadt particle angular acceleration */
    virtual void readParticleKinematics2014_binary( istream& StreamIN,
    	Vector3& dUdt, Vector3& dOmegadt );

    /** @brief Returns the number of bytes of the time integrator data when 
    written in a binary format to an output stream */
    virtual size_t get_numberOfBytes() const ;		    
    //@}


    /** @name Methods Set */
    //@{  
    /** @brief Sets kinematics at time t-2dt from a 1D array of 12 scalars
    (translational velocity, angular velocity, variation of translational
    velocity, variation of angular velocity)
    @param tab 1D array of 4 vectors containing translational velocity, angular
    velocity, variation of translational velocity and variation of angular 
    velocity */
    virtual void setKinematicsNm2( double const* tab );   
    //@}


  protected:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    TimeIntegrator();

    /** @brief Copy constructor
    @param copy copied TimeIntegrator object */
    TimeIntegrator( TimeIntegrator const& copy );
    //@}
};

#endif
