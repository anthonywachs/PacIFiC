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
    @param vtrans translational velocity at time t
    @param dUdt Translational velocity variation dU/dt
    @param transDisplacement translation displacement
    @param dOmegadt Angular velocity variation dom/dt
    @param vrot angular velocity at time t 
    @param meanVRot average angular velocity in interval [t,t+dt]
    @param dt time step magnitude */        
    void Move( Vector3& vtrans, Vector3 const& dUdt,
	Vector3& transDisplacement, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double dt ) ;

    /** @brief Copies kinematics at time t-2dt (translational velocity, angular 
    velocity, variation of translational velocity, variation of angular 
    velocity) in a 1D array 
    @param vit 1D array where kinematics at time t-2dt is copied
    @param i start index to copy in the 1D array */
    void copyKinematicsNm2( double* vit, int i ) const;
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
    Vector3 m_dUdt_nm2; /**< Translational velocity variation at t-dt */
    Vector3 m_dOmegadt_nm2; /**< Angular velocity variation at t-dt */
    //@}      
};

#endif
