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
