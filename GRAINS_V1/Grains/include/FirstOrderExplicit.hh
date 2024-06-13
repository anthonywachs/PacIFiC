#ifndef _FIRSTORDEREXPLICIT_HH_
#define _FIRSTORDEREXPLICIT_HH_

#include "TimeIntegrator.hh"


/** @brief The class FirstOrderExplicit.

    First order explicit integration scheme: x(t+dt)=x(t)+dt*v(t) and 
    v(t+dt)=v(t)+dt*a(t). 

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class FirstOrderExplicit : public TimeIntegrator
{
  public:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Destructor */
    ~FirstOrderExplicit();
  
    /** @brief Default constructor */
    FirstOrderExplicit();  
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
    //@}


  protected:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Copy constructor
    @param copy copied FirstOrderExplicit object */
    FirstOrderExplicit( FirstOrderExplicit const& copy );
    //@}
};

#endif
