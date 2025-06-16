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
