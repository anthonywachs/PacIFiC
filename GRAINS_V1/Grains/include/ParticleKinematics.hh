#ifndef _PARTICLEKINEMATICS_HH_
#define _PARTICLEKINEMATICS_HH_

#include "Kinematics.hh"
#include "WriterXML.hh"
#include "Vector3.hh"
#include "Quaternion.hh"
using namespace solid;

class ParticleKinematicsMemento;
class Particle;
class ParticleKinematics2D;
class ParticleKinematics3D;
class TimeIntegrator;
class Torsor;
class ParticleKinematics;
istream& operator >> ( istream& fileIn, ParticleKinematics& cinematique );
ostream& operator << ( ostream& fileOut, 
	ParticleKinematics const& cinematique );


/** @brief The class ParticleKinematics.

    Manages the particle kinematics.

    @author Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ParticleKinematics : public Kinematics
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Destructor */
    virtual ~ParticleKinematics();
    //@}
  

    /**@name Virtual methods */
    //@{
    /** @brief Returns the class name */
    virtual string className() const = 0;

    /** @brief Returns a clone of ParticleKinematics */
    virtual ParticleKinematics* clone() const = 0;

    /** @brief Computes the angular acceleration in body fixed space
    @param particle particle related to the kinematics 
    @param torque_bf torque in body fixed space
    @param om_bf angular velocity in body-fixed coordinates system 
    @param dOmdt_bf angular acceleration in body fixed space */
    virtual void computeAngularAccelerationBodyFixed( Particle const* particle,
    	Vector3 const& torque_bf, Vector3 const& om_bf, Vector3& dOmdt_bf ) = 0;
			   
    /** @brief Computes explicitly Idw/dt
    @param dw explicit change of angular velocity
    @param inertie particle inertia tensor */
    virtual Vector3 computeExplicitDJomDt( Vector3 const& dw,
	double const* inertie ) const = 0;			   
    //@}


    /** @name Methods */
    //@{
    /** @brief Integrates Newton's law and moves the particle 
    @param particle the particle
    @param dt_particle_vel velocity time step magnitude 
    @param dt_particle_disp motion time step magnitude */
    double Move( Particle* particle, double const& dt_particle_vel, 
    	double const& dt_particle_disp );
	
    /** @brief Advances velocity over dt_particle_vel
    @param particle the particle    
    @param dt_particle_vel velocity time step magnitude */
    void advanceVelocity( Particle* particle, double const& dt_particle_vel ); 	

    /** @brief Returns the total velocity U + om x R given R 
    @param lev lever arm vector */
    Vector3 Velocity( Vector3 const& lev ) const;
   
    /** @brief Writes particle kinematics in an output stream with a high
    precision format
    @param fileOut output stream */
    void writeParticleKinematics( ostream& fileOut ) const;
  
    /** @brief Writes particle kinematics in an output stream with a high
    precision and 2014 format
    @param fileOut output stream 
    @param particle particle related to the kinematics */
    void writeParticleKinematics2014( ostream& fileOut, 
    	Particle const* particle ) const; 
  
    /** @brief Writes particle kinematics in an output stream with a binary 
    and 2014 format
    @param fileOut output stream 
    @param particle particle related to the kinematics */
    void writeParticleKinematics2014_binary( ostream& fileOut, 
    	Particle const* particle );

    /** @brief Reads particle kinematics from a stream in the 2014 format 
    @param StreamIN input stream 
    @param particle particle related to the kinematics */
    void readParticleKinematics2014( istream& StreamIN,
    	Particle* particle ); 
  
    /** @brief Reads particle kinematics from a stream in a binary form in the
    2014 format 
    @param StreamIN input stream 
    @param particle particle related to the kinematics */
    void readParticleKinematics2014_binary( istream& StreamIN,
    	Particle* particle );       
  
    /** @brief Copies kinematics at time t-2dt (translational velocity, angular 
    velocity, variation of translational velocity, variation of angular 
    velocity) in a 1D array 
    @param vit 1D array where kinematics at time t-2dt is copied
    @param i start index to copy in the 1D array */
    void copyKinematicsNm2( double* vit, int i ) const;
    //@}


    /** @name Accessors */
    //@{
    /** @brief Returns the rotation quaternion */
    Quaternion const* getQuaternionRotation() const;

    /** @brief Returns angular velocity */
    Vector3 const* getAngularVelocity() const;

    /** @brief Returns translational velocity */
    Vector3 const* getTranslationalVelocity() const;
    
    /** @brief Returns the number of bytes of the ParticleKinematics when 
    written in a binary format to an output stream */
    size_t get_numberOfBytes() const;    
    //@}
  

    /** @name Set methods */
    //@{
    /** @brief Resets kinematics to 0 */
    void reset();

    /** @brief Sets the rotation quaternion
    @param vecteur0 x component of the quaternion
    @param vecteur1 y component of the quaternion
    @param vecteur2 z component of the quaternion
    @param scalaire scalar component of the quaternion */
    virtual void setQuaternionRotation( double const& vecteur0, 
	double const& vecteur1, 
	double const& vecteur2, 
	double const& scalaire );
	
    /** @brief Sets the rotation quaternion
    @param qrot rotation quaternion */
    virtual void setQuaternionRotation( Quaternion const& qrot );

    /** @brief Sets the angular velocity
    @param omega angular velocity */
    void setAngularVelocity( Vector3 const& omega );
    
    /** @brief Sets the angular velocity
    @param omx x-angular velocity component 
    @param omy y-angular velocity component     
    @param omz z-angular velocity component */
    void setAngularVelocity( double const& omx, double const& omy,
	double const& omz );     

    /** @brief Sets the translation velocity
    @param vtrans translation velocity */
    void setTranslationalVelocity( Vector3 const& vtrans );
    
    /** @brief Sets the translation velocity
    @param vx x-translation velocity component 
    @param vy y-translation velocity component     
    @param vz z-translation velocity component */
    void setTranslationalVelocity( double const& vx, double const& vy,
	double const& vz );    

    /** @brief Sets kinematics at time t-2dt from a 1D array of 12 scalars
    (translational velocity, angular velocity, variation of translational
    velocity, variation of angular velocity)
    @param tab 1D array of 4 vectors containing translational velocity, angular
    velocity, variation of translational velocity and variation of angular 
    velocity */
    void setKinematicsNm2( double const* tab );
    
    /** @brief Sets time integration scheme using the macro variable
    GrainsExec::m_TIScheme */
    void setTimeIntegrationScheme();
    
    /** @brief Sets momentum equation coupling factor
    @param particle_density particle density */
    void setCouplingFactor( double const& particle_density );         
    //@}

  
    /** @name State storing/restoring methods */
    //@{
    /** @brief Saves particle kinematics state */
    void saveState();
  
    /** @brief Creates and returns particle kinematics state */
    ParticleKinematicsMemento* createState();  
  
    /** @brief Restores particle kinematics state */
    void restoreState();
  
    /** @brief Restores particle kinematics state
    @param memento_ kinematics state */
    void restoreState( ParticleKinematicsMemento const* memento_ );  
    //@}


    /**@name Friend methods */
    //@{
    /** @brief Input operator
    @param fileIn input stream
    @param cinematique ParticleKinematics object */
    friend istream& operator >> ( istream& fileIn, 
    	ParticleKinematics& cinematique );
  
    /** @brief Output operator
    @param fileOut output stream
    @param cinematique ParticleKinematics object */
    friend ostream& operator << ( ostream& fileOut, 
	ParticleKinematics const& cinematique );
    //@}
    

  protected:
    /** @name Constructors */
    //@{
    /** @brief Default constructor forbidden if not called by derived class */
    ParticleKinematics();

    /** @brief Copy constructor forbidden if not called by derived class
    @param copy copied ParticleKinematics object */
    ParticleKinematics( ParticleKinematics const& copy );
    //@}


    /** @name Parameters */
    //@{
    Vector3 m_translationalVelocity; /**< Translational velocity */
    Vector3 m_angularVelocity; /**< Angular velocity */
    Vector3 m_angularVelocity_bf; /**< Angular velocity in body-fixed
    	coordinates system */      
    Quaternion m_QuaternionRotation; /**< Rotation quaternion */ 
    Vector3 m_translationalMotionOverDt; /**< Translational motion
  	over dt */
    Vector3 m_averageAngularVelocityOverDt; /**< Average angular velocity over 
  	dt */
    Quaternion m_QuaternionRotationOverDt; /**< Rotation quaternion over dt */
    TimeIntegrator* m_timeIntegrationScheme; /**< Time integration scheme */
    double m_coupling_factor; /**< Pre-factor of momentum change. 
    	(1) purely granular: FluidCorrectedAcceleration is true and fluid 
	density is 0 so coupling factor = 1. (2) coupled to fluid, fluid density
	is not 0: a. FluidCorrectedAcceleration is true so coupling factor = 
	1 - fluid density / particle density or b. FluidCorrectedAcceleration 
	is false so coupling factor = 1 **/
    ParticleKinematicsMemento* m_memento; /**< Object to save the kinematics 
  	state */ 
    //@}
};

#endif
