#ifndef _OBSTACLEKINEMATICSVELOCITY_HH_
#define _OBSTACLEKINEMATICSVELOCITY_HH_

#include "Kinematics.hh"
#include "ObstacleImposedVelocity.hh"
#include "Quaternion.hh"
#include <stdlib.h>

class SimpleObstacle;
class Torsor;
class ObstacleKinematicsMemento;


class ObstacleKinematicsVelocity;
ostream& operator << ( ostream& fileOut, ObstacleKinematicsVelocity const& kine_ );
istream& operator >> ( istream& fileIn, ObstacleKinematicsVelocity& kine_ );


/** @brief The class ObstacleKinematicsVelocity.

    Manages the obstacle kinematics with a prescribed velocity.

    @author Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ObstacleKinematicsVelocity : public Kinematics
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    ObstacleKinematicsVelocity();

    /** @brief Destructor */
    virtual ~ObstacleKinematicsVelocity();
    //@}
    

    /**@name Methods */
    //@{
    /** @brief Adds an imposed velocity motion to the obstacle kinematics
    @param oiv the imposed velocity motion */
    void append( ObstacleImposedVelocity* oiv );

    /** @brief Composes the obstacle kinematics with another "higher level"
    velocity kinematics
    @param other the higher level kinematics
    @param lever lever arm of the higher level kinematics applied to the
    obstacle */
    void Compose( ObstacleKinematicsVelocity const& other, 
    	Vector3 const& lever );

    /** @brief Updates the obstacle translational and angular velocity at time 
    time and translational and angular motion from time - dt to time and 
    returns whether the obstacle moved from time - dt to time
    @param time physical time
    @param dt time step magnitude 
    @param cg center of mass of the obstacle */
    bool ImposedMotion( double time, double dt, Point3 const& cg );
  
    /** @brief Computes the total velocity of the obstacle using the arm lever
    @param om arm lever */
    Vector3 Velocity( Vector3 const& om ) const; 
  
    /** @brief Returns whether there is an active angular motion imposed from
    time - dt to time
    @param time physical time
    @param dt time step magnitude */
    bool activeAngularMotion( double time, double dt ) const; 
    //@}
  

    /** @name Accessors */
    //@{
    /** @brief Returns a pointer to the rotation quaternion over dt */
    Quaternion const* getQuaternionRotationOverDt() const;

    /** @brief Returns a pointer to the translation vector over dt */
    Vector3 const* getTranslation() const;

    /** @brief Returns a pointer to the current angular velocity vector */  
    Vector3 const* getAngularVelocity() const;

    /** @brief Returns a pointer to the current translational velocity vector */
    Vector3 const* getTranslationalVelocity() const;

    /** @brief Returns the list of velocity motions imposed on the obstacle */
    list<ObstacleImposedVelocity*> getChargements() const;
    //@}
  

    /** @name Methods Set */
    //@{
    /** @brief Resets kinematics to 0 */
    void reset();

    /** @brief Sets the velocity and motion using another velocity 
    kinematics 
    @param kine_ the other velocity kinematics */
    void set( ObstacleKinematicsVelocity& kine_ );
  
    /** @brief Sets velocity
    @param vtrans translational velocity 
    @param vrot angular velocity */
    void setVelocity( Vector3 const* vtrans, Vector3 const* vrot );   
    //@}


    /** @name State storing/restoring methods */
    //@{
    /** @brief Saves obstacle kinematics state */
    void saveState();
  
    /** @brief Creates and returns obstacle kinematics state */
    ObstacleKinematicsMemento* createState();  
  
    /** @brief Restores obstacle kinematics state */
    void restoreState();
  
    /** @brief Restores obstacle kinematics state
    @param memento_ kinematics state */
    void restoreState( ObstacleKinematicsMemento const* memento_ );  
    //@}


  private:
    /**@name Methods */
    //@{
    /** @brief Deletes all imposed motions */
    void clearAndDestroy();
    //@}


    /**@name Parameters */
    //@{
    Vector3 m_translationOverTimeStep; /**< Translational motion over 
    	dt */
    Vector3 m_rotationOverTimeStep; /**< Rotation over dt  */
    Quaternion m_QuaternionRotationOverDt; /**< Quaternion describing the 
    	rotation over dt */
    list<ObstacleImposedVelocity*> m_imposedVelocities; /**< list of imposed 
    	velocity motions */  
    Vector3 m_translationalVelocity; /**< obstacle translational velocity */
    Vector3 m_angularVelocity; /**< obstacle angular velocity */
    ObstacleKinematicsMemento *m_memento; /**< to store the kinematics 
  	features */    
    //@}
};

#endif
  
