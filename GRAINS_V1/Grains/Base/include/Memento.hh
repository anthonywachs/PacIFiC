#ifndef _MEMENTO_HH_
#define _MEMENTO_HH_


#include "Quaternion.hh"
#include "Transform.hh"
#include "Vector3.hh"
#include "Particle.hh"
#include "ParticleKinematics.hh"
#include "ParticleKinematics3D.hh"
#include "ParticleKinematics2D.hh"


/** @brief The class memento.

    Memorizes the geometric configuration of a component. 

    @author GRAINS Project - IFP - 2008 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ConfigurationMemento
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Destructor */
    virtual ~ConfigurationMemento() {};
    //@}


  private:
    /** @name Friend classes */
    //@{
    friend class Particle;
    friend class Obstacle;
    friend class CompositeObstacle;  
    friend class SimpleObstacle;
    friend class Component; 
    //@}     
  
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    ConfigurationMemento() {};
    //@}


    /** @name Parameters */
    //@{
    Transform m_position; /**< rigid body's transformation: center of mass 
  	position and body orientation */
    //@}
};




/** @brief The class ParticleKinematicsMemento.

    Memorizes a particle kinematics. 

    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ParticleKinematicsMemento
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Destructor */
    virtual ~ParticleKinematicsMemento() {};
    //@}
  

  private:
    /** @name Friend classes */
    //@{
    friend class ParticleKinematics;
    friend class ParticleKinematics3D;
    friend class ParticleKinematics2D;
    //@}
  
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    ParticleKinematicsMemento() {};
    //@}


    /** @brief Parameters */
    //@{
    Quaternion m_QuaternionRotation; /**< Rotation quaternion */
    solid::Vector3 m_translationalVelocity; /**< Translational velocity */
    Vector3 m_angularVelocity; /**< Angular velocity */
    //@}
};




/** @brief The class ObstacleKinematicsMemento.

    Memorizes an obstacle kinematics. 

    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ObstacleKinematicsMemento
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Destructor */
    virtual ~ObstacleKinematicsMemento() {};
    //@}
  

  private:
    /** @name Friend classes */
    //@{
    friend class ObstacleKinematicsVelocity;
    //@}
  
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    ObstacleKinematicsMemento() {};
    //@}


    /** @brief Parameters */
    //@{
    solid::Vector3 m_translationalVelocity; /**< Translational velocity */
    solid::Vector3 m_angularVelocity; /**< Angular velocity */ 
    //@}
};


#endif
