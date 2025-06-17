#ifndef _OBSTACLEIMPOSEDFORCE_HH_
#define _OBSTACLEIMPOSEDFORCE_HH_

#include "Vector3.hh"
#include "ObstacleImposedVelocity.hh"
using namespace solid;
#include <list>
#include <string>
#include <vector>
#include <iostream>
using namespace std;


class Obstacle;


/** @brief The class ObstacleImposedForce.

    Defines and controls the force imposed on an obstacle.

    @author Institut Francais du Petrole - 2003 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ObstacleImposedForce
{
  public:
    /**@name Constructors & Destructor */
    //@{
    /** @brief Default constructor */
    ObstacleImposedForce();

    /** @brief Constructor with an XML node as input parameter
    @param root XML node
    @param dt time step magnitude
    @param rank MPI rank 
    @param error error in reading the XML node */
    ObstacleImposedForce( DOMNode* root, double dt, int rank, size_t& error );

    /** @brief Destructor */
    ~ObstacleImposedForce();
    //@}


    /**@name Methods */
    //@{
    /** @brief Returns obstacle name */
    string getObstacleName() const;

    /** @brief Returns whether the imposed motion is active over the time
    interval [td,te] and if it is the sub-interval length within [td,te] when it
    is actually active 
    @param td simulation start time
    @param te simulation end time
    @param dt time step magnitude 
    @param subinterval the sub-interval length within [td,te] when it
    is actually active */
    bool isActif( double const& td, double const& te, double const& dt,
    	double& subinterval ) const;

    /** @brief Returns whether the imposed motion is completed at time t
    @param t physical time
    @param dt time step magnitude */
    bool isCompleted( double t, double dt ) const;

    /** @brief Computes the translational velocity at time time and the
    translation over [time-dt,time]
    @param time physical time
    @param dt time step magnitude 
    @param subinterval the sub-interval length within [td,te] when it
    is actually active
    @param vt translational velocity at time time
    @param translation translation over [time-dt,time]
    @param obstacle the obstacle the force is imposed on */
    void translationalVelocity( double time, double dt, 
      double const& subinterval, Vector3& vt, Vector3& translation, 
      const Obstacle* obstacle ); 

    /** @brief Returns the imposed force 
    @param time physical time */
    Vector3 getForce( double time );

    /** @brief Returns the obstacle virtual mass */
    double getMass() const ;

    /** @brief Returns the direction of motion */
    Vector3 const* getDirection() const ;

    /** @brief Returns the imposed force type */
    string getType() const; 
    //@}


  private:
    /**@name Parameters */
    //@{  
    string m_ObstacleName; /**< Obstacle name the force is imposed to */  
    double m_tstart; /**< start time */
    double m_tend; /**< end time */
    Vector3 m_force_amplitude; /**< Imposed force amplitude */
    Vector3 m_force; /**< Imposed force at time t */    
    string m_type; /**< Force type */
    double m_mass; /**< Virtual mass of obstacle */
    bool m_automass; /**< true if the mass is determined by the code */
    Vector3 m_direction; /**< Motion (or force) direction, i.e. unit imposed
    	force vector */
    Vector3 m_translationalVelocity; /**< translational velocity */
    Vector3 m_MultiSin_period; /**< multi-dimensional sinusoidal force period 
    	in each direction */
    Vector3 m_MultiSin_phase_shift; /**< multi-dimensional sinusoidal velocity 
    	phase shift in each direction */    
    double m_vmaxzeroforce; /**< maximum velocity magnitude at early stages if 
    	the force is zero */        
    //@}
    
    
    /**@name Constructors & Destructor */
    //@{
    /** @brief Copy constructor */
    ObstacleImposedForce( ObstacleImposedForce const& copy );
    //@}
    
    
    /**@name Methods */
    //@{
    /** @brief Sets the sinusoidal cyclic force at a given time 
    @param time physical time */
    void MultiSinForce( double time );
    //@}        
};

#endif
