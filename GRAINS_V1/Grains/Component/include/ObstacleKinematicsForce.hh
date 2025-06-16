#ifndef _OBSTACLEKINEMATICSFORCE_HH_
#define _OBSTACLEKINEMATICSFORCE_HH_

#include "Kinematics.hh"
#include "ObstacleImposedForce.hh"
#include "Point3.hh"
#include "Vector3.hh"
using namespace solid;


class Obstacle;


/** @brief The class ObstacleKinematicsForce.

    Manages the obstacle kinematics with an imposed force.

    @author G.FERRER - Aout.2003 - Institut Francais du Petrole - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ObstacleKinematicsForce : public Kinematics
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    ObstacleKinematicsForce();

    /** @brief Destructor */
    ~ObstacleKinematicsForce();
    //@}


    /**@name Methods */
    //@{
    /** @brief Adds an imposed force load to the obstacle kinematics
    @param oif the imposed force load */
    void append( ObstacleImposedForce* oif );

    /** @brief Composes the obstacle kinematics with another "higher level"
    force kinematics
    @param other the higher level kinematics */
    void Compose( ObstacleKinematicsForce const& other );
	
    /** @brief Computes the obstacle velocity and returns whether the obstacle 
    moved from time - dt to time
    @param time physical time
    @param dt time step magnitude 
    @param obstacle obstacle the kinematics is related to */
    bool ImposedMotion( double time, double dt, Obstacle* obstacle );	

    /** @brief Returns translational motion over dt 
    @param dt time step magnitude */
    Vector3 getTranslation( double dt ) const;

    /** @brief Resets kinematics to 0 */
    void reset();

    /** @brief Computes the total velocity of the obstacle using the arm lever
    @param om arm lever */
    Vector3 Velocity( const Vector3 &om ) const; 
    
    /** @brief Returns a pointer to the current translational velocity vector */
    Vector3 const* getTranslationalVelocity() const;    
    //@}


  private:
    /**@name Methods */
    //@{
    /** @brief Deletes all imposed force loads */
    void clearAndDestroy();
    //@}
    
      
    /**@name Parameters */
    //@{
    Vector3 m_translationalVelocity; /**< obstacle translational velocity */
    Vector3 m_angularVelocity; /**< obstacle angular velocity */
    double m_vitesseD; /**< norm of obstacle translational velocity */  
    list<ObstacleImposedForce*> m_imposedForces; /**< list of imposed force 
    	loads */
    Vector3 m_translationOverTimeStep; /**< Translational motion over 
    	dt */
    //@}
};

#endif
