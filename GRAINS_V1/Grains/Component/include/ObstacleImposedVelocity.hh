#ifndef _OBSTACLEIMPOSEDVELOCITY_HH_
#define _OBSTACLEIMPOSEDVELOCITY_HH_

#include "Vector3.hh"
#include "Point3.hh"
using namespace solid;
#include <list>
#include <string>
#include <iostream>
using namespace std;
#include "ReaderXML.hh"


class ObstacleImposedVelocity;
class LinkedCell;
bool operator < ( ObstacleImposedVelocity const& c0,
	ObstacleImposedVelocity const& c1 );
ostream& operator << ( ostream& fileOut, 
	ObstacleImposedVelocity const& motion );
istream& operator >> ( istream& fileIn, ObstacleImposedVelocity& motion );


/** @brief The class ObstacleImposedVelocity.

    Defines and controls the velocity imposed on an obstacle.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ObstacleImposedVelocity
{
  public:
    /**@name Constructors & Destructor */
    //@{
    /** @brief Default constructor */
    ObstacleImposedVelocity();

    /** @brief Constructor with an XML node as input parameter
    @param root XML node
    @param dt time step magnitude
    @param rank MPI rank 
    @param error error in reading the XML node */
    ObstacleImposedVelocity( DOMNode* root, double dt, int rank, 
    	size_t& error );

    /** @brief Destructor */
    ~ObstacleImposedVelocity();
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
  
    /** @brief Returns the translational velocity at time time 
    @param time physical time
    @param dt time step magnitude 
    @param cg center of mass of the obstacle */
    Vector3 const* translationalVelocity( double time, double dt, 
    	Point3 const& cg );

    /** @brief Returns the angular velocity at time time 
    @param time physical time
    @param dt time step magnitude */
    Vector3 const* angularVelocity( double time, double dt );
 
    /** @brief Returns the translational motion over [time-dt,time] 
    @param time physical time
    @param dt time step magnitude 
    @param subinterval the sub-interval length within [time-dt,time] when it
    is actually active
    @param cg center of mass of the obstacle */
    Vector3 translationalMotion( double time, double dt, 
    	double const& subinterval, Point3 const& cg );  
  
    /** @brief Returns the angular motion over [time-dt,time]
    @param time physical time
    @param dt time step magnitude 
    @param subinterval the sub-interval length within [time-dt,time] when it
    is actually active */  
    Vector3 angularMotion( double time, double dt, double const& subinterval ); 

    /** @brief Debug
    @param c debug message */
    void debug( char *c ); 
  
    /** @brief Returns the imposed motion type */
    string getType() const;
    
    /** @brief Updates imposed velocity based on a stress criterion (for 
    cyclic shearing) 
    @param LC linked cell grid */
    void updateImposedVelocity( LinkedCell const* LC );     
    //@}


    /**@name Operators */
    //@{
    /** @brief Operator == based on the object address
    @param other the other ObstacleImposedVelocity */
    bool operator == ( ObstacleImposedVelocity const& other ) const;
    //@}


    /**@name Methods Friends */
    //@{
    /** @brief Operator < based on the start time of the imposed motion.
    Returns true if c0.tdebut < c1.tdebut */
    friend bool operator < ( ObstacleImposedVelocity const& c0,
	ObstacleImposedVelocity const& c1 );
    //@}


  private:
    /**@name Parameters */
    //@{
    string m_ObstacleName; /**< Obstacle name the motion is imposed to */
    string m_type; /**< Motion type */
    double m_tstart; /**< Start time */
    double m_tend; /**< End time */
    Vector3 m_translationalVelocity; /**< translational velocity */
    Vector3 m_previous_translationalVelocity; /**< translational velocity at
    	time - dt */    
    Vector3 m_angularVelocity; /**< angular velocity */
    bool m_rotationCenterIsCenterOfMass; /**< true if the center of rotation
    	is the center of mass of the obstacle. In this case, there is no
    	contribution to the translation motion, otherwise there is */
    Point3 m_rotationCenter; /**< center of rotation */
    double m_amplitude; /**< sinusoidal or step velocity amplitude */
    double m_period; /**< sinusoidal or step velocity period */
    double m_Sin_phase_shift; /**< sinusoidal velocity phase shift */    
    Vector3 m_unit_vitRef; /**< multi-dimensional sinusoidal velocity unit 
    	reference vector */  
    Vector3 m_MultiSin_period; /**< multi-dimensional sinusoidal motion period 
    	in each direction */
    Vector3 m_MultiSin_amplitude; /**< multi-dimensional sinusoidal motion 
    	amplitude in each direction */
    Vector3 m_MultiSin_phase_shift; /**< multi-dimensional sinusoidal velocity 
    	phase shift in each direction */ 
    pair<int,int> m_stressIndices; /**< Stress component indices */
    double m_stress_max; /**< Maximum stress amplitude to revert the motion */
    double m_stress; /**< Current stress */
    double m_previous_stress; /**< Previous stress */ 
    //@}
    

    /**@name Constructors & Destructor */
    //@{
    /** @brief Copy constructor
    @param copy copied ObstacleImposedVelocity object */
    ObstacleImposedVelocity( ObstacleImposedVelocity const& copy );
    //@}    
};

#endif
