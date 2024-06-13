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

    /** @brief Returns the remaining active time interval of the imposed motion
    @param debut simulation start time
    @param fin simulation end time */
    double getTime( double debut, double fin ) const;

    /** @brief Returns whether the imposed motion is activ at time t
    @param t physical time
    @param dt time step magnitude */
    bool isActif( double t, double dt ) const;

    /** @brief Returns whether the imposed motion is completed at time t
    @param t physical time
    @param dt time step magnitude */
    bool isCompleted( double t, double dt ) const;
  
    /** @brief Returns the translational velocity at time t 
    @param time physical time
    @param dt time step magnitude 
    @param cg center of mass of the obstacle */
    Vector3 const* translationalVelocity( double time, double dt, 
    	Point3 const& cg );

    /** @brief Returns the angular velocity at time t 
    @param time physical time
    @param dt time step magnitude */
    Vector3 const* angularVelocity( double time, double dt );
 
    /** @brief Returns the translational motion over dt at time t 
    @param time physical time
    @param dt time step magnitude 
    @param cg center of mass of the obstacle */
    Vector3 translationalMotion( double time, double dt, 
    	Point3 const& cg );  
  
    /** @brief Returns the angular motion over dt at time t 
    @param time physical time
    @param dt time step magnitude */  
    Vector3 angularMotion( double time, double dt ); 

    /** @brief Debug
    @param c debug message */
    void debug( char *c ); 
  
    /** @brief Returns the imposed motion type */
    string getType() const; 
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
    Vector3 m_angularVelocity; /**< angular velocity */
    bool m_rotationCenterIsCenterOfMass; /**< true if the center of rotation
    	is the center of mass of the obstacle. In this case, there is no
    	contribution to the translation motion, otherwise there is */
    Point3 m_rotationCenter; /**< center of rotation */
    double m_Sin_amplitude; /**< sinusoidal velocity amplitude */
    double m_Sin_period; /**< sinusoidal velocity period */
    double m_Sin_phase_shift; /**< sinusoidal velocity phase shift */    
    Vector3 m_unit_vitRef; /**< sinusoidal (cyclic) velocity unit reference 
    	vector */  
    Vector3 m_SinCyclic_period; /**< sinusoidal cyclic motion period in each
    	direction */
    Vector3 m_SinCyclic_amplitude; /**< sinusoidal cyclic motion amplitude in 
    	each direction */
    Vector3 m_SinCyclic_phase_shift; /**< sinusoidal velocity phase shift in 
    	each direction */    
    //@}
    

    /**@name Constructors & Destructor */
    //@{
    /** @brief Copy constructor
    @param copy copied ObstacleImposedVelocity object */
    ObstacleImposedVelocity( ObstacleImposedVelocity const& copy );
    //@}    
};

#endif
