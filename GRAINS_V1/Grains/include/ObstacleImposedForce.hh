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

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation 
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
    string getNom() const;

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
    @param obstacle the obstacle the force is imposed on */
    Vector3 const* translationalVelocity( double time, double dt, 
      Obstacle* obstacle ); 

    /** @brief Returns the imposed force in cyclic mode 
    @param time physical time */
    Vector3 cyclicForce( double time ) const ;

    /** @brief Returns the imposed force */
    Vector3 getForce() const ;

    /** @brief Returns the obstacle virtual mass */
    double getMass() const ;

    /** @brief Returns the direction of displacement */
    Vector3 const* getDirection() const ;

    /** @brief Returns the imposed force type */
    string getType() const; 
    //@}


    /**@name Methods Static */
    //@{
    /** @brief Creates and reads the imposed force features from an input stream
    @param fileIn input stream */
    static ObstacleImposedForce* read( istream& fileIn );
    //@}


  private:
    /**@name Parameters */
    //@{  
    string m_ObstacleName; /**< Obstacle name the force is imposed to */  
    double m_tstart; /**< start time */
    double m_tend; /**< end time */
    Vector3 m_force; /**< Imposed force */
    string m_type; /**< Force type */
    double m_mass; /**< Virtual mass of obstacle */
    Vector3 m_direction; /**< Displacement (or force) direction */
    Vector3 m_translationalVelocity; /**< translational velocity */
    double m_freqX; /**< cyclic motion frequency in x */
    double m_freqY; /**< cyclic motion frequency in y */
    double m_freqZ; /**< cyclic motion frequency in z */
    double m_phase; /**< cyclic motion phase shift */
    Vector3 m_prev; /**< cyclic motion previous position */    
    //@}
    
    
    /**@name Constructors & Destructor */
    //@{
    /** @brief Copy constructor */
    ObstacleImposedForce( ObstacleImposedForce const& copy );
    //@}    
};

#endif
