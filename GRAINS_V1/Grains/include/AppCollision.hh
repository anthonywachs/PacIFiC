#ifndef _APPCOLLISION_HH_
#define _APPCOLLISION_HH_

#include "App.hh"
#include "AllComponents.hh"
#include <iostream>
#include <list>
using namespace std;

class SimpleObstacle;
class Matrix;
class GrainsMPIWrapper;


struct PointForcePostProcessing
{
  Point3 geometricPointOfContact; /**< contact point */
  Vector3 contactForceComp0; /**< contact force exerted on component 0 */
  Point3 PPptComp0; /**< post processing contact point in component 0 */
  Point3 PPptComp1; /**< post processing contact point in component 1 */
}; 


/** @brief The class AppCollision.

    Generic application to detect collisions and compute collision force 
    & torques on rigid bodies.

    @author Institut Francais du Petrole - 2003 - Creation 
    @author A.WACHS - 2021 - Major cleaning & refactoring */
// ============================================================================
class AppCollision : public App
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Destructor */
    virtual ~AppCollision();
    //@}


    /** @name Virtual Methods */
    //@{
    /** @brief Computes forces and torques exerted on rigid bodies
    @param time physical time
    @param dt time step magnitude 
    @param particles active particles */
    virtual void ComputeForces( double time, double dt,
    	list<Particle*> const* particles ) = 0;   	
	
    /** @brief Returns whether a particle is in contact with another component
    using the method Component::isContact
    @param particle particle */
    virtual bool isContact( Particle const* particle ) const;
  
    /** @brief Returns whether a particle is in contact with another component
    using the method Component::isContactWithCrust
    @param particle particle 
    @param BVonly test contact with bounding volume only if true */
    virtual bool isContactWithCrust( Particle const* particle,
    	bool BVonly = false ) const;
  
    /** @brief Returns whether a particle is close to another component
    using the method Component::isClose
    @param particle particle */
    virtual bool isClose( Particle const* particle ) const ; 
  
    /** @brief Returns whether a particle is close to another component
    using the method Component::isCloseWithCrust
    @param particle particle */
    virtual bool isCloseWithCrust( Particle const* particle ) const;
    
    /** @brief Returns whether a point lies inside any particle in the domain
    @param pt point */
    virtual bool isInParticle( Point3 const& pt ) const;    
  
    /** @brief Links a particle with the contact detection algorithm 
    @param particle particle */
    virtual void Link( Particle* particle ) = 0;

    /** @brief Links the root obstacle with the contact detection algorithm at
    the start of the simulation
    @param root_obstacle root obstacle */
    virtual void Link( Obstacle* root_obstacle );

    /** @brief Updates links between particles & obstacles and the contact 
    detection algorithm
    @param time physical time
    @param dt time step magnitude 
    @param particles active particles */
    virtual void LinkUpdate( double time, double dt,
  	list<Particle*>* particles ) = 0;
  
    /** @brief Removes a particle from the contact detection algorithm
    @param particle particle to be removed */
    virtual void remove( Particle* particle ) = 0; 

    /** @brief Removes an obstacle from the contact detection algorithm
    @param obs obstacle to be removed */
    virtual void remove( SimpleObstacle* obs );  	
    //@}


    /** @name Methods */
    //@{  
    /** @brief Computes the average number of particles in the simulation, the
    average is performed over the number of times this method is called 
    @param nbPart number of particles */
    void computeMeanNbParticles( size_t const& nbPart );
  
    /** @brief Adds a contact point to the contact statistics
    @param time physical time 
    @param contactPoint contact point */
    void addToContactsFeatures( double time, PointContact const& contactPoint );
  
    /** @brief Returns a pointer to a simple obstacle based on its ID number
    @param num obstacle number */
    SimpleObstacle* getSimpleObstacle( int const& num ) const;  

    /** @brief Returns maximum overlap */
    double getOverlapMax(); 
  
    /** @brief Returns average overlap */
    double getOverlapMean();   
  
    /** @brief Returns time of maximum overlap */
    double getTimeOverlapMax();
  
    /** @brief Returns average number of iterations of GJK for convergence */
    double getNbIterGJKMean(); 
  
    /** @brief Returns average number of particles */
    double getNbParticlesPerProcMean();      
  
    /** @brief Sets the contact statistics
    @param overlap_max_ maximum overlap
    @param overlap_mean_ average overlap  
    @param time_overlapMax_ time of maximum overlap 
    @param nbIterGJK_ average number of iterations of GJK for convergence */
    void setContactsFeatures( double const& overlap_max_,
	  double const& overlap_mean_,
	  double const& time_overlapMax_,
	  double const& nbIterGJK_ );
	  
    /** @brief Resets postprocessing force index to 0 */
    void resetPPForceIndex();
    
    /** @brief Adds a postprocessing force */
    void addPPForce( Point3 const& pc, Vector3 const& force,
	Component* comp0_, Component* comp1_);
	
    /** @brief Returns the number of postprocessing forces */
    size_t getNbPPForces() const;
    
    /** @brief Returns a pointer to the vector of postprocessing forces */
    vector<struct PointForcePostProcessing> const* getPPForces() const;
    
    /** @brief Resets the list of simple obstacles */
    void resetListSimpleObstacles(); 
    
    /** @brief Computes the macroscopic stress tensor in the whole domain 
    @param wrapper MPI wrapper */
    void computeStressTensor( GrainsMPIWrapper const* wrapper );
	
    /** @brief Sets the parameters to output force statistics
    @param root_ output directory name
    @param freq_ output frequency */
    void setForceStatsParameters( string const& root_,
  	size_t const& freq_ );
	
    /** @brief Returns whether to output force statistics at this time 
    @param enforceOutput force writing   
    @param increaseCounterOnly increases the writing counter only */
    bool outputForceStatsAtThisTime( bool enforceOutput, 
    	bool increaseCounterOnly );
	
    /** @brief Writes force statistics in a file
    @param time physical time
    @param dt time step magnitude
    @param rank process rank 
    @param wrapper MPI wrapper */
    void outputForceStats( double time, double dt, int rank,
    	GrainsMPIWrapper const* wrapper );
	
    /** @brief Initialises output files to write force statistics
    @param rank process rank
    @param coupledFluid whether the simulation is coupled to a fluid solver
    @param time physical time */
    void initialiseForceStatsFiles( int rank,
      	bool coupledFluid, double time );
	
    /** @brief Returns the macroscopic stress tensor in the whole domain */
    Matrix const* getStressTensor() const;
    
    /** @brief Returns a component of the macroscopic stress tensor in the 
    whole domain */
    double getStressTensorComponent( int k, int l ) const;    		
    //@}  

  
  protected:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    AppCollision();
    //@}


    /** @name Parameters */
    //@{  
    Obstacle* m_obstacles; /**< Root obstacle */  
    list<SimpleObstacle*> m_allSimpleObstacles; /**< List of simple obstacles */
    double m_overlap_max; /**< Maximum overlap between 2 colliding
    	components */
    double m_overlap_mean; /**< Average overlap between 2 colliding
    	components */ 
    double m_time_overlapMax; /**< Time of maximum overlap */
    double m_nbIterGJK_mean; /**< Average number of iterations of GJK for
    	convergence */
    double m_nbParticles_mean; /**< Average number of particles */
    vector<struct PointForcePostProcessing> m_allforces;
    size_t m_allforces_index;
    bool m_outputForceStats; /**< whether to to write force statistics */
    string m_outputForceStats_dir; /**< directory name where to write force
  	statistics files, including average macro stress in the whole domain */
    size_t m_outputForceStats_counter; /**< counter for force statistics 
    	output */
    size_t m_outputForceStats_frequency; /**< frequency of force statistics 
    	output */
    static size_t m_allforces_blocksize;
    Matrix m_stressTensor; /**< Macroscopic stress tensor in the whole domain */
  //@}
};

#endif
