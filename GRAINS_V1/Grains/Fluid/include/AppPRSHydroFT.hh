#ifndef _APPPRSHYDROFT_HH_
#define _APPPRSHYDROFT_HH_

#include "App.hh"
#include "Basic.hh"
#include "Error.hh"
#include <list>
using namespace std;


/** @brief The class AppPRSHydroFT.

    Used to add hydrodynamic forces and torques exerted on rigid bodies and 
    computed by a Particle-Resolved external solver.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation 
    @author A.WACHS - 2023 - Creation */
// ============================================================================
class AppPRSHydroFT : public App
{
  public:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    AppPRSHydroFT();

    /** @brief Destructor */
    virtual ~AppPRSHydroFT();
    //@}


    /**@name Virtual methods */
    //@{
    /** @brief Computes forces and torques exerted on rigid bodies
    @param time physical time
    @param dt time step magnitude 
    @param particles active particles */
    virtual void ComputeForces( double time, double dt,
    	list<Particle*> const* particles );  
    //@}  


    /** @name Methods */
    //@{  
    /** @brief Allocates the arrays of hydro force and torque 
    @param nbPart number of particles */
    void allocateHydroFT( size_t const& nbPart );
  
    /** @brief Sets the arrays of hydro force and torque 
    @param hydrovec arrays of hydro force and torque */
    void setHydroFT( vector< vector<double> > const* hydrovec );
    //@} 

  
  private:
    /**@name Parameters */
    //@{  
    vector<Vector3> m_PRSHydroForce; /**< arrays of hydro forces */
    vector<Vector3> m_PRSHydroTorque; /**< arrays of hydro torques */    
    //@}   
};

#endif

     
