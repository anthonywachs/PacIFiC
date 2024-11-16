#ifndef _CONTACTFORCEMODEL_HH_
#define _CONTACTFORCEMODEL_HH_

#include "Basic.hh"
#include "Point3.hh"
#include "PointContact.hh"
#include "Vector3.hh"
#include "ReaderXML.hh"
#include "LinkedCell.hh"
#include <list>
#include <map>
using namespace std;

class Component;
class LinkedCell;


/** @brief The class ContactForceModel.

    Defines the contact forces between two colliding components and computes
    these contact forces.

    @author A.WACHS - IFPEN - 2011 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ContactForceModel
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with a map of contact parameters as inputs
    @param parameters map of parameters */
    ContactForceModel( map<string,double>& parameters );

    /** @brief Destructor */
    virtual ~ContactForceModel();
    //@}


    /**@name Methods */
    //@{  
    /** @brief Returns the name of the contact force model */
    virtual string name() const = 0;   
  
    /** @brief Computes an estimate of the contact time and maximum penetration 
    depth in the case of a gravityless binary collision of spheres, and writes
    the result in an output stream
    @param p0_ first component (Particle)
    @param p1_ second component (Particle or Obstacle)
    @param v0 pre-collisional relative velocity 
    @param dt time step to integrate equations in case an analytical solution is
    not known 
    @param OUT output stream */
    virtual void computeAndWriteEstimates( Component* p0_,  
	Component* p1_,
  	double const& v0, double const& dt,
	ostream& OUT ) const = 0 ;

    /** @brief Computes forces & torques 
    @param p0_ first Component (Particle)
    @param p1_ second Component (Particle ou Obstacle)
    @param contactInfos geometric contact features
    @param LC linked-cell grid
    @param dt time step magnitude
    @param nbContact number of contact points for composite particles */
    virtual bool computeForces( Component* p0_, Component* p1_,
	PointContact const& contactInfos, LinkedCell* LC,
	double dt, int nbContact = 1 ) = 0 ;
    //@}     

  protected:
    /**@name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    ContactForceModel() {};
    //@}  
};

#endif
