#ifndef _HODCCONTACTFORCEMODEL_HH_
#define _HODCCONTACTFORCEMODEL_HH_

#include "ContactForceModel.hh"
#include "Basic.hh"
#include "Point3.hh"
#include "PointContact.hh"
#include "Vector3.hh"
#include "ReaderXML.hh"
#include <map>
using namespace std;

class Component;


/** The class HODCContactForceModel.

    Contact force model involving a normal Hookean spring, a normal Dashpot and
    a tangential Coulomb friction (HO-D-C) to compute the force and torque 
    induced by the contact between two rigid components.

    @author A.WACHS - IFPEN - 2011 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class HODCContactForceModel : public ContactForceModel
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with a map of contact parameters as inputs
    @param parameters map of parameters */
    HODCContactForceModel( map<string,double>& parameters );

    /** @brief Destructor */
    virtual ~HODCContactForceModel();
    //@}


    /** @name Methods */
    //@{  
    /** @brief Returns the name of the contact force model */
    string name() const;
  
    /** @brief Computes an estimate of the contact time and maximum penetration 
    depth in the case of a gravityless binary collision of spheres, and writes
    the result in an output stream
    @param p0_ first component (Particle)
    @param p1_ second component (Particle or Obstacle)
    @param v0 pre-collisional relative velocity 
    @param OUT output stream */
    void computeAndWriteEstimates( Component* p0_,  
	Component* p1_,
  	double const& v0, 
	ostream& OUT ) const ;
  
    /** @brief Computes forces & torques 
    @param p0_ first Component (Particle)
    @param p1_ second Component (Particle ou Obstacle)
    @param contactInfos geometric contact features
    @param LC linked-cell grid
    @param dt time step magnitude
    @param nbContact number of contact points for composite particles */
    bool computeForces( Component* p0_, Component* p1_,
	PointContact const& contactInfos, LinkedCell* LC,
	double dt, int nbContact = 1 ) ;	     		  
    //@}

  
    /** @name Static methods */
    //@{
    /** @brief Reads and returns contact parameter map from an XML node
    @param root XML node */
    static map<string,double> defineParameters( DOMNode* & root );
    //@}


  protected:
    /** @name Parameters */
    //@{
    double stiff; /**< Coefficient d'elasticite */  
    double en; /**< Coefficient de restitution */ 
    double muet; /**< Coefficient de frottement direction tangente */
    double muec; /**< Coefficient de frottement de coulomb */
    double k_m_s; /**< Coefficient */
    //@}


    /**@name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    HODCContactForceModel();
    //@}

  
    /**@name Methods */
    //@{  
    /** @brief Computes maximum penetration depth using a analytical solution
    and a Newton algorithm
    @param theta_ sqrt( omega0*omega0 - mu*mu ) 
    @param mu_ dissipation coefficient
    @param en_ restitution coefficient
    @param tc_ contact time   
    @param v0_ pre-collisional relative velocity */
    double computeDeltaMax( double const& theta_, double const& mu_,
  	double const& en_, double const& tc_, double const& v0_ ) const;
	
    /** @brief Performs forces & torques computation
    @param p0_ first Component (Particle)
    @param p1_ second Component (Particle ou Obstacle)
    @param contactInfos geometric contact features
    @param delFN normal force
    @param delFT tangential force
    @param delM torque */
    void performForcesCalculus( Component* p0_,  Component* p1_,
	PointContact const& contactInfos,
	Vector3& delFN, Vector3& delFT, Vector3& delM );		  
    //@}   
};

#endif
