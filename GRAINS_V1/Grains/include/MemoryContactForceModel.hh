#ifndef _MEMORYCONTACTFORCEMODEL_HH_
#define _MEMORYCONTACTFORCEMODEL_HH_

#include "ContactForceModel.hh"
#include "Basic.hh"
#include "Point3.hh"
#include "PointContact.hh"
#include "Vector3.hh"
#include "ReaderXML.hh"
#include <map>
using namespace std;

class Component;


/** The class MemoryContactForceModel.

    Contact force model involving a Hookean spring and a Dashpot in both the
    normal and tangential directions as well as in the rotational space; and
    a tangential Coulomb friction to compute the force and torque
    induced by the contact between two rigid components.

    @author D.Huet - 2022
// ============================================================================*/
class MemoryContactForceModel : public ContactForceModel
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with a map of contact parameters as inputs
    @param parameters map of parameters */
    MemoryContactForceModel( map<string,double>& parameters );

    /** @brief Destructor */
    virtual ~MemoryContactForceModel();
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

    /** @brief Computation of vector tangent to contact
    @param tij The vector where the result is stored
    @param n_t Eta_t, the tangential damping coefficient
    @param u_t The tangential velocity
    @param kdelta The cumulative displacement */
    void computeTangentialVector(Vector3& tij, double n_t, const Vector3 ut,
  const Vector3 kdelta);


    /** @name Static methods */
    //@{
    /** @brief Reads and returns contact parameter map from an XML node
    @param root XML node */
    static map<string,double> defineParameters( DOMNode* & root );
    //@}


  protected:
    /** @name Parameters */
    //@{
    double stiff; /**< Normal stiffness coefficient */
    double en; /**< Normal restitution coefficient */
    double muet; /**< Tangential damping coefficient */
    double ks; /**< Tangential stiffness coefficient */
    double muec; /**< Tangential Coulomb friction coefficient */
    double eta_r; /**< Viscous rolling damping ratio */
    double mu_r; /**< Rolling resistance coefficient */
    double Jn; /**< dimensionless coefficient used to compute the rotational
                    stiffness coefficient (Ai et al. 2011) */
    double m_f; /**< binary value: if f=0, no viscous rolling damping when the
                   spring rolling friction is saturated. Default value is 0. */
    double epsilon; /**< Criterion on velocity norm for cumulative tangential vector. Default value is 1.e-10. */
    bool rolling_friction; /**< Boolean to switch on/off the rolling friction */
    //@}


    /**@name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    MemoryContactForceModel();
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
    @param delM torque
    @param elementary_id0 ID of elementary particle in case p0_ is a composite
    particle
    @param elementary_id1 ID of elementary particle in case p1_ is a composite
    particle*/
    void performForcesCalculus( Component* p0_,  Component* p1_, double dt,
	PointContact const& contactInfos,
	Vector3& delFN, Vector3& delFT, Vector3& delM, int elementary_id0,
  int elementary_id1 );
    //@}
};

#endif
