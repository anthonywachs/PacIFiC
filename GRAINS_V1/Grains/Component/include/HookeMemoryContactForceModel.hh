#ifndef _HOOKEMEMORYCONTACTFORCEMODEL_HH_
#define _HOOKEMEMORYCONTACTFORCEMODEL_HH_

#include "ContactForceModel.hh"
#include "Basic.hh"
#include "Point3.hh"
#include "PointContact.hh"
#include "Vector3.hh"
#include "ReaderXML.hh"
#include <map>
using namespace std;

class Component;


/** The class HookeMemoryContactForceModel.

    Contact force model involving a Hookean spring and a Dashpot in both the
    normal and tangential directions as well as in the rotational space; and
    a tangential Coulomb friction to compute the force and torque
    induced by the contact between two rigid components.
    
    The force is formulated as follows:
    
    \f{eqnarray*}{
       \mbox{\LARGE $\bm{F}_n$} &=& \mbox{\LARGE $k_n \delta_n \bm{n}
       	- \gamma_n \bm{U}_n$} \\ 
        \mbox{\LARGE $\bm{F}_t$} &=& \mbox{\LARGE $min\left( 
		| -k_t \bm{\delta}_t - \gamma_t \bm{U}_t |, 
		\mu_c | \bm{F}_n | \right)\cdot \bm{t}$}
    \f}
    
    where \f$\delta_n\f$ is the overlap distance in the normal direction, 
    \f$\bm{n}\f$ is the unit normal vector at the contact point, 
    \f$\bm{\delta}_t\f$ is the cumulative tangential displacement at the contact
     point, \f$\bm{U}_n\f$ is the normal relative velocity at the contact point,
     \f$\bm{U}_t\f$ is the tangential relative velocity at the contact point and
    \f$\bm{t}\f$ is the unit tangential vector at the contact point defined as
     \f$\displaystyle -\frac{k_t \bm{\delta}_t + \gamma_t \bm{U}_t}
    {|k_t \bm{\delta}_t + \gamma_t \bm{U}_t|}\f$.
    
    The meaning and expression of the contact force parameters are:
    - \f$k_n\f$ is the normal stiffness coefficient.
    - \f$\gamma_n\f$ is the normal dissipative coefficient defined as
    
    \f[
    \mbox{\LARGE $\displaystyle \gamma_n = - 2\beta\sqrt{k_n m^*} 
    	\:\:\mathrm{with}\:\: 
    	\beta=\frac{\ln(e_n)}{\sqrt{\pi^2+(\ln(e_n))^2}}$}
    \f]
    
    where \f$e_n\f$ is the restitution coefficient, 
    \f$\displaystyle m^*=\frac{1}{\frac{1}{m_0}+\frac{1}{m_1}}\f$ is the reduced
    mass, \f$m_0\f$ is the mass of the first rigid component and \f$m_1\f$ is 
    the mass of the second rigid component. If one of the two components is an
    obstacle, then \f$m_{0\:\mathrm{or}\:1} \rightarrow \infty\f$ such that
    \f$m^*=m_{0\:\mathrm{or}\:1}\f$, and in
    practice we take \f$m_{0\:\mathrm{or}\:1} = 10^{20}\f$.
    
    - \f$\gamma_t\f$ is the tangential dissipative coefficient defined as
    
    \f[
    \mbox{\LARGE $\gamma_t = 2 m^* \eta_t$}
    \f]
    
    where \f$\eta_t\f$ is the input tangential dissipative coefficient. If the
    value specified in the input file is \f$-1\f$, \f$\eta_t\f$ is automatically
    computed such that \f$\gamma_n = \gamma_t\f$, i.e., same damping in the 
    normal and tangential directions. This gives the following expression of 
    \f$\eta_t\f$:
    
    \f[
    \mbox{\LARGE $\eta_t = - 2\beta\sqrt{\frac{k_n}{m^*}}$}
    \f]
    
    - \f$k_t\f$ is the tangential stiffness coefficient.    

    - \f$\mu_c\f$ is the Coulomb tangential friction coefficient. 
    
    The cumulative tangential displacement at the contact point \f$\delta_t\f$
    is defined as follows:
    
    \f[
    \mbox{\LARGE $\displaystyle \bm{\delta}_t = \left\{
    \begin{array}{lcl}
    \bm{q}_{rot} \bm{\delta}_t^{t-\Delta t} \bm{q}_{rot}^{-1} 
    + \int_{t-\Delta t}^t \bm{U}_t(s)ds & \mathrm{if} 
    & |\bm{F}_t| \leq \mu_c | \bm{F}_n | \\
    \frac{-\mu_c| \bm{F}_n |\bm{t} - \gamma_t \bm{U}_t}{k_t} & \mathrm{if} 
    & |\bm{F}_t| > \mu_c | \bm{F}_n |
    \end{array} \right.
    $}\f] 
    
    where \f$\bm{q}_{rot}\f$ is the rotation quaternion from the tangential
    plane at time \f$t-\Delta t\f$ to the tangential plane at time \f$t\f$.

    The rolling resistance torque is formulated as follows:
    
    \f[
    \mbox{\LARGE $\displaystyle \bm{T}_r = \bm{T}_k + \bm{T}_d $}
    \f] 
    
    with    
            
    \f[
    \mbox{\LARGE $\displaystyle \bm{T}_k^t = \left\{
    \begin{array}{lcl}
    \bm{q}_{rot} \bm{T}_k^{t-\Delta t} \bm{q}_{rot}^{-1} 
    - k_r \Delta \bm{\theta} & \mathrm{if} & |\bm{T}_k^t| \leq T_{max} \\
    T_{max}\bm{r} & \mathrm{if} & |\bm{T}_k^t| > T_{max}
    \end{array} \right.
    $}\f] 
    
    where \f$\displaystyle \Delta \bm{\theta} = ( \bm{\omega}_0 -
    \bm{\omega}_1 ) \Delta t\f$ is the rotational increment, 
    \f$\displaystyle T_{max} = \mu_r R^* | \bm{F}_n |\f$ is the saturation
    torque,  and \f$\displaystyle R^*=\frac{1}{\frac{1}{R_0}+\frac{1}{R_1}}\f$ 
    is the reduced radius, \f$R_0\f$ is the radius of the first rigid component 
    and \f$R_1\f$ is the radius of the second rigid component, 
    \f$\displaystyle \bm{r} = \frac{\bm{q}_{rot} \bm{T}_k^{t-\Delta t} 
    \bm{q}_{rot}^{-1} - k_r \Delta \bm{\theta}}{|\bm{q}_{rot} 
    \bm{T}_k^{t-\Delta t} \bm{q}_{rot}^{-1} - k_r \Delta \bm{\theta}|}\f$ is 
    the cumulative rolling direction unit vector.   

    and 
    
    \f[
    \mbox{\LARGE $\displaystyle \bm{T}_d = \left\{
    \begin{array}{lcl}
    -\eta_{pfr} \eta_r ( \bm{\omega}_0 -
    \bm{\omega}_1 ) & \mathrm{if} & |\bm{T}_k^t| < T_{max} \\
    \bm{0} & \mathrm{if} & |\bm{T}_k^t| = T_{max}
    \end{array} \right.
    $}\f]  
    
    In the above, the four coefficients \f$\mu_r\f$, \f$k_r\f$, \f$\eta_r\f$ and
    \f$\eta_{pfr}\f$ are defined as follows:
    
    - \f$\mu_r\f$ is the Coulomb rolling resistance coefficient. 
    
    - \f$k_r = 3 k_n \mu_r^2\ R^{*,2}\f$ is the rolling stiffness coefficient.
    
    - \f$\eta_r = 3 \gamma_n \mu_r^2\ R^{*,2}\f$ is the rolling dissipative
    coefficient. 
    
    - \f$\eta_{pfr}\f$ is a pre-factor to adjust the amount of rolling 
    dissipative resistance.
    

    The input parameters of the model are:
    - \f$k_n\f$: the normal stiffness coefficient \f$\left(N\cdot 
    	m^{-1}\right)\f$
    - \f$e_n\f$: the restitution coefficient (-)
    - \f$\eta_t\f$: the tangential dissipative coefficient 
    	\f$\left(s^{-1}\right)\f$, automatically computed if not specified or
	assigned a value of \f$-1\f$
    - \f$\mu_c\f$: the Coulomb tangential friction coefficient (-)
    - \f$k_t\f$: the tangential stiffness coefficient \f$\left(N\cdot 
    	m^{-1}\right)\f$
    - \f$\mu_r\f$: the Coulomb rolling coefficient (-)
    - \f$\eta_{pfr}\f$: the pre-factor to adjust the amount of rolling 
    dissipative resistance (-)
                 

    @author D.HUET - 2022
    @author A.WACHS - 2024 - Cleaning & documentation */    
// ============================================================================
class HookeMemoryContactForceModel : public ContactForceModel
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with a map of contact parameters as inputs
    @param parameters map of parameters */
    HookeMemoryContactForceModel( map<string,double>& parameters );

    /** @brief Destructor */
    virtual ~HookeMemoryContactForceModel();
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
    @param dt time step to integrate equations in case an analytical solution is
    not known 
    @param OUT output stream */
    void computeAndWriteEstimates( Component* p0_,  
	Component* p1_,
  	double const& v0, double const& dt,
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
    double m_kn; /**< Normal stiffness coefficient */
    double m_en; /**< Normal restitution coefficient */
    double m_etat; /**< Tangential damping coefficient */
    double m_kt; /**< Tangential stiffness coefficient */
    double m_muc; /**< Tangential Coulomb friction coefficient */
    double m_mur; /**< Angular Coulomb-like friction coefficient. Default value
    	is 0. */
    double m_etarpf; /**< Prefactor in [0,1] of the viscous rolling damping 
    	when the spring rolling friction is not saturated. Default value is 
	0. */
    double m_epst; /**< Threshold below which the tangential vector is
    	undefined and arbitrarily set to 0. Default value is 1.e-10. */
    bool m_rolling_friction; /**< Boolean to switch on/off the rolling 
    	resistance model */
    double m_beta; /**< The log(m_en)/sqrt(PI*PI+log(m_en)*log(m_en)) factor */	
    //@}


    /**@name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    HookeMemoryContactForceModel();
    //@}


    /**@name Methods */
    //@{
    /** @brief Computes maximum penetration depth using a analytical solution
    and a Newton algorithm
    @param theta_ sqrt( omega0*omega0 - mu*mu )
    @param eta_ dissipation coefficient
    @param en_ restitution coefficient
    @param tc_ contact time
    @param v0_ pre-collisional relative velocity */
    double computeDeltaMax( double const& theta_, double const& eta_,
  	double const& en_, double const& tc_, double const& v0_ ) const;

    /** @brief Performs forces & torques computation
    @param p0_ first Component (Particle)
    @param p1_ second Component (Particle ou Obstacle)
    @param dt time step magnitude
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
	
    /** @brief Computes the vector tangent at the contact point
    @param tij vector where the result is stored
    @param n_t eta_t, tangential damping coefficient
    @param ut tangential velocity
    @param kdelta cumulative motion */
    void computeTangentialVector( Vector3& tij, double n_t, Vector3 const& ut,
  	Vector3 const& kdelta );	
    //@}
};

#endif
