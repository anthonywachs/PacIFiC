#ifndef _HERTZCONTACTFORCEMODEL_HH_
#define _HERTZCONTACTFORCEMODEL_HH_

#include "ContactForceModel.hh"
#include "Basic.hh"
#include "Point3.hh"
#include "PointContact.hh"
#include "Vector3.hh"
#include "ReaderXML.hh"
#include <map>
using namespace std;

class Component;


/** The class HertzContactForceModel.

    Contact force model involving a normal Hookean spring, a normal Dashpot and
    a tangential Coulomb friction (HO-D-C) supplemented by a relative velocity 
    dependent rolling resistance torque to compute the force and torque 
    induced by the contact between two rigid components.
    
    The force is formulated as follows:
    
    \f{eqnarray*}{
       \mbox{\LARGE $\bm{F}_n$} &=& \mbox{\LARGE $k_n \delta_n \bm{n}
       	- \gamma_n \bm{U}_n$} \\ 
        \mbox{\LARGE $\bm{F}_t$} &=& \mbox{\LARGE $min\left( 
		| - \gamma_t \bm{U}_t |, 
		\mu_c | \bm{F}_n | \right)\cdot \bm{t}$}
    \f}
    
    where \f$\delta_n\f$ is the overlap distance in the normal direction, 
    \f$\bm{n}\f$ is the unit normal vector at the contact point, \f$\bm{U}_n\f$ 
    is the normal relative velocity at the contact point, \f$\bm{U}_t\f$ is 
    the tangential relative velocity at the contact point and \f$\bm{t}\f$ is 
    the unit tangential vector at the contact point defined as
    \f$\displaystyle -\frac{\bm{U}_t}{|\bm{U}_t|}\f$.
    
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

    - \f$\mu_c\f$ is the Coulomb tangential friction coefficient.

    
    The rolling resistance torque is formulated as follows:
      
    \f[
    \mbox{\LARGE $\displaystyle \bm{T}_r = - k_r R^* |\bm{F}_n| 
    	| \bm{\omega}_0 \times \bm{r}_0 
    	- \bm{\omega}_1 \times \bm{r}_1|\frac{\bm{\omega}_0 - \bm{\omega}_1}
	{|\bm{\omega}_0 - \bm{\omega}_1|}$}
    \f]     
      
    where \f$k_r\f$ is the rolling resistance coefficient, \f$\displaystyle
    \frac{\bm{\omega}_0 - \bm{\omega}_1}{|\bm{\omega}_0 - \bm{\omega}_1|}\f$ is 
    the unit vector along the direction of relative angular velocity and 
    \f$\displaystyle | \bm{\omega}_0 \times \bm{r}_0 
    - \bm{\omega}_1 \times \bm{r}_1|\f$ is the norm of the relative tangential
    velocity contributed by the angular velocities at the contact point. 
    
    The input parameters of the model are:
    - \f$k_n\f$: the normal stiffness coefficient \f$\left(N\cdot 
    	m^{-1}\right)\f$
    - \f$e_n\f$: the restitution coefficient (-)
    - \f$\eta_t\f$: the tangential dissipative coefficient 
    	\f$\left(s^{-1}\right)\f$
    - \f$\mu_c\f$: the Coulomb tangential friction coefficient (-)
    - \f$k_r\f$: the rolling resistance coefficient 
    	\f$\left(s\cdot m^{-1}\right)\f$


    @author A.WACHS - 2024 - Creation */
// ============================================================================
class HertzContactForceModel : public ContactForceModel
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with a map of contact parameters as inputs
    @param parameters map of parameters */
    HertzContactForceModel( map<string,double>& parameters );

    /** @brief Destructor */
    virtual ~HertzContactForceModel();
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
    double m_Es; /**< Average Young modulus */  
    double m_en; /**< Normal restitution coefficient */ 
    double m_Gs; /**< Average shear modulus */
    double m_muc; /**< Tangential Coulomb friction coefficient */
    double m_kr; /**< Rolling resistance coefficient */
    double m_beta; /**< The log(m_en)/sqrt(PI*PI+log(m_en)*log(m_en)) factor */
    double m_m2sqrt56; /**< -2*sqrt(5/6) constant */ 
    //@}


    /**@name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    HertzContactForceModel();
    //@}

  
    /**@name Methods */
    //@{  
    /** @brief Computes the sum of the normal forces divided by the effective
    mass
    @param avmass effective mass 
    @param Req effective radius
    @param deltan overlap distance
    @param v relative velocity */
    double computeDvDt( double const& avmass, double const& Req,
  	double const& deltan, double const& v ) const;
	
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
