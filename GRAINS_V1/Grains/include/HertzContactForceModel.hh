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

    Contact force model involving a normal non-linear spring, a normal
    non-linear Dashpot and a tangential Coulomb friction supplemented by a 
    relative velocity dependent rolling resistance torque to compute the force 
    and torque induced by the contact between two rigid components.
    
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
    - \f$k_n\f$ is the normal stiffness coefficient defined as
    
    \f[
    \mbox{\LARGE $\displaystyle k_n = \frac{4}{3}E^*\sqrt{R^* \delta_n}$}
    \f]
    
    where \f$\displaystyle E^* = \frac{1}{\frac{1-\nu_0^2}{E_0}
    +\frac{1-\nu_1^2}{E_1}}\f$ 
    and \f$R^*=\frac{1}{\frac{1}{R_0}+\frac{1}{R_1}}\f$ are the effective Young
    modulus and effective radius of curvature, respectively. \f$R_0\f$ is the 
    radius of curvature of the first rigid component and \f$R_1\f$ is 
    the radius of curvature of the second rigid component. \f$E_0\f$ is the 
    Young modulus of the first rigid component and \f$E_1\f$ is 
    the Young modulus of the second rigid component. Finally, \f$\nu_0\f$ is 
    the Poisson ratio of the first rigid component and \f$\nu_1\f$ is 
    the Poisson ratio of the second rigid component.      
    
    - \f$\gamma_n\f$ is the normal dissipative coefficient defined as
    
    \f[
    \mkern-18mu\mbox{\LARGE $\displaystyle \gamma_n = - 2\sqrt{\frac{5}{6}}
    	\beta\sqrt{S_n m^*} 
    	\:\:\mathrm{with}\:\: 
    	\beta=\frac{\ln(e_n)}{\sqrt{\pi^2+(\ln(e_n))^2}}
	\:\:\mathrm{,}\:\:
	S_n = 2E^*\sqrt{R^* \delta_n}$}
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
    \mbox{\LARGE $\gamma_t = - 2\sqrt{\frac{5}{6}}
    	\beta\sqrt{S_t m^*}\:\:\mathrm{with}\:\:
	S_t = 8G^*\sqrt{R^* \delta_n}$}
    \f]
    
    where \f$\displaystyle G^* = \frac{1}{\frac{2(2-\nu_0)(1+\nu_0)}{E_0}
    +\frac{2(2-\nu_1)(1+\nu_1)}{E_1}}\f$ is the effective shear modulus.     

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

    Note that for most materials, the value of the Poisson ratio belongs to the 
    interval \f$[0.1,0,4]\f$. In this interval, we have \f$1-\nu^2 \approx 2\f$ 
    and \f$2(2-\nu)(1+\nu)\approx 4\f$, such that if the two components in 
    contact are made of the same material, the following simplifying 
    approximations hold:
    
    \f[
    \mbox{\LARGE $E^*\approx\frac{E}{2}
    	\:\:\mathrm{,}\:\:G^*\approx\frac{E}{8}
	\:\:\Rightarrow\:\:G^*\approx 4E^*$}
    \f]     
    
    The input parameters of the model are:
    - \f$E^*\f$: the effective Young modulus \f$\left(N\cdot 
    	m^{-2}=Pa\right)\f$
    - \f$e_n\f$: the restitution coefficient (-)
    - \f$G^*\f$: the effective shear modulus \f$\left(N\cdot 
    	m^{-2}=Pa\right)\f$, automatically computed as \f$4E^*\f$ if not 
	specified 
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
