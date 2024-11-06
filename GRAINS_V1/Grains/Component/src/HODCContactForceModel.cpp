#include "GrainsMPIWrapper.hh"
#include "HODCContactForceModel.hh"
#include "GrainsExec.hh"
#include "Component.hh"
#include "Memento.hh"
#include "LinkedCell.hh"


// ----------------------------------------------------------------------------
// Constructor with a map of contact parameters as inputs
HODCContactForceModel::HODCContactForceModel( map<string,double>& parameters ) 
  : ContactForceModel()
{
  m_kn = parameters["kn"];
  m_en = parameters["en"];
  m_etat = parameters["etat"];
  m_muc = parameters["muc"];
  m_kr = parameters["kr"];
  m_beta = log(m_en) / sqrt( PI * PI + log(m_en) * log(m_en) );
}




// ----------------------------------------------------------------------------
// Destructor
HODCContactForceModel::~HODCContactForceModel()
{}




// ----------------------------------------------------------------------------
string HODCContactForceModel::name() const 
{ 
  return ( "HODCContactForceModel" ); 
}




// ----------------------------------------------------------------------------
// Performs forces & torques computation
void HODCContactForceModel::performForcesCalculus( Component* p0_, 
	Component* p1_, PointContact const& contactInfos,
	Vector3& delFN, Vector3& delFT, Vector3& delM )
{
  Point3 geometricPointOfContact = contactInfos.getContact();
  Vector3 penetration = contactInfos.getOverlapVector();

  // Unit normal vector at contact point
  Vector3 normal( penetration );
  normal /= Norm( normal );
  normal.round();

  // Relative velocity at contact point
  Vector3 tmpV = p0_->getVelocityAtPoint( geometricPointOfContact ) 
  	- p1_->getVelocityAtPoint( geometricPointOfContact );	

  Vector3 v_n = normal * ( tmpV * normal );
  Vector3 v_t = tmpV - v_n;

  // Unit tangential vector along relative velocity at contact point 
  double normv_t = Norm( v_t );
  Vector3 tangent(0.);
  if ( normv_t > EPSILON ) tangent = v_t / normv_t;
  
  // Normal linear elastic force
  delFN = m_kn * penetration;

  // Normal dissipative force  
  double mass0 = p0_->getMass();
  double mass1 = p1_->getMass();  
  double avmass = 1. / ( 1. / mass0 + 1. / mass1 );
  double gamman = - 2. * m_beta * sqrt( avmass * m_kn );
  delFN -= gamman * v_n;  
  double normFN = Norm( delFN );
  
  // Tangential dissipative force
  double gammat = 2. * m_etat * avmass;
  delFT = - gammat * v_t ;  

  // Tangential Coulomg saturation
  double fn = m_muc * normFN;
  double ft = Norm( delFT );
  if ( fn < ft ) delFT = tangent * ( - fn );
  
  // Rolling resistance moment
  if ( m_kr )
  {    
    // Relative angular velocity at contact point
    Vector3 wrel = *p0_->getAngularVelocity() - *p1_->getAngularVelocity();
    double normwrel = Norm( wrel );
    double normwtrel = Norm( ( *p0_->getAngularVelocity() ^ (
    	*p0_->getPosition() - geometricPointOfContact ) )
	- ( *p1_->getAngularVelocity() ^ ( *p1_->getPosition()
    	- geometricPointOfContact ) ) );
    
    // Rolling resistance moment
    if ( normwrel > EPSILON )
    {  
      double Req = 1. / ( 1. / p0_->getEquivalentSphereRadius() 
    	+ 1. / p1_->getEquivalentSphereRadius() ); // Effective radius      
      delM = - ( m_kr * Req * normFN * normwtrel / normwrel ) * wrel ;
    }   
  }
  else delM = 0.0;  
}




// ----------------------------------------------------------------------------
// Computes forces & torques
bool HODCContactForceModel::computeForces( Component* p0_, 
	Component* p1_,
	PointContact const& contactInfos,
	LinkedCell* LC,
	double dt, int nbContact )
{
  // In the case of composite particles, we compute the average of the 
  // contact force over all contact points between two rigid bodies
  // By default, nbContact = 1
  double coef = 1. / double( nbContact );
  
  // In the case of composite particles, we use the composite particle not the
  // elemetary particles
  Component* ref_p0_ = p0_->getMasterComponent() ;
  Component* ref_p1_ = p1_->getMasterComponent() ;

  // Component tags
  int tag_p0_ = ref_p0_->getTag();
  int tag_p1_ = ref_p1_->getTag();   

  Vector3 delFN, delFT, delM; 
  Point3 geometricPointOfContact = contactInfos.getContact();

  // Compute contact force and torque
  performForcesCalculus( ref_p0_, ref_p1_, contactInfos, delFN, delFT, delM );

  // Component p0_
  ref_p0_->addForce( geometricPointOfContact, coef * ( delFN + delFT ), 
  	tag_p1_ );
  if ( m_kr ) ref_p0_->addTorque( delM * coef, tag_p1_ );     
    
  // Component p1_
  ref_p1_->addForce( geometricPointOfContact, coef * ( - delFN - delFT ),
  	 tag_p0_ );
  if ( m_kr ) ref_p1_->addTorque( - delM * coef, tag_p0_ );
    
  // Force postprocessing
  if ( GrainsExec::m_postprocess_forces_at_this_time )
    LC->addPPForce( geometricPointOfContact, coef * ( delFN + delFT ),
	ref_p0_, ref_p1_ );            
  
  return ( true ) ;  
}




// ----------------------------------------------------------------------------
// Reads and returns contact parameter map from an XML node
map<string,double> HODCContactForceModel::defineParameters( DOMNode* & root )
{
  map<string,double> parameters;

  DOMNode* parameter;
  double value;
  
  parameter = ReaderXML::getNode( root, "kn" );
  value = ReaderXML::getNodeValue_Double( parameter );
  parameters["kn"] = value;
  
  parameter = ReaderXML::getNode( root, "muc" );
  value = ReaderXML::getNodeValue_Double( parameter );
  parameters["muc"] = value;
  
  parameter = ReaderXML::getNode( root, "en" );
  value = ReaderXML::getNodeValue_Double( parameter );
  parameters["en"] = value;
  
  parameter = ReaderXML::getNode( root, "etat" );
  value = ReaderXML::getNodeValue_Double( parameter );
  parameters["etat"] = value;
  
  parameter = ReaderXML::getNode( root, "kr" );
  value = ReaderXML::getNodeValue_Double( parameter );
  parameters["kr"] = value;
  
  return ( parameters );
}




// ----------------------------------------------------------------------------
// Computes an estimate of the contact time and maximum penetration 
// depth in the case of a gravityless binary collision of spheres, and writes
// the result in an output stream
void HODCContactForceModel::computeAndWriteEstimates( Component* p0_, 
	Component* p1_, double const& v0, ostream& OUT ) const
{
  double mass0 = p0_->getMass();
  double mass1 = p1_->getMass();
  double avmass = mass0 * mass1 / (mass0 + mass1);
  if (avmass == 0.) avmass = mass1==0. ? 0.5*mass0 : 0.5*mass1;
  double etan = - sqrt(m_kn/avmass) * log(m_en) 
  	/ sqrt( PI*PI+log(m_en)*log(m_en) );

  // Particle/obstacle contact
  if ( p0_->isObstacle() || p1_->isObstacle() )
  {
    Component* particle = NULL;
    Component* obstacle = NULL;
    if ( p0_->isObstacle() )
    {
      particle = p1_;
      obstacle = p0_ ;
    }
    else
    {
      particle = p0_;
      obstacle = p1_ ;
    }    
    mass0 = particle->getMass();
    double delta_allowed = p0_->getCrustThickness()
  	+ p1_->getCrustThickness();
    double omega0 = sqrt(m_kn/mass0); 
    double theta = sqrt(pow(omega0,2.) - pow(0.5*etan,2.));
    double tc = PI / theta;
  
    double delta_max = computeDeltaMax( theta, 0.5*etan, m_en, tc, v0 );

    OUT << "  Particle: material = " << particle->getMaterial()
  	<< " crust thickness = " << particle->getCrustThickness() 
	<< " weight = " << mass0*9.81 << endl;
    OUT << "  Obstacle: material = " << obstacle->getMaterial()
  	<< " crust thickness = " << obstacle->getCrustThickness() 
	<< " weight = " << mass0*9.81 << endl;	
    OUT << "  Maximum overlap allowed by crusts = " 
    	<< delta_allowed << endl;
    OUT << "  Collisional relative velocity = " << v0 << endl;
    OUT << "  Tc = " << tc << "  delta_max = " << delta_max << endl; 
    OUT << "  eta_n = " << etan << endl;        
    OUT << "  Maximum elastic force fel = " << m_kn*delta_allowed << endl;
    OUT << "  fel/weight0 ratio = " << m_kn*delta_allowed/(mass0*9.81) 
  	<< endl;
    if ( delta_max > delta_allowed ) 
    {
      OUT << "  *********************************************" << endl;
      OUT << "  !!!!! WARNING !!!!!" << endl;
      OUT << "  delta_max > maximum overlap allowed by crusts" << endl;
      OUT << "  *********************************************" << endl;
    }	
  }
  // Particle/particle contact 
  else
  {
    double delta_allowed = p0_->getCrustThickness()
  	+ p1_->getCrustThickness();
    double omega0 = sqrt(m_kn/avmass); 
    double theta = sqrt(pow(omega0,2.) - pow(etan,2.));
    double tc = PI / theta;
  
    double delta_max = computeDeltaMax( theta, etan, m_en, tc, v0 );

    OUT << "  Particle 0: material = " << p0_->getMaterial()
  	<< " crust thickness = " << p0_->getCrustThickness() 
	<< " weight = " << mass0*9.81 << endl;
    OUT << "  Particle 1: materiau = " << p1_->getMaterial()
  	<< " crust thickness = " << p1_->getCrustThickness()
	<< " weight = " << mass1*9.81 << endl;
    OUT << "  Maximum overlap allowed by crusts = " << delta_allowed 
    	<< endl;
    OUT << "  Collisional relative velocity = " << v0 << endl;
    OUT << "  Tc = " << tc << "  delta_max = " << delta_max << endl;    
    OUT << "  eta_n = " << etan << endl; 
    OUT << "  Maximum elastic force fel = " << m_kn*delta_allowed << endl;
    OUT << "  fel/weight0 ratio = " << m_kn*delta_allowed/(mass0*9.81) 
  	<< endl;
    OUT << "  fel/weight1 ratio = " << m_kn*delta_allowed/(mass1*9.81)
  	<< endl;
    if ( delta_max > delta_allowed ) 
    {
      OUT << "  *********************************************" << endl;
      OUT << "  !!!!! WARNING !!!!!" << endl;
      OUT << "  delta_max > maximum overlap allowed by crusts" << endl;
      OUT << "  *********************************************" << endl;
    }
  }
  OUT << endl;        
}




// ----------------------------------------------------------------------------
// Computes maximum penetration depth using a analytical solution
// and a Newton algorithm
double HODCContactForceModel::computeDeltaMax( double const& theta_, 
	double const& eta_, double const& en_, double const& tc_, 
	double const& v0_ ) const
{
  double f=1.,df,t0 = tc_ / (en_ > 0.2 ? 2. : 100.);
  while ( fabs(f) > 1e-10 )
  {
    f = ( v0_ / theta_ ) * exp( - eta_ * t0 ) * ( - eta_ * sin ( theta_ * t0 )
    	+ theta_ * cos ( theta_ * t0 ) );
    df = ( v0_ / theta_ ) * exp( - eta_ * t0 ) * 
    	( pow( eta_, 2. ) * sin ( theta_ * t0 )
	- 2. * eta_ * theta_ * cos ( theta_ * t0 )
	- pow( theta_, 2. ) * sin ( theta_ * t0 ) );
    t0 -= f / df ;	
  } 

  double delta_max = ( v0_ / theta_ ) * exp( - eta_ * t0 ) 
    	* sin( theta_ * t0 );
	
  return ( delta_max );
}
