#include "GrainsMPIWrapper.hh"
#include "HertzContactForceModel.hh"
#include "GrainsExec.hh"
#include "Component.hh"
#include "Memento.hh"
#include "LinkedCell.hh"
#include <math.h>


// ----------------------------------------------------------------------------
// Constructor with a map of contact parameters as inputs
HertzContactForceModel::HertzContactForceModel( map<string,double>& parameters ) 
  : ContactForceModel()
{
  m_Es = parameters["Es"];
  m_en = parameters["en"];
  m_Gs = parameters["Gs"];
  m_muc = parameters["muc"];
  m_kr = parameters["kr"];
  m_beta = log(m_en) / sqrt( PI * PI + log(m_en) * log(m_en) );
  m_m2sqrt56 = - 2. * sqrt( 5. / 6. );
}




// ----------------------------------------------------------------------------
// Destructor
HertzContactForceModel::~HertzContactForceModel()
{}




// ----------------------------------------------------------------------------
string HertzContactForceModel::name() const 
{ 
  return ( "HertzContactForceModel" ); 
}




// ----------------------------------------------------------------------------
// Performs forces & torques computation
void HertzContactForceModel::performForcesCalculus( Component* p0_, 
	Component* p1_, PointContact const& contactInfos,
	Vector3& delFN, Vector3& delFT, Vector3& delM )
{
  Point3 geometricPointOfContact = contactInfos.getContact();
  Vector3 penetration = contactInfos.getOverlapVector();
  double deltan = - contactInfos.getOverlapDistance();
  double Req = 1. / ( 1. / p0_->getEquivalentSphereRadius() 
    	+ 1. / p1_->getEquivalentSphereRadius() ); // Effective radius  
  double avmass = 1. / ( 1. / p0_->getMass() 
  	+ 1. / p1_->getMass() ); // Effective mass
  double sqrtReqdeltan = sqrt( Req * deltan );
  double Sn = 2. * m_Es * sqrtReqdeltan;
  double St = 8. * m_Gs * sqrtReqdeltan;  	
	
  // Unit normal vector at contact point
  Vector3 normal( penetration );
  normal /= Norm( normal );
  normal.round();

  // Relative velocity at contact point
  Vector3 tmpV = p0_->getVelocityAtPoint( geometricPointOfContact ) 
  	- p1_->getVelocityAtPoint( geometricPointOfContact );	

  Vector3 v_n = normal * ( tmpV * normal );
  Vector3 v_t = tmpV - v_n;

  // Unit tangential vector in the reverse direction of the relative velocity 
  // at contact point 
  double normv_t = Norm( v_t );
  Vector3 tangent( 0. );
  if ( normv_t > EPSILON ) tangent = - v_t / normv_t;
  
  // Normal linear elastic force
  double kn = ( 4. / 3. ) * m_Es * sqrtReqdeltan;
  delFN = kn * penetration;

  // Normal dissipative force  
  double gamman = m_m2sqrt56 * m_beta * sqrt( avmass * Sn );  
  delFN -= gamman * v_n;  
  double normFN = Norm( delFN );
  
  // Tangential dissipative force
  double gammat = m_m2sqrt56 * m_beta * sqrt( avmass * St );;
  delFT = - gammat * v_t ;  

  // Tangential Coulomg saturation
  double fn = m_muc * normFN;
  double ft = Norm( delFT );
  if ( fn < ft ) delFT = tangent * fn ;
  
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
      delM = - ( m_kr * Req * normFN * normwtrel / normwrel ) * wrel ;  
  }
  else delM = 0.0;  
}




// ----------------------------------------------------------------------------
// Computes forces & torques
bool HertzContactForceModel::computeForces( Component* p0_, 
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
map<string,double> HertzContactForceModel::defineParameters( DOMNode* & root )
{
  map<string,double> parameters;

  DOMNode* parameter;
  double value;
  
  parameter = ReaderXML::getNode( root, "Es" );
  value = ReaderXML::getNodeValue_Double( parameter );
  parameters["Es"] = value;

  parameter = ReaderXML::getNode( root, "en" );
  value = ReaderXML::getNodeValue_Double( parameter );
  parameters["en"] = value;

  parameter = ReaderXML::getNode( root, "muc" );
  value = ReaderXML::getNodeValue_Double( parameter );
  parameters["muc"] = value;
  
  parameter = ReaderXML::getNode( root, "kr" );
  value = ReaderXML::getNodeValue_Double( parameter );
  parameters["kr"] = value;  
  
  parameter = ReaderXML::getNode( root, "Gs" );
  if ( parameter ) value = ReaderXML::getNodeValue_Double( parameter );
  else value = parameters["Es"] / 4.;  
  parameters["Gs"] = value;
  
  return ( parameters );
}




// ----------------------------------------------------------------------------
// Computes an estimate of the contact time and maximum penetration 
// depth in the case of a gravityless binary collision of spheres, and writes
// the result in an output stream
void HertzContactForceModel::computeAndWriteEstimates( Component* p0_, 
	Component* p1_, double const& v0, double const& dt, ostream& OUT ) const
{
  double Req = 1. / ( 1. / p0_->getEquivalentSphereRadius() 
    	+ 1. / p1_->getEquivalentSphereRadius() ); // Effective radius  
  double delta_allowed = p0_->getCrustThickness()
  	+ p1_->getCrustThickness();
  double mass0 = p0_->getMass();
  double mass1 = p1_->getMass();	
  double avmass = 1. / ( 1. / mass0 + 1. / mass1 ); // Effective mass	
  double t = 0, deltan = 0., v = v0, k1d, k1v, k2d, k2v, 
  	delta_max = -1.e20, maxfel = 0.;

  // We integrate the equations of motion with an explicit RK2 scheme
  while ( deltan >= 0. )
  {
    k1d = v * dt;
    k1v = computeDvDt( avmass, Req, deltan, v ) * dt;
    k2d = ( v + k1v / 2. ) * dt;
    k2v = computeDvDt( avmass, Req, deltan + k1d / 2., v + k1v / 2. ) * dt;
    deltan += k2d;
    v += k2v; 
    t += dt;
    delta_max = max( deltan, delta_max );    
  }
  maxfel = ( 4. / 3. ) * m_Es * sqrt( Req ) * pow( delta_max, 1.5 );
   
  OUT << "  Component 0: material = " << p0_->getMaterial()
  	<< " crust thickness = " << p0_->getCrustThickness() 
	<< " weight = " << mass0 * 9.81 << endl;
  OUT << "  Component 1: material = " << p1_->getMaterial()
  	<< " crust thickness = " << p1_->getCrustThickness()
	<< " weight = " << mass1 * 9.81 << endl;
  OUT << "  Maximum overlap allowed by crusts = " << delta_allowed 
    	<< endl;
  OUT << "  Collisional relative velocity = " << v0 << endl;
  OUT << "  Tc = " << t << "  delta_max = " << delta_max << endl;    
  OUT << "  Maximum elastic force fel = " << maxfel << endl;
  if ( !p0_->isObstacle() )
    OUT << "  fel/weight0 ratio = " << maxfel / ( mass0 * 9.81 ) 
  	<< endl;
  if ( !p1_->isObstacle() )
    OUT << "  fel/weight1 ratio = " << maxfel / ( mass1 * 9.81 )
  	<< endl;
  if ( delta_max > delta_allowed ) 
  {
    OUT << "  *********************************************" << endl;
    OUT << "  !!!!! WARNING !!!!!" << endl;
    OUT << "  delta_max > maximum overlap allowed by crusts" << endl;
    OUT << "  *********************************************" << endl;
  }
  OUT << endl;  
  
}




// ----------------------------------------------------------------------------
// Computes the sum of the normal forces divided by the effective mass
double HertzContactForceModel::computeDvDt( double const& avmass, 
	double const& Req, double const& deltan, double const& v ) const
{
  double sqrtReqdeltan = sqrt( Req * deltan );
  double Sn = 2. * m_Es * sqrtReqdeltan;
  double DvDt = ( - ( 4. / 3. ) * m_Es * sqrtReqdeltan * deltan 
  	- m_m2sqrt56 * m_beta * sqrt( avmass * Sn ) * v ) / avmass;
	 
  return ( DvDt );
}
