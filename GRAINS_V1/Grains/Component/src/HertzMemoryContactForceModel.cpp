#include "GrainsMPIWrapper.hh"
#include "HertzMemoryContactForceModel.hh"
#include "GrainsExec.hh"
#include "Component.hh"
#include "Memento.hh"
#include "LinkedCell.hh"
#include <math.h>


// ----------------------------------------------------------------------------
// Constructor with a map of contact parameters as inputs
HertzMemoryContactForceModel::HertzMemoryContactForceModel( 
	map<string,double>& parameters ) 
  : ContactForceModel()
{
  m_Es = parameters["Es"];
  m_en = parameters["en"];
  m_Gs = parameters["Gs"];
  m_muc = parameters["muc"];
  m_mur = parameters["mur"];
  m_etarpf = parameters["etarpf"];
  m_epst = parameters["epst"];

  /** If the user defines a rolling friction coefficient, we will compute the
  rolling friction torque. */
  if ( m_mur ) m_rolling_friction = true;
  else m_rolling_friction = false;

  m_beta = log(m_en) / sqrt( PI * PI + log(m_en) * log(m_en) );
  m_m2sqrt56 = - 2. * sqrt( 5. / 6. );  
}




// ----------------------------------------------------------------------------
// Destructor
HertzMemoryContactForceModel::~HertzMemoryContactForceModel()
{}




// ----------------------------------------------------------------------------
string HertzMemoryContactForceModel::name() const 
{ 
  return ( "HertzMemoryContactForceModel" ); 
}




// ----------------------------------------------------------------------------
// Computes the vector tangent at the contact point
void HertzMemoryContactForceModel::computeTangentialVector( Vector3& tij, 
	double n_t, Vector3 const& ut, Vector3 const& kdelta )
{
  // Definition of the tangential vector (cf Costa et.al., 2015)
  if ( Norm( ut ) > 0.001 )
    tij = - ut / Norm( ut );
  else if ( ( Norm( kdelta ) > m_epst || Norm( n_t * ut ) > m_epst ) &&
      Norm( kdelta + n_t * ut ) > m_epst )
    tij = - ( kdelta + n_t * ut ) / Norm( kdelta + n_t * ut ) ;
  else
    tij = 0.;
}




// ----------------------------------------------------------------------------
// Performs forces & torques computation
void HertzMemoryContactForceModel::performForcesCalculus( Component* p0_,
	Component* p1_, double dt, PointContact const& contactInfos,
	Vector3& delFN, Vector3& delFT, Vector3& delM, int elementary_id0,
  	int elementary_id1 )
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
  Vector3* pkdelta = NULL; // previous vector kt * tangential motion
  Vector3* pprev_normal = NULL; // previous vector normal to the contact plane
  Vector3* pspringRotFriction = NULL; // previous rolling friction spring-torque
  Vector3 tij; // tangential vector (lives in the contact plane)
  Quaternion qrot = 0.; // quaternion from tangential plane at previous time 
  	// step to current tangential plane
  std::tuple<int,int,int> idmap0; // id of this contact for particle 0
  std::tuple<int,int,int> idmap1; // id of this contact for particle 1
  idmap0 = std::make_tuple(elementary_id0,p1_->getID(),elementary_id1);
  idmap1 = std::make_tuple(elementary_id1,p0_->getID(),elementary_id0);   	
	
  // Unit normal vector at contact point
  Vector3 normal( penetration );
  normal /= Norm( normal );
  normal.round();

  // Relative velocity at contact point
  Vector3 tmpV = p0_->getVelocityAtPoint( geometricPointOfContact ) 
  	- p1_->getVelocityAtPoint( geometricPointOfContact );	
  Vector3 v_n = normal * ( tmpV * normal );
  Vector3 v_t = tmpV - v_n;


  // 1) Compute normal force  
  // Normal non-linear elastic force
  double kn = ( 4. / 3. ) * m_Es * sqrtReqdeltan;
  delFN = kn * penetration;

  // Normal non-linear dissipative force  
  double gamman = m_m2sqrt56 * m_beta * sqrt( avmass * Sn );  
  delFN -= gamman * v_n;  
  double normFN = Norm( delFN );


  // 2) Compute tangential force with memory  
  // Tangential dissipative force
  double gammat = m_m2sqrt56 * m_beta * sqrt( avmass * St );;
  Vector3 viscousFT = - gammat * v_t ;  

  // Retrieve the previous cumulative motion (if the contact does not 
  // exist we create it)
  bool contact_existed;
  contact_existed = p0_->getContactMemory( idmap0, pkdelta, pprev_normal,
    pspringRotFriction, true );
  if ( contact_existed ) 
  {
    // Rotate it to the new plane and add the contribution of the current time
    // step
    qrot.setRotFromTwoVectors( *pprev_normal, normal ) ;
    *pkdelta = qrot.rotateVector( *pkdelta );
  }
  double kt = 8 * m_Gs * sqrtReqdeltan;
  *pkdelta += dt * kt * v_t; 
  
  // Update the normal vector in the map
  *pprev_normal = normal;  
  
  // Compute a tentative friction force (including the viscous disspative term)
  computeTangentialVector( tij, gammat, v_t, *pkdelta );  
  double normFT = Norm( - *pkdelta + viscousFT );
  
  // If less than the Coulomb limit, we are done
  if ( normFT <= m_muc * normFN ) 
    delFT = Norm( - *pkdelta + viscousFT ) * tij ;
  else 
  {
    // Otherwise, we modify the cumulative tangential motion and replace
    // the tangential force by the Coulomb limit
    *pkdelta = ( - m_muc * normFN * tij + viscousFT ) / kt;
    delFT = m_muc * normFN * tij ;
  } 
  
  
  // 3) Compute rolling resistance torque with memory (if applicable)
  if ( m_rolling_friction ) 
  {
    // Compute the spring and dashpot coefficients from the particle properties
    // following Jiang et al (2005,2015)
    double kr = 3. * kn * m_mur * m_mur * Req * Req ; // torque spring 
	// stiffness
    double etar = 3. * gamman * m_mur * m_mur * Req * Req ; // torque 
    	// dissipative coefficient
    double max_normMk = m_mur * Req * normFN ; // saturation torque
    
    // Compute the relative angular velocity
    Vector3 wrel = *p0_->getAngularVelocity() - *p1_->getAngularVelocity();

    if ( contact_existed ) 
    {
      // Rotate the cumulative spring torque to the new plane    
      *pspringRotFriction = qrot.rotateVector( *pspringRotFriction );      
    }
    *pspringRotFriction -= kr * dt * wrel ;
    double normMk = Norm( *pspringRotFriction );
    if ( normMk > max_normMk ) 
    {
      // Otherwise, we replace the tangential torque by the Coulomb limit
      *pspringRotFriction *= max_normMk / normMk;
    }
    delM = *pspringRotFriction - m_etarpf * etar * wrel;
  }
  else delM = 0.;

  // Finally, we update the cumulative motions in component p1_
  // If contact_existed was false, it also creates the contact in p1_
  p1_->addDeplContactInMap( idmap1,
    - *pkdelta, - normal, - *pspringRotFriction );   
}




// ----------------------------------------------------------------------------
// Computes forces & torques
bool HertzMemoryContactForceModel::computeForces( Component* p0_, 
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

  // In case of composite particles, we retrieve the ids of the elementary
  // components
  int elementary_id0 = 1;
  int elementary_id1 = 1;
  if ( ref_p0_->isCompositeParticle() ) elementary_id0 = p0_->getID() ;
  if ( ref_p1_->isCompositeParticle() ) elementary_id1 = p1_->getID() ;

  // Component tags
  int tag_p0_ = ref_p0_->getTag();
  int tag_p1_ = ref_p1_->getTag();    

  Vector3 delFN, delFT, delM;
  Point3 geometricPointOfContact = contactInfos.getContact();

  // Compute contact force and torque
  performForcesCalculus( ref_p0_, ref_p1_, dt, contactInfos, delFN, delFT, 
  	delM, elementary_id0, elementary_id1 );

  // Component p0_
  ref_p0_->addForce( geometricPointOfContact, coef * (delFN + delFT), tag_p1_ );
  if ( m_rolling_friction ) ref_p0_->addTorque( delM * coef, tag_p1_ );

  // Component p1_
  ref_p1_->addForce( geometricPointOfContact, coef * ( - delFN - delFT ), 
  	tag_p0_ );
  if ( m_rolling_friction ) ref_p1_->addTorque( - delM * coef, tag_p0_ );

  // Force postprocessing
  if ( GrainsExec::m_postprocess_forces_at_this_time )
    LC->addPPForce( geometricPointOfContact, coef * (delFN + delFT),
	ref_p0_, ref_p1_ );

  return ( true ) ;  
}




// ----------------------------------------------------------------------------
// Reads and returns contact parameter map from an XML node
map<string,double> HertzMemoryContactForceModel::defineParameters( DOMNode* & root )
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
  
  parameter = ReaderXML::getNode( root, "Gs" );
  if ( parameter ) value = ReaderXML::getNodeValue_Double( parameter );
  else value = parameters["Es"] / 4.;  
  parameters["Gs"] = value;
  
  parameter = ReaderXML::getNode( root, "mur" );
  if ( parameter )
  {
    value = ReaderXML::getNodeValue_Double( parameter );
    parameters["mur"] = value;
  }
  else parameters["mur"] = 0.;
  
  parameter = ReaderXML::getNode(root, "etarpf");
  if ( parameter )
  {
    value = ReaderXML::getNodeValue_Double( parameter );
    parameters["etarpf"] = value;
  }
  else parameters["etarpf"] = 0.;
  
  parameter = ReaderXML::getNode( root, "epst" );
  if ( parameter )
  {
    value = ReaderXML::getNodeValue_Double( parameter );
    parameters["epst"] = value;
  }
  else parameters["eps"] = 1.e-10;  
  
  return ( parameters );
}




// ----------------------------------------------------------------------------
// Computes an estimate of the contact time and maximum penetration 
// depth in the case of a gravityless binary collision of spheres, and writes
// the result in an output stream
void HertzMemoryContactForceModel::computeAndWriteEstimates( Component* p0_, 
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
double HertzMemoryContactForceModel::computeDvDt( double const& avmass, 
	double const& Req, double const& deltan, double const& v ) const
{
  double sqrtReqdeltan = sqrt( Req * deltan );
  double Sn = 2. * m_Es * sqrtReqdeltan;
  double DvDt = ( - ( 4. / 3. ) * m_Es * sqrtReqdeltan * deltan 
  	- m_m2sqrt56 * m_beta * sqrt( avmass * Sn ) * v ) / avmass;
	 
  return ( DvDt );
}
