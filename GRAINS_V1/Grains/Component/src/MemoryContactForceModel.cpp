#include "GrainsMPIWrapper.hh"
#include "MemoryContactForceModel.hh"
#include "GrainsExec.hh"
#include "Component.hh"
#include "Memento.hh"
#include "LinkedCell.hh"


// ----------------------------------------------------------------------------
// Constructor with a map of contact parameters as inputs
MemoryContactForceModel::MemoryContactForceModel( map<string,double>& parameters )
  : ContactForceModel()
{
  stiff   = parameters["stiff"];
  en      = parameters["en"   ];
  muet    = parameters["mut"  ];
  muec    = parameters["muc"  ];
  ks      = parameters["ks"   ];
  eta_r   = parameters["nr"   ];
  mu_r    = parameters["mur"  ];
  m_f     = parameters["f"    ];
  epsilon = parameters["eps"  ];

  /** If the user defines a rolling friction coefficient, we will compute the
  rolling friction torque. */
  if ( eta_r || mu_r ) rolling_friction = true;
  else rolling_friction = false;
}




// ----------------------------------------------------------------------------
// Destructor
MemoryContactForceModel::~MemoryContactForceModel()
{}




// ----------------------------------------------------------------------------
// Returns the name of the contact force model
string MemoryContactForceModel::name() const
{
  return ( "MemoryContactForceModel" );
}




// ----------------------------------------------------------------------------
// Computes the vector tangent at the contact point
void MemoryContactForceModel::computeTangentialVector( Vector3& tij, 
	double n_t, const Vector3 ut, const Vector3 kdelta )
{
  // Definition of the tangential vector (cf Costa et.al., 2015)
  if ( Norm(ut) > 0.001 )
    tij = - ut/Norm(ut);
  else if ( ( Norm(kdelta) > epsilon || Norm(n_t*ut) > epsilon ) &&
      Norm(kdelta + n_t*ut) > epsilon )
    tij = -( kdelta + n_t*ut ) / Norm(kdelta + n_t*ut) ;
  else
    tij=Vector3(0.);
}




// ----------------------------------------------------------------------------
// Performs forces & torques computation
void MemoryContactForceModel::performForcesCalculus( Component* p0_,
	Component* p1_, double dt, PointContact const& contactInfos,
	Vector3& delFN, Vector3& delFT, Vector3& delM, int elementary_id0,
  	int elementary_id1 )
{
  Point3 geometricPointOfContact = contactInfos.getContact();
  Vector3 penetration = contactInfos.getOverlapVector();
  Vector3 *pkdelta = NULL; // previous vector ks * tangential displacement
  Vector3 *pprev_normal = NULL; // previous vector normal to the contact plane
  Vector3 *pspringRotFriction = NULL; // previous rolling friction spring-torque
  Vector3 tij = Vector3(0.); // tangential vector (lives in the contact plane)
  Quaternion qrot=0.; // quaternion from tangential plane at previous time step
  // to current tangential plane
  Vector3 w ; // relative angular velocity
  Vector3 wn ; // normal relative angular velocity
  Vector3 wt ; // tangential relative angular velocity
  double radius ; // radius of particle 0 (if any)
  double radius1 ; // radius of particle 1
  double Req = 0. ; // effective radius
  double kr = 0.; // spring torque stiffness
  double Cr = 0.; // critical torque (cf. Ai, 2011)
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

  Vector3 v_n  = normal * ( tmpV * normal );
  Vector3 v_t  = tmpV - v_n;

  // Unit tangential vector along relative velocity at contact point
  double normv_t = Norm( v_t );
  Vector3 tangent(0.);
  if ( normv_t > EPSILON ) tangent = v_t / normv_t;

  // 1) Compute normal force
  // Normal linear elastic force
  delFN = stiff * penetration;

  // Normal dissipative force
  double mass0 = p0_->getMass();
  double mass1 = p1_->getMass();
  double avmass = mass0 * mass1 / ( mass0 + mass1 );
  double omega0 = sqrt( stiff / avmass );
  if ( avmass == 0. )
  {
    avmass = mass1 == 0. ? 0.5 * mass0 : 0.5 * mass1;
    Req = p0_->getEquivalentSphereRadius();
    omega0 = sqrt( 2. * stiff / avmass );
  }
  double muen = - omega0 * log(en) /
  	sqrt( PI * PI + log(en) * log(en) );
  delFN += - muen * 2.0 * avmass * v_n;
  double normFN = Norm( delFN );

  // 2) Compute tangential force with memory
  // Tangential viscous dissipative force
  Vector3 viscousFT = v_t * ( -muet * 2.0 * avmass );

  // Retrieve the previous cumulative displacement (if the contact does not 
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
  *pkdelta += dt*ks*v_t;
  
  // Update the normal vector in the map
  *pprev_normal = normal;  
  
  // Compute a tentative friction force (including the viscous disspative terms)
  computeTangentialVector( tij, muet * 2.0 * avmass, v_t, *pkdelta );
  double normFT = Norm(-*pkdelta + viscousFT);
  
  // If less than the Coulomb limit, we are done
  if ( normFT <= muec * normFN ) 
    delFT = Norm(-*pkdelta + viscousFT) * tij ;
  else 
  {
    // Otherwise, we modify the cumulative tangential displacement and replace
    // the tangential force by the Coulomb limit
    *pkdelta = - muec * normFN * tij + viscousFT;
    delFT = muec * normFN * tij ;
  }

  // 3) Compute rolling resistance moment with memory (if applicable)
  if ( rolling_friction ) 
  {
    // Compute the spring and dashpot coefficients from the particle properties
    radius = p0_->getEquivalentSphereRadius() ;
    if ( !Req )
    {
      radius1 = p1_->getEquivalentSphereRadius() ;
      Req = radius * radius1 / ( radius + radius1 ) ;
    }
    kr = 3. * stiff * mu_r * mu_r * Req * Req ; // c.f. Jiang et al (2005,2015)
    Cr = 3. * (muen * 2.0 * avmass) * mu_r * mu_r * Req * Req ; // c.f. Jiang et
    // al (2005,2015)
    double max_normFT = mu_r * Req * normFN ;
    
    // Compute the relative angular velocity
    w = *p0_->getAngularVelocity() - *p1_->getAngularVelocity();
    wn = ( w * normal ) * normal;
    wt = w - wn ;

    if ( contact_existed ) 
    {
      // Rotate the cumulative spring torque to the new plane
      *pspringRotFriction = qrot.rotateVector( *pspringRotFriction );
    }
    *pspringRotFriction += (-kr) * dt * wt ;
    normFT = Norm(*pspringRotFriction);
    if ( normFT > max_normFT ) 
    {
      // Otherwise, we replace the tangential torque by the Coulomb limit
      *pspringRotFriction *= max_normFT / normFT;
    }
    delM = *pspringRotFriction + (-m_f) * Cr * wt;
  }
  else delM = Vector3(0.);

  // Finally, we update the cumulative displacements in component p1_
  // If contact_existed was false, it also creates the contact in p1_
  p1_->addDeplContactInMap( idmap1,
    -*pkdelta, -normal, -*pspringRotFriction );
}




// ----------------------------------------------------------------------------
// Computes forces & torques
bool MemoryContactForceModel::computeForces( Component* p0_,
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

  Vector3 delFN, delFT, delM;
  Point3 geometricPointOfContact = contactInfos.getContact();

  // Calcul des forces & moments de contact
  performForcesCalculus( ref_p0_, ref_p1_, dt, contactInfos, delFN, delFT, 
  	delM, elementary_id0, elementary_id1 );

  // Component p0_
  ref_p0_->addForce( geometricPointOfContact, coef * (delFN + delFT) );
  if ( rolling_friction ) ref_p0_->addTorque( delM * coef );

  // Component p1_
  ref_p1_->addForce( geometricPointOfContact, coef * ( - delFN - delFT ) );
  if ( rolling_friction ) ref_p1_->addTorque( - delM * coef );

  // Force postprocessing
  if ( GrainsExec::m_output_data_at_this_time )
    LC->addPPForce( geometricPointOfContact, coef * (delFN + delFT),
	ref_p0_, ref_p1_ );

  return ( true ) ;
}




// ----------------------------------------------------------------------------
// Reads and returns contact parameter map from an XML node
map<string,double> MemoryContactForceModel::defineParameters( DOMNode* & root )
{
  map<string,double> parameters;

  DOMNode* parameter;
  double   value;
  parameter = ReaderXML::getNode(root, "stiff");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["stiff"]  = value;
  parameter = ReaderXML::getNode(root, "muc");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["muc"]    = value;
  parameter = ReaderXML::getNode(root, "en");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["en"]    = value;
  parameter = ReaderXML::getNode(root, "mut");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["mut"]    = value;
  // TODO: define default values in the DTD file instead?
  parameter = ReaderXML::getNode(root, "nr");
  if ( parameter )
  {
    value     = ReaderXML::getNodeValue_Double(parameter);
    parameters["nr"]    = value;
  }
  else parameters["nr"] = 0.;
  parameter = ReaderXML::getNode(root, "mur");
  if ( parameter )
  {
    value     = ReaderXML::getNodeValue_Double(parameter);
    parameters["mur"]    = value;
  }
  else parameters["mur"] = 0.;
  parameter = ReaderXML::getNode(root, "f");
  if ( parameter )
  {
    value     = ReaderXML::getNodeValue_Double(parameter);
    parameters["f"]    = value;
  }
  else parameters["f"] = 0.;
  parameter = ReaderXML::getNode(root, "ks");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["ks"]    = value;
  parameter = ReaderXML::getNode(root, "eps");
  if ( parameter )
  {
    value     = ReaderXML::getNodeValue_Double(parameter);
    parameters["eps"]  = value;
  }
  else parameters["eps"] = 1.e-10;

  return ( parameters );
}




// ----------------------------------------------------------------------------
// Computes an estimate of the contact time and maximum penetration
// depth in the case of a gravityless binary collision of spheres, and writes
// the result in an output stream
void MemoryContactForceModel::computeAndWriteEstimates( Component* p0_,
	Component* p1_, double const& v0, ostream& OUT ) const
{
  double mass0 = p0_->getMass();
  double mass1 = p1_->getMass();
  double avmass = mass0 * mass1 / (mass0 + mass1);
  if (avmass == 0.) avmass = mass1==0. ? 0.5*mass0 : 0.5*mass1;
  double muen = - sqrt(stiff/avmass) * log(en) / sqrt( PI*PI+log(en)*log(en) );

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
    double omega0 = sqrt(stiff/mass0);
    double theta = sqrt(pow(omega0,2.) - pow(0.5*muen,2.));
    double tc = PI / theta;

    double delta_max = computeDeltaMax( theta, 0.5*muen, en, tc, v0 );

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
    OUT << "  mu_n = " << muen << endl;
    OUT << "  Maximum elastic force fel = " << stiff*delta_allowed << endl;
    OUT << "  fel/weight0 ratio = " << stiff*delta_allowed/(mass0*9.81)
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
    double omega0 = sqrt(stiff/avmass);
    double theta = sqrt(pow(omega0,2.) - pow(muen,2.));
    double tc = PI / theta;

    double delta_max = computeDeltaMax( theta, muen, en, tc, v0 );

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
    OUT << "  mu_n = " << muen << endl;
    OUT << "  Maximum elastic force fel = " << stiff*delta_allowed << endl;
    OUT << "  fel/weight0 ratio = " << stiff*delta_allowed/(mass0*9.81)
  	<< endl;
    OUT << "  fel/weight1 ratio = " << stiff*delta_allowed/(mass1*9.81)
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
double MemoryContactForceModel::computeDeltaMax( double const& theta_,
	double const& mu_, double const& en_, double const& tc_,
	double const& v0_ ) const
{
  double f=1.,df,t0 = tc_ / (en_ > 0.2 ? 2. : 100.);
  while ( fabs(f) > 1e-10 )
  {
    f = ( v0_ / theta_ ) * exp( - mu_ * t0 ) * ( - mu_ * sin ( theta_ * t0 )
    	+ theta_ * cos ( theta_ * t0 ) );
    df = ( v0_ / theta_ ) * exp( - mu_ * t0 ) *
    	( pow( mu_, 2. ) * sin ( theta_ * t0 )
	- 2. * mu_ * theta_ * cos ( theta_ * t0 )
	- pow( theta_, 2. ) * sin ( theta_ * t0 ) );
    t0 -= f / df ;
  }

  double delta_max = ( v0_ / theta_ ) * exp( - mu_ * t0 )
    	* sin( theta_ * t0 );

  return ( delta_max );
}
