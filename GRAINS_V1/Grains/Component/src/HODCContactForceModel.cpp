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
  stiff = parameters["stiff"];
  en    = parameters["en"   ];
  muet  = parameters["mut"  ];
  muec  = parameters["muc"  ];
  k_m_s = parameters["kms"  ];
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

  Vector3 v_n  = normal * ( tmpV * normal );
  Vector3 v_t  = tmpV - v_n;

  // Unit tangential vector along relative velocity at contact point 
  double normv_t = Norm( v_t );
  Vector3 tangent(0.);
  if ( normv_t > EPSILON ) tangent = v_t / normv_t;
  
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
    omega0 = sqrt( 2. * stiff / avmass );
  }
  double muen = - omega0 * log(en) / 
  	sqrt( PI * PI + log(en) * log(en) );    
  delFN += - muen * 2.0 * avmass * v_n;
  double normFN = Norm( delFN );
  
  // Tangential dissipative force
  delFT = v_t * ( -muet * 2.0 * avmass );  

  // Tangential Coulomg saturation
  double fn = muec * normFN;
  double ft = Norm( delFT );
  if ( fn < ft ) delFT = tangent * (-fn);
  
  // Rolling resistance moment
  if ( k_m_s )
  {
    // Relative angular velocity at contact point
    Vector3 w = *p0_->getAngularVelocity() - *p1_->getAngularVelocity();
    Vector3 wn = ( w * normal ) * normal;
    Vector3 wt = w - wn ;
    double normwt = Norm( wt );

    // Anti-spinning effect along the normal wn
    delM = - k_m_s * normFN * 0.001 * wn ;
    
    // Classical rolling resistance moment
//    if ( normwt > EPSILON )  delM -= k_m_s * normFN * wt / normwt;
    if ( normwt > EPSILON )  delM -= k_m_s * normFN * wt;    
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
  bool compute = true ;

  // Pour les particles composite, on calcule la moyenne des forces appliquees
  // a chaque point de contact. Par defaut nbContact = 1 (particles convexes)
  double coef = 1. / double( nbContact );
  
  // Pour les particles composites, on renvoie la reference du composite
  // et non celle des particles elementaires
  Component* ref_p0_ = p0_->getMasterComponent() ;
  Component* ref_p1_ = p1_->getMasterComponent() ; 

  if ( compute  )
  {
    Vector3 delFN, delFT, delM; 
    Point3 geometricPointOfContact = contactInfos.getContact();

    // Calcul des forces & moments de contact
    performForcesCalculus( ref_p0_, ref_p1_, contactInfos, delFN, delFT, delM );

    // Component p0_
    ref_p0_->addForce( geometricPointOfContact, coef * (delFN + delFT) );
    if ( k_m_s ) ref_p0_->addTorque( delM * coef );     
    
    // Component p1_
    ref_p1_->addForce( geometricPointOfContact, coef * ( - delFN - delFT ) );
    if ( k_m_s ) ref_p1_->addTorque( - delM * coef );
    
    // Force postprocessing
    if ( GrainsExec::m_output_data_at_this_time )
      LC->addPPForce( geometricPointOfContact, coef * (delFN + delFT),
	ref_p0_, ref_p1_ );            
  }
  
  return ( compute ) ;  
}




// ----------------------------------------------------------------------------
// Reads and returns contact parameter map from an XML node
map<string,double> HODCContactForceModel::defineParameters( DOMNode* & root )
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
  parameter = ReaderXML::getNode(root, "kms");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["kms"]    = value;
  
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
double HODCContactForceModel::computeDeltaMax( double const& theta_, 
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
