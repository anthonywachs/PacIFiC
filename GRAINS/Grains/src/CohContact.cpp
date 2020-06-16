#include "MPIWrapperGrains.hh"
#include "CohContact.H"
#include "Grains_Exec.hh"
#include "Composant.H"
#include "Memento.hh"
#include "LinkedCell.H"



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur avec un lot de parametres
CohContact::CohContact( map<string,double>& parameters ) :
  ContactLaw()
{
  stiff = parameters["stiff"];
  stiff_t = parameters["stiff_t"];
  en    = parameters["en"   ];
  muet  = parameters["mut"  ];
  muec  = parameters["muc"  ];
  k_m_s = parameters["kms"  ];
  rupt_dist = parameters["ruptdist"  ];
  initialFmax_dist = parameters["initialmaxdist"  ];
  Grains_Exec::m_withCohesion = true;
}




// ----------------------------------------------------------------------------
// @brief Destructeur
CohContact::~CohContact()
{
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul effectif des forces & moments de contact
// It also sends back the dissipated energy for pp purposes
void CohContact::performForcesCalculusPP( Composant* p0_,
	Composant* p1_, const PointContact &contactInfos,
	Vecteur &delFN, Vecteur &delFT, Vecteur &delM, Scalar &diss )
{
  // Calcul des forces & moments de contact
  // --------------------------------------
  Point geometricPointOfContact = contactInfos.getContact();
  Vecteur penetration = contactInfos.getRecouvrement();

  // Vecteur normal au contact
  Vecteur normal( penetration );
  normal /= Norm( normal );
  normal.round();

  // Calcul de la vitesse relative au point de contact
  Vecteur tmpV = p0_->getVitesse( geometricPointOfContact )
  	- p1_->getVitesse( geometricPointOfContact );

  Vecteur v_n  = normal * ( tmpV * normal );
  Vecteur v_t  = tmpV - v_n;

  // Calcul de la tangente portee par la vitesse relative tangente
  Scalar normv_t = Norm( v_t );
  Vecteur tangent(0.);
  if ( normv_t > EPS ) tangent = v_t / normv_t;

  // Force de restitution elastique lineaire
  delFN = stiff * penetration;

  // Force de restitution elastique normale + Force normale amortie
  double masse0 = p0_->getMasse();
  double masse1 = p1_->getMasse();
  Scalar avmass = masse0 * masse1 / ( masse0 + masse1 );
  Scalar omega0 = sqrt( stiff / avmass );
  if ( avmass == 0. )
  {
    avmass = masse1 == 0. ? 0.5 * masse0 : 0.5 * masse1;
    omega0 = sqrt( 2. * stiff / avmass );
  }
  Scalar muen = - omega0 * log(en) /
  	sqrt( PI * PI + log(en) * log(en) );
  delFN += - muen * 2.0 * avmass * v_n;
  Scalar normFN = Norm( delFN );
  diss =  v_n * (muen * 2.0 * avmass * v_n);

  // Force tangentielle de dissipation
  delFT = v_t * ( -muet * 2.0 * avmass );

  // Force tangentielle par le critere de Coulomb
  Scalar fn = muec * normFN;
  Scalar ft = Norm( delFT );
  if ( fn < ft ) delFT = tangent * (-fn);
  diss +=  - v_t * delFT;

  // Moment de friction de roulement
  if ( k_m_s )
  {
    // Calcul de la vitesse de rotation relative
    Vecteur w = *p0_->getVitesseRotation() - *p1_->getVitesseRotation();
    Vecteur wn = ( w * normal ) * normal;
    Vecteur wt = w - wn ;
    Scalar normwt = Norm( wt );

    // Correction suivant la normale : effet anti-toupie
    delM = - k_m_s * normFN * 0.001 * wn ;

    // Moment de friction de roulement cf. Johnson
    if ( normwt > EPS )  delM -= k_m_s * normFN * wt / normwt;
//    else delM = 0.0;
  }
  else delM = 0.0;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul effectif des forces & moments de contact / cohesion
void CohContact::performForcesCalculus( Composant* p0_,
	Composant* p1_, const PointContact &contactInfos,
	Vecteur &delFN, Vecteur &delFT, Vecteur &delM )
{
  // Get the cohesive properties FmaxDist & Kelast
  Scalar fmax_dist = p0_->get_FmaxDist(p1_->getID());
  Scalar kn_elast = p0_->get_KnElast(p1_->getID());
  Scalar kt_elast = p0_->get_KtElast(p1_->getID());
  Scalar initial_overlap = p0_->get_InitialOverlap(p1_->getID());

  // Calcul des forces & moments de contact
  // --------------------------------------
  Point geometricPointOfContact = contactInfos.getContact();
  Vecteur penetration = contactInfos.getRecouvrement();
  Scalar distance = contactInfos.getDistance();

  // Vecteur normal au contact
  Vecteur normal( penetration );
  normal /= Norm( normal );
  normal.round();

  // Calcul de la vitesse relative au point de contact
  Vecteur tmpV = p0_->getVitesse( geometricPointOfContact )
    - p1_->getVitesse( geometricPointOfContact );

  Vecteur v_n  = normal * ( tmpV * normal );
  Vecteur v_t  = tmpV - v_n;

  // Calcul de la tangente portee par la vitesse relative tangente
  Scalar normv_t = Norm( v_t );
  Vecteur tangent(0.);
  if ( normv_t > EPS )
    tangent = v_t / normv_t;

  // Initialy in contact
  if ( Norm(initial_overlap-distance)<1.e-12 )
  {
    delFN = 0. * normal;
    delFT = 0. * tangent;
  }

  // ER Contact
  else if ( distance < initial_overlap )
  {
   /* cout << "  TEMP ER between P" << p0_->getID() << " & P" << p1_->getID()
    << endl; */

    // Force de restitution elastique lineaire
    delFN = stiff * (distance-initial_overlap) * -normal;

    // Force de restitution elastique normale + Force normale amortie
    double masse0 = p0_->getMasse();
    double masse1 = p1_->getMasse();
    Scalar avmass = masse0 * masse1 / ( masse0 + masse1 );
    Scalar omega0 = sqrt( stiff / avmass );
    if ( avmass == 0. )
    {
      avmass = masse1 == 0. ? 0.5 * masse0 : 0.5 * masse1;
      omega0 = sqrt( 2. * stiff / avmass );
    }
    Scalar muen = - omega0 * log(en) / sqrt( PI * PI + log(en) * log(en) );
    delFN += - muen * 2.0 * avmass * v_n;

    Scalar normFN = Norm( delFN );

    // Force tangentielle de dissipation
    delFT = v_t * ( -muet * 2.0 * avmass );

    // Force tangentielle par le critere de Coulomb
    Scalar fn = muec * normFN;
    Scalar ft = Norm( delFT );
    if ( fn < ft ) delFT = tangent * (-fn);

    // Moment de friction de roulement
    if ( k_m_s )
    {
	  // Calcul de la vitesse de rotation relative
	  Vecteur w = *p0_->getVitesseRotation() - *p1_->getVitesseRotation();
	  Vecteur wn = ( w * normal ) * normal;
	  Vecteur wt = w - wn ;
	  Scalar normwt = Norm( wt );

	  // Correction suivant la normale : effet anti-toupie
	  delM = - k_m_s * normFN * 0.001 * wn ;

	  // Moment de friction de roulement cf. Johnson
	  if ( normwt > EPS )  delM -= k_m_s * normFN * wt / normwt;
	  //    else delM = 0.0;
    }
    else delM = 0.0;
  }

  // Elastic like behavior zone for cohesive force
  else if (distance>initial_overlap && distance<fmax_dist-initial_overlap
           && fmax_dist>0.)
  {
    /*cout << "  TEMP COH ELAST between P" << p0_->getID() << " & P"
      << p1_->getID() << endl; */

    // Force cohesive normale
    delFN = kn_elast * (distance-initial_overlap)* -normal;

    // Force cohesive tangentielle
    delFT = kt_elast * (distance-initial_overlap) * -tangent;

  }

  // Plastic like behavior zone for cohesive force
  else if (distance>fmax_dist-initial_overlap
           && distance<rupt_dist-initial_overlap && fmax_dist>0.)
  {
   /* cout << "  TEMP COH PLAST between P" << p0_->getID() << " & P"
   << p1_->getID() << endl;*/

    Scalar kn_plast = -stiff*(initialFmax_dist-initial_overlap)
                     /(rupt_dist-initialFmax_dist-2.*initial_overlap);
    Scalar bn_plast = stiff*(initialFmax_dist-initial_overlap)
                     * (rupt_dist-initial_overlap)
                     /(rupt_dist-initialFmax_dist-2.*initial_overlap);

    Scalar kt_plast = -stiff_t*(initialFmax_dist-initial_overlap)
                     /(rupt_dist-initialFmax_dist-2.*initial_overlap);
    Scalar bt_plast = stiff_t*(initialFmax_dist-initial_overlap)
                     * (rupt_dist-initial_overlap)
                     /(rupt_dist-initialFmax_dist-2.*initial_overlap);

    // Force cohesive normale
    delFN = (kn_plast*(distance-initial_overlap)+bn_plast)*-normal;
    kn_elast = Norm(delFN)/(distance-initial_overlap);

    // Force cohesive tangentielle
    delFT = (kt_plast*(distance-initial_overlap)+bt_plast)*-tangent;
    kt_elast = Norm(delFT)/(distance-initial_overlap);

    // Endommagement
    p0_->set_CohesiveProperties(p1_->getID(), distance, kn_elast, kt_elast );
    p1_->set_CohesiveProperties(p0_->getID(), distance, kn_elast, kt_elast );
  }

  // Out of range, no forces applied
  else
  {
   /* cout << "  TEMP NO INTERACTION between P" << p0_->getID() << " & P"
   << p1_->getID() << endl; */

    p0_->set_CohesiveProperties(p1_->getID(), 0., 0.,0.);
    p1_->set_CohesiveProperties(p0_->getID(), 0., 0.,0.);
    delFN = 0. * normal;
    delFT = 0. * tangent;
  }
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul des forces & moments de contact
void CohContact::computeForcesPostProcessing( Composant* p0_,
	Composant* p1_, Scalar dt,
	PointContact &contactInfos,
	list<PointForcePostProcessing>* listOfContacts )
{
  Composant* ref_p0_ = p0_->ReferenceComposant() ;
  Composant* ref_p1_ = p1_->ReferenceComposant() ;
  Vecteur delFN, delFT, delM;
  Scalar diss;
  Point geometricPointOfContact = contactInfos.getContact();

  // Calcul des forces & moments de contact
  if ( Grains_Exec::m_ContactDissipation )
    performForcesCalculusPP( ref_p0_, ref_p1_, contactInfos, delFN, delFT,
    delM, diss );
  else  performForcesCalculus( ref_p0_, ref_p1_, contactInfos, delFN, delFT,
    delM);



  // Ajout � la liste des contacts pour post processing
  struct PointForcePostProcessing pfpp;
  pfpp.geometricPointOfContact = geometricPointOfContact;
  pfpp.contactForce = delFN + delFT;
  if ( Grains_Exec::m_ContactDissipation )   pfpp.contactDissip = diss;
  listOfContacts->push_back(pfpp);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul des forces & moments de contact pour un clone periodique
// D. RAKOTONIRINA - Juil 2014 - Modification
bool CohContact::computeForces( Composant* p0_,
	Composant* p1_,
	PointContact &contactInfos,
	LinkedCell *LC,
	Scalar dt, int nbContact )
{
  // cout << " TEMP ENTERING CohContact::computeForces " << endl;
  bool compute = false ;

  // Pour les particules composite, on calcule la moyenne des forces appliquees
  // a chaque point de contact. Par defaut nbContact = 1 (particules convexes)
  double coef = 1./double(nbContact);

  // Pour les particules composites, on renvoie la reference du composite
  // et non celle des particules elementaires
  Composant* ref_p0_ = p0_->ReferenceComposant() ;
  Composant* ref_p1_ = p1_->ReferenceComposant() ;

  // Si le contact concerne 2 clones periodiques et que l'angle des directions
  // periodiques est different de 90,
  // la force n'est pas prise en compte car les particules maitres se touchent
  // �galement et donc le contact est d�j� calcul� une fois
  if( !is_ClonePer_ClonePer( ref_p0_ ,ref_p1_ ) ) compute = true ;
  else if( are_ClonePerDirections_Perp( ref_p0_ ,ref_p1_ ) ) compute = true ;

  if ( compute  )
  {
    Vecteur delFN, delFT, delM;
    Point geometricPointOfContact = contactInfos.getContact();

    // Calcul des forces & moments de contact
    performForcesCalculus( ref_p0_, ref_p1_, contactInfos, delFN, delFT, delM );

    // Composant p0_
    if ( ref_p0_->getID() != -2 )
    {
      // Ajout aux torseurs d'effort au composant p0_
      p0_->addForce( geometricPointOfContact, coef * (delFN + delFT) );
      if ( k_m_s ) p0_->addMoment( delM * coef );
    }
    else
    {
      // Si p0_ est un clone periodique en contact avec une
      // particule p1_ possedant un clone periodique dans la
      // meme direction, alors la particule maitre de p0_ et le clone de
      // p1_ se touchent �galement et donc le contact est d�j� calcul� une fois
      if ( !is_ClonePer_ParticuleWithClonePerSameDirection( ref_p0_,
	  ref_p1_ , LC ) )
      {
        // Ajout aux torseurs d'effort au composant p0_
        p0_->addForce( geometricPointOfContact, (delFN + delFT) * coef );
        if ( k_m_s ) p0_->addMoment( delM * coef );
      }
    }

    // Composant p1_
    if ( ref_p1_->getID() != -2 )
    {
      // Ajout aux torseurs d'effort au composant p1_
      p1_->addForce( geometricPointOfContact, coef * (-delFN - delFT) );
      if ( k_m_s ) p1_->addMoment( - delM * coef );
    }
    else
    {
      // Si p1_ est un clone periodique en contact avec une
      // particule p0_ possedant un clone periodique dans la
      // meme direction, alors la particule maitre de p1_ et le clone de
      // p0_ se touchent �galement et donc le contact est d�j� calcul� une fois
      if( !is_ClonePer_ParticuleWithClonePerSameDirection(
          ref_p1_, ref_p0_, LC ) )
      {
        // Ajout aux torseurs d'effort au composant p1_
        p1_->addForce( geometricPointOfContact, coef * (-delFN - delFT) );
        if ( k_m_s ) p1_->addMoment( - delM * coef );
      }
    }
  }

  return compute ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialization of cohesive properties for particles initially in contact
bool CohContact::InitializeCohesiveObjects( Composant* p0_, Composant* p1_,
    const PointContact &contactInfos, LinkedCell *LC, Scalar dt, int nbContact )
{
  bool compute = false ;

  // Pour les particules composites, on renvoie la reference du composite
  // et non celle des particules elementaires
  Composant* ref_p0_ = p0_->ReferenceComposant() ;
  Composant* ref_p1_ = p1_->ReferenceComposant() ;

  // Si le contact concerne 2 clones periodiques et que l'angle des directions
  // periodiques est different de 90,
  // la force n'est pas prise en compte car les particules maitres se touchent
  // �galement et donc le contact est d�j� calcul� une fois
  if( !is_ClonePer_ClonePer( ref_p0_ ,ref_p1_ ) ) compute = true ;
  else if( are_ClonePerDirections_Perp( ref_p0_ ,ref_p1_ ) ) compute = true ;

  Scalar distance = contactInfos.getDistance();

  if( compute && distance < 0. )
  {
    p0_->initializeCohesiveProperties(
        p1_->getID(), initialFmax_dist, stiff, stiff_t, distance );
    p1_->initializeCohesiveProperties(
        p0_->getID(), initialFmax_dist, stiff, stiff_t, distance );
  }

  return compute ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Define the parameters values for CohContact
map<string,double> CohContact::defineParameters( DOMNode* &root )
{
  map<string,double> parameters;
  DOMNode* parameter;
  double   value;
  parameter = ReaderXML::getNode(root, "stiff");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["stiff"]  = value;
  parameter = ReaderXML::getNode(root, "stiff_t");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["stiff_t"]  = value;
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
  parameter = ReaderXML::getNode(root, "color");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["color"]  = value;
  parameter = ReaderXML::getNode(root, "ruptdist");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["ruptdist"]  = value;
  parameter = ReaderXML::getNode(root, "initialmaxdist");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["initialmaxdist"]  = value;
  return parameters;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul une estimation du temps de contact et de la penetration
// maximale d'un choc elastique binaire avec dissipation d'energie, puis ecrit
// le resultat dans un flux de sortie
void CohContact::computeAndWriteEstimates( Composant* p0_,
	Composant* p1_,
  	const double &v0,
	ostream &OUT ) const
{
  double masse0 = p0_->getMasse();
  double masse1 = p1_->getMasse();
  Scalar avmass = masse0 * masse1 / (masse0 + masse1);
  if (avmass == 0.) avmass = masse1==0. ? 0.5*masse0 : 0.5*masse1;
  Scalar muen = - sqrt(stiff/avmass) * log(en) / sqrt( PI*PI+log(en)*log(en) );

  // Choc particule/obstacle
  if ( p0_->isObstacle() || p1_->isObstacle() )
  {
    Composant *particule=NULL,*obstacle=NULL;
    if ( p0_->isObstacle() )
    {
      particule = p1_;
      obstacle = p0_ ;
    }
    else
    {
      particule = p0_;
      obstacle = p1_ ;
    }
    masse0 = particule->getMasse();
    double delta_allowed = p0_->getRayonInteraction()
  	+ p1_->getRayonInteraction();
    double omega0 = sqrt(stiff/masse0);
    double theta = sqrt(pow(omega0,2.) - pow(0.5*muen,2.));
    double tc = PI / theta;

    double delta_max = computeDeltaMax( theta, 0.5*muen, en, tc, v0 );

    OUT << "  Particule: materiau = " << particule->materiau()
  	<< " croute = " << particule->getRayonInteraction()
	<< " poids = " << masse0*9.81 << endl;
    OUT << "  Obstacle: materiau = " << obstacle->materiau()
  	<< " croute = " << obstacle->getRayonInteraction()
	<< " poids = " << masse0*9.81 << endl;
    OUT << "  Recouvrement autorise par les croutes = "
    	<< delta_allowed << endl;
    OUT << "  Vitesse relative imposee = " << v0 << endl;
    OUT << "  Tc = " << tc << "  delta_max = " << delta_max << endl;
    OUT << "  mu_n = " << muen << endl;
    OUT << "  Force elastique max = " << stiff*delta_allowed << endl;
    OUT << "  Rapport fel/poids0 = " << stiff*delta_allowed/(masse0*9.81)
  	<< endl;
    if ( delta_max > delta_allowed )
    {
      OUT << "  *************************************************" << endl;
      OUT << "  !!!!! WARNING !!!!!" << endl;
      OUT << "  delta_max > recouvrement autorise par les croutes" << endl;
      OUT << "  *************************************************" << endl;
    }
  }
  // Choc particule/particule
  else
  {
    double delta_allowed = p0_->getRayonInteraction()
  	+ p1_->getRayonInteraction();
    double omega0 = sqrt(stiff/avmass);
    double theta = sqrt(pow(omega0,2.) - pow(muen,2.));
    double tc = PI / theta;

    double delta_max = computeDeltaMax( theta, muen, en, tc, v0 );

    OUT << "  Particule 0: materiau = " << p0_->materiau()
  	<< " croute = " << p0_->getRayonInteraction()
	<< " poids = " << masse0*9.81 << endl;
    OUT << "  Particule 1: materiau = " << p1_->materiau()
  	<< " croute = " << p1_->getRayonInteraction()
	<< " poids = " << masse1*9.81 << endl;
    OUT << "  Recouvrement autorise par les croutes = " << delta_allowed
    	<< endl;
    OUT << "  Vitesse relative imposee = " << v0 << endl;
    OUT << "  Tc = " << tc << "  delta_max = " << delta_max << endl;
    OUT << "  mu_n = " << muen << endl;
    OUT << "  Force elastique max = " << stiff*delta_allowed << endl;
    OUT << "  Rapport fel/poids0 = " << stiff*delta_allowed/(masse0*9.81)
  	<< endl;
    OUT << "  Rapport fel/poids1 = " << stiff*delta_allowed/(masse1*9.81)
  	<< endl;
    if ( delta_max > delta_allowed )
    {
      OUT << "  *************************************************" << endl;
      OUT << "  !!!!! WARNING !!!!!" << endl;
      OUT << "  delta_max > recouvrement autorise par les croutes" << endl;
      OUT << "  *************************************************" << endl;
    }
  }
  OUT << endl;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul de la penetration maximale sur la base de la solution
// analytique et d'un algo de Newton
double CohContact::computeDeltaMax( const double &theta_, const double& mu_,
  	const double& en_, const double &tc_, const double &v0_ ) const
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

  return delta_max;
}
