#include "MPIWrapperGrains.hh"
#include "ERHContact.H"
#include "Grains_Exec.hh"
#include "Composant.H"
#include "Memento.hh"
#include "LinkedCell.H"

#ifndef NB_COULOMB_REGIME
#define NB_COULOMB_REGIME 0
#endif
#ifndef NB_STATIC_REGIME
#define NB_STATIC_REGIME 0
#endif

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur avec un lot de parametres
ERHContact::ERHContact( map<string,double>& parameters ) :
  ContactLaw()
{
  stiff = parameters["stiff"];
  en    = parameters["en"   ];
  muet  = parameters["mut"  ];
  muec  = parameters["muc"  ];
  ks    = parameters["ks"   ];
  k_m_s = parameters["kms"  ];
  eta_r = parameters["nr"   ];
  mu_r  = parameters["mur"  ];
  Jn    = parameters["Jn"   ];
  m_f     = parameters["f"  ];
  epsilon = parameters["eps"];

  if ( (!k_m_s) || (!eta_r) || (!mu_r) || (!Jn) || (!m_f) )
  {
    rolling_friction = true;
  }
  else rolling_friction = false;
}




// ----------------------------------------------------------------------------
// @brief Destructeur
ERHContact::~ERHContact()
{
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul effectif des forces & moments de contact
// It also sends back the dissipated energy for pp purposes
void ERHContact::performForcesCalculusPP( Composant* p0_,
	Composant* p1_, PointContact &contactInfos,
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
// Calcul effectif des forces & moments de contact
void ERHContact::performForcesCalculus( Composant* p0_,
	Composant* p1_, Scalar dt, PointContact &contactInfos,
	Vecteur &delFN, Vecteur &delFT, Vecteur &delM)
{
  int *pnbCumulTangent = NULL;
  Vecteur *pvectorTangentFriction = NULL;
  double *ptangentialDepl = NULL;
  Vecteur *pspringRotFriction = NULL;
  double normFTHistory = 0.;
  Vecteur frictionTorque = 0.;
  double tangdepl;
  Vecteur vectorTF;
  Vecteur w ;
  Vecteur wn ;
  Vecteur wt ;
  double radius ;
  double radius1 ;
  double Req = 0. ;
  double kr ;
  double Cr ;
  double max_normFtorque ;

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
    Req = p0_->getRayon();
    omega0 = sqrt( 2. * stiff / avmass );
  }
  Scalar muen = - omega0 * log(en) / sqrt( PI * PI + log(en) * log(en) );
  delFN += - muen * 2.0 * avmass * v_n;
  Scalar normFN = Norm( delFN );

  // Force tangentielle de dissipation
  delFT = v_t * ( -muet * 2.0 * avmass );

  if ( rolling_friction )
  {
    radius = p0_->getRayon() ;
    if (!Req)
    {
      radius1 = p1_->getRayon() ;
      Req = radius * radius1 / (radius + radius1) ;
    }
    //Req to be changed to be able to handle particle-particle contacts
    kr = 2*Jn*Req*normFN ;
    Cr = 2*eta_r*sqrt(3./4 * avmass *kr) * radius;
    // cout << "mur=" << mu_r << ", Req=" << Req << ", normFN=" << normFN << endl;
    max_normFtorque = mu_r * Req * normFN ;
    // cout << "Max torque = " << max_normFtorque << endl ;
    // Calcul de la vitesse de rotation relative
    w = *p0_->getVitesseRotation() - *p1_->getVitesseRotation();
    wn = ( w * normal ) * normal;
    wt = w - wn ;
  }

  // Deplacement cumulatif au point de contact
  tangdepl = Norm(v_t) * dt;
  vectorTF=tangent;

  if ( p0_->ContactInMapIsActive( p1_->getID(), pnbCumulTangent,
    pvectorTangentFriction, ptangentialDepl, pspringRotFriction) )
  {
    // Friction in translation
    *ptangentialDepl += tangdepl;
    if (normv_t>epsilon){
      *pvectorTangentFriction=tangent;
      *pnbCumulTangent = 0.;
    }
    else{
      if (normv_t>EPS){
        Scalar scaleNbCumulTangent=*pnbCumulTangent;
        vectorTF=*pvectorTangentFriction;
        if (vectorTF*v_t<0.){
          //making sure the tangential vector points along -v_t
          scaleNbCumulTangent*=-1;
        }
        *pvectorTangentFriction *= (scaleNbCumulTangent);
        *pvectorTangentFriction+=tangent;
        // Making sure there is no normal component on the tangential vector
        *pvectorTangentFriction-=(*pvectorTangentFriction*normal)*normal;
        *pvectorTangentFriction/=Norm(*pvectorTangentFriction);
        (*pnbCumulTangent)+=1;
      }
      // Nothing to do if normv_t<EPS
    }
    vectorTF=*pvectorTangentFriction;
    normFTHistory = ks * (*ptangentialDepl);

    // Friction in rotation
    if (rolling_friction)
    {
      *pspringRotFriction += (-kr) * dt * wt;
      double normSpringRotFriction = Norm(*pspringRotFriction);
      if (normSpringRotFriction > max_normFtorque)
      {
        *pspringRotFriction *= max_normFtorque / normSpringRotFriction;
        p1_->addDeplContactInMap( p0_->getID(), *pnbCumulTangent,
          *pvectorTangentFriction, *ptangentialDepl, *pspringRotFriction );
        frictionTorque = *pspringRotFriction ;
        frictionTorque += (-m_f) * Cr * wt;
      }
      else
      {
        p1_->addDeplContactInMap( p0_->getID(), *pnbCumulTangent,
          *pvectorTangentFriction, *ptangentialDepl, *pspringRotFriction );
        frictionTorque = *pspringRotFriction ;
        frictionTorque += (-Cr) * wt;
      }
    }
    else
    {
      p1_->addDeplContactInMap( p0_->getID(), *pnbCumulTangent,
        *pvectorTangentFriction, *ptangentialDepl, *pspringRotFriction );
    }



  }
  else
  {
    normFTHistory = ks * tangdepl;
    if (rolling_friction)
    {
      frictionTorque = (-kr) * dt * wt;
      double normSpringRotFriction = Norm(frictionTorque);
      // cout << endl;
      // cout << "------------------------------------------------" << endl ;
      // cout << "Contact has been interrupted" << endl ;
      // cout << "------------------------------------------------" << endl ;
      // cout << endl;
      if (normSpringRotFriction > max_normFtorque)
      {
        // cout << "Saturation of the friction torque from 1st time step" << endl;
        frictionTorque *= max_normFtorque / normSpringRotFriction;
        p1_->addNewContactInMap( p0_->getID(), 0, vectorTF, tangdepl,
          frictionTorque);
        p0_->addNewContactInMap( p1_->getID(), 0, vectorTF, tangdepl,
          frictionTorque);
        frictionTorque += (-m_f) * Cr * wt;
      }
      else
      {
        p1_->addNewContactInMap( p0_->getID(), 0, vectorTF, tangdepl,
          frictionTorque);
        p0_->addNewContactInMap( p1_->getID(), 0, vectorTF, tangdepl,
          frictionTorque);
        frictionTorque += (-Cr) * wt;
      }
    }
  }

  // Friction in translation
  delFT += vectorTF * (-normFTHistory);
  Scalar fn = muec * normFN;
  Scalar ft = Norm( delFT );
  if ( fn < ft )
  {
    delFT = vectorTF * (-fn);
  }

  // Friction in rotation
  if (rolling_friction)
  {
    // Correction suivant la normale : effet anti-toupie
    if (k_m_s) delM = - k_m_s * normFN * 0.001 * wn ;
    else delM = - kr * 0.001 * wn ;


    // Correction with a spring-dashpot model with backwards rolling (see Ai et
    // al., Powder Technology, 2011)
    // cout << "friction torque = " << frictionTorque << endl;
    // delM += frictionTorque;

    // MODEL A : TO REMOVE (just here to reproduce Ai et al. graphs)
    double normwt=Norm(wt);
    if (normwt > EPS) delM += -mu_r * Req * normFN * wt / normwt;

    // // Routine to store friction torque in file
    // ofstream frictionTorqueFile ;
    // frictionTorqueFile.open("Grains/Init/frictionTorque.txt", std::ios_base::app) ;
    // frictionTorqueFile << delM ;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul des forces & moments de contact
void ERHContact::computeForcesPostProcessing( Composant* p0_,
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
  else  performForcesCalculus( ref_p0_, ref_p1_, dt, contactInfos, delFN, delFT,
    delM);

  // Ajout � la liste des contacts pour post processing
  struct PointForcePostProcessing pfpp;
  pfpp.geometricPointOfContact = geometricPointOfContact;
  pfpp.contactForce = delFN + delFT;

  // Ajout des pointeurs des composants en contact pour le reseau de forces
  pfpp.comp0 = ref_p0_;
  pfpp.comp1 = ref_p1_;

  if ( Grains_Exec::m_ContactDissipation )   pfpp.contactDissip = diss;
  listOfContacts->push_back(pfpp);

  // Composant p0_
  // Ajout aux torseurs d'effort au composant p0_
  if ( Grains_Exec::m_ContactforceOutput_instantaneous )
    p0_->addContactForcePP_instantaneous( delFN + delFT );
  // Should add domain limits
  if ( Grains_Exec::m_stressTensor ||
    Grains_Exec::m_particleStressTensor )
      p0_->computeInternalMoments( geometricPointOfContact,
          delFN + delFT );

  // Composant p1_
  // Ajout aux torseurs d'effort au composant p1_
  if ( Grains_Exec::m_ContactforceOutput_instantaneous )
    p1_->addContactForcePP_instantaneous(  -delFN - delFT );
  if ( ( Grains_Exec::m_stressTensor ||
    Grains_Exec::m_particleStressTensor ) && !p1_->isObstacle() )
      p1_->computeInternalMoments( geometricPointOfContact,
        -delFN - delFT );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul des forces & moments de contact pour un clone periodique
// D. RAKOTONIRINA - Juil 2014 - Modification
bool ERHContact::computeForces( Composant* p0_,
	Composant* p1_,
	PointContact &contactInfos,
	LinkedCell *LC,
	Scalar dt, int nbContact )
{
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
  if ( !is_ClonePer_ClonePer( ref_p0_ ,ref_p1_ ) ) compute = true ;
  else if ( are_ClonePerDirections_Perp( ref_p0_ ,ref_p1_ ) ) compute = true ;


  if ( compute  )
  {
    Vecteur delFN, delFT, delM;
    Point geometricPointOfContact = contactInfos.getContact();

    // Calcul des forces & moments de contact
    performForcesCalculus( ref_p0_, ref_p1_, dt, contactInfos, delFN, delFT, delM );

    // Composant p0_
    if ( ref_p0_->getID() != -2 )
    {
      p0_->addForce( geometricPointOfContact,
	      coef * (delFN + delFT) );
      // Ajout aux torseurs d'effort au composant p0_
      if ( Grains_Exec::m_ContactforceOutput )
	p0_->addContactForcePP( coef * (delFN + delFT) );
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
	if ( Grains_Exec::m_ContactforceOutput )
	  p0_->addContactForcePP(  coef * (delFN + delFT) );
        if ( k_m_s ) p0_->addMoment( delM * coef );
      }
    }

    // Composant p1_
    if ( ref_p1_->getID() != -2 )
    {
      // Ajout aux torseurs d'effort au composant p1_
      p1_->addForce( geometricPointOfContact, coef * (-delFN - delFT) );
      if ( Grains_Exec::m_ContactforceOutput )
	p1_->addContactForcePP(  coef * (-delFN - delFT) );
      if ( k_m_s ) p1_->addMoment( - delM * coef );
    }
    else
    {
      // Si p1_ est un clone periodique en contact avec une
      // particule p0_ possedant un clone periodique dans la
      // meme direction, alors la particule maitre de p1_ et le clone de
      // p0_ se touchent �galement et donc le contact est d�j� calcul� une fois
      if ( !is_ClonePer_ParticuleWithClonePerSameDirection( ref_p1_, ref_p0_, LC ) )
      {
        // Ajout aux torseurs d'effort au composant p1_
        p1_->addForce( geometricPointOfContact, coef * (-delFN - delFT) );
	if ( Grains_Exec::m_ContactforceOutput )
	  p1_->addContactForcePP(  coef * (-delFN - delFT) );
        if ( k_m_s ) p1_->addMoment( - delM * coef );
      }
    }
  }

  return compute ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Define the parameters values for ERHContact
map<string,double> ERHContact::defineParameters( DOMNode* &root )
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
  if (parameter)
  {
    value     = ReaderXML::getNodeValue_Double(parameter);
    parameters["kms"]    = value;
  }
  else parameters["kms"] = 0;
  parameter = ReaderXML::getNode(root, "nr");
  if (parameter)
  {
    value     = ReaderXML::getNodeValue_Double(parameter);
    parameters["nr"]    = value;
  }
  else parameters["nr"] = 0;
  parameter = ReaderXML::getNode(root, "mur");
  if (parameter)
  {
    value     = ReaderXML::getNodeValue_Double(parameter);
    parameters["mur"]    = value;
  }
  else parameters["mur"] = 0;
  parameter = ReaderXML::getNode(root, "Jn");
  if (parameter)
  {
    value     = ReaderXML::getNodeValue_Double(parameter);
    parameters["Jn"]    = value;
  }
  else parameters["Jn"] = 0;
  parameter = ReaderXML::getNode(root, "f");
  if (parameter)
  {
    value     = ReaderXML::getNodeValue_Double(parameter);
    parameters["f"]    = value;
  }
  else parameters["f"] = 0;
  parameter = ReaderXML::getNode(root, "ks");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["ks"]    = value;
  parameter = ReaderXML::getNode(root, "eps");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["eps"]  = value;
  parameter = ReaderXML::getNode(root, "color");
  value     = ReaderXML::getNodeValue_Double(parameter);
  parameters["color"]  = value;

  return parameters;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcul une estimation du temps de contact et de la penetration
// maximale d'un choc elastique binaire avec dissipation d'energie, puis ecrit
// le resultat dans un flux de sortie
void ERHContact::computeAndWriteEstimates( Composant* p0_,
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
double ERHContact::computeDeltaMax( const double &theta_, const double& mu_,
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
