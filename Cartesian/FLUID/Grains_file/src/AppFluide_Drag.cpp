#include "Grains_Exec.hh"
#include "AppFluide_Drag.H"
#include "Torseur.H"

// ----------------------------------------------------------------------------
// Constructeur
// M.BERNARD - Janv.2012 - Creation
AppFluide_Drag::AppFluide_Drag() :
  App()
{
}




// ----------------------------------------------------------------------------
// Constructeur
AppFluide_Drag::AppFluide_Drag( istream &fileIn ) :
  App()
{
  read( fileIn );
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
AppFluide_Drag::AppFluide_Drag( DOMNode* root ) :
  App(),
  b_withLiftForces( Grains_Exec::m_withLiftForce ),
  b_withPressureGradient( Grains_Exec::m_withPressureGradient ),
  b_withAddedMass(Grains_Exec::m_addedmass_demcfd )
{
  type_Fhydro = ReaderXML::getNodeAttr_String( root, "Type" );
  if(Grains_Exec::m_withStochasticDrag)
    b_fluctdrag = true;
  else b_fluctdrag = false;
  meanD = EnsComposant::get_SauterMeanDiameter( );
  sum_XY = EnsComposant::get_sumXY( );

}




// ----------------------------------------------------------------------------
// Destructeur
AppFluide_Drag::~AppFluide_Drag()
{
}




// ----------------------------------------------------------------------------
// Description de la force de AppFluide_Drag exercee sur les particules
// M.BERNARD - Janv.2012 - Creation
void AppFluide_Drag::CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules )
{
  double  rhoP=0., mass=0., rp=0., dp=0., particleVolume=0.;
  double  Norm_Vr=0., Rep = 0. ,Cd=0., fd0=0., alpha=0. ;
  double  Ef = 0., chi = 0., normalised_dragForce=0.,Cd_phi = 0.,Cd_Re=0. ;
  double  b=0., f_Re=0., PHI=0., Cd_ReT=0., a_c_fl=0.;
  double  w=0., F0=0., F1=0., F2=0., F3=0. ;
  double  A=0., B=0., Beta=0., Vp=0. , Beta_prime=0.;
  Vecteur dragForce, hydroForce, Vr, stokes_einstein_force, dragForce_fluct;
  Vecteur const* translationVelocity = NULL;
  Scalar rhoF = Particule::getFluideMasseVolumique();
  Scalar muF = Particule::getFluidViscosity();
  Vecteur const* translationVelocity_fluid = NULL;
  Vecteur const* gradp_fluid = NULL;
  Vecteur pressureTerm;
  list<Particule*>::const_iterator il;
  double * GC = new double[3];
  double K0=0., K1=0., K2=0., K3=0., ReSauter=0.;
  Vecteur SaffmanLiftForce, MagnusLiftForce, hydroTorque, AddedMass;
  Vecteur gaussian_rand(0.,0.,0.);
  double sigma;

  for (il=particules->begin() ; il!=particules->end() ; il++)
  {
      // Caracteristiques par particule (uniquement en 3D)
      rhoP = (*il)->getMasseVolumique();
      mass = (*il)->getMasse();
      rp = (*il)->getRayonSphereEquivalente();
  //    rp = (*il)->getRayon();

      dp = 2. * rp;
      Vp = acos(-1)*pow(dp,3)/6. ;
      translationVelocity = (*il)->getVitesseTranslation();
      particleVolume =  mass/rhoP;
      translationVelocity_fluid = (*il)->getVitesseTranslation_fluide();
      Ef = (*il)->get_DEMCFD_volumeFraction() ;
      if ( Ef < 0.25 )
    	  {
	   // cout << " !!ATTENTION!! Saturating the solid volume fraction to eps\n"
	   // << " = 0.75 Real Es = " << 1.-Ef << endl;
	    Ef = 0.25;
	    (*il)->getPosition(GC);
	   // cout << " GC = " << GC[0] << " " << GC[1]<< " "<<  GC[2] << endl;
	  }
      PHI = 1.-Ef;
      Vr.setValue(
		  (*translationVelocity_fluid)[X] - (*translationVelocity)[X],
	 	  (*translationVelocity_fluid)[Y] - (*translationVelocity)[Y],
	 	  (*translationVelocity_fluid)[Z] - (*translationVelocity)[Z]);

      Norm_Vr = Norm(Vr);
      Rep = dp * rhoF * Ef * Norm_Vr / muF;
      Vecteur const* rotU = NULL;
      Vecteur const* rotationVelocity = NULL;
      Vecteur omegaUf, omegaRel;
      double Rer=0.;
      // Added mass force computation
      if ( b_withAddedMass )
      {
        AddedMass = 0.5*Vp*rhoF*(Vr -
	  (*il)->getTranslationalSlipVelocityPreviousTime())/dt;
        (*il)->setVelocitySlipPreviousTime(Vr);
      }
      else
        AddedMass.setValue(0., 0., 0.);
      if( b_withLiftForces )
      {
	rotU = (*il)->getVorticity_fluide();
	rotationVelocity = (*il)->getVitesseRotation();
	omegaUf = 0.5*(*rotU);
	omegaRel = omegaUf - *rotationVelocity;
	Rer = rhoF*pow(dp,2)*Norm(omegaRel)/muF;

	double SaffmanLiftCoef = 0., alpha_Mei=0.;
	double MagnusLiftCoef = 0.;
	double rotationCoef = 0.;

	if( Norm(omegaUf) < 1e-20 ) // because rotU = 0 au premier pas de temps.
	{
          SaffmanLiftForce.setValue(0., 0., 0.);
          MagnusLiftForce.setValue(0., 0., 0.);
	}
	else
	{
          SaffmanLiftForce.setValue(
                    Vr[1]*omegaUf[Z]-Vr[2]*omegaUf[Y],
                    Vr[2]*omegaUf[X]-Vr[0]*omegaUf[Z],
                    Vr[0]*omegaUf[Y]-Vr[1]*omegaUf[X] );
          alpha_Mei = rp*Norm(omegaUf)/Norm_Vr;
          if( Rep > 40 ) SaffmanLiftCoef = 0.0524 * pow(alpha_Mei*Rep,0.5);
          else SaffmanLiftCoef = (1-0.3314)*pow(alpha_Mei,0.5)*exp(-Rep/10) +
              0.3314*pow(alpha_Mei,0.5);
          SaffmanLiftForce = 1.615*4.*pow(rp,2)*pow(rhoF*muF,0.5)*
        	SaffmanLiftCoef * SaffmanLiftForce / pow(Norm(omegaUf),0.5);

          // ========================
          // ======== MAGNUS ========
          MagnusLiftForce.setValue(
                    omegaRel[1]*Vr[2]-omegaRel[2]*Vr[1],
                    omegaRel[2]*Vr[0]-omegaRel[0]*Vr[2],
                    omegaRel[0]*Vr[1]-omegaRel[1]*Vr[0] );
              // On utilise la formulation de Oesterle
              //  - valide pour Rep<140 !!!
              //  - ~= Rer/Rep pour Rep<1
              // Formulation de  Lun&Liu1997 douteuse
              // ne fit pas du tout avec l'autre
              // MagnusLiftCoef = 2*rp*Norm(omegaRel)/Norm_Vr*
              //    	(0.178+0.822*pow(Rep,-0.522));
          MagnusLiftCoef = 0.45+(Rer/Rep-0.45)*
            exp(-0.05684*pow(Rer,0.4)*pow(Rep,0.3));
          MagnusLiftForce = 0.5*acos(-1)*pow(rp,2)*rhoF*Norm_Vr*
            MagnusLiftCoef * MagnusLiftForce / Norm(omegaRel);
	}

	// ==============================
	// ======== HYDRO TORQUE ========
	if( Rer < 1.e-20 )
          rotationCoef = 0.;
	else if( Rer < 32 )
          rotationCoef = 64*acos(-1)/Rer;
	else
          rotationCoef = 12.9/sqrt(Rer)+128.4/Rer;
	hydroTorque = 0.5*rhoF*pow(rp,5)*rotationCoef*
      		  Norm(omegaRel)*omegaRel;
      }
      else
      {
	SaffmanLiftForce.setValue(0., 0., 0.);
	MagnusLiftForce.setValue(0., 0., 0.);
	hydroTorque.setValue(0., 0., 0.);
      }

      //////////////////////////////////////////////////////////////
      if( type_Fhydro == "None" )
      {
	dragForce.setValue(0.,0.,0.) ;
      }

      //////////////////////////////////////////////////////////////
      else if( type_Fhydro == "Stokes" )
      {
	dragForce = 3.* PI * muF * dp * Vr ;
      }

      //////////////////////////////////////////////////////////////
      else if( type_Fhydro == "Stokes_WithUp" )
      {
	dragForce = - 3.* PI * muF * dp * (*translationVelocity) ;
      }

      //////////////////////////////////////////////////////////////
      // Di Felice 1994. Plutot pour du liquide/solid
      else if ( type_Fhydro == "DiFelice" )
      {
	// Cd = pow((0.63 + 4.8*pow(Rep,-0.5)),2);

	// Drag Force from Schiller-Naumann
	if(Rep < 1.e-12) Cd = (24/1.e-12);
	else if(Rep < 1000.) Cd = (24/Rep)*(1+0.15*pow(Rep,0.687));
	else            Cd = 0.44;

	chi = 3.7 - 0.65 * exp(-0.5*pow((1.5-log10(Rep)),2));
	fd0 =  0.5 * Cd * rhoF * PI * pow(rp, 2) *
      	   pow( Ef, 2) * Norm_Vr ;
	dragForce = fd0 * Vr * pow(Ef, -(chi+2.));
	// Following 1 line is only for PostProcessing purpose
      //  Beta = Cd * pow(Ef, -(chi+2.));
      }

      //////////////////////////////////////////////////////////////
      // ATTENTION: You need to be cautious about the range of dimless
      // numbers: 4.5<Re<7 ==> for Liquid/Solid
      //          16.54 <Re<28 ==> for Gas/Solid
      //In either case you should comment/uncomment below
      //
        // Di Felice 1994. Plutot pour du liquide/solid
      else if ( type_Fhydro == "DiFelice_Fixed_PlusFluc" )
      {
        double beta=-1;
        double alphaFd=0.260;
       	double randNum=0.,Fdev=0.;
        Vecteur dragForce0;
	Vecteur GRN(0.,0.,0.);

	// Cd = pow((0.63 + 4.8*pow(Rep,-0.5)),2);

	// Drag Force from Schiller-Naumann
	if(Rep < 1.e-12) Cd = (24/1.e-12);
	else if(Rep < 1000.) Cd = (24/Rep)*(1+0.15*pow(Rep,0.687));
	else            Cd = 0.44;

	chi = 3.7 - 0.65 * exp(-0.5*pow((1.5-log10(Rep)),2));
	fd0 =  0.5 * Cd * rhoF * PI * pow(rp, 2) *
      	   pow( Ef, 2) * Norm_Vr ;
	dragForce0 = fd0 * Vr * pow(Ef, -(chi+2.));

        if ( time == dt  )
        {
         randNum=((double) rand() / (RAND_MAX));
        (*il)->set_rnd(randNum,0,0);
        }
        else
        {
         GRN = *((*il)->get_rnd());
         randNum = GRN[0];
        }
         // Compute the inverse function
         double sgn,a,x,term1,term2,pi,erf_inv;
         pi=3.1415926535897931;
         a=0.147;
         x=2*randNum-1.;
         sgn = (x < 0) ? -1. : 1.;
         term1=2/(pi*a);
         term2=log(1.-x*x);
         erf_inv=sgn*pow(-term1-term2/2.+pow(pow(term1+term2/2.,2)-term2/a,0.5),0.5);
         // Compute the Nusselt number
         Fdev=beta+exp(erf_inv*pow(2.,0.5)*alphaFd);
         dragForce = dragForce0*(1.+Fdev); //      dragForce = fd0 * Vr * pow(Ef, -(chi));
      }

      else if ( type_Fhydro == "DiFelice_fluct" )
      {
	// Drag Force from Schiller-Naumann
	if(Rep < 1.e-12) Cd = (24/1.e-12);
	else if(Rep < 1000.) Cd = (24/Rep)*(1+0.15*pow(Rep,0.687));
	else            Cd = 0.44;

	chi = 3.7 - 0.65 * exp(-0.5*pow((1.5-log10(Rep)),2));
	fd0 =  0.5 * Cd * rhoF * PI * pow(rp, 2) *
      	   pow( Ef, 2) * Norm_Vr ;
	dragForce = fd0 * Vr * pow(Ef, -(chi+2.));
        Beta = Cd * pow(Ef, -(chi+2.));
          if ( int(time/dt)%int(b_fluidTimeStep/dt)==0 )
	  {
	 // PRS data are given here: 4.5<Re<7 ==> for Liquid/Solid
	 //                          16.54 <Re<28 ==> for Gas/Solid
	 //   if (Rep > 7.)
	 //     sigma = 13.;
	 //   else if (Rep <4.5)
	 //     sigma = 52.;
         //   else
	 //   sigma = (Rep-4.5)*(13.-52.)/(7.5-4.5) + 52. ;
	    if (Rep > 28.)
	      sigma = 16.54;
	    else if (Rep <16.54)
	      sigma =  43.9;
            else
	      sigma = (Rep-16.54)*(16.-43.9)/(28.-16.54) + 43.9;

	    gaussian_rand =*((*il)->get_rnd());
	    gaussian_rand[0] = gaussian_rand[0]*a_c_fl +
	        pow(1.-pow(a_c_fl,2),0.5) * gaussian_rand[1];
	    gaussian_rand[1] = rand_normal(0.,1.);
            if ( gaussian_rand[1]>2. )
                gaussian_rand[1] = 2.;
            if ( gaussian_rand[1]<-2. )
                gaussian_rand[1] = -2.;
	    gaussian_rand[1] = gaussian_rand[1]*sigma;
	    (*il)->set_rnd(gaussian_rand[0],gaussian_rand[1],
	        gaussian_rand[2]);
	  }
	  else gaussian_rand = *((*il)->get_rnd());
          Beta_prime = gaussian_rand[0]*a_c_fl +
	        pow(1.-pow(a_c_fl,2),0.5) * gaussian_rand[1];
          if ( fabs(Beta_prime) > fabs(Beta) )
            Beta_prime = Beta * Beta_prime / fabs(Beta_prime);
	  dragForce_fluct = rhoF*PI*rp*rp*Vr*Norm(Vr)*Ef*Ef*Beta_prime/2.;
      }


     //////////////////////////////////////////////////////////////
      // From Beetstra2007
      else if ( type_Fhydro == "bidisperse_Beetstra" )
      {
	A = (10*PHI)/pow(Ef,2) +
      		  (pow(Ef, 2))*(1+1.5*sqrt(PHI));

	B = (0.413*Rep)/(24*pow(Ef,2)) *
      		  (pow(Ef,-1) + 3*PHI*Ef+8.4*pow(Rep,-0.343))/
		  (1+pow(10,(3*PHI))*pow(Rep,-0.5-2*PHI)) ;

	Beta = 18*muF*Ef*PHI*(A+B)/(dp*dp);

	dragForce = particleVolume*Beta*Vr/(Ef*PHI);

	double yi = dp / meanD;
	dragForce = ( Ef*yi + PHI*pow(yi,2) + 0.064*Ef*pow(yi,3) )
      						  * dragForce;
      }

      //////////////////////////////////////////////////////////////
      // From Cello & DiRenzo 2010
      else if ( type_Fhydro == "bidisperse_DiRenzo" )
      {
	ReSauter = meanD * rhoF * Ef * Norm_Vr / muF;
	stokes_einstein_force = 3*PI * muF * dp * Ef * Vr ;

	K0 = PHI/(1+3*Ef);
	K1 = (1+128*K0+715*K0*K0)/(Ef*Ef*(1+49.5*K0));
	K2 = (1+0.13*ReSauter*6.66e-4*ReSauter*ReSauter)/
      		  (1+3.42e-2*ReSauter+6.92e-6*ReSauter*ReSauter)-1;
	K3 = (2*ReSauter*ReSauter/(1+ReSauter))*
      	  ((-410*Ef+9.2e7*ReSauter*pow(K0,20)+1.9e3*Ef*Ef-6.6e-2*ReSauter)/
	  (6.6e3*Ef+4.92e-4*ReSauter-4.3e4*Ef*Ef-1.31e-4*ReSauter*ReSauter+
		  7.38e4*pow(Ef,3)));
	fd0 = K1 + K2*pow(Ef,4) + K3*(1-pow(Ef,4));

	double yi = dp / meanD;

	Beta = yi+(PHI/Ef)*((PHI-0.27)/(1-0.27))*((yi*yi-yi)/sum_XY);

  //      dragForce = fd0 * stokes_einstein_force / Ef;
	dragForce = fd0 * stokes_einstein_force;

      }

      //////////////////////////////////////////////////////////////
      // Ergun pour milieu dense + Wen & Yu pour milieu dilue
      else if (type_Fhydro == "Gidaspow")
      {
	if(Ef<=0.8) // ERGUN
	{
          // Muller2009 -> *Ef - Model A
          // Deb2014 -> *Ef - Model B (inspire de Muller2009)
	  // Sakai2014 -> *Ef - Model A
	  // Deen2007 -> Gauche *Ef ; Droite *Ef**2 - Model A
	  // Jajcevic2013 -> *Ef - Model A
          Beta = ( 150. * pow(PHI,2) * muF ) /
		( pow(Ef,2) * pow(dp,2)) +
		(1.75 * rhoF * PHI * Norm_Vr) /
		(Ef * dp) ;
	}
	else // WEN & YU (Ef > 0,8)
	{
          if(Rep > 1000.)  Cd = 0.44 ;
          else             Cd = 24.*(1.+0.15*pow(Rep,0.687))/Rep ;

          // Same as Muller2009 - Model A
          //         Deb2014 - Model B (pourquoi Ergun est *E ? ->rho*g ?)
          //         Sakai2014 - Model A
          Beta = 0.75* Cd * rhoF * Ef * PHI * Norm_Vr *
		  pow(Ef,-2.65) / dp ;
	}
	dragForce = particleVolume*Beta*Vr/PHI ;
      }

      //////////////////////////////////////////////////////////////
      //  Wen&Yu
      //
      else if (type_Fhydro == "WenAndYu")
      // modification due to gradP
      {
         double beta_WenYu=0.;
         if(Rep > 1000.)  Cd = 0.44 ;
         else             Cd = 24.*(1.+0.15*pow(Rep,0.687))/Rep ;

         beta_WenYu = 0.75* Cd * rhoF * PHI * Ef * Norm_Vr *
                            pow(Ef,-2.7) / dp ;

         dragForce = particleVolume*beta_WenYu*Vr/(Ef*PHI);
    //      dragForce = particleVolume*Beta*Vr/PHI;
      }
      //
      //////////////////////////////////////////////////////////////
      // Ergun + Wen&Yu + smoothing function
      else if (type_Fhydro == "Huilin")
      // modification due to gradP
      {
	double beta_Ergun=0., beta_WenYu=0., phi_Huilin=0.;

	beta_Ergun =
            ( 150.*pow(PHI,2)*muF ) / ( Ef*dp*dp ) +
	    ( 1.75*rhoF*PHI*Norm_Vr ) / dp ;

	if(Rep > 1000.)  Cd = 0.44 ;
	else             Cd = 24.*(1.+0.15*pow(Rep,0.687))/Rep ;

	beta_WenYu = 0.75* Cd * rhoF * PHI * Ef * Norm_Vr *
		  pow(Ef,-2.7) / dp ;

	phi_Huilin = (atan(150*1.75*(0.2-PHI)))/(acos(-1)) + 0.5;
	Beta = (1-phi_Huilin)*beta_Ergun + phi_Huilin*beta_WenYu;

	dragForce = particleVolume*Beta*Vr/(Ef*PHI);
  //      dragForce = particleVolume*Beta*Vr/PHI;
      }

      //////////////////////////////////////////////////////////////
      // BEESTRA 2007
      // /E one more time because we use model B instead of model A
      // modification due to dradP
      else if (type_Fhydro == "Beetstra_2007")
      {
	stokes_einstein_force = 3*PI * muF * dp * Ef * Vr ;

	fd0 = (10*PHI)/pow(Ef,2) +
      		  (pow(Ef, 2))*(1+1.5*sqrt(PHI));

	alpha = (0.413*Rep)/(24*pow(Ef,2)) *
      		  (pow(Ef,-1) + 3*PHI*Ef+8.4*pow(Rep,-0.343))/
		  (1+pow(10,(3*PHI))*pow(Rep,-0.5-2*PHI)) ;

  //      dragForce = (fd0 + alpha) * stokes_einstein_force / Ef;
	dragForce = (fd0 + alpha) * stokes_einstein_force;

      }

      //////////////////////////////////////////////////////////////
      // BEESTRA 2007
      // /E one more time because we use model B instead of model A
      // modification due to dradP
      else if (type_Fhydro == "myBeetstra")
      {
	A = (10*PHI)/pow(Ef,2) +
      		  (pow(Ef, 2))*(1+1.5*sqrt(PHI));

	B = (0.413*Rep)/(24*pow(Ef,2)) *
     		  (pow(Ef,-1) + 3*PHI*Ef+8.4*pow(Rep,-0.343))/
		  (1+pow(10,(3*PHI))*pow(Rep,-0.5-2*PHI)) ;
	Beta = 18*muF*Ef*PHI*(A+B)/(dp*dp);
	dragForce = particleVolume*Beta*Vr/(Ef*PHI);

//	 Following 4 lines are only for PostProcessing purpose
//	if (Rep > 0.01)
//	  Beta = (A+B)*24./(Rep*Ef);
//	else
//	  Beta = (A+B)*24./(0.01*Ef);
      }
      //////////////////////////////////////////////////////////////
      // ATTENTION: You need to be cautious about the range of dimless
      // numbers: 3.<Re<10.5 ==> for Liquid/Solid
      //          16.54 <Re<28 ==> for Gas/Solid
      //In either case you should comment/uncomment below

      else if (type_Fhydro == "myBeetstra_fluct")
      {
	A = (10*PHI)/pow(Ef,2) +
      		  (pow(Ef, 2))*(1+1.5*sqrt(PHI));

	B = (0.413*Rep)/(24*pow(Ef,2)) *
      		  (pow(Ef,-1) + 3*PHI*Ef+8.4*pow(Rep,-0.343))/
		  (1+pow(10,(3*PHI))*pow(Rep,-0.5-2*PHI)) ;
	if (Rep > 0.001)
	{
          Beta = (A+B)*24./(Rep*Ef);
          dragForce = rhoF*PI*rp*rp*Vr*Norm(Vr)*Ef*Ef*Beta/2.;
//        We do not want to add randomness during solid timestepping
          if ( int(time/dt)%int(b_fluidTimeStep/dt)==0 )
	  {
	// PRS data are given here: 4.5<Re<7 ==> for Liquid/Solid
	//                          16.54 <Re<28 ==> for Gas/Solid
	//    if (Rep > 28.)
	//      sigma = 16.54;
	//    else if (Rep <16.54)
	//      sigma =  43.9;
        //    else
	//      sigma = (Rep-16.54)*(16.-43.9)/(28.-16.54) + 43.9;

        // 2017: extended liquid/solid with a new fit:
        // cd_prime = ax**b; a = 512.875, b = -1.54919
	    if (Rep <3.)
            {
              sigma = 93.5;
              //tau = 0.0210687 ;
              tau = 0.030058;
            }
            else if (Rep > 10.5)
            {
              sigma = 13.5;
              //tau = 0.01137045;
              tau = 0.0094056835972865;
            }
            else
            {
              sigma = 512.875 * pow(Rep,-1.54919);
              //tau = -0.0012931*Rep + 0.0249483;
              tau = 0.0832629*pow(Rep,-0.92741);
            }

            a_c_fl = exp(-b_fluidTimeStep/tau);
	    gaussian_rand =*((*il)->get_rnd());
	    gaussian_rand[0] = gaussian_rand[0]*a_c_fl +
	        pow(1.-pow(a_c_fl,2),0.5) * gaussian_rand[1];
	    gaussian_rand[1] = rand_normal(0.,1.);
            if ( gaussian_rand[1]>2. )
                gaussian_rand[1] = 2.;
            if ( gaussian_rand[1]<-2. )
                gaussian_rand[1] = -2.;
	    gaussian_rand[1] = gaussian_rand[1]*sigma;
            gaussian_rand[2] = a_c_fl;
	    (*il)->set_rnd(gaussian_rand[0],gaussian_rand[1],
	        gaussian_rand[2]);
	  }
	  else gaussian_rand = *((*il)->get_rnd());
          Beta_prime = gaussian_rand[0]*gaussian_rand[2] +
	        pow(1.-pow(gaussian_rand[2],2),0.5) * gaussian_rand[1];
          if ( fabs(Beta_prime) > fabs(Beta) )
            Beta_prime = Beta * Beta_prime / fabs(Beta_prime);
	  dragForce_fluct = rhoF*PI*rp*rp*Vr*Norm(Vr)*Ef*Ef*Beta_prime/2.;
	}
	else
	{
	  Beta = 18*muF*Ef*PHI*(A+B)/(dp*dp);
	  dragForce = particleVolume*Beta*Vr/(Ef*PHI);
	  // Following 4 lines are only for PostProcessing purpose
	  if (Rep > 0.01)
	    Beta = (A+B)*24./(Rep*Ef);
	  else
	    Beta = (A+B)*24./(0.01*Ef);
	}
      }

      else if (type_Fhydro == "Beestra_from_PepiotDesjardins")
      {
	// Probleme de Rep/Ef (peut etre un oubli de leur par)
	Rep = Rep/Ef;
	stokes_einstein_force = 3.*PI * muF * dp * Ef * Vr ;

	A = (10*PHI)/pow(Ef,2) +
      		  (pow(Ef, 2))*(1+1.5*sqrt(PHI));

	B = (0.413*Rep)/(24*pow(Ef,2)) *
      		  (pow(Ef,-1) + 3*PHI*Ef+8.4*pow(Rep,-0.343))/
		  (1+pow(10,(3*PHI))*pow(Rep,-0.5-2*PHI)) ;

	dragForce = (A + B) * stokes_einstein_force;
      }

      //////////////////////////////////////////////////////////////
      // Beestra from Muller2009 /!\ formulation A
      else if (type_Fhydro == "Beestra_from_Muller")
      {
	A=180.+(18.*pow(Ef,4)*(1.+1.5*sqrt(PHI)))/PHI;
	B=(0.31*(pow(Ef,-1)+3.*PHI*Ef+8.4*pow(Rep,-0.343)))/
      	  (1.+pow(10.,3.*PHI)*pow(Rep,(2.*Ef-2.5)));

	Beta = A*muF*pow(1.-Ef,2)/(dp*dp*Ef) + B*muF*PHI*Rep/(dp*dp);

	dragForce = Vp*Beta*Vr/PHI ;
      }


      //////////////////////////////////////////////////////////////
      // cf Beestra_2007
      else if (type_Fhydro == "Ergun_Carman")
      {
	stokes_einstein_force = 6. * muF * rp * PI * Vr ;
	fd0 = (10*PHI)/pow(Ef, 2) +
      		  pow(Ef, 2)*(1.+1.5*sqrt(PHI)) ;
	b = 180.;
	alpha = b / (18.*pow(Ef,2)) ;

	normalised_dragForce = fd0 + alpha ;

	dragForce = normalised_dragForce * stokes_einstein_force ;
      }

      //////////////////////////////////////////////////////////////
      // cf Beestra_2007
      else if (type_Fhydro == "Hill_Beestra")
      {
  // F_E_0_HILL    = 10.*(1.-E)/E**2 + E**2 * (1.+1.5*sqrt(1.-E))
  // ALPHA_E_HILL(x) = 0.03365*E + 0.106*E*(1.-E) + 0.0116/E**4 #+(6.*(1.-E)-10.*(1.-E)**2)/(Re(x)*E**2)
  // F_normalised_HILL(x)  = F_E_0_HILL + ALPHA_E_HILL(x) * Re(x)
  // Fd_HILL(x)   = F_normalised_HILL(x) * F_stokes_einstein(x)

	stokes_einstein_force = 6. * muF * rp * PI * Vr ;
	fd0 = 10.*PHI/pow(Ef,2) +
      		  pow(Ef, 2)*(1.+1.5*sqrt(PHI)) ;

	alpha = 0.03365*Ef + 0.106*Ef*PHI +
      		  0.0116/pow(Ef,4) +
		  (6.*PHI - 10.*pow(1.-Ef,2))/(Rep*pow(Ef,2)) ;

	normalised_dragForce = fd0 + alpha ;

	dragForce = normalised_dragForce * stokes_einstein_force ;
      }

      //////////////////////////////////////////////////////////////
      // cf Alobaid 2012 --> probleme avec B, 0.6057*E**2
      else if (type_Fhydro == "Alobaid") //Koch et al 2001 moderate Re
      {
	if ( Ef < 0.6 )
          A = 180. ;
	else
          A = (18.*pow(Ef,3))/PHI *
	    (1.+3.*sqrt(PHI/2)+(135.*PHI*log(PHI))/64. +
	    16.14*PHI ) / (1.+0.681*PHI-8.48*pow(PHI,2)+
	    8.16*pow(PHI,3) );

	B = 0.6057*pow(PHI,2) + 1.908*PHI*pow(Ef,2)+
          0.209*pow(Ef,-3) ;

	Beta = (muF/ pow(dp,2)) *
            (A*pow(PHI,2)/Ef +
	    B*PHI*Rep ) ;

	dragForce = particleVolume*Beta*Vr / (Ef*PHI) ;
      }


      //////////////////////////////////////////////////////////////
      // cf VanDerHoef2006 and link2009comparison
      else if (type_Fhydro == "Hill")
      {
	if ( Ef < 0.6 )
          A = 180. ;
	else
          A = (18.*pow(Ef,3))/PHI *
	    (1.+3.*sqrt(PHI/2)+(135.*PHI*log(PHI))/64. +
	    16.14*PHI ) / (1.+0.681*PHI-8.48*pow(PHI,2)+
	    8.16*pow(PHI,3) );

	B = 0.6057*pow(Ef,2) + 1.908*PHI*pow(Ef,2)+
          0.209*pow(Ef,-3) ;

	Beta = (muF/ pow(dp,2)) *
            (A*pow(PHI,2)/Ef +
	    B*PHI*Rep ) ;

	dragForce = particleVolume*Beta*Vr / (Ef*PHI) ;
      }


      //////////////////////////////////////////////////////////////
      // cf Beestra_2007
      else if (type_Fhydro == "Gibilaro")
      {
  // F_Re_Gibilaro(x)          = 17.3/18. + 0.336*Re(x) / 18.
  // F_normalised_Gibilaro(x)  = F_Re_Gibilaro(x) * E**(-3.8)
  // Fd_Gibilaro(x)            = F_normalised_Gibilaro(x) * F_stokes_einstein(x)

	stokes_einstein_force = 6. * muF * rp * PI * Vr ;
	f_Re = 17.3/18. + 0.336*Rep/18. ;

	normalised_dragForce = f_Re * pow(Ef,-3.8) ;

	dragForce = normalised_dragForce * stokes_einstein_force ;
      }

      //////////////////////////////////////////////////////////////
      else if ( type_Fhydro == "Benyahia" ) // Voir Benyahia2006
      {
	Rep=0.5*Rep;
	stokes_einstein_force = 6. * muF * rp * PI * Vr ;

	w  = exp(-10.*(0.4-PHI)/PHI) ;

	if( (Ef>0.6) && (Ef>0.99) )
          F0 = (1-w)*( (1. + 3.*sqrt(PHI/2.) + (135./64.)*
		    PHI*log(PHI) + 17.14*PHI ) /
		    ( 1.+0.681*PHI - 8.48*pow(1-Ef,2) +
		    8.16*pow(1-Ef,3) ) ) +
		    w*(10.*PHI/pow(Ef,3)) ;
	else if ( Ef <= 0.6 )
          F0 = 10.*PHI/pow(Ef,3) ;
	else
          F0 = 0. ;

	if( Ef>0.9 && Ef<0.99 )
          F1 = (sqrt(2/PHI)) / 40. ;
	else if ( Ef < 0.9 )
          F1 = 0.11 + 0.00051*exp(11.6*PHI) ;

	if( Ef > 0.6 )
          F2 = (1-w)*( (1+3*sqrt(PHI/2) +
	    (135./64.)*PHI*log(PHI) + 17.89*(1-Ef) ) /
	    ( 1.+0.681*(1-Ef)-11.03*pow(1-Ef,2) +
	    15.41*pow(1-Ef,3) ) ) +
	    w*(10.*PHI/pow(Ef,3)) ;
	else
          F2 = 10.*(1-Ef)/pow(Ef,3) ;

	if( Ef > 0.9047 )
          F3 = 0.9351*(1-Ef) + 0.03667 ;
	else
          F3 = 0.0673 + 0.212*(1-Ef) + 0.0232/pow(Ef,5) ;

	if( (Ef > 0.99) && (Rep<=((F2-1.)/(3./8.-F3))) )
          normalised_dragForce = 1. + 3./8. * Rep ;
	else if( (Ef<0.99) &&
      		  (Rep<=((F3+sqrt(pow(F3,2)-4.*F1*(F0-F2)))/(2.*F1))) )
          normalised_dragForce = F0 + F1*pow(Rep,2) ;
	else if( ( (Ef<0.99) &&
      		  (Rep>((F3+sqrt(pow(F3,2)-4.*F1*(F0-F2)))/(2.*F1))) ) ||
		  (Ef>=0.99 && Rep>((F2-1.)/(3./8.-F3))) )
          normalised_dragForce = F2 + (F3)*Rep ;
	else
	{
          cout <<"ERROR while computing Benyahia drag force"<<endl<<
		 "Rep     = " << Rep <<endl<<
		 "Ef = " << Ef << endl;
          exit(0) ;
	}
	dragForce = normalised_dragForce * stokes_einstein_force ;
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Fhydro == "Tenneti_2011" )
      {
       	// Drag Force from Schiller-Naumann (in a stokes-einstein-
	// based drag formulation)
	if(Rep < 1.e-12) Cd = 1.;
	else if(Rep < 1000.) Cd = 1+0.15*pow(Rep,0.687);
	else            Cd = 0.44*Rep/24.;

	stokes_einstein_force = 3*PI * muF * dp * Ef * Vr ;

	Cd_phi = (5.81*PHI)/pow(Ef,3) + 0.48*pow(PHI,1./3.)/pow(Ef,4);
	Cd_Re = pow(PHI,3)*Rep*(0.95+ 0.61*pow(PHI,3)/pow(Ef,2) );
	normalised_dragForce = Cd/pow(Ef,3) + Cd_Re + Cd_phi;
	dragForce = normalised_dragForce * stokes_einstein_force ;
	// Following 4 lines are only for PostProcessing purpose
	// We do not use Eps in the denominator because "normalised_dragForce
	// \equiv (A+B)/eps". compare with the definition of drag in myBeetstra
	if (Rep > 0.01)
	  Beta = (normalised_dragForce)*24./(Rep);
	else
	  Beta = (normalised_dragForce)*24./(0.01);
      }
      //////////////////////////////////////////////////////////////
      else if (type_Fhydro == "Bogner_2015")
      {
        stokes_einstein_force = 3*PI * muF * dp * Ef * Vr ;
        Cd = pow(Ef,-5.726)*(1.751 + 0.151*pow(Rep,0.684)
             -0.445*pow(1.+Rep,1.04*PHI) - 0.16*pow(1.+Rep,0.0003*PHI));
        dragForce = Cd * stokes_einstein_force;
      }
      //////////////////////////////////////////////////////////////
      else if (type_Fhydro == "Schiller-Naumann" )
      {
        // Drag Force from Schiller-Naumann
	if(Rep < 1.e-12) Cd = (24/1.e-12);
	else if(Rep < 1000.) Cd = (24/Rep)*(1+0.15*pow(Rep,0.687));
	else            Cd = 0.44;
	fd0 =  0.5 * Cd * rhoF * PI * pow(rp, 2) *
      	   pow( Ef, 2) * Norm_Vr ;
	dragForce = fd0 * Vr;
      }
      //////////////////////////////////////////////////////////////
      else if (type_Fhydro == "Tang_2015")
      {
	stokes_einstein_force = 3*PI * muF * dp * Ef * Vr ;
	Cd_phi = (10.*PHI)/pow(Ef,2.) + pow(Ef,2)*(1.+1.5*pow(PHI,0.5) ) ;
	Cd_Re = 0.11*PHI*(1.+PHI) - 0.00456/pow(Ef,4.) +
		(0.169*Ef + 0.0644/pow(Ef,4))*pow(Rep,-0.343) ;
	normalised_dragForce = Cd_phi + Cd_Re * Rep ;
	dragForce = stokes_einstein_force * normalised_dragForce /Ef;
      }
      //////////////////////////////////////////////////////////////
      else if (type_Fhydro == "Tang_2016_mobile")
      {
	stokes_einstein_force = 3*PI * muF * dp * Ef * Vr ;
	Cd_phi = (10.*PHI)/pow(Ef,2) + pow(Ef,2)*(1.+1.5*pow(PHI,0.5) ) ;
	Cd_Re = 0.11*PHI*(1.+PHI) - 0.00456/pow(Ef,4.) +
		(0.169*Ef + 0.0644/pow(Ef,4))*pow(Rep,-0.343) ;
	Cd_ReT = 2.98 * (2.108*pow(Rep,0.85)*pow(rhoP/rhoF,-0.5)) *
		PHI/pow(Ef,2) ;

	normalised_dragForce = Cd_phi + Cd_Re * Rep + Cd_ReT;

	dragForce = stokes_einstein_force * normalised_dragForce /Ef;
      }
      else
      {
	cout <<"ERROR with type_Fhydro"<<endl;
	exit(0) ;
      }

      // If we add the pressure term to the total hydroForce,
      // then we multiply the dragForce by Ef !! (model A)
      if( b_withPressureGradient )
      {
	gradp_fluid = (*il)->getGradientPression_fluide();
	pressureTerm = - particleVolume*(*gradp_fluid) ;
	double ratio = fabs(Norm(pressureTerm)/Norm(dragForce));
	if (ratio >10000000.)
	    pressureTerm = 0.;
        hydroForce = Ef * dragForce + SaffmanLiftForce + MagnusLiftForce
		+AddedMass + pressureTerm;
      }
      else
      {
	if (b_fluctdrag)
	  hydroForce = dragForce + SaffmanLiftForce + MagnusLiftForce+
	  AddedMass + dragForce_fluct ;
        else
	  hydroForce = dragForce + SaffmanLiftForce + MagnusLiftForce+
	               AddedMass;
      }

      // -----------------------------------------------------------
      // Ajout de la force au torseur des efforts dans Composant.cpp
      (*il)->addBodyForce( hydroForce ) ;

      // ---------------------------------------------------------
      // Ajout du moment au torseur des efforts dans Composant.cpp
      (*il)->addMoment( hydroTorque ) ;

      // -----------------------------------------------------------
      // Affectation a la structure m_fluidInfos dans Particule.cpp
      (*il)->setHydroForce_PP( &hydroForce ) ;
      if (b_fluctdrag)
      {
        hydroForce = hydroForce - dragForce_fluct;
        (*il)->setHydroForce( &hydroForce ) ;
      }
      else  (*il)->setHydroForce( &hydroForce ) ;

    // -----------------------------------------------------------
      // Affectation a la structure m_fluidInfos dans Particule.cpp
      if( b_slipveloutput )    (*il)->setSlipVel(&Vr) ;

    } // end if particles are not in halozone
//  }//  end loop on particles

  delete[] GC;
}




// ----------------------------------------------------------------------------
// Lecture des constantes de AppFluide_Drag
void AppFluide_Drag::read( istream &fileIn )
{
  fileIn >> type_Fhydro ;
}




// ----------------------------------------------------------------------------
// Set if slip velocity is on output
void AppFluide_Drag::set_slipveloutput( bool slipveloutput )
{
  b_slipveloutput = slipveloutput;
}




//----------------------------------------------------------------------------
// Generate random number using Box Muller distribution
double AppFluide_Drag::rand_normal(double mean, double stddev)
{
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}




// ----------------------------------------------------------------------------
// Set fluid timestep
void AppFluide_Drag::set_simultime( double dt )
{
  b_fluidTimeStep = dt;
}
