#include "Grains_Exec.hh"
#include "AppFluide_Temperature.H"
//#include "Torseur.H"

double AppFluide_Temperature::m_heatCapacityS = 0.;

// ----------------------------------------------------------------------------
// Constructeur
// M.BERNARD - July.2015 - Creation
// A. HAMMOUTI - mars 2016 - Modification
AppFluide_Temperature::AppFluide_Temperature() :
  App()
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
AppFluide_Temperature::AppFluide_Temperature( DOMNode* nTemperature ) :
  App()
{
  string interSolid="";

  DOMNode* nSolidTemp = ReaderXML::getNode( nTemperature, "Solid" );
  if( nSolidTemp )
  {
    // ========================================================================
    if( ReaderXML::hasNodeAttr_Double( nSolidTemp, "initTempS" ) )
      m_initialTempS = ReaderXML::getNodeAttr_Double( nSolidTemp, "initTempS" );
    else
      cout << "initTempS attribute requested in Solid Temperature module"
           << endl;

    // ========================================================================
    if( ReaderXML::hasNodeAttr_Double( nSolidTemp, "CpS" ) )
      m_heatCapacityS = ReaderXML::getNodeAttr_Double( nSolidTemp, "CpS" );
    else
      cout << "CpS attribute requested in Solid Temperature module"
           << endl;

    // ========================================================================
    if( ReaderXML::hasNodeAttr_String( nSolidTemp, "interSolid" ) )
    {
      interSolid = ReaderXML::getNodeAttr_String( nSolidTemp, "interSolid" );
      if( interSolid == "yes" )
      {
        Grains_Exec::m_withSolidTemperature = true;
	cout << "  Solid Temperature" << endl;

        if( ReaderXML::hasNodeAttr_Double( nSolidTemp, "thermalConductS" ) )
          m_diffCoefS = ReaderXML::getNodeAttr_Double( nSolidTemp, "thermalConductS" );
        else
          cout << "thermalConductS attribute requested in Solid Temperature module"
               << endl;
      }
    }
  }
  else
  {
    cout << " 'Solid' module requested in Temperature module !" << endl;
    exit(0);
  }

  DOMNode* nFluidTemp = ReaderXML::getNode( nTemperature, "Fluid" );
  if( nFluidTemp )
  {
    if( Grains_Exec::m_withHydroForce )
    {
      Grains_Exec::m_withFluidTemperature = true;
      cout << "  Fluid Temperature" << endl;
      // ========================================================================
      if( ReaderXML::hasNodeAttr_Double( nFluidTemp, "initTempF" ) )
      {
        if( !Grains_Exec::m_withdemcfd )
          m_initialTempF = ReaderXML::getNodeAttr_Double( nFluidTemp, "initTempF" );
        else
          cout << "WARNING : initial temperature written in simul.xml will" << endl
               << "          be neglected and the computed one will be" << endl
               << "          considered instead" << endl;
      }
      // ========================================================================
      if( ReaderXML::hasNodeAttr_Double( nFluidTemp, "thermalConductF" ) )
      {
        if( !Grains_Exec::m_withdemcfd )
          m_diffCoefF = ReaderXML::getNodeAttr_Double( nFluidTemp, "thermalConductF" );
        else
	  cout << "WARNING : Fluid thermal conductivity  writen in simul.xml"
               << endl
               << "          will be neclected and the computed one will be"
               << endl
               << "          considered instead" << endl;
      }
      // ========================================================================
      if( ReaderXML::hasNodeAttr_Double( nFluidTemp, "CpF" ) )
      {
        if( !Grains_Exec::m_withdemcfd )
          m_heatCapacityF = ReaderXML::getNodeAttr_Double( nFluidTemp, "CpF" );
        else
	  cout << "WARNING : Fluid heat capacity  writen in simul.xml"
               << endl
               << "          will be neclected and the computed one will be"
               << endl
               << "          considered instead" << endl;
      }
    }
    else
      cout << "WARNING : Fluid temperature coupling impossible without "
           << endl
           << "          HydroForce module, please add it to Forces module"
           << endl;
  }

  DOMNode* nNusselt = ReaderXML::getNode( nTemperature, "Nusselt" );
  if( nNusselt )
  {
    string isStochasticNusselt;
    if(  ReaderXML::hasNodeAttr_String( nNusselt, "Type" ) )
    {
      type_Nusselt  = ReaderXML::getNodeAttr_String( nNusselt, "Type" );
      cout << "TypeNusselt: " << type_Nusselt << endl;
      if ( type_Nusselt == "LocalAnalysis" )
      {
        if (ReaderXML::hasNodeAttr_Double( nNusselt, "Lbox" ))
	{
          m_cellDiamRatio = ReaderXML::getNodeAttr_Double( nNusselt, "Lbox" );
	}
        else
        {
          cout << "WARNING : Using 'LocalAnalysis' closure law without"
	  << endl
	  <<"                specifying cell to particle diameter ratio"
	  << endl
	  <<"                in simul.xml !!" << endl;
	  exit(0);
        }
      }
      else
      {
       m_cellDiamRatio = 10.;
      }
    }
    if(  ReaderXML::hasNodeAttr_String( nNusselt, "WithStochasticNusselt" ) )
    {
     isStochasticNusselt=ReaderXML::getNodeAttr_String( nNusselt, "WithStochasticNusselt");
    }
    if ( isStochasticNusselt == "yes"  )
    {
     Grains_Exec::m_withStochasticNusselt = true;
     b_fluctNu = true;
     cout << "With Stochastic Nusselt" << endl;
    if (type_Nusselt == "myGunn_fluct"|| type_Nusselt == "Sun_fluct")
    {
      tau = ReaderXML::getNodeAttr_Double( nNusselt, "Tau" );
      cout << "  Fluctuation time scale = "
	<<tau<<endl;
    }
    if (type_Nusselt == "myGunn_Dynamic_PlusFluc_LogNormal")
    {
      tau = ReaderXML::getNodeAttr_Double( nNusselt, "Tau" );
      cout << "  Fluctuation time scale = "
	<<tau<<endl;
    }
    }
    else
      cout << "WARNING : Specify Nusselt modele in simul.xml will" << endl;
  }
  else
  {
    cout << " 'Nusselt' module requested in Temperature module !" << endl;
    exit(0);
  }
  Particule::setFluidInitialTemperature( m_initialTempF ) ;
//  Particule::setSolidInitialTemperature( m_initialTempS ) ;

}




// ----------------------------------------------------------------------------
// Destructeur
AppFluide_Temperature::~AppFluide_Temperature()
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialization step
void AppFluide_Temperature::InitializeTemperature( Scalar time, Scalar dt,
    list<Particule*> const* particules,
    const Scalar fluidThermalConductivity_ )
{
//  cout << "TEMPORARY: in AppFluide_Temperature::InitializeTemperature" << endl;

  if( Grains_Exec::m_withdemcfd )
    m_diffCoefF = fluidThermalConductivity_;

  list<Particule*>::const_iterator il;
  for (il=particules->begin() ; il!=particules->end() ; il++)
  {
    (*il)->set_solidTemperature( m_initialTempS );
//    (*il)->set_heatCapacity( m_heatCapacityS );
    if( Grains_Exec::m_withFluidTemperature )
      (*il)->set_DEMCFD_fluidTemperature( m_initialTempF );
  }

  if( Grains_Exec::m_withFluidTemperature )
  {
    double muF = Particule::getFluidViscosity();
    m_Prandtl = m_heatCapacityF * muF / m_diffCoefF;
  }
}




// ----------------------------------------------------------------------------
// Description de la force de AppFluide_Drag exercee sur les particules
// M.BERNARD - June 2015 - Creation
void AppFluide_Temperature::CalculerForces( Scalar time, Scalar dt,
    list<Particule*> const* particules )
{
  if( Grains_Exec::m_withFluidTemperature )
  {
    Scalar dp=0., Rep=0.;
    Scalar Ef=0., Norm_Vr=0.;
    Scalar Nu=0., heatFlux=0.;
    Scalar Nulam=0., Nuturb=0., Nuiso=0., A=0., B=0., n=0.;
    double const* Tp;
	double const* Tf;
    Scalar rhoF = Particule::getFluideMasseVolumique();
    Scalar muF = Particule::getFluidViscosity();
    Vecteur const* translationVelocity_fluid = NULL;
    Vecteur const* translationVelocity = NULL;
    Vecteur Vr;
    double * GC = new double[3];

    list<Particule*>::const_iterator il;

    for (il=particules->begin() ; il!=particules->end() ; il++)
    {
      (*il)->getPosition(GC);
      dp = 2. * (*il)->getRayonSphereEquivalente();
      //mp = (*il)->getMasse();
      translationVelocity = (*il)->getVitesseTranslation();
      translationVelocity_fluid = (*il)->getVitesseTranslation_fluide();
      Ef = (*il)->get_DEMCFD_volumeFraction() ;
      Tf = (*il)->get_DEMCFD_fluidTemperature() ;

      Vr.setValue(
          (*translationVelocity_fluid)[X] - (*translationVelocity)[X],
          (*translationVelocity_fluid)[Y] - (*translationVelocity)[Y],
          (*translationVelocity_fluid)[Z] - (*translationVelocity)[Z]);

      Norm_Vr = Norm(Vr);
      Rep = dp * rhoF * Ef * Norm_Vr / muF;
      Tp = (*il)->get_solidTemperature();

      if( type_Nusselt == "None" )
      {
	Nu = 0.0;
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "Diffusion")
      {
	Nu = 2.0;
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "RanzMarshall")
      {
	Nu = 2. + 0.6 * pow(Rep,0.5)*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "Feng")
      {
	Nu = 0.992 + pow(Rep*m_Prandtl,1./3.)+0.1*pow(Rep,1./3.)*pow(Rep*m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "Gunn")
      {
	Nu = (7.-10.*Ef+5.*pow(Ef,2.))*(1.+0.7*pow(Ef*Rep,0.2)*pow(m_Prandtl,1./3.)) +
           (1.33-2.4*Ef+1.2*pow(Ef,2.))*pow(Ef*Rep,0.7)*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
       else if ( type_Nusselt == "myGunn_FixedBed_PlusFluc")
      {
        double beta=-1.;
        double alphaNu=0.336;
	double randNum=0.,Nu0=0.,Nudev=0.;
	Vecteur GRN(0.,0.,0.);

        Nu0 = (7.-10.*Ef+5.*pow(Ef,2.))*(1.+0.7*pow(Rep,0.2)*pow(m_Prandtl,1./3.))+(1.33-2.4*Ef+1.2*pow(Ef,2.))*pow(Rep,0.7)*pow(m_Prandtl,1./3.);

        if ( time == dt  )
        {
         randNum=((double) rand() / (RAND_MAX));
        (*il)->set_rnd_Nu(randNum,0.,0.);
        }
        else
        {
         GRN = *((*il)->get_rnd_Nu());
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
         Nudev=beta+exp(erf_inv*pow(2.,0.5)*alphaNu);
         Nu = Nu0*(1.+Nudev);
      }
      //////////////////////////////////////////////////////////////
       else if ( type_Nusselt == "myGunn_Dynamic_PlusFluc_WhiteNoise")
      {
	double Nu0=0.,Nu_prime=0.;
        double sigma=0.,a_c_fl=0.;
	Vecteur gaussian_rand(1.,1.,0.);

        Nu0 = (7.-10.*Ef+5.*pow(Ef,2.))*(1.+0.7*pow(Rep,0.2)*pow(m_Prandtl,1./3.))+(1.33-2.4*Ef+1.2*pow(Ef,2.))*pow(Rep,0.7)*pow(m_Prandtl,1./3.);

          if ( int(time/dt)%int(b_fluidTimeStep/dt)==0 )
          {
           // PRS data are given here: 4.5<Re<7 ==> for Liquid/Solid
           // PRS data are given here: 4.5<Re<7 ==> for Liquid/Solid
            if (Rep<3.)
            {
              tau=0.02403;
            }
            else if (Rep>10.5)
            {
              tau=0.0144;
            }
            else
            {
             tau=0.047*pow(Rep,-0.513);
            }
            a_c_fl=exp(-b_fluidTimeStep/tau);
            gaussian_rand =*((*il)->get_rnd_Nu());
            gaussian_rand[0] = gaussian_rand[0]*a_c_fl +
                pow(1.-pow(a_c_fl,2),0.5) * gaussian_rand[1];
            gaussian_rand[1] = rand_normal(0.,1.);
            if ( gaussian_rand[1]>2. )
                gaussian_rand[1] = 2.;
            if ( gaussian_rand[1]<-2. )
                gaussian_rand[1] = -2.;
            gaussian_rand[1] = gaussian_rand[1]*sigma;
            gaussian_rand[2] = a_c_fl;
            (*il)->set_rnd_Nu(gaussian_rand[0],gaussian_rand[1],
                gaussian_rand[2]);
          }
          else gaussian_rand = *((*il)->get_rnd_Nu());
          Nu_prime = gaussian_rand[0]*gaussian_rand[2]+
                pow(1.-pow(gaussian_rand[2],2),0.5) * gaussian_rand[1];
          Nu = Nu0+Nu_prime;

      }
     //////////////////////////////////////////////////////////////
       else if ( type_Nusselt == "myGunn_Dynamic_PlusFluc_LogNormal")
      {
        double beta=-1.;
        double alphaNu=0.336;
	double Nu0=0.,Nudev=1.,Nu_prime=1.,mX=0.,sX=0.;
        double a_c_fl=0.,a_c_fl_grn=0.;
	Vecteur gaussian_rand(1.,1.,0.);
        double mu=0.;
        // median of the shifted deviation distribution
        mu=exp(0.);
        Nu0 = (7.-10.*Ef+5.*pow(Ef,2.))*(1.+0.7*pow(Rep,0.2)*pow(m_Prandtl,1./3.))+(1.33-2.4*Ef+1.2*pow(Ef,2.))*pow(Rep,0.7)*pow(m_Prandtl,1./3.);
        // Compute the parameter a
        a_c_fl_grn=exp(-b_fluidTimeStep/tau);
        mX=exp(0.5*pow(alphaNu,2.));
        sX=pow(exp(2*pow(alphaNu,2.))-exp(pow(alphaNu,2.)),0.5);
        a_c_fl=mX*mX/(sX*sX)*(exp(a_c_fl_grn*log(pow(sX,2)/pow(mX,2)+1))-1);
          if ( int(time/dt)%int(b_fluidTimeStep/dt)==0 )
          {
           // PRS data are given here: 4.5<Re<7 ==> for Liquid/Solid
            if (Rep<3.)
            {
              tau=0.02403;
            }
            else if (Rep>10.5)
            {
              tau=0.0144;
            }
            else
            {
             tau=0.047*pow(Rep,-0.513);
            }
            //a_c_fl=pow(0.4,b_fluidTimeStep/tau);
            //a_c_fl=0.4;
            // Compute the parameter a
            a_c_fl_grn=exp(-b_fluidTimeStep/tau);
            mX=exp(0.5*pow(alphaNu,2.));
            sX=pow(exp(2*pow(alphaNu,2.))-exp(pow(alphaNu,2.)),0.5);
            a_c_fl=mX*mX/(sX*sX)*(exp(a_c_fl_grn*log(pow(sX,2)/pow(mX,2)+1))-1);
            // Get the old data
            gaussian_rand =*((*il)->get_rnd_Nu());
            // Compute log-normal distribution of the deviation
            Nudev=exp(log(mu)+a_c_fl*(log(gaussian_rand[0])-log(mu))+log(exp(gaussian_rand[1])));
            // Compute the fluctuating term
            gaussian_rand[0]=Nudev;
            gaussian_rand[1] = rand_normal(0.,alphaNu*pow(1.-pow(a_c_fl,2.),0.5));
            if ( gaussian_rand[1]>2.*alphaNu*pow(1.-pow(a_c_fl,2.),0.5) )
                gaussian_rand[1] = 2.*alphaNu*pow(1.-pow(a_c_fl,2.),0.5);
            if ( gaussian_rand[1]<-2.*alphaNu*pow(1.-pow(a_c_fl,2.),0.5) )
                gaussian_rand[1] = -2.*alphaNu*pow(1.-pow(a_c_fl,2.),0.5);
            gaussian_rand[2] = a_c_fl;
            (*il)->set_rnd_Nu(gaussian_rand[0],gaussian_rand[1],
                gaussian_rand[2]);
          }
          else gaussian_rand = *((*il)->get_rnd_Nu());
          Nudev=exp(log(mu)+a_c_fl*(log(gaussian_rand[0])-log(mu))+log(exp(gaussian_rand[1])));
          Nu_prime = beta+Nudev;
          Nu = Nu0*(1.+Nu_prime);
      }
     //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "Deen")
      {
	Nu = (7.-10.*Ef+5.*pow(Ef,2.))*(1.+0.17*pow(Rep,0.2)*pow(m_Prandtl,1./3.)) +
           (1.33-2.31*Ef+1.16*pow(Ef,2.))*pow(Rep,0.7)*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "Gnielinski")
      {
        Nulam = 0.664*pow(Rep,0.5)*pow(m_Prandtl,1./3.);
	Nuturb = 0.037*pow(Rep,0.8)*m_Prandtl/(1+2.443*pow(Rep,-0.1)*(pow(m_Prandtl,2./3.)-1));
	Nuiso = 2.+pow(pow(Nulam,2.)+pow(Nuturb,2.),0.5);
	Nu = (1.+1.5*(1.-Ef))*Nuiso;
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "Whitaker")
      {
        Nu = (1.-Ef)/Ef*(0.5*pow(Rep,0.5)+0.2*pow(Rep,2./3.))*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "Rowe")
      {
        A = 2./(1.-pow(1.-Ef,1./3.));
	B = 2./(3.*Ef);
	n = (2.+0.65*pow(Ef*Rep,-0.28))/(3.*(1.+4.65*pow(Ef*Rep,-0.28)));
        Nu = A+B*pow(Ef*Rep,n)*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "Jeschar")
      {
        Nu = 2.+1.12*pow((1.-Ef)/Ef*Rep,0.5)*pow(m_Prandtl,1./3.)+0.0056*Rep*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "OurE03_2Dp")
      {
	Nu = dp/m_diffCoefF*2.500*pow(Rep,0.334)*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if ( type_Nusselt == "OurE03_1Dp")
      {
	Nu = dp/m_diffCoefF*3.568*pow(Rep,0.172)*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if (type_Nusselt == "LittDataFit")
      {
        Nu =2.*(1.+pow(Ef,2./3.-pow(Ef,2.)))+(0.552-4.5908*(1.-Ef)+17.392*pow(1.-Ef,2.)-
	    14.623*pow(1.-Ef,3.))*pow(Rep,0.6)*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if (type_Nusselt == "MyDataFit")
      {
        A=14.105*pow(1.-Ef,2.)+3.121*(1.-Ef)+1.60;
	B=-0.602*pow(1.-Ef,2.)+0.131*(1.-Ef)+0.3804;
	Nu=A*pow(Ef*Rep*m_Prandtl,B);
      }
      //////////////////////////////////////////////////////////////
      else if (type_Nusselt == "Subramaniam2016")
      {
        A=(-0.46+1.77*Ef+0.69*pow(Ef,2.))/pow(Ef,3.);
	B=1.37-2.4*Ef+1.2*pow(Ef,2.);
	Nu=A+B*pow(Rep,0.7)*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if (type_Nusselt == "Deen2014")
      {
        Nu = (7.-10.*Ef+5.*pow(Ef,2.))*(1.+0.17*pow(Ef*Rep,0.2)*pow(m_Prandtl,1./3.)) +
           (1.33-2.31*Ef+1.16*pow(Ef,2.))*pow(Ef*Rep,0.7)*pow(m_Prandtl,1./3.);
      }
      //////////////////////////////////////////////////////////////
      else if (type_Nusselt == "LocalAnalysis")
      {
        A=(1.0395*Ef+8.24905)*pow(m_cellDiamRatio,-2.795*Ef+1.2295);
	B=(0.5375*Ef-0.26695)*log(m_cellDiamRatio)-0.207*Ef+0.2765;
	Nu=A*pow(Ef*Rep*m_Prandtl,B);
        //cout << "Particle : " << (*il) << endl;
        //cout << "DX/dp : " << m_cellDiamRatio << endl;
        //cout << "A : " << A << endl;
        //cout << "B : " << B << endl;
        //cout << "Nu : " << Nu << endl;

        //int fi = rand() % 100 + 1;
        //double f = (double)fi / 100. - 0.5;
        //int sign=(f > 0) ? 1 : ((f < 0) ? -1 : 0);
        //Nu+=double(sign)*30.*f*exp(-pow(1.-double(sign)*f,2));

        // FLO
        double mean=0.;
        double stddev=4.;
        static double n2 = 0.0;
        static int n2_cached = 0;
        double result=0.;
        Vecteur GRN(0.,0.,0.);
        //cout << "Time : " << time << endl;
        if ( time == dt )
        {
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
                 result = n1*stddev + mean;
                 n2_cached = 1;
                 //cout << "First res : " << result << endl;
             }
         }
         else
         {
             n2_cached = 0;
             result = n2*stddev + mean;
             //cout << "First res : " << result << endl;
         }

        (*il)->set_rnd(result,result,result);
       }
       else
       {
         GRN = *((*il)->get_rnd());
         result = GRN[0];
         //cout <<  "res : " << result << endl;
       }
       Nu+=result;
      }
      //////////////////////////////////////////////////////////////
      else
      {
	cout <<"ERROR with type_Nusselt"<<endl;
	exit(0) ;
      }

      heatFlux = PI * dp * m_diffCoefF * Nu * (*Tf-*Tp);
/*cout.precision(17);
cout << "AppFluid_Temperature line 568 DEBUG \n";
cout << " heatFlux = " << heatFlux << endl;
cout << " *Tf = " << *Tf << endl;
cout << " *Tp = " << *Tp << endl;
*/
      // Add it to the sum of heat fluxes on particle
      (*il)->add_heatFlux( heatFlux );

      // Set it as a property of the particle
      (*il)->set_fluidSolidHeatFlux( heatFlux );
      (*il)->set_solidNusselt( Nu );

    } // end loop on particles
  }
}

//----------------------------------------------------------------------------
// Generate random number using Box Muller distribution
double AppFluide_Temperature::rand_normal(double mean, double stddev)
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
void AppFluide_Temperature::set_simultime( double dt )
{
  b_fluidTimeStep = dt;
}
