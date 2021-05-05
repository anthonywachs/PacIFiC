/**
# Coupling interface for Grains3D and Basilisk 
*/

#include "Grains_BuilderFactory.H"
#include "InterfaceGrains.h"
#include <iostream>
#include <cstring>
#include <string>

#ifdef __cplusplus
extern "C" {
#endif
  
  static GrainsCoupledWithFluid * grains = NULL;

  
  void Init_Grains ( char const* inputfile, 
  	double fluid_density, const bool b_restart,
        const bool b_initializeClonePer,
        const double grid_size,
	const bool is_solidsolver_parallel,
	const int my_rank, const int nb_procs) 
  {
    string simulation_file( inputfile );
    
    ReaderXML::initialize();
        
    string simulation_file_exe = Grains_BuilderFactory::init( simulation_file, 
    	my_rank, is_solidsolver_parallel ? nb_procs : 1 );
  
    DOMElement* rootNode = ReaderXML::getRoot (simulation_file_exe);
    
    grains = Grains_BuilderFactory::createCoupledWithFluid( rootNode, 
    	fluid_density, grid_size );

    if ( b_restart ) grains->setReloadSame();
    grains->Construction (rootNode);
    grains->Forces (rootNode);
    grains->Chargement (rootNode);
    if ( b_initializeClonePer ) grains->initializeClonesPeriodiques();

    ReaderXML::terminate();
     
    cout << "Construction of Grains completed" << endl;
  }


  void Simu_Grains( const bool predictor, const bool isPredictorCorrector, 
  	const bool explicit_added_mass) 
  {
    grains->Simulation( predictor, isPredictorCorrector, explicit_added_mass );	
  }

  
  char* GrainsToBasilisk( int* pstrsize )
  {
    // We use the interface function of PeliGRIFF
    istringstream iss;
    grains->WriteParticulesInFluid( iss );

    // Temporary: we remove the formatting and separate each entry by " "
    string buff;
    *pstrsize = int(iss.str().size()) + 1;    
    char* pstr = new char [*pstrsize];
    int pos = 0;
    while ( iss >> buff )
    {
      std::strcpy( &pstr[pos], buff.c_str() );
      pos += int(buff.size());
      std::strcpy( &pstr[pos], " " ); 
      pos += 1;     
    }    
    
    // Return the pointer to the char
    return pstr;
  }


  void Setdt_Grains( double const dtfluid ) 
  {
    grains->set_timeStep( dtfluid );
  }

  
  void Setdt_AddedMassGrains( double dtfluid ) 
  {
    grains->set_ExplicitAddedMasstimeStep( dtfluid );
  }  


  void Setviscosity_Grains( double const viscosity ) 
  {
    grains->setFluidViscosity( viscosity );
  }


  void SaveResults_Grains() 
  {
    static unsigned int ppcounter = 0;
    if ( !ppcounter ) 
      grains->InitialPostProcessing( 6 );
    else
      grains->doPostProcessing( 6 );
    ++ppcounter; 
  }  


  void checkParaviewPostProcessing_Grains( char * solid_resDir )
  {    
    grains->checkParaviewPostProcessing( "grains", solid_resDir, true );
  }

  
  void UpdateVelocityGrains( double arrayv[][6], const int m, 
  	bool explicit_added_mass ) 
  {
    // Transfer into a vector< vector<double> >
    vector<double> buf( 6, 0.);
    vector< vector<double> > vecv( m, buf );
    for (size_t i=0;i<size_t(m);++i)
      for (size_t j=0;j<6;++j)
        vecv[i][j] = arrayv[i][j];

    // We use the interface function of PeliGRIFF
    grains->UpdateParticulesVelocities( vecv, explicit_added_mass );
  }  

  
  void ActivateExplicitAddedMass( bool restart ) 
  {
    grains->AddExplicitAddedMass( restart, "" );
  }

  
  void InitializeExplicitAddedMass( bool restart, char * rootfilename ) 
  {
    grains->InitializeExplicitAddedMassRestart( restart, rootfilename,
    	"Basilisk" );
  }

  
  void SetInitialCycleNumber( int cycle0 ) 
  {
    grains->setInitialCycleNumber( cycle0 );
  }

  
  void SetInitialTime( double tinit ) 
  {
    grains->setInitialTime( tinit );
  }    
  
#ifdef __cplusplus
}
#endif
