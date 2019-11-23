/**
# Coupling interface for Grains3D and Basilisk 
*/

#include "Grains_BuilderFactory.H"
#include "InterfaceGrains.h"

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


  void Data_GrainsToCstruct( struct BasiliskDataStructure * b, const int m ) 
  {
    grains->WriteParticulesInFluid( b );
  }


  void Data_CstructToGrains( struct BasiliskDataStructure * b ) {}


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


  void Update_Velocity_Grains( struct BasiliskDataStructure * b, 
  	bool explicit_added_mass ) 
  {
    grains->UpdateParticulesVelocities( b, explicit_added_mass );
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
