/** 
# Interface functions for Grains/Basilisk
*/

#ifndef _INTERFACEGRAINS3DBASILISK_H_ 
#define _INTERFACEGRAINS3DBASILISK_H_ 

#ifdef __cplusplus
extern "C" {
#endif

  void Init_Grains ( char const* inputfile,
  	double fluid_density, const bool b_restart,
        const bool b_initializeClonePer,
        const double grid_size,
	const bool is_solidsolver_parallel,
	const int my_rank, const int nb_procs );
  
  void Simu_Grains( bool predictor, const bool isPredictorCorrector, 
  	const bool explicit_added_mass );
  
  char* GrainsToBasilisk( int* pstrsize );

  void SetInitialTime( double tinit );
  
  void Setdt_Grains( double const dtfluid );  
  
  void Setdt_AddedMassGrains( double dtfluid );

  void Setviscosity_Grains( double const viscosity );

  void SaveResults_Grains();

  void checkParaviewPostProcessing_Grains( char* solid_resDir );
	
  void UpdateVelocityGrains( double arrayv[][6], const int m, 
  	bool explicit_added_mass );	
  
  void ActivateExplicitAddedMass( bool restart );
  
  void InitializeExplicitAddedMass( bool b_restart, char* rootfilename );
  
  void SetInitialCycleNumber( int cycle0 );

#ifdef __cplusplus
}
#endif

#endif
