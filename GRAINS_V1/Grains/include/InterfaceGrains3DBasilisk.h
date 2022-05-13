/** 
# Interface functions for Grains/Basilisk
*/

#ifndef _INTERFACEGRAINS3DBASILISK_H_ 
#define _INTERFACEGRAINS3DBASILISK_H_ 

#ifdef __cplusplus
extern "C" {
#endif

  void Init_Grains ( char const* inputfile,
  	double fluid_density, const bool b_restart );
  
  void Simu_Grains( const double dt_fluid, const bool explicit_added_mass );
  
  char* GrainsToBasilisk( int* pstrsize );

  void SetInitialTime( double tinit );
  
  void Setdt_AddedMassGrains( double dtfluid );

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
