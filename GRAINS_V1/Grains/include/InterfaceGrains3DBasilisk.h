/** 
# Interface functions for Grains/Basilisk
*/

#ifndef _INTERFACEGRAINS3DBASILISK_H_ 
#define _INTERFACEGRAINS3DBASILISK_H_ 

#ifdef __cplusplus
extern "C" {
#endif

  void Init_Grains ( char const* inputfile,
  	double fluid_rho, const bool b_restart,
	const bool b_fluidcorrectedacc );
  
  void Simu_Grains( const double dt_fluid );
  
  char* GrainsToBasilisk( int* pstrsize );

  void SetInitialTime( double tinit );

  void SaveResults_Grains();

  void checkParaviewPostProcessing_Grains( char* solid_resDir );

  void UpdateVelocityGrains( double** arrayv, const int m );
  
  void SetInitialCycleNumber( int cycle0 );
  
  void NumberOfRigidBodiesInBasilisk( size_t* nparticles, size_t* nobstacles );

#ifdef __cplusplus
}
#endif

#endif
