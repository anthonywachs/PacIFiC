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
  
  void Simu_Grains( const double dt_fluid );
  
  char* GrainsToBasilisk( int* pstrsize );

  void SetInitialTime( double tinit );

  void SaveResults_Grains();

  void checkParaviewPostProcessing_Grains( char* solid_resDir );

  void UpdateVelocityGrains( double arrayv[][6], const int m,
  	bool bsplit_explicit_acceleration );
  
  void SetInitialCycleNumber( int cycle0 );
  
  size_t NumberOfRigidBodiesInBasilisk();

#ifdef __cplusplus
}
#endif

#endif
