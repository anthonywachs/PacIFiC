/* The Grains3D plugin */ 
# include "dlmfd-plugin.h"

/* File names definition */
# ifndef grains_result_dir
# define grains_result_dir "Grains/Simu"
# endif
# ifndef grains_inputfile
# define grains_inputfile "Grains/Simu/simul.xml"
# endif

/* Coupling Interface for Grains3D */
# include "InterfaceGrains.h"

/* Additional helper functions for the coupling with Grains3D */
# include "BasiliskGrainsCouplingFunctions.h"

struct BasiliskDataStructure BasiliskData[NPARTICLES];

/* Here we overload the generic events defined in the general DLMFD plugin
   DLMFD.h such that it uses Grains3D as a granular solver */


/** Overloading of the granular solver init event */
// -------------------------------------------------
event GranularSolver_init (t < -1.)
{
  // Initialize Grains with its parameters 
  bool b_intitializeClonePer = false;
  double grid_size = 0.;
  bool is_solidsolver_parallel = false;
  int my_rank = pid();
  int nb_procs = npe();
  char* pstr = NULL;
  int pstrsize = 0;  

  // Grains runs in sequential 
  if ( pid() == 0 )
  {
    // Output the call to Grains3D
    printf( "Grains3D\n" );
    
    // Initialize Grains
    Init_Grains( grains_inputfile, rhoval, restarted_simu, 
    	b_intitializeClonePer, grid_size,
	is_solidsolver_parallel, my_rank, nb_procs );

    // Set initial time
    SetInitialTime( trestart );

    // Activate explicit added mass if solid density < fluid density
    if ( b_explicit_added_mass )
    {
      ActivateExplicitAddedMass( restarted_simu ); 
      char rootfilename[80] = "";
      strcpy( rootfilename, result_dir );
      strcat( rootfilename, "/" );
      strcat( rootfilename, result_particle_vp_rootfilename );      
      InitializeExplicitAddedMass( restarted_simu, rootfilename );
    }
    
    // Transfer the data to the common C structure
//    Data_GrainsToCstruct( &BasiliskData[0], NPARTICLES );
    
    pstr = GrainsToBasilisk( &pstrsize ); 
//     printf("Length of string = %d %lu\n",pstrsize,strlen(pstr));    
//     printf("%s\n",pstr);    

    // Check that Paraview writer is activated
    checkParaviewPostProcessing_Grains( grains_result_dir );
    
    // Set the initial cycle number and do initial post-processing
    if ( restarted_simu )
    { 
      SetInitialCycleNumber( init_cycle_number - 1 );
      
      // In Grains in reload mode the initial post-processing does not
      // output any result but simply open existing files and recover the 
      // initial cycle number
      SaveResults_Grains();
    }
    else SetInitialCycleNumber( init_cycle_number );    
  }
// 
//   // Update Basilisk particle structure
//   UpdateParticlesBasilisk( &BasiliskData[0], particles, NPARTICLES,
//   	b_explicit_added_mass, rhoval );

  // Unallocate the BasiliskDataStructure used for Grains-Basilisk
  // communication.  At this point Basilisk has all the particle data
  // in the structure particles 
//  unallocateBasiliskDataStructure( &BasiliskData[0], NPARTICLES );

  pstr = UpdateParticlesBasilisk2( pstr, pstrsize, particles, NPARTICLES,
  	b_explicit_added_mass, rhoval );  
  free( pstr ); 
  
  if ( pid() == 0 ) print_all_particles( particles );     
} 




/** Overloading of the granular solver predictor event */
// ------------------------------------------------------
event GranularSolver_predictor (t < -1.)
{
  char* pstr = NULL;
  int pstrsize = 0;
  
  // Predictor step: pure granular problem solved by Grains (Grains works
  // in serial only )
  if ( pid() == 0 ) 
  {
    // Output the call to Grains3D
    printf ("run Grains3D\n");

    // Set the fluid time step magnitude in Grains3D
    Setdt_Grains( dt );
    
    // Run the granular simulation
    Simu_Grains( true, false, b_explicit_added_mass );

    // Transfer the data to the common C structure
//    Data_GrainsToCstruct( &BasiliskData[0], NPARTICLES );
    
    pstr = GrainsToBasilisk( &pstrsize ); 
//     printf("Length of string = %d %lu\n",pstrsize,strlen(pstr));    
//     printf("%s\n",pstr);    
    
    // Set dt for explicit mass calculation in Grains3D at next time step
    if ( b_explicit_added_mass ) Setdt_AddedMassGrains( dt ) ;
  }
    
  // Update Basilisk particle structure
//    UpdateParticlesBasilisk( &BasiliskData[0], particles, NPARTICLES, 
//    	b_explicit_added_mass, rhoval );
	
    pstr = UpdateParticlesBasilisk2( pstr, pstrsize, particles, NPARTICLES,
    	b_explicit_added_mass, rhoval );  
    free( pstr ); 
}




/** Overloading of the granular solver velocity update event */
// ------------------------------------------------------------
event GranularSolver_updateVelocity (t < -1.)
{
  // Output the call to Grains3D
  if ( pid() == 0 ) printf ("Grains3D\n");

  // Update Basilisk particle structure  
//  UpdateBasiliskStructure( &BasiliskData[0], particles, NPARTICLES );

  // Update particles velocity on the granular side
//  if ( pid() == 0 )
//    Update_Velocity_Grains( &BasiliskData[0], b_explicit_added_mass );
    
  if ( pid() == 0 )
  {
    UpdateDLMFDtoGS_vel( DLMFDtoGS_vel, particles, NPARTICLES );  
    UpdateVelocityGrains( DLMFDtoGS_vel, NPARTICLES, b_explicit_added_mass );
  }

  // Unallocate the BasiliskDataStructure used for Grains-Basilisk communication
//  unallocateBasiliskDataStructure( &BasiliskData[0], NPARTICLES );
}
