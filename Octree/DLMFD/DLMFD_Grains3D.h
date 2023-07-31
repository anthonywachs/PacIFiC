/** 
# The Grains3D plugin 
*/
 
# include "DLMFD_Plugin.h"

/** File names definition */
# ifndef grains_result_dir
#   define grains_result_dir "Grains/Simu"
# endif
# ifndef grains_inputfile
#   define grains_inputfile "Grains/Simu/simul.xml"
# endif

/** Split explicit acceleration treatment in case of particles are lighter 
than the fluid */
# ifndef b_split_explicit_acceleration
#   define b_split_explicit_acceleration false
# endif

/** Coupling Interface for Grains3D */
# include "InterfaceGrains3DBasilisk.h"

/** Additional functions for the coupling with Grains3D */
# include "BasiliskGrains3DCouplingFunctions.h"


/** Here we overload the generic events defined in the general DLMFD plugin 
DLMFD_Plugin.h such that it uses Grains3D as a granular solver */


/** Overloading of the granular solver init event */
//----------------------------------------------------------------------------
event GranularSolver_init (t < -1.)
//----------------------------------------------------------------------------
{
  // Initialize Grains with its parameters 
  char* pstr = NULL;
  int pstrsize = 0;  

  // Grains runs in sequential 
  if ( pid() == 0 )
  {
    // Output the call to Grains3D
    printf( "Grains3D\n" );
    
    // Initialize Grains
    Init_Grains( grains_inputfile, rhoval, restarted_simu );

    // Set initial time
    SetInitialTime( trestart );

    // In case part of the particle acceleration is computed explicitly
    // when particles are lighter than the fluid
    if ( b_split_explicit_acceleration )
    {
      // TO DO
    }
    
    // Transfer the data from Grains to an array of characters
    pstr = GrainsToBasilisk( &pstrsize ); 

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

  // Update all particle data 
  pstr = UpdateParticlesBasilisk( pstr, pstrsize, particles, NPARTICLES,
  	!b_split_explicit_acceleration, rhoval );  
  free( pstr );      
} 




/** Overloading of the granular solver predictor event */
//----------------------------------------------------------------------------
event GranularSolver_predictor (t < -1.)
//----------------------------------------------------------------------------
{
  char* pstr = NULL;
  int pstrsize = 0;
  
  // Predictor step: pure granular problem solved by Grains (Grains works
  // in serial only )
  if ( pid() == 0 ) 
  {
    // Output the call to Grains3D
    printf ("run Grains3D\n");
    
    // Run the granular simulation
    Simu_Grains( dt );

    // Transfer the data from Grains to an array of characters
    pstr = GrainsToBasilisk( &pstrsize );     
    
    // Set dt for split explicit acceleration computation in Grains3D at 
    // next time step
    if ( b_split_explicit_acceleration ) 
    {
      // TO DO
    }    
  }
    
  // Update all particle data 
  pstr = UpdateParticlesBasilisk( pstr, pstrsize, particles, NPARTICLES,
    	!b_split_explicit_acceleration, rhoval );  
  free( pstr ); 
}




/** Overloading of the granular solver velocity update event */
//----------------------------------------------------------------------------
event GranularSolver_updateVelocity (t < -1.)
//----------------------------------------------------------------------------
{    
  if ( pid() == 0 )
  {
    // Output the call to Grains3D
    printf ("Grains3D\n");

    // Copy velocity in a 2D array
    UpdateDLMFDtoGS_vel( DLMFDtoGS_vel, particles, NPARTICLES );  

    // Update particle velocity in Grains using the 2D array
    UpdateVelocityGrains( DLMFDtoGS_vel, NPARTICLES, 
    	b_split_explicit_acceleration );
  }
}
