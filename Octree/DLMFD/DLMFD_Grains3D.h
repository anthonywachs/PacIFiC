/** 
# The Grains3D plugin 
*/
 
# include "DLMFD_Plugin.h"

/** File names definition */
# ifndef GRAINS_RESULT_DIR
#   define GRAINS_RESULT_DIR "Grains/Simu"
# endif
# ifndef GRAINS_INPUTFILE
#   define GRAINS_INPUTFILE "Grains/Simu/simul.xml"
# endif

/** Split explicit acceleration treatment in case of particles are lighter 
than the fluid */
# ifndef B_SPLIT_EXPLICIT_ACCELERATION
#   define B_SPLIT_EXPLICIT_ACCELERATION 0
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
    Init_Grains( GRAINS_INPUTFILE, FLUID_DENSITY, restarted_simu );

    // Set initial time
    SetInitialTime( trestart );

    // In case part of the particle acceleration is computed explicitly
    // when particles are lighter than the fluid
    if ( B_SPLIT_EXPLICIT_ACCELERATION )
    {
      // TO DO
    }

    // Get the number of rigid bodies sent by Grains to Basilisk
    // Note: this number must be constant over the simulation
    NumberOfRigidBodiesInBasilisk( &nbParticles, &NbObstacles );
    
    // Transfer the data from Grains to an array of characters
    pstr = GrainsToBasilisk( &pstrsize ); 

    // Check that Paraview writer is activated
    checkParaviewPostProcessing_Grains( GRAINS_RESULT_DIR );
    
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
 
# if _MPI
    // Broadcast the number of rigid bodies
    MPI_Bcast( &nbParticles, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &NbObstacles, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );    
# endif
  nbRigidBodies = nbParticles + NbObstacles;  
  if ( RIGIDBODIES_AS_FIXED_OBSTACLES )
  {
    nbParticles = 0;
    NbObstacles = nbRigidBodies;
  }    
  
  // Allocate number of rigid bodies dependent arrays
  allocate_np_dep_arrays( nbRigidBodies, nbParticles, NRBDATA, &allRigidBodies, 
  	&DLMFDtoGS_vel, &vpartbuf, &pdata, !RIGIDBODIES_AS_FIXED_OBSTACLES,
	&fdata );

  // Update all rigid body data 
  pstr = UpdateParticlesBasilisk( pstr, pstrsize, allRigidBodies, nbRigidBodies,
  	!B_SPLIT_EXPLICIT_ACCELERATION, FLUID_DENSITY ); 	 
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
  // in serial only)
  if ( pid() == 0 ) 
  {
    if ( RIGIDBODIES_AS_FIXED_OBSTACLES )
    {  
      // Output the call to Grains3D
      printf ("Grains3D sends RB data\n");
    }
    else
    {    
      // Output the call to Grains3D
      printf ("Grains3D runs and sends RB data\n");      
      
      // Run the granular simulation
      Simu_Grains( dt );
    }

    // Transfer the data from Grains to an array of characters
    pstr = GrainsToBasilisk( &pstrsize );     
    
    // Set dt for split explicit acceleration computation in Grains3D at 
    // next time step
    if ( B_SPLIT_EXPLICIT_ACCELERATION ) 
    {
      // TO DO
    }    
  }
    
  // Update all rigid body data 
  pstr = UpdateParticlesBasilisk( pstr, pstrsize, allRigidBodies, nbRigidBodies,
    	!B_SPLIT_EXPLICIT_ACCELERATION, FLUID_DENSITY );  
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
    UpdateDLMFDtoGS_vel( DLMFDtoGS_vel, allRigidBodies, nbRigidBodies );  

    // Update particle velocity in Grains using the 2D array
    UpdateVelocityGrains( DLMFDtoGS_vel, nbRigidBodies, 
    	B_SPLIT_EXPLICIT_ACCELERATION );
  }
}
