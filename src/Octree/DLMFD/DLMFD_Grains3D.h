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

/** Coupling Interface for Grains3D */
# include "DLMFD_InterfaceGrains3DBasilisk.h"

/** Additional functions for the coupling with Grains3D */
# include "DLMFD_BasiliskGrains3DCouplingFunctions.h"


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

  // Grains runs in serial 
  if ( pid() == 0 )
  {
    // Output the call to Grains3D
    printf( "Grains3D\n" );
    
    // Initialize Grains
    Init_Grains( GRAINS_INPUTFILE, FLUID_DENSITY, restarted_simu,
    	!B_SPLIT_EXPLICIT_ACCELERATION );

    // Set initial time
    SetInitialTime( trestart );
    
    // Get the number of reference rigid bodies sent by Grains to Basilisk
    // Note: this number must be constant over the simulation
    NumberOfReferenceRigidBodiesInBasilisk( &nbReferenceRigidBodies );      
    
    // Transfer the data from Grains to an array of characters
    pstr = GrainsToBasiliskReference( &pstrsize );
         
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
    MPI_Bcast( &nbReferenceRigidBodies, 1, MPI_UNSIGNED_LONG, 0, 
    	MPI_COMM_WORLD );    
# endif

  // Create the array of reference rigid bodies
  ReferenceRigidBodies = (RigidBody*) calloc( nbReferenceRigidBodies, 
  	sizeof(RigidBody) );
  pstr = CreateReferenceRBBasilisk( pstr, pstrsize, ReferenceRigidBodies,
  	nbReferenceRigidBodies, FLUID_DENSITY );
  free( pstr );
 
  if ( pid() == 0 )
  {
    // Get the number of rigid bodies sent by Grains to Basilisk
    // Note: this number must be constant over the simulation
    NumberOfRigidBodiesInBasilisk( &nbParticles, &NbObstacles );
    
    // Transfer the data from Grains to an array of characters
    pstr = GrainsToBasilisk( &pstrsize );       
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
  	&RBnumToIndex, &DLMFDtoGS_vel, &vpartbuf, &pdata, 
	!RIGIDBODIES_AS_FIXED_OBSTACLES, &fdata );

  // Update all rigid body data 
  pstr = UpdateParticlesBasilisk( pstr, pstrsize, allRigidBodies, nbRigidBodies,
  	ReferenceRigidBodies, FLUID_DENSITY, true );	 	 
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
  }
    
  // Update all rigid body data 
  pstr = UpdateParticlesBasilisk( pstr, pstrsize, allRigidBodies, nbRigidBodies,
    	ReferenceRigidBodies, FLUID_DENSITY, false );  
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
    UpdateVelocityGrains( DLMFDtoGS_vel, nbRigidBodies );
  }
}




/** Overloading of the granular solver result saving event */
//----------------------------------------------------------------------------
event GranularSolver_saveResults (t < -1.)
//----------------------------------------------------------------------------
{
  if ( pid() == 0 ) SaveResults_Grains();  
}
