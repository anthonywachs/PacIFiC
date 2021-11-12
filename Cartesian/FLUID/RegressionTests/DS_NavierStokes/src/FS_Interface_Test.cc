#include <FS_Interface_Test.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DomainBuilder.hh>
#include <FV_DiscreteField.hh>
#include <FV_DiscreteField_Centered.hh>
#include <FV_SystemNumbering.hh>
#include <FV_Mesh.hh>
#include <FV.hh>
#include <FV_TimeIterator.hh>
#include <MAC.hh>
#include <MAC_Root.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_Application.hh>
#include <FS_SolidPlugIn_BuilderFactory.hh>
#include <FS_SolidPlugIn.hh>
#include <FS_Grains3DPlugIn.hh>
#include <DS_AllRigidBodies.hh>
#include <math.h>
#include <cstdlib>
#include <iostream>
using std::endl;


FS_Interface_Test const* FS_Interface_Test::PROTOTYPE
			= new FS_Interface_Test() ;


//---------------------------------------------------------------------------
FS_Interface_Test:: FS_Interface_Test( void )
//--------------------------------------------------------------------------
   : FV_OneStepIteration( "FS_Interface_Test" )
{
   MAC_LABEL( "FS_Interface_Test:: FS_Interface_Test" ) ;
}




//---------------------------------------------------------------------------
FS_Interface_Test*
FS_Interface_Test:: create_replica( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FS_Interface_Test:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FS_Interface_Test* result =
                        new FS_Interface_Test( a_owner, dom, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;

}




//---------------------------------------------------------------------------
FS_Interface_Test:: FS_Interface_Test( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp )
   , PP ( dom->discrete_field( "pressure" ) )
   , UU ( dom->discrete_field( "velocity" ) )   
   , solidSolver( NULL )
   , solidFluid_transferStream( NULL )
   , allrigidbodies( NULL )
   , b_particles_as_fixed_obstacles( false )
{
   MAC_LABEL( "FS_Interface_Test:: FS_Interface_Test" ) ;

   // Call of MAC_Communicator routine to set the rank of each proces and
   // the number of processes during execution of FS_Interface_Test
   macCOMM = MAC_Exec::communicator();
   my_rank = macCOMM->rank();
   nb_ranks = macCOMM->nb_ranks();
   is_master = 0;
   
   // Get space dimension
   dimension = PP->primary_grid()->nb_space_dimensions() ;

   // Treat all particles as fixed obstacles
   if ( exp->has_entry( "Particles_as_FixedObstacles" ) ) 
     b_particles_as_fixed_obstacles = exp->bool_data( 
     	"Particles_as_FixedObstacles" ) ;
   
   // Solid solver type is Grains3D
   solidSolverType = "Grains3D";
   b_solidSolver_parallel = false;
   solidSolver_insertionFile = "Grains/Init/insert.xml";
   solidSolver_simulationFile = "Grains/Res/simul.xml";
   int error = 0;   
   solidSolver = FS_SolidPlugIn_BuilderFactory:: create( solidSolverType,
	solidSolver_insertionFile, solidSolver_simulationFile,
        1., false, b_particles_as_fixed_obstacles, 1., b_solidSolver_parallel, 
	error );
	
}




//---------------------------------------------------------------------------
FS_Interface_Test:: ~FS_Interface_Test( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FS_Interface_Test:: ~FS_Interface_Test" ) ;

   if ( solidSolver ) delete solidSolver;
   if ( solidFluid_transferStream ) delete solidFluid_transferStream;
   if ( allrigidbodies ) delete allrigidbodies;

}




//---------------------------------------------------------------------------
void
FS_Interface_Test:: do_one_inner_iteration( 
	FV_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FS_Interface_Test:: do_one_inner_iteration" ) ;
   MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "FS_Interface_Test:: do_one_inner_iteration" ) ;
   start_solving_timer() ;

   if ( !b_particles_as_fixed_obstacles )
   {
     solidSolver->getSolidBodyFeatures( solidFluid_transferStream );
     allrigidbodies->update( *solidFluid_transferStream );

     // Display the geometric features of all rigid bodies
     string space( 3, ' ' ) ;
     for (size_t i = 0; i < nb_ranks; ++i)
     {
       if ( i == my_rank )
       {       
         MAC::out() << space << "Rank " << my_rank << endl ;
         allrigidbodies->display_geometric( MAC::out(), 3 );
         MAC::out() << endl;      
       }
       macCOMM->barrier();
     }

     // Display the features of all rigid bodies
     for (size_t i = 0; i < nb_ranks; ++i)
     {
       if ( i == my_rank )
       {       
         MAC::out() << space << "Rank " << my_rank << endl ;
         allrigidbodies->display( MAC::out(), 3 );
         MAC::out() << endl;      
       }
       macCOMM->barrier();
     }
   }
      
   stop_solving_timer() ;
   stop_total_timer() ;   	
   
}




//---------------------------------------------------------------------------
void
FS_Interface_Test:: do_before_time_stepping( 
	FV_TimeIterator const* t_it, 
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FS_Interface_Test:: do_before_time_stepping" ) ;
   
   start_total_timer( "FS_Interface_Test:: do_before_time_stepping" ) ; 

   solidFluid_transferStream = NULL;
   solidSolver->getSolidBodyFeatures( solidFluid_transferStream );

   allrigidbodies = new DS_AllRigidBodies( dimension, 
   	*solidFluid_transferStream, b_particles_as_fixed_obstacles );
	
   // Display the geometric features of all rigid bodies
   string space( 3, ' ' ) ;
   for (size_t i = 0; i < nb_ranks; ++i)
   {
     if ( i == my_rank )
     {       
       MAC::out() << space << "Rank " << my_rank << endl ;
       allrigidbodies->display_geometric( MAC::out(), 3 );
       MAC::out() << endl;      
     }
     macCOMM->barrier();
   }

   // Display the features of all rigid bodies
   for (size_t i = 0; i < nb_ranks; ++i)
   {
     if ( i == my_rank )
     {       
       MAC::out() << space << "Rank " << my_rank << endl ;
       allrigidbodies->display( MAC::out(), 3 );
       MAC::out() << endl;      
     }
     macCOMM->barrier();
   }  
        
   stop_total_timer() ;  
      
}




//---------------------------------------------------------------------------
void
FS_Interface_Test:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FS_Interface_Test:: do_after_time_stepping" ) ;  

   start_total_timer( 
   	"FS_Interface_Test:: do_after_time_stepping" ) ;

   // Compute the hydro force and torque on all rigid bodies
   allrigidbodies->compute_hydro_force_torque( PP, UU );
   
   stop_total_timer() ;     

}




//---------------------------------------------------------------------------
void
FS_Interface_Test:: do_before_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( 
   	"FS_Interface_Test:: do_before_inner_iterations_stage" ) ;

   start_total_timer( 
   	"FS_Interface_Test:: do_before_inner_iterations_stage" ) ;
   
   stop_total_timer() ;
      
}




//---------------------------------------------------------------------------
void
FS_Interface_Test:: do_after_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FS_Interface_Test:: do_after_inner_iterations_stage" ) ;
   
   start_total_timer( 
   	"FS_Interface_Test:: do_after_inner_iterations_stage" ) ;	
   
   stop_total_timer() ;
   
}
  
