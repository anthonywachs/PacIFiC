#include <DS_DirectionSplitting.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DomainBuilder.hh>
#include <FV_DiscreteField.hh>
#include <FV_SystemNumbering.hh>
#include <FV_Mesh.hh>
#include <FV_TimeIterator.hh>
#include <MAC.hh>
#include <MAC_Root.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Vector.hh>
#include <MAC_BoolArray2D.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_Application.hh>
#include <PAC_Misc.hh>
#include <time.h>
#include <sys/time.h>
#include <math.h>


DS_DirectionSplitting const* DS_DirectionSplitting::PROTOTYPE
                                                 = new DS_DirectionSplitting() ;


//---------------------------------------------------------------------------
DS_DirectionSplitting:: DS_DirectionSplitting( void )
//--------------------------------------------------------------------------
   : FV_OneStepIteration( "DS_DirectionSplitting" )
   , ComputingTime("Solver")
{
   MAC_LABEL( "DS_DirectionSplitting:: DS_DirectionSplitting" ) ;

}




//---------------------------------------------------------------------------
DS_DirectionSplitting*
DS_DirectionSplitting:: create_replica( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   DS_DirectionSplitting* result =
                        new DS_DirectionSplitting( a_owner, dom, exp );

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;

}

//---------------------------------------------------------------------------
DS_DirectionSplitting:: DS_DirectionSplitting( MAC_Object* a_owner,
		                                   FV_DomainAndFields const* dom,
                                         MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp )
   , ComputingTime("Solver")
   , rho( 1. )
   , mu( 1. )
   , kai( 1. )
   , AdvectionScheme( "TVD" )
   , AdvectionTimeAccuracy( 1 )
   , b_restart( false )
   , is_solids( false )
   , insertion_type ( "Grains3D" )
   , is_stressCal ( false )
   , ViscousStressOrder ( "second" )
   , surface_cell_scale ( 1. )
   , is_surfacestressOUT ( false )
   , stressCalFreq ( 1 )
   , is_par_motion ( false )
   , grid_check_for_solid ( 3.0 )
   , FlowSolver ( 0 )
{
   MAC_LABEL( "DS_DirectionSplitting:: DS_DirectionSplitting" ) ;

   // Is the run a follow up of a previous job
   b_restart = MAC_Application::is_follow();

   // Read Density
   if ( exp->has_entry( "Density" ) ) {
     rho = exp->double_data( "Density" ) ;
     exp->test_data( "Density", "Density>0." ) ;
   }

   // Read Viscosity
   if ( exp->has_entry( "Viscosity" ) ) {
     mu = exp->double_data( "Viscosity" ) ;
     exp->test_data( "Viscosity", "Viscosity>0." ) ;
   }

   // Read the presence of particles
   if ( exp->has_entry( "Particles" ) )
     is_solids = exp->bool_data( "Particles" ) ;

   // Read Kai
   if ( exp->has_entry( "Kai" ) ) {
     kai = exp->double_data( "Kai" ) ;
     exp->test_data( "Kai", "Kai>=0." ) ;
   }

   // Advection scheme
   if ( exp->has_entry( "AdvectionScheme" ) )
     AdvectionScheme = exp->string_data( "AdvectionScheme" );
   if ( AdvectionScheme != "Upwind"
     && AdvectionScheme != "TVD"
     && AdvectionScheme != "Centered" )
   {
     string error_message="   - Upwind\n   - TVD\n   - Centered";
     MAC_Error::object()->raise_bad_data_value( exp,
        "AdvectionScheme", error_message );
   }

   // Advection term time accuracy
   if ( exp->has_entry( "AdvectionTimeAccuracy" ) )
     AdvectionTimeAccuracy = exp->int_data( "AdvectionTimeAccuracy" );
   if ( AdvectionTimeAccuracy != 1 && AdvectionTimeAccuracy != 2 )
   {
     string error_message ="   - 1\n   - 2\n   ";
     MAC_Error::object()->raise_bad_data_value( exp,
	"AdvectionTimeAccuracy", error_message );
   }

   if (is_solids) {
      insertion_type = exp->string_data( "InsertionType" ) ;
      MAC_ASSERT( insertion_type == "Grains3D" ) ;

      // Read weather the sress calculation on particle is ON/OFF
      if ( exp->has_entry( "Stress_calculation" ) )
        is_stressCal = exp->bool_data( "Stress_calculation" ) ;

      if (is_stressCal) {
         if ( exp->has_entry( "ViscousStressOrder" ) ) {
           ViscousStressOrder = exp->string_data( "ViscousStressOrder" );
           if ( ViscousStressOrder != "first"
           && ViscousStressOrder != "second") {
              string error_message="- first\n   - second";
              MAC_Error::object()->raise_bad_data_value( exp,
                 "ViscousStressOrder", error_message );
           }
         }
         surface_cell_scale = exp->double_data( "SurfaceCellScale" ) ;
      }
      // Read if the discretized surface force output os ON/OFF
      if (is_stressCal) {
         if ( exp->has_entry( "Surface_Stress_Output" ))
            is_surfacestressOUT = exp->bool_data( "Surface_Stress_Output" ) ;
         if ( exp->has_entry( "Stress_calculation_frequency" ))
            stressCalFreq = exp->int_data("Stress_calculation_frequency");
      }

      // Read weather the particle motion is ON/OFF
      if ( exp->has_entry( "Particle_motion" ) )
        is_par_motion = exp->bool_data( "Particle_motion" ) ;

      if ( exp->has_entry( "GridCheckforSolid" ) )
         grid_check_for_solid = exp->double_data( "GridCheckforSolid" );

   }

   // Create structure to input in the solver system
   struct DS2NS inputDataNS;
   inputDataNS.rho_ = rho ;
   inputDataNS.mu_ = mu ;
   inputDataNS.kai_ = kai ;
   inputDataNS.AdvectionScheme_ = AdvectionScheme ;
   inputDataNS.AdvectionTimeAccuracy_ = AdvectionTimeAccuracy ;
   inputDataNS.b_restart_ = b_restart ;
   inputDataNS.is_solids_ = is_solids ;
   inputDataNS.insertion_type_ = insertion_type;
   inputDataNS.is_stressCal_ = is_stressCal;
   inputDataNS.ViscousStressOrder_ = ViscousStressOrder;
   inputDataNS.surface_cell_scale_ = surface_cell_scale;
   inputDataNS.is_surfacestressOUT_ = is_surfacestressOUT;
   inputDataNS.stressCalFreq_ = stressCalFreq;
   inputDataNS.is_par_motion_ = is_par_motion;
   inputDataNS.grid_check_for_solid_ = grid_check_for_solid;
   inputDataNS.dom_ = dom ;

   MAC_ModuleExplorer* set = exp->create_subexplorer( 0, "DS_NavierStokes" ) ;
   FlowSolver = DS_NavierStokes::create( this, set, inputDataNS ) ;
   set->destroy() ;

//    // Create the temperature solver
//    struct NavierStokes2Temperature inputData;
//    inputData.rho_ = rho ;
//    inputData.b_restart_ = b_restart ;
//    inputData.dom_ = dom ;
//    inputData.UF_ = UF ;
//    inputData.AdvectionScheme_ = AdvectionScheme ;
//    inputData.ViscousStressOrder_ = ViscousStressOrder;
//    inputData.AdvectionTimeAccuracy_ = AdvectionTimeAccuracy ;
//    inputData.is_solids_ = is_solids ;
//    inputData.is_stressCal_ = is_stressCal ;
//    inputData.Npart_ = Npart ;
//    inputData.loc_thres_ = loc_thres ;
//    inputData.level_set_type_ = level_set_type ;
//    inputData.Npoints_ = Npoints ;
//    inputData.Pmin_ = Pmin ;
//    inputData.ar_ = ar ;
//    inputData.pole_loc_ = pole_loc ;
//    inputData.particle_information_ = &particle_information ;
//    inputData.insertion_type_ = insertion_type ;
//    inputData.is_par_motion_ = is_par_motion ;
//    inputData.IntersectionMethod_ = IntersectionMethod ;
//    inputData.tolerance_ = tolerance ;
//
//    MAC_ModuleExplorer* set = exp->create_subexplorer( 0, "DDS_HeatTransfer" ) ;
//    Solver_Temperature = DDS_HeatTransfer::create( this, set, inputData ) ;
//    set->destroy() ;
//
//    // Timing routines
//    if ( my_rank == is_master )
//    {
//      SCT_insert_app("Matrix_Assembly&Initialization");
//      SCT_insert_app("Pressure predictor");
// //     SCT_insert_app("Explicit Velocity step");
// //     SCT_insert_app("Velocity x update");
// //     SCT_insert_app("Velocity y update");
// //     SCT_insert_app("Velocity z update");
//      SCT_insert_app("Velocity update");
//      SCT_insert_app("Penalty Step");
//      SCT_insert_app("Pressure Update");
//      SCT_insert_app("Writing CSV");
//      SCT_get_elapsed_time("Objects_Creation");
//    }
}

//---------------------------------------------------------------------------
DS_DirectionSplitting:: ~DS_DirectionSplitting( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: ~DS_DirectionSplitting" ) ;

}

//---------------------------------------------------------------------------
void
DS_DirectionSplitting:: do_one_inner_iteration( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: do_one_inner_iteration" ) ;
   MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;


   // Flow solver
   start_total_timer( "DS_NavierStokes:: do_one_inner_iteration" ) ;
   start_solving_timer() ;
   FlowSolver->do_one_inner_iteration( t_it ) ;
   stop_solving_timer() ;
   stop_total_timer() ;

   // Temperature solver
   // Solver_Temperature->do_one_inner_iteration( t_it ) ;


}

//---------------------------------------------------------------------------
void
DS_DirectionSplitting:: do_before_time_stepping( FV_TimeIterator const* t_it,
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: do_before_time_stepping" ) ;


   FV_OneStepIteration::do_before_time_stepping( t_it, basename ) ;

   // Flow solver
   start_total_timer( "DS_NavierStokes:: do_before_time_stepping" ) ;
   FlowSolver->do_before_time_stepping( t_it, basename ) ;
   stop_total_timer() ;

   // Temperature solver



}




//---------------------------------------------------------------------------
void
DS_DirectionSplitting:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: do_after_time_stepping" ) ;


   // Flow solver
   start_total_timer( "DS_NavierStokes:: do_after_time_stepping" ) ;
   FlowSolver->do_after_time_stepping() ;
   stop_total_timer() ;
   // Temperature solver
   // Solver_Temperature->do_after_time_stepping() ;


}

//---------------------------------------------------------------------------
void
DS_DirectionSplitting:: do_before_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: do_before_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_before_inner_iterations_stage( t_it ) ;

   // Flow solver
   start_total_timer( "DS_NavierStokes:: do_before_inner_iterations_stage" ) ;
   FlowSolver->do_before_inner_iterations_stage( t_it ) ;
   stop_total_timer() ;
   // Temperature solver
   // Solver_Temperature->do_before_inner_iterations_stage( t_it ) ;


}

//---------------------------------------------------------------------------
void
DS_DirectionSplitting:: do_after_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: do_after_inner_iterations_stage" ) ;


   FV_OneStepIteration::do_after_inner_iterations_stage( t_it ) ;

   // Flow solver
   start_total_timer( "DS_NavierStokes:: do_after_inner_iterations_stage" ) ;
   FlowSolver->do_after_inner_iterations_stage( t_it ) ;
   stop_total_timer() ;
   // Temperature solver
   // Solver_Temperature->do_after_inner_iterations_stage( t_it ) ;


}

//---------------------------------------------------------------------------
void
DS_DirectionSplitting:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: do_additional_savings" ) ;


   // Flow solver
   start_total_timer( "DS_NavierStokes:: do_additional_savings" ) ;
   FlowSolver->do_additional_savings( t_it, cycleNumber ) ;
   stop_total_timer() ;
   // Temperature solver
   // Solver_Temperature->do_additional_savings( t_it, cycleNumber ) ;


}
