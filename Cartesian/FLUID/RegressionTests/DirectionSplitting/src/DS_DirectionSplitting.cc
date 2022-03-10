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
#include <FS_SolidPlugIn_BuilderFactory.hh>
#include <FS_SolidPlugIn.hh>
#include <FS_Grains3DPlugIn.hh>
#include <DS_AllRigidBodies.hh>


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
   , is_HE( false )
   , is_NS( false )
   , is_NSwithHE( false )
   , insertion_type ( "Grains3D" )
   , is_stressCal ( false )
   , ViscousStressOrder ( "second" )
   , surface_cell_scale ( 1. )
   , is_surfacestressOUT ( false )
   , stressCalFreq ( 1 )
   , is_par_motion ( false )
   , FlowSolver ( 0 )
   , HeatSolver ( 0 )
   , allrigidbodies ( 0 )
   , b_particles_as_fixed_obstacles( true )
{
   MAC_LABEL( "DS_DirectionSplitting:: DS_DirectionSplitting" ) ;

   // Is the run a follow up of a previous job
   b_restart = MAC_Application::is_follow();
   macCOMM = MAC_Exec::communicator();

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
   if ( AdvectionTimeAccuracy != 1 && AdvectionTimeAccuracy != 2 ) {
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

      // Critical distance
      if ( dom->primary_grid()->is_translation_active() ) {
        if ( exp->has_entry( "Critical_Distance_Translation" ) )
           critical_distance_translation= exp->double_data(
             										"Critical_Distance_Translation" );
        else {
           string error_message=" Projection-Translation is active but ";
           error_message+="Critical_Distance_Translation is NOT defined.";
           MAC_Error::object()->raise_bad_data_value( exp,
              							 "Projection_Translation", error_message );
        }
      }

   }

   // Read the solids filename
   if (is_solids && (insertion_type == "Grains3D")) {
      solidSolverType = "Grains3D";
      b_solidSolver_parallel = false;
      solidSolver_insertionFile = "Grains/Init/insert.xml";
      solidSolver_simulationFile = "Grains/Res/simul.xml";
   }

   if (dom->discrete_field( "temperature" ) &&
       dom->discrete_field( "velocity" )) {
      is_NSwithHE = true;
   } else if (dom->discrete_field( "temperature" )) {
      is_HE = true;
   } else if (dom->discrete_field( "velocity" )) {
      is_NS = true;
   }

   // Create Grains3D if solidSolverType is Grains3D;
   if (is_solids) {
      int error = 0;
      solidSolver = FS_SolidPlugIn_BuilderFactory:: create( solidSolverType,
         solidSolver_insertionFile,
         solidSolver_simulationFile,
         1., false,
         b_particles_as_fixed_obstacles,
         1., b_solidSolver_parallel,
         error );

      solidFluid_transferStream = NULL;
      solidSolver->getSolidBodyFeatures( solidFluid_transferStream );
   }

   size_t dim = (is_HE) ? dom->discrete_field( "temperature" )->primary_grid()
                                                       ->nb_space_dimensions()
                        : dom->discrete_field( "velocity" )->primary_grid()
                                                       ->nb_space_dimensions();

   // Read the gravity vector or direction of enforced motion
   doubleVector gg( dim, 0 );
   if ( exp->has_entry( "Gravity_vector" ) )
      gg = exp->doubleVector_data( "Gravity_vector" );
   gravity_vector = MAC_DoubleVector::create( this, gg );

   // Create rigid bodies objects depending on which PDE to solve
   if (is_solids) {
      if (is_NS) {
         allrigidbodies = new DS_AllRigidBodies( dim
                          , *solidFluid_transferStream
                          , b_particles_as_fixed_obstacles
                          , dom->discrete_field( "velocity" )
                          , dom->discrete_field( "pressure" )
                          , rho
                          , gravity_vector
                          , surface_cell_scale
                          , macCOMM
                          , mu );
      } else if (is_HE) {
         allrigidbodies = new DS_AllRigidBodies( dim
                          , *solidFluid_transferStream
                          , b_particles_as_fixed_obstacles
                          , dom->discrete_field( "temperature" )
                          , surface_cell_scale
                          , macCOMM
                          , mu );
      } else if (is_NSwithHE) {
         allrigidbodies = new DS_AllRigidBodies( dim
                          , *solidFluid_transferStream
                          , b_particles_as_fixed_obstacles
                          , dom->discrete_field( "velocity" )
                          , dom->discrete_field( "pressure" )
                          , dom->discrete_field( "temperature" )
                          , surface_cell_scale
                          , macCOMM
                          , mu );
      }
   }


   // Create structure to input in the NS solver
   if (is_NS || is_NSwithHE) {
      struct DS2NS inputDataNS;
      inputDataNS.rho_ = rho ;
      inputDataNS.mu_ = mu ;
      inputDataNS.kai_ = kai ;
      inputDataNS.AdvectionScheme_ = AdvectionScheme ;
      inputDataNS.AdvectionTimeAccuracy_ = AdvectionTimeAccuracy ;
      inputDataNS.b_restart_ = b_restart ;
      inputDataNS.is_solids_ = is_solids ;
      inputDataNS.is_stressCal_ = is_stressCal;
      inputDataNS.ViscousStressOrder_ = ViscousStressOrder;
      inputDataNS.surface_cell_scale_ = surface_cell_scale;
      inputDataNS.is_surfacestressOUT_ = is_surfacestressOUT;
      inputDataNS.stressCalFreq_ = stressCalFreq;
      inputDataNS.is_par_motion_ = is_par_motion;
      inputDataNS.dom_ = dom;
      inputDataNS.allrigidbodies_ = allrigidbodies;
      inputDataNS.critical_distance_translation_ = critical_distance_translation;

      MAC_ModuleExplorer* set = exp->create_subexplorer( 0, "DS_NavierStokes" ) ;
      FlowSolver = DS_NavierStokes::create( this, set, inputDataNS ) ;
      set->destroy() ;
   }

   // Create the temperature solver
   if (is_HE || is_NSwithHE) {
      struct DS2HE inputDataHE;
      inputDataHE.rho_ = rho ;
      inputDataHE.b_restart_ = b_restart ;
      inputDataHE.AdvectionScheme_ = AdvectionScheme ;
      inputDataHE.ViscousStressOrder_ = ViscousStressOrder;
      inputDataHE.AdvectionTimeAccuracy_ = AdvectionTimeAccuracy ;
      inputDataHE.is_solids_ = is_solids ;
      inputDataHE.is_NSwithHE_ = is_NSwithHE ;
      inputDataHE.is_stressCal_ = is_stressCal ;
      inputDataHE.is_par_motion_ = is_par_motion ;
      inputDataHE.dom_ = dom ;
      inputDataHE.allrigidbodies_ = allrigidbodies ;

      MAC_ModuleExplorer* set = exp->create_subexplorer( 0, "DS_HeatTransfer" ) ;
      HeatSolver = DS_HeatTransfer::create( this, set, inputDataHE ) ;
      set->destroy() ;
   }


}




//---------------------------------------------------------------------------
DS_DirectionSplitting:: ~DS_DirectionSplitting( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: ~DS_DirectionSplitting" ) ;

   if ( is_solids ) {
      if ( solidSolver ) delete solidSolver;
      if ( solidFluid_transferStream ) delete solidFluid_transferStream;
      if ( allrigidbodies ) delete allrigidbodies;
   }

}




//---------------------------------------------------------------------------
void
DS_DirectionSplitting:: do_one_inner_iteration( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: do_one_inner_iteration" ) ;
   MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   // Flow solver
   if (is_NS || is_NSwithHE) {
      start_total_timer( "DS_NavierStokes:: do_one_inner_iteration" ) ;
      start_solving_timer() ;
      FlowSolver->do_one_inner_iteration( t_it ) ;
      stop_solving_timer() ;
      stop_total_timer() ;
   }

   // Temperature solver
   if (is_HE || is_NSwithHE) {
      start_total_timer( "DS_HeatTransfer:: do_one_inner_iteration" ) ;
      start_solving_timer() ;
      HeatSolver->do_one_inner_iteration( t_it ) ;
      stop_solving_timer() ;
      stop_total_timer() ;
   }

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
   if (is_NS || is_NSwithHE) {
      start_total_timer( "DS_NavierStokes:: do_before_time_stepping" ) ;
      FlowSolver->do_before_time_stepping( t_it, basename ) ;
      stop_total_timer() ;
   }

   // Temperature solver
   if (is_HE || is_NSwithHE) {
      start_total_timer( "DS_HeatTransfer:: do_before_time_stepping" ) ;
      HeatSolver->do_before_time_stepping( t_it, basename ) ;
      stop_total_timer() ;
   }


}




//---------------------------------------------------------------------------
void
DS_DirectionSplitting:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: do_after_time_stepping" ) ;

   // Flow solver
   if (is_NS || is_NSwithHE) {
      start_total_timer( "DS_NavierStokes:: do_after_time_stepping" ) ;
      FlowSolver->do_after_time_stepping() ;
      stop_total_timer() ;
   }

   // Temperature solver
   if (is_HE || is_NSwithHE) {
      start_total_timer( "DS_HeatTransfer:: do_after_time_stepping" ) ;
      HeatSolver->do_after_time_stepping() ;
      stop_total_timer() ;
   }

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
   if (is_NS || is_NSwithHE) {
      start_total_timer( "DS_NavierStokes:: do_before_inner_iterations_stage" ) ;
      FlowSolver->do_before_inner_iterations_stage( t_it ) ;
      stop_total_timer() ;
   }
   // Temperature solver
   if (is_HE || is_NSwithHE) {
      start_total_timer( "DS_HeatTransfer:: do_before_inner_iterations_stage" ) ;
      HeatSolver->do_before_inner_iterations_stage( t_it ) ;
      stop_total_timer() ;
   }

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
   if (is_NS || is_NSwithHE) {
      start_total_timer( "DS_NavierStokes:: do_after_inner_iterations_stage" ) ;
      FlowSolver->do_after_inner_iterations_stage( t_it ) ;
      stop_total_timer() ;
   }
   // Temperature solver
   if (is_HE || is_NSwithHE) {
      start_total_timer( "DS_HeatTransfer:: do_after_inner_iterations_stage" ) ;
      HeatSolver->do_after_inner_iterations_stage( t_it ) ;
      stop_total_timer() ;
   }


}




//---------------------------------------------------------------------------
void
DS_DirectionSplitting:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_DirectionSplitting:: do_additional_savings" ) ;


   // Flow solver
   if (is_NS || is_NSwithHE) {
      start_total_timer( "DS_NavierStokes:: do_additional_savings" ) ;
      FlowSolver->do_additional_savings( t_it, cycleNumber ) ;
      stop_total_timer() ;
   }
   // Temperature solver
   if (is_HE || is_NSwithHE) {
      start_total_timer( "DS_HeatTransfer:: do_additional_savings" ) ;
      HeatSolver->do_additional_savings( t_it, cycleNumber ) ;
      stop_total_timer() ;
   }


}
