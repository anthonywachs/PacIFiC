#include <DDS_NavierStokes.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DomainBuilder.hh>
#include <FV_DiscreteField.hh>
#include <DDS_NavierStokesSystem.hh>
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
#include <intVector.hh>
#include <LA_Vector.hh>
#include <LA_SeqMatrix.hh>
#include <LA_MatrixIterator.hh>
#include <PAC_Misc.hh>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <FS_SolidPlugIn_BuilderFactory.hh>
#include <FS_SolidPlugIn.hh>
#include <FS_Grains3DPlugIn.hh>
#include <DS_AllRigidBodies.hh>



DDS_NavierStokes const* DDS_NavierStokes::PROTOTYPE
                                                 = new DDS_NavierStokes() ;


//---------------------------------------------------------------------------
DDS_NavierStokes:: DDS_NavierStokes( void )
//--------------------------------------------------------------------------
   : FV_OneStepIteration( "DDS_NavierStokes" )
   , ComputingTime("Solver")
{
   MAC_LABEL( "DDS_NavierStokes:: DDS_NavierStokes" ) ;

}




//---------------------------------------------------------------------------
DDS_NavierStokes*
DDS_NavierStokes:: create_replica( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   DDS_NavierStokes* result =
                        new DDS_NavierStokes( a_owner, dom, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;

}

//---------------------------------------------------------------------------
DDS_NavierStokes:: DDS_NavierStokes( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp )
   , ComputingTime("Solver")
   , UF ( dom->discrete_field( "velocity" ) )
   , PF ( dom->discrete_field( "pressure" ) )
   , GLOBAL_EQ( 0 )
   , peclet( 1. )
   , mu( 1. )
   , kai( 1. )
   , AdvectionScheme( "TVD" )
   , AdvectionTimeAccuracy( 1 )
   , DivRelax( 10 )
   , rho( 1. )
   , is_solids( false )
   , is_par_motion( false )
   , is_stressCal( false )
   , is_surfacestressOUT( false )
   , IntersectionMethod ( "Bisection" )
   , ViscousStressOrder ( "second" )
   , PressureStressOrder ( "first" )
   , tolerance ( 1.e-6 )
   , grid_check_for_solid ( 3.0 )
   , b_particles_as_fixed_obstacles( true )
   , b_projection_translation( dom->primary_grid()->is_translation_active() )
   , b_grid_has_been_translated_since_last_output( false )
   , b_grid_has_been_translated_at_previous_time( false )
   , critical_distance_translation( 0. )
   , translation_direction( 0 )
   , bottom_coordinate( 0. )
   , translated_distance( 0. )
   , gravity_vector( 0 )
{
   MAC_LABEL( "DDS_NavierStokes:: DDS_NavierStokes" ) ;
   MAC_ASSERT( UF->discretization_type() == "staggered" ) ;
   MAC_ASSERT( PF->discretization_type() == "centered" ) ;
   if (b_projection_translation) {
      MAC_ASSERT( PF->storage_depth() == 3 ) ;
      MAC_ASSERT( UF->storage_depth() == 6 ) ;
   } else {
      MAC_ASSERT( PF->storage_depth() == 2 ) ;
      MAC_ASSERT( UF->storage_depth() == 5 ) ;
   }

   // Call of MAC_Communicator routine to set the rank of each proces and
   // the number of processes during execution
   macCOMM = MAC_Exec::communicator();
   my_rank = macCOMM->rank();
   nb_procs = macCOMM->nb_ranks();
   is_master = 0;

   is_periodic[0][0] = false;
   is_periodic[0][1] = false;
   is_periodic[0][2] = false;
   is_periodic[1][0] = false;
   is_periodic[1][1] = false;
   is_periodic[1][2] = false;

   // Timing routines
   if ( my_rank == is_master )
   {
     CT_set_start();
     SCT_insert_app("Objects_Creation");
     SCT_set_start("Objects_Creation");
   }

   // Is the run a follow up of a previous job
   b_restart = MAC_Application::is_follow();


   // Clear results directory in case of a new run
   if( !b_restart ) PAC_Misc::clearAllFiles( "Res", "Savings", my_rank ) ;

   // Get space dimension
   dim = UF->primary_grid()->nb_space_dimensions() ;
   nb_comps[0] = PF->nb_components() ;
   nb_comps[1] = UF->nb_components() ;

   if ( dim == 1 )
   {
     string error_message="Space dimension should either 2 or 3";
     MAC_Error::object()->raise_bad_data_value(exp,
	"nb_space_dimensions",
	error_message );
   }

   // Create the Direction Splitting subcommunicators
   create_DDS_subcommunicators();

   // Read Density
   if ( exp->has_entry( "Density" ) )
   {
     rho = exp->double_data( "Density" ) ;
     exp->test_data( "Density", "Density>0." ) ;
   }

   // Read Viscosity
   if ( exp->has_entry( "Viscosity" ) )
   {
     mu = exp->double_data( "Viscosity" ) ;
     exp->test_data( "Viscosity", "Viscosity>0." ) ;
   }

   // Read Intersection tolerance
   if ( exp->has_entry( "IntersectionTolerance" ) )
   {
     tolerance = exp->double_data( "IntersectionTolerance" ) ;
     exp->test_data( "IntersectionTolerance", "IntersectionTolerance>0." ) ;
   }

   // Read the presence of particles
   if ( exp->has_entry( "Particles" ) )
     is_solids = exp->bool_data( "Particles" ) ;

   if (is_solids) {
      Npart = exp->int_data( "NParticles" ) ;
      insertion_type = exp->string_data( "InsertionType" ) ;
      MAC_ASSERT( insertion_type == "file" || insertion_type == "Grains3D" ) ;
      loc_thres = exp->double_data( "Local_threshold" ) ;
      if ( exp->has_entry( "LevelSetType" ) )
         level_set_type = exp->string_data( "LevelSetType" );
      if ( level_set_type != "Cube" && level_set_type != "Cylinder" &&
           level_set_type != "Sphere" && level_set_type != "Ellipsoid" &&
           level_set_type != "PipeX" && level_set_type != "Superquadric") {
         string error_message="- Cube\n   - Sphere\n   - Cylinder\n   "
                              "- Superquadric\n   - Ellipsoid\n   - PipeX";
         MAC_Error::object()->raise_bad_data_value( exp,
                                       "LevelSetType", error_message );
      }

      // Read the solids filename
      if (insertion_type == "file") {
         solidSolverType = "Grains3D";
         b_solidSolver_parallel = false;
         solidSolver_insertionFile = "Grains/Init/insert.xml";
         solidSolver_simulationFile = "Grains/Res/simul.xml";
         solid_filename = exp->string_data( "Particle_FileName" );
      }

      // Read weather the sress calculation on particle is ON/OFF
      if ( exp->has_entry( "Stress_calculation" ) )
        is_stressCal = exp->bool_data( "Stress_calculation" ) ;

      // Read if the discretized surface force output os ON/OFF
      if (is_stressCal) {
         if ( exp->has_entry( "Surface_Stress_Output" ))
            is_surfacestressOUT = exp->bool_data( "Stress_calculation" ) ;
      }

      // Read weather the particle motion is ON/OFF
      if ( exp->has_entry( "Particle_motion" ) )
        is_par_motion = exp->bool_data( "Particle_motion" ) ;

      if ( exp->has_entry( "GridCheckforSolid" ) )
         grid_check_for_solid = exp->double_data( "GridCheckforSolid" );


      if (is_par_motion && (insertion_type=="file")) {
         // Read the type for particle motion
         if ( exp->has_entry( "Motion_type" ) ) {
            motion_type = exp->string_data( "Motion_type" ) ;
            MAC_ASSERT( motion_type == "Sine" || motion_type == "Hydro" ) ;
         }
         // Read the gravity vector or direction of enforced motion
         doubleVector gg( dim, 0 );
         if ( exp->has_entry( "Gravity_vector" ) )
            gg = exp->doubleVector_data( "Gravity_vector" );
         gravity_vector = MAC_DoubleVector::create( this, gg );

         if (motion_type == "Sine") {
            Amp = exp->double_data( "Amplitude" ) ;
            freq = exp->double_data( "Frequency" ) ;
         } else if (motion_type == "Hydro") {
            rho_s = exp->double_data( "Solid_Density" );
         }
      }

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

         if ( exp->has_entry( "PressureStressOrder" ) ) {
            PressureStressOrder = exp->string_data( "PressureStressOrder" );
            if ( PressureStressOrder != "first"
              && PressureStressOrder != "second"
              && PressureStressOrder != "second_withNeumannBC") {
               string error_message="- first\n   "
                                    "- second\n   "
                                    "- second_withNeumannBC";
               MAC_Error::object()->raise_bad_data_value( exp,
                                    "PressureStressOrder", error_message );
            }
         }

         Npoints = exp->double_data( "Npoints" ) ;

         if (dim == 3) {
            if ((level_set_type == "Sphere") || (level_set_type == "Cylinder")) {
               Pmin = exp->int_data( "Pmin" ) ;
               ar = exp->double_data( "aspect_ratio" ) ;
               pole_loc = exp->int_data( "pole_loc" ) ;
            }
         }
      }
   }

   // Read Kai
   if ( exp->has_entry( "Kai" ) )
   {
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
   if ( AdvectionScheme == "TVD"
   	&& UF->primary_grid()->get_security_bandwidth() < 2 )
   {
     string error_message="   >= 2 with TVD scheme";
     MAC_Error::object()->raise_bad_data_value( exp,
        "security_bandwidth", error_message );
   }

   if ( is_stressCal == true &&
        UF->primary_grid()->get_security_bandwidth() < 4 )
   {
     string error_message="   >= 4 for correct stress calculations on solids";
     MAC_Error::object()->raise_bad_data_value( exp,
        "security_bandwidth", error_message );
   }

   // Critical distance
   if( b_projection_translation )
   {
     if( exp->has_entry( "Critical_Distance_Translation" ) )
       critical_distance_translation=exp->double_data(
           "Critical_Distance_Translation" );
     else
     {
       string error_message=" Projection-Translation is active but ";
       error_message+="Critical_Distance_Translation is NOT defined.";
       MAC_Error::object()->raise_bad_data_value( exp,
           "Projection_Translation", error_message );
     }
   }

   // Method for calculating intersections
   if ( exp->has_entry( "IntersectionMethod" ) )
   {
     IntersectionMethod = exp->string_data( "IntersectionMethod" ) ;
     if ( IntersectionMethod != "Bisection" && IntersectionMethod != "Newton") {
        string error_message="   - Bisection\n   - Newton";
        MAC_Error::object()->raise_bad_data_value( exp,
           "IntersectionMethod", error_message );
     }
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

   // Number of iteration to relax the change in divergence stencil
   if ( exp->has_entry( "DivRelax" ) ) {
     DivRelax = exp->int_data( "DivRelax" );
     exp->test_data( "DivRelax", "DivRelax>=0" ) ;
   }

   // Periodic boundary condition check for velocity
   U_periodic_comp = UF->primary_grid()->get_periodic_directions();
   is_periodic[1][0] = U_periodic_comp->operator()( 0 );
   is_periodic[1][1] = U_periodic_comp->operator()( 1 );
   if(dim >2)
      is_periodic[1][2] = U_periodic_comp->operator()( 2 );

   // Periodic boundary condition check for pressure
   P_periodic_comp = PF->primary_grid()->get_periodic_directions();
   is_periodic[0][0] = P_periodic_comp->operator()( 0 );
   is_periodic[0][1] = P_periodic_comp->operator()( 1 );
   if(dim >2)
      is_periodic[0][2] = P_periodic_comp->operator()( 2 );

   // Create structure to input in the solver system
   struct NavierStokes2System inputData;
   inputData.is_solids_ = is_solids ;
   inputData.is_stressCal_ = is_stressCal ;
   inputData.Npart_ = Npart ;
   inputData.level_set_type_ = level_set_type ;
   inputData.Npoints_ = Npoints ;
   inputData.ar_ = ar ;

   // Build the matrix system
   MAC_ModuleExplorer* se = exp->create_subexplorer( 0,"DDS_NavierStokesSystem" ) ;
   GLOBAL_EQ = DDS_NavierStokesSystem::create( this, se, UF, PF, inputData ) ;
   se->destroy() ;

   // Create Grains3D if solidSolverType is Grains3D;
   int error = 0;
   solidSolver = FS_SolidPlugIn_BuilderFactory:: create( solidSolverType,
	                                            solidSolver_insertionFile,
                                               solidSolver_simulationFile,
                                               1., false,
                                               b_particles_as_fixed_obstacles,
                                               1., b_solidSolver_parallel,
                                               error );

   // Timing routines
   if ( my_rank == is_master )
   {
     SCT_insert_app("Matrix_Assembly&Initialization");
     SCT_insert_app("Pressure predictor");
     SCT_insert_app("Velocity update");
     SCT_insert_app("Penalty Step");
     SCT_insert_app("Pressure Update");
     SCT_get_elapsed_time("Objects_Creation");
   }
}

//---------------------------------------------------------------------------
DDS_NavierStokes:: ~DDS_NavierStokes( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: ~DDS_NavierStokes" ) ;

   free_DDS_subcommunicators() ;

   if ( solidSolver ) delete solidSolver;
   if ( solidFluid_transferStream ) delete solidFluid_transferStream;
   if ( allrigidbodies ) delete allrigidbodies;

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_one_inner_iteration( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_one_inner_iteration" ) ;
   MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "DDS_NavierStokes:: do_one_inner_iteration" ) ;
   start_solving_timer() ;

   if ( my_rank == is_master ) SCT_set_start("Pressure predictor");
   NS_first_step(t_it);
   if ( my_rank == is_master ) SCT_get_elapsed_time( "Pressure predictor" );

   if ( my_rank == is_master ) SCT_set_start( "Velocity update" );
   NS_velocity_update(t_it);
   if ( my_rank == is_master ) SCT_get_elapsed_time( "Velocity update" );

   if ( my_rank == is_master ) SCT_set_start( "Penalty Step" );
   NS_pressure_update(t_it);
   if ( my_rank == is_master ) SCT_get_elapsed_time( "Penalty Step" );

   if ( my_rank == is_master ) SCT_set_start( "Pressure Update" );
   NS_final_step(t_it);
   if ( my_rank == is_master ) SCT_get_elapsed_time( "Pressure Update" );

   UF->copy_DOFs_value( 0, 1 );

   stop_solving_timer() ;
   stop_total_timer() ;

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_before_time_stepping( FV_TimeIterator const* t_it,
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_before_time_stepping" ) ;

   start_total_timer( "DDS_NavierStokes:: do_before_time_stepping" ) ;

   if ( my_rank == is_master ) SCT_set_start("Matrix_Assembly&Initialization");

   FV_OneStepIteration::do_before_time_stepping( t_it, basename ) ;

   allocate_mpi_variables(PF);
   allocate_mpi_variables(UF);

   // Initialize velocity vector at the matrix level
   GLOBAL_EQ->initialize_DS_velocity();
   GLOBAL_EQ->initialize_DS_pressure();

   // Setting ugradu as zero at start of simulation
   if (b_restart == false) ugradu_initialization ( );

   // Calculate row index for each field in each direction`
   calculate_row_indexes ( PF );
   calculate_row_indexes ( UF );

   // Projection-Translation
   if ( b_projection_translation )
   {
     set_translation_vector() ;

     if ( MVQ_translation_vector(translation_direction) < 0. )
       bottom_coordinate = (*UF->primary_grid()->get_global_main_coordinates())
                [translation_direction](0) ;
     else
       bottom_coordinate = (*UF->primary_grid()->get_global_main_coordinates())
           [translation_direction]((*UF->primary_grid()->get_global_max_index())
                                (translation_direction)) ;

     build_links_translation() ;
   }

   // Generate solid particles if required
   if (is_solids) {

      solidFluid_transferStream = NULL;
      solidSolver->getSolidBodyFeatures( solidFluid_transferStream );

      allrigidbodies = new DS_AllRigidBodies( dim,
      	*solidFluid_transferStream, b_particles_as_fixed_obstacles, UF, PF );

      Solids_generation( );

      if (my_rank == 0)
         cout << "Finished particle generation... \n" << endl;

      if (my_rank == 0)
         cout << "Finished intersection calculations... \n" << endl;

      vector<size_t> vec{ 0, 1, 3};
      if (dim == 3) vec.push_back(4);

      initialize_grid_nodes_on_rigidbody(vec);

      if (my_rank == 0)
         cout << "Finished field initializations... \n" << endl;

      if (is_stressCal) {
         // Generate discretization of surface in approximate equal area
         generate_surface_discretization ();
      }
      if (my_rank == 0)
         cout << "Finished particle surface discretizations... \n" << endl;
   }

   // Direction splitting
   // Assemble 1D tridiagonal matrices
   assemble_1D_matrices(PF,t_it);
   assemble_1D_matrices(UF,t_it);

   if (my_rank == 0)
      cout << "Finished assembling pre-coefficient matrix... \n" << endl;

   if ( my_rank == is_master )
      SCT_get_elapsed_time( "Matrix_Assembly&Initialization" );

   stop_total_timer() ;

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_after_time_stepping" ) ;

   // Elapsed time by sub-problems

   // SCT_set_start( "Writing CSV" );
//   write_output_field(PF,0);
//   write_output_field(UF,1);
   // SCT_get_elapsed_time( "Writing CSV" );

   output_L2norm_velocity(0);
   output_L2norm_pressure(0);
// //   error_with_analytical_solution_poiseuille();
// //   error_with_analytical_solution_couette(PF,0);
// //   error_with_analytical_solution_couette(UF,1);

   if ( my_rank == is_master )
   {
     double cputime = CT_get_elapsed_time();
     cout << endl << "Full problem" << endl;
     write_elapsed_time_smhd(cout,cputime,"Computation time");
     SCT_get_summary(cout,cputime);
   }

   deallocate_mpi_variables();
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_before_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_before_inner_iterations_stage" ) ;

   start_total_timer( "DDS_NavierStokes:: do_before_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_before_inner_iterations_stage( t_it ) ;

   //  Projection_Translation
   if ( b_projection_translation )
     b_grid_has_been_translated_at_previous_time = false;

   // if ((is_par_motion) && (is_solids)) {
   //    update_particle_system(t_it);
   //    node_property_calculation(PF);
   //    node_property_calculation(UF);
   //    nodes_field_initialization(0);
   //    nodes_field_initialization(1);
   //    nodes_field_initialization(3);
   //    if (dim == 3) nodes_field_initialization(4);
   //
   //    // Direction splitting
   //    // Assemble 1D tridiagonal matrices
   //    assemble_1D_matrices(PF,t_it);
   //    assemble_1D_matrices(UF,t_it);
   // }

   // Perform matrix level operations before each time step
   GLOBAL_EQ->at_each_time_step( );

   stop_total_timer() ;

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_after_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_after_inner_iterations_stage" ) ;

   start_total_timer( "DDS_NavierStokes:: do_after_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_after_inner_iterations_stage( t_it ) ;

   // Compute velocity change over the time step
   double velocity_time_change = GLOBAL_EQ->compute_DS_velocity_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master ) cout << "velocity change = " <<
     	MAC::doubleToString( ios::scientific, 5, velocity_time_change ) << endl;

   // Compute hydrodynamic forces by surface viscous stress
   // if (is_stressCal) {
   //    compute_fluid_velocity_particle_interaction(t_it);
   // }
   //
   // if (level_set_type == "PipeX") error_with_analytical_solution_poiseuille3D();
/*
   double vel_divergence = get_velocity_divergence(t_it);

   if ( my_rank == is_master ) {
      string fileName = "./DS_results/max_divergence.csv" ;
      ofstream MyFile( fileName.c_str(), ios::app ) ;
      MyFile << t_it -> time() << "," << vel_divergence << endl;
      MyFile.close( ) ;
   }
*/
//   write_output_field(UF,t_it);

   double cfl = UF->compute_CFL( t_it, 0 );
   if ( my_rank == is_master )
      MAC::out() << "CFL: "<< cfl <<endl;

   // Projection translation
   if ( b_projection_translation ) {

      PartInput solid = GLOBAL_EQ->get_solid(0);

      double distance_to_bottom =
       MAC::abs(solid.coord[translation_direction]->item(0)-bottom_coordinate);

      if ( distance_to_bottom < critical_distance_translation ) {

         if ( my_rank == is_master )
            MAC::out() << "         -> -> -> -> -> -> -> -> -> -> -> ->"
               << endl << "         !!!     Domain Translation      !!!"
               << endl << "         -> -> -> -> -> -> -> -> -> -> -> ->"
               << endl;

         b_grid_has_been_translated_at_previous_time = true;

         translated_distance += MVQ_translation_vector( translation_direction );
         if ( my_rank == is_master )
            MAC::out() << "         Translated distance = " <<
                translated_distance << endl;

         fields_projection();

         if ( MVQ_translation_vector(translation_direction) < 0. )
            bottom_coordinate = (*UF->primary_grid()->get_global_main_coordinates())
                                [translation_direction](0) ;
         else
            bottom_coordinate = (*UF->primary_grid()->get_global_main_coordinates())
               [translation_direction]((*UF->primary_grid()->get_global_max_index())
                                (translation_direction)) ;

         b_grid_has_been_translated_since_last_output = true;
      }
   }

   stop_total_timer() ;

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_additional_savings" ) ;

   start_total_timer( "DDS_NavierStokes:: do_additional_savings" ) ;

   stop_total_timer() ;

   GLOBAL_EQ->display_debug();
}

// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: error_with_analytical_solution_couette (FV_DiscreteField const* FF, size_t const& field)
// //---------------------------------------------------------------------------
// {
//    MAC_LABEL( "DDS_NavierStokes:: error_with_analytical_solution_couette" ) ;
//
//    size_t_vector min_unknown_index(dim,0);
//    size_t_vector max_unknown_index(dim,0);
//
//    PartInput solid = GLOBAL_EQ->get_solid(0);
//    NodeProp node = GLOBAL_EQ->get_node_property(field,0);
//
//    for (size_t comp=0;comp<nb_comps[field];++comp) {
//       // Compute error
//       double computed_field = 0., analytical_solution = 0., error_L2 = 0.;
//       // Get local min and max indices
//       for (size_t l=0;l<dim;++l) {
//          min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
//          max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
//       }
//
//       double xp = solid.coord[0]->item(0);
//       double yp = solid.coord[1]->item(0);
//       double r1 = solid.size->item(0);
//       double o1 = solid.ang_vel[2]->item(0);
//       double r2 = solid.size->item(1);
//       double o2 = solid.ang_vel[2]->item(1);
//       double a = (o2*pow(r2,2.)-o1*pow(r1,2.))/(pow(r2,2.)-pow(r1,2.));
//       double b = (o1-o2)*(pow(r1,2.)*pow(r2,2.))/(pow(r2,2.)-pow(r1,2.));
//       //double mean_press = 1./(r2-r1)*(pow(a,2.)*(pow(r2,3.)-pow(r1,3.))/6. + 2.*a*b*((r2*log(r2)-r2)-(r1*log(r1)-r1)) + pow(b,2.)/2./r2 - pow(b,2.)/2./r1);
//       double mean_press = 1./(pow(r2,2)-pow(r1,2))*(pow(a,2.)/4.*(pow(r2,4.)-pow(r1,4.)) + a*b*(2.*pow(r2,2.)*log(r2)-pow(r2,2.)) - a*b*(2.*pow(r1,2.)*log(r1)-pow(r1,2.)) - pow(b,2.)*log(r2/r1));
//
//       for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
//          double x = FF->get_DOF_coordinate( i, comp, 0 ) ;
//          for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
//             double y = FF->get_DOF_coordinate( j, comp, 1 ) ;
//             if (dim == 2) {
//                size_t k=0;
//                computed_field = FF->DOF_value( i, j, k, comp, 0 ) ;
//
//                size_t p = return_node_index(FF,comp,i,j,k);
//                if (node.void_frac[comp]->item(p) != 1.) {
//                   double r = pow(pow(x-xp,2.)+pow(y-yp,2.),0.5);
//                   if (field == 1) {
//                      double v_theta = a*r + b/r;
//                      if (comp == 0) {
//                         analytical_solution = -v_theta*(y-yp)/r;
// //                        analytical_solution = sin(x)*sin(y);
//                      } else if (comp == 1) {
//                         analytical_solution = v_theta*(x-xp)/r;
// //                        analytical_solution = cos(x)*cos(y);
//                      }
//                   } else if (field == 0) {
// //                     analytical_solution = 0.;
//                      analytical_solution = rho*(pow(a,2.)*pow(r,2.)/2. + 2.*a*b*log(r) - pow(b,2.)/2./pow(r,2.));
//                      analytical_solution -= rho*mean_press;
// //                     analytical_solution = sin(x+y);
//                   }
//
//                   if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
// //                     error_L2 += MAC::sqr( analytical_solution)
//                      error_L2 += MAC::sqr( computed_field - analytical_solution)
//                                             * FF->get_cell_measure( i, j, k, comp ) ;
//                   }
//                }
//             }
//          }
//       }
//
//       error_L2 = macCOMM->sum( error_L2 );
//       error_L2 = MAC::sqrt(error_L2);
//
//       if ( my_rank == 0 )
//          cout << "L2 Error with analytical solution field " << FF->name() << ", component " << comp << " = " << std::fixed << std::setprecision(16) << error_L2 << endl;
//    }
// }

// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: error_with_analytical_solution_poiseuille3D ( )
// //---------------------------------------------------------------------------
// {
//    MAC_LABEL("DDS_NavierStokes:: error_with_analytical_solution_poiseuille3D" ) ;
//
//    // Parameters
//    size_t cpp=10;
//    double bodyterm=0.;
//
//    // Periodic pressure gradient
//    if ( UF->primary_grid()->is_periodic_flow() ) {
//       cpp = UF->primary_grid()->get_periodic_flow_direction() ;
//       bodyterm = UF->primary_grid()->get_periodic_pressure_drop() /
//                ( UF->primary_grid()->get_main_domain_max_coordinate( cpp )
//                - UF->primary_grid()->get_main_domain_min_coordinate( cpp ) ) ;
//    }
//
//    for (size_t comp=0;comp<nb_comps[1];++comp) {
//       // Get nb of local dof
//       size_t_vector local_dof_number( dim, 0 );
//       for (size_t l=0;l<dim;++l)
//          local_dof_number(l) = UF->get_local_nb_dof( comp, l ) ;
//
//       double computed_field = 0., analytical_solution = 0.;
//       double error_L2 = 0.;
//
//       if (is_solids) {
//          PartInput solid = GLOBAL_EQ->get_solid(0);
//          NodeProp node = GLOBAL_EQ->get_node_property(1,0);
//          double yp = solid.coord[1]->item(0);
//          double zp = solid.coord[2]->item(0);
//          double rp = solid.size->item(0);
//
//          // Compute error
//          for (size_t i=0;i<local_dof_number(0);++i) {
//             for (size_t j=0;j<local_dof_number(1);++j) {
//                double y = UF->get_DOF_coordinate( j, comp, 1 ) - yp ;
//                for (size_t k=0;k<local_dof_number(2);++k) {
//                   double z = UF->get_DOF_coordinate( k, comp, 2 ) - zp ;
//                   computed_field = UF->DOF_value( i, j, k, comp, 0 ) ;
//
//                   double ri = pow(pow(y,2)+pow(z,2),0.5);
//
//                   size_t p = return_node_index(UF,comp,i,j,k);
//                   if (node.void_frac[comp]->item(p) == 0.) {
//                      if (comp == cpp) {
//                         analytical_solution = (bodyterm/(4.*mu))*(ri*ri-rp*rp);
//                      }
//                      if ( UF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) )
//                         error_L2 += MAC::sqr( computed_field - analytical_solution ) * UF->get_cell_measure( i, j, k, comp ) ;
//                   }
//                }
//             }
//          }
//       }
//
//       error_L2 = macCOMM->sum( error_L2 );
//       error_L2 = MAC::sqrt(error_L2);
//
//       if ( my_rank == 0 )
//          cout << "L2 Error with analytical solution field " << UF->name() << ", component " << comp << " = " << std::fixed << std::setprecision(16) << error_L2 << endl;
//    }
//
// }
//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: error_with_analytical_solution_poiseuille ( )
// //---------------------------------------------------------------------------
// {
//    MAC_LABEL(
//    "DDS_NavierStokes:: error_with_analytical_solution_poiseuille" ) ;
//
//    // Parameters
//    size_t cpp=10;
//    double bodyterm=0., height;
//
//    // Periodic pressure gradient
//    if ( UF->primary_grid()->is_periodic_flow() ) {
//       cpp = UF->primary_grid()->get_periodic_flow_direction() ;
//       bodyterm = UF->primary_grid()->get_periodic_pressure_drop() /
//                ( UF->primary_grid()->get_main_domain_max_coordinate( cpp )
//                - UF->primary_grid()->get_main_domain_min_coordinate( cpp ) ) ;
//    }
//
//    height = ( UF->primary_grid()->get_main_domain_max_coordinate(1) - UF->primary_grid()->get_main_domain_min_coordinate(1) ) /2. ;
//
//    for (size_t comp=0;comp<nb_comps[1];++comp) {
//       // Get nb of local dof
//       size_t_vector local_dof_number( dim, 0 );
//       for (size_t l=0;l<dim;++l)
//          local_dof_number(l) = UF->get_local_nb_dof( comp, l ) ;
//
//       // Compute error
//       double computed_field = 0., analytical_solution = 0.;
//       double error_L2 = 0.;
//       for (size_t i=0;i<local_dof_number(0);++i) {
//          for (size_t j=0;j<local_dof_number(1);++j) {
//             double y = UF->get_DOF_coordinate( j, comp, 1 ) ;
//
//             if ( dim == 2 ) {
//                size_t k = 0 ;
//                computed_field = UF->DOF_value( i, j, k, comp, 0 ) ;
//
//                if (comp == cpp) {
//                   analytical_solution = (bodyterm/(2.*mu))*((y-height)*(y-height)-height*height);
//                }
//
// //               analytical_solution = MAC::sin(MAC::pi()*x)*MAC::sin(MAC::pi()*y);
//
//                if ( UF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) )
// 	          error_L2 += MAC::sqr( computed_field - analytical_solution ) * UF->get_cell_measure( i, j, k, comp ) ;
// 	    } else {
//                for (size_t k=0;k<local_dof_number(2);++k) {
//                   computed_field = UF->DOF_value( i, j, k, comp, 0 ) ;
//
//                   if (comp == cpp) {
//                      analytical_solution = (bodyterm/(2.*mu))*(y*y-height*height);
//                   }
//
//                   if ( UF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) )
//                      error_L2 += MAC::sqr( computed_field - analytical_solution ) * UF->get_cell_measure( i, j, k, comp ) ;
// 	       }
// 	    }
//          }
//       }
//
//       error_L2 = macCOMM->sum( error_L2 );
//       error_L2 = MAC::sqrt(error_L2);
//
//       if ( my_rank == 0 )
//          cout << "L2 Error with analytical solution field " << UF->name() << ", component " << comp << " = " << std::fixed << std::setprecision(16) << error_L2 << endl;
//    }
//
// }

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: ugradu_initialization ( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: ugradu_initialization" ) ;

  for (size_t comp=0;comp<nb_comps[1];comp++) UF->set_DOFs_value( comp, 2, 0.);

}



//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: update_particle_system(FV_TimeIterator const* t_it)
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL( "DDS_NavierStokes:: update_particle_system" ) ;
//
//   // Structure of particle input data
//   PartForces hydro_forces = GLOBAL_EQ->get_forces(0);
//   PartForces hydro_torque = GLOBAL_EQ->get_torque(0);
//
//   doubleVector const& gg = gravity_vector->to_double_vector();
//
//   if (insertion_type == "file") {
//
//      // Structure of particle input data
//      PartInput solid = GLOBAL_EQ->get_solid(0);
//
//      for (size_t i=0;i<Npart;i++) {
//         double rp = solid.size->item(i);
//         double mass_p = (dim == 3) ? rho_s*(4./3.)*MAC::pi()*pow(rp,3.) :
//                                      rho_s*MAC::pi()*pow(rp,2.)*1.;
//         double moi = 0.;
// 	if (level_set_type == "Sphere")
//            moi = (dim == 3) ? (2./5.)*mass_p*rp*rp :
// 		              (1./2.)*mass_p*rp*rp ;
//
//         doubleVector pos(dim,0);
//         doubleVector vel(dim,0);
//         doubleVector acc(dim,0);
//         doubleVector ang_vel(dim,0);
//         doubleVector ang_acc(dim,0);
//
//         for (size_t dir=0;dir<dim;dir++) {
//            pos(dir) = solid.coord[dir]->item(i);
//            vel(dir) = solid.vel[dir]->item(i);
//
//            if (motion_type == "Sine") {
//               vel(dir) = gg(dir)*Amp*MAC::cos(2.*MAC::pi()*freq*t_it->time());
//               pos(dir) = pos(dir) + vel(dir)*t_it->time_step();
//            } else if (motion_type == "Hydro") {
//               acc(dir) = gg(dir)*(1-rho/rho_s) + (hydro_forces.press[dir]->item(i)
// 	                                       +  hydro_forces.vel[dir]->item(i)) / mass_p ;
//               vel(dir) = vel(dir) + acc(dir)*t_it->time_step();
//               pos(dir) = pos(dir) + vel(dir)*t_it->time_step();
//               ang_acc(dir) = (hydro_torque.press[dir]->item(i)
//                            +  hydro_torque.vel[dir]->item(i)) / moi ;
//               ang_vel(dir) = ang_vel(dir) + ang_acc(dir)*t_it->time_step();
//            }
//
//            solid.coord[dir]->set_item(i,pos(dir));
//            solid.vel[dir]->set_item(i,vel(dir));
//            solid.ang_vel[dir]->set_item(i,ang_vel(dir));
//         }
//      }
//   }
// }
//



//---------------------------------------------------------------------------
void
DDS_NavierStokes:: calculate_row_indexes ( FV_DiscreteField const* FF)
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: calculate_row_indexes" ) ;

   size_t field = (FF == PF) ? 0 : 1;

   for (size_t comp = 0; comp < nb_comps[field]; comp++) {
      // Get local min and max indices
      size_t_vector min_unknown_index(3,0);
      size_t_vector max_unknown_index(3,0);

      for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) =
                           FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) =
                           FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (size_t dir = 0; dir < dim; dir++) {
         LA_SeqMatrix* row_index = GLOBAL_EQ->get_row_indexes(field,dir,comp);
         switch (dir) {
            case 0:
             for (size_t j=min_unknown_index(1); j<=max_unknown_index(1); j++) {
              for (size_t k=min_unknown_index(2); k<=max_unknown_index(2); k++) {
                  size_t p = (j - min_unknown_index(1))
                           + (1 + max_unknown_index(1) - min_unknown_index(1))
                           * (k - min_unknown_index(2));
                  row_index->set_item(j,k,(double)p);
              }
             }
             break;
            case 1:
             for (size_t i=min_unknown_index(0); i<=max_unknown_index(0); i++) {
              for (size_t k=min_unknown_index(2); k<=max_unknown_index(2); k++) {
                  size_t p = (i - min_unknown_index(0))
                           + (1 + max_unknown_index(0) - min_unknown_index(0))
                           * (k - min_unknown_index(2));
                  row_index->set_item(i,k,(double)p);
              }
             }
             break;
            case 2:
             for (size_t i=min_unknown_index(0); i<=max_unknown_index(0); i++) {
              for (size_t j=min_unknown_index(1); j<=max_unknown_index(1); j++) {
                  size_t p = (i - min_unknown_index(0))
                           + (1 + max_unknown_index(0) - min_unknown_index(0))
                           * (j - min_unknown_index(1));
                  row_index->set_item(i,j,(double)p);
              }
             }
             break;
         }
      }
   }
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: trans_rotation_matrix (size_t const& m
                                        , class doubleVector& delta
                                        , size_t const& comp
                                        , size_t const& field)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: trans_rotation_matrix" ) ;

  PartInput solid = GLOBAL_EQ->get_solid(0);

  // yaw along z-axis; pitch along y-axis; roll along x-axis
  doubleArray2D rot_matrix(3,3,0);

  // Rotation matrix assemble
  if (insertion_type == "file") {
     double roll = (MAC::pi()/180.)*solid.thetap->item(m,0);
     double pitch = (MAC::pi()/180.)*solid.thetap->item(m,1);
     double yaw = (MAC::pi()/180.)*solid.thetap->item(m,2);
     rot_matrix(0,0) = MAC::cos(yaw)*MAC::cos(pitch);
     rot_matrix(1,0) = MAC::cos(yaw)*MAC::sin(pitch)*MAC::sin(roll)
                                    - MAC::sin(yaw)*MAC::cos(roll);
     rot_matrix(2,0) = MAC::cos(yaw)*MAC::sin(pitch)*MAC::cos(roll)
                                    + MAC::sin(yaw)*MAC::sin(roll);
     rot_matrix(0,1) = MAC::sin(yaw)*MAC::cos(pitch);
     rot_matrix(1,1) = MAC::sin(yaw)*MAC::sin(pitch)*MAC::sin(roll)
                                    + MAC::cos(yaw)*MAC::cos(roll);
     rot_matrix(2,1) = MAC::sin(yaw)*MAC::sin(pitch)*MAC::cos(roll)
                                    - MAC::cos(yaw)*MAC::sin(roll);
     rot_matrix(0,2) = -MAC::sin(pitch);
     rot_matrix(1,2) = MAC::cos(pitch)*MAC::sin(roll);
     rot_matrix(2,2) = MAC::cos(pitch)*MAC::cos(roll);
  } else if (insertion_type == "Grains3D") {
     rot_matrix(0,0) = solid.thetap->item(m,0);
     rot_matrix(1,0) = solid.thetap->item(m,1);
     rot_matrix(2,0) = solid.thetap->item(m,2);
     rot_matrix(0,1) = solid.thetap->item(m,3);
     rot_matrix(1,1) = solid.thetap->item(m,4);
     rot_matrix(2,1) = solid.thetap->item(m,5);
     rot_matrix(0,2) = solid.thetap->item(m,6);
     rot_matrix(1,2) = solid.thetap->item(m,7);
     rot_matrix(2,2) = solid.thetap->item(m,8);
  }

  double delta_x = delta(0)*rot_matrix(0,0)
                 + delta(1)*rot_matrix(0,1)
                 + delta(2)*rot_matrix(0,2);
  double delta_y = delta(0)*rot_matrix(1,0)
                 + delta(1)*rot_matrix(1,1)
                 + delta(2)*rot_matrix(1,2);
  double delta_z = delta(0)*rot_matrix(2,0)
                 + delta(1)*rot_matrix(2,1)
                 + delta(2)*rot_matrix(2,2);

  delta(0) = delta_x;
  delta(1) = delta_y;
  delta(2) = delta_z;
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: rotation_matrix (size_t const& m
                                  , class doubleVector& delta
                                  , size_t const& comp
                                  , size_t const& field)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: rotation_matrix" ) ;

  PartInput solid = GLOBAL_EQ->get_solid(0);

  // yaw along z-axis; pitch along y-axis; roll along x-axis
  doubleArray2D rot_matrix(3,3,0);

  // Rotation matrix assemble
  if (insertion_type == "file") {
     double roll = (MAC::pi()/180.)*solid.thetap->item(m,0);
     double pitch = (MAC::pi()/180.)*solid.thetap->item(m,1);
     double yaw = (MAC::pi()/180.)*solid.thetap->item(m,2);
     rot_matrix(0,0) = MAC::cos(yaw)*MAC::cos(pitch);
     rot_matrix(0,1) = MAC::cos(yaw)*MAC::sin(pitch)*MAC::sin(roll)
                                    - MAC::sin(yaw)*MAC::cos(roll);
     rot_matrix(0,2) = MAC::cos(yaw)*MAC::sin(pitch)*MAC::cos(roll)
                                    + MAC::sin(yaw)*MAC::sin(roll);
     rot_matrix(1,0) = MAC::sin(yaw)*MAC::cos(pitch);
     rot_matrix(1,1) = MAC::sin(yaw)*MAC::sin(pitch)*MAC::sin(roll)
                                    + MAC::cos(yaw)*MAC::cos(roll);
     rot_matrix(1,2) = MAC::sin(yaw)*MAC::sin(pitch)*MAC::cos(roll)
                                    - MAC::cos(yaw)*MAC::sin(roll);
     rot_matrix(2,0) = -MAC::sin(pitch);
     rot_matrix(2,1) = MAC::cos(pitch)*MAC::sin(roll);
     rot_matrix(2,2) = MAC::cos(pitch)*MAC::cos(roll);
  } else if (insertion_type == "Grains3D") {
     rot_matrix(0,0) = solid.thetap->item(m,0);
     rot_matrix(0,1) = solid.thetap->item(m,1);
     rot_matrix(0,2) = solid.thetap->item(m,2);
     rot_matrix(1,0) = solid.thetap->item(m,3);
     rot_matrix(1,1) = solid.thetap->item(m,4);
     rot_matrix(1,2) = solid.thetap->item(m,5);
     rot_matrix(2,0) = solid.thetap->item(m,6);
     rot_matrix(2,1) = solid.thetap->item(m,7);
     rot_matrix(2,2) = solid.thetap->item(m,8);
  }

  double delta_x = delta(0)*rot_matrix(0,0)
                 + delta(1)*rot_matrix(0,1)
                 + delta(2)*rot_matrix(0,2);
  double delta_y = delta(0)*rot_matrix(1,0)
                 + delta(1)*rot_matrix(1,1)
                 + delta(2)*rot_matrix(1,2);
  double delta_z = delta(0)*rot_matrix(2,0)
                 + delta(1)*rot_matrix(2,1)
                 + delta(2)*rot_matrix(2,2);

  delta(0) = delta_x;
  delta(1) = delta_y;
  delta(2) = delta_z;
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: Solids_generation ()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: Solids_generation" ) ;

  // Structure of particle input data
  PartInput solid = GLOBAL_EQ->get_solid(0);

  double xp,yp,zp,Rp,tx,ty,tz,vx,vy,vz,wx,wy,wz,Tp,off;

  ifstream inFile;
  std::ostringstream os2;
  os2 << "./InputFiles/" << solid_filename;
  std::string filename = os2.str();

  inFile.open(filename.c_str());
  string line;
  getline(inFile,line);
  for (size_t i=0;i<Npart;i++) {
     inFile >> xp >> yp >> zp >> Rp >> tx >> ty >> tz >> vx >> vy >> vz >> wx >> wy >> wz >> Tp >> off;
     solid.coord[0]->set_item(i,xp);
     solid.coord[1]->set_item(i,yp);
     solid.coord[2]->set_item(i,zp);
     solid.size->set_item(i,Rp);
     solid.thetap->set_item(i,0,tx);
     solid.thetap->set_item(i,1,ty);
     solid.thetap->set_item(i,2,tz);
     solid.vel[0]->set_item(i,vx);
     solid.vel[1]->set_item(i,vy);
     solid.vel[2]->set_item(i,vz);
     solid.ang_vel[0]->set_item(i,wx);
     solid.ang_vel[1]->set_item(i,wy);
     solid.ang_vel[2]->set_item(i,wz);
     solid.temp->set_item(i,Tp);
     solid.inside->set_item(i,off);
  }
  inFile.close();
}




//---------------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_field_matrix ( FV_DiscreteField const* FF,
                                           FV_TimeIterator const* t_it,
                                           double const& gamma,
                                           size_t const& comp,
                                           size_t const& dir,
                                           size_t const& j,
                                           size_t const& k,
                                           size_t const& r_index )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: assemble_field_matrix" ) ;

   // Parameters
   double dxr,dxl,dx,xR,xL,xC,right=0.,left=0.,center=0.;

   size_t field = (FF == PF) ? 0 : 1 ;

   bool r_bound = false;
   bool l_bound = false;
   // All the proc will have open right bound,
   // except last proc for non periodic systems
   if ((is_periodic[field][dir] != 1)
    && (rank_in_i[dir] == nb_ranks_comm_i[dir]-1))
      r_bound = true;
   // All the proc will have open left bound,
   // except first proc for non periodic systems
   if ((is_periodic[field][dir] != 1)
    && (rank_in_i[dir] == 0))
      l_bound = true;

   // Get local min and max indices
   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);

   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) =
                     FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      max_unknown_index(l) =
                     FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   // Perform assembling
   size_t m, i;
   TDMatrix* A = GLOBAL_EQ-> get_A(field);

   double Aee_diagcoef=0.;

   for (m=0,i=min_unknown_index(dir);i<=max_unknown_index(dir);++i,++m) {
      int i_temp=0; int j_temp=0; int k_temp=0;
      if (dir == 0) {
         i_temp = (int) i;
         j_temp = (int) min_unknown_index(1);
         k_temp = (int) min_unknown_index(2);
      } else if (dir == 1) {
         i_temp = (int) min_unknown_index(0);
         j_temp = (int) i;
         k_temp = (int) min_unknown_index(2);
      } else if (dir == 2) {
         i_temp = (int) min_unknown_index(0);
         j_temp = (int) min_unknown_index(1);
         k_temp = (int) i;
      }

      xC= FF->get_DOF_coordinate( i, comp, dir) ;

      // Check if the index is at right domain
      // boundary with neumann or dirichlet BC
      if ( (i==max_unknown_index(dir))
        && r_bound
        && FF->DOF_on_BC(i_temp,j_temp,k_temp,comp)) {
         xR = 0.;
      } else {
         xR= FF->get_DOF_coordinate( i+1,comp, dir) ;
      }

      // Check if the index is at left domain
      // boundary with neumann or dirichlet BC
      if ( (i==min_unknown_index(dir))
        && l_bound
        && FF->DOF_on_BC(i_temp,j_temp,k_temp,comp)) {
         xL = 0.;
      } else {
         xL= FF->get_DOF_coordinate( i-1, comp, dir) ;
      }

      dx = FF->get_cell_size( i,comp, dir);

      dxr= xR - xC;
      dxl= xC - xL;

      double value=0., unsteady_term=0.;

      if (field == 0) {
         right = -1.0/(dxr);
         left = -1.0/(dxl);
         // add unsteady term for pressure field
         unsteady_term = 1.0*dx;
      } else if (field == 1) {
         right = -gamma/(dxr);
         left = -gamma/(dxl);
         // add unsteady term for velocity field
         unsteady_term = rho*FF->get_cell_size(i,comp,dir)/t_it->time_step();
      }

      if ((is_solids) && (field == 1)) {
         size_t p=0;
         if (dir == 0) {
            p = FF->DOF_local_number(i,j,k,comp);
         } else if (dir == 1) {
            p = FF->DOF_local_number(j,i,k,comp);
         } else if (dir == 2) {
            p = FF->DOF_local_number(j,k,i,comp);
         }

         size_t_array2D* intersect_vector =
                           allrigidbodies->get_intersect_vector_on_grid(FF);
         doubleArray2D* intersect_distance =
                           allrigidbodies->get_intersect_distance_on_grid(FF);
         size_t_vector* void_frac =
                           allrigidbodies->get_void_fraction_on_grid(FF);

         if (void_frac->operator()(p) == 0) {
            // if left node is inside the solid particle
            if (intersect_vector->operator()(p,2*dir+0) == 1) {
               left = -gamma/intersect_distance->operator()(p,2*dir+0);
            }
            // if right node is inside the solid particle
            if (intersect_vector->operator()(p,2*dir+1) == 1) {
               right = -gamma/intersect_distance->operator()(p,2*dir+1);
            }
         } else if (void_frac->operator()(p) == 1) {
            // if center node is inside the solid particle
            left = 0.;
            right = 0.;
         }

         center = -(right+left);

         if (intersect_vector->operator()(p,2*dir+0) == 1) left = 0.;
         if (intersect_vector->operator()(p,2*dir+1) == 1) right = 0.;
      } else {
         center = - (right+left);
      }

      // Since, this function is used in all directions;
      // ii, jj, and kk are used to convert the passed
      // arguments corresponding to correct direction
      int ii=0,jj=0,kk=0;

      // Condition for handling the pressure neumann conditions at wall
      if (i == min_unknown_index(dir) && l_bound) {
         if (dir == 0) {
            ii = (int) i-1;
            jj = (int) min_unknown_index(1);
            kk = (int) min_unknown_index(2);
         } else if (dir == 1) {
            ii = (int) min_unknown_index(0);
            jj = (int) i-1;
            kk = (int) min_unknown_index(2);
         } else if (dir == 2) {
            ii = (int) min_unknown_index(0);
            jj = (int) min_unknown_index(1);
            kk = (int) i-1;
         }
         if (FF->DOF_in_domain(ii,jj,kk,comp)
          && FF->DOF_has_imposed_Dirichlet_value
                        ((size_t)ii,(size_t)jj,(size_t)kk,comp)) {
            // For Dirichlet boundary condition
            value = center;
         } else {
            // For Neumann homogeneous boundary condition
            value = -right;
         }
      } else if (i==max_unknown_index(dir) && r_bound) {
         if (dir == 0) {
            ii = (int) i+1;
            jj = (int) max_unknown_index(1);
            kk = (int) max_unknown_index(2);
         } else if (dir == 1) {
            ii = (int) max_unknown_index(0);
            jj = (int) i+1;
            kk = (int) max_unknown_index(2);
         } else if (dir == 2) {
            ii = (int) max_unknown_index(0);
            jj = (int) max_unknown_index(1);
            kk = (int) i+1;
         }
         if (FF->DOF_in_domain(ii,jj,kk,comp)
          && FF->DOF_has_imposed_Dirichlet_value
                        ((size_t)ii,(size_t)jj,(size_t)kk,comp)) {
            // For Dirichlet boundary condition
            value = center;
         } else {
            // For Neumann homogeneous boundary condition
            value = -left;
         }
      } else {
         value = center;
      }

      value = value + unsteady_term;

      // Set Aie, Aei and Ae
      if ((!l_bound) && (i == min_unknown_index(dir))) {
         // Periodic boundary condition at minimum unknown index
         // First proc has non zero value in Aie,Aei for first & last index
         if (rank_in_i[dir] == 0) {
            A[dir].ie[comp][r_index]->set_item(m,nb_ranks_comm_i[dir]-1,left);
            A[dir].ei[comp][r_index]->set_item(nb_ranks_comm_i[dir]-1,m,left);
         } else {
            A[dir].ie[comp][r_index]->set_item(m,rank_in_i[dir]-1,left);
            A[dir].ei[comp][r_index]->set_item(rank_in_i[dir]-1,m,left);
         }
      }

      if ((!r_bound) && (i == max_unknown_index(dir))) {
         // Periodic boundary condition at maximum unknown index
         // For last index, Aee comes from this proc as it
         // is interface unknown wrt this proc
         A[dir].ie[comp][r_index]->set_item(m-1,rank_in_i[dir],left);
         Aee_diagcoef = value;
         A[dir].ei[comp][r_index]->set_item(rank_in_i[dir],m-1,left);
      }

      // Set Aii_sub_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1)
       && (is_periodic[field][dir] != 1)) {
         if (i > min_unknown_index(dir))
            A[dir].ii_sub[comp][r_index]->set_item(m-1,left);
      } else {
         if (i < max_unknown_index(dir))
            if (i > min_unknown_index(dir))
               A[dir].ii_sub[comp][r_index]->set_item(m-1,left);
      }

      // Set Aii_super_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1)
       && (is_periodic[field][dir] != 1)) {
         if (i < max_unknown_index(dir))
            A[dir].ii_super[comp][r_index]->set_item(m,right);
      } else {
         if (i < max_unknown_index(dir)-1)
            A[dir].ii_super[comp][r_index]->set_item(m,right);
      }

      // Set Aii_main_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1)
       && (is_periodic[field][dir] != 1)) {
         A[dir].ii_main[comp][r_index]->set_item(m,value);
      } else {
         if (i<max_unknown_index(dir))
            A[dir].ii_main[comp][r_index]->set_item(m,value);
      }
   } // End of for loop

   GLOBAL_EQ->pre_thomas_treatment(comp,dir,A,r_index);

   return (Aee_diagcoef);
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_field_schur_matrix (size_t const& comp
                                              , size_t const& dir
                                              , double const& Aee_diagcoef
                                              , size_t const& field
                                              , size_t const& r_index )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: assemble_field_schur_matrix" ) ;
   // Compute the product matrix for each proc

   TDMatrix* A = GLOBAL_EQ-> get_A(field);

   if (nb_ranks_comm_i[dir]>1) {

      ProdMatrix* Ap = GLOBAL_EQ->get_Ap(field);
      ProdMatrix* Ap_proc0 = GLOBAL_EQ->get_Ap_proc0(field);

      GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir,field,r_index);

      LA_SeqMatrix* product_matrix = Ap[dir].ei_ii_ie[comp];

      size_t nbrow = product_matrix->nb_rows();
      // Create a copy of product matrix to receive matrix,
      // this will eliminate the memory leak issue
      // which caused by "create_copy" command
      for (size_t k=0;k<nbrow;k++) {
         for (size_t j=0;j<nbrow;j++) {
            Ap_proc0[dir].ei_ii_ie[comp]
                              ->set_item(k,j,product_matrix->item(k,j));
         }
      }

      LA_SeqMatrix* receive_matrix = Ap_proc0[dir].ei_ii_ie[comp];

      if ( rank_in_i[dir] == 0 ) {
         A[dir].ee[comp][r_index]->set_item(0,0,Aee_diagcoef);
         for (size_t i=1;i<(size_t)nb_ranks_comm_i[dir];++i) {

            // Create the container to receive
            size_t nbrows = product_matrix->nb_rows();
            size_t nb_received_data = (size_t)pow(nbrows,2)+1;
            double * received_data = new double [nb_received_data];

            // Receive the data
            static MPI_Status status ;
            MPI_Recv( received_data, (int)nb_received_data,
                     MPI_DOUBLE, (int)i, 0, DDS_Comm_i[dir], &status ) ;

            // Transfer the received data to the receive matrix
            for (int k=0;k<(int)nbrows;k++) {
               for (int j=0;j<(int)nbrows;j++) {
                  // Assemble the global product matrix by
                  // adding contributions from all the procs
                  receive_matrix->add_to_item(k,j,received_data[k*(nbrows)+j]);
               }
            }

  	         if (is_periodic[field][dir] == 0) {
               if (i<(size_t)nb_ranks_comm_i[dir]-1) {
                  // Assemble the global Aee matrix
                  // No periodic condition in x.
                  // So no fe contribution from last proc
                  A[dir].ee[comp][r_index]
                           ->set_item(i,i,received_data[nb_received_data-1]);
               }
            } else {
               // Assemble the global Aee matrix
               // Periodic condition in x.
               // So there is fe contribution from last proc
               A[dir].ee[comp][r_index]
                           ->set_item(i,i,received_data[nb_received_data-1]);
            }
            delete [] received_data;
         }
      } else {
         // Create the packed data container
         size_t nbrows = product_matrix->nb_rows();
         size_t nb_send_data = (size_t)pow(nbrows,2)+1;
         double * packed_data = new double [nb_send_data];

         // Fill the packed data container with Aie
         // Iterator only fetches the values present. Zeros are not fetched.

         for (size_t i=0 ; i<nbrows ; i++ ) {
            for ( size_t j=0 ; j<nbrows ; j++ ) {
               // Packing rule
               // Pack the product matrix into a vector
               packed_data[i*nbrows+j]=product_matrix->item(i,j);
            }
         }

         // Fill the last element of packed data
         // with the diagonal coefficient Aee
         packed_data[nb_send_data-1] = Aee_diagcoef;

         // Send the data
         MPI_Send( packed_data, (int)nb_send_data,
                                    MPI_DOUBLE, 0, 0, DDS_Comm_i[dir] ) ;

         delete [] packed_data;
      }

      // Assemble the schlur complement in the master proc

      if (rank_in_i[dir] == 0){
         TDMatrix* Schur = GLOBAL_EQ-> get_Schur(field);
         size_t nb_row = Schur[dir].ii_main[comp][r_index]->nb_rows();
         for (int p = 0; p < (int)nb_row; p++) {
            Schur[dir].ii_main[comp][r_index]
                              ->set_item(p,A[dir].ee[comp][r_index]->item(p,p)
                                          -receive_matrix->item(p,p));
            if (p < (int)nb_row-1)
               Schur[dir].ii_super[comp][r_index]
                              ->set_item(p,-receive_matrix->item(p,p+1));
            if (p > 0)
               Schur[dir].ii_sub[comp][r_index]
                              ->set_item(p-1,-receive_matrix->item(p,p-1));
            // In case of periodic and multi-processor,
            // there will be a variant of Tridiagonal matrix
            // instead of normal format
            if (is_periodic[field][dir] == 1) {
               Schur[dir].ie[comp][r_index]
                              ->set_item(p,0,-receive_matrix->item(p,nb_row));
               Schur[dir].ei[comp][r_index]
                              ->set_item(0,p,-receive_matrix->item(nb_row,p));
            }
         }
         // Pre-thomas treatment on Schur complement
         GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur,r_index);

         // In case of periodic and multi-processor, there will be a variant
         // of Tridiagonal matrix instead of normal format
         // So, Schur complement of Schur complement is calculated
         if (is_periodic[field][dir] == 1) {
            Schur[dir].ee[comp][r_index]
                  ->set_item(0,0,A[dir].ee[comp][r_index]->item(nb_row,nb_row)
                                 -receive_matrix->item(nb_row,nb_row));

            ProdMatrix* SchurP = GLOBAL_EQ->get_SchurP(field);
            GLOBAL_EQ->compute_product_matrix_interior(Schur
                                                      ,SchurP
                                                      ,comp
                                                      ,0
                                                      ,dir
                                                      ,r_index);

            TDMatrix* DoubleSchur = GLOBAL_EQ-> get_DoubleSchur(field);
            nb_row = DoubleSchur[dir].ii_main[comp][r_index]->nb_rows();
            DoubleSchur[dir].ii_main[comp][r_index]
                  ->set_item(0,Schur[dir].ee[comp][r_index]->item(0,0)
                              -SchurP[dir].ei_ii_ie[comp]->item(0,0));
         }
      }
   } else if (is_periodic[field][dir] == 1) {
      // Condition for single processor in any
      // direction with periodic boundary conditions
      ProdMatrix* Ap = GLOBAL_EQ->get_Ap(field);
      GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir, field,r_index);

      LA_SeqMatrix* product_matrix = Ap[dir].ei_ii_ie[comp];
      LA_SeqMatrix* receive_matrix =
                    product_matrix->create_copy(this,product_matrix);

      A[dir].ee[comp][r_index]->set_item(0,0,Aee_diagcoef);

      TDMatrix* Schur = GLOBAL_EQ-> get_Schur(field);
      size_t nb_row = Schur[dir].ii_main[comp][r_index]->nb_rows();
      for (int p = 0; p < (int)nb_row; p++) {
         Schur[dir].ii_main[comp][r_index]
                  ->set_item(p,A[dir].ee[comp][r_index]->item(p,p)
                              -receive_matrix->item(p,p));
         if (p < (int)nb_row-1)
            Schur[dir].ii_super[comp][r_index]
                  ->set_item(p,-receive_matrix->item(p,p+1));
         if (p > 0)
            Schur[dir].ii_sub[comp][r_index]
                  ->set_item(p-1,-receive_matrix->item(p,p-1));
      }
      GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur,r_index);
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_1D_matrices ( FV_DiscreteField const* FF
                                        , FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: assemble_1D_matrices" ) ;

   double gamma = mu/2.0;

   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   // Assemble the matrices for pressure field(0) and velocity(1) field
   size_t field = (FF == PF) ? 0 : 1;

   for (size_t comp=0;comp<nb_comps[field];comp++) {
      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) =
                        FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) =
                        FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (size_t dir=0;dir<dim;dir++) {
         size_t dir_j, dir_k;

         if (dir == 0) {
            dir_j = 1; dir_k = 2;
         } else if (dir == 1) {
            dir_j = 0; dir_k = 2;
         } else if (dir == 2) {
            dir_j = 0; dir_k = 1;
         }

         size_t local_min_k = (dim == 2) ? 0 : min_unknown_index(dir_k);
         size_t local_max_k = (dim == 2) ? 0 : max_unknown_index(dir_k);

         LA_SeqMatrix* row_index = GLOBAL_EQ->get_row_indexes(field,dir,comp);

         for (size_t j = min_unknown_index(dir_j);
                    j <= max_unknown_index(dir_j);++j) {
            for (size_t k = local_min_k; k <= local_max_k; ++k) {
               size_t r_index = row_index->item(j,k);
               double Aee_diagcoef =
                     assemble_field_matrix (FF,t_it,gamma,comp,dir,j,k,r_index);

               assemble_field_schur_matrix(comp,dir,Aee_diagcoef,field,r_index);
            }
         }
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: NS_first_step ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: NS_first_step" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);
  // First Equation

  // Get local min and max indices
  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
         for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
            // Set P*_n as sum of P_(n-1/2)+phi_(n-1/2)
            double value = PF->DOF_value( i, j, k, 0, 0 )
                         + PF->DOF_value( i, j, k, 0, 1 );
            GLOBAL_EQ->update_global_P_vector(i,j,k,value);
         }
     }
  }

  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec_P();
  // Tranfer back to field
  PF->update_free_DOFs_value( 1, GLOBAL_EQ->get_solution_DS_pressure() ) ;

  PF->set_neumann_DOF_values();

  // Calculate pressure forces on the solid particles
  // if (is_stressCal) {
  //    compute_fluid_pressure_particle_interaction(t_it);
  // }

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_velocity_diffusion_terms ( )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: assemble_velocity_diffusion_terms" ) ;

   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);

   vector<doubleVector*> vel_diffusion = GLOBAL_EQ->get_velocity_diffusion();

   for (size_t comp=0;comp<nb_comps[1];comp++) {
      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) =
                        UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) =
                        UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
         for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
            for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
               size_t p = UF->DOF_local_number(i,j,k,comp);
               // dxx of level3
               vel_diffusion[0]->operator()(p) =
                                    compute_un_component(comp,i,j,k,0,3);
               // dyy of level1(2D) or level4(3D)
               if (dim == 2) {
                  vel_diffusion[1]->operator()(p) =
                                    compute_un_component(comp,i,j,k,1,1);
               } else {
                  vel_diffusion[1]->operator()(p) =
                                    compute_un_component(comp,i,j,k,1,4);
               }
               // dzz of level3
               if (dim == 3)
                  vel_diffusion[2]->operator()(p) =
                                    compute_un_component(comp,i,j,k,2,1);

            }
         }
      }
   }
}




//---------------------------------------------------------------------------
double
DDS_NavierStokes:: compute_un_component ( size_t const& comp,
                                          size_t const& i,
                                          size_t const& j,
                                          size_t const& k,
                                          size_t const& dir,
                                          size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: compute_un_component" ) ;

   double xhr=1.,xhl=1.,xright=0.,xleft=0.,yhr=1.,yhl=1.,yright=0.,yleft=0.;
   double zhr=1.,zhl=1.,zright=0.,zleft=0., value=0.;

   size_t_array2D* intersect_vector =
                           allrigidbodies->get_intersect_vector_on_grid(UF);
   doubleArray2D* intersect_distance =
                           allrigidbodies->get_intersect_distance_on_grid(UF);
   doubleArray2D* intersect_fieldVal =
                           allrigidbodies->get_intersect_fieldValue_on_grid(UF);

   size_t p = UF->DOF_local_number(i,j,k,comp);

   switch (dir) {
      case 0:
         if (UF->DOF_in_domain( (int)i+1, (int)j, (int)k, comp)) {
            xright = UF->DOF_value( i+1, j, k, comp, level )
                   - UF->DOF_value( i, j, k, comp, level ) ;
            xhr = UF->get_DOF_coordinate( i+1,comp, 0 )
                - UF->get_DOF_coordinate( i, comp, 0 ) ;
         }

         if (UF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp)) {
            xleft = UF->DOF_value( i, j, k, comp, level )
                  - UF->DOF_value( i-1, j, k, comp, level ) ;
            xhl = UF->get_DOF_coordinate( i, comp, 0 )
                - UF->get_DOF_coordinate( i-1, comp, 0 ) ;
         }

         if (is_solids) {
            size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
            if (void_frac->operator()(p) == 0) {
               if (intersect_vector->operator()(p,2*dir+0) == 1) {
                  xleft = UF->DOF_value( i, j, k, comp, level )
                        - intersect_fieldVal->operator()(p,2*dir+0);
                  xhl = intersect_distance->operator()(p,2*dir+0);
               }
               if (intersect_vector->operator()(p,2*dir+1) == 1) {
                  xright = intersect_fieldVal->operator()(p,2*dir+1)
                         - UF->DOF_value( i, j, k, comp, level );
                  xhr = intersect_distance->operator()(p,2*dir+1);
               }
            } else {
               xright = 0.; xleft = 0.;
            }
         }

         //xvalue = xright/xhr - xleft/xhl;
         if (UF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp)
            && UF->DOF_in_domain( (int)i+1, (int)j, (int)k, comp))
            value = xright/xhr - xleft/xhl;
         else if (UF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp))
            value = - xleft/xhl;
         else
            value = xright/xhr;
         break;
      case 1:
         if (UF->DOF_in_domain((int)i, (int)j+1, (int)k, comp)) {
            yright = UF->DOF_value( i, j+1, k, comp, level )
                   - UF->DOF_value( i, j, k, comp, level ) ;
            yhr = UF->get_DOF_coordinate( j+1,comp, 1 )
                - UF->get_DOF_coordinate( j, comp, 1 ) ;
         }

         if (UF->DOF_in_domain((int)i, (int)j-1, (int)k, comp)) {
            yleft = UF->DOF_value( i, j, k, comp, level )
                  - UF->DOF_value( i, j-1, k, comp, level ) ;
            yhl = UF->get_DOF_coordinate( j, comp, 1 )
                - UF->get_DOF_coordinate( j-1, comp, 1 ) ;
         }

         if (is_solids) {
            size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
            if (void_frac->operator()(p) == 0) {
               if (intersect_vector->operator()(p,2*dir+0) == 1) {
                  yleft = UF->DOF_value( i, j, k, comp, level )
                        - intersect_fieldVal->operator()(p,2*dir+0);
                  yhl = intersect_distance->operator()(p,2*dir+0);
               }
               if (intersect_vector->operator()(p,2*dir+1) == 1) {
                  yright = intersect_fieldVal->operator()(p,2*dir+1)
                         - UF->DOF_value( i, j, k, comp, level );
                  yhr = intersect_distance->operator()(p,2*dir+1);
               }
            } else {
               yleft = 0.; yright = 0.;
            }
         }

         //yvalue = yright/yhr - yleft/yhl;
         if (UF->DOF_in_domain((int)i, (int)j-1, (int)k, comp)
            && UF->DOF_in_domain((int)i, (int)j+1, (int)k, comp))
            value = yright/yhr - yleft/yhl;
         else if(UF->DOF_in_domain((int)i, (int)j-1, (int)k, comp))
            value = - yleft/yhl;
         else
            value = yright/yhr;
         break;
      case 2:
         if (UF->DOF_in_domain((int)i, (int)j, (int)k+1, comp)) {
            zright = UF->DOF_value( i, j, k+1, comp, level )
                   - UF->DOF_value( i, j, k, comp, level ) ;
            zhr = UF->get_DOF_coordinate( k+1,comp, 2 )
                - UF->get_DOF_coordinate( k, comp, 2 ) ;
         }

         if (UF->DOF_in_domain((int)i, (int)j, (int)k-1, comp)) {
            zleft = UF->DOF_value( i, j, k, comp, level )
                  - UF->DOF_value( i, j, k-1, comp, level ) ;
            zhl = UF->get_DOF_coordinate( k, comp, 2 )
                - UF->get_DOF_coordinate( k-1, comp, 2 ) ;
         }

         if (is_solids) {
            size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
            if (void_frac->operator()(p) == 0) {
               if (intersect_vector->operator()(p,2*dir+0) == 1) {
                  zleft = UF->DOF_value( i, j, k, comp, level )
                        - intersect_fieldVal->operator()(p,2*dir+0);
                  zhl = intersect_distance->operator()(p,2*dir+0);
               }
               if (intersect_vector->operator()(p,2*dir+1) == 1) {
                  zright = intersect_fieldVal->operator()(p,2*dir+1)
                         - UF->DOF_value( i, j, k, comp, level );
                  zhr = intersect_distance->operator()(p,2*dir+1);
               }
            } else {
               zleft = 0.; zright = 0.;
            }
         }

         //zvalue = zright/zhr - zleft/zhl;
         if (UF->DOF_in_domain((int)i, (int)j, (int)k-1, comp)
            && UF->DOF_in_domain((int)i, (int)j, (int)k+1, comp))
            value = zright/zhr - zleft/zhl;
         else if(UF->DOF_in_domain((int)i, (int)j, (int)k-1, comp))
            value = - zleft/zhl;
         else
            value = zright/zhr;
         break;
   }

   return(value);

}




//---------------------------------------------------------------------------
double
DDS_NavierStokes:: velocity_local_rhs ( size_t const& j
                                      , size_t const& k
                                      , double const& gamma
                                      , FV_TimeIterator const* t_it
                                      , size_t const& comp
                                      , size_t const& dir)
//---------------------------------------------------------------------------
{

   MAC_LABEL("DDS_NavierStokes:: velocity_local_rhs" ) ;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = UF->get_min_index_unknown_handled_by_proc(comp,l) ;
     max_unknown_index(l) = UF->get_max_index_unknown_handled_by_proc(comp,l) ;
   }

   // Compute VEC_rhs_x = rhs in x
   double fe=0.;

   // Since, this function is used in all directions;
   // ii, jj, and kk are used to convert the passed
   // arguments corresponding to correct direction
   size_t ii=0,jj=0,kk=0;
   size_t level=0;

   // Vector for fi
   LocalVector* VEC = GLOBAL_EQ->get_VEC(1);

   vector<doubleVector*> vel_diffusion = GLOBAL_EQ->get_velocity_diffusion();
   size_t_array2D* intersect_vector =
                           allrigidbodies->get_intersect_vector_on_grid(UF);
   doubleArray2D* intersect_distance =
                           allrigidbodies->get_intersect_distance_on_grid(UF);
   doubleArray2D* intersect_fieldVal =
                           allrigidbodies->get_intersect_fieldValue_on_grid(UF);

   for (size_t i=min_unknown_index(dir);i<=max_unknown_index(dir);++i) {
      if (dir == 0) {
         ii = i; jj = j; kk = k; level = 0;
      } else if (dir == 1) {
         ii = j; jj = i; kk = k; level = 3;
      } else if (dir == 2) {
         ii = j; jj = k; kk = i; level = 4;
      }

      size_t pos = i - min_unknown_index(dir);
      size_t p = UF->DOF_local_number(ii,jj,kk,comp);

      double value= vel_diffusion[dir]->operator()(p);

      if (is_solids) {
         if (intersect_vector->operator()(p,2*dir+0) == 1) {
            value = value - intersect_fieldVal->operator()(p,2*dir+0)
                           /intersect_distance->operator()(p,2*dir+0);
         }
         if (intersect_vector->operator()(p,2*dir+1) == 1) {
            value = value - intersect_fieldVal->operator()(p,2*dir+1)
                           /intersect_distance->operator()(p,2*dir+1);
         }
      }

      double dC = UF->get_cell_size(i,comp,dir);

      double temp_val = UF->DOF_value(ii,jj,kk,comp,level)
                        *dC*rho/t_it->time_step() - gamma*value;

      if (is_periodic[1][dir] == 0) {
        if (rank_in_i[dir] == nb_ranks_comm_i[dir]-1) {
           VEC[dir].local_T[comp]->set_item( pos,temp_val);
        } else {
           if (i == max_unknown_index(dir))
              fe = temp_val;
           else
              VEC[dir].local_T[comp]->set_item( pos,temp_val);
        }
      } else {
           if (i == max_unknown_index(dir))
              fe = temp_val;
           else
              VEC[dir].local_T[comp]->set_item( pos,temp_val);
      }

   }

   // Effect of boundary conditions in case of non-periodic direction
   int m = int(min_unknown_index(dir)) - 1;

   if (dir == 0) {
      ii = m; jj = j; kk = k;
   } else if (dir == 1) {
      ii = j; jj = m; kk = k;
   } else if (dir == 2) {
      ii = j; jj = k; kk = m;
   }

   if ( UF->DOF_in_domain((int)ii,(int)jj,(int)kk,comp))
      if ( UF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
         double ai = 1. / (UF->get_DOF_coordinate(m+1,comp,dir)
                         - UF->get_DOF_coordinate(m,comp,dir));
         double dirichlet_value = UF->DOF_value(ii,jj,kk,comp,1) ;
         VEC[dir].local_T[comp]->add_to_item( 0, + gamma*ai*dirichlet_value );
      }

   m = int(max_unknown_index(dir)) + 1;

   if (dir == 0) {
      ii = m; jj = j; kk = k;
   } else if (dir == 1) {
      ii = j; jj = m; kk = k;
   } else if (dir == 2) {
      ii = j; jj = k; kk = m;
   }

   if ( UF->DOF_in_domain((int)ii,(int)jj,(int)kk,comp))
      if ( UF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
         double ai = 1. / (UF->get_DOF_coordinate(m,comp,dir)
                         - UF->get_DOF_coordinate(m-1,comp,dir));
         double dirichlet_value = UF->DOF_value(ii,jj,kk,comp,1) ;
        VEC[dir].local_T[comp]->add_to_item(VEC[dir].local_T[comp]->nb_rows()-1,
                                                 + gamma*ai*dirichlet_value );
      }

   return fe;
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: unpack_compute_ue_pack(size_t const& comp
                                        , size_t const& dir
                                        , size_t const& p
                                        , size_t const& field)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: unpack_compute_ue_pack" ) ;

   LocalVector* VEC = GLOBAL_EQ->get_VEC(field);

   size_t nb_interface_unknowns = VEC[dir].T[comp]->nb_rows();

   for (size_t i=0;i<nb_interface_unknowns;i++) {
      VEC[dir].T[comp]->set_item(i,0);
      VEC[dir].interface_T[comp]->set_item(i,0);
   }

   if (is_periodic[field][dir])
      VEC[dir].T[comp]->set_item(nb_ranks_comm_i[dir]-1,
                     first_pass[field][dir].send[comp][rank_in_i[dir]][3*p]);

   VEC[dir].T[comp]->set_item(0,
                     first_pass[field][dir].send[comp][rank_in_i[dir]][3*p+1]);
   VEC[dir].interface_T[comp]->set_item(0,
                     first_pass[field][dir].send[comp][rank_in_i[dir]][3*p+2]);


   // Vec_temp might contain previous values
   for (size_t i=1;i<(size_t)nb_ranks_comm_i[dir];i++) {
      if (i!=(size_t)nb_ranks_comm_i[dir]-1) {
         VEC[dir].T[comp]->add_to_item(i-1,
                           first_pass[field][dir].receive[comp][i][3*p]);
         VEC[dir].T[comp]->add_to_item(i,
                           first_pass[field][dir].receive[comp][i][3*p+1]);
         // Assemble the interface rhs fe
         VEC[dir].interface_T[comp]->set_item(i,
                           first_pass[field][dir].receive[comp][i][3*p+2]);
      } else {
         if (is_periodic[field][dir] ==0) {
            VEC[dir].T[comp]->add_to_item(i-1,
                           first_pass[field][dir].receive[comp][i][3*p]);
         } else{
            VEC[dir].T[comp]->add_to_item(i-1,
                           first_pass[field][dir].receive[comp][i][3*p]);
            // If periodic in x, last proc has an interface unknown
            VEC[dir].T[comp]->add_to_item(i,
                           first_pass[field][dir].receive[comp][i][3*p+1]);
            VEC[dir].interface_T[comp]->set_item(i,
                           first_pass[field][dir].receive[comp][i][3*p+2]);
         }
      }
   }

   // Get fe - Aei*xi to solve for ue
   for (size_t i=0;i<nb_interface_unknowns;i++) {
      VEC[dir].interface_T[comp]->set_item(i,
                           VEC[dir].interface_T[comp]->item(i)
                           -VEC[dir].T[comp]->item(i));
   }

   // Solve for ue (interface unknowns) in the master proc
   DS_interface_unknown_solver(VEC[dir].interface_T[comp],comp,dir,field,p);

   for (size_t i=1;i<(size_t)nb_ranks_comm_i[dir];++i) {
      if (i != (size_t)nb_ranks_comm_i[dir]-1) {
         second_pass[field][dir].send[comp][i][2*p+0] =
                                          VEC[dir].interface_T[comp]->item(i-1);
         second_pass[field][dir].send[comp][i][2*p+1] =
                                          VEC[dir].interface_T[comp]->item(i);
      } else {
         second_pass[field][dir].send[comp][i][2*p+0] =
                                          VEC[dir].interface_T[comp]->item(i-1);
         if (is_periodic[field][dir])
            second_pass[field][dir].send[comp][i][2*p+1] =
                                          VEC[dir].interface_T[comp]->item(i);
         else
            second_pass[field][dir].send[comp][i][2*p+1] = 0;
      }
   }
}

//----------------------------------------------------------------------
void
DDS_NavierStokes::DS_interface_unknown_solver( LA_SeqVector* interface_rhs
                                             , size_t const& comp
                                             , size_t const& dir
                                             , size_t const& field
                                             , size_t const& r_index )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_interface_unknown_solver" ) ;

   TDMatrix* Schur = GLOBAL_EQ->get_Schur(field);

   // Condition for variant of Tridiagonal Schur
   // complement in Perioidic direction with multi-processor
   if ((is_periodic[field][dir] == 1) && (nb_ranks_comm_i[dir] != 1)) {
      LocalVector* Schur_VEC = GLOBAL_EQ->get_Schur_VEC(field);
      TDMatrix* DoubleSchur = GLOBAL_EQ->get_DoubleSchur(field);

      // Transfer interface_rhs to Schur VEC (i.e. S_fi and S_fe)
      size_t nrows = Schur_VEC[dir].local_T[comp]->nb_rows();
      for (size_t i = 0; i < nrows; i++) {
          Schur_VEC[dir].local_T[comp]->set_item(i,interface_rhs->item(i));
      }
      Schur_VEC[dir].interface_T[comp]->set_item(0,interface_rhs->item(nrows));

      // Calculate Sei*(Sii)-1*S_fi
      compute_Aei_ui(Schur,Schur_VEC,comp,dir,r_index);

      // Calculate S_fe - Sei*(Sii)-1*S_fi
      Schur_VEC[dir].interface_T[comp]->set_item(0,
                                 Schur_VEC[dir].interface_T[comp]->item(0)
                                -Schur_VEC[dir].T[comp]->item(0));

      // Calculate S_ue, using Schur complement of Schur complement
      GLOBAL_EQ->mod_thomas_algorithm(DoubleSchur,
                                      Schur_VEC[dir].interface_T[comp],
                                      comp,
                                      dir,
                                      r_index);

      // Calculate S_fi-Sie*S_ue
      Schur[dir].ie[comp][r_index]
               ->multiply_vec_then_add(Schur_VEC[dir].interface_T[comp]
                                      ,Schur_VEC[dir].local_T[comp],-1.0,1.0);

      // Calculate S_ui
      GLOBAL_EQ->mod_thomas_algorithm(Schur,
                                      Schur_VEC[dir].local_T[comp],
                                      comp,
                                      dir,
                                      r_index);

      // Transfer back the solution to interface_rhs
      for (size_t i = 0; i < nrows; i++) {
          interface_rhs->set_item(i,Schur_VEC[dir].local_T[comp]->item(i));
      }
      interface_rhs->set_item(nrows,Schur_VEC[dir].interface_T[comp]->item(0));
   } else {
      GLOBAL_EQ->mod_thomas_algorithm(Schur, interface_rhs, comp, dir,r_index);
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: unpack_ue(size_t const& comp
                           , double * received_data
                           , size_t const& dir
                           , size_t const& p
                           , size_t const& field)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: unpack_ue" ) ;

   LocalVector* VEC = GLOBAL_EQ->get_VEC(field);

   if (rank_in_i[dir] != nb_ranks_comm_i[dir]-1) {
      VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,
                                                   received_data[2*p]);
      VEC[dir].interface_T[comp]->set_item(rank_in_i[dir],
                                                   received_data[2*p+1]);
   } else {
      if (is_periodic[field][dir] ==0) {
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,
                                                   received_data[2*p]);
      } else {
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,
                                                   received_data[2*p]);
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir],
                                                   received_data[2*p+1]);
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: solve_interface_unknowns ( FV_DiscreteField* FF
                                            , double const& gamma
                                            , FV_TimeIterator const* t_it
                                            , size_t const& comp
                                            , size_t const& dir)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: solve_interface_unknowns" ) ;

   size_t field = (FF == PF) ? 0 : 1 ;

   // Get local min and max indices
   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);
   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) =
                        FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      max_unknown_index(l) =
                        FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   TDMatrix* A = GLOBAL_EQ-> get_A(field);
   LocalVector* VEC = GLOBAL_EQ->get_VEC(field);

   // Array declaration for sending data from master to all slaves
   size_t local_min_j=0, local_max_j=0;
   size_t local_min_k=0, local_max_k=0;

   if (dir == 0) {
      local_min_j = min_unknown_index(1);
      local_max_j = max_unknown_index(1);
      if (dim == 3) {
         local_min_k = min_unknown_index(2);
         local_max_k = max_unknown_index(2);
      }
   } else if (dir == 1) {
      local_min_j = min_unknown_index(0);
      local_max_j = max_unknown_index(0);
      if (dim == 3) {
         local_min_k = min_unknown_index(2);
         local_max_k = max_unknown_index(2);
      }
   } else if (dir == 2) {
      local_min_j = min_unknown_index(0);
      local_max_j = max_unknown_index(0);
      local_min_k = min_unknown_index(1);
      local_max_k = max_unknown_index(1);
   }

   LA_SeqMatrix* row_index = GLOBAL_EQ->get_row_indexes(field,dir,comp);

   // Send and receive the data first pass
   if ( rank_in_i[dir] == 0 ) {
      if (nb_ranks_comm_i[dir] != 1) {
         for (size_t i = 1; i < (size_t)nb_ranks_comm_i[dir]; ++i) {
            // Receive the data
            static MPI_Status status;
            MPI_Recv( first_pass[field][dir].receive[comp][i],
                (int) first_pass[field][dir].size[comp],
                MPI_DOUBLE, (int) i, 0, DDS_Comm_i[dir], &status ) ;
         }
      }

      for (size_t j = local_min_j; j <= local_max_j; j++) {
         for (size_t k = local_min_k; k <= local_max_k; k++) {

            size_t p = row_index->item(j,k);

            unpack_compute_ue_pack(comp,dir,p,field);

  	         // Need to have the original rhs function
            // assembled for corrosponding j,k pair
            assemble_local_rhs(j,k,gamma,t_it,comp,dir,field);

            // Setup RHS = fi - Aie*xe for solving ui
            A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp]
                                             ,VEC[dir].local_T[comp],-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_NavierStokes_solver(FF,j,k,min_unknown_index(dir)
                                                            ,comp,dir,p);
         }
      }

   } else {
      // Send the packed data to master
      MPI_Send( first_pass[field][dir].send[comp][rank_in_i[dir]],
          (int) first_pass[field][dir].size[comp],
          MPI_DOUBLE, 0, 0, DDS_Comm_i[dir] ) ;
   }

   // Send the data from master iff multi processor are used
   if (nb_ranks_comm_i[dir] != 1) {
      if ( rank_in_i[dir] == 0 ) {
         for (size_t i = 1; i < (size_t)nb_ranks_comm_i[dir]; ++i) {
            MPI_Send( second_pass[field][dir].send[comp][i],
                (int) second_pass[field][dir].size[comp],
                MPI_DOUBLE,(int) i, 0, DDS_Comm_i[dir] ) ;
         }
      } else {
         // Receive the data
         static MPI_Status status ;
         MPI_Recv( second_pass[field][dir].send[comp][rank_in_i[dir]],
             (int) first_pass[field][dir].size[comp],
             MPI_DOUBLE, 0, 0, DDS_Comm_i[dir], &status ) ;

         // Solve the system of equations in each proc
         for (size_t j = local_min_j; j <= local_max_j; j++) {
            for (size_t k = local_min_k; k <= local_max_k; k++) {
               size_t p = row_index->item(j,k);

               unpack_ue(comp,second_pass[field][dir].send[comp][rank_in_i[dir]]
                             ,dir,p,field);

               // Need to have the original rhs function
               // assembled for corrosponding j,k pair
               assemble_local_rhs(j,k,gamma,t_it,comp,dir,field);

               // Setup RHS = fi - Aie*xe for solving ui
               A[dir].ie[comp][p]
                       ->multiply_vec_then_add(VEC[dir].interface_T[comp]
                                              ,VEC[dir].local_T[comp],-1.0,1.0);

               // Solve ui and transfer solution into distributed vector
               GLOBAL_EQ->DS_NavierStokes_solver(FF,j,k,min_unknown_index(dir)
                                                            ,comp,dir,p);
            }
         }
      }
   }
}
//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: second_order_pressure_stress(size_t const& parID, size_t const& Np, FV_TimeIterator const* t_it )
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: second_order_pressure_stress" ) ;
//
//   size_t i0_temp;
//   size_t comp = 0;
//
//   // Structure of particle input data
//   PartInput solid = GLOBAL_EQ->get_solid(0);
//   // Structure of particle surface input data
//   SurfaceDiscretize surface = GLOBAL_EQ->get_surface();
//
//   // comp won't matter as the particle position is independent of comp
//   double xp = solid.coord[0]->item(parID);
//   double yp = solid.coord[1]->item(parID);
//   double zp = solid.coord[2]->item(parID);
//   double ri = solid.size->item(parID);
//
//   doubleArray2D ipoint(3,3,0.);
//   doubleVector fini(3,0);
//   doubleVector stress(Np,0);
//   boolVector in_domain(2,true);        //true if ghost point in the computational domain
//   size_t_vector in_parID(2,0);         //Store particle ID if level_set becomes negative
//   boolArray2D found(3,dim,false);
//   size_t_array2D i0(3,dim,0);
//
//   size_t_vector min_unknown_index(dim,0);
//   size_t_vector max_unknown_index(dim,0);
//
//   doubleVector Dmin(dim,0);
//   doubleVector Dmax(dim,0);
//   doubleVector rotated_coord(3,0);
//   doubleVector rotated_normal(3,0);
//   doubleArray2D surface_force(Np,3,0);
//   doubleArray2D surface_point(Np,3,0);
//
//   // Structure of particle input data
//   PartForces hydro_forces = GLOBAL_EQ->get_forces(0);
//   PartForces hydro_torque = GLOBAL_EQ->get_torque(0);
//
//   // Force vector
//   doubleVector force(3,0);
//   doubleVector torque(3,0);
//
//   for (size_t i=0;i<Np;i++) {
//      // Get local min and max indices
//      for (size_t l=0;l<dim;++l) {
//         min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( comp, l );
//         max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( comp, l );
//         Dmin(l) = PF->primary_grid()->get_min_coordinate_on_current_processor(l);
//         Dmax(l) = PF->primary_grid()->get_max_coordinate_on_current_processor(l);
//      }
//
//      // Rotating surface points
//      rotated_coord(0) = ri*surface.coordinate[0]->item(i);
//      rotated_coord(1) = ri*surface.coordinate[1]->item(i);
//      rotated_coord(2) = ri*surface.coordinate[2]->item(i);
//
//      rotation_matrix(parID,rotated_coord,comp,0);
//
//      ipoint(0,0) = xp + rotated_coord(0);
//      ipoint(0,1) = yp + rotated_coord(1);
//      ipoint(0,2) = zp + rotated_coord(2);
//
//      // Rotating surface normal
//      rotated_normal(0) = surface.normal[0]->item(i);
//      rotated_normal(1) = surface.normal[1]->item(i);
//      rotated_normal(2) = surface.normal[2]->item(i);
//
//      rotation_matrix(parID,rotated_normal,comp,0);
//
//      // Correction in case of periodic boundary condition in any direction
//      for (size_t dir=0;dir<dim;dir++) {
//         // PBC on rotated surface points
//         if (is_periodic[0][dir]) {
//            double isize = PF->primary_grid()->get_main_domain_max_coordinate(dir) - PF->primary_grid()->get_main_domain_min_coordinate(dir);
//            double imin = PF->primary_grid()->get_main_domain_min_coordinate(dir);
//            ipoint(0,dir) = ipoint(0,dir) - MAC::floor((ipoint(0,dir)-imin)/isize)*isize;
//         }
//         // Finding the grid indexes next to ghost points
//         found(0,dir) = FV_Mesh::between(PF->get_DOF_coordinates_vector(comp,dir), ipoint(0,dir), i0_temp) ;
//         if (found(0,dir) == 1) i0(0,dir) = i0_temp ;
//      }
//
//      // Accessing the smallest grid size in domain
//      double dh = PF->primary_grid()->get_smallest_grid_size();
//
//      bool status = (dim==2) ? ((ipoint(0,0) > Dmin(0)) && (ipoint(0,0) <= Dmax(0)) && (ipoint(0,1) > Dmin(1)) && (ipoint(0,1) <= Dmax(1))) :
//                               ((ipoint(0,0) > Dmin(0)) && (ipoint(0,0) <= Dmax(0)) && (ipoint(0,1) > Dmin(1)) && (ipoint(0,1) <= Dmax(1))
//                                                                                    && (ipoint(0,2) > Dmin(2)) && (ipoint(0,2) <= Dmax(2)));
//
//      if (status) {
//         for (size_t l=0;l<dim;++l) {
//            // Ghost points in normal direction at the particle surface
//            ipoint(1,l) = ipoint(0,l) + dh*rotated_normal(l);
//            ipoint(2,l) = ipoint(0,l) + 2.*dh*rotated_normal(l);
//
//            // Periocid boundary conditions
//            if (is_periodic[0][l]) {
//               double isize = PF->primary_grid()->get_main_domain_max_coordinate(l) - PF->primary_grid()->get_main_domain_min_coordinate(l);
//               double imin = PF->primary_grid()->get_main_domain_min_coordinate(l);
//               ipoint(1,l) = ipoint(1,l) - MAC::floor((ipoint(1,l)-imin)/isize)*isize;
//               ipoint(2,l) = ipoint(2,l) - MAC::floor((ipoint(2,l)-imin)/isize)*isize;
//            }
//
//            // Finding the grid indexes next to ghost points
//            found(1,l) = FV_Mesh::between(PF->get_DOF_coordinates_vector(comp,l), ipoint(1,l), i0_temp);
//            if (found(1,l) == 1) i0(1,l) = i0_temp;
//            found(2,l) = FV_Mesh::between(PF->get_DOF_coordinates_vector(comp,l), ipoint(2,l), i0_temp);
//            if (found(2,l) == 1) i0(2,l) = i0_temp;
//         }
//
// 	// Check the points in doamin or not
//         in_domain(0) = (dim == 2) ? found(1,0) && found(1,1) :
//                                     found(1,0) && found(1,1) && found(1,2) ;
//         in_domain(1) = (dim == 2) ? found(2,0) && found(2,1) :
//                                     found(2,0) && found(2,1) && found(2,2) ;
//
// 	// Calculation of field variable on ghost point(0,0)
//         for (size_t level=0; level<2;level++) {
//            // Calculation of field variable on ghost point(1)
//            if (in_domain(0))
//               fini(1) = (dim == 2) ? ghost_field_estimate_on_face (PF,comp,i0(1,0),i0(1,1),0, ipoint(1,0), ipoint(1,1),0, 0.,2,level) :
//                         ghost_field_estimate_in_box (PF,comp,i0(1,0),i0(1,1),i0(1,2),ipoint(1,0),ipoint(1,1),ipoint(1,2),0.,level,parID) ;
//            // Calculation of field variable on ghost point(2)
//            if (in_domain(1))
//               fini(2) = (dim == 2) ? ghost_field_estimate_on_face (PF,comp,i0(2,0),i0(2,1),0, ipoint(2,0), ipoint(2,1),0, 0.,2,level) :
//                         ghost_field_estimate_in_box (PF,comp,i0(2,0),i0(2,1),i0(2,2),ipoint(2,0),ipoint(2,1),ipoint(2,2),0.,level,parID) ;
//
//            // Both points 1 and 2 in the computational domain
//            if (in_domain(0) && in_domain(1)) {
//               fini(0) = 2.*fini(1) - fini(2);
//            // Point 1 in the computational domain
//            } else if (in_domain(0)) {
//               fini(0) = fini(1);
//            // Point 2 in the computational domain
//            } else if (in_domain(1)) {
//               fini(0) = fini(2);
//            }
//            stress(i) = stress(i) - fini(0)/2.;
// 	}
//      }
//
//      double scale = (dim == 2) ? ri : ri*ri;
//
//      // Ref: Keating thesis Pg-85
//      // point_coord*(area) --> Component of area in particular direction
//      double fx = stress(i)*rotated_normal(0)*(surface.area->item(i)*scale);
//      double fy = stress(i)*rotated_normal(1)*(surface.area->item(i)*scale);
//      double fz = stress(i)*rotated_normal(2)*(surface.area->item(i)*scale);
//
//      force(0) = force(0) + fx ;
//      force(1) = force(1) + fy ;
//      force(2) = force(2) + fz ;
//
//      surface_point(i,0) = ipoint(0,0);
//      surface_point(i,1) = ipoint(0,1);
//      surface_point(i,2) = ipoint(0,2);
//
//      surface_force(i,0) = fx;
//      surface_force(i,1) = fy;
//      surface_force(i,2) = fz;
//
//      torque(0) = torque(0) + fz*rotated_coord(1) - fy*rotated_coord(2);
//      torque(1) = torque(1) + fx*rotated_coord(2) - fz*rotated_coord(0);
//      torque(2) = torque(2) + fy*rotated_coord(0) - fx*rotated_coord(1);
//   }
//
//   write_surface_discretized_forces(PF,Np,parID,surface_point,surface_force,t_it);
//
//   hydro_forces.press[0]->set_item(parID,force(0)) ;
//   hydro_forces.press[1]->set_item(parID,force(1)) ;
//   hydro_forces.press[2]->set_item(parID,force(2)) ;
//
//   hydro_torque.press[0]->set_item(parID,torque(0)) ;
//   hydro_torque.press[1]->set_item(parID,torque(1)) ;
//   hydro_torque.press[2]->set_item(parID,torque(2)) ;
// }
//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: second_order_pressure_stress_withNeumannBC(size_t const& parID, size_t const& Np, FV_TimeIterator const* t_it)
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: second_order_pressure_stress_withNeumannBC" ) ;
//
//   size_t i0_temp;
//   size_t comp = 0;
//
//   // Structure of particle input data
//   PartInput solid = GLOBAL_EQ->get_solid(0);
//   // Structure of particle surface input data
//   SurfaceDiscretize surface = GLOBAL_EQ->get_surface();
//
//   // comp won't matter as the particle position is independent of comp
//   double xp = solid.coord[0]->item(parID);
//   double yp = solid.coord[1]->item(parID);
//   double zp = solid.coord[2]->item(parID);
//   double ri = solid.size->item(parID);
//
//   doubleArray2D point(3,3,0);
//   doubleVector fini(3,0);
//   doubleVector level_set(2,1.);
//   boolArray2D point_in_domain(1,2,true);
//   size_t_vector in_parID(2,0);         //Store particle ID if level_set becomes negative
//   boolArray2D found(dim,3,false);
//   size_t_array2D i0(3,3,0);
//   doubleVector stress(Np,0);
//   doubleArray2D surface_force(Np,3,0);
//   doubleArray2D surface_point(Np,3,0);
//
//   size_t_vector min_unknown_index(dim,0);
//   size_t_vector max_unknown_index(dim,0);
//
//   doubleVector Dmin(dim,0);
//   doubleVector Dmax(dim,0);
//   doubleVector rotated_coord(3,0);
//   doubleVector rotated_normal(3,0);
//   intVector sign(dim,0);
//
//   // Structure of particle input data
//   PartForces hydro_forces = GLOBAL_EQ->get_forces(0);
//   PartForces hydro_torque = GLOBAL_EQ->get_torque(0);
//
//   // Force vector
//   doubleVector force(3,0);
//   doubleVector torque(3,0);
//
//   for (size_t i=0;i<Np;i++) {
//      // Get local min and max indices
//      // Get min and max coordinates in the current processor
//      for (size_t l=0;l<dim;++l) {
//          min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( comp, l );
//          max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( comp, l );
//          Dmin(l) = PF->primary_grid()->get_min_coordinate_on_current_processor(l);
//          Dmax(l) = PF->primary_grid()->get_max_coordinate_on_current_processor(l);
//      }
//
//      // Rotating surface points
//      rotated_coord(0) = ri*surface.coordinate[0]->item(i);
//      rotated_coord(1) = ri*surface.coordinate[1]->item(i);
//      rotated_coord(2) = ri*surface.coordinate[2]->item(i);
//
//      rotation_matrix(parID,rotated_coord,comp,0);
//
//      point(0,0) = xp + rotated_coord(0);
//      point(0,1) = yp + rotated_coord(1);
//      point(0,2) = zp + rotated_coord(2);
//
//      // Rotating surface normal
//      rotated_normal(0) = surface.normal[0]->item(i);
//      rotated_normal(1) = surface.normal[1]->item(i);
//      rotated_normal(2) = surface.normal[2]->item(i);
//
//      rotation_matrix(parID,rotated_normal,comp,0);
//
//      for (size_t dir=0;dir<dim;dir++) {
//         // PBC on rotated surface points
//         if (is_periodic[0][dir]) {
//            double isize = PF->primary_grid()->get_main_domain_max_coordinate(dir) - PF->primary_grid()->get_main_domain_min_coordinate(dir);
//            double imin = PF->primary_grid()->get_main_domain_min_coordinate(dir);
//            point(0,dir) = point(0,dir) - MAC::floor((point(0,dir)-imin)/isize)*isize;
//         }
//         // Finding the grid indexes next to ghost points
//         found(dir,0) = FV_Mesh::between(PF->get_DOF_coordinates_vector(comp,dir), point(0,dir), i0_temp);
//         if (found(dir,0) == 1) i0(0,dir) = i0_temp;
//      }
//
//      // Accessing the smallest grid size in domain
//      double dh = PF->primary_grid()->get_smallest_grid_size();
//
//      bool status = (dim==2) ? ((point(0,0) > Dmin(0)) && (point(0,0) <= Dmax(0)) && (point(0,1) > Dmin(1)) && (point(0,1) <= Dmax(1))) :
//                               ((point(0,0) > Dmin(0)) && (point(0,0) <= Dmax(0)) && (point(0,1) > Dmin(1)) && (point(0,1) <= Dmax(1))
//                                                                                  && (point(0,2) > Dmin(2)) && (point(0,2) <= Dmax(2)));
//      double threshold = pow(loc_thres,0.5)*dh;
//      size_t major_dir=4;
//
//      if (status) {
//         for (size_t dir=0;dir<dim;dir++) {
//            sign(dir) = (rotated_normal(dir) > 0.) ? 1 : -1;
// 	   if (dim == 2) {
//               if (MAC::abs(rotated_normal(dir)) == MAC::max(MAC::abs(rotated_normal(0)),MAC::abs(rotated_normal(1))))
//                  major_dir = dir;
// 	   } else if (dim == 3) {
// 	      if (MAC::abs(rotated_normal(dir)) == MAC::max(MAC::abs(rotated_normal(0)),
// 				                   MAC::max(MAC::abs(rotated_normal(1)),MAC::abs(rotated_normal(2)))))
// 		 major_dir = dir;
// 	   }
// 	}
//
// 	// Ghost points in i for the calculation of i-derivative of field
//         ghost_points_generation( PF, point, i0, sign(major_dir), comp, major_dir, point_in_domain, rotated_normal);
//
//         // Assuming all ghost points are in fluid
//         level_set(0) = 1.; level_set(1) = 1.;
//
//         // Checking all the ghost points in the solid/fluid, and storing the parID if present in solid
//         for (size_t m=0;m<Npart;m++) {
//            // In Normal direction
//            if (level_set(0) > threshold) {
//               level_set(0) = level_set_function(PF,m,comp,point(1,0),point(1,1),point(1,2),level_set_type,0);
//               level_set(0) *= solid.inside->item(m);
//               if (level_set(0) < threshold) in_parID(0) = m;
//            }
//            if (level_set(1) > threshold) {
//               level_set(1) = level_set_function(PF,m,comp,point(2,0),point(2,1),point(2,2),level_set_type,0);
//               level_set(1) *= solid.inside->item(m);
//               if (level_set(1) < threshold) in_parID(1) = m;
//            }
//         }
//
//         double dx1 = pow(pow(point(1,0)-point(0,0),2) + pow(point(1,1)-point(0,1),2) + pow(point(1,2)-point(0,2),2),0.5);
//         double dx2 = pow(pow(point(2,0)-point(0,0),2) + pow(point(2,1)-point(0,1),2) + pow(point(2,2)-point(0,2),2),0.5);
//
//         // Calculation of field variable on ghost point(0,0)
//         for (size_t level=0; level<2;level++) {
//            // Calculation of field variable on ghost point(1)
//            if ((level_set(0) > threshold) && point_in_domain(0,0)) fini(1) = third_order_ghost_field_estimate(PF, comp, point(1,0), point(1,1), point(1,2), i0(1,0), i0(1,1), i0(1,2), major_dir, sign, level);
//            // Calculation of field variable on ghost point(2)
//            if ((level_set(1) > threshold) && point_in_domain(0,1)) fini(2) = third_order_ghost_field_estimate(PF, comp, point(2,0), point(2,1), point(2,2), i0(2,0), i0(2,1), i0(2,2), major_dir, sign, level);
//            // Extrapolation of point on surface using point(1) and point(2), assuming Neumann BC on particle boundary
//            fini(0) = ((dx2/dx1)*fini(1) - (dx1/dx2)*fini(2))/(dx2/dx1 - dx1/dx2);
//            // Stress accumulation
//            stress(i) = stress(i) - fini(0)/2.;
//         }
//      }
//
//      double scale = (dim == 2) ? ri : ri*ri;
//
//      // Ref: Keating thesis Pg-85
//      // point_coord*(area) --> Component of area in particular direction
//      double fx = stress(i)*rotated_normal(0)*(surface.area->item(i)*scale);
//      double fy = stress(i)*rotated_normal(1)*(surface.area->item(i)*scale);
//      double fz = stress(i)*rotated_normal(2)*(surface.area->item(i)*scale);
//
//      force(0) = force(0) + fx ;
//      force(1) = force(1) + fy ;
//      force(2) = force(2) + fz ;
//
//      surface_point(i,0) = point(0,0);
//      surface_point(i,1) = point(0,1);
//      surface_point(i,2) = point(0,2);
//
//      surface_force(i,0) = fx;
//      surface_force(i,1) = fy;
//      surface_force(i,2) = fz;
//
//      torque(0) = torque(0) + fz*rotated_coord(1) - fy*rotated_coord(2);
//      torque(1) = torque(1) + fx*rotated_coord(2) - fz*rotated_coord(0);
//      torque(2) = torque(2) + fy*rotated_coord(0) - fx*rotated_coord(1);
//   }
//
//   write_surface_discretized_forces(PF,Np,parID,surface_point,surface_force,t_it);
//
//   hydro_forces.press[0]->set_item(parID,force(0)) ;
//   hydro_forces.press[1]->set_item(parID,force(1)) ;
//   hydro_forces.press[2]->set_item(parID,force(2)) ;
//
//   hydro_torque.press[0]->set_item(parID,torque(0)) ;
//   hydro_torque.press[1]->set_item(parID,torque(1)) ;
//   hydro_torque.press[2]->set_item(parID,torque(2)) ;
// }
//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: first_order_pressure_stress(size_t const& parID, size_t const& Np, FV_TimeIterator const* t_it)
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: first_order_pressure_stress" ) ;
//
//   size_t i0_temp;
//   double ri=0.;
//   bool found = 0;
//
//   doubleVector point(3,0);
//   doubleVector stress(Np,0);
//   size_t i0, j0, k0=0;
//
//   size_t_vector min_unknown_index(dim,0);
//   size_t_vector max_unknown_index(dim,0);
//
//   doubleVector Dmin(dim,0);
//   doubleVector Dmax(dim,0);
//   doubleVector rotated_coord(3,0);
//   doubleVector rotated_normal(3,0);
//   doubleArray2D surface_force(Np,3,0);
//   doubleArray2D surface_point(Np,3,0);
//
//   // Structure of particle input data
//   PartInput solid = GLOBAL_EQ->get_solid(0);
//   // Structure of particle surface data
//   SurfaceDiscretize surface = GLOBAL_EQ->get_surface();
//   // Structure of particle input data
//   PartForces hydro_forces = GLOBAL_EQ->get_forces(0);
//   PartForces hydro_torque = GLOBAL_EQ->get_torque(0);
//
//   // Force vector
//   doubleVector force(3,0);
//   doubleVector torque(3,0);
//
//   for (size_t i=0;i<Np;i++) {
//      double sx = surface.coordinate[0]->item(i);
//      double sy = surface.coordinate[1]->item(i);
//      double sz = surface.coordinate[2]->item(i);
//
//      double s_nx = surface.normal[0]->item(i);
//      double s_ny = surface.normal[1]->item(i);
//      double s_nz = surface.normal[2]->item(i);
//
//      double s_area = surface.area->item(i);
//
//      for (size_t comp=0;comp<nb_comps[0];comp++) {
//         // Get local min and max indices
//         // One extra grid cell needs to considered, since ghost points can be
//         // located in between the min/max index handled by the proc
//         for (size_t l=0;l<dim;++l) {
//            min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( comp, l );
//            max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( comp, l );
// 	   Dmin(l) = PF->primary_grid()->get_min_coordinate_on_current_processor(l);
// 	   Dmax(l) = PF->primary_grid()->get_max_coordinate_on_current_processor(l);
//         }
//
//         double xp = solid.coord[0]->item(parID);
//         double yp = solid.coord[1]->item(parID);
//         double zp = solid.coord[2]->item(parID);
//         ri = solid.size->item(parID);
//
// 	// Rotating surface points
// 	rotated_coord(0) = ri*sx;
// 	rotated_coord(1) = ri*sy;
// 	rotated_coord(2) = ri*sz;
//
//         rotation_matrix(parID,rotated_coord,comp,0);
//
//         point(0) = xp + rotated_coord(0);
//         point(1) = yp + rotated_coord(1);
//         point(2) = zp + rotated_coord(2);
//
// 	// Rotating surface vectors
// 	rotated_normal(0) = s_nx;
// 	rotated_normal(1) = s_ny;
// 	rotated_normal(2) = s_nz;
//         rotation_matrix(parID,rotated_normal,comp,0);
//
//       	// Displacement correction in case of periodic boundary condition in any or all directions
//         for (size_t dir=0;dir<dim;dir++) {
//            if (is_periodic[0][dir]) {
//               double isize = PF->primary_grid()->get_main_domain_max_coordinate(dir) - PF->primary_grid()->get_main_domain_min_coordinate(dir);
//               double imin = PF->primary_grid()->get_main_domain_min_coordinate(dir);
//               point(dir) = point(dir) - MAC::floor((point(dir)-imin)/isize)*isize;
//            }
//         }
//
//         bool status = (dim==2) ? ((point(0) > Dmin(0)) && (point(0) <= Dmax(0)) && (point(1) > Dmin(1)) && (point(1) <= Dmax(1))) :
//                                  ((point(0) > Dmin(0)) && (point(0) <= Dmax(0)) && (point(1) > Dmin(1)) && (point(1) <= Dmax(1))
//                                                                                 && (point(2) > Dmin(2)) && (point(2) <= Dmax(2)));
//
//         if (status) {
//            // Finding the grid indexes next to ghost points
//            found = FV_Mesh::between(PF->get_DOF_coordinates_vector(comp,0), point(0), i0_temp);
//            if (found == 1) i0 = i0_temp;
//
//            found = FV_Mesh::between(PF->get_DOF_coordinates_vector(comp,1), point(1), i0_temp);
//            if (found == 1) j0 = i0_temp;
//
//            if (dim == 3) {
//               found = FV_Mesh::between(PF->get_DOF_coordinates_vector(comp,2), point(2), i0_temp);
//               if (found == 1) k0 = i0_temp;
//            }
//
//            // Calculation of field variable on ghost point(0,0)
//            for (size_t level=0; level<2;level++) {
//               if (dim == 2) {
//                  double press0 = ghost_field_estimate_on_face (PF,comp,i0,j0,0,point(0),point(1),0,0.,2,level);
//                  stress(i) = stress(i) - press0/2.;
//               } else if (dim == 3) {
//                  double press0 = ghost_field_estimate_in_box (PF,comp,i0,j0,k0,point(0),point(1),point(2),0.,level,0);
//                  stress(i) = stress(i) - press0/2.;
//               }
//            }
//         }
//      }
//
//      double scale = (dim == 2) ? ri : ri*ri;
//
//      // Ref: Keating thesis Pg-85
//      // point_coord*(area) --> Component of area in particular direction
//      double fx = stress(i)*rotated_normal(0)*(s_area*scale);
//      double fy = stress(i)*rotated_normal(1)*(s_area*scale);
//      double fz = stress(i)*rotated_normal(2)*(s_area*scale);
//
//      surface_point(i,0) = point(0);
//      surface_point(i,1) = point(1);
//      surface_point(i,2) = point(2);
//
//      surface_force(i,0) = fx;
//      surface_force(i,1) = fy;
//      surface_force(i,2) = fz;
//
//      force(0) = force(0) + fx ;
//      force(1) = force(1) + fy ;
//      force(2) = force(2) + fz ;
//
//      torque(0) = torque(0) + fz*rotated_coord(1) - fy*rotated_coord(2);
//      torque(1) = torque(1) + fx*rotated_coord(2) - fz*rotated_coord(0);
//      torque(2) = torque(2) + fy*rotated_coord(0) - fx*rotated_coord(1);
//
//   }
//
//   write_surface_discretized_forces(PF,Np,parID,surface_point,surface_force,t_it);
//
//   hydro_forces.press[0]->set_item(parID,force(0)) ;
//   hydro_forces.press[1]->set_item(parID,force(1)) ;
//   hydro_forces.press[2]->set_item(parID,force(2)) ;
//
//   hydro_torque.press[0]->set_item(parID,torque(0)) ;
//   hydro_torque.press[1]->set_item(parID,torque(1)) ;
//   hydro_torque.press[2]->set_item(parID,torque(2)) ;
//
// }
//
//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: compute_fluid_pressure_particle_interaction( FV_TimeIterator const* t_it)
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: compute_fluid_pressure_particle_interaction" ) ;
//
//   // Structure of particle input data
//   PartForces hydro_forces = GLOBAL_EQ->get_forces(0);
//   PartForces hydro_torque = GLOBAL_EQ->get_torque(0);
//
//   size_t Nmax = 0.;
//   if (level_set_type == "Sphere") {
//      Nmax = (dim == 2) ? ((size_t) Npoints) : ((size_t) (2*Npoints)) ;
//   } else if (level_set_type == "Cube") {
//      Nmax = (dim == 2) ? ((size_t) (4*Npoints)) : ((size_t) (6*pow(Npoints,2))) ;
//   } else if (level_set_type == "Cylinder") {
//      double Npm1 = round(pow(MAC::sqrt(Npoints) - MAC::sqrt(MAC::pi()/ar),2.));
//      double Ncyl = (Npoints - Npm1);
//      double dh = 1. - MAC::sqrt(Npm1/Npoints);
//      double Nr = round(2./dh);
//      Nmax = (dim == 3) ? ((size_t)(2*Npoints + Nr*Ncyl)) : 0 ;
//   }
//
//   for (size_t parID = 0; parID < Npart; parID++) {
//
//      // Contribution due to pressure tensor
//      compute_pressure_force_on_particle(parID, Nmax, t_it);
//      // Gathering information from all procs
//      hydro_forces.press[0]->set_item(parID,
// 		     macCOMM->sum(hydro_forces.press[0]->item(parID))) ;
//      hydro_forces.press[1]->set_item(parID,
// 		     macCOMM->sum(hydro_forces.press[1]->item(parID))) ;
//      hydro_forces.press[2]->set_item(parID,
// 		     macCOMM->sum(hydro_forces.press[2]->item(parID))) ;
//      hydro_torque.press[0]->set_item(parID,
// 		     macCOMM->sum(hydro_torque.press[0]->item(parID))) ;
//      hydro_torque.press[1]->set_item(parID,
// 		     macCOMM->sum(hydro_torque.press[1]->item(parID))) ;
//      hydro_torque.press[2]->set_item(parID,
// 		     macCOMM->sum(hydro_torque.press[2]->item(parID))) ;
//
//   }
// }
//
//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: compute_fluid_velocity_particle_interaction( FV_TimeIterator const* t_it)
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: compute_fluid_velocity_particle_interaction" ) ;
//
//   string fileName = "./DS_results/particle_forces.csv" ;
//
//   ofstream MyFile( fileName.c_str(), ios::app ) ;
//   if ((t_it->time() == t_it->time_step()) && (my_rank == 0))
//      MyFile << "Time,parID,Npoints,x,y,z,vx,vy,vz,Fpx,Fpy,Fpz,Fvx,Fvy,Fvz,Tpx,Tpy,Tpz,Tvx,Tvy,Tvz" << endl;
//
//   doubleArray2D avg_force(3,2,0);
//   doubleArray2D avg_torque(3,2,0);
//
//   // Structure of hydrodynamic forces
//   PartForces hydro_forces = GLOBAL_EQ->get_forces(0);
//   PartForces hydro_torque = GLOBAL_EQ->get_torque(0);
//   // Structure of particle input data
//   PartInput solid = GLOBAL_EQ->get_solid(0);
//
//   size_t Nmax = 0.;
//   if (level_set_type == "Sphere") {
//      Nmax = (dim == 2) ? ((size_t) Npoints) : ((size_t) (2*Npoints)) ;
//   } else if (level_set_type == "Cube") {
//      Nmax = (dim == 2) ? ((size_t) (4*Npoints)) : ((size_t) (6*pow(Npoints,2))) ;
//   } else if (level_set_type == "Cylinder") {
//      double Npm1 = round(pow(MAC::sqrt(Npoints) - MAC::sqrt(MAC::pi()/ar),2.));
//      double Ncyl = (Npoints - Npm1);
//      double dh = 1. - MAC::sqrt(Npm1/Npoints);
//      double Nr = round(2./dh);
//      Nmax = (dim == 3) ? ((size_t)(2*Npoints + Nr*Ncyl)) : 0 ;
//   }
//
//   for (size_t parID = 0; parID < Npart; parID++) {
//      // comp won't matter as the particle position is independent of comp
//      double xp = solid.coord[0]->item(parID);
//      double yp = solid.coord[1]->item(parID);
//      double zp = solid.coord[2]->item(parID);
//      double vx = solid.vel[0]->item(parID);
//      double vy = solid.vel[1]->item(parID);
//      double vz = solid.vel[2]->item(parID);
//
//      // Contribution of stress tensor
//      compute_velocity_force_on_particle(parID, Nmax, t_it);
//      // Gathering information from all procs
//      hydro_forces.vel[0]->set_item(parID,
// 		     macCOMM->sum(hydro_forces.vel[0]->item(parID))) ;
//      hydro_forces.vel[1]->set_item(parID,
// 		     macCOMM->sum(hydro_forces.vel[1]->item(parID))) ;
//      hydro_forces.vel[2]->set_item(parID,
// 		     macCOMM->sum(hydro_forces.vel[2]->item(parID))) ;
//      hydro_torque.vel[0]->set_item(parID,
// 		     macCOMM->sum(hydro_torque.vel[0]->item(parID))) ;
//      hydro_torque.vel[1]->set_item(parID,
// 		     macCOMM->sum(hydro_torque.vel[1]->item(parID))) ;
//      hydro_torque.vel[2]->set_item(parID,
// 		     macCOMM->sum(hydro_torque.vel[2]->item(parID))) ;
//
//
//      if (my_rank == 0) {
// 	avg_force(0,0) += hydro_forces.press[0]->item(parID);
// 	avg_force(1,0) += hydro_forces.press[1]->item(parID);
// 	avg_force(2,0) += hydro_forces.press[2]->item(parID);
// 	avg_force(0,1) += hydro_forces.vel[0]->item(parID);
// 	avg_force(1,1) += hydro_forces.vel[1]->item(parID);
// 	avg_force(2,1) += hydro_forces.vel[2]->item(parID);
//
// 	avg_torque(0,0) += hydro_torque.press[0]->item(parID);
// 	avg_torque(1,0) += hydro_torque.press[1]->item(parID);
// 	avg_torque(2,0) += hydro_torque.press[2]->item(parID);
// 	avg_torque(0,1) += hydro_torque.vel[0]->item(parID);
// 	avg_torque(1,1) += hydro_torque.vel[1]->item(parID);
// 	avg_torque(2,1) += hydro_torque.vel[2]->item(parID);
//
//         MyFile << t_it -> time() << "," << parID << "," << Nmax << "," << xp << "," << yp << "," << zp
// 	                                                        << "," << vx << "," << vy << "," << vz
// 					                        << "," << hydro_forces.press[0]->item(parID)
//                                                                 << "," << hydro_forces.press[1]->item(parID)
//                                                                 << "," << hydro_forces.press[2]->item(parID)
//                                                                 << "," << hydro_forces.vel[0]->item(parID)
//                                                                 << "," << hydro_forces.vel[1]->item(parID)
//                                                                 << "," << hydro_forces.vel[2]->item(parID)
// 					                        << "," << hydro_torque.press[0]->item(parID)
//                                                                 << "," << hydro_torque.press[1]->item(parID)
//                                                                 << "," << hydro_torque.press[2]->item(parID)
//                                                                 << "," << hydro_torque.vel[0]->item(parID)
//                                                                 << "," << hydro_torque.vel[1]->item(parID)
//                                                                 << "," << hydro_torque.vel[2]->item(parID) << endl;
//      }
//   }
//   if (my_rank == 0) cout << "Average force on particles: " << (avg_force(0,0) + avg_force(0,1))/double(Npart)
//                                                     << "," << (avg_force(1,0) + avg_force(1,1))/double(Npart)
//                                                     << "," << (avg_force(2,0) + avg_force(2,1))/double(Npart) <<endl;
//   if (my_rank == 0) cout << "Average torque on particles: " << (avg_torque(0,0) + avg_torque(0,1))/double(Npart)
//                                                      << "," << (avg_torque(1,0) + avg_torque(1,1))/double(Npart)
//                                                      << "," << (avg_torque(2,0) + avg_torque(2,1))/double(Npart) <<endl;
//   if (my_rank == 0) MyFile.close( ) ;
// }
//
//---------------------------------------------------------------------------
void
DDS_NavierStokes:: compute_surface_points_on_cube(size_t const& Np)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_NavierStokes:: compute_surface_points_on_cube" ) ;

  // Structure of particle input data
  SurfaceDiscretize surface = GLOBAL_EQ->get_surface();

  double dp = 2./((double)Np);

  if (dim == 3) {
     // Generating discretization on surface
     doubleVector lsp(Np,0.);
     for (size_t i=0; i<Np; i++) lsp(i) = -1. + dp*(double(i)+0.5);

     size_t cntr = 0;

     for (size_t i=0; i<Np; i++) {
        for (size_t j=0; j<Np; j++) {
           //Front
           surface.coordinate[0]->set_item(cntr,lsp(i));
           surface.coordinate[1]->set_item(cntr,lsp(j));
           surface.coordinate[2]->set_item(cntr,1.);
           surface.area->set_item(cntr,dp*dp);
           surface.normal[2]->set_item(cntr,1.);
           //Behind
           surface.coordinate[0]->set_item(Np*Np+cntr,lsp(i));
           surface.coordinate[1]->set_item(Np*Np+cntr,lsp(j));
           surface.coordinate[2]->set_item(Np*Np+cntr,-1.);
           surface.area->set_item(Np*Np+cntr,dp*dp);
           surface.normal[2]->set_item(Np*Np+cntr,-1.);
           //Top
           surface.coordinate[0]->set_item(2*Np*Np+cntr,lsp(i));
           surface.coordinate[2]->set_item(2*Np*Np+cntr,lsp(j));
           surface.coordinate[1]->set_item(2*Np*Np+cntr,1.);
           surface.area->set_item(2*Np*Np+cntr,dp*dp);
           surface.normal[1]->set_item(2*Np*Np+cntr,1.);
           //Bottom
           surface.coordinate[0]->set_item(3*Np*Np+cntr,lsp(i));
           surface.coordinate[2]->set_item(3*Np*Np+cntr,lsp(j));
           surface.coordinate[1]->set_item(3*Np*Np+cntr,-1.);
           surface.area->set_item(3*Np*Np+cntr,dp*dp);
           surface.normal[1]->set_item(3*Np*Np+cntr,-1.);
           //Right
           surface.coordinate[1]->set_item(4*Np*Np+cntr,lsp(i));
           surface.coordinate[2]->set_item(4*Np*Np+cntr,lsp(j));
           surface.coordinate[0]->set_item(4*Np*Np+cntr,1.);
           surface.area->set_item(4*Np*Np+cntr,dp*dp);
           surface.normal[0]->set_item(4*Np*Np+cntr,1.);
           //Left
           surface.coordinate[1]->set_item(5*Np*Np+cntr,lsp(i));
           surface.coordinate[2]->set_item(5*Np*Np+cntr,lsp(j));
           surface.coordinate[0]->set_item(5*Np*Np+cntr,-1.);
           surface.area->set_item(5*Np*Np+cntr,dp*dp);
           surface.normal[0]->set_item(5*Np*Np+cntr,-1.);

           cntr++;
        }
     }
  } else if (dim == 2) {
     // Generating discretization on surface
     double lsp=0.;
     for (size_t i=0; i<Np; i++) {
         lsp = -1. + dp*((double)i+0.5);
      	//Bottom
      	surface.coordinate[0]->set_item(i,lsp);
      	surface.coordinate[1]->set_item(i,-1.);
      	surface.area->set_item(i,dp);
      	surface.normal[1]->set_item(i,-1.);
      	//Top
      	surface.coordinate[0]->set_item(Np+i,lsp);
      	surface.coordinate[1]->set_item(Np+i,1.);
      	surface.area->set_item(Np+i,dp);
      	surface.normal[1]->set_item(Np+i,1.);
      	//Left
      	surface.coordinate[0]->set_item(2*Np+i,-1.);
      	surface.coordinate[1]->set_item(2*Np+i,lsp);
      	surface.area->set_item(2*Np+i,dp);
      	surface.normal[0]->set_item(2*Np+i,-1.);
      	//Right
      	surface.coordinate[0]->set_item(3*Np+i,1.);
      	surface.coordinate[1]->set_item(3*Np+i,lsp);
      	surface.area->set_item(3*Np+i,dp);
      	surface.normal[0]->set_item(3*Np+i,1.);
     }
  }
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: compute_surface_points_on_cylinder(class doubleVector& k
                                                    , size_t const& Nring)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_NavierStokes:: compute_surface_points_on_cylinder" ) ;

  // Structure of particle input data
  SurfaceDiscretize surface = GLOBAL_EQ->get_surface();

  // Radius of the rings in lamber projection plane
  doubleVector Rring(Nring,0.);

  Rring(Nring-1) = 1.;

  size_t maxby2 = (size_t) k(Nring-1);

  if (dim == 3) {
     // Calculation for all rings except at the pole
     for (int i=(int)Nring-1; i>0; --i) {
        double Ri = Rring(i);
        Rring(i-1) = MAC::sqrt(k(i-1)/k(i))*Rring(i);
        Rring(i) = (Rring(i) + Rring(i-1))/2.;
        double d_theta = 2.*MAC::pi()/(k(i)-k(i-1));
        // Theta initialize as 1% of the d_theta,
        // so there would be no chance of point overlap with mesh gridlines
        double theta = 0.01*d_theta;
        for (int j=(int)k(i-1); j<k(i); j++) {
           // For top disk
           theta = theta + d_theta;
           surface.coordinate[0]->set_item(j,Rring(i)*MAC::cos(theta));
           surface.coordinate[1]->set_item(j,Rring(i)*MAC::sin(theta));
           surface.coordinate[2]->set_item(j,1.);
           surface.area->set_item(j,0.5*d_theta*(pow(Ri,2)-pow(Rring(i-1),2)));
      	  // For bottom disk
      	  surface.coordinate[0]->set_item(maxby2+j,Rring(i)*MAC::cos(theta));
           surface.coordinate[1]->set_item(maxby2+j,Rring(i)*MAC::sin(theta));
           surface.coordinate[2]->set_item(maxby2+j,-1.);
           surface.area->set_item(maxby2+j,0.5*d_theta*(pow(Ri,2)-pow(Rring(i-1),2)));
           // Create surface normal vectors
           surface.normal[0]->set_item(j,0.);
      	  surface.normal[1]->set_item(j,0.);
      	  surface.normal[2]->set_item(j,1.);
           surface.normal[0]->set_item(maxby2+j,0.);
      	  surface.normal[1]->set_item(maxby2+j,0.);
      	  surface.normal[2]->set_item(maxby2+j,-1.);
        }
     }

     // Calculation at the ring on pole (i=0)
     double Ri = Rring(0);
     Rring(0) = Rring(0)/2.;
     double d_theta = 2.*MAC::pi()/(k(0));
     // Theta initialize as 1% of the d_theta,
     // so there would be no chance of point overlap with mesh gridlines
     double theta = 0.01*d_theta;
     if (k(0)>1) {
        for (int j=0; j < k(0); j++) {
           // For top disk
           theta = theta + d_theta;
           surface.coordinate[0]->set_item(j,Rring(0)*MAC::cos(theta));
           surface.coordinate[1]->set_item(j,Rring(0)*MAC::sin(theta));
           surface.coordinate[2]->set_item(j,1.);
           surface.area->set_item(j,0.5*d_theta*pow(Ri,2));
           // For bottom disk
           surface.coordinate[0]->set_item(maxby2+j,Rring(0)*MAC::cos(theta));
           surface.coordinate[1]->set_item(maxby2+j,Rring(0)*MAC::sin(theta));
           surface.coordinate[2]->set_item(maxby2+j,-1.);
           surface.area->set_item(maxby2+j,0.5*d_theta*pow(Ri,2.));
      	  // Create surface normal vectors
      	  surface.normal[0]->set_item(j,0.);
      	  surface.normal[1]->set_item(j,0.);
      	  surface.normal[2]->set_item(j,1.);
      	  surface.normal[0]->set_item(maxby2+j,0.);
      	  surface.normal[1]->set_item(maxby2+j,0.);
      	  surface.normal[2]->set_item(maxby2+j,-1.);
        }
     } else {
        // For top disk
        surface.coordinate[0]->set_item(0,0.);
        surface.coordinate[1]->set_item(0,0.);
        surface.coordinate[2]->set_item(0,1.);
        surface.area->set_item(0,0.5*d_theta*pow(Ri,2.));
        // For bottom disk
        surface.coordinate[0]->set_item(maxby2,0.);
        surface.coordinate[1]->set_item(maxby2,0.);
        surface.coordinate[2]->set_item(maxby2,-1.);
        surface.area->set_item(maxby2,0.5*d_theta*pow(Ri,2.));
        // Create surface normal vectors
        surface.normal[0]->set_item(0,0.);
        surface.normal[1]->set_item(0,0.);
        surface.normal[2]->set_item(0,1.);
        surface.normal[0]->set_item(maxby2,0.);
        surface.normal[1]->set_item(maxby2,0.);
        surface.normal[2]->set_item(maxby2,-1.);
     }

     // Generating one ring of points on cylindrical surface
     // Can be used to calculate stress on whole surface by
     // a constant shift of points

     // Estimating number of points on cylindrical surface
     int pts_1_ring = (int)(k(Nring-1) - k(Nring-2));
     int cyl_rings = (int)round(2./(1-MAC::sqrt(k(Nring-2)/k(Nring-1))));
     double cell_area = 2.*MAC::pi()/(double(pts_1_ring))*(2./(double(cyl_rings)));

     d_theta = 2.*MAC::pi()/pts_1_ring;
     for (int j=0; j<cyl_rings; j++) {
        theta = 0.01*d_theta;
        for (int ij=0; ij<pts_1_ring; ij++) {
           theta = theta + d_theta;
           int n = 2*(int)k(Nring-1) + j*pts_1_ring + ij;
           surface.coordinate[0]->set_item(n,MAC::cos(theta));
           surface.coordinate[1]->set_item(n,MAC::sin(theta));
           surface.coordinate[2]->set_item(n,-1.+ 2.*(j+0.5)/(double(cyl_rings)));
           surface.area->set_item(n,cell_area);
           surface.normal[0]->set_item(n,MAC::cos(theta));
           surface.normal[1]->set_item(n,MAC::sin(theta));
           surface.normal[2]->set_item(n,0.);
        }
     }
     // write_surface_discretization(2*maxby2 + pts_1_ring*cyl_rings);
  }
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: compute_surface_points_on_sphere(class doubleVector& eta,
                                                    class doubleVector& k,
                                                    class doubleVector& Rring,
                                                    size_t const& Nring)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_NavierStokes:: compute_surface_points_on_sphere" ) ;

  // Structure of particle input data
  SurfaceDiscretize surface = GLOBAL_EQ->get_surface();

  if (dim == 3) {
     size_t maxby2 = (size_t) k(Nring-1);
     // Calculation for all rings except at the pole
     for (int i=(int)Nring-1; i>0; --i) {
        double Ri = Rring(i);
        Rring(i) = (Rring(i) + Rring(i-1))/2.;
        eta(i) = (eta(i) + eta(i-1))/2.;
        double d_theta = 2.*MAC::pi()/(k(i)-k(i-1));
        // Theta initialize as 1% of the d_theta,
        // so there would be no chance of point overlap with mesh gridlines
        double theta = 0.01*d_theta;
        for (int j=(int)k(i-1); j<k(i); j++) {
           theta = theta + d_theta;
           if (pole_loc == 2) {
              surface.coordinate[0]->set_item(j,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate[1]->set_item(j,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate[2]->set_item(j,MAC::cos(eta(i)));
              surface.area->set_item(j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
              // For second half of sphere
              surface.coordinate[0]->set_item(maxby2+j,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate[1]->set_item(maxby2+j,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate[2]->set_item(maxby2+j,-MAC::cos(eta(i)));
              surface.area->set_item(maxby2+j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
           } else if (pole_loc == 1) {
              surface.coordinate[2]->set_item(j,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate[0]->set_item(j,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate[1]->set_item(j,MAC::cos(eta(i)));
              surface.area->set_item(j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
              // For second half of sphere
              surface.coordinate[2]->set_item(maxby2+j,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate[0]->set_item(maxby2+j,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate[1]->set_item(maxby2+j,-MAC::cos(eta(i)));
              surface.area->set_item(maxby2+j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
           } else if (pole_loc == 0) {
              surface.coordinate[1]->set_item(j,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate[2]->set_item(j,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate[0]->set_item(j,MAC::cos(eta(i)));
              surface.area->set_item(j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
              // For second half of sphere
              surface.coordinate[1]->set_item(maxby2+j,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate[2]->set_item(maxby2+j,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate[0]->set_item(maxby2+j,-MAC::cos(eta(i)));
              surface.area->set_item(maxby2+j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
           }
	   // Create surface normal vectors
	   surface.normal[0]->set_item(j,surface.coordinate[0]->item(j));
	   surface.normal[1]->set_item(j,surface.coordinate[1]->item(j));
	   surface.normal[2]->set_item(j,surface.coordinate[2]->item(j));
	   surface.normal[0]->set_item(maxby2+j,surface.coordinate[0]->item(maxby2+j));
	   surface.normal[1]->set_item(maxby2+j,surface.coordinate[1]->item(maxby2+j));
	   surface.normal[2]->set_item(maxby2+j,surface.coordinate[2]->item(maxby2+j));

        }
     }

     // Calculation at the ring on pole (i=0)
     double Ri = Rring(0);
     Rring(0) = Rring(0)/2.;
     eta(0) = eta(0)/2.;
     double d_theta = 2.*MAC::pi()/(k(0));
     // Theta initialize as 1% of the d_theta, so there would be no chance of point overlap with mesh gridlines
     double theta = 0.01*d_theta;
     if (k(0)>1) {
        for (int j=0; j < k(0); j++) {
           theta = theta + d_theta;
           if (pole_loc == 2) {
              surface.coordinate[0]->set_item(j,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate[1]->set_item(j,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate[2]->set_item(j,MAC::cos(eta(0)));
              surface.area->set_item(j,0.5*d_theta*pow(Ri,2.));
              // For second half of sphere
              surface.coordinate[0]->set_item(maxby2+j,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate[1]->set_item(maxby2+j,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate[2]->set_item(maxby2+j,-MAC::cos(eta(0)));
              surface.area->set_item(maxby2+j,0.5*d_theta*pow(Ri,2.));
           } else if (pole_loc == 1) {
              surface.coordinate[2]->set_item(j,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate[0]->set_item(j,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate[1]->set_item(j,MAC::cos(eta(0)));
              surface.area->set_item(j,0.5*d_theta*pow(Ri,2.));
              // For second half of sphere
              surface.coordinate[2]->set_item(maxby2+j,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate[0]->set_item(maxby2+j,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate[1]->set_item(maxby2+j,-MAC::cos(eta(0)));
              surface.area->set_item(maxby2+j,0.5*d_theta*pow(Ri,2.));
           } else if (pole_loc == 0) {
              surface.coordinate[1]->set_item(j,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate[2]->set_item(j,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate[0]->set_item(j,MAC::cos(eta(0)));
              surface.area->set_item(j,0.5*d_theta*pow(Ri,2.));
              // For second half of sphere
              surface.coordinate[1]->set_item(maxby2+j,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate[2]->set_item(maxby2+j,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate[0]->set_item(maxby2+j,-MAC::cos(eta(0)));
              surface.area->set_item(maxby2+j,0.5*d_theta*pow(Ri,2.));
           }
	   // Create surface normal vectors
	   surface.normal[0]->set_item(j,surface.coordinate[0]->item(j));
	   surface.normal[1]->set_item(j,surface.coordinate[1]->item(j));
	   surface.normal[2]->set_item(j,surface.coordinate[2]->item(j));
	   surface.normal[0]->set_item(maxby2+j,surface.coordinate[0]->item(maxby2+j));
	   surface.normal[1]->set_item(maxby2+j,surface.coordinate[1]->item(maxby2+j));
	   surface.normal[2]->set_item(maxby2+j,surface.coordinate[2]->item(maxby2+j));
        }
     } else {
        if (pole_loc == 2) {
           surface.coordinate[0]->set_item(0,0.);
           surface.coordinate[1]->set_item(0,0.);
           surface.coordinate[2]->set_item(0,1.);
           surface.area->set_item(0,0.5*d_theta*pow(Ri,2.));
           // For second half of sphere
           surface.coordinate[0]->set_item(maxby2,0.);
           surface.coordinate[1]->set_item(maxby2,0.);
           surface.coordinate[2]->set_item(maxby2,-1.);
           surface.area->set_item(maxby2,0.5*d_theta*pow(Ri,2.));
        } else if (pole_loc == 1) {
           surface.coordinate[2]->set_item(0,0.);
           surface.coordinate[0]->set_item(0,0.);
           surface.coordinate[1]->set_item(0,1.);
           surface.area->set_item(0,0.5*d_theta*pow(Ri,2.));
           // For second half of sphere
           surface.coordinate[2]->set_item(maxby2,0.);
           surface.coordinate[0]->set_item(maxby2,0.);
           surface.coordinate[1]->set_item(maxby2,-1.);
           surface.area->set_item(maxby2,0.5*d_theta*pow(Ri,2.));
        } else if (pole_loc == 0) {
           surface.coordinate[1]->set_item(0,0.);
           surface.coordinate[2]->set_item(0,0.);
           surface.coordinate[0]->set_item(0,1.);
           surface.area->set_item(0,0.5*d_theta*pow(Ri,2.));
           // For second half of sphere
           surface.coordinate[1]->set_item(maxby2,0.);
           surface.coordinate[2]->set_item(maxby2,0.);
           surface.coordinate[0]->set_item(maxby2,-1.);
           surface.area->set_item(maxby2,0.5*d_theta*pow(Ri,2.));
        }
        // Create surface normal vectors
        surface.normal[0]->set_item(0,surface.coordinate[0]->item(0));
        surface.normal[1]->set_item(0,surface.coordinate[1]->item(0));
        surface.normal[2]->set_item(0,surface.coordinate[2]->item(0));
        surface.normal[0]->set_item(maxby2,surface.coordinate[0]->item(maxby2));
        surface.normal[1]->set_item(maxby2,surface.coordinate[1]->item(maxby2));
        surface.normal[2]->set_item(maxby2,surface.coordinate[2]->item(maxby2));
     }
  } else if (dim == 2) {
     double d_theta = 2.*MAC::pi()/(double(Nring));
     double theta = 0.01*d_theta;
     for (int j=0; j < (int) Nring; j++) {
        theta = theta + d_theta;
        surface.coordinate[0]->set_item(j,MAC::cos(theta));
        surface.coordinate[1]->set_item(j,MAC::sin(theta));
        surface.area->set_item(j,d_theta);
        // Create surface normal vectors
        surface.normal[0]->set_item(j,surface.coordinate[0]->item(j));
        surface.normal[1]->set_item(j,surface.coordinate[1]->item(j));
     }
  }
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: generate_surface_discretization()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_NavierStokes:: generate_surface_discretization" ) ;

  size_t kmax = (int) Npoints;

  if ( level_set_type == "Sphere" ) {
     // Reference paper: Becker and Becker, A general rule for disk and hemisphere partition into
     // equal-area cells, Computational Geometry 45 (2012) 275-283
     double eta_temp = MAC::pi()/2.;
     double k_temp = (double) kmax;
     double Ro_temp = MAC::sqrt(2);
     double Rn_temp = MAC::sqrt(2);
     size_t cntr = 0;

     // Estimating the number of rings on the hemisphere
     while (k_temp > double(Pmin+2)) {
        eta_temp = eta_temp - 2./ar*MAC::sqrt(MAC::pi()/k_temp)*MAC::sin(eta_temp/2.);
        Rn_temp = 2.*MAC::sin(eta_temp/2.);
        k_temp = round(k_temp*pow(Rn_temp/Ro_temp,2.));
        Ro_temp = Rn_temp;
        cntr++;
     }

     size_t Nrings = cntr+1;

     // Summation of total discretized points with increase in number of rings radially
     doubleVector k(Nrings,0.);
     // Zenithal angle for the sphere
     doubleVector eta(Nrings,0.);
     // Radius of the rings in lamber projection plane
     doubleVector Rring(Nrings,0.);

     // Assigning the maximum number of discretized points to the last element of the array
     k(Nrings-1) = (double) kmax;
     // Zenithal angle for the last must be pi/2.
     eta(Nrings-1) = MAC::pi()/2.;
     // Radius of last ring in lamber projection plane
     Rring(Nrings-1) = MAC::sqrt(2.);

     for (int i=int(Nrings)-2; i>=0; --i) {
        eta(i) = eta(i+1) - 2./ar*MAC::sqrt(MAC::pi()/k(i+1))*MAC::sin(eta(i+1)/2.);
        Rring(i) = 2.*MAC::sin(eta(i)/2.);
        k(i) = round(k(i+1)*pow(Rring(i)/Rring(i+1),2.));
        if (i==0) k(0) = (double) Pmin;
     }

     // Discretize the particle surface into approximate equal area cells
     if (dim == 3) {
        compute_surface_points_on_sphere(eta, k, Rring, Nrings);
        //	write_surface_discretization(2*kmax);
     } else {
        compute_surface_points_on_sphere(eta, k, Rring, kmax);
        //	write_surface_discretization(kmax);
     }
  } else if (level_set_type == "Cube") {
     compute_surface_points_on_cube(kmax);
     //  write_surface_discretization((size_t)(6*pow(kmax,2)));
  } else if (level_set_type == "Cylinder") {
     // Reference paper: Becker and Becker, A general rule for disk and hemisphere partition into
     // equal-area cells, Computational Geometry 45 (2012) 275-283

     double p = MAC::pi()/ar;
     double k_temp = (double) kmax;
     size_t cntr = 0;

     // Estimating the number of rings on either of the disc
     while (k_temp > double(Pmin+2)) {
        k_temp = round(pow(MAC::sqrt(k_temp) - MAC::sqrt(p),2.));
        cntr++;
     }

     size_t Nrings = cntr+1;

     // Summation of total discretized points with increase in number of rings radially
     doubleVector k(Nrings,0.);
     // Assigning the maximum number of discretized points to the last element of the array
     k(Nrings-1) = (double) kmax;

     for (int i=(int)Nrings-2; i>=0; --i) {
        k(i) = round(pow(MAC::sqrt(k(i+1)) - MAC::sqrt(p),2.));
        if (i==0) k(0) = (double) Pmin;
     }

     if (dim == 3) {
        compute_surface_points_on_cylinder(k, Nrings);
     }
  }
}

// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: ghost_points_generation( FV_DiscreteField* FF, class doubleArray2D& point, class size_t_array2D& i0, int const& sign, size_t const& comp,size_t const& major_dir, class boolArray2D& point_in_domain, class doubleVector& rotated_vector )
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: ghost_points_generation" ) ;
//
//   intVector i0_temp(2,0);
//   size_t i0_t;
//
//   // Ghost points in i for the calculation of i-derivative of field
//   i0_temp(0) = (sign == 1) ? (int(i0(0,major_dir)) + 1*sign) : (int(i0(0,major_dir)) + 0*sign);
//   i0_temp(1) = (sign == 1) ? (int(i0(0,major_dir)) + 2*sign) : (int(i0(0,major_dir)) + 1*sign);
//
//   if ((i0_temp(0) >= 0) && (i0_temp(0) < (int)FF->get_local_nb_dof(comp,major_dir)))
//      point(1,major_dir) = FF->get_DOF_coordinate(i0_temp(0), comp, major_dir);
//   if ((i0_temp(1) >= 0) && (i0_temp(1) < (int)FF->get_local_nb_dof(comp,major_dir)))
//      point(2,major_dir) = FF->get_DOF_coordinate(i0_temp(1), comp, major_dir);
//
//   if (MAC::abs(point(0,major_dir)-point(1,major_dir)) < MAC::abs(point(1,major_dir)-point(2,major_dir))) {
//      i0_temp(0) = (sign == 1) ? (int(i0(0,major_dir)) + 2*sign) : (int(i0(0,major_dir)) + 1*sign);
//      i0_temp(1) = (sign == 1) ? (int(i0(0,major_dir)) + 3*sign) : (int(i0(0,major_dir)) + 2*sign);
//
//      if ((i0_temp(0) >= 0) && (i0_temp(0) < (int)FF->get_local_nb_dof(comp,major_dir)))
//         point(1,major_dir) = FF->get_DOF_coordinate(i0_temp(0), comp, major_dir);
//      if ((i0_temp(1) >= 0) && (i0_temp(1) < (int)FF->get_local_nb_dof(comp,major_dir)))
//         point(2,major_dir) = FF->get_DOF_coordinate(i0_temp(1), comp, major_dir);
//   }
//
//   // Velocity field
//   if (FF == UF) {
//      if ((i0_temp(0) < 0) || (i0_temp(0) >= (int)FF->get_local_nb_dof(comp,major_dir))) {
//         point_in_domain(0,major_dir) = 0;
//      } else {
//         point_in_domain(0,major_dir) = 1;
//      }
//
//      if ((i0_temp(1) < 0) || (i0_temp(1) >= (int)FF->get_local_nb_dof(comp,major_dir))) {
//         point_in_domain(1,major_dir) = 0;
//      } else {
//         point_in_domain(1,major_dir) = 1;
//      }
//
//      i0(1,major_dir) = i0_temp(0);
//      i0(2,major_dir) = i0_temp(1);
//   // Pressure field
//   } else {
//      double t1 = (point(1,major_dir) - point(0,major_dir))/rotated_vector(major_dir);
//      double t2 = (point(2,major_dir) - point(0,major_dir))/rotated_vector(major_dir);
//
//      for (size_t dir=0; dir<dim; dir++) {
//         if (dir != major_dir) {
//            point(1,dir) = point(0,dir) + rotated_vector(dir)*t1;
//            point(2,dir) = point(0,dir) + rotated_vector(dir)*t2;
//         }
//      }
//
//      boolArray2D found(2,dim,true);
//
//      for (size_t dir=0; dir<dim; dir++) {
//         if (dir == major_dir) {
//            // Checking the points in domain or not
//            if ((i0_temp(0) < 0) || (i0_temp(0) >= (int)FF->get_local_nb_dof(0,dir))) {
//               found(0,dir) = 0;
//            }
//            if ((i0_temp(1) < 0) || (i0_temp(1) >= (int)FF->get_local_nb_dof(0,dir))) {
//               found(1,dir) = 0;
//            }
//            i0(1,dir) = i0_temp(0);
//            i0(2,dir) = i0_temp(1);
//         } else {
//            found(0,dir) = FV_Mesh::between(FF->get_DOF_coordinates_vector(0,dir), point(1,dir), i0_t);
//            if (found(0,dir)) i0(1,dir) = i0_t;
//            found(1,dir) = FV_Mesh::between(FF->get_DOF_coordinates_vector(0,dir), point(2,dir), i0_t);
//            if (found(1,dir)) i0(2,dir) = i0_t;
//         }
//      }
//
//      // Checking the ghost points in domain or not
//      point_in_domain(0,0) = (dim == 2) ? found(0,0) && found(0,1) :
//                                          found(0,0) && found(0,1) && found(0,2) ;
//      point_in_domain(0,1) = (dim == 2) ? found(1,0) && found(1,1) :
//                                          found(1,0) && found(1,1) && found(1,2) ;
//   }
// }
//
//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: compute_pressure_force_on_particle(size_t const& parID, size_t const& Np, FV_TimeIterator const* t_it )
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: compute_pressure_force_on_particle" ) ;
//
//   if (PressureStressOrder == "first") {
//      first_order_pressure_stress(parID, Np, t_it );
//   } else if (PressureStressOrder == "second") {
//      second_order_pressure_stress(parID, Np, t_it );
//   } else if (PressureStressOrder == "second_withNeumannBC") {
//      second_order_pressure_stress_withNeumannBC(parID, Np, t_it );
//   }
// }
//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: compute_velocity_force_on_particle(size_t const& parID, size_t const& Np, FV_TimeIterator const* t_it )
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: compute_velocity_force_on_particle" ) ;
//
//   if (ViscousStressOrder == "first") {
//      first_order_viscous_stress(parID, Np, t_it );
//   } else if (ViscousStressOrder == "second") {
//      second_order_viscous_stress(parID, Np, t_it );
//   }
// }
//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: second_order_viscous_stress(size_t const& parID, size_t const& Np, FV_TimeIterator const* t_it )
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: second_order_viscous_stress" ) ;
//
//   size_t i0_temp;
//   double dfdx=0.,dfdy=0., dfdz=0.;
//
//   // Structure of particle input data
//   PartInput solid = GLOBAL_EQ->get_solid(0);
//   // Structure of particle surface input data
//   SurfaceDiscretize surface = GLOBAL_EQ->get_surface();
//
//   // comp won't matter as the particle position is independent of comp
//   double xp = solid.coord[0]->item(parID);
//   double yp = solid.coord[1]->item(parID);
//   double zp = solid.coord[2]->item(parID);
//   double ri = solid.size->item(parID);
//
//   doubleArray2D point(3,3,0);
//   doubleArray2D fini(3,3,0);
//   doubleArray2D stress(Np,6,0);         //xx,yy,zz,xy,yz,zx
//   doubleArray2D level_set(dim,2,1.);
//   size_t_array2D in_parID(dim,2,0);         //Store particle ID if level_set becomes negative
//   boolArray2D found(dim,3,false);
//   size_t_array2D i0(3,3,0);
//   vector<double> net_vel(3,0.);
//   doubleArray2D surface_force(Np,3,0);
//   doubleArray2D surface_point(Np,3,0);
//
//   size_t_vector min_unknown_index(dim,0);
//   size_t_vector max_unknown_index(dim,0);
//
//   doubleVector Dmin(dim,0);
//   doubleVector Dmax(dim,0);
//   doubleVector rotated_coord(3,0);
//   doubleVector rotated_normal(3,0);
//   intVector sign(dim,0);
//   boolArray2D point_in_domain(2,dim,1);
//
//   // Structure of particle input data
//   PartForces hydro_forces = GLOBAL_EQ->get_forces(0);
//   PartForces hydro_torque = GLOBAL_EQ->get_torque(0);
//
//   // Force vector
//   doubleVector force(3,0);
//   doubleVector torque(3,0);
//
//   for (size_t i=0;i<Np;i++) {
//      for (size_t comp=0;comp<nb_comps[1];comp++) {
//         // Get local min and max indices
// 	// Get min and max coordinates in the current processor
//         for (size_t l=0;l<dim;++l) {
//            min_unknown_index(l) = UF->get_min_index_unknown_handled_by_proc( comp, l );
//            max_unknown_index(l) = UF->get_max_index_unknown_handled_by_proc( comp, l );
// 	   Dmin(l) = UF->primary_grid()->get_min_coordinate_on_current_processor(l);
// 	   Dmax(l) = UF->primary_grid()->get_max_coordinate_on_current_processor(l);
//         }
//
// 	// Rotating surface points
// 	rotated_coord(0) = ri*surface.coordinate[0]->item(i);
// 	rotated_coord(1) = ri*surface.coordinate[1]->item(i);
// 	rotated_coord(2) = ri*surface.coordinate[2]->item(i);
//
//         rotation_matrix(parID,rotated_coord,comp,1);
//
//         point(0,0) = xp + rotated_coord(0);
//         point(0,1) = yp + rotated_coord(1);
//         point(0,2) = zp + rotated_coord(2);
//
//         // Rotating surface normal
// 	rotated_normal(0) = surface.normal[0]->item(i);
// 	rotated_normal(1) = surface.normal[1]->item(i);
// 	rotated_normal(2) = surface.normal[2]->item(i);
//
// 	rotation_matrix(parID,rotated_normal,comp,1);
//
//         for (size_t dir=0;dir<dim;dir++) {
// 	   // PBC on rotated surface points
// 	   if (is_periodic[1][dir]) {
//               double isize = UF->primary_grid()->get_main_domain_max_coordinate(dir) - UF->primary_grid()->get_main_domain_min_coordinate(dir);
//               double imin = UF->primary_grid()->get_main_domain_min_coordinate(dir);
//               point(0,dir) = point(0,dir) - MAC::floor((point(0,dir)-imin)/isize)*isize;
//            }
//            // Finding the grid indexes next to ghost points
//            found(dir,0) = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,dir), point(0,dir), i0_temp);
//            if (found(dir,0) == 1) i0(0,dir) = i0_temp;
// 	}
//
// 	// Accessing the smallest grid size in domain
//         double dh = UF->primary_grid()->get_smallest_grid_size();
//
//         bool status = (dim==2) ? ((point(0,0) > Dmin(0)) && (point(0,0) <= Dmax(0)) && (point(0,1) > Dmin(1)) && (point(0,1) <= Dmax(1))) :
//                                  ((point(0,0) > Dmin(0)) && (point(0,0) <= Dmax(0)) && (point(0,1) > Dmin(1)) && (point(0,1) <= Dmax(1))
//                                                                                     && (point(0,2) > Dmin(2)) && (point(0,2) <= Dmax(2)));
//         double threshold = pow(loc_thres,0.5)*dh;
//
//         if (status) {
//            for (size_t dir=0;dir<dim;dir++) {
//               sign(dir) = (rotated_normal(dir) > 0.) ? 1 : -1;
//
//               // Ghost points in i for the calculation of i-derivative of field
//               ghost_points_generation( UF, point, i0, sign(dir), comp, dir, point_in_domain, rotated_normal);
// 	      // Assuming all ghost points are in fluid
//               level_set(dir,0) = 1.; level_set(dir,1) = 1.;
// 	   }
//
//            // Checking all the ghost points in the solid/fluid, and storing the parID if present in solid
//            for (size_t m=0;m<Npart;m++) {
// 	      // In X
//               if (level_set(0,0) > threshold) {
//                  level_set(0,0) = level_set_function(UF,m,comp,point(1,0),point(0,1),point(0,2),level_set_type,1);
//                  level_set(0,0) *= solid.inside->item(m);
//                  if (level_set(0,0) < threshold) in_parID(0,0) = m;
//               }
//               if (level_set(0,1) > threshold) {
//                  level_set(0,1) = level_set_function(UF,m,comp,point(2,0),point(0,1),point(0,2),level_set_type,1);
//                  level_set(0,1) *= solid.inside->item(m);
//                  if (level_set(0,1) < threshold) in_parID(0,1) = m;
//               }
// 	      // In Y
//               if (level_set(1,0) > threshold) {
//                  level_set(1,0) = level_set_function(UF,m,comp,point(0,0),point(1,1),point(0,2),level_set_type,1);
//                  level_set(1,0) *= solid.inside->item(m);
//                  if (level_set(1,0) < threshold) in_parID(1,0) = m;
//               }
//               if (level_set(1,1) > threshold) {
//                  level_set(1,1) = level_set_function(UF,m,comp,point(0,0),point(2,1),point(0,2),level_set_type,1);
//                  level_set(1,1) *= solid.inside->item(m);
//                  if (level_set(1,1) < threshold) in_parID(1,1) = m;
//               }
// 	      // In Z
// 	      if (dim == 3) {
//                  if (level_set(2,0) > threshold) {
//                     level_set(2,0) = level_set_function(UF,m,comp,point(0,0),point(0,1),point(1,2),level_set_type,1);
//                     level_set(2,0) *= solid.inside->item(m);
//                     if (level_set(2,0) < threshold) in_parID(2,0) = m;
//                  }
//                  if (level_set(2,1) > threshold) {
//                     level_set(2,1) = level_set_function(UF,m,comp,point(0,0),point(0,1),point(2,2),level_set_type,1);
//                     level_set(2,1) *= solid.inside->item(m);
//                     if (level_set(2,1) < threshold) in_parID(2,1) = m;
//                  }
// 	      }
//            }
//
//            // Calculation of field variable on ghost point(0,0)
//            impose_solid_velocity_for_ghost(net_vel,comp,point(0,0),point(0,1),point(0,2),parID);
//            fini(0,0) = net_vel[comp];
//            fini(0,1) = net_vel[comp];
//            fini(0,2) = net_vel[comp];
//
//            // Calculation of field variable on ghost point(1,0)
//            if ((level_set(0,0) > threshold) && point_in_domain(0,0)) {
//               fini(1,0) = third_order_ghost_field_estimate(UF, comp, point(1,0), point(0,1), point(0,2), i0(1,0), i0(0,1), i0(0,2), 0, sign,0);
//            } else if (level_set(0,0) <= threshold) {
//               impose_solid_velocity_for_ghost(net_vel,comp,point(1,0),point(0,1),point(0,2),in_parID(0,0));
//               fini(1,0) = net_vel[comp];
//            }
//            // Calculation of field variable on ghost point(2,0)
//            if ((level_set(0,1) > threshold) && point_in_domain(1,0)) {
//               fini(2,0) = third_order_ghost_field_estimate(UF, comp, point(2,0), point(0,1), point(0,2), i0(2,0), i0(0,1), i0(0,2), 0, sign,0);
//            } else if (level_set(0,1) <= threshold) {
//               impose_solid_velocity_for_ghost(net_vel,comp,point(2,0),point(0,1),point(0,2),in_parID(0,1));
//               fini(2,0) = net_vel[comp];
//            }
//            // Calculation of field variable on ghost point(1,1)
//            if ((level_set(1,0) > threshold) && point_in_domain(0,1)) {
// 	      fini(1,1) = third_order_ghost_field_estimate(UF, comp, point(0,0), point(1,1), point(0,2), i0(0,0), i0(1,1), i0(0,2), 1, sign,0);
//            } else if (level_set(1,0) <= threshold) {
//               impose_solid_velocity_for_ghost(net_vel,comp,point(0,0),point(1,1),point(0,2),in_parID(1,0));
//               fini(1,1) = net_vel[comp];
// 	   }
//            // Calculation of field variable on ghost point(2,1)
//            if ((level_set(1,1) > threshold) && point_in_domain(1,1)) {
//               fini(2,1) = third_order_ghost_field_estimate(UF, comp, point(0,0), point(2,1), point(0,2), i0(0,0), i0(2,1), i0(0,2), 1, sign,0);
// 	   } else if (level_set(1,1) <= threshold) {
//               impose_solid_velocity_for_ghost(net_vel,comp,point(0,0),point(2,1),point(0,2),in_parID(1,1));
// 	      fini(2,1) = net_vel[comp];
//            }
//
// 	   if (dim == 3) {
//               // Calculation of field variable on ghost point(1,2)
//               if ((level_set(2,0) > threshold) && point_in_domain(0,2)) {
//                 fini(1,2) = third_order_ghost_field_estimate(UF, comp, point(0,0), point(0,1), point(1,2), i0(0,0), i0(0,1), i0(1,2), 2, sign, 0);
//               } else if (level_set(2,0) <= threshold) {
//                  impose_solid_velocity_for_ghost(net_vel,comp,point(0,0),point(0,1),point(1,2),in_parID(2,0));
//                  fini(1,2) = net_vel[comp];
//               }
//               // Calculation of field variable on ghost point(2,2)
//               if ((level_set(2,1) > threshold) && point_in_domain(1,2)) {
//                 fini(2,2) = third_order_ghost_field_estimate(UF, comp, point(0,0), point(0,1), point(2,2), i0(0,0), i0(0,1), i0(2,2), 2, sign, 0);
// 	      } else if (level_set(2,1) <= threshold) {
//                  impose_solid_velocity_for_ghost(net_vel,comp,point(0,0),point(0,1),point(2,2),in_parID(2,1));
// 	         fini(2,2) = net_vel[comp];
//               }
//  	   }
//
//            // Derivative in x
//            // Point 1 and 2 in computational domain
// 	   if (point_in_domain(0,0) && point_in_domain(1,0)) {
//               if ((level_set(0,0) > threshold) && (level_set(0,1) > threshold)) {
// 	         double dx1 = (point(1,0)-point(0,0));
//                  double dx2 = (point(2,0)-point(0,0));
//                  dfdx = mu*((fini(1,0) - fini(0,0))*dx2/dx1 - (fini(2,0) - fini(0,0))*dx1/dx2)/(dx2-dx1);
//               // Point 1 in fluid and 2 in the solid
//               } else if ((level_set(0,0) > threshold) && (level_set(0,1) <= threshold)) {
//                  double dx1 = (point(1,0)-point(0,0));
//                  dfdx = mu*(fini(1,0) - fini(0,0))/dx1;
//               // Point 1 is present in solid
//               } else if (level_set(0,0) <= threshold) {
// 	         double dx1 = (point(1,0)-point(0,0));
//                  dfdx = mu*(fini(1,0) - fini(0,0))/dx1;
// 	      }
//            // Point 1 in computational domain
// 	   } else if (point_in_domain(0,0) && !point_in_domain(1,0)) {
//               double dx1 = (point(1,0)-point(0,0));
//               dfdx = mu*(fini(1,0) - fini(0,0))/dx1;
//            // Particle close to wall
//            } else if (!point_in_domain(0,0)) {
//               i0(1,0) = (sign(0) == 1) ? (i0(0,0) + 1*sign(0)) : (i0(0,0) + 0*sign(0));
//               point(1,0) = UF->get_DOF_coordinate(i0(1,0), comp, 0);
//               fini(1,0) = third_order_ghost_field_estimate(UF, comp, point(1,0), point(0,1), point(0,2), i0(1,0), i0(0,1), i0(0,2), 0, sign,0);
//               double dx1 = (point(1,0)-point(0,0));
//               dfdx = mu*(fini(1,0) - fini(0,0))/dx1;
// 	   }
//
//            // Derivative in y
//            // Point 1 and 2 in computational domain
//            if (point_in_domain(0,1) && point_in_domain(1,1)) {
//               if ((level_set(1,0) > threshold) && (level_set(1,1) > threshold)) {
//                  double dy1 = (point(1,1)-point(0,1));
// 	         double dy2 = (point(2,1)-point(0,1));
//                  dfdy = mu*((fini(1,1) - fini(0,1))*dy2/dy1 - (fini(2,1) - fini(0,1))*dy1/dy2)/(dy2-dy1);
//               // Point 1 in fluid and 2 in the solid
//               } else if ((level_set(1,0) > threshold) && (level_set(1,1) <= threshold)) {
// 	         double dy1 = (point(1,1)-point(0,1));
//                  dfdy = mu*(fini(1,1) - fini(0,1))/dy1;
//               // Point 1 is present in solid
//               } else if (level_set(1,0) <= threshold) {
//                  double dy1 = (point(1,1)-point(0,1));
//                  dfdy = mu*(fini(1,1) - fini(0,1))/dy1;
// 	      }
//            // Point 1 in computational domain
// 	   } else if (point_in_domain(0,1) && !point_in_domain(1,1)) {
//               double dy1 = (point(1,1)-point(0,1));
//               dfdy = mu*(fini(1,1) - fini(0,1))/dy1;
//            // Particle close to wall
//            } else if (!point_in_domain(0,1)) {
//               i0(1,1) = (sign(1) == 1) ? (i0(0,1) + 1*sign(1)) : (i0(0,1) + 0*sign(1));
//               point(1,1) = UF->get_DOF_coordinate(i0(1,1), comp, 1);
// 	      fini(1,1) = third_order_ghost_field_estimate(UF, comp, point(0,0), point(1,1), point(0,2), i0(0,0), i0(1,1), i0(0,2), 1, sign,0);
//               double dy1 = (point(1,1)-point(0,1));
//               dfdy = mu*(fini(1,1) - fini(0,1))/dy1;
//            }
//
//            // Derivative in z
// 	   if (dim == 3) {
//               // Point 1 and 2 in computational domain
//               if (point_in_domain(0,2) && point_in_domain(1,2)) {
//                  if ((level_set(2,0) > threshold) && (level_set(2,1) > threshold)) {
//                     double dz1 = (point(1,2)-point(0,2));
//                     double dz2 = (point(2,2)-point(0,2));
//                     dfdz = mu*((fini(1,2) - fini(0,2))*dz2/dz1 - (fini(2,2) - fini(0,2))*dz1/dz2)/(dz2-dz1);
//                  // Point 1 in fluid and 2 in solid
//                  } else if ((level_set(2,0) > threshold) && (level_set(2,1) <= threshold)) {
// 	            double dz1 = (point(1,2)-point(0,2));
//                     dfdz = mu*(fini(1,2) - fini(0,2))/dz1;
//                  // Point 1 is present in solid
//                  } else if (level_set(2,0) <= threshold) {
//                     double dz1 = (point(1,2)-point(0,2));
//                     dfdz = mu*(fini(1,2) - fini(0,2))/dz1;
//                  }
//               // Point 1 in computational domain
// 	      } else if (point_in_domain(0,2) && !point_in_domain(1,2)) {
//                  double dz1 = (point(1,2)-point(0,2));
//                  dfdz = mu*(fini(1,2) - fini(0,2))/dz1;
//               // Particle close to wall
//               } else if (!point_in_domain(0,2)) {
//                  i0(1,2) = (sign(2) == 1) ? (i0(0,2) + 1*sign(2)) : (i0(0,2) + 0*sign(2));
//                  point(1,2) = UF->get_DOF_coordinate(i0(1,2), comp, 2);
//                 fini(1,2) = third_order_ghost_field_estimate(UF, comp, point(0,0), point(0,1), point(1,2), i0(0,0), i0(0,1), i0(1,2), 2, sign, 0);
//                  double dz1 = (point(1,2)-point(0,2));
//                  dfdz = mu*(fini(1,2) - fini(0,2))/dz1;
// 	      }
// 	   }
//
//            if (comp == 0) {
//               stress(i,0) = 2.*dfdx;
//               stress(i,3) = stress(i,3) + dfdy;
//               stress(i,5) = stress(i,5) + dfdz;
//            } else if (comp == 1) {
//               stress(i,1) = 2.*dfdy;
//               stress(i,3) = stress(i,3) + dfdx;
//               stress(i,4) = stress(i,4) + dfdz;
//            } else if (comp == 2) {
//               stress(i,2) = 2.*dfdz;
//               stress(i,4) = stress(i,4) + dfdy;
//               stress(i,5) = stress(i,5) + dfdx;
//            }
// 	}
//      }
//
//      double scale = (dim == 2) ? ri : ri*ri;
//
//      // Ref: Keating thesis Pg-85
//      // point_coord*(area) --> Component of area in particular direction
//      double fx = stress(i,0)*rotated_normal(0)*(surface.area->item(i)*scale)
//                + stress(i,3)*rotated_normal(1)*(surface.area->item(i)*scale)
//                + stress(i,5)*rotated_normal(2)*(surface.area->item(i)*scale);
//      double fy = stress(i,3)*rotated_normal(0)*(surface.area->item(i)*scale)
//                + stress(i,1)*rotated_normal(1)*(surface.area->item(i)*scale)
//                + stress(i,4)*rotated_normal(2)*(surface.area->item(i)*scale);
//      double fz = stress(i,5)*rotated_normal(0)*(surface.area->item(i)*scale)
//                + stress(i,4)*rotated_normal(1)*(surface.area->item(i)*scale)
//                + stress(i,2)*rotated_normal(2)*(surface.area->item(i)*scale);
//
//      surface_point(i,0) = point(0,0);
//      surface_point(i,1) = point(0,1);
//      surface_point(i,2) = point(0,2);
//
//      surface_force(i,0) = fx;
//      surface_force(i,1) = fy;
//      surface_force(i,2) = fz;
//
//      force(0) = force(0) + fx ;
//      force(1) = force(1) + fy ;
//      force(2) = force(2) + fz ;
//
//      torque(0) = torque(0) + fz*rotated_coord(1) - fy*rotated_coord(2);
//      torque(1) = torque(1) + fx*rotated_coord(2) - fz*rotated_coord(0);
//      torque(2) = torque(2) + fy*rotated_coord(0) - fx*rotated_coord(1);
//   }
//
//   write_surface_discretized_forces(UF,Np,parID,surface_point,surface_force,t_it);
//
//   hydro_forces.vel[0]->set_item(parID,force(0)) ;
//   hydro_forces.vel[1]->set_item(parID,force(1)) ;
//   hydro_forces.vel[2]->set_item(parID,force(2)) ;
//
//   hydro_torque.vel[0]->set_item(parID,torque(0)) ;
//   hydro_torque.vel[1]->set_item(parID,torque(1)) ;
//   hydro_torque.vel[2]->set_item(parID,torque(2)) ;
// }
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: first_order_viscous_stress(size_t const& parID, size_t const& Np, FV_TimeIterator const* t_it )
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: first_order_viscous_stress" ) ;
//
//   size_t i0_temp;
//   double dfdx=0.,dfdy=0., dfdz=0., dzh=0.;
//
//   // Structure of particle input data
//   PartInput solid = GLOBAL_EQ->get_solid(0);
//   // Structure of particle surface input data
//   SurfaceDiscretize surface = GLOBAL_EQ->get_surface();
//
//   // comp won't matter as the particle position is independent of comp
//   double xp = solid.coord[0]->item(parID);
//   double yp = solid.coord[1]->item(parID);
//   double zp = solid.coord[2]->item(parID);
//   double ri = solid.size->item(parID);
//
//   doubleVector xpoint(3,0);
//   doubleVector ypoint(3,0);
//   doubleVector zpoint(3,0);
//   doubleVector finx(3,0);
//   doubleVector finy(3,0);
//   doubleVector finz(3,0);
//   doubleArray2D stress(Np,6,0);         //xx,yy,zz,xy,yz,zx
//   doubleArray2D level_set(dim,2,1.);
//   boolArray2D in_domain(dim,2,true);        //true if ghost point in the computational domain
//   size_t_array2D in_parID(dim,2,0);         //Store particle ID if level_set becomes negative
//   boolArray2D found(dim,3,false);
//   size_t_vector i0_x(3,0);
//   size_t_vector i0_y(3,0);
//   size_t_vector i0_z(3,0);
//   vector<double> net_vel(3,0.);
//   doubleArray2D surface_force(Np,3,0);
//   doubleArray2D surface_point(Np,3,0);
//
//   size_t_vector min_unknown_index(dim,0);
//   size_t_vector max_unknown_index(dim,0);
//
//   doubleVector Dmin(dim,0);
//   doubleVector Dmax(dim,0);
//   doubleVector rotated_coord(3,0);
//   doubleVector rotated_normal(3,0);
//
//   // Structure of particle input data
//   PartForces hydro_forces = GLOBAL_EQ->get_forces(0);
//   PartForces hydro_torque = GLOBAL_EQ->get_torque(0);
//
//   // Force vector
//   doubleVector force(3,0);
//   doubleVector torque(3,0);
//
//   for (size_t i=0;i<Np;i++) {
//      for (size_t comp=0;comp<nb_comps[1];comp++) {
//         // Get local min and max indices
//         // One extra grid cell needs to considered, since ghost points can be
//         // located in between the min/max index handled by the proc
//         for (size_t l=0;l<dim;++l) {
//            min_unknown_index(l) = UF->get_min_index_unknown_handled_by_proc( comp, l );
//            max_unknown_index(l) = UF->get_max_index_unknown_handled_by_proc( comp, l );
//            if (rank_in_i[l] == 0) {
//               Dmin(l) = UF->get_DOF_coordinate( min_unknown_index(l), comp, l ) - UF->get_cell_size(min_unknown_index(l),comp,l);
//               Dmax(l) = UF->get_DOF_coordinate( max_unknown_index(l), comp, l ) + UF->get_cell_size(max_unknown_index(l),comp,l);
//            } else  {
//               Dmin(l) = UF->get_DOF_coordinate( min_unknown_index(l), comp, l );
//               Dmax(l) = UF->get_DOF_coordinate( max_unknown_index(l), comp, l ) + UF->get_cell_size(max_unknown_index(l),comp,l);
//            }
//         }
//
// 	// Rotating surface points
// 	rotated_coord(0) = ri*surface.coordinate[0]->item(i);
// 	rotated_coord(1) = ri*surface.coordinate[1]->item(i);
// 	rotated_coord(2) = ri*surface.coordinate[2]->item(i);
//
//         rotation_matrix(parID,rotated_coord,comp,1);
//
//         xpoint(0) = xp + rotated_coord(0);
//         ypoint(0) = yp + rotated_coord(1);
//         zpoint(0) = zp + rotated_coord(2);
//
//         // Rotating surface normal
// 	rotated_normal(0) = surface.normal[0]->item(i);
// 	rotated_normal(1) = surface.normal[1]->item(i);
// 	rotated_normal(2) = surface.normal[2]->item(i);
//
// 	rotation_matrix(parID,rotated_normal,comp,1);
//
//         if (is_periodic[1][0]) {
//            double isize = UF->primary_grid()->get_main_domain_max_coordinate(0) - UF->primary_grid()->get_main_domain_min_coordinate(0);
//            double imin = UF->primary_grid()->get_main_domain_min_coordinate(0);
//            xpoint(0) = xpoint(0) - MAC::floor((xpoint(0)-imin)/isize)*isize;
//         }
//         if (is_periodic[1][1]) {
//            double isize = UF->primary_grid()->get_main_domain_max_coordinate(1) - UF->primary_grid()->get_main_domain_min_coordinate(1);
//            double imin = UF->primary_grid()->get_main_domain_min_coordinate(1);
//            ypoint(0) = ypoint(0) - MAC::floor((ypoint(0)-imin)/isize)*isize;
//         }
//         if (is_periodic[1][2]) {
//            double isize = UF->primary_grid()->get_main_domain_max_coordinate(2) - UF->primary_grid()->get_main_domain_min_coordinate(2);
//            double imin = UF->primary_grid()->get_main_domain_min_coordinate(2);
//            zpoint(0) = zpoint(0) - MAC::floor((zpoint(0)-imin)/isize)*isize;
//         }
//         // Displacement correction in case of periodic boundary condition in any or all directions
//
//         // Finding the grid indexes next to ghost points
//         found(0,0) = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,0), xpoint(0), i0_temp);
//         if (found(0,0) == 1) i0_x(0) = i0_temp;
//
//         found(1,0) = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,1), ypoint(0), i0_temp);
//         if (found(1,0) == 1) i0_y(0) = i0_temp;
//
//         double dxh = UF->get_cell_size(i0_x(0),comp,0) ;
//         double dyh = UF->get_cell_size(i0_y(0),comp,1) ;
//
//         if (dim == 3) {
//            found(2,0) = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,2), zpoint(0), i0_temp);
//            if (found(2,0) == 1) i0_z(0) = i0_temp;
//            dzh = UF->get_cell_size(i0_z(0),comp,2) ;
//         }
//
//         double dh = (dim == 2) ? (dxh+dyh)/2. : (dxh+dyh+dzh)/3.;
//
//         bool status = (dim==2) ? ((xpoint(0) > Dmin(0)) && (xpoint(0) <= Dmax(0)) && (ypoint(0) > Dmin(1)) && (ypoint(0) <= Dmax(1))) :
//                                  ((xpoint(0) > Dmin(0)) && (xpoint(0) <= Dmax(0)) && (ypoint(0) > Dmin(1)) && (ypoint(0) <= Dmax(1))
//                                                                                   && (zpoint(0) > Dmin(2)) && (zpoint(0) <= Dmax(2)));
//         double threshold = pow(loc_thres,0.5)*dh;
//
//         if (status) {
//            double sign_x = (rotated_normal(0) > 0.) ? 1. : -1.;
//            double sign_y = (rotated_normal(1) > 0.) ? 1. : -1.;
//            double sign_z = 1.;
//
//            // Ghost points in x for the calculation of x-derivative of field
//            xpoint(1) = xpoint(0) + sign_x*dh;
//            xpoint(2) = xpoint(1) + sign_x*dh;
//
//            if (is_periodic[1][0]) {
//               double isize = UF->primary_grid()->get_main_domain_max_coordinate(0) - UF->primary_grid()->get_main_domain_min_coordinate(0);
//               double imin = UF->primary_grid()->get_main_domain_min_coordinate(0);
//               xpoint(1) = xpoint(1) - MAC::floor((xpoint(1)-imin)/isize)*isize;
//               xpoint(2) = xpoint(2) - MAC::floor((xpoint(2)-imin)/isize)*isize;
//            }
//
//            // Ghost points in y for the calculation of y-derivative of field
//            ypoint(1) = ypoint(0) + sign_y*dh;
//            ypoint(2) = ypoint(1) + sign_y*dh;
//
//            if (is_periodic[1][1]) {
//               double isize = UF->primary_grid()->get_main_domain_max_coordinate(1) - UF->primary_grid()->get_main_domain_min_coordinate(1);
//               double imin = UF->primary_grid()->get_main_domain_min_coordinate(1);
//               ypoint(1) = ypoint(1) - MAC::floor((ypoint(1)-imin)/isize)*isize;
//               ypoint(2) = ypoint(2) - MAC::floor((ypoint(2)-imin)/isize)*isize;
//            }
//
//            if (dim == 3) {
//               sign_z = (rotated_normal(2) > 0.) ? 1. : -1.;
//               // Ghost points in z for the calculation of z-derivative of field
//               zpoint(1) = zpoint(0) + sign_z*dh;
//               zpoint(2) = zpoint(1) + sign_z*dh;
//
//               if (is_periodic[1][2]) {
//                  double isize = UF->primary_grid()->get_main_domain_max_coordinate(2) - UF->primary_grid()->get_main_domain_min_coordinate(2);
//                  double imin = UF->primary_grid()->get_main_domain_min_coordinate(2);
//                  zpoint(1) = zpoint(1) - MAC::floor((zpoint(1)-imin)/isize)*isize;
//                  zpoint(2) = zpoint(2) - MAC::floor((zpoint(2)-imin)/isize)*isize;
//               }
//            }
//
//            // Assuming all ghost points are in fluid
//            level_set(0,0) = 1.; level_set(0,1) = 1.;
//            level_set(1,0) = 1.; level_set(1,1) = 1.;
//            if (dim == 3) {level_set(2,0) = 1.; level_set(2,1) = 1.;}
//
//            // Checking all the ghost points in the solid/fluid, and storing the parID if present in solid
//            for (size_t m=0;m<Npart;m++) {
//               if (level_set(0,0) > threshold) {
//                  level_set(0,0) = level_set_function(UF,m,comp,xpoint(1),ypoint(0),zpoint(0),level_set_type,1);
//                  level_set(0,0) *= solid.inside->item(m);
//                  if (level_set(0,0) < threshold) in_parID(0,0) = m;
//               }
//               if (level_set(0,1) > threshold) {
//                  level_set(0,1) = level_set_function(UF,m,comp,xpoint(2),ypoint(0),zpoint(0),level_set_type,1);
//                  level_set(0,1) *= solid.inside->item(m);
//                  if (level_set(0,1) < threshold) in_parID(0,1) = m;
//               }
//               if (level_set(1,0) > threshold) {
//                  level_set(1,0) = level_set_function(UF,m,comp,xpoint(0),ypoint(1),zpoint(0),level_set_type,1);
//                  level_set(1,0) *= solid.inside->item(m);
//                  if (level_set(1,0) < threshold) in_parID(1,0) = m;
//               }
//               if (level_set(1,1) > threshold) {
//                  level_set(1,1) = level_set_function(UF,m,comp,xpoint(0),ypoint(2),zpoint(0),level_set_type,1);
//                  level_set(1,1) *= solid.inside->item(m);
//                  if (level_set(1,1) < threshold) in_parID(1,1) = m;
//               }
//               if (dim == 3) {
//                  if (level_set(2,0) > threshold) {
//                     level_set(2,0) = level_set_function(UF,m,comp,xpoint(0),ypoint(0),zpoint(1),level_set_type,1);
//                     level_set(2,0) *= solid.inside->item(m);
//                     if (level_set(2,0) < threshold) in_parID(2,0) = m;
//                  }
//                  if (level_set(2,1) > threshold) {
//                     level_set(2,1) = level_set_function(UF,m,comp,xpoint(0),ypoint(0),zpoint(2),level_set_type,1);
//                     level_set(2,1) *= solid.inside->item(m);
//                     if (level_set(2,1) < threshold) in_parID(2,1) = m;
//                  }
//               }
//            }
//
//            // Finding the grid indexes next to ghost points
//            found(0,1) = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,0), xpoint(1), i0_temp);
//            if (found(0,1) == 1) i0_x(1) = i0_temp;
//
//            found(0,2) = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,0), xpoint(2), i0_temp);
//            if (found(0,2) == 1) i0_x(2) = i0_temp;
//
//            found(1,1) = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,1), ypoint(1), i0_temp);
//            if (found(1,1) == 1) i0_y(1) = i0_temp;
//
//            found(1,2) = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,1), ypoint(2), i0_temp);
//            if (found(1,2) == 1) i0_y(2) = i0_temp;
//
//            // Calculation of field variable on ghost point(0,0)
//            impose_solid_velocity_for_ghost(net_vel,comp,xpoint(0),ypoint(0),zpoint(0),parID);
//            finx(0) = net_vel[comp];
//            finy(0) = net_vel[comp];
//            finz(0) = net_vel[comp];
//
//            if (dim == 2) {
//               in_domain(0,0) = found(0,1) && found(1,0);
//               in_domain(0,1) = found(0,2) && found(1,0);
//               in_domain(1,0) = found(0,0) && found(1,1);
//               in_domain(1,1) = found(0,0) && found(1,2);
// 	      // Calculation of field variable on ghost point(1,0)
//               if ((level_set(0,0) > threshold) && in_domain(0,0)) {
//                   finx(1) = ghost_field_estimate_on_face (UF,comp,i0_x(1),i0_y(0),0, xpoint(1), ypoint(0),0, dh,2,0);
// 	      } else if ((level_set(0,0) <= threshold) && in_domain(0,0)) {
// 		  finx(1) = net_vel[comp];
// 	      }
//               // Calculation of field variable on ghost point(2,0)
//               if ((level_set(0,1) > threshold) && in_domain(0,1)) {
//                   finx(2) = ghost_field_estimate_on_face (UF,comp,i0_x(2),i0_y(0),0, xpoint(2), ypoint(0),0, dh,2,0);
// 	      } else if ((level_set(0,1) <= threshold) && in_domain(0,1)) {
//                   finx(2) = net_vel[comp];
// 	      }
//               // Calculation of field variable on ghost point(0,1)
//               if ((level_set(1,0) > threshold) && in_domain(1,0)) {
//                   finy(1) = ghost_field_estimate_on_face (UF,comp,i0_x(0),i0_y(1),0, xpoint(0), ypoint(1),0, dh,2,0);
// 	      } else if ((level_set(1,0) <= threshold) && in_domain(1,0)) {
// 		  finy(1) = net_vel[comp];
// 	      }
//               // Calculation of field variable on ghost point(0,2)
//               if ((level_set(1,1) > threshold) && in_domain(1,1)) {
//                   finy(2) = ghost_field_estimate_on_face (UF,comp,i0_x(0),i0_y(2),0, xpoint(0), ypoint(2),0, dh,2,0);
// 	      } else if ((level_set(1,1) <= threshold) && in_domain(1,1)) {
// 		  finy(2) = net_vel[comp];
// 	      }
//            } else if (dim == 3) {
//               found(2,1) = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,2), zpoint(1), i0_temp);
//               if (found(2,1) == 1) i0_z(1) = i0_temp;
//
//               found(2,2) = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,2), zpoint(2), i0_temp);
//               if (found(2,2) == 1) i0_z(2) = i0_temp;
//
//               in_domain(0,0) = found(0,1) && found(1,0) && found(2,0);
//               in_domain(0,1) = found(0,2) && found(1,0) && found(2,0);
//               in_domain(1,0) = found(1,1) && found(0,0) && found(2,0);
//               in_domain(1,1) = found(1,2) && found(0,0) && found(2,0);
//               in_domain(2,0) = found(2,1) && found(0,0) && found(1,0);
//               in_domain(2,1) = found(2,2) && found(0,0) && found(1,0);
//
//               // Calculation of field variable on ghost point(1,0,0)
//               if ((level_set(0,0) > threshold) && in_domain(0,0)) {
//                  finx(1) = ghost_field_estimate_in_box (UF,comp,i0_x(1),i0_y(0),i0_z(0),xpoint(1),ypoint(0),zpoint(0),dh,0,parID);
// 	      } else if ((level_set(0,0) <= threshold) && in_domain(0,0)) {
// 	         finx(1) = net_vel[comp];
// 	      }
//               // Calculation of field variable on ghost point(2,0,0)
//               if ((level_set(0,1) > threshold) && in_domain(0,1)) {
//                  finx(2) = ghost_field_estimate_in_box (UF,comp,i0_x(2),i0_y(0),i0_z(0),xpoint(2),ypoint(0),zpoint(0),dh,0,parID);
//               } else if ((level_set(0,1) <= threshold) && in_domain(0,1)) {
// 	         finx(2) = net_vel[comp];
// 	      }
//               // Calculation of field variable on ghost point(0,1,0)
//               if ((level_set(1,0) > threshold) && in_domain(1,0)) {
//                  finy(1) = ghost_field_estimate_in_box (UF,comp,i0_x(0),i0_y(1),i0_z(0),xpoint(0),ypoint(1),zpoint(0),dh,0,parID);
// 	      } else if ((level_set(1,0) <= threshold) && in_domain(1,0)) {
// 		 finy(1) = net_vel[comp];
// 	      }
//               // Calculation of field variable on ghost point(0,2,0)
//               if ((level_set(1,1) > threshold) && in_domain(1,1)) {
//                  finy(2) = ghost_field_estimate_in_box (UF,comp,i0_x(0),i0_y(2),i0_z(0),xpoint(0),ypoint(2),zpoint(0),dh,0,parID);
// 	      } else if ((level_set(1,1) <= threshold) && in_domain(1,1)) {
// 		 finy(2) = net_vel[comp];
// 	      }
//               // Calculation of field variable on ghost point(0,0,1)
//               if ((level_set(2,0) > threshold) && in_domain(2,0)) {
//                  finz(1) = ghost_field_estimate_in_box (UF,comp,i0_x(0),i0_y(0),i0_z(1),xpoint(0),ypoint(0),zpoint(1),dh,0,parID);
//               } else if ((level_set(2,0) <= threshold) && in_domain(2,0)) {
// 		 finz(1) = net_vel[comp];
// 	      }
//               // Calculation of field variable on ghost point(0,0,2)
//               if ((level_set(2,1) > threshold) && in_domain(2,1)) {
//                  finz(2) = ghost_field_estimate_in_box (UF,comp,i0_x(0),i0_y(0),i0_z(2),xpoint(0),ypoint(0),zpoint(2),dh,0,parID);
// 	      } else if ((level_set(2,1) <= threshold) && in_domain(2,1)) {
// 		 finz(2) = net_vel[comp];
// 	      }
//
//               // Derivative in z
//               // Both points 1 and 2 are in fluid, and both in the computational domain
//               if ((level_set(2,0) > threshold) && (level_set(2,1) > threshold) && (in_domain(2,0) && in_domain(2,1))) {
//                  dfdz = mu*(-finz(2) + 4.*finz(1) - 3.*finz(0))/2./dh;
//               // Point 1 in fluid and 2 is either in the solid or out of the computational domain
//               } else if ((level_set(2,0) > threshold) && ((level_set(2,1) <= threshold) || ((in_domain(2,1) == 0) && (in_domain(2,0) == 1)))) {
//                  dfdz = mu*(finz(1) - finz(0))/dh;
//               // Point 1 is present in solid
//               } else if (level_set(2,0) <= threshold) {
//                  impose_solid_velocity_for_ghost(net_vel,comp,xpoint(0),ypoint(0),zpoint(1),in_parID(2,0));
//                  dfdz = mu*(net_vel[comp] - finz(0))/dh;
//               // Point 1 is out of the computational domain
//               } else if (in_domain(2,0) == 0) {
//                  double dh_wall = (sign_z > 0.) ? MAC::abs(zpoint(0)-UF->primary_grid()->get_main_domain_max_coordinate(2)) :
//                                                   MAC::abs(zpoint(0)-UF->primary_grid()->get_main_domain_min_coordinate(2)) ;
//                  size_t ix,iy,iz;
//                  bool found_x = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,0), xpoint(0), i0_temp);
//                  if (found_x == 1) ix = i0_temp;
//                  bool found_y = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,1), ypoint(0), i0_temp);
//                  if (found_y == 1) iy = i0_temp;
//                  bool found_z = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,2), zpoint(0)+sign_z*dh_wall, i0_temp);
//                  if (found_z == 1) iz = i0_temp;
//                  finz(1) = ghost_field_estimate_in_box (UF,comp,ix,iy,iz,xpoint(0),ypoint(0),zpoint(0)+sign_z*dh_wall,dh_wall,0,parID);
//                  dfdz = mu*(finz(1) - finz(0))/dh_wall;
//               }
//               dfdz *= sign_z;
//            }
//
//            // Derivative in x
//            // Both points 1 and 2 are in fluid, and both in the computational domain
//            if ((level_set(0,0) > threshold) && (level_set(0,1) > threshold) && (in_domain(0,0) && in_domain(0,1))) {
//               dfdx = mu*(-finx(2) + 4.*finx(1) - 3.*finx(0))/2./dh;
//            // Point 1 in fluid and 2 is either in the solid or out of the computational domain
//            } else if ((level_set(0,0) > threshold) && ((level_set(0,1) <= threshold) || ((in_domain(0,1) == 0) && (in_domain(0,0) == 1)))) {
//               dfdx = mu*(finx(1) - finx(0))/dh;
//            // Point 1 is present in solid
//            } else if (level_set(0,0) <= threshold) {
//               impose_solid_velocity_for_ghost(net_vel,comp,xpoint(1),ypoint(0),zpoint(0),in_parID(0,0));
//               dfdx = mu*(net_vel[comp] - finx(0))/dh;
//            // Point 1 is out of the computational domain
//            } else if (in_domain(0,0) == 0) {
//               double dh_wall = (sign_x > 0.) ? MAC::abs(xpoint(0)-UF->primary_grid()->get_main_domain_max_coordinate(0)) :
//                                                MAC::abs(xpoint(0)-UF->primary_grid()->get_main_domain_min_coordinate(0)) ;
//               size_t ix,iy,iz=0;
//               bool found_x = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,0), xpoint(0)+sign_x*dh_wall, i0_temp);
//               if (found_x == 1) ix = i0_temp;
//               bool found_y = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,1), ypoint(0), i0_temp);
//               if (found_y == 1) iy = i0_temp;
//               bool found_z = (dim == 3) ? FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,2), zpoint(0), i0_temp) : 0;
//               if (found_z == 1) iz = i0_temp;
//               finx(1) = (dim == 2) ? ghost_field_estimate_on_face (UF,comp,ix,iy,0, xpoint(0)+sign_x*dh_wall, ypoint(0),0, dh_wall,2,0) :
//                                      ghost_field_estimate_in_box (UF,comp,ix,iy,iz, xpoint(0)+sign_x*dh_wall, ypoint(0),zpoint(0),dh_wall,0,parID);
//               dfdx = mu*(finx(1) - finx(0))/dh_wall;
//            }
//            dfdx *= sign_x;
//
//            // Derivative in y
//            // Both points 1 and 2 are in fluid, and both in the computational domain
//            if ((level_set(1,0) > threshold) && (level_set(1,1) > threshold) && (in_domain(1,0) && in_domain(1,1))) {
//               dfdy = mu*(-finy(2) + 4.*finy(1) - 3.*finy(0))/2./dh;
//            // Point 1 in fluid and 2 is either in the solid or out of the computational domain
//            } else if ((level_set(1,0) > threshold) && ((level_set(1,1) <= threshold) || ((in_domain(1,1) == 0) && (in_domain(1,0) == 1)))) {
//               dfdy = mu*(finy(1) - finy(0))/dh;
//            // Point 1 is present in solid
//            } else if (level_set(1,0) <= threshold) {
//               impose_solid_velocity_for_ghost(net_vel,comp,xpoint(0),ypoint(1),zpoint(0),in_parID(1,0));
//               dfdy = mu*(net_vel[comp] - finy(0))/dh;
//            // Point 1 is out of the computational domain
//            } else if (in_domain(1,0) == 0) {
//               double dh_wall = (sign_y > 0.) ? MAC::abs(ypoint(0)-UF->primary_grid()->get_main_domain_max_coordinate(1)) :
//                                                MAC::abs(ypoint(0)-UF->primary_grid()->get_main_domain_min_coordinate(1)) ;
//               size_t ix,iy,iz=0;
//               bool found_x = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,0), xpoint(0), i0_temp);
//               if (found_x == 1) ix = i0_temp;
//               bool found_y = FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,1), ypoint(0)+sign_y*dh_wall, i0_temp);
//               if (found_y == 1) iy = i0_temp;
//               bool found_z = (dim == 3) ? FV_Mesh::between(UF->get_DOF_coordinates_vector(comp,2), zpoint(0), i0_temp) : 0;
//               if (found_z == 1) iz = i0_temp;
//               finy(1) = (dim == 2) ? ghost_field_estimate_on_face (UF,comp,ix,iy,0, xpoint(0), ypoint(0)+sign_y*dh_wall,0, dh_wall,2,0) :
//                                      ghost_field_estimate_in_box (UF,comp,ix,iy,iz, xpoint(0), ypoint(0)+sign_y*dh_wall,zpoint(0),dh_wall,0,parID);
//               dfdy = mu*(finy(1) - finy(0))/dh_wall;
//            }
//            dfdy *= sign_y;
//
//            if (comp == 0) {
//               stress(i,0) = 2.*dfdx;
//               stress(i,3) = stress(i,3) + dfdy;
//               stress(i,5) = stress(i,5) + dfdz;
//            } else if (comp == 1) {
//               stress(i,1) = 2.*dfdy;
//               stress(i,3) = stress(i,3) + dfdx;
//               stress(i,4) = stress(i,4) + dfdz;
//            } else if (comp == 2) {
//               stress(i,2) = 2.*dfdz;
//               stress(i,4) = stress(i,4) + dfdy;
//               stress(i,5) = stress(i,5) + dfdx;
//            }
//
// 	}
//      }
//
//      double scale = (dim == 2) ? ri : ri*ri;
//
//      // Ref: Keating thesis Pg-85
//      // point_coord*(area) --> Component of area in particular direction
//      double fx = stress(i,0)*rotated_normal(0)*(surface.area->item(i)*scale)
//                + stress(i,3)*rotated_normal(1)*(surface.area->item(i)*scale)
//                + stress(i,5)*rotated_normal(2)*(surface.area->item(i)*scale);
//
//      double fy = stress(i,3)*rotated_normal(0)*(surface.area->item(i)*scale)
//                + stress(i,1)*rotated_normal(1)*(surface.area->item(i)*scale)
//                + stress(i,4)*rotated_normal(2)*(surface.area->item(i)*scale);
//
//      double fz = stress(i,5)*rotated_normal(0)*(surface.area->item(i)*scale)
//                + stress(i,4)*rotated_normal(1)*(surface.area->item(i)*scale)
//                + stress(i,2)*rotated_normal(2)*(surface.area->item(i)*scale);
//
//      surface_force(i,0) = fx;
//      surface_force(i,1) = fy;
//      surface_force(i,2) = fz;
//
//      surface_point(i,0) = xpoint(0);
//      surface_point(i,1) = ypoint(0);
//      surface_point(i,2) = zpoint(0);
//
//      force(0) = force(0) + fx ;
//      force(1) = force(1) + fy ;
//      force(2) = force(2) + fz ;
//
//      torque(0) = torque(0) + fz*rotated_coord(1) - fy*rotated_coord(2);
//      torque(1) = torque(1) + fx*rotated_coord(2) - fz*rotated_coord(0);
//      torque(2) = torque(2) + fy*rotated_coord(0) - fx*rotated_coord(1);
//   }
//
//   write_surface_discretized_forces(UF,Np,parID,surface_point,surface_force,t_it);
//
//   hydro_forces.vel[0]->set_item(parID,force(0)) ;
//   hydro_forces.vel[1]->set_item(parID,force(1)) ;
//   hydro_forces.vel[2]->set_item(parID,force(2)) ;
//
//   hydro_torque.vel[0]->set_item(parID,torque(0)) ;
//   hydro_torque.vel[1]->set_item(parID,torque(1)) ;
//   hydro_torque.vel[2]->set_item(parID,torque(2)) ;
// }
//
// //---------------------------------------------------------------------------
// double
// DDS_NavierStokes:: quadratic_interpolation3D ( FV_DiscreteField* FF, size_t const& comp, double const& xp, double const& yp, double const& zp, size_t const& ii, size_t const& ji, size_t const& ki, size_t const& ghost_points_dir, class intVector& sign, size_t const& level)
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: quadratic_interpolation3D" ) ;
//
//   // Structure of particle input data
//   PartInput solid = GLOBAL_EQ->get_solid(0);
//
//   doubleArray2D point(3,3,0);
//   // Directional indexes of ghost points
//   size_t_vector index(3,0);
//   // Directional indexes of secondary ghost points
//   size_t_array2D index_g(3,3,0);
//   // Coordinates of secondary ghost points
//   doubleArray2D coord_g(3,3,0.);
//   // level set value of the secondary ghost points
//   doubleVector level_set(3,1.);
//   // Store particle ID if level_set becomes negative
//   size_t_vector in_parID(3,0);
//   // Presence in solid or not
//   boolVector point_in_solid(3,0);
//   // Presence in domain or not
//   boolVector point_in_domain(3,1);
//   vector<double> net_vel(3,0.);
//   // Decide which scheme to use
//   string scheme = "quadratic";
//
//   size_t sec_ghost_dir = 0;
//   size_t sec_interpol_dir = 0;
//
//   point(0,0) = xp; point(0,1) = yp; point(0,2) = zp;
//   index(0) = ii; index(1) = ji; index(2) = ki;
//
//   // Ghost points generated in y and then quadratic interpolation in z will generate the same
//   // stencil if ghost points are generated in z and the quadratic interpolation done in y
//   if (ghost_points_dir == 0) {
//      sec_ghost_dir = 1;
//      sec_interpol_dir = 2;
//   } else if (ghost_points_dir == 1) {
//      sec_ghost_dir = 0;
//      sec_interpol_dir = 2;
//   } else if (ghost_points_dir == 2) {
//      sec_ghost_dir = 0;
//      sec_interpol_dir = 1;
//   }
//
//   // Generate indexes of secondary ghost points in sec_ghost_dir direction
//   gen_dir_index_of_secondary_ghost_points(FF, index, sign, sec_ghost_dir, index_g, point_in_domain, comp);
//
//   // Assume all secondary ghost points in fluid
//   double x0 = (point_in_domain(0)) ? FF->get_DOF_coordinate(index_g(0,sec_ghost_dir), comp, sec_ghost_dir) : 0. ;
//   double x1 = (point_in_domain(1)) ? FF->get_DOF_coordinate(index_g(1,sec_ghost_dir), comp, sec_ghost_dir) : 0. ;
//   double x2 = (point_in_domain(2)) ? FF->get_DOF_coordinate(index_g(2,sec_ghost_dir), comp, sec_ghost_dir) : 0. ;
//
//   if (sec_ghost_dir == 0) {
//      coord_g(0,0) = x0; coord_g(0,1) = point(0,1); coord_g(0,2) = point(0,2);
//      coord_g(1,0) = x1; coord_g(1,1) = point(0,1); coord_g(1,2) = point(0,2);
//      coord_g(2,0) = x2; coord_g(2,1) = point(0,1); coord_g(2,2) = point(0,2);
//   } else if (sec_ghost_dir == 1) {
//      coord_g(0,0) = point(0,0); coord_g(0,1) = x0; coord_g(0,2) = point(0,2);
//      coord_g(1,0) = point(0,0); coord_g(1,1) = x1; coord_g(1,2) = point(0,2);
//      coord_g(2,0) = point(0,0); coord_g(2,1) = x2; coord_g(2,2) = point(0,2);
//   } else if (sec_ghost_dir == 2) {
//      coord_g(0,0) = point(0,0); coord_g(0,1) = point(0,1); coord_g(0,2) = x0;
//      coord_g(1,0) = point(0,0); coord_g(1,1) = point(0,1); coord_g(1,2) = x1;
//      coord_g(2,0) = point(0,0); coord_g(2,1) = point(0,1); coord_g(2,2) = x2;
//   }
//
//   double dh = FF->primary_grid()->get_smallest_grid_size();
//
//   double threshold = pow(loc_thres,0.5)*dh;
//   // Checking the secondary ghost points in the solid/fluid, and storing the parID if present in solid
//   for (size_t m=0;m<Npart;m++) {
//      // x0
//      if (level_set(0) > threshold) {
//         level_set(0) = level_set_function(FF,m,comp,coord_g(0,0),coord_g(0,1),coord_g(0,2),level_set_type,1);
//         level_set(0) *= solid.inside->item(m);
//         if (level_set(0) < threshold) { in_parID(0) = m; point_in_solid(0) = 1; }
//      }
//      // x1
//      if (level_set(1) > threshold) {
//         level_set(1) = level_set_function(FF,m,comp,coord_g(1,0),coord_g(1,1),coord_g(1,2),level_set_type,1);
//         level_set(1) *= solid.inside->item(m);
//         if (level_set(1) < threshold) { in_parID(1) = m; point_in_solid(1) = 1; }
//      }
//      // x2
//      if (level_set(2) > threshold) {
//         level_set(2) = level_set_function(FF,m,comp,coord_g(2,0),coord_g(2,1),coord_g(2,2),level_set_type,1);
//         level_set(2) *= solid.inside->item(m);
//         if (level_set(2) < threshold) { in_parID(2) = m; point_in_solid(2) = 1; }
//      }
//   }
//
//   double f0 = 0., f1 = 0., f2 = 0., del = 0.;
//
//   // Estimate the field values at the secondary ghost points
//   f0 = (point_in_domain(0)) ? quadratic_interpolation2D(FF,comp,coord_g(0,0),coord_g(0,1),coord_g(0,2),index_g(0,0),index_g(0,1),index_g(0,2),sec_interpol_dir,sign,level) : 0.;
//   f1 = (point_in_domain(1)) ? quadratic_interpolation2D(FF,comp,coord_g(1,0),coord_g(1,1),coord_g(1,2),index_g(1,0),index_g(1,1),index_g(1,2),sec_interpol_dir,sign,level) : 0.;
//   f2 = (point_in_domain(2)) ? quadratic_interpolation2D(FF,comp,coord_g(2,0),coord_g(2,1),coord_g(2,2),index_g(2,0),index_g(2,1),index_g(2,2),sec_interpol_dir,sign,level) : 0.;
//
//   // Ghost points corrections
//   if (point_in_domain(0) && point_in_domain(1) && point_in_domain(2)) {
//      // 0 in solid, rest in fluid
//      if (point_in_solid(0) && !point_in_solid(1) && !point_in_solid(2)) {
// 	if (FF == UF) {
//            if (sec_ghost_dir == 0) {
//               del = find_intersection_for_ghost(FF, x0, x1, yp, zp, in_parID(0), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,x1-del,yp,zp,in_parID(0));
//            } else if (sec_ghost_dir == 1) {
//               del = find_intersection_for_ghost(FF, x0, x1, xp, zp, in_parID(0), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,x1-del,zp,in_parID(0));
//            } else if (sec_ghost_dir == 2) {
//               del = find_intersection_for_ghost(FF, x0, x1, xp, yp, in_parID(0), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,yp,x1-del,in_parID(0));
//            }
//            x0 = x1 - del;
//            f0 = net_vel[comp];
// 	} else {
//            scheme = "linear12";
// 	}
//      // 2 in solid, rest in fluid
//      } else if (!point_in_solid(0) && !point_in_solid(1) && point_in_solid(2)) {
// 	if (FF == UF) {
//            if (sec_ghost_dir == 0) {
//               del = find_intersection_for_ghost(FF, x1, x2, yp, zp, in_parID(2), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,x1+del,yp,zp,in_parID(2));
//            } else if (sec_ghost_dir == 1) {
//               del = find_intersection_for_ghost(FF, x1, x2, xp, zp, in_parID(2), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,x1+del,zp,in_parID(2));
//            } else if (sec_ghost_dir == 2) {
//               del = find_intersection_for_ghost(FF, x1, x2, xp, yp, in_parID(2), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,yp,x1+del,in_parID(2));
//            }
//            x2 = x1 + del;
//            f2 = net_vel[comp];
// 	} else {
//            scheme = "linear01";
// 	}
//      // 0, 2 in solid; 1 in fluid
//      } else if (point_in_solid(0) && !point_in_solid(1) && point_in_solid(2)) {
// 	if (FF == UF) {
//            if (sec_ghost_dir == 0) {
//               del = find_intersection_for_ghost(FF, x0, x1, yp, zp, in_parID(0), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,x1-del,yp,zp,in_parID(0));
//            } else if (sec_ghost_dir == 1) {
//               del = find_intersection_for_ghost(FF, x0, x1, xp, zp, in_parID(0), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,x1-del,zp,in_parID(0));
//            } else if (sec_ghost_dir == 2) {
//               del = find_intersection_for_ghost(FF, x0, x1, xp, yp, in_parID(0), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,yp,x1-del,in_parID(0));
//            }
//            x0 = x1 - del;
//            f0 = net_vel[comp];
//
//            if (sec_ghost_dir == 0) {
//               del = find_intersection_for_ghost(FF, x1, x2, yp, zp, in_parID(2), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,x1+del,yp,zp,in_parID(2));
//            } else if (sec_ghost_dir == 1) {
//               del = find_intersection_for_ghost(FF, x1, x2, xp, zp, in_parID(2), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,x1+del,zp,in_parID(2));
//            } else if (sec_ghost_dir == 2) {
//               del = find_intersection_for_ghost(FF, x1, x2, xp, yp, in_parID(2), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,yp,x1+del,in_parID(2));
//            }
//            x2 = x1 + del;
//            f2 = net_vel[comp];
// 	} else {
// 	   scheme = "linear1";
// 	}
//      // 0, 1 in solid; 2 in fluid
//      } else if (point_in_solid(0) && point_in_solid(1) && !point_in_solid(2)) {
// 	if (FF == UF) {
//            if (sec_ghost_dir == 0) {
//               del = find_intersection_for_ghost(FF, x1, x2, yp, zp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,x2-del,yp,zp,in_parID(1));
//            } else if (sec_ghost_dir == 1) {
//               del = find_intersection_for_ghost(FF, x1, x2, xp, zp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,x2-del,zp,in_parID(1));
//            } else if (sec_ghost_dir == 2) {
//               del = find_intersection_for_ghost(FF, x1, x2, xp, yp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,yp,x2-del,in_parID(1));
//            }
//            x1 = x2 - del;
//            f1 = net_vel[comp];
//            scheme = "linear12";
// 	} else {
// 	   scheme = "linear2";
// 	}
//      // 1, 2 in solid; 0 in fluid
//      } else if (!point_in_solid(0) && point_in_solid(1) && point_in_solid(2)) {
// 	if (FF == UF) {
//            if (sec_ghost_dir == 0) {
//               del = find_intersection_for_ghost(FF, x0, x1, yp, zp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,x0+del,yp,zp,in_parID(1));
//            } else if (sec_ghost_dir == 1) {
//               del = find_intersection_for_ghost(FF, x0, x1, xp, zp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,x0+del,zp,in_parID(1));
//            } else if (sec_ghost_dir == 2) {
//               del = find_intersection_for_ghost(FF, x0, x1, xp, yp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,yp,x0+del,in_parID(1));
//            }
//            x1 = x0 + del;
//            f1 = net_vel[comp];
//            scheme = "linear01";
// 	} else {
// 	   scheme = "linear0";
// 	}
//      }
//   // Point 0 and 1 are in domain, 2 not in domain
//   } else if (point_in_domain(0) && point_in_domain(1) && !point_in_domain(2)) {
//      scheme = "linear01";
//      // 0 in fluid; 1 in solid
//      if (!point_in_solid(0) && point_in_solid(1)) {
// 	if (FF == UF) {
//            if (sec_ghost_dir == 0) {
//               del = find_intersection_for_ghost(FF, x0, x1, yp, zp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,x0+del,yp,zp,in_parID(1));
//            } else if (sec_ghost_dir == 1) {
//               del = find_intersection_for_ghost(FF, x0, x1, xp, zp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,x0+del,zp,in_parID(1));
//            } else if (sec_ghost_dir == 2) {
//               del = find_intersection_for_ghost(FF, x0, x1, xp, yp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,yp,x0+del,in_parID(1));
//            }
//            x1 = x0 + del;
//            f1 = net_vel[comp];
// 	} else {
// 	   scheme = "linear0";
// 	}
//      // 0 in solid, 1 in fluid
//      } else if (point_in_solid(0) && !point_in_solid(1)) {
// 	if (FF == UF) {
//            if (sec_ghost_dir == 0) {
//               del = find_intersection_for_ghost(FF, x0, x1, yp, zp, in_parID(0), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,x1-del,yp,zp,in_parID(0));
//            } else if (sec_ghost_dir == 1) {
//               del = find_intersection_for_ghost(FF, x0, x1, xp, zp, in_parID(0), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,x1-del,zp,in_parID(0));
//            } else if (sec_ghost_dir == 2) {
//               del = find_intersection_for_ghost(FF, x0, x1, xp, yp, in_parID(0), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,yp,x1-del,in_parID(0));
//            }
//            x0 = x1 - del;
//            f0 = net_vel[comp];
// 	} else {
// 	   scheme = "linear1";
// 	}
//      }
//   // Point 1 and 2 are in domain, 0 not in domain
//   } else if (!point_in_domain(0) && point_in_domain(1) && point_in_domain(2)) {
//      scheme = "linear12";
//      // 1 in fluid; 2 in solid
//      if (!point_in_solid(1) && point_in_solid(2)) {
// 	if (FF == UF) {
//            if (sec_ghost_dir == 0) {
//               del = find_intersection_for_ghost(FF, x1, x2, yp, zp, in_parID(2), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,x1+del,yp,zp,in_parID(2));
//            } else if (sec_ghost_dir == 1) {
//               del = find_intersection_for_ghost(FF, x1, x2, xp, zp, in_parID(2), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,x1+del,zp,in_parID(2));
//            } else if (sec_ghost_dir == 2) {
//               del = find_intersection_for_ghost(FF, x1, x2, xp, yp, in_parID(2), comp, sec_ghost_dir, dh, 1, 0, 1);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,yp,x1+del,in_parID(2));
//            }
//            x2 = x1 + del;
//            f2 = net_vel[comp];
// 	} else {
// 	   scheme = "linear1";
// 	}
//      // 1 in solid, 2 in fluid
//      } else if (point_in_solid(1) && !point_in_solid(2)) {
// 	if (FF == UF) {
//            if (sec_ghost_dir == 0) {
//               del = find_intersection_for_ghost(FF, x1, x2, yp, zp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,x2-del,yp,zp,in_parID(1));
//            } else if (sec_ghost_dir == 1) {
//               del = find_intersection_for_ghost(FF, x1, x2, xp, zp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,x2-del,zp,in_parID(1));
//            } else if (sec_ghost_dir == 2) {
//               del = find_intersection_for_ghost(FF, x1, x2, xp, yp, in_parID(1), comp, sec_ghost_dir, dh, 1, 0, 0);
//               impose_solid_velocity_for_ghost(net_vel,comp,xp,yp,x2-del,in_parID(1));
//            }
//            x1 = x2 - del;
//            f1 = net_vel[comp];
// 	} else {
//            scheme = "linear2";
// 	}
//      }
//   }
//
//   double l0 = 0., l1 = 0., l2 = 0.;
//   double result = 0.;
//
//   if (scheme == "quadratic") {
//      l0 = (point(0,sec_ghost_dir) - x1)*(point(0,sec_ghost_dir) - x2)/(x0 - x1)/(x0 - x2);
//      l1 = (point(0,sec_ghost_dir) - x0)*(point(0,sec_ghost_dir) - x2)/(x1 - x0)/(x1 - x2);
//      l2 = (point(0,sec_ghost_dir) - x0)*(point(0,sec_ghost_dir) - x1)/(x2 - x0)/(x2 - x1);
//      result = f0*l0 + f1*l1 + f2*l2;
//   } else if (scheme == "linear01") {
//      l0 = (point(0,sec_ghost_dir) - x1)/(x0 - x1);
//      l1 = (point(0,sec_ghost_dir) - x0)/(x1 - x0);
//      result = f0*l0 + f1*l1;
//   } else if (scheme == "linear12") {
//      l1 = (point(0,sec_ghost_dir) - x2)/(x1 - x2);
//      l2 = (point(0,sec_ghost_dir) - x1)/(x2 - x1);
//      result = f1*l1 + f2*l2;
//   } else if (scheme == "linear0") {
//      result = f0;
//   } else if (scheme == "linear1") {
//      result = f1;
//   } else if (scheme == "linear2") {
//      result = f2;
//   }
//
//   return(result);
// }
//
// //---------------------------------------------------------------------------
// double
// DDS_NavierStokes:: third_order_ghost_field_estimate ( FV_DiscreteField* FF, size_t const& comp, double const& xp, double const& yp, double const& zp, size_t const& ii, size_t const& ji, size_t const& ki, size_t const& ghost_points_dir, class intVector& sign, size_t const& level)
// //---------------------------------------------------------------------------
// {
//    MAC_LABEL("DDS_NavierStokes:: third_order_ghost_field_estimate" ) ;
//
// // Call respective functions based on 2D or 3D domain
//
//    double result = 0.;
//
//    if (dim == 2) {
//       result = (ghost_points_dir == 0) ? quadratic_interpolation2D(FF, comp, xp, yp, zp, ii, ji, ki, 1, sign, level) :
//                                          quadratic_interpolation2D(FF, comp, xp, yp, zp, ii, ji, ki, 0, sign, level) ;
//    } else if (dim == 3) {
//       result = quadratic_interpolation3D(FF, comp, xp, yp, zp, ii, ji, ki, ghost_points_dir, sign, level) ;
//    }
//
//    return(result);
// }
//
// //---------------------------------------------------------------------------
// void
// DDS_NavierStokes:: gen_dir_index_of_secondary_ghost_points (FV_DiscreteField* FF, class size_t_vector& index, class intVector& sign, size_t const& interpol_dir, class size_t_array2D& index_g, class boolVector& point_in_domain, size_t const& comp)
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL("DDS_NavierStokes:: gen_dir_index_of_secondary_ghost_points" ) ;
//
//   intVector i0_temp(3,0);
//
//   for (size_t dir=0;dir<dim;dir++) {
//      index_g(0,dir) = index(dir);
//      index_g(1,dir) = index(dir);
//      index_g(2,dir) = index(dir);
//   }
//
//   if (sign(interpol_dir) > 0.) {
//      i0_temp(0) = int(index(interpol_dir));
//      i0_temp(1) = int(index(interpol_dir)) + 1;
//      i0_temp(2) = int(index(interpol_dir)) + 2;
//   } else if (sign(interpol_dir) <= 0.) {
//      i0_temp(0) = int(index(interpol_dir)) - 1;
//      i0_temp(1) = int(index(interpol_dir));
//      i0_temp(2) = int(index(interpol_dir)) + 1;
//   }
//
//   // Checking the ghost points in domain or not
//   if ((i0_temp(0) < 0) || (i0_temp(0) >= (int)FF->get_local_nb_dof(comp,interpol_dir))) {
//      point_in_domain(0) = 0;
//   } else {
//      point_in_domain(0) = 1;
//   }
//
//   if ((i0_temp(1) < 0) || (i0_temp(1) >= (int)FF->get_local_nb_dof(comp,interpol_dir))) {
//      point_in_domain(1) = 0;
//   } else {
//      point_in_domain(1) = 1;
//   }
//
//   if ((i0_temp(2) < 0) || (i0_temp(2) >= (int)FF->get_local_nb_dof(comp,interpol_dir))) {
//      point_in_domain(2) = 0;
//   } else {
//      point_in_domain(2) = 1;
//   }
//
//   index_g(0,interpol_dir) = i0_temp(0);
//   index_g(1,interpol_dir) = i0_temp(1);
//   index_g(2,interpol_dir) = i0_temp(2);
//
// }
//
//
//
//
// //---------------------------------------------------------------------------
// double
// DDS_NavierStokes:: quadratic_interpolation2D ( FV_DiscreteField* FF, size_t const& comp, double const& xp, double const& yp, double const& zp, size_t const& i0, size_t const& j0, size_t const& k0, size_t const& interpol_dir, class intVector& sign, size_t const& level)
// //---------------------------------------------------------------------------
// {
//    MAC_LABEL("DDS_NavierStokes:: quadratic_interpolation2D" ) ;
//
// // Calculates the field value at the ghost points
// // near the particle boundary using the quadratic interpolation
// // inspired from Johansen 1998;
// // xp,yp,zp are the ghost point coordinated; interpol_dir is the direction
// // in which the additional points will be used for quadratic interpolation
// /*
//    ofstream outputFile ;
//    std::ostringstream os2;
//    os2 << "./DS_results/velocity_drag_extras.csv";
//    std::string filename = os2.str();
//    outputFile.open(filename.c_str(), ios::app);
// */
//    size_t field = (FF == UF) ? 1 : 0;
//    // Node information for field(1)
//    NodeProp node = GLOBAL_EQ->get_node_property(field,0);
//    // Intersection information for field in fluid(0)
//    BoundaryBisec* bf_intersect = GLOBAL_EQ->get_b_intersect(field,0);
//
//    // Directional index of point
//    size_t_vector index(3,0);
//    boolVector point_in_solid(3,0);
//    boolVector point_in_domain(3,1);
//    doubleVector xi(3,0.);
//    // Directional indexes of ghost points
//    size_t_array2D index_g(3,3,0);
//    // Local node index of ghost points
//    size_t_vector node_index(3,0);
//    // Decide which scheme to use
//    string scheme = "quadratic";
//
//    xi(0) = xp; xi(1) = yp; xi(2) = zp;
//    index(0) = i0; index(1) = j0; index(2) = k0;
//
//    double f0=0.,f1=0.,f2=0.;
//    double x0=0.,x1=0.,x2=0.;
//    double l0=0.,l1=0.,l2=0.;
//
//    // Creating ghost points for quadratic interpolation
//    gen_dir_index_of_secondary_ghost_points(FF, index, sign, interpol_dir, index_g, point_in_domain, comp);
//
//    // Assume all the ghost points in fluid
//    // Storing the field values assuming all ghost points in fluid and domain
//    // Check weather the ghost points are in solid or not; TRUE if they are
//    if (point_in_domain(0)) {
//       x0 = FF->get_DOF_coordinate(index_g(0,interpol_dir), comp, interpol_dir);
//       f0 = FF->DOF_value( index_g(0,0), index_g(0,1), index_g(0,2), comp, level );
//       node_index(0) = return_node_index(FF,comp,index_g(0,0),index_g(0,1),index_g(0,2));
//       point_in_solid(0) = node.void_frac[comp]->item(node_index(0));
//    }
//
//    if (point_in_domain(1)) {
//       x1 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir);
//       f1 = FF->DOF_value( index_g(1,0), index_g(1,1), index_g(1,2), comp, level );
//       node_index(1) = return_node_index(FF,comp,index_g(1,0),index_g(1,1),index_g(1,2));
//       point_in_solid(1) = node.void_frac[comp]->item(node_index(1));
//    }
//
//    if (point_in_domain(2)) {
//       x2 = FF->get_DOF_coordinate(index_g(2,interpol_dir), comp, interpol_dir);
//       f2 = FF->DOF_value( index_g(2,0), index_g(2,1), index_g(2,2), comp, level );
//       node_index(2) = return_node_index(FF,comp,index_g(2,0),index_g(2,1),index_g(2,2));
//       point_in_solid(2) = node.void_frac[comp]->item(node_index(2));
//    }
//
//    // Ghost points corrections
//    // All points in domain
//    if (point_in_domain(0) && point_in_domain(1) && point_in_domain(2)) {
//       // 0 in solid, rest in fluid
//       if (point_in_solid(0) && !point_in_solid(1) && !point_in_solid(2)) {
// 	 if (FF == UF) {
//             x0 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) - bf_intersect[interpol_dir].value[comp]->item(node_index(1),0);
//             f0 = bf_intersect[interpol_dir].field_var[comp]->item(node_index(1),0);
// 	 } else {
//             scheme = "linear12";
// 	 }
//       // 2 in solid, rest in fluid
//       } else if (!point_in_solid(0) && !point_in_solid(1) && point_in_solid(2)) {
// 	 if (FF == UF) {
//             x2 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) + bf_intersect[interpol_dir].value[comp]->item(node_index(1),1);
//             f2 = bf_intersect[interpol_dir].field_var[comp]->item(node_index(1),1);
// 	 } else {
//             scheme = "linear01";
// 	 }
//       // 0, 2 in solid; 1 in fluid
//       } else if (point_in_solid(0) && !point_in_solid(1) && point_in_solid(2)) {
// 	 if (FF == UF) {
//             x0 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) - bf_intersect[interpol_dir].value[comp]->item(node_index(1),0);
//             f0 = bf_intersect[interpol_dir].field_var[comp]->item(node_index(1),0);
//             x2 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) + bf_intersect[interpol_dir].value[comp]->item(node_index(1),1);
//             f2 = bf_intersect[interpol_dir].field_var[comp]->item(node_index(1),1);
// 	 } else {
//             scheme = "linear1";
// 	 }
//       // 0, 1 in solid; 2 in fluid
//       } else if (point_in_solid(0) && point_in_solid(1) && !point_in_solid(2)) {
//          if (FF == UF) {
//             x1 = FF->get_DOF_coordinate(index_g(2,interpol_dir), comp, interpol_dir) - bf_intersect[interpol_dir].value[comp]->item(node_index(2),0);
//             f1 = bf_intersect[interpol_dir].field_var[comp]->item(node_index(2),0);
//             scheme = "linear12";
// 	 } else {
//             scheme = "linear2";
// 	 }
//       // 1, 2 in solid; 0 in fluid
//       } else if (!point_in_solid(0) && point_in_solid(1) && point_in_solid(2)) {
// 	 if (FF == UF) {
//             x1 = FF->get_DOF_coordinate(index_g(0,interpol_dir), comp, interpol_dir) + bf_intersect[interpol_dir].value[comp]->item(node_index(0),1);
//             f1 = bf_intersect[interpol_dir].field_var[comp]->item(node_index(0),1);
//             scheme = "linear01";
// 	 } else {
//             scheme = "linear0";
// 	 }
//       }
//    // Point 0 and 1 are in domain, 2 not in domain
//    } else if (point_in_domain(0) && point_in_domain(1) && !point_in_domain(2)) {
//       scheme = "linear01";
//       // 0 in fluid; 1 in solid
//       if (!point_in_solid(0) && point_in_solid(1)) {
// 	 if (FF == UF) {
//             x1 = FF->get_DOF_coordinate(index_g(0,interpol_dir), comp, interpol_dir) + bf_intersect[interpol_dir].value[comp]->item(node_index(0),1);
//             f1 = bf_intersect[interpol_dir].field_var[comp]->item(node_index(0),1);
// 	 } else {
//             scheme = "linear0";
// 	 }
//       // 0 in solid, 1 in fluid
//       } else if (point_in_solid(0) && !point_in_solid(1)) {
// 	 if (FF == UF) {
//             x0 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) - bf_intersect[interpol_dir].value[comp]->item(node_index(1),0);
//             f0 = bf_intersect[interpol_dir].field_var[comp]->item(node_index(1),0);
// 	 } else {
//             scheme = "linear1";
// 	 }
//       }
//    // Point 1 and 2 are in domain, 0 not in domain
//    } else if (!point_in_domain(0) && point_in_domain(1) && point_in_domain(2)) {
//       scheme = "linear12";
//       // 1 in fluid; 2 in solid
//       if (!point_in_solid(1) && point_in_solid(2)) {
// 	 if (FF == UF) {
//             x2 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) + bf_intersect[interpol_dir].value[comp]->item(node_index(1),1);
//             f2 = bf_intersect[interpol_dir].field_var[comp]->item(node_index(1),1);
// 	 } else {
//             scheme = "linear1";
// 	 }
//       // 1 in solid, 2 in fluid
//       } else if (point_in_solid(1) && !point_in_solid(2)) {
// 	 if (FF == UF) {
//             x1 = FF->get_DOF_coordinate(index_g(2,interpol_dir), comp, interpol_dir) - bf_intersect[interpol_dir].value[comp]->item(node_index(2),0);
//             f1 = bf_intersect[interpol_dir].field_var[comp]->item(node_index(2),0);
// 	 } else {
//             scheme = "linear2";
// 	 }
//       }
//    }
// /*
//    if (comp == 0) {
//       if (interpol_dir == 0) {
//          double yt = FF->get_DOF_coordinate(index(1), comp, 1);
//          outputFile << x0 << "," << yt << "," << 0. << "," << f0 << endl;
//          outputFile << x1 << "," << yt << "," << 0. << "," << f1 << endl;
//          outputFile << x2 << "," << yt << "," << 0. << "," << f2 << endl;
//       } else if (interpol_dir == 1) {
//          double xt = FF->get_DOF_coordinate(index(0), comp, 0);
//          outputFile << xt << "," << x0 << "," << 0. << "," << f0 << endl;
//          outputFile << xt << "," << x1 << "," << 0. << "," << f1 << endl;
//          outputFile << xt << "," << x2 << "," << 0. << "," << f2 << endl;
//       }
//    }
// */
//    double result = 0.;
//
//    if (scheme == "quadratic") {
//       l0 = (xi(interpol_dir) - x1)*(xi(interpol_dir) - x2)/(x0 - x1)/(x0 - x2);
//       l1 = (xi(interpol_dir) - x0)*(xi(interpol_dir) - x2)/(x1 - x0)/(x1 - x2);
//       l2 = (xi(interpol_dir) - x0)*(xi(interpol_dir) - x1)/(x2 - x0)/(x2 - x1);
//       result = f0*l0 + f1*l1 + f2*l2;
//    } else if (scheme == "linear01") {
//       l0 = (xi(interpol_dir) - x1)/(x0 - x1);
//       l1 = (xi(interpol_dir) - x0)/(x1 - x0);
//       result = f0*l0 + f1*l1;
//    } else if (scheme == "linear12") {
//       l1 = (xi(interpol_dir) - x2)/(x1 - x2);
//       l2 = (xi(interpol_dir) - x1)/(x2 - x1);
//       result = f1*l1 + f2*l2;
//    } else if (scheme == "linear0") {
//       result = f0;
//    } else if (scheme == "linear1") {
//       result = f1;
//    } else if (scheme == "linear2") {
//       result = f2;
//    }
//
//    return(result);
// }
//
// //---------------------------------------------------------------------------
// double
// DDS_NavierStokes:: ghost_field_estimate_in_box ( FV_DiscreteField* FF, size_t const& comp, size_t const& i0, size_t const& j0, size_t const& k0, double const& x0, double const& y0, double const& z0, double const& dh, size_t const& level, size_t const& parID)
// //---------------------------------------------------------------------------
// {
//    MAC_LABEL("DDS_NavierStokes:: ghost_field_estimate_in_box" ) ;
//
// // Calculates the field value at the ghost points in the box
// // near the particle boundary considering boundary affects;
// // x0,y0,z0 are the ghost point coordinated; i0,j0,k0 is the
// // bottom left index of the grid coordinate
//
//    doubleArray2D vel(dim,2,0);
//    doubleArray2D del(dim,2,0);
//    double temp = 0.;
//
//    if (FF == UF) {
//       vector<double> net_vel(3,0.);
//
//       // Behind face
//       temp = FF->get_DOF_coordinate(k0,comp , 2);
//       double face_solid = level_set_function (FF,parID,comp,x0,y0,z0,level_set_type,1)*
//                           level_set_function (FF,parID,comp,x0,y0,temp,level_set_type,1);
//
//       if (face_solid > 0) {
//          vel(2,0) = ghost_field_estimate_on_face (FF,comp,i0,j0,k0,x0,y0,temp,dh,2,0);
//          del(2,0) = MAC::abs(temp - z0);
//       } else {
//          del(2,0) = find_intersection_for_ghost(FF, temp, z0, x0, y0, parID, comp, 2, dh, 1, 0, 0);
//          impose_solid_velocity_for_ghost(net_vel,comp,x0,y0,z0-del(2,0),parID);
//          vel(2,0) = net_vel[comp];
//       }
//
//       // Front face
//       temp = FF->get_DOF_coordinate(k0+1,comp , 2);
//       face_solid = level_set_function (FF,parID,comp,x0,y0,z0,level_set_type,1)*
//                    level_set_function (FF,parID,comp,x0,y0,temp,level_set_type,1);
//
//       if (face_solid > 0) {
//          vel(2,1) = ghost_field_estimate_on_face (FF,comp,i0,j0,k0+1,x0,y0,temp,dh,2,0);
//          del(2,1) = MAC::abs(temp - z0);
//       } else {
//          del(2,1) = find_intersection_for_ghost(FF, z0, temp, x0, y0, parID, comp, 2, dh, 1, 0, 1);
//          impose_solid_velocity_for_ghost(net_vel,comp,x0,y0,z0+del(2,1),parID);
//          vel(2,1) = net_vel[comp];
//       }
//
//       // Left face
//       temp = FF->get_DOF_coordinate(i0,comp, 0);
//       face_solid = level_set_function (FF,parID,comp,x0,y0,z0,level_set_type,1)*
//                    level_set_function (FF,parID,comp,temp,y0,z0,level_set_type,1);
//
//       if (face_solid > 0) {
//          vel(0,0) = ghost_field_estimate_on_face (FF,comp,i0,j0,k0,temp,y0,z0,dh,0,0);
//          del(0,0) = MAC::abs(temp - x0);
//       } else {
//          del(0,0) = find_intersection_for_ghost(FF, temp, x0, y0, z0, parID, comp, 0, dh, 1, 0, 0);
//          impose_solid_velocity_for_ghost(net_vel,comp,x0-del(0,0),y0,z0,parID);
//          vel(0,0) = net_vel[comp];
//       }
//
//       // Right face
//       temp = FF->get_DOF_coordinate(i0+1,comp, 0);
//       face_solid = level_set_function (FF,parID,comp,x0,y0,z0,level_set_type,1)*
//                    level_set_function (FF,parID,comp,temp,y0,z0,level_set_type,1);
//
//       if (face_solid > 0) {
//          vel(0,1) = ghost_field_estimate_on_face (FF,comp,i0+1,j0,k0,temp,y0,z0,dh,0,0);
//          del(0,1) = MAC::abs(temp - x0);
//       } else {
//          del(0,1) = find_intersection_for_ghost(FF, x0, temp, y0, z0, parID, comp, 0, dh, 1, 0, 1);
//          impose_solid_velocity_for_ghost(net_vel,comp,x0+del(0,1),y0,z0,parID);
//          vel(0,1) = net_vel[comp];
//       }
//
//       // Bottom face
//       temp = FF->get_DOF_coordinate(j0,comp, 1);
//       face_solid = level_set_function (FF,parID,comp,x0,y0,z0,level_set_type,1)*
//                    level_set_function (FF,parID,comp,x0,temp,z0,level_set_type,1);
//
//       if (face_solid > 0) {
//          vel(1,0) = ghost_field_estimate_on_face (FF,comp,i0,j0,k0,x0,temp,z0,dh,1,0);
//          del(1,0) = MAC::abs(temp - y0);
//       } else {
//          del(1,0) = find_intersection_for_ghost(FF, temp, y0, x0, z0, parID, comp, 1, dh, 1, 0, 0);
//          impose_solid_velocity_for_ghost(net_vel,comp,x0,y0-del(1,0),z0,parID);
//          vel(1,0) = net_vel[comp];
//       }
//
//       // Top face
//       temp = FF->get_DOF_coordinate(j0+1,comp, 1);
//       face_solid = level_set_function (FF,parID,comp,x0,y0,z0,level_set_type,1)*
//                    level_set_function (FF,parID,comp,x0,temp,z0,level_set_type,1);
//
//       if (face_solid > 0) {
//          vel(1,1) = ghost_field_estimate_on_face (FF,comp,i0,j0+1,k0,x0,temp,z0,dh,1,0);
//          del(1,1) = MAC::abs(temp - y0);
//       } else {
//          del(1,1) = find_intersection_for_ghost(FF, y0, temp, x0, z0, parID, comp, 1, dh, 1, 0, 1);
//          impose_solid_velocity_for_ghost(net_vel,comp,x0,y0+del(1,1),z0,parID);
//          vel(1,1) = net_vel[comp];
//       }
//    } else if (FF == PF) {
//       // Behind face
//       temp = FF->get_DOF_coordinate(k0, comp, 2);
//       vel(2,0) = ghost_field_estimate_on_face (FF,comp,i0,j0,k0,x0,y0,temp,0.,2,level);
//       del(2,0) = MAC::abs(temp - z0);
//       // Front face
//       temp = FF->get_DOF_coordinate(k0+1, comp, 2);
//       vel(2,1) = ghost_field_estimate_on_face (FF,comp,i0,j0,k0+1,x0,y0,temp,0.,2,level);
//       del(2,1) = MAC::abs(temp - z0);
//       // Left face
//       temp = FF->get_DOF_coordinate(i0, comp, 0);
//       vel(0,0) = ghost_field_estimate_on_face (FF,comp,i0,j0,k0,temp,y0,z0,0.,0,level);
//       del(0,0) = MAC::abs(temp - x0);
//       // Right face
//       temp = FF->get_DOF_coordinate(i0+1, comp, 0);
//       vel(0,1) = ghost_field_estimate_on_face (FF,comp,i0+1,j0,k0,temp,y0,z0,0.,0,level);
//       del(0,1) = MAC::abs(temp - x0);
//       // Bottom face
//       temp = FF->get_DOF_coordinate(j0, comp, 1);
//       vel(1,0) = ghost_field_estimate_on_face (FF,comp,i0,j0,k0,x0,temp,z0,0.,1,level);
//       del(1,0) = MAC::abs(temp - y0);
//       // Bottom face
//       temp = FF->get_DOF_coordinate(j0+1, comp, 1);
//       vel(1,1) = ghost_field_estimate_on_face (FF,comp,i0,j0+1,k0,x0,temp,z0,0.,1,level);
//       del(1,1) = MAC::abs(temp - y0);
//    }
//
//    double value = (1./3.)*((vel(0,1)*del(0,0)+vel(0,0)*del(0,1))/(del(0,0)+del(0,1)) +
//                            (vel(1,1)*del(1,0)+vel(1,0)*del(1,1))/(del(1,0)+del(1,1)) +
//                            (vel(2,1)*del(2,0)+vel(2,0)*del(2,1))/(del(2,0)+del(2,1)));
//
//    return(value);
// }
//
// //---------------------------------------------------------------------------
// double
// DDS_NavierStokes:: ghost_field_estimate_on_face ( FV_DiscreteField* FF, size_t const& comp, size_t const& i0, size_t const& j0, size_t const& k0, double const& x0, double const& y0, double const& z0, double const& dh, size_t const& face_vec, size_t const& level)
// //---------------------------------------------------------------------------
// {
//    MAC_LABEL("DDS_NavierStokes:: ghost_field_estimate_on_face" ) ;
//
// // Calculates the field value on a face at the ghost points
// // near the particle boundary considering boundary affects;
// // x0,y0,z0 are the ghost point coordinated; i0,j0,k0 is the
// // bottom left index of the face cell; face_vec is the
// // normal vector of the face (i.e. 0 is x,1 is y, 2 is z)
//
//    size_t field = (FF==UF) ? 1 : 0;
//
//    vector<double> net_vel(3,0.);
//
//    BoundaryBisec* bf_intersect = GLOBAL_EQ->get_b_intersect(field,0);    // intersect information for field(1) in fluid(0)
//    NodeProp node = GLOBAL_EQ->get_node_property(field,0);                 // node information for field(1)
//
//    size_t_array2D p(dim,dim,0);
//    // arrays of vertex indexes of the cube/square
//    size_t_array2D ix(dim,dim,0);
//    size_t_array2D iy(dim,dim,0);
//    size_t_array2D iz(dim,dim,0);
//    // Min and max of the cell containing ghost point
//    doubleArray2D extents(dim,2,0);
//    // Field value at vertex of face
//    doubleArray2D f(2,2,0);
//    // Interpolated values at the walls
//    doubleArray2D fwall(2,2,0);
//    // Distance of grid/particle walls from the ghost point
//    doubleArray2D del_wall(2,2,0);
//    // Ghost point coordinate
//    doubleVector xghost(3,0);
//
//    // Direction in the plane of face
//    size_t dir1=0, dir2=0;
//
//    if (face_vec == 0) {
//       ix(0,0) = i0,     iy(0,0) = j0,    iz(0,0) = k0;
//       ix(1,0) = i0,     iy(1,0) = j0,    iz(1,0) = k0+1;
//       ix(0,1) = i0,     iy(0,1) = j0+1,  iz(0,1) = k0;
//       ix(1,1) = i0,     iy(1,1) = j0+1,  iz(1,1) = k0+1;
//       dir1 = 2, dir2 = 1;
//    } else if (face_vec == 1) {
//       ix(0,0) = i0,     iy(0,0) = j0,    iz(0,0) = k0;
//       ix(1,0) = i0+1,   iy(1,0) = j0,    iz(1,0) = k0;
//       ix(0,1) = i0,     iy(0,1) = j0,    iz(0,1) = k0+1;
//       ix(1,1) = i0+1,   iy(1,1) = j0,    iz(1,1) = k0+1;
//       dir1 = 0, dir2 = 2;
//    } else if (face_vec == 2) {
//       ix(0,0) = i0,     iy(0,0) = j0,    iz(0,0) = k0;
//       ix(1,0) = i0+1,   iy(1,0) = j0,    iz(1,0) = k0;
//       ix(0,1) = i0,     iy(0,1) = j0+1,  iz(0,1) = k0;
//       ix(1,1) = i0+1,   iy(1,1) = j0+1,  iz(1,1) = k0;
//       dir1 = 0, dir2 = 1;
//    }
//
//    xghost(0) = x0;
//    xghost(1) = y0;
//    xghost(2) = z0;
//
//    for (size_t i=0; i < 2; i++) {
//       for (size_t j=0; j < 2; j++) {
//          // Face vertex index
//          p(i,j) = return_node_index(FF,comp,ix(i,j),iy(i,j),iz(i,j));
//          // Vertex field values
//          f(i,j) = FF->DOF_value( ix(i,j), iy(i,j), iz(i,j), comp, level );
//       }
//    }
//
//    // Min x-coordinate in the grid cell
//    extents(0,0) = FF->get_DOF_coordinate( ix(0,0), comp, 0 ) ;
//    // Max x-coordinate in the grid cell
//    extents(0,1) = FF->get_DOF_coordinate( ix(1,1), comp, 0 ) ;
//    // Min y-coordinate in the grid cell
//    extents(1,0) = FF->get_DOF_coordinate( iy(0,0), comp, 1 ) ;
//    // Max y-coordinate in the grid cell
//    extents(1,1) = FF->get_DOF_coordinate( iy(1,1), comp, 1 ) ;
//
//    if (dim == 3) {
//       // Min z-coordinate in the grid cell
//       extents(2,0) = FF->get_DOF_coordinate( iz(0,0), comp, 2 ) ;
//       // Max z-coordinate in the grid cell
//       extents(2,1) = FF->get_DOF_coordinate( iz(1,1), comp, 2 ) ;
//    }
//
//
//    // Contribution from left and right wall
//    for (size_t i = 0; i < 2; i++) {     // 0 --> left; 1 --> right
//       if ((field == 0) || ((node.void_frac[comp]->item(p(i,0)) == 0) && (node.void_frac[comp]->item(p(i,1)) == 0))) {
//          fwall(0,i) = ((extents(dir2,1) - xghost(dir2))*f(i,0) + (xghost(dir2) - extents(dir2,0))*f(i,1))/(extents(dir2,1)-extents(dir2,0));
//          del_wall(0,i) = MAC::abs(extents(dir1,i) - xghost(dir1));
//       // if bottom vertex is in fluid domain
//       } else if ((node.void_frac[comp]->item(p(i,0)) == 0) && (bf_intersect[dir2].offset[comp]->item(p(i,0),1) == 1)) {
//          double yint = bf_intersect[dir2].value[comp]->item(p(i,0),1);
//          // Condition where intersection distance is more than ghost point distance, it means that the ghost
//          // point can be projected on the particle wall
//          if (yint >= (xghost(dir2)-extents(dir2,0))) {
//             fwall(0,i) = ((extents(dir2,0)+yint-xghost(dir2))*f(i,0)+(xghost(dir2)-extents(dir2,0))*bf_intersect[dir2].field_var[comp]->item(p(i,0),1))/yint;
//             del_wall(0,i) = MAC::abs(extents(dir1,i) - xghost(dir1));
//          // Ghost point cannot be projected on the particle wall, as the solid surface come first
//          } else {
//             size_t id = (size_t) node.parID[comp]->item(p(i,1));
//             if (face_vec > dir2) {
//                del_wall(0,i) = (i==1) ?
//                    find_intersection_for_ghost(FF, xghost(dir1), extents(dir1,i), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, field, level, i) :
//                    find_intersection_for_ghost(FF, extents(dir1,i), xghost(dir1), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, field, level, i) ;
//             } else {
//                del_wall(0,i) = (i==1) ?
//                    find_intersection_for_ghost(FF, xghost(dir1), extents(dir1,i), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, field, level, i) :
//                    find_intersection_for_ghost(FF, extents(dir1,i), xghost(dir1), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, field, level, i) ;
//             }
//
//             if (dir1 == 0) {
//                (i == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)+del_wall(0,i),xghost(1),xghost(2),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)-del_wall(0,i),xghost(1),xghost(2),id) ;
//             } else if (dir1 == 1) {
//                (i == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)+del_wall(0,i),xghost(2),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)-del_wall(0,i),xghost(2),id) ;
//             } else if (dir1 == 2) {
//                (i == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)+del_wall(0,i),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)-del_wall(0,i),id) ;
//             }
//             fwall(0,i) = net_vel[comp];
//          }
//       // if top vertex is in fluid domain
//       } else if ((node.void_frac[comp]->item(p(i,1)) == 0) && (bf_intersect[dir2].offset[comp]->item(p(i,1),0) == 1)) {
//          double yint = bf_intersect[dir2].value[comp]->item(p(i,1),0);
//          // Condition where intersection distance is more than ghost point distance, it means that the ghost
//          // point can be projected on the wall
//          if (yint >= (extents(dir2,1)-xghost(dir2))) {
//             fwall(0,i) = ((xghost(dir2)+yint-extents(dir2,1))*f(i,1)+(extents(dir2,1)-xghost(dir2))*bf_intersect[dir2].field_var[comp]->item(p(i,1),0))/yint;
//             del_wall(0,i) = MAC::abs(extents(dir1,i) - xghost(dir1));
//          // Ghost point cannot be projected on the wall, as the solid surface come first
//          } else {
//             size_t id = (size_t) node.parID[comp]->item(p(i,0));
//             if (face_vec > dir2) {
//                del_wall(0,i) = (i==1) ?
//                    find_intersection_for_ghost(FF, xghost(dir1), extents(dir1,i), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, field, level, i) :
//                    find_intersection_for_ghost(FF, extents(dir1,i), xghost(dir1), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, field, level, i) ;
//             } else {
//                del_wall(0,i) = (i==1) ?
//                    find_intersection_for_ghost(FF, xghost(dir1), extents(dir1,i), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, field, level, i) :
//                    find_intersection_for_ghost(FF, extents(dir1,i), xghost(dir1), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, field, level, i) ;
//             }
//
//             if (dir1 == 0) {
//                (i == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)+del_wall(0,i),xghost(1),xghost(2),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)-del_wall(0,i),xghost(1),xghost(2),id) ;
//             } else if (dir1 == 1) {
//                (i == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)+del_wall(0,i),xghost(2),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)-del_wall(0,i),xghost(2),id) ;
//             } else if (dir1 == 2) {
//                (i == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)+del_wall(0,i),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)-del_wall(0,i),id) ;
//             }
//             fwall(0,i) = net_vel[comp];
//          }
//       // if both vertex's are in solid domain
//       } else if ((node.void_frac[comp]->item(p(i,0)) == 1) && (node.void_frac[comp]->item(p(i,1)) == 1)) {
//          size_t id = (size_t) node.parID[comp]->item(p(i,0));
//
//          if (face_vec > dir2) {
//             del_wall(0,i) = (i==1) ?
//                 find_intersection_for_ghost(FF, xghost(dir1), extents(dir1,i), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, field, level, i) :
//                 find_intersection_for_ghost(FF, extents(dir1,i), xghost(dir1), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, field, level, i) ;
//          } else {
//             del_wall(0,i) = (i==1) ?
//                 find_intersection_for_ghost(FF, xghost(dir1), extents(dir1,i), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, field, level, i) :
//                 find_intersection_for_ghost(FF, extents(dir1,i), xghost(dir1), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, field, level, i) ;
//          }
//
//          if (dir1 == 0) {
//             (i == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)+del_wall(0,i),xghost(1),xghost(2),id) :
//                        impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)-del_wall(0,i),xghost(1),xghost(2),id) ;
//          } else if (dir1 == 1) {
//             (i == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)+del_wall(0,i),xghost(2),id) :
//                        impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)-del_wall(0,i),xghost(2),id) ;
//          } else if (dir1 == 2) {
//             (i == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)+del_wall(0,i),id) :
//                        impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)-del_wall(0,i),id) ;
//          }
//          fwall(0,i) = net_vel[comp];
//       }
//    }
//
//    // Contribution from top and bottom wall
//    for (size_t j = 0; j < 2; j++) {         // 0 --> bottom; 1 --> top
//       if ((field == 0) || ((node.void_frac[comp]->item(p(0,j)) == 0) && (node.void_frac[comp]->item(p(1,j)) == 0))) {
//          fwall(1,j) = ((extents(dir1,1) - xghost(dir1))*f(0,j) + (xghost(dir1) - extents(dir1,0))*f(1,j))/(extents(dir1,1)-extents(dir1,0));
//          del_wall(1,j) = MAC::abs(extents(dir2,j) - xghost(dir2));
//       // if left vertex is in fluid domain
//       } else if ((node.void_frac[comp]->item(p(0,j)) == 0) && (bf_intersect[dir1].offset[comp]->item(p(0,j),1) == 1)) {
//          double xint = bf_intersect[dir1].value[comp]->item(p(0,j),1);
//          if (xint >= (xghost(dir1)-extents(dir1,0))) {
//             fwall(1,j) = ((extents(dir1,0)+xint-xghost(dir1))*f(0,j)+(xghost(dir1)-extents(dir1,0))*bf_intersect[dir1].field_var[comp]->item(p(0,j),1))/xint;
//             del_wall(1,j) = MAC::abs(extents(dir2,j) - xghost(dir2));
//          } else {
//             size_t id = (size_t) node.parID[comp]->item(p(1,j));
//             if (face_vec > dir1) {
//                del_wall(1,j) = (j==1) ?
//                    find_intersection_for_ghost(FF, xghost(dir2), extents(dir2,j), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, field, level, j) :
//                    find_intersection_for_ghost(FF, extents(dir2,j), xghost(dir2), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, field, level, j) ;
//             } else {
//                del_wall(1,j) = (j==1) ?
//                    find_intersection_for_ghost(FF, xghost(dir2), extents(dir2,j), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, field, level, j) :
//                    find_intersection_for_ghost(FF, extents(dir2,j), xghost(dir2), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, field, level, j) ;
//             }
//             if (dir2 == 0) {
//                (j == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)+del_wall(1,j),xghost(1),xghost(2),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)-del_wall(1,j),xghost(1),xghost(2),id) ;
//             } else if (dir2 == 1) {
//                (j == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)+del_wall(1,j),xghost(2),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)-del_wall(1,j),xghost(2),id) ;
//             } else if (dir2 == 2) {
//                (j == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)+del_wall(1,j),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)-del_wall(1,j),id) ;
//             }
//             fwall(1,j) = net_vel[comp];
//          }
//       // if right vertex is in fluid domain
//       } else if ((node.void_frac[comp]->item(p(1,j)) == 0) && (bf_intersect[dir1].offset[comp]->item(p(1,j),0) == 1)) {
//          double xint = bf_intersect[dir1].value[comp]->item(p(1,j),0);
//          if (xint >= (extents(dir1,1)-xghost(dir1))) {
//             fwall(1,j) = ((xghost(dir1)+xint-extents(dir1,1))*f(1,j)+(extents(dir1,1)-xghost(dir1))*bf_intersect[dir1].field_var[comp]->item(p(1,j),0))/xint;
//             del_wall(1,j) = MAC::abs(extents(dir2,j) - xghost(dir2));
//          } else {
//             size_t id = (size_t) node.parID[comp]->item(p(0,j));
//             if (face_vec > dir1) {
//                del_wall(1,j) = (j==1) ?
//                    find_intersection_for_ghost(FF, xghost(dir2), extents(dir2,j), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, field, level, j) :
//                    find_intersection_for_ghost(FF, extents(dir2,j), xghost(dir2), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, field, level, j) ;
//             } else {
//                del_wall(1,j) = (j==1) ?
//                    find_intersection_for_ghost(FF, xghost(dir2), extents(dir2,j), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, field, level, j) :
//                    find_intersection_for_ghost(FF, extents(dir2,j), xghost(dir2), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, field, level, j) ;
//             }
//             if (dir2 == 0) {
//                (j == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)+del_wall(1,j),xghost(1),xghost(2),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)-del_wall(1,j),xghost(1),xghost(2),id) ;
//             } else if (dir2 == 1) {
//                (j == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)+del_wall(1,j),xghost(2),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)-del_wall(1,j),xghost(2),id) ;
//             } else if (dir2 == 2) {
//                (j == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)+del_wall(1,j),id) :
//                           impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)-del_wall(1,j),id) ;
//             }
//             fwall(1,j) = net_vel[comp];
//          }
//       // if both vertex's are in solid domain
//       } else if ((node.void_frac[comp]->item(p(0,j)) == 1) && (node.void_frac[comp]->item(p(1,j)) == 1)) {
//          size_t id = (size_t) node.parID[comp]->item(p(0,j));
//          if (face_vec > dir1) {
//             del_wall(1,j) = (j==1) ?
//                 find_intersection_for_ghost(FF, xghost(dir2), extents(dir2,j), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, field, level, j) :
//                 find_intersection_for_ghost(FF, extents(dir2,j), xghost(dir2), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, field, level, j) ;
//          } else {
//             del_wall(1,j) = (j==1) ?
//                 find_intersection_for_ghost(FF, xghost(dir2), extents(dir2,j), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, field, level, j) :
//                 find_intersection_for_ghost(FF, extents(dir2,j), xghost(dir2), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, field, level, j) ;
//          }
//          if (dir2 == 0) {
//             (j == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)+del_wall(1,j),xghost(1),xghost(2),id) :
//                        impose_solid_velocity_for_ghost(net_vel,comp,xghost(0)-del_wall(1,j),xghost(1),xghost(2),id) ;
//          } else if (dir2 == 1) {
//             (j == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)+del_wall(1,j),xghost(2),id) :
//                        impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1)-del_wall(1,j),xghost(2),id) ;
//          } else if (dir2 == 2) {
//             (j == 1) ? impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)+del_wall(1,j),id) :
//                        impose_solid_velocity_for_ghost(net_vel,comp,xghost(0),xghost(1),xghost(2)-del_wall(1,j),id) ;
//          }
//          fwall(1,j) = net_vel[comp];
//       }
//    }
//
//    double field_value = (1./2.)*((del_wall(0,1)*fwall(0,0) + del_wall(0,0)*fwall(0,1))/(del_wall(0,1)+del_wall(0,0)) +
//                                  (del_wall(1,0)*fwall(1,1) + del_wall(1,1)*fwall(1,0))/(del_wall(1,0)+del_wall(1,1)));
//
//    return (field_value);
//
// }

//---------------------------------------------------------------------------
double
DDS_NavierStokes:: compute_p_component ( size_t const& comp
                                       , size_t const& i
                                       , size_t const& j
                                       , size_t const& k)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: compute_p_component" ) ;
   FV_SHIFT_TRIPLET shift ;
   double value=0.;

   // Here we use shift_staggeredToStaggered to shitf centered to staggered
   // for the following reason: the ith component of UU is staggered with the
   // centered field only in the ith direction, which shares its ith location
   // with the non-ith components of the velocity field
   // Example:
   // For ux, it is staggered with the centered field tf in the x
   // direction only, in the y & z direction, ux and tf have the same location
   // Now, in the x direction, tf is located at the same x position as uy and
   // uz and hence the shift_staggeredToStaggered can be used for tf
   // When interpolating a centered field to a staggered field, we use
   // shift_staggeredToStaggered for each ith component and consider the ith
   // shift in the ith direction only, i.e.:
   // * for ux, use shift_staggeredToStaggered(0) and shift in the x direction
   // with shift.i (the xth component of the shift) only
   // * for uy, use shift_staggeredToStaggered(1) and shift in the y direction
   // with shift.j (the yth component of the shift) only
   // * for uz, use shift_staggeredToStaggered(2) and shift in the z direction
   // with shift.k (the zth component of the shift) only
   shift = UF->shift_staggeredToStaggered( comp ) ;

   double dxC = UF->get_cell_size( i, comp, 0 ) ;
   double dyC = UF->get_cell_size( j, comp, 1 ) ;
   double dzC = (dim == 3) ? UF->get_cell_size( k, comp, 2 ) : 1 ;

   switch (comp) {
     case 0:
        value = (PF->DOF_value( shift.i+i, j, k, 0, 1 )
               - PF->DOF_value( shift.i+i-1, j, k, 0, 1 ))*dyC*dzC;
        break;
     case 1:
        value = (PF->DOF_value( i, shift.j+j, k, 0, 1 )
               - PF->DOF_value( i, shift.j+j-1, k, 0, 1 ))*dxC*dzC;
        break;
     case 2:
        value = (PF->DOF_value( i, j, shift.k+k, 0, 1 )
               - PF->DOF_value( i, j, shift.k+k-1, 0, 1 ))*dxC*dyC;
        break;
   }

   return(value);
}

//---------------------------------------------------------------------------
double
DDS_NavierStokes:: compute_adv_component ( size_t const& comp,
                                           size_t const& i,
                                           size_t const& j,
                                           size_t const& k)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: compute_adv_component" ) ;
   double ugradu = 0., value = 0.;

   if ( AdvectionScheme == "TVD" ) {
      ugradu = assemble_advection_TVD(1,rho,1,i,j,k,comp)
             - rho*UF->DOF_value(i,j,k,comp,1)*divergence_of_U(i,j,k,comp,1);
   } else if ( AdvectionScheme == "Upwind" ) {
      ugradu = assemble_advection_Upwind(1,rho,1,i,j,k,comp)
             - rho*UF->DOF_value(i,j,k,comp,1)*divergence_of_U(i,j,k,comp,1);
   } else if ( AdvectionScheme == "Centered" ) {
      ugradu = assemble_advection_Centered(1,rho,1,i,j,k,comp)
             - rho*UF->DOF_value(i,j,k,comp,1)*divergence_of_U(i,j,k,comp,1);
   }

   if ( AdvectionTimeAccuracy == 1 ) {
      value = ugradu;
   } else {
      value = 1.5*ugradu - 0.5*UF->DOF_value(i,j,k,comp,2);
      UF->set_DOF_value(i,j,k,comp,2,ugradu);
   }

   return(value);
}

//----------------------------------------------------------------------
double
DDS_NavierStokes:: divergence_of_U( size_t const& i
                                  , size_t const& j
                                  , size_t const& k
                                  , size_t const& component
                                  , size_t const& level)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: divergence_of_U" );

   // Parameters
   double flux = 0.;
   double AdvectorValueC = 0., AdvectorValueRi = 0.,
    AdvectorValueLe = 0., AdvectorValueTo = 0., AdvectorValueBo = 0.,
    AdvectorValueFr = 0., AdvectorValueBe = 0, AdvectorValueToLe = 0.,
    AdvectorValueToRi = 0., AdvectorValueBoLe = 0., AdvectorValueBoRi = 0.,
    AdvectorValueFrLe = 0., AdvectorValueFrRi = 0., AdvectorValueBeLe = 0.,
    AdvectorValueBeRi = 0., AdvectorValueFrTo = 0., AdvectorValueFrBo = 0.,
    AdvectorValueBeTo = 0., AdvectorValueBeBo = 0.,
    ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.;

   // Comment: staggered unknowns always have a defined value at +1/-1
   // indices in directions different from their component number,
   // i.e. u in x, v in y and w in z.
   // For instance, if u on the right or left boundary is an unknown with
   // homogeneous Neumann BC, then the flux on the right or left needs special
   // treatment using the center value.
   // Otherwise, whether one of the +1/-1 DOF values is on a
   // boundary or not, and whether that boundary has a Dirichlet or Neumann
   // condition is irrelevant, this +1/-1 DOF always has the right value.
   // For Neumann, this is guaranted by
   // FV_BoundaryCondition:: set_free_DOF_values in
   // FV_DiscreteField:: update_free_DOFs_value or
   // FV_DiscreteField:: add_to_free_DOFs_value
   FV_SHIFT_TRIPLET shift = UF->shift_staggeredToStaggered( component );

   double dxC = UF->get_cell_size(i,component,0) ;
   double dyC = UF->get_cell_size(j,component,1) ;
   double dzC = (dim == 3) ? UF->get_cell_size(k,component,2) : 0;

   AdvectorValueC = UF->DOF_value( i, j, k, component, level );

   // The First Component (u)
   if ( component == 0 ) {
      // Right (U_X)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
         ur = AdvectorValueC ;
      else {
         AdvectorValueRi = UF->DOF_value(i+1, j, k, component, level );
         ur = 0.5 * ( AdvectorValueC + AdvectorValueRi );
      }

      // Left (U_X)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
         ul = AdvectorValueC;
      else {
         AdvectorValueLe = UF->DOF_value(i-1, j, k, component, level );
         ul = 0.5 * ( AdvectorValueC + AdvectorValueLe );
      }

      // Top (U_Y)
      AdvectorValueToLe = UF->DOF_value(i+shift.i-1, j+shift.j, k, 1, level );
      AdvectorValueToRi = UF->DOF_value(i+shift.i, j+shift.j, k, 1, level );
      vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );

      // Bottom (U_Y)
      AdvectorValueBoLe = UF->DOF_value(i+shift.i-1, j+shift.j-1, k, 1, level );
      AdvectorValueBoRi = UF->DOF_value(i+shift.i, j+shift.j-1, k, 1, level );
      vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );

      if (dim == 3) {
         // Front (U_Z)
         AdvectorValueFrLe = UF->DOF_value(i+shift.i-1, j, k+shift.k, 2, level );
         AdvectorValueFrRi = UF->DOF_value(i+shift.i, j, k+shift.k, 2, level );
         wf = 0.5 * ( AdvectorValueFrLe + AdvectorValueFrRi );

         // Behind (U_Z)
         AdvectorValueBeLe = UF->DOF_value(i+shift.i-1, j, k+shift.k-1, 2, level );
         AdvectorValueBeRi = UF->DOF_value(i+shift.i, j, k+shift.k-1, 2, level );
         wb = 0.5 * ( AdvectorValueBeLe + AdvectorValueBeRi );
      }
   } else if (component == 1) {
      // The second Component (v)
      // Right (V_X)
      AdvectorValueToRi = UF->DOF_value(i+shift.i, j+shift.j, k, 0, level );
      AdvectorValueBoRi = UF->DOF_value(i+shift.i, j+shift.j-1, k, 0, level );
      ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );

      // Left (V_X)
      AdvectorValueToLe = UF->DOF_value(i+shift.i-1, j+shift.j, k, 0, level );
      AdvectorValueBoLe = UF->DOF_value(i+shift.i-1, j+shift.j-1, k, 0, level );
      ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );

      // Top (V_Y)
      if ( UF->DOF_color(i, j, k, component ) == FV_BC_TOP )
         vt = AdvectorValueC;
      else {
         AdvectorValueTo = UF->DOF_value(i, j+1, k, component, level );
         vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
      }

      // Bottom (V_Y)
      if ( UF->DOF_color(i, j, k, component ) == FV_BC_BOTTOM )
         vb = AdvectorValueC;
      else {
         AdvectorValueBo = UF->DOF_value(i, j-1, k, component, level );
         vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
      }

      if (dim == 3) {
         // Front (V_Z)
         AdvectorValueFrTo = UF->DOF_value(i, j+shift.j, k+shift.k, 2, level );
         AdvectorValueFrBo = UF->DOF_value(i, j+shift.j-1, k+shift.k, 2, level );
         wf = 0.5 * ( AdvectorValueFrTo + AdvectorValueFrBo );

         // Behind (V_Z)
         AdvectorValueBeTo = UF->DOF_value(i, j+shift.j, k+shift.k-1, 2, level );
         AdvectorValueBeBo = UF->DOF_value(i, j+shift.j-1, k+shift.k-1, 2, level );
         wb = 0.5 * ( AdvectorValueBeTo + AdvectorValueBeBo );
      }
   } else {
      // The Third Component (w)
      // Right (W_X)
      AdvectorValueFrRi = UF->DOF_value(i+shift.i, j, k+shift.k, 0, level );
      AdvectorValueBeRi = UF->DOF_value(i+shift.i, j, k+shift.k-1, 0, level );
      ur = 0.5 * ( AdvectorValueFrRi + AdvectorValueBeRi );

      // Left (W_X)
      AdvectorValueFrLe = UF->DOF_value(i+shift.i-1, j, k+shift.k, 0, level );
      AdvectorValueBeLe = UF->DOF_value(i+shift.i-1, j, k+shift.k-1, 0, level );
      ul = 0.5 * ( AdvectorValueFrLe + AdvectorValueBeLe );

      // Top (W_Y)
      AdvectorValueFrTo = UF->DOF_value(i, j+shift.j, k+shift.k, 1, level );
      AdvectorValueBeTo = UF->DOF_value(i, j+shift.j, k+shift.k-1, 1, level );
      vt = 0.5 * ( AdvectorValueFrTo + AdvectorValueBeTo );

      // Bottom (W_Y)
      AdvectorValueFrBo = UF->DOF_value(i, j+shift.j-1, k+shift.k, 1, level );
      AdvectorValueBeBo = UF->DOF_value(i, j+shift.j-1, k+shift.k-1, 1, level );
      vb = 0.5 * ( AdvectorValueFrBo + AdvectorValueBeBo );

      // Front (W_Z)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
         wf = AdvectorValueC;
      else {
         AdvectorValueFr = UF->DOF_value(i, j, k+1, component, level );
         wf = 0.5 * ( AdvectorValueFr + AdvectorValueC );
      }

      // Behind (W_Z)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND )
         wb = AdvectorValueC;
      else {
         AdvectorValueBe = UF->DOF_value(i, j, k-1, component, level );
         wb = 0.5 * ( AdvectorValueBe + AdvectorValueC );
      }
   }

   if (dim == 2) {
      flux = ((vt - vb) * dxC + (ur - ul) * dyC);
   } else if (dim == 3) {
      flux = (vt - vb) * dxC * dzC + (ur - ul) * dyC * dzC + (wf - wb) * dxC * dyC;
   }
   return ( flux );
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_DS_un_at_rhs ( FV_TimeIterator const* t_it,
                                           double const& gamma)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: assemble_DS_un_at_rhs" ) ;

   // Assemble the diffusive components of velocity once in each iteration
   assemble_velocity_diffusion_terms ( );

   double bodyterm=0.;
   size_t cpp = 10;

   // Periodic pressure gradient
   if ( UF->primary_grid()->is_periodic_flow() ) {
      cpp = UF->primary_grid()->get_periodic_flow_direction() ;
      bodyterm = UF->primary_grid()->get_periodic_pressure_drop() /
               ( UF->primary_grid()->get_main_domain_max_coordinate( cpp )
               - UF->primary_grid()->get_main_domain_min_coordinate( cpp ) ) ;
   }

   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);

   min_unknown_index(2) = (dim == 3) ? 0 : 0;
   max_unknown_index(2) = (dim == 3) ? 0 : 1;

   vector<doubleVector*> vel_diffusion = GLOBAL_EQ->get_velocity_diffusion();
   size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);

   for (size_t comp=0;comp<nb_comps[1];comp++) {
      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) =
                        UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) =
                        UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
         double dxC = UF->get_cell_size( i, comp, 0 ) ;
         for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
            double dyC = UF->get_cell_size( j, comp, 1 ) ;
            for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
               double dzC = (dim == 3) ? UF->get_cell_size( k, comp, 2 ) : 1 ;

               size_t p = UF->DOF_local_number(i,j,k,comp);
               // Dxx for un
               double xvalue = vel_diffusion[0]->operator()(p);
               // Dyy for un
               double yvalue = vel_diffusion[1]->operator()(p);
               // Dzz for un
               double zvalue = (dim == 3) ? vel_diffusion[2]->operator()(p) : 0;
               // Pressure contribution
               double pvalue = compute_p_component(comp,i,j,k);
               // Advection contribution
               double adv_value = compute_adv_component(comp,i,j,k);
               //
               if (is_solids) {
                  if (void_frac->operator()(p) == 1) {
                     pvalue = 0.; adv_value = 0.;
                  }
               }

               double rhs = gamma*(xvalue*dyC*dzC
                                 + yvalue*dxC*dzC
                                 + zvalue*dxC*dyC)
                          - pvalue - adv_value
                          + (UF->DOF_value( i, j, k, comp, 1 )*dxC*dyC*dzC*rho)
                                                         /(t_it -> time_step());


               if ( cpp==comp ) rhs += - bodyterm*dxC*dyC*dzC;

               if (is_solids) {
                  if (void_frac->operator()(p) == 1) {
                     if ( cpp==comp ) rhs += bodyterm*dxC*dyC*dzC;
                  }
               }

               UF->set_DOF_value( i, j, k, comp, 0,
                                  rhs*(t_it -> time_step())/(dxC*dyC*dzC*rho));
            }
         }
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: Solve_i_in_jk ( FV_DiscreteField* FF
                                 , FV_TimeIterator const* t_it
                                 , size_t const& dir_i
                                 , size_t const& dir_j
                                 , size_t const& dir_k
                                 , double const& gamma)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_NavierStokes:: Solve_i_in_jk" ) ;

  size_t field = (FF == PF) ? 0 : 1 ;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  for (size_t comp=0;comp<nb_comps[field];comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) =
                        FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) =
                        FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_unknown_index(dir_k);
        local_max_k = max_unknown_index(dir_k);
     }

     LocalVector* VEC = GLOBAL_EQ->get_VEC(field) ;
     TDMatrix* A = GLOBAL_EQ->get_A(field);
     LA_SeqMatrix* row_index = GLOBAL_EQ->get_row_indexes(field,dir_i,comp);

     // Solve in i
     if ((nb_ranks_comm_i[dir_i]>1)||(is_periodic[field][dir_i] == 1)) {
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j){
           for (size_t k=local_min_k; k <= local_max_k; ++k) {
              size_t r_index = row_index->item(j,k);
              // Assemble fi and return fe for each proc locally
              double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir_i,field);
              // Calculate Aei*ui in each proc locally
              compute_Aei_ui(A,VEC,comp,dir_i,r_index);
              // Pack Aei_ui and fe for sending it to master
              data_packing (FF,j,k,fe,comp,dir_i);
           }
        }
        solve_interface_unknowns ( FF, gamma, t_it, comp, dir_i);

     } else if (is_periodic[field][dir_i] == 0) {
        // Serial mode with non-periodic condition
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j){
           for (size_t k=local_min_k; k <= local_max_k; ++k) {
              size_t r_index = row_index->item(j,k);
              assemble_local_rhs(j,k,gamma,t_it,comp,dir_i,field);
              GLOBAL_EQ->DS_NavierStokes_solver(FF,j,k,min_unknown_index(dir_i)
                                                           ,comp,dir_i,r_index);
           }
        }
     }
  }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: data_packing ( FV_DiscreteField const* FF
                                , size_t const& j
                                , size_t const& k
                                , double const& fe
                                , size_t const& comp
                                , size_t const& dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: data_packing" ) ;

   size_t field = (FF == PF) ? 0 : 1 ;

   LocalVector* VEC = GLOBAL_EQ->get_VEC(field) ;
   LA_SeqMatrix* row_index = GLOBAL_EQ->get_row_indexes(field,dir,comp);

   double *packed_data = first_pass[field][dir].send[comp][rank_in_i[dir]];

   // Pack the data
   size_t vec_pos = row_index->item(j,k);

   if (rank_in_i[dir] == 0) {
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_periodic[field][dir])
          packed_data[3*vec_pos+0] =
                        VEC[dir].T[comp]->item(nb_ranks_comm_i[dir]-1);
      else
          packed_data[3*vec_pos+0] = 0;

      packed_data[3*vec_pos+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);

   } else if (rank_in_i[dir] == nb_ranks_comm_i[dir]-1) {
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_periodic[field][dir])
          packed_data[3*vec_pos+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);
      else
          packed_data[3*vec_pos+1]=0;

      packed_data[3*vec_pos+0] = VEC[dir].T[comp]->item(rank_in_i[dir]-1);

   } else {
      packed_data[3*vec_pos+0] = VEC[dir].T[comp]->item(rank_in_i[dir]-1);
      packed_data[3*vec_pos+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);
   }

   // Send the fe values and 0 for last proc
   packed_data[3*vec_pos+2] = fe;

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: compute_Aei_ui (struct TDMatrix* arr
                                 , struct LocalVector* VEC
                                 , size_t const& comp
                                 , size_t const& dir
                                 , size_t const& r_index)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: compute_Aei_ui" ) ;
   // create a replica of local rhs vector in local solution vector
   for (size_t i = 0; i < VEC[dir].local_T[comp]->nb_rows(); i++){
      VEC[dir].local_solution_T[comp]
                  ->set_item(i,VEC[dir].local_T[comp]->item(i));
   }

   // Solve for ui locally and put it in local solution vector
   GLOBAL_EQ->mod_thomas_algorithm(arr, VEC[dir].local_solution_T[comp]
                                                      , comp, dir,r_index);

   for (size_t i = 0; i < VEC[dir].T[comp]->nb_rows(); i++){
      VEC[dir].T[comp]->set_item(i,0);
   }

   // Calculate Aei*ui in each proc locally and put it in T vector
   arr[dir].ei[comp][r_index]
                  ->multiply_vec_then_add(VEC[dir].local_solution_T[comp]
                                                         ,VEC[dir].T[comp]);

}

//---------------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_local_rhs ( size_t const& j
                                      , size_t const& k
                                      , double const& gamma
                                      , FV_TimeIterator const* t_it
                                      , size_t const& comp
                                      , size_t const& dir
                                      , size_t const& field )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: assemble_local_rhs" ) ;
   double fe = 0.;
   if (field == 0) {
      fe = pressure_local_rhs(j,k,t_it,dir);
   } else if (field == 1) {
      fe = velocity_local_rhs(j,k,gamma,t_it,comp,dir);
   }
   return(fe);
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: NS_velocity_update ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: NS_velocity_update" ) ;

   double gamma=mu;

   assemble_DS_un_at_rhs(t_it,gamma);
   // Update gamma based for invidual direction
   gamma = mu/2.0;

   Solve_i_in_jk(UF,t_it,0,1,2,gamma);
   // Synchronize the distributed DS solution vector
   GLOBAL_EQ->synchronize_DS_solution_vec();
   // Tranfer back to field
   UF->update_free_DOFs_value( 3, GLOBAL_EQ->get_solution_DS_velocity() ) ;
   if (is_solids) initialize_grid_nodes_on_rigidbody({3});

   Solve_i_in_jk(UF,t_it,1,0,2,gamma);
   // Synchronize the distributed DS solution vector
   GLOBAL_EQ->synchronize_DS_solution_vec();
   // Tranfer back to field
   if (dim == 2) {
      UF->update_free_DOFs_value( 0 , GLOBAL_EQ->get_solution_DS_velocity() ) ;
      if (is_solids) initialize_grid_nodes_on_rigidbody({0});
   } else if (dim == 3) {
      UF->update_free_DOFs_value( 4 , GLOBAL_EQ->get_solution_DS_velocity() ) ;
      if (is_solids) initialize_grid_nodes_on_rigidbody({4});
   }

   if (dim == 3) {
      Solve_i_in_jk(UF,t_it,2,0,1,gamma);
      // Synchronize the distributed DS solution vector
      GLOBAL_EQ->synchronize_DS_solution_vec();
      // Tranfer back to field
      UF->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_DS_velocity() ) ;
      if (is_solids) initialize_grid_nodes_on_rigidbody({0});
   }

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes::initialize_grid_nodes_on_rigidbody( vector<size_t> const& list )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes::initialize_grid_nodes_on_rigidbody" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  // Vector for solid presence
  size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
  size_t_vector* parID = allrigidbodies->get_rigidbodyIDs_on_grid(UF);

  for (size_t comp = 0; comp < nb_comps[1]; comp++) {
     // Get local min and max indices
     for (size_t dir = 0; dir < dim; dir++) {
        if (is_periodic[1][dir]) {
           min_unknown_index(dir) =
                     UF->get_min_index_unknown_handled_by_proc( comp, dir ) - 1;
           max_unknown_index(dir) =
                     UF->get_max_index_unknown_handled_by_proc( comp, dir ) + 1;
        } else {
           min_unknown_index(dir) =
                     UF->get_min_index_unknown_handled_by_proc( comp, dir );
           max_unknown_index(dir) =
                     UF->get_max_index_unknown_handled_by_proc( comp, dir );
        }
     }

     for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
        double xC = UF->get_DOF_coordinate( i, comp, 0 ) ;
        for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
           double yC = UF->get_DOF_coordinate( j, comp, 1 ) ;
           for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
              double zC = (dim == 2) ? 0 : UF->get_DOF_coordinate( k, comp, 2 );
              geomVector pt(xC,yC,zC);
              size_t p = UF->DOF_local_number(i,j,k,comp);
              if (void_frac->operator()(p) == 1) {
                 size_t par_id = parID->operator()(p);
                 geomVector rb_vel = allrigidbodies->rigid_body_velocity(par_id,pt);
                 for (size_t level : list)
                  UF->set_DOF_value( i, j, k, comp, level,rb_vel(comp));
              }
           }
        }
     }
  }
}




//---------------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_velocity_gradients (class doubleVector& grad
                                              , size_t const& i
                                              , size_t const& j
                                              , size_t const& k
                                              , size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: assemble_velocity_gradients" ) ;

   FV_SHIFT_TRIPLET shift = PF->shift_staggeredToCentered() ;

   size_t comp = 0;

   size_t_array2D* intersect_vector =
                           allrigidbodies->get_intersect_vector_on_grid(PF);
   doubleArray2D* intersect_distance =
                           allrigidbodies->get_intersect_distance_on_grid(PF);
   doubleArray2D* intersect_fieldVal =
                           allrigidbodies->get_intersect_fieldValue_on_grid(PF);
   size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(PF);

   size_t p = PF->DOF_local_number(i,j,k,comp);

   // Dxx for un
   double xh= UF->get_DOF_coordinate( shift.i+i,0, 0 )
            - UF->get_DOF_coordinate( shift.i+i-1, 0, 0 ) ;
   double xvalue = UF->DOF_value( shift.i+i, j, k, 0, level )
                 - UF->DOF_value( shift.i+i-1, j, k, 0, level) ;
   // Dyy for un
   double yvalue = UF->DOF_value( i, shift.j+j, k, 1, level)
                 - UF->DOF_value( i, shift.j+j-1, k, 1, level) ;
   double yh= UF->get_DOF_coordinate( shift.j+j,1, 1 )
            - UF->get_DOF_coordinate( shift.j+j-1, 1, 1 ) ;

   double bx = xh;
   double by = yh;

   if (is_solids) {
      if (void_frac->operator()(p) == 0) {
         if (intersect_vector->operator()(p,2*0+0) == 1) {
            xvalue = UF->DOF_value( shift.i+i, j, k, 0, level)
                   - intersect_fieldVal->operator()(p,2*0+0);
            xh = intersect_distance->operator()(p,2*0+0)
               + PF->get_cell_size( i, 0, 0 )/2.;
         }
         if (intersect_vector->operator()(p,2*0+1) == 1) {
            xvalue = intersect_fieldVal->operator()(p,2*0+1)
                   - UF->DOF_value( shift.i+i-1, j, k, 0, level);
            xh = intersect_distance->operator()(p,2*0+1)
               + PF->get_cell_size( i, 0, 0 )/2.;
         }
         if ((intersect_vector->operator()(p,2*0+1) == 1)
          && (intersect_vector->operator()(p,2*0+0) == 1)) {
            xvalue = intersect_fieldVal->operator()(p,2*0+1)
                   - intersect_fieldVal->operator()(p,2*0+0);
            xh = intersect_distance->operator()(p,2*0+1)
               + intersect_distance->operator()(p,2*0+0);
         }
      } else {
         xvalue = 0.;
      }
   }

   if (is_solids) {
      if (void_frac->operator()(p) == 0) {
         if (intersect_vector->operator()(p,2*1+0) == 1) {
            yvalue = UF->DOF_value( i, shift.j+j, k, 1, level)
                   - intersect_fieldVal->operator()(p,2*1+0);
            yh = intersect_distance->operator()(p,2*1+0)
               + PF->get_cell_size( j, 0, 1 )/2.;
         }
         if (intersect_vector->operator()(p,2*1+1) == 1) {
            yvalue = intersect_fieldVal->operator()(p,2*1+1)
                   - UF->DOF_value( i,shift.j+j-1, k, 1, level);
            yh = intersect_distance->operator()(p,2*1+1)
               + PF->get_cell_size( j, 0, 1 )/2.;
         }
         if ((intersect_vector->operator()(p,2*1+1) == 1)
          && (intersect_vector->operator()(p,2*1+0) == 1)) {
            yvalue = intersect_fieldVal->operator()(p,2*1+1)
                   - intersect_fieldVal->operator()(p,2*1+0);
            yh = intersect_distance->operator()(p,2*1+1)
               + intersect_distance->operator()(p,2*1+0);
         }
      } else {
         yvalue = 0.;
      }
   }

   bx = xh/bx;
   by = yh/by;

   double beta = min(1.,min(bx,by));

   if (dim == 3) {
      // Dzz for un
      double zh = UF->get_DOF_coordinate( shift.k+k,2, 2 )
                - UF->get_DOF_coordinate( shift.k+k-1, 2, 2 ) ;
      double zvalue = UF->DOF_value( i, j, shift.k+k, 2, level)
                    - UF->DOF_value( i, j, shift.k+k-1, 2, level) ;

      double bz = zh;

      if (is_solids) {
         if (void_frac->operator()(p) == 0) {
            if (intersect_vector->operator()(p,2*2+0) == 1) {
               zvalue = UF->DOF_value( i, j, shift.k+k, 2, level)
                      - intersect_fieldVal->operator()(p,2*2+0);
               zh = intersect_distance->operator()(p,2*2+0)
                  + PF->get_cell_size( k, 0, 2 )/2.;
            }
            if (intersect_vector->operator()(p,2*2+1) == 1) {
               zvalue = intersect_fieldVal->operator()(p,2*2+1)
                      - UF->DOF_value( i, j,shift.k+k-1, 2, level);
               zh = intersect_distance->operator()(p,2*2+1)
                  + PF->get_cell_size( k, 0, 2 )/2.;
            }
            if ((intersect_vector->operator()(p,2*2+1) == 1)
             && (intersect_vector->operator()(p,2*2+0) == 1)) {
               zvalue = intersect_fieldVal->operator()(p,2*2+1)
                      - intersect_fieldVal->operator()(p,2*2+0);
               zh = intersect_distance->operator()(p,2*2+1)
                  + intersect_distance->operator()(p,2*2+0);
            }
         } else {
            zvalue = 0.;
         }
      }

      bz = zh/bz;

      grad(2) = zvalue/zh;

      beta = min(1.,min(bx,min(by,bz)));
   }

   grad(0) = xvalue/xh;
   grad(1) = yvalue/yh;

   return(beta);

}


//---------------------------------------------------------------------------
double
DDS_NavierStokes:: calculate_velocity_divergence ( size_t const& i, size_t const& j, size_t const& k, size_t const& level, FV_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: calculate_velocity_divergence" ) ;

   doubleVector grad(3,0);

   double beta = assemble_velocity_gradients(grad,i,j,k,level);

   double value = beta*(grad(0) + grad(1) + grad(2));

   return(value);
}

//---------------------------------------------------------------------------
double
DDS_NavierStokes:: pressure_local_rhs ( size_t const& j
                                      , size_t const& k
                                      , FV_TimeIterator const* t_it
                                      , size_t const& dir)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: pressure_local_rhs" ) ;
   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
      max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
   }

   size_t pos;

   // Compute VEC_rhs_x = rhs in x
   double fe=0.;
   double value=0.;

   // Vector for fi
   LocalVector* VEC = GLOBAL_EQ->get_VEC(0);
   doubleVector* divergence = GLOBAL_EQ->get_node_divergence(0);

   for (size_t i=min_unknown_index(dir);i<=max_unknown_index(dir);++i) {
      double dx = PF->get_cell_size( i, 0, dir );
      if (dir == 0) {
         double vel_div = calculate_velocity_divergence(i,j,k,0,t_it);
         value = -(rho*vel_div*dx)/(t_it -> time_step());
         size_t p = PF->DOF_local_number(i,j,k,0);
         divergence->operator()(p) = vel_div;
      } else if (dir == 1) {
         value = PF->DOF_value( j, i, k, 0, 1 )*dx;
      } else if (dir == 2) {
         value = PF->DOF_value( j, k, i, 0, 1 )*dx;
      }

      pos = i - min_unknown_index(dir);

      if (is_periodic[0][dir] == 0) {
         if (rank_in_i[dir] == nb_ranks_comm_i[dir]-1) {
            VEC[dir].local_T[0]->set_item( pos, value);
         } else {
            if (i == max_unknown_index(dir))
               fe = value;
            else
               VEC[dir].local_T[0]->set_item( pos, value);
         }
      } else {
         if (i == max_unknown_index(dir))
            fe = value;
         else
            VEC[dir].local_T[0]->set_item( pos, value);
      }
   }

   return fe;
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: correct_pressure_1st_layer_solid (size_t const& level )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: correct_pressure_1st_layer_solid" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  NodeProp node = GLOBAL_EQ->get_node_property(0,0);
  size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(PF);

  size_t comp = 0;

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
           size_t p = PF->DOF_local_number(i,j,k,comp);
           node.bound_cell->set_item(p,0.);
           if (void_frac->operator()(p) == 1) {
              double value = 0., count = 0.;
              size_t p1 = PF->DOF_local_number(i+1,j,k,comp);
              size_t p2 = PF->DOF_local_number(i-1,j,k,comp);
              size_t p3 = PF->DOF_local_number(i,j+1,k,comp);
              size_t p4 = PF->DOF_local_number(i,j-1,k,comp);
              if (void_frac->operator()(p1) == 0) {
                 node.bound_cell->set_item(p,1.);
                 value += PF->DOF_value( i+1, j, k, comp, level );
                 count += 1.;
              }
              if (void_frac->operator()(p2) == 0) {
                 node.bound_cell->set_item(p,1.);
                 value += PF->DOF_value( i-1, j, k, comp, level );
                 count += 1.;
              }
              if (void_frac->operator()(p3) == 0) {
                 node.bound_cell->set_item(p,1.);
                 value += PF->DOF_value( i, j+1, k, comp, level );
                 count += 1.;
              }
              if (void_frac->operator()(p4) == 0) {
                 node.bound_cell->set_item(p,1.);
                 value += PF->DOF_value( i, j-1, k, comp, level );
                 count += 1.;
              }

              if (dim == 3) {
                 size_t p5 = PF->DOF_local_number(i,j,k+1,comp);
                 size_t p6 = PF->DOF_local_number(i,j,k-1,comp);

                 if (void_frac->operator()(p5) == 0) {
                    node.bound_cell->set_item(p,1.);
                    value += PF->DOF_value( i, j, k+1, comp, level );
                    count += 1.;
                 }
                 if (void_frac->operator()(p6) == 0) {
                    node.bound_cell->set_item(p,1.);
                    value += PF->DOF_value( i, j, k-1, comp, level );
                    count += 1.;
                 }
              }

              if (count != 0.) value = value/count;
              GLOBAL_EQ->update_global_P_vector(i,j,k,value);
           }
        }
     }
  }

  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec_P();
  // Tranfer back to field
  PF->update_free_DOFs_value( level, GLOBAL_EQ->get_solution_DS_pressure() ) ;
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: correct_pressure_2nd_layer_solid (size_t const& level )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: correct_pressure_2nd_layer_solid" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  NodeProp node = GLOBAL_EQ->get_node_property(0,0);
  size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(PF);

  size_t comp = 0;

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1);j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2);k <= max_unknown_index(2); ++k) {
           size_t p = PF->DOF_local_number(i,j,k,comp);
           if ((void_frac->operator()(p) == 1)
            && (node.bound_cell->item(p) == 0)) {
              double value = 0., count = 0.;
              size_t p1 = PF->DOF_local_number(i+1,j,k,comp);
              size_t p2 = PF->DOF_local_number(i-1,j,k,comp);
              size_t p3 = PF->DOF_local_number(i,j+1,k,comp);
              size_t p4 = PF->DOF_local_number(i,j-1,k,comp);
              if (node.bound_cell->item(p1) == 1.) {
                 node.bound_cell->set_item(p,2.);
                 value += PF->DOF_value( i+1, j, k, comp, level );
                 count += 1.;
              }
              if (node.bound_cell->item(p2) == 1.) {
                 node.bound_cell->set_item(p,2.);
                 value += PF->DOF_value( i-1, j, k, comp, level );
                 count += 1.;
              }
              if (node.bound_cell->item(p3) == 1.) {
                 node.bound_cell->set_item(p,2.);
                 value += PF->DOF_value( i, j+1, k, comp, level );
                 count += 1.;
              }
              if (node.bound_cell->item(p4) == 1.) {
                 node.bound_cell->set_item(p,2.);
                 value += PF->DOF_value( i, j-1, k, comp, level );
                 count += 1.;
              }

              if (dim == 3) {
                 size_t p5 = PF->DOF_local_number(i,j,k+1,comp);
                 size_t p6 = PF->DOF_local_number(i,j,k-1,comp);
                 if (node.bound_cell->item(p5) == 1.) {
                    node.bound_cell->set_item(p,2.);
                    value += PF->DOF_value( i, j, k+1, comp, level );
                    count += 1.;
                 }
                 if (node.bound_cell->item(p6) == 1.) {
                    node.bound_cell->set_item(p,2.);
                    value += PF->DOF_value( i, j, k-1, comp, level );
                    count += 1.;
                 }
              }

              if (count != 0.) value = value/count;
              GLOBAL_EQ->update_global_P_vector(i,j,k,value);
           }
        }
     }
  }
  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec_P();
  // Tranfer back to field
  PF->update_free_DOFs_value( level, GLOBAL_EQ->get_solution_DS_pressure() ) ;
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: correct_mean_pressure (size_t const& level )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: correct_mean_pressure" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(PF);
  size_t comp = 0;

  double mean=0.,nb_global_unknown=0.;

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
           if (is_solids) {
              size_t p = PF->DOF_local_number(i,j,k,comp);
              if (void_frac->operator()(p) == 0) {
                 mean += PF->DOF_value( i, j, k, comp, level );
                 nb_global_unknown += 1.;
              }
           }
        }
     }
  }

  mean = macCOMM->sum( mean ) ;
  nb_global_unknown = macCOMM->sum( nb_global_unknown ) ;

  mean = mean/nb_global_unknown;

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
           if (is_solids) {
              size_t p = PF->DOF_local_number(i,j,k,comp);
              if (void_frac->operator()(p) == 0) {
                 double value = PF->DOF_value( i, j, k, comp, level );
                 GLOBAL_EQ->update_global_P_vector(i,j,k,value-mean);
              }
           }
        }
     }
  }

  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec_P();
  // Tranfer back to field
  PF->update_free_DOFs_value( level, GLOBAL_EQ->get_solution_DS_pressure() ) ;

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: NS_pressure_update ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: NS_pressure_update" ) ;

  double gamma=mu/2.0;

  Solve_i_in_jk (PF,t_it,0,1,2,gamma);
  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec_P();
  // Tranfer back to field
  PF->update_free_DOFs_value( 1, GLOBAL_EQ->get_solution_DS_pressure() ) ;

  Solve_i_in_jk (PF,t_it,1,0,2,gamma);
  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec_P();
  // Tranfer back to field
  PF->update_free_DOFs_value( 1, GLOBAL_EQ->get_solution_DS_pressure() ) ;

  if (dim == 3) {
     Solve_i_in_jk (PF,t_it,2,0,1,gamma);
     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec_P();
     // Tranfer back to field
     PF->update_free_DOFs_value( 1, GLOBAL_EQ->get_solution_DS_pressure() ) ;
  }

  if (PF->all_BCs_nonDirichlet(0)) {
     correct_mean_pressure(1);
  }

  if (is_solids) {
     correct_pressure_1st_layer_solid(1);
     correct_pressure_2nd_layer_solid(1);
  }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: NS_final_step ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: NS_final_step" ) ;

   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);

   doubleVector* divergence = GLOBAL_EQ->get_node_divergence(0);
   doubleVector* divergence_old = GLOBAL_EQ->get_node_divergence(1);

   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
      max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
   }

   for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
      for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
         for (size_t k = min_unknown_index(2);k <= max_unknown_index(2); ++k) {
            size_t p = PF->DOF_local_number(i,j,k,0);
            double vel_div0 = divergence->operator()(p);
            double vel_div1 = divergence_old->operator()(p);

            // Assemble the bodyterm
            double value = PF->DOF_value( i, j, k, 0, 0 )
                         + PF->DOF_value( i, j, k, 0, 1 )
                         - 0.5*kai*mu*(vel_div0+vel_div1);
            GLOBAL_EQ->update_global_P_vector(i,j,k,value);
         }
      }
   }

   // Synchronize the distributed DS solution vector
   GLOBAL_EQ->synchronize_DS_solution_vec_P();
   // Tranfer back to field
   PF->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_DS_pressure() ) ;

   if (PF->all_BCs_nonDirichlet(0)) {
      correct_mean_pressure(0);
   }

   // Store the divergence to be used in the next time iteration
   for (size_t i = 0; i < divergence->size(); i++)
      divergence_old->operator()(i) = divergence->operator()(i);

   if (is_solids) {
      correct_pressure_1st_layer_solid(0);
      correct_pressure_2nd_layer_solid(0);
   }
   // Propagate values to the boundaries depending on BC conditions
   PF->set_neumann_DOF_values();
}

//
// //----------------------------------------------------------------------
// void
// DDS_NavierStokes::write_surface_discretized_forces(
// 		FV_DiscreteField const* FF,
// 		size_t const& Np,
// 		size_t const& parID,
// 		class doubleArray2D& point,
// 		class doubleArray2D& force,
// 		FV_TimeIterator const* t_it)
// //----------------------------------------------------------------------
// {
//
//   if ((is_surfacestressOUT) && (t_it->iteration_number()%100 == 0)) {
//      // Collect particle data from all the procs, if any.
//      macCOMM->sum_array(force);
//
//      if (my_rank == 0) {
//         ofstream outputFile ;
//
//         std::ostringstream os2;
//         if (FF == PF) {
//            os2 << "./DS_results/surface_press_force_parID_" << parID << ".csv";
//         } else {
//            os2 << "./DS_results/surface_visco_force_parID_" << parID << ".csv";
//         }
//         std::string filename = os2.str();
//         outputFile.open(filename.c_str());
//
//         outputFile << "Sid,x,y,z,Fx,Fy,Fz" << endl;
//
//         for (size_t i=0;i<Np;i++) {
//            outputFile << i << "," << point(i,0) << ","
//                                   << point(i,1) << ","
//                                   << point(i,2) << ","
// 				  << force(i,0) << ","
//                                   << force(i,1) << ","
//                                   << force(i,2) << endl;
//         }
//         outputFile.close();
//      }
//   }
//
// }
//
//
// //----------------------------------------------------------------------
// void
// DDS_NavierStokes::write_surface_discretization(size_t const& Np)
// //----------------------------------------------------------------------
// {
//   ofstream outputFile ;
//
//   std::ostringstream os2;
//   os2 << "./DS_results/surface_discretization.csv";
//   std::string filename = os2.str();
//   outputFile.open(filename.c_str());
//
//   outputFile << "x,y,z" << endl;
//
//   // Structure of particle surface discretization
//   SurfaceDiscretize surface = GLOBAL_EQ->get_surface();
//
//   for (size_t i=0;i<Np;i++) {
//      outputFile << surface.coordinate[0]->item(i) << ","
// 		<< surface.coordinate[1]->item(i) << ","
// 		<< surface.coordinate[2]->item(i) << endl;
//   }
//   outputFile.close();
// }
//
//
// //----------------------------------------------------------------------
// void
// DDS_NavierStokes::write_output_field(FV_DiscreteField const* FF, FV_TimeIterator const* t_it)
// //----------------------------------------------------------------------
// {
//   ofstream outputFile ;
//
//   std::ostringstream os2;
//   os2 << "./DS_results/outputNS_" << t_it->iteration_number() << ".csv";
//   std::string filename = os2.str();
//   outputFile.open(filename.c_str());
//
//   size_t field = (FF == UF) ? 1 : 0 ;
//
//   size_t i,j,k;
//   outputFile << "x,y,z,par_ID,void_frac,left,lv,right,rv,bottom,bov,top,tv" << endl;//,behind,bev,front,fv" << endl;
// //  outputFile << "x,y,z,id,lambdax,lambday,void_frac,fresh,neighx,neighy,count_fresh,counter_neigh,div,vel" << endl;//,behind,bev,front,fv" << endl;
//
//   size_t_vector min_index(dim,0);
//   size_t_vector max_index(dim,0);
//
//   NodeProp node = GLOBAL_EQ->get_node_property(field,0);
// //  FreshNode* fresh = GLOBAL_EQ->get_fresh_node();
// //  DivNode* divergence = GLOBAL_EQ->get_node_divergence();
//   BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(field,0);
//
//   for (size_t comp=0;comp<nb_comps[field];comp++) {
//      // Get local min and max indices
//      for (size_t l=0;l<dim;++l) {
// //        min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l );
// //        max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l );
//         min_index(l) = 0;
//         max_index(l) = FF->get_local_nb_dof( comp, l );
//      }
//
//      size_t local_min_k = 0;
//      size_t local_max_k = 1;
//
//      if (dim == 3) {
//         local_min_k = min_index(2);
//         local_max_k = max_index(2);
//      }
//
//      for (i=min_index(0);i<max_index(0);++i) {
//         double xC = FF->get_DOF_coordinate( i, comp, 0 ) ;
//         for (j=min_index(1);j<max_index(1);++j) {
//            double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;
//            for (k=local_min_k;k<local_max_k;++k) {
//               double zC = 0.;
//               if (dim == 3) zC = FF->get_DOF_coordinate( k, comp, 2 ) ;
//               size_t p = return_node_index(FF,comp,i,j,k);
//               double id = node.parID[comp]->item(p);
//               double voidf = node.void_frac[comp]->item(p);
//
// //              double lambdax = divergence[0].lambda[0]->item(p);
// //              double lambday = divergence[0].lambda[1]->item(p);
// //	      double div = divergence[0].div->item(p);
//               outputFile << xC << "," << yC << "," << zC << "," << id << "," << voidf;
// //              outputFile << xC << "," << yC << "," << zC << "," << p << "," << lambdax << "," << lambday << "," << voidf << "," << fresh[0].flag->item(p) << "," << fresh[0].neigh[0]->item(p) << "," << fresh[0].neigh[1]->item(p) << "," << fresh[0].flag_count->item(p) << "," << fresh[0].neigh_count->item(p) << "," << div << "," << fresh[0].sep_vel->item(p) << endl;
//               for (size_t dir = 0; dir < dim; dir++) {
//                   for (size_t off = 0; off < 2; off++) {
//                       outputFile << "," << b_intersect[dir].offset[comp]->item(p,off) << "," << b_intersect[dir].value[comp]->item(p,off);
//                   }
//               }
//               outputFile << endl;
//            }
//         }
//      }
//   }
//   outputFile.close();
// }
//
// //----------------------------------------------------------------------
// double
// DDS_NavierStokes::get_velocity_divergence(FV_TimeIterator const* t_it)
// //----------------------------------------------------------------------
// {
//   size_t_vector min_unknown_index(dim,0);
//   size_t_vector max_unknown_index(dim,0);
//   double div_velocity = 0.;
//   double max_divu=0.;
//
//   DivNode* divergence = GLOBAL_EQ->get_node_divergence();
//
//   for (size_t l=0;l<dim;++l) {
//     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
//     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
//   }
//
//   ofstream outputFile ;
// /*
//   std::ostringstream os2;
//   os2 << "./DS_results/divergence_" << t_it->iteration_number() << ".csv";
//   std::string filename = os2.str();
//   outputFile.open(filename.c_str());
//   outputFile << "i,j,div" << endl;
// */
//   for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
//      for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
//         double dx = PF->get_cell_size( i, 0, 0 );
//         double dy = PF->get_cell_size( j, 0, 1 );
//         if (dim == 2) {
//            size_t p = return_node_index(PF,0,i,j,0);
//            double vel_div = divergence[0].div->item(p);
//            max_divu = MAC::max( MAC::abs(vel_div), max_divu );
//            div_velocity += vel_div * ( dx * dy );
//
// //           outputFile << PF->get_DOF_coordinate( i, 0, 0 ) << "," << PF->get_DOF_coordinate( j, 0, 1 ) << "," << vel_div << endl;
//         } else {
//            for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
//               double dz = PF->get_cell_size( k, 0, 2 );
//               size_t p = return_node_index(PF,0,i,j,k);
//               double vel_div = divergence[0].div->item(p);
//               max_divu = MAC::max( MAC::abs(vel_div), max_divu );
//               div_velocity += vel_div * ( dx * dy * dz );
//            }
//         }
//      }
//   }
//
//   div_velocity = macCOMM->sum( div_velocity ) ;
// //  div_velocity = MAC::sqrt( div_velocity );
//   max_divu = macCOMM->max( max_divu ) ;
//   if ( my_rank == is_master )
//     MAC::out() << "Norm L2 div(u) = "<< MAC::doubleToString( ios::scientific, 12, div_velocity ) << " Max div(u) = " << MAC::doubleToString( ios::scientific, 12, max_divu ) << endl;
//
// //  outputFile.close( ) ;
//   return(max_divu);
// }
//
//
//
//
//----------------------------------------------------------------------
void
DDS_NavierStokes::output_L2norm_pressure( size_t const& level )
//----------------------------------------------------------------------
{
  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  double L2normP = 0.;
  double cell_P=0.,max_P=0.;

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(PF);

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
           double dx = PF->get_cell_size( i, 0, 0 );
           double dy = PF->get_cell_size( j, 0, 1 );
           double dz = (dim == 3) ? PF->get_cell_size( k, 0, 2 ) : 1;
           cell_P = PF->DOF_value( i, j, k, 0, level );
           max_P = MAC::max( MAC::abs(cell_P), max_P );
           if (is_solids) {
              size_t p = PF->DOF_local_number(i,j,k,0);
              if (void_frac->operator()(p) == 0) {
                 L2normP += cell_P * cell_P * dx * dy * dz;
              }
           } else {
              L2normP += cell_P * cell_P * dx * dy * dz;
           }
        }
     }
  }

  L2normP = macCOMM->sum( L2normP ) ;
  L2normP = MAC::sqrt( L2normP );
  max_P = macCOMM->max( max_P ) ;
  if ( my_rank == is_master )
     MAC::out() << "Norm L2 P = "
                << MAC::doubleToString( ios::scientific, 12, L2normP )
                << " Max P = "
                << MAC::doubleToString( ios::scientific, 12, max_P ) << endl;

}




//----------------------------------------------------------------------
void
DDS_NavierStokes::output_L2norm_velocity( size_t const& level )
//----------------------------------------------------------------------
{
  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);

  for (size_t comp = 0; comp < nb_comps[1]; comp++) {
     for (size_t l = 0; l < dim; ++l) {
        min_unknown_index(l) =
                     UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) =
                     UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     double L2normU = 0.;
     double cell_U=0.,max_U=0.;

     for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
        for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
           for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
              double dx = UF->get_cell_size( i,comp, 0 );
              double dy = UF->get_cell_size( j,comp, 1 );
              double dz = (dim == 3) ? UF->get_cell_size( k,comp, 2 ) : 1 ;
              cell_U = UF->DOF_value( i, j, k, comp, level );
              max_U = MAC::max( MAC::abs(cell_U), max_U );
              if (is_solids) {
                 size_t p = UF->DOF_local_number(i,j,k,comp);
                 if (void_frac->operator()(p) == 0) {
                    L2normU += cell_U * cell_U * dx * dy * dz;
                 }
              } else {
                 L2normU += cell_U * cell_U * dx * dy * dz;
              }
           }
        }
     }

     L2normU = macCOMM->sum( L2normU ) ;
     L2normU = MAC::sqrt( L2normU );
     max_U = macCOMM->max( max_U ) ;
     if ( my_rank == is_master )
        MAC::out() << "Component: " << comp
                   << " Norm L2 U = "
                   << MAC::doubleToString( ios::scientific, 12, L2normU )
                   << " Max U = "
                   << MAC::doubleToString( ios::scientific, 12, max_U ) << endl;
  }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: create_DDS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: create_DDS_subcommunicators" ) ;

   int color = 0, key = 0;
   int const* MPI_coordinates_world =
                              UF->primary_grid()->get_MPI_coordinates() ;
   int const* MPI_number_of_coordinates =
                              UF->primary_grid()->get_domain_decomposition() ;

   if (dim == 2) {
      // Assign color and key for splitting in x
      color = MPI_coordinates_world[1];
      key = MPI_coordinates_world[0];
      // Split by direction in x
      processor_splitting (color,key,0);

      // Assign color and key for splitting in y
      color = MPI_coordinates_world[0];
      key = MPI_coordinates_world[1];
      // Split by direction in y
      processor_splitting (color,key,1);
   } else {
      // Assign color and key for splitting in x
      color = MPI_coordinates_world[1]
            + MPI_coordinates_world[2]*MPI_number_of_coordinates[1] ;
      key = MPI_coordinates_world[0];
      // Split by direction in x
      processor_splitting (color,key,0);

      // Assign color and key for splitting in y
      color = MPI_coordinates_world[2]
            + MPI_coordinates_world[0]*MPI_number_of_coordinates[2];
      key = MPI_coordinates_world[1];
      // Split by direction in y
      processor_splitting (color,key,1);

      // Assign color and key for splitting in y
      color = MPI_coordinates_world[0]
            + MPI_coordinates_world[1]*MPI_number_of_coordinates[0];;
      key = MPI_coordinates_world[2];

      // Split by direction in y
      processor_splitting (color,key,2);
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: processor_splitting ( int const& color
                                       , int const& key
                                       , size_t const& dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: processor_splitting" ) ;

   MPI_Comm_split(MPI_COMM_WORLD, color, key, &DDS_Comm_i[dir]);
   MPI_Comm_size( DDS_Comm_i[dir], &nb_ranks_comm_i[dir] ) ;
   MPI_Comm_rank( DDS_Comm_i[dir], &rank_in_i[dir] ) ;

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: allocate_mpi_variables (FV_DiscreteField const* FF)
//---------------------------------------------------------------------------
{

   size_t field = (FF == PF) ? 0 : 1 ;

   for (size_t dir = 0; dir < dim; dir++) {
      first_pass[field][dir].size = new size_t [nb_comps[field]];
      second_pass[field][dir].size = new size_t [nb_comps[field]];
      for (size_t comp = 0; comp < nb_comps[field]; comp++) {
         size_t local_min_j=0, local_max_j=0;
         size_t local_min_k=0, local_max_k=0;

         // Get local min and max indices
         size_t_vector min_unknown_index(dim,0);
         size_t_vector max_unknown_index(dim,0);
         for (size_t l=0;l<dim;++l) {
            min_unknown_index(l) =
                        FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
            max_unknown_index(l) =
                        FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
         }

         if (dir == 0) {
            local_min_j = min_unknown_index(1);
            local_max_j = max_unknown_index(1);
            if (dim == 3) {
               local_min_k = min_unknown_index(2);
               local_max_k = max_unknown_index(2);
            }
         } else if (dir == 1) {
            local_min_j = min_unknown_index(0);
            local_max_j = max_unknown_index(0);
            if (dim == 3) {
               local_min_k = min_unknown_index(2);
               local_max_k = max_unknown_index(2);
            }
         } else if (dir == 2) {
            local_min_j = min_unknown_index(0);
            local_max_j = max_unknown_index(0);
            local_min_k = min_unknown_index(1);
            local_max_k = max_unknown_index(1);
         }

         size_t local_length_j = (local_max_j-local_min_j+1);
         size_t local_length_k = (local_max_k-local_min_k+1);

         if (dim != 3) {
            first_pass[field][dir].size[comp] = 3*local_length_j;
            second_pass[field][dir].size[comp] = 2*local_length_j;
         } else if (dim == 3) {
            first_pass[field][dir].size[comp] =
                                    3*local_length_j*local_length_k;
            second_pass[field][dir].size[comp] =
                                    2*local_length_j*local_length_k;
         }
      }
   }

   // Array declarations
   for (size_t dir = 0; dir < dim; dir++) {
      first_pass[field][dir].send =
                     new double** [nb_comps[field]];
      first_pass[field][dir].receive =
                     new double** [nb_comps[field]];
      second_pass[field][dir].send =
                     new double** [nb_comps[field]];
      second_pass[field][dir].receive =
                     new double** [nb_comps[field]];
      for (size_t comp = 0; comp < nb_comps[field]; comp++) {
         first_pass[field][dir].send[comp] =
                        new double* [nb_ranks_comm_i[dir]];
         first_pass[field][dir].receive[comp] =
                        new double* [nb_ranks_comm_i[dir]];
         second_pass[field][dir].send[comp] =
                        new double* [nb_ranks_comm_i[dir]];
         second_pass[field][dir].receive[comp] =
                        new double* [nb_ranks_comm_i[dir]];
         for (size_t i = 0; i < (size_t) nb_ranks_comm_i[dir]; i++) {
            first_pass[field][dir].send[comp][i] =
                           new double[first_pass[field][dir].size[comp]];
            first_pass[field][dir].receive[comp][i] =
                           new double[first_pass[field][dir].size[comp]];
            second_pass[field][dir].send[comp][i] =
                           new double[second_pass[field][dir].size[comp]];
            second_pass[field][dir].receive[comp][i] =
                           new double[second_pass[field][dir].size[comp]];
         }
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: deallocate_mpi_variables ()
//---------------------------------------------------------------------------
{
   // Array declarations
   for (size_t field = 0; field < 2; field++) {
      for (size_t dir = 0; dir < dim; dir++) {
         for (size_t comp = 0; comp < nb_comps[field]; comp++) {
            for (size_t i = 0; i < (size_t) nb_ranks_comm_i[dir]; i++) {
               delete [] first_pass[field][dir].send[comp][i];
               delete [] first_pass[field][dir].receive[comp][i];
               delete [] second_pass[field][dir].send[comp][i];
               delete [] second_pass[field][dir].receive[comp][i];
            }
            delete [] first_pass[field][dir].send[comp];
            delete [] first_pass[field][dir].receive[comp];
            delete [] second_pass[field][dir].send[comp];
            delete [] second_pass[field][dir].receive[comp];
         }
         delete [] first_pass[field][dir].send;
         delete [] first_pass[field][dir].receive;
         delete [] second_pass[field][dir].send;
         delete [] second_pass[field][dir].receive;
         delete [] first_pass[field][dir].size;
         delete [] second_pass[field][dir].size;
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: free_DDS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: free_DDS_subcommunicators" ) ;
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: set_translation_vector()
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: set_translation_vector" ) ;

   MVQ_translation_vector.resize( dim );
   translation_direction = UF->primary_grid()->get_translation_direction() ;
   MVQ_translation_vector( translation_direction ) =
        UF->primary_grid()->get_translation_magnitude() ;
}


//---------------------------------------------------------------------------
void
DDS_NavierStokes::build_links_translation()
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: build_links_translation" ) ;

   UF->create_transproj_interpolation() ;
   PF->create_transproj_interpolation() ;

}



//---------------------------------------------------------------------------
void
DDS_NavierStokes::fields_projection()
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: fields_projection" ) ;

   FV_Mesh* pmesh = const_cast<FV_Mesh*>(UF->primary_grid()) ;

   pmesh->translation() ;

   UF->translation_projection( 0, 5, 0 ) ;
   synchronize_velocity_field ( 0 ) ;

   UF->translation_projection( 1, 5, 0 ) ;
   synchronize_velocity_field ( 1 ) ;

   UF->translation_projection( 2, 5, 0 ) ;
   synchronize_velocity_field ( 2 ) ;

   UF->translation_projection( 3, 5, 0 ) ;
   synchronize_velocity_field ( 3 ) ;

   UF->translation_projection( 4, 5 ) ;
   synchronize_velocity_field ( 4 ) ;

   PF->translation_projection( 0, 2, 0 ) ;
   synchronize_pressure_field( 0 ) ;

   PF->translation_projection( 1, 2 ) ;
   synchronize_pressure_field( 1 ) ;

}


//---------------------------------------------------------------------------
void
DDS_NavierStokes:: synchronize_velocity_field ( size_t level )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: synchronize_velocity_field" ) ;

   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);

   for (size_t comp=0;comp<nb_comps[1];++comp) {

      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) =
                     UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) =
                     UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
       for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
          GLOBAL_EQ->update_global_U_vector(i,j,k,comp,
                                    UF->DOF_value(i, j, k, comp, level));
        }
       }
      }
   }

   // Synchronize the distributed DS solution vector
   GLOBAL_EQ->synchronize_DS_solution_vec();
   // Tranfer back to field
   UF->update_free_DOFs_value( level, GLOBAL_EQ->get_solution_DS_velocity() ) ;

}


//---------------------------------------------------------------------------
void
DDS_NavierStokes:: synchronize_pressure_field ( size_t level )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: synchronize_pressure_field" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  // Get local min and max indices
  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
   for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
    for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
       GLOBAL_EQ->update_global_P_vector(i,j,k,
                                 PF->DOF_value(i, j, k, 0, level));
    }
   }
  }

  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec_P();
  // Tranfer back to field
  PF->update_free_DOFs_value( level, GLOBAL_EQ->get_solution_DS_pressure() ) ;

}


//----------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_advection_Centered( size_t const& advecting_level,
                                                double const& coef,
                                                size_t const& advected_level,
                                                size_t const& i,
                                                size_t const& j,
                                                size_t const& k,
                                                size_t const& component )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: assemble_advection_Centered" );

   // Parameters
   double dxC = 0., dyC = 0., dzC = 0.;
   double AdvectedValueC = 0., AdvectedValueRi = 0., AdvectedValueLe = 0.,
    AdvectedValueTo = 0., AdvectedValueBo = 0., AdvectedValueFr = 0.,
    AdvectedValueBe = 0, AdvectorValueC = 0., AdvectorValueRi = 0.,
    AdvectorValueLe = 0., AdvectorValueTo = 0., AdvectorValueBo = 0.,
    AdvectorValueFr = 0., AdvectorValueBe = 0, AdvectorValueToLe = 0.,
    AdvectorValueToRi = 0., AdvectorValueBoLe = 0., AdvectorValueBoRi = 0.,
    AdvectorValueFrLe = 0., AdvectorValueFrRi = 0., AdvectorValueBeLe = 0.,
    AdvectorValueBeRi = 0., AdvectorValueFrTo = 0., AdvectorValueFrBo = 0.,
    AdvectorValueBeTo = 0., AdvectorValueBeBo = 0.,
    ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
    fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;

   // Comment: staggered unknowns always have a defined value at +1/-1
   // indices in directions different from their component number,
   // i.e. u in x, v in y and w in z.
   // For instance, if u on the right or left boundary is an unknown with
   // homogeneous Neumann BC, then the flux on the right or left needs special
   // treatment using the center value.
   // Otherwise, whether one of the +1/-1 DOF values is on a
   // boundary or not, and whether that boundary has a Dirichlet or Neumann
   // condition is irrelevant, this +1/-1 DOF always has the right value.
   // For Neumann, this is guaranted by
   // FV_BoundaryCondition:: set_free_DOF_values in
   // FV_DiscreteField:: update_free_DOFs_value or
   // FV_DiscreteField:: add_to_free_DOFs_value
   FV_SHIFT_TRIPLET shift = UF->shift_staggeredToStaggered( component );

   dxC = UF->get_cell_size(i,component,0) ;
   dyC = UF->get_cell_size(j,component,1) ;
   if (dim == 3) dzC = UF->get_cell_size(k,component,2) ;

   AdvectedValueC = UF->DOF_value( i, j, k, component, advected_level );
   AdvectorValueC = UF->DOF_value( i, j, k, component, advecting_level );

   // The First Component (u)
   if ( component == 0 ) {
      // Right (U_X)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
         fri = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueRi = UF->DOF_value(i+1, j, k, component, advected_level );
         AdvectorValueRi = UF->DOF_value(i+1, j, k, component, advecting_level );
         ur = 0.5 * ( AdvectorValueC + AdvectorValueRi );
         fri = ur * ur;
      }

      // Left (U_X)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
         fle = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueLe = UF->DOF_value(i-1, j, k, component, advected_level );
         AdvectorValueLe = UF->DOF_value(i-1, j, k, component, advecting_level );
         ul = 0.5 * ( AdvectorValueC + AdvectorValueLe );
         fle = ul * ul;
      }

      // Top (U_Y)
      AdvectedValueTo = UF->DOF_value(i, j+1, k, component, advected_level );
      AdvectorValueToLe = UF->DOF_value(i+shift.i-1, j+shift.j, k, 1, advecting_level );
      AdvectorValueToRi = UF->DOF_value(i+shift.i, j+shift.j, k, 1, advecting_level );
      vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
      fto = vt * 0.5 * ( AdvectedValueC + AdvectedValueTo );

      // Bottom (U_Y)
      AdvectedValueBo = UF->DOF_value(i, j-1, k, component, advected_level );
      AdvectorValueBoLe = UF->DOF_value(i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
      AdvectorValueBoRi = UF->DOF_value(i+shift.i, j+shift.j-1, k, 1, advecting_level );
      vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
      fbo = vb * 0.5 * ( AdvectedValueBo + AdvectedValueC );

      if (dim == 3) {
         // Front (U_Z)
         AdvectedValueFr = UF->DOF_value(i, j, k+1, component, advected_level );
         AdvectorValueFrLe = UF->DOF_value(i+shift.i-1, j, k+shift.k, 2, advecting_level );
         AdvectorValueFrRi = UF->DOF_value(i+shift.i, j, k+shift.k, 2, advecting_level );
         wf = 0.5 * ( AdvectorValueFrLe + AdvectorValueFrRi );
         ffr = wf * 0.5 * ( AdvectedValueC + AdvectedValueFr);

         // Behind (U_Z)
         AdvectedValueBe = UF->DOF_value(i, j, k-1, component, advected_level );
         AdvectorValueBeLe = UF->DOF_value(i+shift.i-1, j, k+shift.k-1, 2, advecting_level );
         AdvectorValueBeRi = UF->DOF_value(i+shift.i, j, k+shift.k-1, 2, advecting_level );
         wb = 0.5 * ( AdvectorValueBeLe + AdvectorValueBeRi );
         fbe = wb * 0.5 * ( AdvectedValueBe + AdvectedValueC );
      }
   } else if (component == 1) {
      // The second Component (v)
      // Right (V_X)
      AdvectedValueRi = UF->DOF_value(i+1, j, k, component, advected_level );
      AdvectorValueToRi = UF->DOF_value(i+shift.i, j+shift.j, k, 0, advecting_level );
      AdvectorValueBoRi = UF->DOF_value(i+shift.i, j+shift.j-1, k, 0, advecting_level );
      ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
      fri = ur * 0.5 * ( AdvectedValueC + AdvectedValueRi );

      // Left (V_X)
      AdvectedValueLe = UF->DOF_value(i-1, j, k, component, advected_level );
      AdvectorValueToLe = UF->DOF_value(i+shift.i-1, j+shift.j, k, 0, advecting_level );
      AdvectorValueBoLe = UF->DOF_value(i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
      ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
      fle = ul * 0.5 * ( AdvectedValueLe + AdvectedValueC );

      // Top (V_Y)
      if ( UF->DOF_color(i, j, k, component ) == FV_BC_TOP )
         fto = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueTo = UF->DOF_value(i, j+1, k, component, advected_level );
         AdvectorValueTo = UF->DOF_value(i, j+1, k, component, advecting_level );
         vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
         fto = vt * 0.5 * ( AdvectedValueC + AdvectedValueTo );
      }

      // Bottom (V_Y)
      if ( UF->DOF_color(i, j, k, component ) == FV_BC_BOTTOM )
         fbo = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueBo = UF->DOF_value(i, j-1, k, component, advected_level );
         AdvectorValueBo = UF->DOF_value(i, j-1, k, component, advecting_level );
         vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
         fbo = vb * 0.5 * ( AdvectedValueBo + AdvectedValueC );
      }

      if (dim == 3) {
         // Front (V_Z)
         AdvectedValueFr = UF->DOF_value(i, j, k+1, component, advected_level );
         AdvectorValueFrTo = UF->DOF_value(i, j+shift.j, k+shift.k, 2, advecting_level );
         AdvectorValueFrBo = UF->DOF_value(i, j+shift.j-1, k+shift.k, 2, advecting_level );
         wf = 0.5 * ( AdvectorValueFrTo + AdvectorValueFrBo );
         ffr = wf * 0.5 * ( AdvectedValueC + AdvectedValueFr );

         // Behind (V_Z)
         AdvectedValueBe = UF->DOF_value(i, j, k-1, component, advected_level );
         AdvectorValueBeTo = UF->DOF_value(i, j+shift.j, k+shift.k-1, 2, advecting_level );
         AdvectorValueBeBo = UF->DOF_value(i, j+shift.j-1, k+shift.k-1, 2, advecting_level );
         wb = 0.5 * ( AdvectorValueBeTo + AdvectorValueBeBo );
         fbe = wb * 0.5 * ( AdvectedValueBe + AdvectedValueC );
      }
   } else {
      // The Third Component (w)
      // Right (W_X)
      AdvectedValueRi = UF->DOF_value(i+1, j, k, component, advected_level );
      AdvectorValueFrRi = UF->DOF_value(i+shift.i, j, k+shift.k, 0, advecting_level );
      AdvectorValueBeRi = UF->DOF_value(i+shift.i, j, k+shift.k-1, 0, advecting_level );
      ur = 0.5 * ( AdvectorValueFrRi + AdvectorValueBeRi );
      fri = ur * 0.5 * ( AdvectedValueC + AdvectedValueRi );

      // Left (W_X)
      AdvectedValueLe = UF->DOF_value(i-1, j, k, component, advected_level );
      AdvectorValueFrLe = UF->DOF_value(i+shift.i-1, j, k+shift.k, 0, advecting_level );
      AdvectorValueBeLe = UF->DOF_value(i+shift.i-1, j, k+shift.k-1, 0, advecting_level );
      ul = 0.5 * ( AdvectorValueFrLe + AdvectorValueBeLe );
      fle = ul * 0.5 * ( AdvectedValueLe + AdvectedValueC );

      // Top (W_Y)
      AdvectedValueTo = UF->DOF_value(i, j+1, k, component, advected_level );
      AdvectorValueFrTo = UF->DOF_value(i, j+shift.j, k+shift.k, 1, advecting_level );
      AdvectorValueBeTo = UF->DOF_value(i, j+shift.j, k+shift.k-1, 1, advecting_level );
      vt = 0.5 * ( AdvectorValueFrTo + AdvectorValueBeTo );
      fto = vt * 0.5 * ( AdvectedValueC + AdvectedValueTo );

      // Bottom (W_Y)
      AdvectedValueBo = UF->DOF_value(i, j-1, k, component, advected_level );
      AdvectorValueFrBo = UF->DOF_value(i, j+shift.j-1, k+shift.k, 1, advecting_level );
      AdvectorValueBeBo = UF->DOF_value(i, j+shift.j-1, k+shift.k-1, 1, advecting_level );
      vb = 0.5 * ( AdvectorValueFrBo + AdvectorValueBeBo );
      fbo = vb * 0.5 * ( AdvectedValueBo + AdvectedValueC );

      // Front (W_Z)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
         ffr = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueFr = UF->DOF_value(i, j, k+1, component, advected_level );
         AdvectorValueFr = UF->DOF_value(i, j, k+1, component, advecting_level );
         wf = 0.5 * ( AdvectorValueFr + AdvectorValueC );
         ffr = wf * 0.5 * ( AdvectedValueC + AdvectedValueFr );
      }

      // Behind (W_Z)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND )
         fbe = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueBe = UF->DOF_value(i, j, k-1, component, advected_level );
         AdvectorValueBe = UF->DOF_value(i, j, k-1, component, advecting_level );
         wb = 0.5 * ( AdvectorValueBe + AdvectorValueC );
         fbe = wb * 0.5 * ( AdvectedValueBe + AdvectedValueC );
      }
   }

   if (dim == 2) {
      flux = ((fto - fbo) * dxC + (fri - fle) * dyC);
   } else if (dim == 3) {
      flux = (fto - fbo) * dxC * dzC + (fri - fle) * dyC * dzC + (ffr - fbe) * dxC * dyC;
   }
   return ( coef * flux );
}



//----------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_advection_Upwind(
  size_t const& advecting_level, double const& coef, size_t const& advected_level,
  size_t const& i, size_t const& j, size_t const& k, size_t const& component )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: assemble_advection_Upwind" );

   // Parameters
   double dxC = 0., dyC = 0., dzC = 0.;
   double AdvectedValueC = 0., AdvectedValueRi = 0., AdvectedValueLe = 0.,
    AdvectedValueTo = 0., AdvectedValueBo = 0., AdvectedValueFr = 0.,
    AdvectedValueBe = 0, AdvectorValueC = 0., AdvectorValueRi = 0.,
    AdvectorValueLe = 0., AdvectorValueTo = 0., AdvectorValueBo = 0.,
    AdvectorValueFr = 0., AdvectorValueBe = 0, AdvectorValueToLe = 0.,
    AdvectorValueToRi = 0., AdvectorValueBoLe = 0., AdvectorValueBoRi = 0.,
    AdvectorValueFrLe = 0., AdvectorValueFrRi = 0., AdvectorValueBeLe = 0.,
    AdvectorValueBeRi = 0., AdvectorValueFrTo = 0., AdvectorValueFrBo = 0.,
    AdvectorValueBeTo = 0., AdvectorValueBeBo = 0.,
    ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
    fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;

   // Comment: staggered unknowns always have a defined value at +1/-1
   // indices in directions different from their component number,
   // i.e. u in x, v in y and w in z.
   // For instance, if u on the right or left boundary is an unknown with
   // homogeneous Neumann BC, then the flux on the right or left needs special
   // treatment using the center value.
   // Otherwise, whether one of the +1/-1 DOF values is on a
   // boundary or not, and whether that boundary has a Dirichlet or Neumann
   // condition is irrelevant, this +1/-1 DOF always has the right value.
   // For Neumann, this is guaranted by
   // FV_BoundaryCondition:: set_free_DOF_values in
   // FV_DiscreteField:: update_free_DOFs_value or
   // FV_DiscreteField:: add_to_free_DOFs_value
   FV_SHIFT_TRIPLET shift = UF->shift_staggeredToStaggered( component );

   dxC = UF->get_cell_size(i,component,0) ;
   dyC = UF->get_cell_size(j,component,1) ;
   if (dim == 3) dzC = UF->get_cell_size(k,component,2) ;

   AdvectedValueC = UF->DOF_value( i, j, k, component, advected_level );
   AdvectorValueC = UF->DOF_value( i, j, k, component, advecting_level );

   // The First Component (u)
   if ( component == 0 ) {
      // Right (U_X)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
         fri = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueRi = UF->DOF_value(i+1, j, k, component, advected_level );
         AdvectorValueRi = UF->DOF_value(i+1, j, k, component, advecting_level );
         ur = 0.5 * ( AdvectorValueC + AdvectorValueRi );
         if ( ur > 0. ) fri = ur * AdvectedValueC;
         else fri = ur * AdvectedValueRi;
      }

      // Left (U_X)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
         fle = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueLe = UF->DOF_value(i-1, j, k, component, advected_level );
         AdvectorValueLe = UF->DOF_value(i-1, j, k, component, advecting_level );
         ul = 0.5 * ( AdvectorValueC + AdvectorValueLe );
         if ( ul > 0. ) fle = ul * AdvectedValueLe;
         else fle = ul * AdvectedValueC;
      }

      // Top (U_Y)
      AdvectedValueTo = UF->DOF_value(i, j+1, k, component, advected_level );
      AdvectorValueToLe = UF->DOF_value(i+shift.i-1, j+shift.j, k, 1, advecting_level );
      AdvectorValueToRi = UF->DOF_value(i+shift.i, j+shift.j, k, 1, advecting_level );
      vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
      if ( vt > 0. ) fto = vt * AdvectedValueC;
      else fto = vt * AdvectedValueTo;

      // Bottom (U_Y)
      AdvectedValueBo = UF->DOF_value(i, j-1, k, component, advected_level );
      AdvectorValueBoLe = UF->DOF_value(i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
      AdvectorValueBoRi = UF->DOF_value(i+shift.i, j+shift.j-1, k, 1, advecting_level );
      vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
      if ( vb > 0. ) fbo = vb * AdvectedValueBo;
      else fbo = vb * AdvectedValueC;

      if (dim == 3) {
         // Front (U_Z)
         AdvectedValueFr = UF->DOF_value(i, j, k+1, component, advected_level );
         AdvectorValueFrLe = UF->DOF_value(i+shift.i-1, j, k+shift.k, 2, advecting_level );
         AdvectorValueFrRi = UF->DOF_value(i+shift.i, j, k+shift.k, 2, advecting_level );
         wf = 0.5 * ( AdvectorValueFrLe + AdvectorValueFrRi );
         if ( wf > 0. ) ffr = wf * AdvectedValueC;
         else ffr = wf * AdvectedValueFr;

         // Behind (U_Z)
         AdvectedValueBe = UF->DOF_value(i, j, k-1, component, advected_level );
         AdvectorValueBeLe = UF->DOF_value(i+shift.i-1, j, k+shift.k-1, 2, advecting_level );
         AdvectorValueBeRi = UF->DOF_value(i+shift.i, j, k+shift.k-1, 2, advecting_level );
         wb = 0.5 * ( AdvectorValueBeLe + AdvectorValueBeRi );
         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
         else fbe = wb * AdvectedValueC;
      }
   } else if (component == 1) {
      // The second Component (v)
      // Right (V_X)
      AdvectedValueRi = UF->DOF_value(i+1, j, k, component, advected_level );
      AdvectorValueToRi = UF->DOF_value(i+shift.i, j+shift.j, k, 0, advecting_level );
      AdvectorValueBoRi = UF->DOF_value(i+shift.i, j+shift.j-1, k, 0, advecting_level );
      ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
      if ( ur > 0. ) fri = ur * AdvectedValueC;
      else fri = ur * AdvectedValueRi;

      // Left (V_X)
      AdvectedValueLe = UF->DOF_value(i-1, j, k, component, advected_level );
      AdvectorValueToLe = UF->DOF_value(i+shift.i-1, j+shift.j, k, 0, advecting_level );
      AdvectorValueBoLe = UF->DOF_value(i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
      ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
      if ( ul > 0. ) fle = ul * AdvectedValueLe;
      else fle = ul * AdvectedValueC;

      // Top (V_Y)
      if ( UF->DOF_color(i, j, k, component ) == FV_BC_TOP )
         fto = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueTo = UF->DOF_value(i, j+1, k, component, advected_level );
         AdvectorValueTo = UF->DOF_value(i, j+1, k, component, advecting_level );
         vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
         if ( vt > 0. ) fto = vt * AdvectedValueC;
         else fto = vt * AdvectedValueTo;
      }

      // Bottom (V_Y)
      if ( UF->DOF_color(i, j, k, component ) == FV_BC_BOTTOM )
         fbo = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueBo = UF->DOF_value(i, j-1, k, component, advected_level );
         AdvectorValueBo = UF->DOF_value(i, j-1, k, component, advecting_level );
         vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
         else fbo = vb * AdvectedValueC;
      }

      if (dim == 3) {
         // Front (V_Z)
         AdvectedValueFr = UF->DOF_value(i, j, k+1, component, advected_level );
         AdvectorValueFrTo = UF->DOF_value(i, j+shift.j, k+shift.k, 2, advecting_level );
         AdvectorValueFrBo = UF->DOF_value(i, j+shift.j-1, k+shift.k, 2, advecting_level );
         wf = 0.5 * ( AdvectorValueFrTo + AdvectorValueFrBo );
         if ( wf > 0. ) ffr = wf * AdvectedValueC;
         else ffr = wf * AdvectedValueFr;

         // Behind (V_Z)
         AdvectedValueBe = UF->DOF_value(i, j, k-1, component, advected_level );
         AdvectorValueBeTo = UF->DOF_value(i, j+shift.j, k+shift.k-1, 2, advecting_level );
         AdvectorValueBeBo = UF->DOF_value(i, j+shift.j-1, k+shift.k-1, 2, advecting_level );
         wb = 0.5 * ( AdvectorValueBeTo + AdvectorValueBeBo );
         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
         else fbe = wb * AdvectedValueC;
      }
   } else {
      // The Third Component (w)
      // Right (W_X)
      AdvectedValueRi = UF->DOF_value(i+1, j, k, component, advected_level );
      AdvectorValueFrRi = UF->DOF_value(i+shift.i, j, k+shift.k, 0, advecting_level );
      AdvectorValueBeRi = UF->DOF_value(i+shift.i, j, k+shift.k-1, 0, advecting_level );
      ur = 0.5 * ( AdvectorValueFrRi + AdvectorValueBeRi );
      if ( ur > 0. ) fri = ur * AdvectedValueC;
      else fri = ur * AdvectedValueRi;

      // Left (W_X)
      AdvectedValueLe = UF->DOF_value(i-1, j, k, component, advected_level );
      AdvectorValueFrLe = UF->DOF_value(i+shift.i-1, j, k+shift.k, 0, advecting_level );
      AdvectorValueBeLe = UF->DOF_value(i+shift.i-1, j, k+shift.k-1, 0, advecting_level );
      ul = 0.5 * ( AdvectorValueFrLe + AdvectorValueBeLe );
      if ( ul > 0. ) fle = ul * AdvectedValueLe;
      else fle = ul * AdvectedValueC;

      // Top (W_Y)
      AdvectedValueTo = UF->DOF_value(i, j+1, k, component, advected_level );
      AdvectorValueFrTo = UF->DOF_value(i, j+shift.j, k+shift.k, 1, advecting_level );
      AdvectorValueBeTo = UF->DOF_value(i, j+shift.j, k+shift.k-1, 1, advecting_level );
      vt = 0.5 * ( AdvectorValueFrTo + AdvectorValueBeTo );
      if ( vt > 0. ) fto = vt * AdvectedValueC;
      else fto = vt * AdvectedValueTo;

      // Bottom (W_Y)
      AdvectedValueBo = UF->DOF_value(i, j-1, k, component, advected_level );
      AdvectorValueFrBo = UF->DOF_value(i, j+shift.j-1, k+shift.k, 1, advecting_level );
      AdvectorValueBeBo = UF->DOF_value(i, j+shift.j-1, k+shift.k-1, 1, advecting_level );
      vb = 0.5 * ( AdvectorValueFrBo + AdvectorValueBeBo );
      if ( vb > 0. ) fbo = vb * AdvectedValueBo;
      else fbo = vb * AdvectedValueC;

      // Front (W_Z)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
         ffr = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueFr = UF->DOF_value(i, j, k+1, component, advected_level );
         AdvectorValueFr = UF->DOF_value(i, j, k+1, component, advecting_level );
         wf = 0.5 * ( AdvectorValueFr + AdvectorValueC );
         if ( wf > 0. ) ffr = wf * AdvectedValueC;
         else ffr = wf * AdvectedValueFr;
      }

      // Behind (W_Z)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND )
         fbe = AdvectorValueC * AdvectedValueC;
      else {
         AdvectedValueBe = UF->DOF_value(i, j, k-1, component, advected_level );
         AdvectorValueBe = UF->DOF_value(i, j, k-1, component, advecting_level );
         wb = 0.5 * ( AdvectorValueBe + AdvectorValueC );
         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
         else fbe = wb * AdvectedValueC;
      }
   }

   if (dim == 2) {
      flux = ((fto - fbo) * dxC + (fri - fle) * dyC);
   } else if (dim == 3) {
      flux = (fto - fbo) * dxC * dzC + (fri - fle) * dyC * dzC + (ffr - fbe) * dxC * dyC;
   }
   return ( coef * flux );
}

//----------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_advection_TVD(
  size_t const& advecting_level, double const& coef, size_t const& advected_level,
         size_t const& i, size_t const& j, size_t const& k, size_t const& component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: assemble_advection_TVD" );

   // Parameters
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   double xC = 0., yC = 0., zC = 0., xr = 0., xR = 0., xl = 0., xL = 0.,
    yt = 0., yT = 0., yb = 0., yB = 0.,
    zf = 0., zF = 0., zb = 0., zB = 0.;
   double dxC = 0., dyC = 0., dzC = 0., dxr = 0., dxl = 0., dxCr = 0.,
    dxCl = 0., dxRr = 0., dxR = 0., dxLl = 0., dyt = 0., dyb = 0.,
    dyCt = 0., dyCb = 0., dyTt = 0., dyT = 0., dyBb = 0., dzf = 0.,
    dzb = 0., dzCf = 0., dzCb = 0., dzFf = 0., dzF = 0., dzBb = 0.;

   double AdvectedValueC = 0., AdvectedValueRi = 0., AdvectedValueLe = 0.,
    AdvectedValueTo = 0., AdvectedValueBo = 0., AdvectedValueFr = 0.,
    AdvectedValueBe = 0, AdvectedValueLeLe=0., AdvectedValueRiRi=0.,
    AdvectedValueBoBo=0., AdvectedValueToTo=0., AdvectedValueBeBe=0.,
    AdvectedValueFrFr=0., AdvectorValueC = 0., AdvectorValueRi = 0.,
    AdvectorValueLe = 0., AdvectorValueTo = 0., AdvectorValueBo = 0.,
    AdvectorValueFr = 0., AdvectorValueBe = 0, AdvectorValueToLe = 0.,
    AdvectorValueToRi = 0., AdvectorValueBoLe = 0., AdvectorValueBoRi = 0.,
    AdvectorValueFrLe = 0., AdvectorValueFrRi = 0., AdvectorValueBeLe = 0.,
    AdvectorValueBeRi = 0., AdvectorValueFrTo = 0., AdvectorValueFrBo = 0.,
    AdvectorValueBeTo = 0., AdvectorValueBeBo = 0.;
   double ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
    fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
   double cRip12 = 0., cLip12 = 0., cRim12 = 0., cLim12 = 0., thetaC = 0.,
    thetaRi = 0., thetaLe = 0., thetaTo = 0., thetaBo = 0., thetaFr = 0.,
    thetaBe = 0.;

   FV_SHIFT_TRIPLET shift = UF->shift_staggeredToStaggered( component ) ;

   // Perform assembling
   xC = UF->get_DOF_coordinate( i, component, 0 );
   dxC = UF->get_cell_size( i, component, 0 ) ;
   yC = UF->get_DOF_coordinate( j, component, 1 );
   dyC = UF->get_cell_size( j, component, 1 ) ;
   if (dim == 3) {
      zC =UF->get_DOF_coordinate( k, component, 2 ) ;
      dzC =UF->get_cell_size( k, component, 2 ) ;
   }

   AdvectorValueC = UF->DOF_value( i, j, k, component, advecting_level );
   AdvectedValueC = UF->DOF_value( i, j, k, component, advected_level );

   // The First component (u)
   if (component == 0) {
      // Right and Left
      // --------------
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT ) {
         AdvectorValueRi = AdvectorValueC;
         AdvectedValueRi = AdvectedValueC;
      } else {
         AdvectorValueRi = UF->DOF_value( i+1, j, k, component, advecting_level );
         AdvectedValueRi =UF->DOF_value( i+1, j, k, component, advected_level );
      }

      if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT ) {
         AdvectorValueLe = AdvectorValueC;
         AdvectedValueLe = AdvectedValueC;
      } else {
         AdvectorValueLe = UF->DOF_value( i-1, j, k, component, advecting_level );
         AdvectedValueLe = UF->DOF_value( i-1, j, k, component, advected_level );
      }

      thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ?
      ( AdvectedValueC - AdvectedValueLe ) / ( AdvectedValueRi - AdvectedValueC ) : 1.e20;

      // Right (X)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
         fri = AdvectorValueC * AdvectedValueC;
      else {
         ur = 0.5 * ( AdvectorValueRi + AdvectorValueC );
         if ( UF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT ) {
            if ( ur > 0. ) fri = ur * AdvectedValueC;
            else fri = ur * AdvectedValueRi;
         } else {
            xr =UF->get_DOF_coordinate( i+shift.i, 1, 0 );
            xR =UF->get_DOF_coordinate( i+1, component, 0 );
            dxCr = xr - xC;
            dxr  = xR - xC;
            cLip12 = AdvectedValueC + ( dxCr / dxr ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                     * ( AdvectedValueRi - AdvectedValueC );

            dxRr = xR - xr;
            dxR = UF->get_cell_size( i+1, component, 0 );
            AdvectedValueRiRi =UF->DOF_value( i+2, j, k, component, advected_level );

            thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 1.e-20 ?
            ( AdvectedValueRi - AdvectedValueC ) / ( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
            cRip12 = AdvectedValueRi - ( dxRr / dxR ) * FV_DiscreteField::SuperBee_phi(thetaRi)
                                                      * ( AdvectedValueRiRi - AdvectedValueRi );
            fri = 0.5 * ( ur * ( cRip12 + cLip12 ) - fabs(ur) * ( cRip12 - cLip12 ) );
         }
      }

      // Left (X)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
         fle = AdvectorValueC * AdvectedValueC;
      else {
         ul = 0.5 * ( AdvectorValueLe + AdvectorValueC );
         if ( UF->DOF_color(i-1, j, k, component ) == FV_BC_LEFT ) {
            if ( ul > 0. ) fle = ul * AdvectedValueLe;
            else fle = ul * AdvectedValueC;
         } else {
            xl =UF->get_DOF_coordinate( i+shift.i-1, 1, 0 );
            xL =UF->get_DOF_coordinate( i-1, component, 0 );
            dxl  = xC - xL;
            dxLl = xl - xL;

            AdvectedValueLeLe =UF->DOF_value( i-2, j, k, component, advected_level );

            thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
            ( AdvectedValueLe - AdvectedValueLeLe ) / ( AdvectedValueC - AdvectedValueLe ) : 1.e20;
            cLim12 = AdvectedValueLe + ( dxLl / dxl ) * FV_DiscreteField::SuperBee_phi(thetaLe)
                                                      * ( AdvectedValueC - AdvectedValueLe );
            if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
               cRim12 = AdvectedValueC;
            else {
               xR =UF->get_DOF_coordinate( i+1, component, 0 );
               dxr  = xR - xC;
               dxCl = xC - xl;

               cRim12 = AdvectedValueC - ( dxCl / dxr ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                        * ( AdvectedValueRi - AdvectedValueC );
            }

            fle = 0.5 * ( ul * ( cRim12 + cLim12 ) - fabs(ul) * ( cRim12 - cLim12 ) );
         }
      }

      // Top and Bottom
      // --------------
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP ) {
         AdvectorValueTo = AdvectorValueC;
         AdvectedValueTo = AdvectedValueC;
      } else {
         AdvectorValueTo = UF->DOF_value( i, j+1, k, component, advecting_level );
         AdvectedValueTo =UF->DOF_value( i, j+1, k, component, advected_level );
      }

      if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM ) {
         AdvectorValueBo = AdvectorValueC;
         AdvectedValueBo = AdvectedValueC;
      } else {
         AdvectorValueBo = UF->DOF_value( i, j-1, k, component, advecting_level );
         AdvectedValueBo =UF->DOF_value( i, j-1, k, component, advected_level );
      }

      thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ?
      ( AdvectedValueC - AdvectedValueBo ) / ( AdvectedValueTo - AdvectedValueC ) : 1.e20;

      // Top (Y)
      AdvectorValueToLe = UF->DOF_value( i+shift.i-1, j+shift.j, k, 1, advecting_level );
      AdvectorValueToRi = UF->DOF_value( i+shift.i, j+shift.j, k, 1, advecting_level );
      vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
      if ( UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP
        || UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP_LEFT
        || UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP_RIGHT ) {
         if ( vt > 0. ) fto = vt * AdvectedValueC;
         else fto = vt * AdvectedValueTo;
      } else {
         yt = UF->get_DOF_coordinate( j+shift.j, 1, 1 );
         yT = UF->get_DOF_coordinate( j+1, component, 1 );
         dyCt = yt - yC;
         dyt  = yT - yC;

         cLip12 = AdvectedValueC + ( dyCt / dyt ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                   * ( AdvectedValueTo - AdvectedValueC );
         dyTt = yT - yt;
         dyT =UF->get_cell_size( j+1, component, 1 );

         AdvectedValueToTo =UF->DOF_value( i, j+2, k, component, advected_level );

         thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 1.e-20 ?
         ( AdvectedValueTo - AdvectedValueC ) / ( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;

         cRip12 = AdvectedValueTo - ( dyTt / dyT ) * FV_DiscreteField::SuperBee_phi(thetaTo)
                                                 * ( AdvectedValueToTo - AdvectedValueTo );

         fto = 0.5 * ( vt * ( cRip12 + cLip12 ) - fabs(vt) * ( cRip12 - cLip12 ) );
      }

      // Bottom (Y)
      AdvectorValueBoLe = UF->DOF_value( i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
      AdvectorValueBoRi = UF->DOF_value( i+shift.i, j+shift.j-1, k, 1, advecting_level );
      vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
      if ( UF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
        || UF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM_LEFT
        || UF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM_RIGHT ) {
         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
         else fbo = vb * AdvectedValueC;
      } else {
         yb =UF->get_DOF_coordinate( j+shift.j-1, 1, 1 );
         yB =UF->get_DOF_coordinate( j-1, component, 1 );
         dyb  = yC - yB;
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP)
            cRim12 = AdvectedValueC;
         else {
            yT = UF->get_DOF_coordinate( j+1, component, 1 );
            dyt  = yT - yC;
            dyCb = yC - yb;
            cRim12 = AdvectedValueC - ( dyCb / dyt ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                      * ( AdvectedValueTo - AdvectedValueC );
         }
         dyBb = yb - yB;
         AdvectedValueBoBo =UF->DOF_value( i, j-2, k, component, advected_level );

         thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
         ( AdvectedValueBo - AdvectedValueBoBo ) / ( AdvectedValueC - AdvectedValueBo ) : 1.e20;
         cLim12 = AdvectedValueBo + ( dyBb / dyb ) * FV_DiscreteField::SuperBee_phi(thetaBo)
                                                      * ( AdvectedValueC - AdvectedValueBo );
         fbo = 0.5 * ( vb * ( cRim12 + cLim12 ) - fabs(vb) * ( cRim12 - cLim12 ) );
      }

      if (dim == 3) {
         // Front and Behind
         // ----------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT ) {
            AdvectorValueFr = AdvectorValueC;
            AdvectedValueFr = AdvectedValueC;
         } else {
            AdvectorValueFr = UF->DOF_value( i, j, k+1, component, advecting_level );
            AdvectedValueFr =UF->DOF_value( i, j, k+1, component, advected_level );
         }

         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND ) {
            AdvectorValueBe = AdvectorValueC;
            AdvectedValueBe = AdvectedValueC;
         } else {
            AdvectorValueBe = UF->DOF_value( i, j, k-1, component, advecting_level );
            AdvectedValueBe =UF->DOF_value( i, j, k-1, component, advected_level );
         }

         thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ?
         ( AdvectedValueC - AdvectedValueBe ) / ( AdvectedValueFr - AdvectedValueC ) : 1.e20;

         // Front (Z)
         AdvectorValueFrLe = UF->DOF_value( i+shift.i-1, j, k+shift.k, 2, advecting_level );
         AdvectorValueFrRi = UF->DOF_value( i+shift.i, j, k+shift.k, 2, advecting_level );
         wf = 0.5 * (AdvectorValueFrLe + AdvectorValueFrRi);
         if ( UF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT
           || UF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT_LEFT
           || UF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT_RIGHT ) {
            if ( wf > 0. ) ffr = wf * AdvectedValueC;
            else ffr = wf * AdvectedValueFr;
         } else {
            zf =UF->get_DOF_coordinate( k+shift.k, 2, 2 );
            zF =UF->get_DOF_coordinate( k+1, component, 2 );
            dzCf = zf - zC;
            dzf  = zF - zC;
            cLip12 = AdvectedValueC + ( dzCf / dzf ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                      * ( AdvectedValueFr - AdvectedValueC );
            dzFf = zF - zf;
            dzF =UF->get_cell_size( k+1, component, 2 );
            AdvectedValueFrFr =UF->DOF_value( i, j, k+2, component, advected_level );

            thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 1.e-20 ?
            ( AdvectedValueFr - AdvectedValueC ) / ( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
            cRip12 = AdvectedValueFr - ( dzFf / dzF ) * FV_DiscreteField::SuperBee_phi(thetaFr)
                                                     * ( AdvectedValueFrFr - AdvectedValueFr );
            ffr = 0.5 * ( wf * ( cRip12 + cLip12 ) - fabs(wf) * ( cRip12 - cLip12 ) );
         }

         // Behind (Z)
         AdvectorValueBeLe = UF->DOF_value( i+shift.i-1, j, k+shift.k-1, 2, advecting_level );
         AdvectorValueBeRi = UF->DOF_value( i+shift.i, j, k+shift.k-1, 2, advecting_level );
         wb = 0.5 * ( AdvectorValueBeLe + AdvectorValueBeRi );
         if ( UF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND
           || UF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND_LEFT
           || UF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND_RIGHT ) {
            if ( wb > 0. ) fbe = wb * AdvectedValueBe;
            else fbe = wb * AdvectedValueC;
         } else {
            zb =UF->get_DOF_coordinate( k+shift.k-1, 2, 2 );
            zB =UF->get_DOF_coordinate( k-1, component, 2 );
            dzb = zC - zB;
            if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
               cRim12 = AdvectedValueC;
            else {
               zF =UF->get_DOF_coordinate( k+1, component, 2 );
               dzf  = zF - zC;
               dzCb = zC - zb;
               cRim12 = AdvectedValueC - ( dzCb / dzf ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                         * ( AdvectedValueFr - AdvectedValueC );
            }
            dzBb = zb - zB;
            AdvectedValueBeBe =UF->DOF_value( i, j, k-2, component, advected_level );

            thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
            ( AdvectedValueBe - AdvectedValueBeBe ) / ( AdvectedValueC - AdvectedValueBe ) : 1.e20;
            cLim12 = AdvectedValueBe + ( dzBb / dzb ) * FV_DiscreteField::SuperBee_phi(thetaBe)
                                                        * ( AdvectedValueC - AdvectedValueBe );
            fbe = 0.5 * ( wb * ( cRim12 + cLim12 ) - fabs(wb) * ( cRim12 - cLim12 ) );
         }
      }
   } else if (component == 1) {
      // The second component (v)
      // Right and Left
      // --------------
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT ) {
         AdvectorValueRi = AdvectorValueC;
         AdvectedValueRi = AdvectedValueC;
      } else {
         AdvectorValueRi = UF->DOF_value( i+1, j, k, component, advecting_level );
         AdvectedValueRi =UF->DOF_value( i+1, j, k, component, advected_level );
      }

      if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT ) {
         AdvectorValueLe = AdvectorValueC;
         AdvectedValueLe = AdvectedValueC;
      } else {
         AdvectorValueLe = UF->DOF_value( i-1, j, k, component, advecting_level );
         AdvectedValueLe =UF->DOF_value( i-1, j, k, component, advected_level );
      }

      thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ?
      ( AdvectedValueC - AdvectedValueLe ) / ( AdvectedValueRi - AdvectedValueC ) : 1.e20;

      // Right (X)
      AdvectorValueToRi = UF->DOF_value( i+shift.i, j+shift.j, k, 0, advecting_level );
      AdvectorValueBoRi = UF->DOF_value( i+shift.i, j+shift.j-1, k, 0, advecting_level );
      ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
      if ( UF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
        || UF->DOF_color( i+1, j, k, component ) == FV_BC_BOTTOM_RIGHT
        || UF->DOF_color( i+1, j, k, component ) == FV_BC_TOP_RIGHT ) {
         if ( ur > 0. ) fri = ur * AdvectedValueC;
         else fri = ur * AdvectedValueRi;
      } else {
         xr =UF->get_DOF_coordinate( i+shift.i, 0, 0 );
         xR =UF->get_DOF_coordinate( i+1, component, 0 );
         dxCr = xr - xC;
         dxr  = xR - xC;

         cLip12 = AdvectedValueC + ( dxCr / dxr ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                   * ( AdvectedValueRi - AdvectedValueC );

         dxRr = xR - xr;
         dxR =UF->get_cell_size( i+1, component, 0 );
         AdvectedValueRiRi =UF->DOF_value( i+2, j, k, component, advected_level );

         thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 1.e-20 ?
         ( AdvectedValueRi - AdvectedValueC ) / ( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
         cRip12 = AdvectedValueRi - ( dxRr / dxR ) * FV_DiscreteField::SuperBee_phi(thetaRi)
                                                  * ( AdvectedValueRiRi - AdvectedValueRi );
         fri = 0.5 * ( ur * ( cRip12 + cLip12 ) - fabs(ur) * ( cRip12 - cLip12 ) );
      }

      // Left (X)
      AdvectorValueToLe = UF->DOF_value( i+shift.i-1, j+shift.j, k, 0, advecting_level );
      AdvectorValueBoLe = UF->DOF_value( i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
      ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
      if ( UF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT
        || UF->DOF_color( i-1, j, k, component ) == FV_BC_BOTTOM_LEFT
        || UF->DOF_color( i-1, j, k, component ) == FV_BC_TOP_LEFT) {
         if ( ul > 0. ) fle = ul * AdvectedValueLe;
         else fle = ul * AdvectedValueC;
      } else {
         xl =UF->get_DOF_coordinate( i+shift.i-1, 0, 0 );
         xL =UF->get_DOF_coordinate( i-1, component, 0 );
         dxl  = xC - xL;
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
            cRim12 = AdvectedValueC;
         else {
            xR =UF->get_DOF_coordinate( i+1, component, 0 );
            dxr  = xR - xC;
            dxCl = xC - xl;
            cRim12 = AdvectedValueC - ( dxCl / dxr ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                      * ( AdvectedValueRi - AdvectedValueC );
         }
         dxLl = xl - xL;
         AdvectedValueLeLe =UF->DOF_value( i-2, j, k, component, advected_level );

         thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
         ( AdvectedValueLe - AdvectedValueLeLe ) / ( AdvectedValueC - AdvectedValueLe ) : 1.e20;
         cLim12 = AdvectedValueLe + ( dxLl / dxl ) * FV_DiscreteField::SuperBee_phi(thetaLe)
                                                     * ( AdvectedValueC - AdvectedValueLe );
         fle = 0.5 * ( ul * ( cRim12 + cLim12 ) - fabs(ul) * ( cRim12 - cLim12 ) );
      }

      // Top and Bottom
      // --------------
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP ) {
         AdvectorValueTo = AdvectorValueC;
         AdvectedValueTo = AdvectedValueC;
      } else {
         AdvectorValueTo = UF->DOF_value( i, j+1, k, component, advecting_level );
         AdvectedValueTo =UF->DOF_value( i, j+1, k, component, advected_level );
      }

      if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM ) {
         AdvectorValueBo = AdvectorValueC;
         AdvectedValueBo = AdvectedValueC;
      } else {
         AdvectorValueBo = UF->DOF_value( i, j-1, k, component, advecting_level );
         AdvectedValueBo =UF->DOF_value( i, j-1, k, component, advected_level );
      }

      thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ?
      ( AdvectedValueC - AdvectedValueBo ) / ( AdvectedValueTo - AdvectedValueC ) : 1.e20;

      // Top (Y)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
         fto = AdvectorValueC * AdvectedValueC;
      else {
         vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
         if ( UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP ) {
            if ( vt > 0. ) fto = vt * AdvectedValueC;
            else fto = vt * AdvectedValueTo;
         } else {
            yt =UF->get_DOF_coordinate( j+shift.j, 0, 1 );
            yT =UF->get_DOF_coordinate( j+1, component, 1 );
            dyCt = yt - yC;
            dyt  = yT - yC;
            cLip12 = AdvectedValueC + ( dyCt / dyt ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                      * ( AdvectedValueTo - AdvectedValueC );

            dyTt = yT - yt;
            dyT =UF->get_cell_size( j+1, component, 1 );
            AdvectedValueToTo =UF->DOF_value( i, j+2, k, component, advected_level );

            thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 1.e-20 ?
            ( AdvectedValueTo - AdvectedValueC ) / ( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
            cRip12 = AdvectedValueTo - ( dyTt / dyT ) * FV_DiscreteField::SuperBee_phi(thetaTo)
                                                     * ( AdvectedValueToTo - AdvectedValueTo );
            fto = 0.5 * ( vt * ( cRip12 + cLip12 ) - fabs(vt) * ( cRip12 - cLip12 ) );
         }
      }

      // Bottom (Y)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
         fbo = AdvectorValueC * AdvectedValueC;
      else {
         vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
         if ( UF->DOF_color(i,j-1,k,component) == FV_BC_BOTTOM ) {
            if ( vb > 0. ) fbo = vb * AdvectedValueBo;
            else fbo = vb * AdvectedValueC;
         } else {
            yb =UF->get_DOF_coordinate( j+shift.j-1, 0, 1 );
            yB =UF->get_DOF_coordinate( j-1, component, 1 );
            dyb  = yC - yB;

            dyBb = yb - yB;
            AdvectedValueBoBo =UF->DOF_value( i, j-2, k, component, advected_level );

            thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
            ( AdvectedValueBo - AdvectedValueBoBo ) / ( AdvectedValueC - AdvectedValueBo ) : 1.e20;
            cLim12 = AdvectedValueBo + ( dyBb / dyb ) * FV_DiscreteField::SuperBee_phi(thetaBo)
                                                        * ( AdvectedValueC - AdvectedValueBo );

            if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
               cRim12 = AdvectedValueC;
            else {
               yT =UF->get_DOF_coordinate( j+1, component, 1 );
               dyt  = yT - yC;
               dyCb = yC - yb;
               cRim12 = AdvectedValueC - ( dyCb / dyt ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                         * ( AdvectedValueTo - AdvectedValueC );
            }
            fbo = 0.5 * ( vb * ( cRim12 + cLim12 ) - fabs(vb) * ( cRim12 - cLim12 ) );
         }
      }

      if (dim == 3) {
         // Front and Behind
         // ----------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT ) {
            AdvectorValueFr = AdvectorValueC;
            AdvectedValueFr = AdvectedValueC;
         } else {
            AdvectorValueFr = UF->DOF_value( i, j, k+1, component, advecting_level );
            AdvectedValueFr =UF->DOF_value( i, j, k+1, component, advected_level );
         }

         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND ) {
            AdvectorValueBe = AdvectorValueC;
            AdvectedValueBe = AdvectedValueC;
         } else {
            AdvectorValueBe = UF->DOF_value( i, j, k-1, component, advecting_level );
            AdvectedValueBe =UF->DOF_value( i, j, k-1, component, advected_level );
         }

         thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ?
         ( AdvectedValueC - AdvectedValueBe ) / ( AdvectedValueFr - AdvectedValueC ) : 1.e20;

         // Front (Z)
         AdvectorValueFrBo = UF->DOF_value( i, j+shift.j-1, k+shift.k, 2, advecting_level );
         AdvectorValueFrTo = UF->DOF_value( i, j+shift.j, k+shift.k, 2, advecting_level );
         wf = 0.5 * ( AdvectorValueFrBo + AdvectorValueFrTo );
         if ( UF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT
           || UF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT_BOTTOM
           || UF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT_TOP ) {
            if ( wf > 0. ) ffr = wf * AdvectedValueC;
            else ffr = wf * AdvectedValueFr;
         } else {
            zf =UF->get_DOF_coordinate( k+shift.k, 2, 2 );
            zF =UF->get_DOF_coordinate( k+1, component, 2 );
            dzCf = zf - zC;
            dzf  = zF - zC;
            cLip12 = AdvectedValueC + ( dzCf / dzf ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                      * ( AdvectedValueFr - AdvectedValueC );
            dzFf = zF - zf;
            dzF =UF->get_cell_size( k+1, component, 2 );
            AdvectedValueFrFr =UF->DOF_value( i, j, k+2, component, advected_level );

            thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 1.e-20 ?
            ( AdvectedValueFr - AdvectedValueC ) / ( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
            cRip12 = AdvectedValueFr - ( dzFf / dzF ) * FV_DiscreteField::SuperBee_phi(thetaFr)
                                                     * ( AdvectedValueFrFr - AdvectedValueFr );
            ffr = 0.5 * ( wf * ( cRip12 + cLip12 ) - fabs(wf) * ( cRip12 - cLip12 ) );
         }

         // Behind (Z)
         AdvectorValueBeBo = UF->DOF_value( i, j+shift.j-1, k+shift.k-1, 2, advecting_level );
         AdvectorValueBeTo = UF->DOF_value( i, j+shift.j, k+shift.k-1, 2, advecting_level );
         wb = 0.5 * ( AdvectorValueBeBo + AdvectorValueBeTo );
         if ( UF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND
           || UF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND_BOTTOM
           || UF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND_TOP ) {
            if ( wb > 0. ) fbe = wb * AdvectedValueBe;
            else fbe = wb * AdvectedValueC;
         } else {
            zb =UF->get_DOF_coordinate( k+shift.k-1, 2, 2 );
            zB =UF->get_DOF_coordinate( k-1, component, 2 );
            dzb  = zC - zB;
            if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
               cRim12 = AdvectedValueC;
            else {
               zF =UF->get_DOF_coordinate( k+1, component, 2 );
               dzf  = zF - zC;
               dzCb = zC - zb;
               cRim12 = AdvectedValueC - ( dzCb / dzf ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                         * ( AdvectedValueFr - AdvectedValueC );
            }
            dzBb = zb - zB;
            AdvectedValueBeBe =UF->DOF_value( i, j, k-2, component, advected_level );

            thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
            ( AdvectedValueBe - AdvectedValueBeBe ) / ( AdvectedValueC - AdvectedValueBe ) : 1.e20;
            cLim12 = AdvectedValueBe + ( dzBb / dzb ) * FV_DiscreteField::SuperBee_phi(thetaBe)
                                                        * ( AdvectedValueC - AdvectedValueBe );
            fbe = 0.5 * ( wb * ( cRim12 + cLim12 ) - fabs(wb) * ( cRim12 - cLim12 ) );
         }
      }
   } else if (component == 2) {
      // The Third component (w)
      // Right and Left
      // --------------
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT ) {
         AdvectorValueRi = AdvectorValueC;
         AdvectedValueRi = AdvectedValueC;
      } else {
         AdvectorValueRi = UF->DOF_value( i+1, j, k, component, advecting_level );
         AdvectedValueRi =UF->DOF_value( i+1, j, k, component, advected_level );
      }

      if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT ) {
         AdvectorValueLe = AdvectorValueC;
         AdvectedValueLe = AdvectedValueC;
      } else {
         AdvectorValueLe = UF->DOF_value( i-1, j, k, component, advecting_level );
         AdvectedValueLe =UF->DOF_value( i-1, j, k, component, advected_level );
      }

      thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ?
      ( AdvectedValueC - AdvectedValueLe ) / ( AdvectedValueRi - AdvectedValueC ) : 1.e20;

      // Right (X)
      AdvectorValueFrRi = UF->DOF_value( i+shift.i, j, k+shift.k, 0, advecting_level );
      AdvectorValueBeRi = UF->DOF_value( i+shift.i, j, k+shift.k-1, 0, advecting_level );
      ur = 0.5 * (AdvectorValueFrRi + AdvectorValueBeRi);
      if ( UF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
        || UF->DOF_color( i+1, j, k, component ) == FV_BC_BEHIND_RIGHT
        || UF->DOF_color( i+1, j, k, component ) == FV_BC_FRONT_RIGHT ) {
         if ( ur > 0. ) fri = ur * AdvectedValueC;
         else fri = ur * AdvectedValueRi;
      } else {
         xr =UF->get_DOF_coordinate( i+shift.i, 0, 0 );
         xR =UF->get_DOF_coordinate( i+1, component, 0 );
         dxCr = xr - xC;
         dxr  = xR - xC;
         cLip12 = AdvectedValueC + ( dxCr / dxr ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                   * ( AdvectedValueRi - AdvectedValueC );

         dxRr = xR - xr;
         dxR =UF->get_cell_size( i+1, component, 0 );
         AdvectedValueRiRi =UF->DOF_value( i+2, j, k, component, advected_level );

         thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 1.e-20 ?
         ( AdvectedValueRi - AdvectedValueC ) / ( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
         cRip12 = AdvectedValueRi - ( dxRr / dxR ) * FV_DiscreteField::SuperBee_phi(thetaRi)
                                                  * ( AdvectedValueRiRi - AdvectedValueRi );
         fri = 0.5 * ( ur * ( cRip12 + cLip12 ) - fabs(ur) * ( cRip12 - cLip12 ) );
      }

      // Left (X)
      AdvectorValueFrLe = UF->DOF_value(i+shift.i-1, j, k+shift.k, 0, advecting_level );
      AdvectorValueBeLe = UF->DOF_value(i+shift.i-1, j, k+shift.k-1, 0, advecting_level );
      ul = 0.5 * ( AdvectorValueFrLe + AdvectorValueBeLe );
      if ( UF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT
        || UF->DOF_color( i-1, j, k, component ) == FV_BC_BEHIND_LEFT
        || UF->DOF_color( i-1, j, k, component ) == FV_BC_FRONT_LEFT ) {
         if ( ul > 0. ) fle = ul * AdvectedValueLe;
         else fle = ul * AdvectedValueC;
      } else {
         xl =UF->get_DOF_coordinate( i+shift.i-1, 0, 0 );
         xL =UF->get_DOF_coordinate( i-1, component, 0 );
         dxl  = xC - xL;
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
            cRim12 = AdvectedValueC;
         else {
            xR =UF->get_DOF_coordinate( i+1, component, 0 );
            dxr  = xR - xC;
            dxCl = xC - xl;
            cRim12 = AdvectedValueC - ( dxCl / dxr ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                      * ( AdvectedValueRi - AdvectedValueC );
         }

         dxLl = xl - xL;
         AdvectedValueLeLe =UF->DOF_value( i-2, j, k, component, advected_level );

         thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
         ( AdvectedValueLe - AdvectedValueLeLe ) / ( AdvectedValueC - AdvectedValueLe ) : 1.e20;
         cLim12 = AdvectedValueLe + ( dxLl / dxl ) * FV_DiscreteField::SuperBee_phi(thetaLe)
                                                    * ( AdvectedValueC - AdvectedValueLe );
         fle = 0.5 * ( ul * ( cRim12 + cLim12 ) - fabs(ul) * ( cRim12 - cLim12 ) );
      }

      // Top and Bottom
      // --------------
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP ) {
         AdvectorValueTo = AdvectorValueC;
         AdvectedValueTo = AdvectedValueC;
      } else {
         AdvectorValueTo = UF->DOF_value( i, j+1, k, component, advecting_level );
         AdvectedValueTo =UF->DOF_value( i, j+1, k, component, advected_level );
      }

      if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM ) {
         AdvectorValueBo = AdvectorValueC;
         AdvectedValueBo = AdvectedValueC;
      } else {
         AdvectorValueBo = UF->DOF_value( i, j-1, k, component, advecting_level );
         AdvectedValueBo =UF->DOF_value( i, j-1, k, component, advected_level );
      }

      thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ?
      ( AdvectedValueC - AdvectedValueBo ) / ( AdvectedValueTo - AdvectedValueC ) : 1.e20;

      // Top (Y)
      AdvectorValueBeTo = UF->DOF_value( i, j+shift.j, k+shift.k-1, 1, advecting_level );
      AdvectorValueFrTo = UF->DOF_value( i, j+shift.j, k+shift.k, 1, advecting_level );
      vt = 0.5 * ( AdvectorValueBeTo + AdvectorValueFrTo );
      if ( UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP
        || UF->DOF_color( i, j+1, k, component ) == FV_BC_BEHIND_TOP
        || UF->DOF_color( i, j+1, k, component ) == FV_BC_FRONT_TOP ) {
         if ( vt > 0. ) fto = vt * AdvectedValueC;
         else fto = vt * AdvectedValueTo;
      } else {
         yt =UF->get_DOF_coordinate( j+shift.j, 1, 1 );
         yT =UF->get_DOF_coordinate( j+1, component, 1 );
         dyCt = yt - yC;
         dyt  = yT - yC;
         cLip12 = AdvectedValueC + ( dyCt / dyt ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                   * ( AdvectedValueTo - AdvectedValueC );
         dyTt = yT - yt;
         dyT =UF->get_cell_size( j+1, component, 1 );
         AdvectedValueToTo =UF->DOF_value( i, j+2, k, component, advected_level );

         thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 1.e-20 ?
         ( AdvectedValueTo - AdvectedValueC ) / ( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
         cRip12 = AdvectedValueTo - ( dyTt / dyT ) * FV_DiscreteField::SuperBee_phi(thetaTo)
                                                  * ( AdvectedValueToTo - AdvectedValueTo );
         fto = 0.5 * ( vt * ( cRip12 + cLip12 ) - fabs(vt) * ( cRip12 - cLip12 ) );
      }

      // Bottom (Y)
      AdvectorValueBeBo = UF->DOF_value( i, j+shift.j-1, k+shift.k-1, 1, advecting_level );
      AdvectorValueFrBo = UF->DOF_value( i, j+shift.j-1, k+shift.k, 1, advecting_level );
      vb = 0.5 * ( AdvectorValueBeBo + AdvectorValueFrBo );
      if ( UF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
        || UF->DOF_color( i, j-1, k, component ) == FV_BC_BEHIND_BOTTOM
        || UF->DOF_color( i, j-1, k, component ) == FV_BC_FRONT_BOTTOM ) {
         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
         else fbo = vb * AdvectedValueC;
      } else {
         yb =UF->get_DOF_coordinate( j+shift.j-1, 1, 1 );
         yB =UF->get_DOF_coordinate( j-1, component, 1 );
         dyb  = yC - yB;
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
            cRim12 = AdvectedValueC;
         else {
            yT =UF->get_DOF_coordinate( j+1, component, 1 );
            dyt  = yT - yC;
            dyCb = yC - yb;
            cRim12 = AdvectedValueC - ( dyCb / dyt ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                      * ( AdvectedValueTo - AdvectedValueC );
         }
         dyBb = yb - yB;
         AdvectedValueBoBo =UF->DOF_value( i, j-2, k, component, advected_level );

         thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
         ( AdvectedValueBo - AdvectedValueBoBo ) / ( AdvectedValueC - AdvectedValueBo ) : 1.e20;
         cLim12 = AdvectedValueBo + ( dyBb / dyb ) * FV_DiscreteField::SuperBee_phi(thetaBo)
                                                     * ( AdvectedValueC - AdvectedValueBo );
         fbo = 0.5 * ( vb * ( cRim12 + cLim12 ) - fabs(vb) * ( cRim12 - cLim12 ) );
      }

      // Front and Behind
      // ----------------
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT ) {
         AdvectorValueFr = AdvectorValueC;
         AdvectedValueFr = AdvectedValueC;
      } else {
         AdvectorValueFr = UF->DOF_value( i, j, k+1, component, advecting_level );
         AdvectedValueFr =UF->DOF_value( i, j, k+1, component, advected_level );
      }

      if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND ) {
         AdvectorValueBe = AdvectorValueC;
         AdvectedValueBe = AdvectedValueC;
      } else {
         AdvectorValueBe = UF->DOF_value( i, j, k-1, component, advecting_level );
         AdvectedValueBe =UF->DOF_value( i, j, k-1, component, advected_level );
      }

      thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ?
      ( AdvectedValueC - AdvectedValueBe ) / ( AdvectedValueFr - AdvectedValueC ) : 1.e20;

      // Front (Z)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
         ffr = AdvectorValueC * AdvectedValueC;
      else {
         wf = 0.5 * ( AdvectorValueFr + AdvectorValueC );
         if ( UF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT ) {
            if ( wf > 0. ) ffr = wf * AdvectedValueC;
            else ffr = wf * AdvectedValueFr;
         } else {
            zf =UF->get_DOF_coordinate( k+shift.k, 0, 2 );
            zF =UF->get_DOF_coordinate( k+1, component, 2 );
            dzCf = zf - zC;
            dzf  = zF - zC;
            cLip12 = AdvectedValueC + ( dzCf / dzf ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                      * ( AdvectedValueFr - AdvectedValueC );

            dzFf = zF - zf;
            dzF =UF->get_cell_size( k+1, component, 2 );
            AdvectedValueFrFr =UF->DOF_value( i, j, k+2, component, advected_level );

            thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 1.e-20 ?
            ( AdvectedValueFr - AdvectedValueC ) / ( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
            cRip12 = AdvectedValueFr - ( dzFf / dzF ) * FV_DiscreteField::SuperBee_phi(thetaFr)
                                                     * ( AdvectedValueFrFr - AdvectedValueFr );
            ffr = 0.5 * ( wf * ( cRip12 + cLip12 ) - fabs(wf) * ( cRip12 - cLip12 ) );
         }
      }

      // Behind (Z)
      if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND )
         fbe = AdvectorValueC * AdvectedValueC;
      else {
         wb = 0.5 * (AdvectorValueBe + AdvectorValueC);
         if ( UF->DOF_color(i, j, k-1, component ) == FV_BC_BEHIND ) {
            if (wb > 0.) fbe = wb * AdvectedValueBe;
            else fbe = wb * AdvectedValueC;
         } else {
            zb =UF->get_DOF_coordinate( k+shift.k-1, 0, 2 );
            zB =UF->get_DOF_coordinate( k-1, component, 2 );
            dzb  = zC - zB;
            if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
               cRim12 = AdvectedValueC;
            else {
               zF =UF->get_DOF_coordinate( k+1, component, 2 );
               dzf  = zF - zC;
               dzCb = zC - zb;
               cRim12 = AdvectedValueC - ( dzCb / dzf ) * FV_DiscreteField::SuperBee_phi(thetaC)
                                                         * ( AdvectedValueFr - AdvectedValueC );
            }
            dzBb = zb - zB;
            AdvectedValueBeBe =UF->DOF_value( i, j, k-2, component, advected_level );

            thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
            ( AdvectedValueBe - AdvectedValueBeBe ) / ( AdvectedValueC - AdvectedValueBe ) : 1.e20;
            cLim12 = AdvectedValueBe + ( dzBb / dzb ) * FV_DiscreteField::SuperBee_phi(thetaBe)
                                                        * ( AdvectedValueC - AdvectedValueBe );
            fbe = 0.5 * ( wb * ( cRim12 + cLim12 ) - fabs(wb) * ( cRim12 - cLim12 ) );
         }
      }
   }

   if (dim == 2) {
      flux = ( fto - fbo ) * dxC + ( fri - fle ) * dyC;
   } else if (dim == 3) {
      flux = ( fto - fbo ) * dxC * dzC + ( fri - fle ) * dyC * dzC + ( ffr - fbe ) * dxC * dyC;
   }
   return ( coef * flux );
}
