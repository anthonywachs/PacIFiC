#include <DS_NavierStokes.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DomainBuilder.hh>
#include <FV_DiscreteField.hh>
#include <DS_NavierStokesSystem.hh>
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
#include <algorithm>
#include <DS_AllRigidBodies.hh>


//---------------------------------------------------------------------------
DS_NavierStokes*
DS_NavierStokes:: create( MAC_Object* a_owner,
		                    MAC_ModuleExplorer const* exp,
                          struct DS2NS const& transfer )
//---------------------------------------------------------------------------
{
 MAC_LABEL( "DS_NavierStokes:: create" ) ;
 MAC_CHECK_PRE( exp != 0 ) ;

 DS_NavierStokes* result =
                      new DS_NavierStokes( a_owner, exp, transfer ) ;

 MAC_CHECK_POST( result != 0 ) ;
 MAC_CHECK_POST( result->owner() == a_owner ) ;

 return( result ) ;

}




//---------------------------------------------------------------------------
DS_NavierStokes:: DS_NavierStokes( MAC_Object* a_owner,
                              	  MAC_ModuleExplorer const* exp,
                                   struct DS2NS const& fromDS )
//---------------------------------------------------------------------------
   : MAC_Object( a_owner )
	, ComputingTime("Solver")
   , UF ( fromDS.dom_->discrete_field( "velocity" ) )
   , PF ( fromDS.dom_->discrete_field( "pressure" ) )
   , GLOBAL_EQ( 0 )
   , mu( fromDS.mu_ )
   , kai( fromDS.kai_ )
   , AdvectionScheme( fromDS.AdvectionScheme_ )
	, StencilCorrection( fromDS.StencilCorrection_ )
	, FluxRedistThres( fromDS.FluxRedistThres_ )
   , AdvectionTimeAccuracy( fromDS.AdvectionTimeAccuracy_ )
   , rho( fromDS.rho_ )
	, b_restart ( fromDS.b_restart_ )
   , is_solids( fromDS.is_solids_ )
	, is_stressCal ( fromDS.is_stressCal_ )
	, ViscousStressOrder ( fromDS.ViscousStressOrder_ )
	, PressureStressOrder ( fromDS.PressureStressOrder_ )
	, stressCalFreq ( fromDS.stressCalFreq_ )
	, is_par_motion ( fromDS.is_par_motion_ )
	, allrigidbodies ( fromDS.allrigidbodies_ )
   , b_projection_translation( fromDS.dom_->primary_grid()
														->is_translation_active() )
   , b_grid_has_been_translated_since_last_output( false )
   , b_grid_has_been_translated_at_previous_time( false )
   , critical_distance_translation( fromDS.critical_distance_translation_ )
   , translation_direction( 0 )
   , bottom_coordinate( 0. )
   , translated_distance( 0. )
   , gravity_vector( 0 )
	, exceed (false)
	, turn (false)
{
   MAC_LABEL( "DS_NavierStokes:: DS_NavierStokes" ) ;
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
   if ( my_rank == is_master ) {
     CT_set_start();
     SCT_insert_app("Objects_Creation");
     SCT_set_start("Objects_Creation");
   }

	if ( AdvectionScheme == "TVD"
     && UF->primary_grid()->get_security_bandwidth() < 2 ) {
     string error_message="   >= 2 with TVD scheme";
     MAC_Error::object()->raise_bad_data_value( exp,
        "security_bandwidth", error_message );
   }

   // Get space dimension
   dim = UF->primary_grid()->nb_space_dimensions() ;
   nb_comps[0] = PF->nb_components() ;
   nb_comps[1] = UF->nb_components() ;

   if ( dim == 1 ) {
     string error_message="Space dimension should either 2 or 3";
     MAC_Error::object()->raise_bad_data_value(exp,
		  													  "nb_space_dimensions",
															  error_message );
   }

   // Create the Direction Splitting subcommunicators
   create_DS_subcommunicators();

   if ( is_stressCal == true &&
        UF->primary_grid()->get_security_bandwidth() < 4 ) {
     string error_message="   >= 4 for correct stress calculations on solids";
     MAC_Error::object()->raise_bad_data_value( exp,
        "security_bandwidth", error_message );
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
   struct NS2System inputData;
   inputData.is_solids_ = is_solids ;
   inputData.is_stressCal_ = is_stressCal ;

   // Build the matrix system
   MAC_ModuleExplorer* se = exp->create_subexplorer( 0,"DS_NavierStokesSystem" ) ;
   GLOBAL_EQ = DS_NavierStokesSystem::create( this, se, UF, PF, inputData ) ;
   se->destroy() ;

	if (UF->primary_grid()->is_periodic_flow_rate()) {
		// controller = DS_PID::create(0.5,1.,0);
		// Set the initial pressure drop, required in case of given flow rate
		if (!b_restart) {
			pressure_drop = -100.;
			string fileName = "./DS_results/flow_pressure_history.csv" ;
			std::ofstream MyFile;
			if (macCOMM->rank() == 0) {
				MyFile.open( fileName.c_str()) ;
				MyFile << "t Qc dP exceed turn" << endl;
			}
			MyFile.close();
		} else if (b_restart) {
			string fileName = "./DS_results/flow_pressure_history.csv" ;
			std::ifstream MyFile;
			MyFile.open( fileName.c_str() ) ;
			MyFile.seekg(-31,ios_base::end);
	      MyFile >> Qold >> pressure_drop >> exceed >> turn;
			MyFile.close();
		}
	}

   // Timing routines
   if ( my_rank == is_master ) {
     SCT_insert_app("Matrix_Assembly&Initialization");
	  SCT_insert_app("Matrix_RE_Assembly&Initialization");
	  // SCT_insert_app("Stencil");
	  // SCT_insert_app("Schur");
     SCT_insert_app("Pressure predictor");
     SCT_insert_app("Velocity update");
     SCT_insert_app("Penalty Step");
     SCT_insert_app("Pressure Update");
     SCT_insert_app("Viscous stress");
     SCT_insert_app("Pressure stress");
     SCT_get_elapsed_time("Objects_Creation");
   }
}

//---------------------------------------------------------------------------
DS_NavierStokes:: ~DS_NavierStokes( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: ~DS_NavierStokes" ) ;

   free_DS_subcommunicators() ;

}

//---------------------------------------------------------------------------
void
DS_NavierStokes:: do_one_inner_iteration( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: do_one_inner_iteration" ) ;

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

}

//---------------------------------------------------------------------------
void
DS_NavierStokes:: do_before_time_stepping( FV_TimeIterator const* t_it,
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: do_before_time_stepping" ) ;

	if ( my_rank == is_master ) SCT_set_start("Matrix_Assembly&Initialization");

   allocate_mpi_variables(PF);
   allocate_mpi_variables(UF);

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

   if (is_solids) {
		// Build void frac and intersection variable
		allrigidbodies->build_solid_variables_on_fluid_grid(PF);
		allrigidbodies->build_solid_variables_on_fluid_grid(UF);
		// Compute void fraction for pressure and velocity field
		allrigidbodies->compute_void_fraction_on_grid(PF, false);
		allrigidbodies->compute_void_fraction_on_grid(UF, false);
		// Compute intersection with RB for pressure and velocity field
		allrigidbodies->compute_grid_intersection_with_rigidbody(PF, false);
		allrigidbodies->compute_grid_intersection_with_rigidbody(UF, false);
		// Compute the face fractions of CutCell is active
		if (StencilCorrection == "CutCell") {
			allrigidbodies->compute_cutCell_geometric_parameters(PF);
			allrigidbodies->compute_cutCell_geometric_parameters(UF);
		}

		if (my_rank == 0)
         cout << "Finished void fraction and grid intersection... \n" << endl;

		// Field initialization
      vector<size_t> vec{ 0, 1, 3};
      if (dim == 3) vec.push_back(4);
      initialize_grid_nodes_on_rigidbody(vec);

      if (my_rank == 0)
         cout << "Finished field initializations... \n" << endl;

   }

   // Direction splitting
   // Assemble 1D tridiagonal matrices
   assemble_1D_matrices(PF,t_it);
   assemble_1D_matrices(UF,t_it);

   if (my_rank == 0)
      cout << "Finished assembling pre-coefficient matrix... \n" << endl;

   if ( my_rank == is_master )
      SCT_get_elapsed_time( "Matrix_Assembly&Initialization" );

}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: do_after_time_stepping" ) ;

   // Elapsed time by sub-problems

   // SCT_set_start( "Writing CSV" );
   // write_output_field(PF);
   // write_output_field(UF);
   // SCT_get_elapsed_time( "Writing CSV" );

   output_L2norm_velocity(0);
   output_L2norm_pressure(0);

   if ( my_rank == is_master )
   {
     double cputime = CT_get_elapsed_time();
     cout << endl
	  		 << "========================================================" << endl
			 << "                Navier Stokes Problem                   " << endl
			 << "========================================================" << endl;
     write_elapsed_time_smhd(cout,cputime,"Computation time");
     SCT_get_summary(cout,cputime);
   }

   deallocate_mpi_variables();
}

//---------------------------------------------------------------------------
void
DS_NavierStokes:: do_before_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: do_before_inner_iterations_stage" ) ;

	if ( my_rank == is_master )
		SCT_set_start( "Matrix_RE_Assembly&Initialization" );

   if ((is_par_motion) && (is_solids)) {
		// Solve equation of motion for all RB and update pos,vel
		allrigidbodies->solve_RB_equation_of_motion(t_it);
		allrigidbodies->generate_list_of_local_RB();
		allrigidbodies->initialize_surface_variables_for_all_RB();
		allrigidbodies->compute_surface_variables_for_all_RB();
		allrigidbodies->compute_halo_zones_for_all_rigid_body();
		allrigidbodies->create_neighbour_list_for_AllRB();
		// Clear void fraction and intersection data
		allrigidbodies->clear_GrainsRB_data_on_grid(PF);
		allrigidbodies->clear_GrainsRB_data_on_grid(UF);
		// Compute void fraction for pressure and velocity field
		allrigidbodies->compute_void_fraction_on_grid(PF, true);
		allrigidbodies->compute_void_fraction_on_grid(UF, true);
		// Compute intersection with RB for pressure and velocity field
		allrigidbodies->compute_grid_intersection_with_rigidbody(PF, true);
		allrigidbodies->compute_grid_intersection_with_rigidbody(UF, true);
		// Compute the face fractions of CutCell is active
		if (StencilCorrection == "CutCell") {
			allrigidbodies->compute_cutCell_geometric_parameters(PF);
			allrigidbodies->compute_cutCell_geometric_parameters(UF);
		}
		// Compute the fresh node coming in fluid in this step
		allrigidbodies->compute_fresh_nodes(PF);
		allrigidbodies->compute_fresh_nodes(UF);

		// Field initialization
		vector<size_t> vec{ 0, 1, 3};
		if (dim == 3) vec.push_back(4);
		initialize_grid_nodes_on_rigidbody(vec);

		// Extrapolate advecvtion term on fresh nodes
		allrigidbodies->extrapolate_scalar_on_fresh_nodes(UF,2);
		// Interpolate velocity on the fresh nodes in fluid domain
		allrigidbodies->interpolate_vector_on_fresh_nodes(UF,0);
		allrigidbodies->interpolate_vector_on_fresh_nodes(UF,1);
		allrigidbodies->interpolate_vector_on_fresh_nodes(UF,3);
		if (dim == 3)
			allrigidbodies->interpolate_vector_on_fresh_nodes(UF,4);
		// allrigidbodies->extrapolate_scalar_on_fresh_nodes(PF,0);
		// allrigidbodies->extrapolate_scalar_on_fresh_nodes(PF,1);

	   // Direction splitting
      // Assemble 1D tridiagonal matrices
      assemble_1D_matrices(PF,t_it);
      assemble_1D_matrices(UF,t_it);
   }

	//  Projection_Translation
	if ( b_projection_translation )
		b_grid_has_been_translated_at_previous_time = false;

	if ( my_rank == is_master )
		SCT_get_elapsed_time( "Matrix_RE_Assembly&Initialization" );

}

//---------------------------------------------------------------------------
void
DS_NavierStokes:: do_after_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: do_after_inner_iterations_stage" ) ;

   // Compute hydrodynamic forces by surface viscous stress
   if ( my_rank == is_master ) SCT_set_start( "Viscous stress" );
   if (is_stressCal && (t_it->iteration_number() % stressCalFreq == 0)) {
      allrigidbodies->compute_viscous_force_and_torque_for_allRB(ViscousStressOrder);
   }
   if ( my_rank == is_master ) SCT_get_elapsed_time( "Viscous stress" );

   output_L2norm_divergence();

/*
   double vel_divergence = get_velocity_divergence(t_it);

   if ( my_rank == is_master ) {
      string fileName = "./DS_results/max_divergence.csv" ;
      ofstream MyFile( fileName.c_str(), ios::app ) ;
      MyFile << t_it -> time() << "," << vel_divergence << endl;
      MyFile.close( ) ;
   }
*/

   double cfl = UF->compute_CFL( t_it, 0 );
   if ( my_rank == is_master )
      MAC::out() << "CFL: "<< cfl <<endl;

	// allrigidbodies->write_volume_conservation(t_it);

   // Projection translation
   if ( b_projection_translation ) {

      double min_coord = allrigidbodies->get_min_RB_coord(translation_direction);

      double distance_to_bottom = MAC::abs(min_coord-bottom_coordinate);

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

}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: do_additional_savings" ) ;

   GLOBAL_EQ->display_debug();
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: ugradu_initialization ( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_NavierStokes:: ugradu_initialization" ) ;

  for (size_t comp=0;comp<nb_comps[1];comp++) UF->set_DOFs_value( comp, 2, 0.);

}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: calculate_row_indexes ( FV_DiscreteField const* FF)
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: calculate_row_indexes" ) ;

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
         size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(field,dir,comp);
         switch (dir) {
            case 0:
             for (size_t j=min_unknown_index(1); j<=max_unknown_index(1); j++) {
              for (size_t k=min_unknown_index(2); k<=max_unknown_index(2); k++) {
                  size_t p = (j - min_unknown_index(1))
                           + (1 + max_unknown_index(1) - min_unknown_index(1))
                           * (k - min_unknown_index(2));
                  row_index->operator()(j,k) = p;
              }
             }
             break;
            case 1:
             for (size_t i=min_unknown_index(0); i<=max_unknown_index(0); i++) {
              for (size_t k=min_unknown_index(2); k<=max_unknown_index(2); k++) {
                  size_t p = (i - min_unknown_index(0))
                           + (1 + max_unknown_index(0) - min_unknown_index(0))
                           * (k - min_unknown_index(2));
						row_index->operator()(i,k) = p;
              }
             }
             break;
            case 2:
             for (size_t i=min_unknown_index(0); i<=max_unknown_index(0); i++) {
              for (size_t j=min_unknown_index(1); j<=max_unknown_index(1); j++) {
                  size_t p = (i - min_unknown_index(0))
                           + (1 + max_unknown_index(0) - min_unknown_index(0))
                           * (j - min_unknown_index(1));
						row_index->operator()(i,j) = p;
              }
             }
             break;
         }
      }
   }
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: assemble_field_matrix ( FV_DiscreteField const* FF
                                        , FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: assemble_field_matrix" ) ;

	// if ( my_rank == is_master )
	// 	SCT_set_start("Stencil");

	// Assemble the matrices for pressure field(0) and velocity(1) field
	size_t field = (FF == PF) ? 0 : 1 ;
	double gamma = mu/2.0;
	TDMatrix* A = GLOBAL_EQ-> get_A(field);

	size_t_vector min_unknown_index(3,0);
	size_t_vector max_unknown_index(3,0);

	// doubleArray2D* CC_face_frac = allrigidbodies->get_CC_face_fraction(FF);
	// doubleArray2D* CC_vol = allrigidbodies->get_CC_cell_volume(FF);

	for (size_t comp=0;comp<nb_comps[field];comp++) {
		// Get local min and max indices
		for (size_t l=0;l<dim;++l) {
			min_unknown_index(l) =
								FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
			max_unknown_index(l) =
								FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
		}

		for (size_t dir = 0; dir < dim; dir++) {

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

			size_t dir_j = (dir == 0) ? 1 : 0;
			size_t dir_k = (dir == 2) ? 1 : 2;

			size_t local_min_k = (dim == 2) ? 0 : min_unknown_index(dir_k);
			size_t local_max_k = (dim == 2) ? 0 : max_unknown_index(dir_k);

			size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(field,dir,comp);

			for (size_t j = min_unknown_index(dir_j);
						  j <= max_unknown_index(dir_j);++j) {
				for (size_t k = local_min_k; k <= local_max_k; ++k) {
					size_t r_index = row_index->operator()(j,k);
   				// Perform assembling
					double Aee_diagcoef = 0.;
   				for (size_t m = 0, i = min_unknown_index(dir);
											 i <= max_unknown_index(dir); ++i,++m) {
						int i_temp = (int)((dir == 0) ? i : min_unknown_index(0));
						int j_temp = (int)((dir == 1) ? i : min_unknown_index(1));
						int k_temp = (int)((dir == 2) ? i : min_unknown_index(2));

      				double xC = FF->get_DOF_coordinate( i, comp, dir) ;

				      // Check if the index is at right domain
				      // boundary with neumann or dirichlet BC
						double xR = FF->get_DOF_coordinate( i+1,comp, dir) ;
				      if ( (i==max_unknown_index(dir))
				        && r_bound
				        && FF->DOF_on_BC(i_temp,j_temp,k_temp,comp)) {
				         xR = 0.;
				      }

				      // Check if the index is at left domain
				      // boundary with neumann or dirichlet BC
						double xL = FF->get_DOF_coordinate( i-1, comp, dir) ;
				      if ( (i==min_unknown_index(dir))
				        && l_bound
				        && FF->DOF_on_BC(i_temp,j_temp,k_temp,comp)) {
				         xL = 0.;
				      }

      				double dxr = xR - xC;
      				double dxl = xC - xL;
						double dx = FF->get_cell_size( i, comp, dir);
						// double pre_factor = (2. * dx) / (dxr + dxl);

						double right = (FF == PF) ? -1.0/dxr : -gamma/dxr ;
						double left = (FF == PF) ? -1.0/dxl : -gamma/dxl ;
						double unsteady_term = (FF == PF) ? 1.0*dx
																	 : rho*dx/t_it->time_step();

						double center = - (right+left);

				      if ((is_solids) && (FF == UF)) {
				         size_t_array2D* intersect_vector = allrigidbodies
														->get_intersect_vector_on_grid(FF);
				         doubleArray2D* intersect_distance = allrigidbodies
														->get_intersect_distance_on_grid(FF);
				         size_t_array2D* void_frac = allrigidbodies
														->get_void_fraction_on_grid(FF);

							size_t p = 0;
							if (dir == 0) {
								p = FF->DOF_local_number(i,j,k,comp);
							} else if (dir == 1) {
								p = FF->DOF_local_number(j,i,k,comp);
							} else if (dir == 2) {
								p = FF->DOF_local_number(j,k,i,comp);
							}

				         if (void_frac->operator()(p,0) == 0) {
				            // if left node is inside the solid particle
				            if (intersect_vector->operator()(p,2*dir+0) == 1) {
				               left = -gamma/intersect_distance->operator()(p,2*dir+0);
									// pre_factor = (2. * dx) / (dxr + intersect_distance->operator()(p,2*dir+0));
								}
				            // if right node is inside the solid particle
				            if (intersect_vector->operator()(p,2*dir+1) == 1) {
				               right = -gamma/intersect_distance->operator()(p,2*dir+1);
									// pre_factor = (2. * dx) / (dxl + intersect_distance->operator()(p,2*dir+1));
								}
				         } else if (void_frac->operator()(p,0) != 0) {
				            // if center node is inside the solid particle
				            left = 0.;
				            right = 0.;
				         }

         				center = -(right+left);

         				if (intersect_vector->operator()(p,2*dir+0) == 1)
								left = 0.;
         				if (intersect_vector->operator()(p,2*dir+1) == 1)
								right = 0.;
      				}

						// left *= pre_factor;
						// right *= pre_factor;
						// center *= pre_factor;

						double value = center;

				      // Condition for handling the pressure
						// neumann conditions at wall
				      if (i == min_unknown_index(dir) && l_bound) {
							int ii = (int)((dir == 0) ? i-1 : min_unknown_index(0));
							int jj = (int)((dir == 1) ? i-1 : min_unknown_index(1));
							int kk = (int)((dir == 2) ? i-1 : min_unknown_index(2));

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
							int ii = (int)((dir == 0) ? i+1 : max_unknown_index(0));
							int jj = (int)((dir == 1) ? i+1 : max_unknown_index(1));
							int kk = (int)((dir == 2) ? i+1 : max_unknown_index(2));

				         if (FF->DOF_in_domain(ii,jj,kk,comp)
				          && FF->DOF_has_imposed_Dirichlet_value
				                        ((size_t)ii,(size_t)jj,(size_t)kk,comp)) {
				            // For Dirichlet boundary condition
				            value = center;
				         } else {
				            // For Neumann homogeneous boundary condition
				            value = -left;
				         }
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


					// Storing Aee for MPI communication
					ProdMatrix* Ap = GLOBAL_EQ->get_Ap(field);
					double* local_coeff = data_for_S[field][dir].send[comp][0];
					size_t nbrow = Ap[dir].ei_ii_ie[comp]->nb_rows();
					size_t ii = (nbrow*nbrow + 1)*r_index;
					local_coeff[ii] = Aee_diagcoef;
				}
			}

		}
	}

	// if ( my_rank == is_master )
	// 	SCT_get_elapsed_time("Stencil");

}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: assemble_field_schur_matrix ( FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: assemble_field_schur_matrix" ) ;
   // Compute the product matrix for each proc

	// if ( my_rank == is_master )
	// 	SCT_set_start("Schur");

	size_t field = (FF == PF) ? 0 : 1;

   TDMatrix* A = GLOBAL_EQ-> get_A(field);
	ProdMatrix* Ap = GLOBAL_EQ->get_Ap(field);

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

      for (size_t dir = 0; dir < dim; dir++) {
         size_t dir_j = (dir == 0) ? 1 : 0;
			size_t dir_k = (dir == 2) ? 1 : 2;

			size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(field,dir,comp);

         size_t local_min_k = (dim == 2) ? 0 : min_unknown_index(dir_k);
         size_t local_max_k = (dim == 2) ? 0 : max_unknown_index(dir_k);

			if (nb_ranks_comm_i[dir] > 1) {

				size_t nbrow = Ap[dir].ei_ii_ie[comp]->nb_rows();
				double* local_packet = data_for_S[field][dir].send[comp][0];

				// Calculating the product matrix and creating the container for MPI
	         for (size_t j = min_unknown_index(dir_j);
	                    j <= max_unknown_index(dir_j);++j) {
	            for (size_t k = local_min_k; k <= local_max_k; ++k) {
						size_t r_index = row_index->operator()(j,k);

						GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir,field,r_index);

						// Create the data container to send and in the master as well
						for (size_t kk = 0; kk < nbrow; kk++) {
							for (size_t jj = 0; jj < nbrow; jj++) {
								size_t ii = (nbrow*nbrow + 1)*r_index
											 + nbrow*kk + jj + 1;
								local_packet[ii] = Ap[dir].ei_ii_ie[comp]->item(kk,jj);
							}
						}
					}
				}

				if (rank_in_i[dir] != 0 ) {
					// Send the packed data to master
			      MPI_Send( data_for_S[field][dir].send[comp][0],
			          (int) data_for_S[field][dir].size[comp],
			          	MPI_DOUBLE, 0, 0, DS_Comm_i[dir] ) ;
				}

				// Assemble the global product matrix by adding contribution from
				// all procs
				if (rank_in_i[dir] == 0 ) {
					for (size_t i = 1; i < (size_t)nb_ranks_comm_i[dir]; ++i) {
						// Recieve the data packet sent by master
						static MPI_Status status;
						MPI_Recv( data_for_S[field][dir].receive[comp][0],
							(int) data_for_S[field][dir].size[comp],
							MPI_DOUBLE, (int) i, 0, DS_Comm_i[dir], &status ) ;

		         	for (size_t j = min_unknown_index(dir_j);
		                    		j <= max_unknown_index(dir_j);++j) {
		            	for (size_t k = local_min_k; k <= local_max_k; ++k) {
								size_t r_index = row_index->operator()(j,k);
								size_t p = (nbrow*nbrow + 1)*r_index;
								double ee_proc0 = local_packet[p];
								A[dir].ee[comp][r_index]->set_item(0,0,ee_proc0);

								for (size_t kk = 0;kk < nbrow; kk++) {
									for (size_t jj = 0;jj < nbrow; jj++) {
										size_t ii = (nbrow*nbrow + 1)*r_index
													 + nbrow*kk + jj + 1;
										double value =
												data_for_S[field][dir].receive[comp][0][ii];
										local_packet[ii] += value;
									}
								}

								size_t ii = (nbrow*nbrow + 1)*r_index;
								double value_Aee =
												data_for_S[field][dir].receive[comp][0][ii];

								// Assemble the global Aee matrix
								// Only for (nb_proc-1) in case no PBC
								// Otherwise for (nb_proc)
								if (!is_periodic[field][dir]) {
									if (i<(size_t)nb_ranks_comm_i[dir]-1) {
										A[dir].ee[comp][r_index]->set_item(i,i,value_Aee);
									}
								} else {
               				A[dir].ee[comp][r_index]->set_item(i,i,value_Aee);
            				}
							}
						}
					}
				}

				// Assemble the schlur complement in the master proc
				if (rank_in_i[dir] == 0) {
					TDMatrix* Schur = GLOBAL_EQ-> get_Schur(field);
					for (size_t j = min_unknown_index(dir_j);
		                    j <= max_unknown_index(dir_j);++j) {
		            for (size_t k = local_min_k; k <= local_max_k; ++k) {
							size_t r_index = row_index->operator()(j,k);
							size_t schur_size = Schur[dir].ii_main[comp][r_index]
													  ->nb_rows();
							for (int p = 0; p < (int)schur_size; p++) {
								size_t ii = (nbrow*nbrow + 1)*r_index
											 + nbrow*p + p + 1;
								Schur[dir].ii_main[comp][r_index]
										->set_item(p,A[dir].ee[comp][r_index]->item(p,p)
											     								-local_packet[ii]);
								if (p < (int)schur_size-1)
									Schur[dir].ii_super[comp][r_index]
										->set_item(p,-local_packet[ii+1]);
								if (p > 0)
									Schur[dir].ii_sub[comp][r_index]
										->set_item(p-1,-local_packet[ii-1]);
								// In case of periodic and multi-processor,
								// there will be a variant of Tridiagonal matrix
								// instead of normal format
								if (is_periodic[field][dir] == 1) {
									ii = (nbrow*nbrow + 1)*r_index
										+ nbrow*p + schur_size + 1;
									Schur[dir].ie[comp][r_index]->set_item(p,0,
																				-local_packet[ii]);
									ii = (nbrow*nbrow + 1)*r_index
										+ nbrow*schur_size + p + 1;
									Schur[dir].ei[comp][r_index]->set_item(0,p,
																				-local_packet[ii]);
								}
							}
							// Pre-thomas treatment on Schur complement
							GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur,r_index);

							// In case of periodic and multi-processor, there will be
							// a variant of Tridiagonal matrix instead of normal format
							// So, Schur complement of Schur complement is calculated
							if (is_periodic[field][dir] == 1) {
								size_t ii = ((size_t)pow(nbrow,2) + 1)*r_index
											 + nbrow*schur_size
											 + schur_size + 1;
								Schur[dir].ee[comp][r_index]->set_item(0,0,
										A[dir].ee[comp][r_index]->item(schur_size,schur_size)
									 										- local_packet[ii]);

								ProdMatrix* SchurP = GLOBAL_EQ->get_SchurP(field);
								GLOBAL_EQ->compute_product_matrix_interior(Schur
															 ,SchurP,comp,0,dir,r_index);

								TDMatrix* DoubleSchur = GLOBAL_EQ-> get_DoubleSchur(field);
								DoubleSchur[dir].ii_main[comp][r_index]
								->set_item(0,Schur[dir].ee[comp][r_index]->item(0,0)
								-SchurP[dir].ei_ii_ie[comp]->item(0,0));
							}
						}
					}
				}
			// Condition for single processor in any
			// direction with periodic boundary conditions
			} else if (is_periodic[field][dir] == 1) {
				for (size_t j = min_unknown_index(dir_j);
								j <= max_unknown_index(dir_j);++j) {
					for (size_t k = local_min_k; k <= local_max_k; ++k) {
						size_t r_index = row_index->operator()(j,k);

						GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir,field,r_index);

						LA_SeqMatrix* product_matrix = Ap[dir].ei_ii_ie[comp];

						A[dir].ee[comp][r_index]->set_item(0,0
								,data_for_S[field][dir].send[comp][0][0]);

						TDMatrix* Schur = GLOBAL_EQ-> get_Schur(field);
						size_t nb_row = Schur[dir].ii_main[comp][r_index]->nb_rows();
						for (int p = 0; p < (int)nb_row; p++) {
							Schur[dir].ii_main[comp][r_index]
										->set_item(p,A[dir].ee[comp][r_index]->item(p,p)
														-product_matrix->item(p,p));
							if (p < (int)nb_row-1)
								Schur[dir].ii_super[comp][r_index]
										->set_item(p,-product_matrix->item(p,p+1));
							if (p > 0)
							Schur[dir].ii_sub[comp][r_index]
							->set_item(p-1,-product_matrix->item(p,p-1));
						}
						GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur,r_index);
					}
				}
			}
      }
   }

	// if ( my_rank == is_master )
	// 	SCT_get_elapsed_time("Schur");
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: assemble_1D_matrices ( FV_DiscreteField const* FF
                                        , FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: assemble_1D_matrices" ) ;

   // Assemble field matrix
   assemble_field_matrix (FF,t_it);
   // Calculate and assemble Schur complement
	assemble_field_schur_matrix(FF);
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: NS_first_step ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_NavierStokes:: NS_first_step" ) ;

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
				PF->set_DOF_value( i, j, k, 0, 1, value);
         }
     }
  }

  // Synchronize pressure field
  PF->synchronize(1);

  if (is_solids) {
	  allrigidbodies->extrapolate_pressure_inside_RB(PF,1);
  	// correct_pressure_1st_layer_solid(1);
  	// correct_pressure_2nd_layer_solid(1);
  }

  PF->set_neumann_DOF_values();

  // Calculate pressure forces on the solid particles
  if ( my_rank == is_master ) SCT_set_start( "Pressure stress" );
  if (is_stressCal && (t_it->iteration_number() % stressCalFreq == 0)) {
     allrigidbodies->compute_pressure_force_and_torque_for_allRB(PressureStressOrder);
  }
  if ( my_rank == is_master ) SCT_get_elapsed_time( "Pressure stress" );

}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: assemble_velocity_diffusion_terms ( )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: assemble_velocity_diffusion_terms" ) ;

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

		// Flux redistribution
		if (is_solids && (StencilCorrection == "CutCell")) {
			size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
			intVector* ownerID = allrigidbodies->get_CC_ownerID(UF);
			for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
				for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
					for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
						size_t p = UF->DOF_local_number(i,j,k,comp);
						if ((ownerID->operator()(p) != -1) && (void_frac->operator()(p,0) != 0)) {
							for (size_t dir = 0; dir < dim; dir++) {
								size_t_vector i0_r(3), i0_l(3);
								i0_l(0) = i; i0_l(1) = j; i0_l(2) = k;
								i0_r = i0_l;
								i0_l(dir) = i0_l(dir) - 1;
								i0_r(dir) = i0_r(dir) + 1;
								size_t p_lft = UF->DOF_local_number(i0_l(0),i0_l(1),i0_l(2),comp);
								size_t p_rht = UF->DOF_local_number(i0_r(0),i0_r(1),i0_r(2),comp);
								if (vel_diffusion[dir]->operator()(p) != 0.) {
									if (void_frac->operator()(p_lft,0) == 0) {
										vel_diffusion[dir]->operator()(p_lft) += vel_diffusion[dir]->operator()(p);
										vel_diffusion[dir]->operator()(p) = 0.;
									} else if (vel_diffusion[dir]->operator()(p_rht) == 0) {
										vel_diffusion[dir]->operator()(p_rht) += vel_diffusion[dir]->operator()(p);
										vel_diffusion[dir]->operator()(p) = 0.;
									}
								}
							}
						}
					}
				}
			}
		}
   }
}




//---------------------------------------------------------------------------
double
DS_NavierStokes:: compute_un_component ( size_t const& comp,
                                          size_t const& i,
                                          size_t const& j,
                                          size_t const& k,
                                          size_t const& dir,
                                          size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: compute_un_component" ) ;

	double value = 0.;

	if (StencilCorrection == "CutCell") {
		value = compute_un_component_FV(comp,i,j,k,dir,level);
	} else {
		value = compute_un_component_FD(comp,i,j,k,dir,level);
	}

	return(value);

}




//---------------------------------------------------------------------------
std::tuple<double,double>
DS_NavierStokes:: compute_first_derivative ( size_t const& comp,
	                                          size_t const& i,
	                                          size_t const& j,
	                                          size_t const& k,
	                                          size_t const& dir,
	                                          size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: compute_first_derivative" ) ;

   size_t_array2D* intersect_vector = (is_solids) ?
                           allrigidbodies->get_intersect_vector_on_grid(UF)
									: 0;
   doubleArray2D* intersect_distance = (is_solids) ?
                           allrigidbodies->get_intersect_distance_on_grid(UF)
									: 0;
   doubleArray2D* intersect_fieldVal = (is_solids) ?
                           allrigidbodies->get_intersect_fieldValue_on_grid(UF)
									: 0;

   size_t p = UF->DOF_local_number(i,j,k,comp);
	double value = 0., dC = 1.;
	double fc = UF->DOF_value( i, j, k, comp, level );
	double dxl_act = 1., dxr_act = 1.;
	double dxl = 1., dxr = 1.;

   if (dir == 0) {
		dC = UF->get_cell_size(i, comp, 0);
		double xc = UF->get_DOF_coordinate( i,comp, 0 );
		double xl = xc, xr = xc;
		double fr = fc, fl = fc;

      if (UF->DOF_in_domain( (int)i+1, (int)j, (int)k, comp)) {
         fr = UF->DOF_value( i+1, j, k, comp, level );
         xr = UF->get_DOF_coordinate( i+1,comp, 0 );
      }

      if (UF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp)) {
         fl = UF->DOF_value( i-1, j, k, comp, level ) ;
         xl = UF->get_DOF_coordinate( i-1, comp, 0 ) ;
      }

		dxl_act = xc - xl;
		dxr_act = xr - xc;
		dxl = dxl_act;
		dxr = dxr_act;

      if (is_solids) {
         size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
         if (void_frac->operator()(p,0) == 0) {
            if (intersect_vector->operator()(p,2*dir+0) == 1) {
               fl = intersect_fieldVal->operator()(p,2*dir+0);
               xl = xc - intersect_distance->operator()(p,2*dir+0);
					dC = ((xc - xl) < 0.5*dC ) ? 0.5 * dC + (xc - xl) : dC;
					dxl = xc - xl;
            }
            if (intersect_vector->operator()(p,2*dir+1) == 1) {
               fr = intersect_fieldVal->operator()(p,2*dir+1);
               xr = xc + intersect_distance->operator()(p,2*dir+1);
					dC = ((xr - xc) < 0.5*dC ) ? 0.5 * dC + (xr - xc) : dC;
					dxr = xr - xc;
            }
         }
      }

		// Derivative
      if (UF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp)
         && UF->DOF_in_domain( (int)i+1, (int)j, (int)k, comp)) {
			if (dxl != dxl_act) {
				if (dxl < 0.5 * dxl_act) {
					value = (fr - fl) / (dxr + dxl);
				} else {
					value = 1./(dxl_act + dxr_act)
							* ((dxl_act/dxr_act)*(fr-fc) - (dxr_act/dxl)*(fl-fc));
				}
			} else if (dxr != dxr_act) {
				if (dxr < 0.5 * dxr_act) {
					value = (fr - fl) / (dxr + dxl);
				} else {
					value = 1./(dxl_act + dxr_act)
							* ((dxl_act/dxr)*(fr-fc) - (dxr_act/dxl_act)*(fl-fc));
				}
			} else {
				value = 1./(dxl_act + dxr_act)
						* ((dxl_act/dxr_act)*(fr-fc) - (dxr_act/dxl_act)*(fl-fc));
			}
      } else if (UF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp)) {
         value = - (fl - fc)/dxl;
      } else if (UF->DOF_in_domain( (int)i+1, (int)j, (int)k, comp)) {
         value = (fr - fc)/dxr;
		} else {
			value = 0.;
		}

	} else if (dir == 1) {
		dC = UF->get_cell_size(j, comp, 1);
		double xc = UF->get_DOF_coordinate( j,comp, 1 );
		double xl = xc, xr = xc;
		double fr = fc, fl = fc;

      if (UF->DOF_in_domain((int)i, (int)j+1, (int)k, comp)) {
         fr = UF->DOF_value( i, j+1, k, comp, level );
         xr = UF->get_DOF_coordinate( j+1,comp, 1 );
      }

      if (UF->DOF_in_domain((int)i, (int)j-1, (int)k, comp)) {
         fl = UF->DOF_value( i, j-1, k, comp, level ) ;
         xl = UF->get_DOF_coordinate( j-1, comp, 1 ) ;
      }

		dxl_act = xc - xl;
		dxr_act = xr - xc;
		dxl = dxl_act;
		dxr = dxr_act;

      if (is_solids) {
         size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
         if (void_frac->operator()(p,0) == 0) {
            if (intersect_vector->operator()(p,2*dir+0) == 1) {
               fl = intersect_fieldVal->operator()(p,2*dir+0);
               xl = xc - intersect_distance->operator()(p,2*dir+0);
					dC = ((xc - xl) < 0.5*dC ) ? 0.5 * dC + (xc - xl) : dC;
					dxl = xc - xl;
            }
            if (intersect_vector->operator()(p,2*dir+1) == 1) {
               fr = intersect_fieldVal->operator()(p,2*dir+1);
               xr = xc + intersect_distance->operator()(p,2*dir+1);
					dC = ((xr - xc) < 0.5*dC ) ? 0.5 * dC + (xr - xc) : dC;
					dxr = xr - xc;
            }
         }
      }

		// Derivative
      if (UF->DOF_in_domain((int)i, (int)j-1, (int)k, comp)
         && UF->DOF_in_domain((int)i, (int)j+1, (int)k, comp)) {
				if (dxl != dxl_act) {
					if (dxl < 0.5 * dxl_act) {
						value = (fr - fl) / (dxr + dxl);
					} else {
						value = 1./(dxl_act + dxr_act)
								* ((dxl_act/dxr_act)*(fr-fc) - (dxr_act/dxl)*(fl-fc));
					}
				} else if (dxr != dxr_act) {
					if (dxr < 0.5 * dxr_act) {
						value = (fr - fl) / (dxr + dxl);
					} else {
						value = 1./(dxl_act + dxr_act)
								* ((dxl_act/dxr)*(fr-fc) - (dxr_act/dxl_act)*(fl-fc));
					}
				} else {
					value = 1./(dxl_act + dxr_act)
							* ((dxl_act/dxr_act)*(fr-fc) - (dxr_act/dxl_act)*(fl-fc));
				}
      } else if(UF->DOF_in_domain((int)i, (int)j-1, (int)k, comp)) {
         value = - (fl - fc)/dxl;
      } else if(UF->DOF_in_domain((int)i, (int)j+1, (int)k, comp)) {
         value = (fr - fc)/dxr;
		} else {
			value = 0.;
		}

	} else if (dir == 2) {
		dC = UF->get_cell_size(k, comp, 2);
		double xc = UF->get_DOF_coordinate( k,comp, 2 );
		double xl = xc, xr = xc;
		double fr = fc, fl = fc;

      if (UF->DOF_in_domain((int)i, (int)j, (int)k+1, comp)) {
         fr = UF->DOF_value( i, j, k+1, comp, level );
         xr = UF->get_DOF_coordinate( k+1,comp, 2 );
      }

      if (UF->DOF_in_domain((int)i, (int)j, (int)k-1, comp)) {
         fl = UF->DOF_value( i, j, k-1, comp, level ) ;
         xl = UF->get_DOF_coordinate( k-1, comp, 2 ) ;
      }

		dxl_act = xc - xl;
		dxr_act = xr - xc;
		dxl = dxl_act;
		dxr = dxr_act;

      if (is_solids) {
         size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
         if (void_frac->operator()(p,0) == 0) {
            if (intersect_vector->operator()(p,2*dir+0) == 1) {
               fl = intersect_fieldVal->operator()(p,2*dir+0);
               xl = xc - intersect_distance->operator()(p,2*dir+0);
					dC = ((xc - xl) < 0.5*dC ) ? 0.5 * dC + (xc - xl) : dC;
					dxl = xc - xl;
            }
            if (intersect_vector->operator()(p,2*dir+1) == 1) {
               fr = intersect_fieldVal->operator()(p,2*dir+1);
               xr = xc + intersect_distance->operator()(p,2*dir+1);
					dC = ((xr - xc) < 0.5*dC ) ? 0.5 * dC + (xr - xc) : dC;
					dxr = xr - xc;
            }
         }
      }

		// Derivative
      if (UF->DOF_in_domain((int)i, (int)j, (int)k-1, comp)
         && UF->DOF_in_domain((int)i, (int)j, (int)k+1, comp)) {
				if (dxl != dxl_act) {
					if (dxl < 0.5 * dxl_act) {
						value = (fr - fl) / (dxr + dxl);
					} else {
						value = 1./(dxl_act + dxr_act)
								* ((dxl_act/dxr_act)*(fr-fc) - (dxr_act/dxl)*(fl-fc));
					}
				} else if (dxr != dxr_act) {
					if (dxr < 0.5 * dxr_act) {
						value = (fr - fl) / (dxr + dxl);
					} else {
						value = 1./(dxl_act + dxr_act)
								* ((dxl_act/dxr)*(fr-fc) - (dxr_act/dxl_act)*(fl-fc));
					}
				} else {
					value = 1./(dxl_act + dxr_act)
							* ((dxl_act/dxr_act)*(fr-fc) - (dxr_act/dxl_act)*(fl-fc));
				}
      } else if(UF->DOF_in_domain((int)i, (int)j, (int)k-1, comp)) {
         value = - (fl - fc)/dxl;
      } else if(UF->DOF_in_domain((int)i, (int)j, (int)k+1, comp)) {
         value = (fr - fc)/dxr;
		} else {
			value = 0.;
		}
	}

   return(std::make_tuple(value,dC));

}




//---------------------------------------------------------------------------
double
DS_NavierStokes:: compute_un_component_FD ( size_t const& comp,
                                          size_t const& i,
                                          size_t const& j,
                                          size_t const& k,
                                          size_t const& dir,
                                          size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: compute_un_component_FD" ) ;

   double xhr=1.,xhl=1.,xright=0.,xleft=0.,yhr=1.,yhl=1.,yright=0.,yleft=0.;
   double zhr=1.,zhl=1.,zright=0.,zleft=0., value=0.;

   size_t_array2D* intersect_vector = (is_solids) ?
                           allrigidbodies->get_intersect_vector_on_grid(UF)
									: 0;
   doubleArray2D* intersect_distance = (is_solids) ?
                           allrigidbodies->get_intersect_distance_on_grid(UF)
									: 0;
   doubleArray2D* intersect_fieldVal = (is_solids) ?
                           allrigidbodies->get_intersect_fieldValue_on_grid(UF)
									: 0;

   size_t p = UF->DOF_local_number(i,j,k,comp);
	double dxC = UF->get_cell_size(i, comp, 0);
	double dyC = UF->get_cell_size(j, comp, 1);
	double dzC = (dim == 3) ? UF->get_cell_size(k, comp, 2) : 1.;
	// double pre_factor = 1.;

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
            size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
            if (void_frac->operator()(p,0) == 0) {
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
            && UF->DOF_in_domain( (int)i+1, (int)j, (int)k, comp)) {
            value = xright/xhr - xleft/xhl;
				// pre_factor = (2.* dxC) / (xhl + xhr);
         } else if (UF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp)) {
            value = - xleft/xhl;
				// pre_factor = dxC / xhl;
         } else {
            value = xright/xhr;
				// pre_factor = dxC/ xhr;
			}
			// value *= pre_factor * dyC * dzC;
			value *= dyC * dzC;
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
            size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
            if (void_frac->operator()(p,0) == 0) {
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
            && UF->DOF_in_domain((int)i, (int)j+1, (int)k, comp)) {
            value = yright/yhr - yleft/yhl;
				// pre_factor = (2. * dyC) / (yhr + yhl);
         } else if(UF->DOF_in_domain((int)i, (int)j-1, (int)k, comp)) {
            value = - yleft/yhl;
				// pre_factor = dyC / yhl;
         } else {
            value = yright/yhr;
				// pre_factor = dyC / yhr;
			}
			// value *= pre_factor * dxC * dzC;
			value *= dxC * dzC;
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
            size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
            if (void_frac->operator()(p,0) == 0) {
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
            && UF->DOF_in_domain((int)i, (int)j, (int)k+1, comp)) {
            value = zright/zhr - zleft/zhl;
				// pre_factor = (2. * dzC) / (zhl + zhr);
         } else if(UF->DOF_in_domain((int)i, (int)j, (int)k-1, comp)) {
            value = - zleft/zhl;
				// pre_factor = dzC / zhl;
         } else {
            value = zright/zhr;
				// pre_factor = dzC / zhr;
			}
			// value *= pre_factor * dxC * dyC;
			value *= dxC * dyC;
         break;
   }

   return(value);

}




//---------------------------------------------------------------------------
double
DS_NavierStokes:: compute_un_component_FV ( size_t const& comp,
                                       	  size_t const& i,
                                            size_t const& j,
                                            size_t const& k,
                                            size_t const& dir,
                                            size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: compute_un_component_FV" ) ;

	intVector* ownerID = allrigidbodies->get_CC_ownerID(UF);
	doubleArray3D* CC_face_cen = allrigidbodies->get_CC_face_centroid(UF);
	doubleArray2D* CC_face_frac = allrigidbodies->get_CC_face_fraction(UF);
	// size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);

	size_t p = UF->DOF_local_number(i,j,k,comp);

	geomVector xL(3), xR(3);
	xL(0) = CC_face_cen->operator()(p,2*dir+0,0);
	xL(1) = CC_face_cen->operator()(p,2*dir+0,1);
	xL(2) = CC_face_cen->operator()(p,2*dir+0,2);

	xR(0) = CC_face_cen->operator()(p,2*dir+1,0);
	xR(1) = CC_face_cen->operator()(p,2*dir+1,1);
	xR(2) = CC_face_cen->operator()(p,2*dir+1,2);

	double value = 0.;

	// if ((ownerID->operator()(p) != -1) && (void_frac->operator()(p) == 0)) {
	if (ownerID->operator()(p) != -1) {
		double fleft = 0.;
		if (CC_face_frac->operator()(p,2*dir+0) != 0.) {
			fleft = allrigidbodies->calculate_diffusive_flux(p, comp,dir, xL,ownerID->operator()(p),level)
					* CC_face_frac->operator()(p,2*dir+0);
		}
		double fright = 0.;
		if (CC_face_frac->operator()(p,2*dir+1) != 0.) {
		   fright = allrigidbodies->calculate_diffusive_flux(p, comp,dir, xR,ownerID->operator()(p),level)
					 * CC_face_frac->operator()(p,2*dir+1);
		}

		double RBflux = allrigidbodies->calculate_diffusive_flux_fromRB(UF,p,comp,dir,level);

		value = (fright - fleft + RBflux);
	} else {
		value = compute_un_component_FD(comp, i, j, k, dir, level);
	}

   return(value);

}




//---------------------------------------------------------------------------
double
DS_NavierStokes:: velocity_local_rhs ( size_t const& j
                                      , size_t const& k
                                      , double const& gamma
                                      , FV_TimeIterator const* t_it
                                      , size_t const& comp
                                      , size_t const& dir)
//---------------------------------------------------------------------------
{

   MAC_LABEL("DS_NavierStokes:: velocity_local_rhs" ) ;

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
   size_t_array2D* intersect_vector = (is_solids) ?
                           allrigidbodies->get_intersect_vector_on_grid(UF)
									: 0;
   doubleArray2D* intersect_distance = (is_solids) ?
                           allrigidbodies->get_intersect_distance_on_grid(UF)
									: 0;
   doubleArray2D* intersect_fieldVal = (is_solids) ?
                           allrigidbodies->get_intersect_fieldValue_on_grid(UF)
									: 0;

   for (size_t i=min_unknown_index(dir);i<=max_unknown_index(dir);++i) {
      if (dir == 0) {
         ii = i; jj = j; kk = k; level = 0;
      } else if (dir == 1) {
         ii = j; jj = i; kk = k; level = 3;
      } else if (dir == 2) {
         ii = j; jj = k; kk = i; level = 4;
      }

		double dxC = UF->get_cell_size(ii, comp, 0);
		double dyC = UF->get_cell_size(jj, comp, 1);
		double dzC = (dim == 3) ? UF->get_cell_size(kk, comp, 2) : 1.;
		double area = 0.;

		if (dir == 0) {
			area = dyC * dzC;
      } else if (dir == 1) {
			area = dxC * dzC;
      } else if (dir == 2) {
			area = dxC * dyC;
      }

      size_t pos = i - min_unknown_index(dir);
      size_t p = UF->DOF_local_number(ii,jj,kk,comp);
		double dC = UF->get_cell_size(i,comp,dir);

      double value = vel_diffusion[dir]->operator()(p)/area;

      if (is_solids) {
			if ((intersect_vector->operator()(p,2*dir+0) == 1)
			 && (intersect_vector->operator()(p,2*dir+1) == 1)) {
				 // double pre_factor = 2 * dC / (intersect_distance->operator()(p,2*dir+0)
			 		// 								  + intersect_distance->operator()(p,2*dir+1));
				 value = value - (intersect_fieldVal->operator()(p,2*dir+0)
                            /intersect_distance->operator()(p,2*dir+0)
								   + intersect_fieldVal->operator()(p,2*dir+1)
                            /intersect_distance->operator()(p,2*dir+1));// * pre_factor;
			} else if (intersect_vector->operator()(p,2*dir+0) == 1) {
				// double pre_factor = 2 * dC / (dC + intersect_distance->operator()(p,2*dir+0));
            value = value - intersect_fieldVal->operator()(p,2*dir+0)
                           /intersect_distance->operator()(p,2*dir+0);// * pre_factor;
         } else if (intersect_vector->operator()(p,2*dir+1) == 1) {
				// double pre_factor = 2 * dC / (dC + intersect_distance->operator()(p,2*dir+1));
            value = value - intersect_fieldVal->operator()(p,2*dir+1)
                           /intersect_distance->operator()(p,2*dir+1);// * pre_factor;
         }
      }


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
DS_NavierStokes:: unpack_compute_ue_pack(size_t const& comp
                                        , size_t const& dir
                                        , size_t const& p
                                        , size_t const& field)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: unpack_compute_ue_pack" ) ;

   LocalVector* VEC = GLOBAL_EQ->get_VEC(field);

   size_t nb_interface_unknowns = VEC[dir].T[comp]->nb_rows();

   VEC[dir].T[comp]->set(0.);
   VEC[dir].interface_T[comp]->set(0.);

   if (is_periodic[field][dir])
      VEC[dir].T[comp]->set_item(nb_ranks_comm_i[dir]-1,
                     first_pass[field][dir].send[comp][0][3*p]);

   VEC[dir].T[comp]->set_item(0,
                     first_pass[field][dir].send[comp][0][3*p+1]);
   VEC[dir].interface_T[comp]->set_item(0,
                     first_pass[field][dir].send[comp][0][3*p+2]);


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

   for (size_t i = 1; i < (size_t)nb_ranks_comm_i[dir]; ++i) {
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
DS_NavierStokes::DS_interface_unknown_solver( LA_SeqVector* interface_rhs
                                             , size_t const& comp
                                             , size_t const& dir
                                             , size_t const& field
                                             , size_t const& r_index )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokesSystem:: DS_interface_unknown_solver" ) ;

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
DS_NavierStokes:: unpack_ue(size_t const& comp
                           , double * received_data
                           , size_t const& dir
                           , size_t const& p
                           , size_t const& field)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: unpack_ue" ) ;

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
DS_NavierStokes:: solve_interface_unknowns ( FV_DiscreteField* FF
                                            , double const& gamma
                                            , FV_TimeIterator const* t_it
                                            , size_t const& comp
                                            , size_t const& dir
													  	  , size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: solve_interface_unknowns" ) ;

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
   size_t local_min_j = (dir == 0) ? min_unknown_index(1) : min_unknown_index(0);
	size_t local_max_j = (dir == 0) ? max_unknown_index(1) : max_unknown_index(0);
   size_t local_min_k = (dim == 2) ? 0 :
							  ((dir == 2) ? min_unknown_index(1) : min_unknown_index(2));
   size_t local_max_k = (dim == 2) ? 0 :
							  ((dir == 2) ? max_unknown_index(1) : max_unknown_index(2));

   size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(field,dir,comp);

   // Send and receive the data first pass
   if ( rank_in_i[dir] == 0 ) {
      if (nb_ranks_comm_i[dir] != 1) {
         for (size_t i = 1; i < (size_t)nb_ranks_comm_i[dir]; ++i) {
            // Receive the data
            static MPI_Status status;
            MPI_Recv( first_pass[field][dir].receive[comp][i],
                (int) first_pass[field][dir].size[comp],
                MPI_DOUBLE, (int) i, 0, DS_Comm_i[dir], &status ) ;
         }
      }

      for (size_t j = local_min_j; j <= local_max_j; j++) {
         for (size_t k = local_min_k; k <= local_max_k; k++) {

            size_t p = row_index->operator()(j,k);

            unpack_compute_ue_pack(comp,dir,p,field);

  	         // Need to have the original rhs function
            // assembled for corrosponding j,k pair
            assemble_local_rhs(j,k,gamma,t_it,comp,dir,field);

            // Setup RHS = fi - Aie*xe for solving ui
            A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp]
                                             ,VEC[dir].local_T[comp],-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_NavierStokes_solver(FF,j,k,min_unknown_index(dir)
                                                            ,comp,dir,p,level);
         }
      }

   } else {
      // Send the packed data to master
      MPI_Send( first_pass[field][dir].send[comp][0],
          (int) first_pass[field][dir].size[comp],
          MPI_DOUBLE, 0, 0, DS_Comm_i[dir] ) ;
   }

   // Send the data from master iff multi processor are used
   if (nb_ranks_comm_i[dir] != 1) {
      if ( rank_in_i[dir] == 0 ) {
         for (size_t i = 1; i < (size_t)nb_ranks_comm_i[dir]; ++i) {
            MPI_Send( second_pass[field][dir].send[comp][i],
                (int) second_pass[field][dir].size[comp],
                MPI_DOUBLE,(int) i, 0, DS_Comm_i[dir] ) ;
         }
      } else {
         // Receive the data
         static MPI_Status status ;
         MPI_Recv( second_pass[field][dir].receive[comp][0],
             (int) second_pass[field][dir].size[comp],
             MPI_DOUBLE, 0, 0, DS_Comm_i[dir], &status ) ;

         // Solve the system of equations in each proc
         for (size_t j = local_min_j; j <= local_max_j; j++) {
            for (size_t k = local_min_k; k <= local_max_k; k++) {
               size_t p = row_index->operator()(j,k);

               unpack_ue(comp,second_pass[field][dir].receive[comp][0]
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
                                                            ,comp,dir,p,level);
            }
         }
      }
   }
}




//---------------------------------------------------------------------------
double
DS_NavierStokes:: compute_p_component ( size_t const& comp
                                       , size_t const& i
                                       , size_t const& j
                                       , size_t const& k)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: compute_p_component" ) ;
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
DS_NavierStokes:: compute_adv_component ( FV_TimeIterator const* t_it,
														 size_t const& comp,
                                           size_t const& i,
                                           size_t const& j,
                                           size_t const& k)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: compute_adv_component" ) ;
   double ugradu = 0., value = 0.;

	double dxC = UF->get_cell_size(i,comp,0) ;
	double dyC = UF->get_cell_size(j,comp,1) ;
	double dzC = (dim == 3) ? UF->get_cell_size(k,comp,2) : 1.;
	doubleArray2D* divergence = GLOBAL_EQ->get_node_divergence(1);
	size_t p = UF->DOF_local_number(i, j, k, comp);

   if ( AdvectionScheme == "TVD" ) {
      ugradu = assemble_advection_TVD(1,rho,1,i,j,k,comp)
             - rho * UF->DOF_value(i,j,k,comp,1)
				  		 * divergence->operator()(p,0)
						 * dxC * dyC * dzC;
   } else if ( AdvectionScheme == "Upwind" ) {
      ugradu = assemble_advection_Upwind(1,rho,1,i,j,k,comp)
				 - rho * UF->DOF_value(i,j,k,comp,1)
				       * divergence->operator()(p,0)
					    * dxC * dyC * dzC;
	} else if ( AdvectionScheme == "CenteredFD" ) {
		// No divergence correction requried as we calculate the u div(u) instead of div(uu)
		ugradu = assemble_advection_Centered_FD(rho,i,j,k,comp,1);
	} else if ( AdvectionScheme == "Centered" && StencilCorrection == "FD" ) {
 		ugradu = assemble_advection_Centered(1,rho,1,i,j,k,comp)
				 - rho * UF->DOF_value(i,j,k,comp,1)
				       * divergence->operator()(p,0)
					    * dxC * dyC * dzC;
	} else if ( AdvectionScheme == "Centered" && StencilCorrection == "CutCell" ) {
		intVector* ownerID = allrigidbodies->get_CC_ownerID(UF);
		doubleArray2D* CC_vol = allrigidbodies->get_CC_cell_volume(UF);
		if (ownerID->operator()(p) != -1) {
			ugradu = assemble_advection_Centered_CutCell(t_it,rho,i,j,k,comp,1);
		} else {
			ugradu = assemble_advection_Centered(1,rho,1,i,j,k,comp);
		}
		ugradu -= rho * UF->DOF_value(i,j,k,comp,1)
						  * divergence->operator()(p,0)
						  * CC_vol->operator()(p,0);
   }

   if (( AdvectionTimeAccuracy == 1 )) {
      value = ugradu;
   } else {
      value = 1.5*ugradu - 0.5*UF->DOF_value(i,j,k,comp,2);
      UF->set_DOF_value(i,j,k,comp,2,ugradu);
   }

   return(value);
}

//----------------------------------------------------------------------
double
DS_NavierStokes:: divergence_of_U_noCorrection( size_t const& i
				                                  , size_t const& j
				                                  , size_t const& k
				                                  , size_t const& component
				                                  , size_t const& level)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: divergence_of_U_noCorrection" );

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

	if (is_solids) {
		size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
		size_t p = UF->DOF_local_number(i, j, k, component);
		if (void_frac->operator()(p,0) != 0) return(flux);
	}

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
      flux = ((vt - vb) * dxC + (ur - ul) * dyC) / dxC / dyC;
   } else if (dim == 3) {
      flux = ((vt - vb) * dxC * dzC
			  + (ur - ul) * dyC * dzC
			  + (wf - wb) * dxC * dyC) / dxC / dyC / dzC;
	}

   return ( flux );
}




//---------------------------------------------------------------------------
double
DS_NavierStokes:: get_current_mean_flow_speed()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokesSystem:: get_current_mean_flow_speed" ) ;

   size_t comp = UF->primary_grid()->get_periodic_flow_direction();

   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   // Get local min and max indices
   for (size_t l = 0; l < dim; ++l) {
     min_unknown_index(l) = UF->get_min_index_unknown_handled_by_proc( comp,l );
     max_unknown_index(l) = UF->get_max_index_unknown_handled_by_proc( comp,l );
   }

   double value = 0.;
   double volume = 0.;

   for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
      double dx = UF->get_cell_size( i, comp, 0 );
      for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
         double dy = UF->get_cell_size( j, comp, 1 );
         for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
            double dz = (dim == 3) ? UF->get_cell_size( k, comp, 2 ) : 1.;
            value += UF->DOF_value(i, j, k, comp, 0) * dx * dy * dz;
            volume += dx * dy * dz;
         }
      }
   }

   value = macCOMM->sum( value ) ;
   volume = macCOMM->sum( volume ) ;

   return(value/volume);
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: predicted_pressure_drop (FV_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: predicted_pressure_drop" ) ;

	size_t p_flow_dir = UF->primary_grid()->get_periodic_flow_direction();
	double Lx = UF->primary_grid()->get_main_domain_max_coordinate( 0 )
				 - UF->primary_grid()->get_main_domain_min_coordinate( 0 );
	double Ly = UF->primary_grid()->get_main_domain_max_coordinate( 1 )
				 - UF->primary_grid()->get_main_domain_min_coordinate( 1 );
   double Lz = (dim == 3) ? UF->primary_grid()->get_main_domain_max_coordinate( 2 )
				 				  - UF->primary_grid()->get_main_domain_min_coordinate( 2 )
								  : 1.;

   double cross_sec_area = Lx * Ly * Lz;

	if (p_flow_dir == 0) {
		cross_sec_area = Ly * Lz;
	} else if (p_flow_dir == 1) {
		cross_sec_area = Lx * Lz;
	} else if (p_flow_dir == 2) {
		cross_sec_area = Lx * Ly;
	}

	double Um = get_current_mean_flow_speed();
	double Qc = cross_sec_area * Um;
	double Qset = UF->primary_grid()->get_periodic_flow_rate();

	if ((fabs(Qc - Qset) / Qset > 1e-5) && (t_it->iteration_number() % 1 == 0)) {
      if ((Qc / Qset) > 1.) {
         exceed = true;
			// pressure_drop -= controller->calculate(Qset,Qc,t_it->time_step());
	      if (turn == false) {
	         pressure_drop *= 0.5;
	         if ((Qc - Qold) < 0.) turn = true;
         } else {
            if ((Qc / Qold) > 1.) {
               pressure_drop *= (Qset / Qc);
            } else {
					// Threshold controls the oscillation
               if (fabs(Qc - Qold) / Qc < 1e-6)
                  pressure_drop *= (Qset / Qc);
            }
         }
      } else {
         if (exceed == true) {
				// pressure_drop -= controller->calculate(Qset,Qc,t_it->time_step());
            if ((Qc / Qold) < 1.) {
               pressure_drop *= (Qset / Qc);
            } else {
					// Threshold controls the oscillation
               if (fabs(Qc - Qold) / Qc < 1e-6)
                  pressure_drop *= (Qset / Qc);
            }
         }
      }
	}
   Qold = Qc;

	string fileName = "./DS_results/flow_pressure_history.csv" ;
	std::ofstream MyFile;

	if (macCOMM->rank() == 0) {
		MyFile.open( fileName.c_str(), std::ios::app ) ;
      MyFile << t_it->time() << " " << MAC::doubleToString( ios::scientific, 6, Qc)
									  << " " << MAC::doubleToString( ios::scientific, 6, pressure_drop)
									  << " " << exceed
									  << " " << turn << endl;
   }

}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: assemble_DS_un_at_rhs ( FV_TimeIterator const* t_it,
                                           double const& gamma)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: assemble_DS_un_at_rhs" ) ;

   // Assemble the diffusive components of velocity once in each iteration
   assemble_velocity_diffusion_terms ( );
	compute_velocity_divergence(UF);
	assemble_velocity_advection_terms ( t_it );

   double bodyterm=0.;
   size_t cpp = 10;

   // Periodic pressure gradient
   if ( UF->primary_grid()->is_periodic_flow() ) {
      cpp = UF->primary_grid()->get_periodic_flow_direction() ;
		if (UF->primary_grid()->is_periodic_pressure_drop() &&
	      !UF->primary_grid()->is_periodic_flow_rate()) {
	      bodyterm = UF->primary_grid()->get_periodic_pressure_drop() /
	               ( UF->primary_grid()->get_main_domain_max_coordinate( cpp )
	               - UF->primary_grid()->get_main_domain_min_coordinate( cpp ) ) ;
		} else if (UF->primary_grid()->is_periodic_flow_rate()) {
			predicted_pressure_drop(t_it);
			bodyterm = pressure_drop /
	               ( UF->primary_grid()->get_main_domain_max_coordinate( cpp )
	               - UF->primary_grid()->get_main_domain_min_coordinate( cpp ) ) ;
		}
   }

   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);

	doubleArray2D* CC_vol = (is_solids) ?
							allrigidbodies->get_CC_cell_volume(UF) : 0;
	vector<doubleVector*> advection = GLOBAL_EQ->get_velocity_advection();
   vector<doubleVector*> vel_diffusion = GLOBAL_EQ->get_velocity_diffusion();
   size_t_array2D* void_frac = (is_solids) ?
							allrigidbodies->get_void_fraction_on_grid(UF) : 0;

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
					// Cell volume
					double cellV = dxC * dyC * dzC;
					if ((StencilCorrection == "CutCell")
					 && (CC_vol->operator()(p,0) >= FluxRedistThres*cellV)) {
					 	cellV = CC_vol->operator()(p,0);
					}

               // Dxx for un
               double xvalue = vel_diffusion[0]->operator()(p);
               // Dyy for un
               double yvalue = vel_diffusion[1]->operator()(p);
               // Dzz for un
               double zvalue = (dim == 3) ? vel_diffusion[2]->operator()(p) : 0;
               // Pressure contribution
               double pvalue = compute_p_component(comp,i,j,k);
               // Advection contribution
               double adv_value = advection[0]->operator()(p);
               //
               if (is_solids) {
                  if (void_frac->operator()(p,0) != 0) {
                     pvalue = 0.;
                  }
               }

               double rhs = gamma*(xvalue + yvalue + zvalue)
                          - pvalue - adv_value
                          + (UF->DOF_value( i, j, k, comp, 1 )*cellV*rho)
                                                         /(t_it -> time_step());


               if ( cpp==comp ) rhs += - bodyterm*cellV;

               if (is_solids) {
                  if (void_frac->operator()(p,0) != 0) {
                     if ( cpp==comp ) rhs += bodyterm*cellV;
                  }
               }

               UF->set_DOF_value( i, j, k, comp, 0,
                                  rhs*(t_it -> time_step())/(cellV*rho));
            }
         }
      }
   }
}

//---------------------------------------------------------------------------
void
DS_NavierStokes:: Solve_i_in_jk ( FV_DiscreteField* FF
                                 , FV_TimeIterator const* t_it
                                 , size_t const& dir_i
                                 , size_t const& dir_j
                                 , size_t const& dir_k
                                 , double const& gamma
											, size_t const& level)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_NavierStokes:: Solve_i_in_jk" ) ;

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
     size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(field,dir_i,comp);

     // Solve in i
     if ((nb_ranks_comm_i[dir_i]>1)||(is_periodic[field][dir_i] == 1)) {
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j){
           for (size_t k=local_min_k; k <= local_max_k; ++k) {
              size_t r_index = row_index->operator()(j,k);
              // Assemble fi and return fe for each proc locally
              double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir_i,field);
              // Calculate Aei*ui in each proc locally
              compute_Aei_ui(A,VEC,comp,dir_i,r_index);
              // Pack Aei_ui and fe for sending it to master
              data_packing (FF,r_index,fe,comp,dir_i);
           }
        }
        solve_interface_unknowns ( FF, gamma, t_it, comp, dir_i,level);

     } else if (is_periodic[field][dir_i] == 0) {
        // Serial mode with non-periodic condition
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j){
           for (size_t k=local_min_k; k <= local_max_k; ++k) {
              size_t r_index = row_index->operator()(j,k);
              assemble_local_rhs(j,k,gamma,t_it,comp,dir_i,field);
              GLOBAL_EQ->DS_NavierStokes_solver(FF,j,k,min_unknown_index(dir_i)
                                                  ,comp,dir_i,r_index,level);
           }
        }
     }
  }
}

//---------------------------------------------------------------------------
void
DS_NavierStokes:: data_packing ( FV_DiscreteField const* FF
                                , size_t const& p
                                , double const& fe
                                , size_t const& comp
                                , size_t const& dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: data_packing" ) ;

   size_t field = (FF == PF) ? 0 : 1 ;

   LocalVector* VEC = GLOBAL_EQ->get_VEC(field) ;
   double *packed_data = first_pass[field][dir].send[comp][0];

   if (rank_in_i[dir] == 0) {
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_periodic[field][dir])
          packed_data[3*p+0] =
                        VEC[dir].T[comp]->item(nb_ranks_comm_i[dir]-1);
      else
          packed_data[3*p+0] = 0;

      packed_data[3*p+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);

   } else if (rank_in_i[dir] == nb_ranks_comm_i[dir]-1) {
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_periodic[field][dir])
          packed_data[3*p+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);
      else
          packed_data[3*p+1]=0;

      packed_data[3*p+0] = VEC[dir].T[comp]->item(rank_in_i[dir]-1);

   } else {
      packed_data[3*p+0] = VEC[dir].T[comp]->item(rank_in_i[dir]-1);
      packed_data[3*p+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);
   }

   // Send the fe values and 0 for last proc
   packed_data[3*p+2] = fe;

}

//---------------------------------------------------------------------------
void
DS_NavierStokes:: compute_Aei_ui (struct TDMatrix* arr
                                 , struct LocalVector* VEC
                                 , size_t const& comp
                                 , size_t const& dir
                                 , size_t const& r_index)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: compute_Aei_ui" ) ;
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
DS_NavierStokes:: assemble_local_rhs ( size_t const& j
                                     , size_t const& k
                                     , double const& gamma
                                     , FV_TimeIterator const* t_it
                                     , size_t const& comp
                                     , size_t const& dir
                                     , size_t const& field )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: assemble_local_rhs" ) ;
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
DS_NavierStokes:: NS_velocity_update ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: NS_velocity_update" ) ;

   double gamma=mu;

   assemble_DS_un_at_rhs(t_it,gamma);
   // Update gamma based for invidual direction
   gamma = mu/2.0;

   Solve_i_in_jk(UF,t_it,0,1,2,gamma,3);
   // Synchronize the velocity field
	UF->synchronize( 3 );
   if (is_solids) initialize_grid_nodes_on_rigidbody({3});

	size_t level = (dim == 2) ? 0 : 4 ;
   Solve_i_in_jk(UF,t_it,1,0,2,gamma,level);
	// Synchronize the velocity field
	UF->synchronize( level );
   if (is_solids) initialize_grid_nodes_on_rigidbody({level});

   if (dim == 3) {
      Solve_i_in_jk(UF,t_it,2,0,1,gamma,0);
		// Synchronize the velocity field
		UF->synchronize( 0 );
      if (is_solids) initialize_grid_nodes_on_rigidbody({0});
   }

	// Compute velocity change over the time step
	double velocity_time_change = compute_DS_velocity_change()/t_it->time_step();
	if ( my_rank == is_master ) cout << "velocity change = " <<
		MAC::doubleToString( ios::scientific, 5, velocity_time_change ) << endl;

	compute_velocity_divergence(PF);

}




//---------------------------------------------------------------------------
void
DS_NavierStokes::initialize_grid_nodes_on_rigidbody( vector<size_t> const& list )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_NavierStokes::initialize_grid_nodes_on_rigidbody" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  // Vector for solid presence
  size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);

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
              if (void_frac->operator()(p,0) != 0) {
                 size_t par_id = void_frac->operator()(p,0) - 1;
                 geomVector rb_vel = allrigidbodies->rigid_body_velocity(par_id,pt);
					  // geomVector rb_vel = allrigidbodies->rigid_body_GC_velocity(par_id);
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
DS_NavierStokes:: assemble_velocity_gradients (class doubleVector& grad
                                              , size_t const& i
                                              , size_t const& j
                                              , size_t const& k
                                              , size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: assemble_velocity_gradients" ) ;

   FV_SHIFT_TRIPLET shift = PF->shift_staggeredToCentered() ;

   size_t comp = 0;

   size_t_array2D* intersect_vector = (is_solids) ?
                     allrigidbodies->get_intersect_vector_on_grid(PF) : 0;
   doubleArray2D* intersect_distance = (is_solids) ?
                     allrigidbodies->get_intersect_distance_on_grid(PF) : 0;
   doubleArray2D* intersect_fieldVal = (is_solids) ?
                     allrigidbodies->get_intersect_fieldValue_on_grid(PF) : 0;
   size_t_array2D* void_frac = (is_solids) ?
							allrigidbodies->get_void_fraction_on_grid(PF) : 0;

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

   if (is_solids && StencilCorrection == "FD") {
      if (void_frac->operator()(p,0) == 0) {
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

   if (is_solids && StencilCorrection == "FD") {
      if (void_frac->operator()(p,0) == 0) {
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

      if (is_solids && StencilCorrection == "FD") {
         if (void_frac->operator()(p,0) == 0) {
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
DS_NavierStokes:: calculate_velocity_divergence_FD ( size_t const& i,
                                                   size_t const& j,
                                                   size_t const& k,
                                                   size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: calculate_velocity_divergence" ) ;

   doubleVector grad(3,0);

   double beta = assemble_velocity_gradients(grad,i,j,k,level);

   double value = beta*(grad(0) + grad(1) + grad(2));

   return(value);
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: assemble_velocity_advection_terms ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: assemble_velocity_advection_terms" ) ;

	vector<doubleVector*> advection = GLOBAL_EQ->get_velocity_advection();
	size_t_array2D* void_frac = (is_solids) ?
							allrigidbodies->get_void_fraction_on_grid(UF) : 0;

	size_t_vector min_unknown_index(3,0);
	size_t_vector max_unknown_index(3,0);

	for (size_t comp=0;comp<nb_comps[1];comp++) {
		// Get local min and max indices
		for (size_t l=0;l<dim;++l) {
			min_unknown_index(l) =
								UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
			max_unknown_index(l) =
								UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
		}

		// Compute the advection value and store
		for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
			for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
				for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
					size_t p = UF->DOF_local_number(i,j,k,comp);
					// Advection contribution
					advection[0]->operator()(p) = compute_adv_component(t_it,comp,i,j,k);
				}
			}
		}

		if (is_solids && StencilCorrection == "CutCell") {
			size_t UF_UNK_MAX = UF->nb_local_unknowns();
			intVector* ownerID = allrigidbodies->get_CC_ownerID(UF);
			// Flux redistribution
			for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
				for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
					for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
						size_t p = UF->DOF_local_number(i,j,k,comp);
						if (ownerID->operator()(p) != -1) {
							size_t p_lft = UF->DOF_local_number(i-1,j,k,comp);
							size_t p_rht = UF->DOF_local_number(i+1,j,k,comp);
							size_t p_bot = UF->DOF_local_number(i,j-1,k,comp);
							size_t p_top = UF->DOF_local_number(i,j+1,k,comp);
							size_t p_bhd = (dim == 2) ? 0 : UF->DOF_local_number(i,j,k-1,comp);
							size_t p_frt = (dim == 2) ? 0 : UF->DOF_local_number(i,j,k+1,comp);

							vector<double> wht = allrigidbodies->
							flux_redistribution_factor(UF,i,j,k,comp,FluxRedistThres);

						   // Flux redistribution
						   double sum = wht[0] + wht[1] + wht[2] + wht[3] + wht[4] + wht[5];

							if (sum > 0.) {
								if (p_lft <= UF_UNK_MAX)
									advection[0]->operator()(p_lft) +=
													wht[0]/sum * advection[0]->operator()(p);
								if (p_rht <= UF_UNK_MAX)
									advection[0]->operator()(p_rht) +=
													wht[1]/sum * advection[0]->operator()(p);
								if (p_bot <= UF_UNK_MAX)
									advection[0]->operator()(p_bot) +=
													wht[2]/sum * advection[0]->operator()(p);
								if (p_top <= UF_UNK_MAX)
									advection[0]->operator()(p_top) +=
													wht[3]/sum * advection[0]->operator()(p);
								if (dim == 3) {
									if (p_bhd <= UF_UNK_MAX)
										advection[0]->operator()(p_bhd) +=
													wht[4]/sum * advection[0]->operator()(p);
									if (p_frt <= UF_UNK_MAX)
										advection[0]->operator()(p_frt) +=
													wht[5]/sum * advection[0]->operator()(p);
								}
								advection[0]->operator()(p) = 0.;
							} else if ((sum == 0.) && (void_frac->operator()(p,0) != 0)) {
								advection[0]->operator()(p) = 0.;
							}
						}
					}
				}
			}
		}

	}

}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: compute_velocity_divergence ( FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: compute_velocity_divergence" ) ;

	size_t field = (FF == PF) ? 0 : 1;
	size_t UNK_MAX = FF->nb_local_unknowns();
	doubleArray2D* divergence = GLOBAL_EQ->get_node_divergence(field);
	doubleArray2D* CC_vol = (is_solids) ? allrigidbodies->get_CC_cell_volume(FF) : 0;
	intVector* ownerID = (is_solids) ? allrigidbodies->get_CC_ownerID(FF) : 0;

	for (size_t comp = 0; comp < nb_comps[field]; comp++) {
		// Get local min and max indices
		size_t_vector min_unknown_index(3,0);
		size_t_vector max_unknown_index(3,0);
		for (size_t l = 0; l < dim; ++l) {
			min_unknown_index(l) =
					FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
			max_unknown_index(l) =
					FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
		}

		for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
			double dx = FF->get_cell_size( i, comp, 0 );
			for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
				double dy = FF->get_cell_size( j, comp, 1 );
				for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
					double dz = (dim == 3) ? FF->get_cell_size( k, comp, 2 ) : 1.;
					size_t p = FF->DOF_local_number(i,j,k,comp);
					if (is_solids) {
						if (ownerID->operator()(p) == -1) {
							double vel_div = (FF == UF) ?
												divergence_of_U_noCorrection(i,j,k,comp,0)
											 : calculate_velocity_divergence_FD(i,j,k,0);
				         divergence->operator()(p,0) = vel_div * dx * dy * dz;
						} else if (ownerID->operator()(p) != -1) {
						   double rht_flux = allrigidbodies->velocity_flux(FF,p,0,0,1,0);
						   double lft_flux = allrigidbodies->velocity_flux(FF,p,0,0,0,0);
						   double top_flux = allrigidbodies->velocity_flux(FF,p,1,1,1,0);
						   double bot_flux = allrigidbodies->velocity_flux(FF,p,1,1,0,0);
						   double fnt_flux = (dim == 2) ? 0.
						 	       			 : allrigidbodies->velocity_flux(FF,p,2,2,1,0);
						   double bhd_flux = (dim == 2) ? 0.
						 	    				 : allrigidbodies->velocity_flux(FF,p,2,2,0,0);

							double xvalue = (rht_flux - lft_flux);
							double yvalue = (top_flux - bot_flux);
							double zvalue = (fnt_flux - bhd_flux);
							double RB_flux = allrigidbodies->calculate_velocity_flux_fromRB(FF,p);
							double dFV = 0.;//(CC_vol->operator()(p,0) - CC_vol->operator()(p,1))/0.001;

							double value = xvalue + yvalue + zvalue + RB_flux + dFV;

							// In case of nodes at the boundary of proc, especially
							// nodes not handled by the proc
							value = isnan(value) ? 0. : value;

				         divergence->operator()(p,0) = value;
						}
					} else {
						double vel_div = (FF == UF) ?
											divergence_of_U_noCorrection(i,j,k,comp,0)
										 : calculate_velocity_divergence_FD(i,j,k,0);
						divergence->operator()(p,0) = vel_div * dx * dy * dz;
					}
				}
			}
		}

		size_t_array2D* void_frac = (is_solids) ? allrigidbodies->get_void_fraction_on_grid(FF) : 0;

		// Flux redistribution
		if ((StencilCorrection == "CutCell") && (FF == UF)) {
			for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
				for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
					for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
						size_t p = FF->DOF_local_number(i,j,k,comp);
						if (ownerID->operator()(p) != -1) {
							size_t p_lft = FF->DOF_local_number(i-1,j,k,comp);
							size_t p_rht = FF->DOF_local_number(i+1,j,k,comp);
							size_t p_bot = FF->DOF_local_number(i,j-1,k,comp);
							size_t p_top = FF->DOF_local_number(i,j+1,k,comp);
							size_t p_bhd = (dim == 2) ? 0 : FF->DOF_local_number(i,j,k-1,comp);
							size_t p_frt = (dim == 2) ? 0 : FF->DOF_local_number(i,j,k+1,comp);

							vector<double> wht = allrigidbodies
							->flux_redistribution_factor(FF,i,j,k,comp,FluxRedistThres);

							// Flux redistribution
							double sum = wht[0] + wht[1] + wht[2] + wht[3] + wht[4] + wht[5];

							if (sum > 0.) {
								if (p_lft <= UNK_MAX)
									divergence->operator()(p_lft,0) +=
													wht[0]/sum * divergence->operator()(p,0);
								if (p_rht <= UNK_MAX)
									divergence->operator()(p_rht,0) +=
													wht[1]/sum * divergence->operator()(p,0);
								if (p_bot <= UNK_MAX)
									divergence->operator()(p_bot,0) +=
													wht[2]/sum * divergence->operator()(p,0);
								if (p_top <= UNK_MAX)
									divergence->operator()(p_top,0) +=
													wht[3]/sum * divergence->operator()(p,0);
								if (dim == 3) {
									if (p_bhd <= UNK_MAX)
										divergence->operator()(p_bhd,0) +=
													wht[4]/sum * divergence->operator()(p,0);
									if (p_frt <= UNK_MAX)
										divergence->operator()(p_frt,0) +=
													wht[5]/sum * divergence->operator()(p,0);
								}
								divergence->operator()(p,0) = 0.;
							} else if ((sum == 0.) && (void_frac->operator()(p,0) != 0)) {
								divergence->operator()(p,0) = 0.;
							}
						}
					}
				}
			}
		}

		// Converting flux to divergence
		for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
			double dx = FF->get_cell_size( i, comp, 0 );
			for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
				double dy = FF->get_cell_size( j, comp, 1 );
				for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
					double dz = (dim == 3) ? FF->get_cell_size( k, comp, 2 ) : 1.;
					size_t p = FF->DOF_local_number(i, j, k, comp);
					double volume = dx * dy * dz;
					if (is_solids && (FF == UF) && (StencilCorrection == "CutCell")
					 && (CC_vol->operator()(p,0) >= FluxRedistThres*volume)) {
					 	volume = CC_vol->operator()(p,0);
					}
					// if (is_solids && (void_frac->operator()(p) == 0) && (StencilCorrection == "CutCell"))
					// 	volume = CC_vol->operator()(p,0);
					divergence->operator()(p,0) = divergence->operator()(p,0)/volume;
				}
			}
		}
	}
}




//---------------------------------------------------------------------------
double
DS_NavierStokes:: pressure_local_rhs ( size_t const& j
                                      , size_t const& k
                                      , FV_TimeIterator const* t_it
                                      , size_t const& dir)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: pressure_local_rhs" ) ;
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
   doubleArray2D* divergence = GLOBAL_EQ->get_node_divergence(0);

   for (size_t i=min_unknown_index(dir);i<=max_unknown_index(dir);++i) {
      double dx = PF->get_cell_size( i, 0, dir );
      if (dir == 0) {
			size_t p = PF->DOF_local_number(i,j,k,0);
         double vel_div = divergence->operator()(p,0);
         value = -(rho*vel_div*dx)/(t_it -> time_step());
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
DS_NavierStokes:: correct_pressure_1st_layer_solid (size_t const& level )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_NavierStokes:: correct_pressure_1st_layer_solid" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  NodeProp node = GLOBAL_EQ->get_node_property(0,0);
  size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(PF);

  size_t comp = 0;

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
           size_t p = PF->DOF_local_number(i,j,k,comp);
           node.bound_cell->set_item(p,0.);
           if (void_frac->operator()(p,0) != 0) {
              double value = 0., count = 0.;
              size_t p1 = PF->DOF_local_number(i+1,j,k,comp);
              size_t p2 = PF->DOF_local_number(i-1,j,k,comp);
              size_t p3 = PF->DOF_local_number(i,j+1,k,comp);
              size_t p4 = PF->DOF_local_number(i,j-1,k,comp);
              if (void_frac->operator()(p1,0) == 0) {
                 node.bound_cell->set_item(p,1.);
                 value += PF->DOF_value( i+1, j, k, comp, level );
                 count += 1.;
              }
              if (void_frac->operator()(p2,0) == 0) {
                 node.bound_cell->set_item(p,1.);
                 value += PF->DOF_value( i-1, j, k, comp, level );
                 count += 1.;
              }
              if (void_frac->operator()(p3,0) == 0) {
                 node.bound_cell->set_item(p,1.);
                 value += PF->DOF_value( i, j+1, k, comp, level );
                 count += 1.;
              }
              if (void_frac->operator()(p4,0) == 0) {
                 node.bound_cell->set_item(p,1.);
                 value += PF->DOF_value( i, j-1, k, comp, level );
                 count += 1.;
              }

              if (dim == 3) {
                 size_t p5 = PF->DOF_local_number(i,j,k+1,comp);
                 size_t p6 = PF->DOF_local_number(i,j,k-1,comp);

                 if (void_frac->operator()(p5,0) == 0) {
                    node.bound_cell->set_item(p,1.);
                    value += PF->DOF_value( i, j, k+1, comp, level );
                    count += 1.;
                 }
                 if (void_frac->operator()(p6,0) == 0) {
                    node.bound_cell->set_item(p,1.);
                    value += PF->DOF_value( i, j, k-1, comp, level );
                    count += 1.;
                 }
              }

              if (count != 0.) value = value/count;
				  PF->set_DOF_value( i, j, k, comp, level, value);
           }
        }
     }
  }

  // Synchronize the pressure field
  PF->synchronize( level );
}

//---------------------------------------------------------------------------
void
DS_NavierStokes:: correct_pressure_2nd_layer_solid (size_t const& level )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_NavierStokes:: correct_pressure_2nd_layer_solid" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  NodeProp node = GLOBAL_EQ->get_node_property(0,0);
  size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(PF);

  size_t comp = 0;

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1);j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2);k <= max_unknown_index(2); ++k) {
           size_t p = PF->DOF_local_number(i,j,k,comp);
           if ((void_frac->operator()(p,0) != 0)
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
				  PF->set_DOF_value( i, j, k, comp, level, value);
           }
        }
     }
  }

  // Synchronize the pressure field
  PF->synchronize( level );
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: correct_mean_pressure (size_t const& level )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_NavierStokes:: correct_mean_pressure" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  size_t_array2D* void_frac = (is_solids) ?
  						allrigidbodies->get_void_fraction_on_grid(PF) : 0;
  size_t comp = 0;

  double mean=0.,nb_global_unknown=0.;

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
           if (is_solids) {
              size_t p = PF->DOF_local_number(i,j,k,comp);
              if (void_frac->operator()(p,0) == 0) {
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
              if (void_frac->operator()(p,0) == 0) {
                 double value = PF->DOF_value( i, j, k, comp, level );
					  PF->set_DOF_value( i, j, k, comp, level, value-mean);
              }
           }
        }
     }
  }

  // Synchronize the pressure field
  PF->synchronize( level );

}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: NS_pressure_update ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_NavierStokes:: NS_pressure_update" ) ;

  double gamma=mu/2.0;

  Solve_i_in_jk (PF,t_it,0,1,2,gamma,1);
  // Synchronize the pressure field
  PF->synchronize( 1 );

  Solve_i_in_jk (PF,t_it,1,0,2,gamma,1);
  // Synchronize the pressure field
  PF->synchronize( 1 );

  if (dim == 3) {
     Solve_i_in_jk (PF,t_it,2,0,1,gamma,1);
	  // Synchronize the pressure field
	  PF->synchronize( 1 );
  }

  if (PF->all_BCs_nonDirichlet(0)) {
     correct_mean_pressure(1);
  }

  if (is_solids) {
	  allrigidbodies->extrapolate_pressure_inside_RB(PF,1);
	  // correct_pressure_1st_layer_solid(1);
	  // correct_pressure_2nd_layer_solid(1);
  }
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: NS_final_step ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: NS_final_step" ) ;

   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);

   doubleArray2D* divergencePF = GLOBAL_EQ->get_node_divergence(0);
	doubleArray2D* divergenceUF = GLOBAL_EQ->get_node_divergence(1);
	doubleArray2D* CC_volPF = (is_solids) ? allrigidbodies->get_CC_cell_volume(PF) : 0;
	doubleArray2D* CC_volUF = (is_solids) ? allrigidbodies->get_CC_cell_volume(UF) : 0;
	size_t_array2D* void_fracPF = (is_solids) ? allrigidbodies->get_void_fraction_on_grid(PF) : 0;
	size_t_array2D* void_fracUF = (is_solids) ? allrigidbodies->get_void_fraction_on_grid(UF) : 0;

   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
      max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
   }

   for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
      for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
         for (size_t k = min_unknown_index(2);k <= max_unknown_index(2); ++k) {
            size_t p = PF->DOF_local_number(i,j,k,0);
            double vel_div0 = divergencePF->operator()(p,0);
            double vel_div1 = divergencePF->operator()(p,1);

            // Assemble the bodyterm
            double value = PF->DOF_value( i, j, k, 0, 0 )
                         + PF->DOF_value( i, j, k, 0, 1 )
                         - 0.5*kai*mu*(vel_div0+vel_div1);
				PF->set_DOF_value( i, j, k, 0, 0, value);
         }
      }
   }

   // Synchronize the pressure field
	PF->synchronize( 0 );

   if (PF->all_BCs_nonDirichlet(0)) {
      correct_mean_pressure(0);
   }

   // Store the divergence to be used in the next time iteration
   for (size_t i = 0; i < PF->nb_local_unknowns() ; i++) {
      divergencePF->operator()(i,1) = divergencePF->operator()(i,0);
		if (is_solids) {
			CC_volPF->operator()(i,1) = CC_volPF->operator()(i,0);
			void_fracPF->operator()(i,1) = void_fracPF->operator()(i,0);
		}
	}

	// Store the divergence to be used in the next time iteration
	for (size_t i = 0; i < UF->nb_local_unknowns() ; i++) {
		divergenceUF->operator()(i,1) = divergenceUF->operator()(i,0);
		if (is_solids) {
			CC_volUF->operator()(i,1) = CC_volUF->operator()(i,0);
			void_fracUF->operator()(i,1) = void_fracUF->operator()(i,0);
		}
	}

   if (is_solids) {
		allrigidbodies->extrapolate_pressure_inside_RB(PF,0);
		// correct_pressure_1st_layer_solid(0);
		// correct_pressure_2nd_layer_solid(0);
   }
   // Propagate values to the boundaries depending on BC conditions
   PF->set_neumann_DOF_values();
}




//----------------------------------------------------------------------
void
DS_NavierStokes::write_output_field(FV_DiscreteField const* FF)
//----------------------------------------------------------------------
{
  ofstream outputFile ;

  std::ostringstream os2;
  os2 << "./DS_results/intersection_data_" << my_rank << ".csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());

  size_t field = (FF == UF) ? 1 : 0 ;

  size_t i,j,k;
  outputFile << "x,y,z,id,void_frac,div,left,lv"
             << ",right,rv"
             << ",bottom,bov"
             << ",top,tv"
             << ",behind,bev"
             << ",front,fv"
				 << ",fresh"
				 << ",field"
				 << endl;

  size_t_vector min_index(dim,0);
  size_t_vector max_index(dim,0);

  boolVector* fresh_node = (is_solids) ?
                    allrigidbodies->get_fresh_nodes(FF) : 0;
  size_t_array2D* void_frac = (is_solids) ?
                    allrigidbodies->get_void_fraction_on_grid(FF) : 0;
  size_t_array2D* intersect_vector = (is_solids) ?
                    allrigidbodies->get_intersect_vector_on_grid(FF) : 0;
  doubleArray2D* intersect_distance = (is_solids) ?
                    allrigidbodies->get_intersect_distance_on_grid(FF) : 0;

  doubleArray2D* divergence = GLOBAL_EQ->get_node_divergence(field);

  for (size_t comp=0;comp<nb_comps[field];comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_index(2);
        local_max_k = max_index(2);
     }

     for (i=min_index(0);i<=max_index(0);++i) {
        double xC = FF->get_DOF_coordinate( i, comp, 0 ) ;
        for (j=min_index(1);j<=max_index(1);++j) {
           double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;
           for (k=local_min_k;k<=local_max_k;++k) {
              double zC = 0.;
              if (dim == 3) zC = FF->get_DOF_coordinate( k, comp, 2 ) ;
              size_t p = FF->DOF_local_number(i,j,k,comp);

              outputFile << xC << "," << yC << "," << zC
				  << "," << p
              << "," << void_frac->operator()(p,0)
				  << "," << divergence->operator()(p,0)
              << "," << intersect_vector->operator()(p,0)
              << "," << intersect_distance->operator()(p,0)
              << "," << intersect_vector->operator()(p,1)
              << "," << intersect_distance->operator()(p,1)
              << "," << intersect_vector->operator()(p,2)
              << "," << intersect_distance->operator()(p,2)
              << "," << intersect_vector->operator()(p,3)
              << "," << intersect_distance->operator()(p,3)
              << "," << intersect_vector->operator()(p,4)
              << "," << intersect_distance->operator()(p,4)
              << "," << intersect_vector->operator()(p,5)
              << "," << intersect_distance->operator()(p,5)
				  << "," << fresh_node->operator()(p)
				  << "," << FF->DOF_value(i,j,k,comp,0)
              << endl;
           }
        }
     }
  }
  outputFile.close();
}




//----------------------------------------------------------------------
void
DS_NavierStokes::output_L2norm_divergence( )
//----------------------------------------------------------------------
{
  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  double div_velocity = 0.;
  double max_divu = 0.;

  doubleArray2D* divergence = GLOBAL_EQ->get_node_divergence(0);

  for (size_t l=0;l<dim;++l) {
    min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
    max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
     for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
        double dx = PF->get_cell_size( i, 0, 0 );
        double dy = PF->get_cell_size( j, 0, 1 );
        if (dim == 2) {
           size_t p = PF->DOF_local_number(i,j,0,0);
           double vel_div = divergence->operator()(p,0);
           max_divu = MAC::max( MAC::abs(vel_div), max_divu );
           div_velocity += vel_div * vel_div * ( dx * dy );
        } else {
           for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
              double dz = PF->get_cell_size( k, 0, 2 );
              size_t p = PF->DOF_local_number(i,j,k,0);
              double vel_div = divergence->operator()(p,0);
              max_divu = MAC::max( MAC::abs(vel_div), max_divu );
              div_velocity += vel_div * vel_div * ( dx * dy * dz );
           }
        }
     }
  }

  div_velocity = macCOMM->sum( div_velocity ) ;
  div_velocity = MAC::sqrt( div_velocity );
  max_divu = macCOMM->max( max_divu ) ;
  if ( my_rank == is_master )
    MAC::out() << "Norm L2 div(u) = "
               << MAC::doubleToString( ios::scientific, 12, div_velocity )
               << " Max div(u) = "
               << MAC::doubleToString( ios::scientific, 12, max_divu )
               << endl;

  if ( my_rank == is_master ) {
     string fileName = "./DS_results/max_divergence.csv" ;
     ofstream MyFile( fileName.c_str(), ios::app ) ;
     MyFile << div_velocity << "," << max_divu << endl;
     MyFile.close( ) ;
  }
}




//----------------------------------------------------------------------
void
DS_NavierStokes::output_L2norm_pressure( size_t const& level )
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

  size_t_array2D* void_frac = (is_solids) ?
  									  allrigidbodies->get_void_fraction_on_grid(PF) : 0;

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
              if (void_frac->operator()(p,0) == 0) {
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
double
DS_NavierStokes:: compute_DS_velocity_change( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: compute_DS_velocity_change" ) ;

	size_t_vector min_unknown_index(3,0);
	size_t_vector max_unknown_index(3,0);

	double sum_sq_U=0.,sum_sq_dU=0.;

	for (size_t comp = 0; comp < nb_comps[1]; comp++) {
		for (size_t l = 0; l < dim; ++l) {
			min_unknown_index(l) =
							 UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
			max_unknown_index(l) =
							 UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
		}


		for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
			for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
				for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
					sum_sq_U += pow(UF->DOF_value( i, j, k, comp, 0 ),2.);
				   sum_sq_dU += pow(UF->DOF_value( i, j, k, comp, 0 )
									   - UF->DOF_value( i, j, k, comp, 1 ),2.);

				}
			}
		}
	}

	sum_sq_U = macCOMM->sum(sum_sq_U);
	sum_sq_dU = macCOMM->sum(sum_sq_dU);

   return ( MAC::sqrt(sum_sq_dU/sum_sq_U) ) ;
}




//----------------------------------------------------------------------
void
DS_NavierStokes::output_L2norm_velocity( size_t const& level )
//----------------------------------------------------------------------
{
  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  size_t_array2D* void_frac = (is_solids) ?
  									allrigidbodies->get_void_fraction_on_grid(UF) : 0;

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
                 if (void_frac->operator()(p,0) == 0) {
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
DS_NavierStokes:: create_DS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: create_DS_subcommunicators" ) ;

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
DS_NavierStokes:: processor_splitting ( int const& color
                                       , int const& key
                                       , size_t const& dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: processor_splitting" ) ;

   MPI_Comm_split(MPI_COMM_WORLD, color, key, &DS_Comm_i[dir]);
   MPI_Comm_size( DS_Comm_i[dir], &nb_ranks_comm_i[dir] ) ;
   MPI_Comm_rank( DS_Comm_i[dir], &rank_in_i[dir] ) ;

}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: allocate_mpi_variables (FV_DiscreteField const* FF)
//---------------------------------------------------------------------------
{

   size_t field = (FF == PF) ? 0 : 1 ;

   for (size_t dir = 0; dir < dim; dir++) {
      first_pass[field][dir].size = new size_t [nb_comps[field]];
      second_pass[field][dir].size = new size_t [nb_comps[field]];
		data_for_S[field][dir].size = new size_t [nb_comps[field]];
      for (size_t comp = 0; comp < nb_comps[field]; comp++) {
         // Get local min and max indices
         size_t_vector min_unknown_index(dim,0);
         size_t_vector max_unknown_index(dim,0);
         for (size_t l=0;l<dim;++l) {
            min_unknown_index(l) =
                        FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
            max_unknown_index(l) =
                        FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
         }

			size_t local_min_j = (dir == 0) ? min_unknown_index(1)
													  : min_unknown_index(0);
			size_t local_max_j = (dir == 0) ? max_unknown_index(1)
													  : max_unknown_index(0);
			size_t local_min_k = (dim == 2) ? 0 :
									  ((dir == 2) ? min_unknown_index(1)
									  				  : min_unknown_index(2));
			size_t local_max_k = (dim == 2) ? 0 :
									  ((dir == 2) ? max_unknown_index(1)
									  				  : max_unknown_index(2));

         size_t local_length_j = (local_max_j-local_min_j+1);
         size_t local_length_k = (local_max_k-local_min_k+1);

         if (dim != 3) {
            first_pass[field][dir].size[comp] = 3*local_length_j;
            second_pass[field][dir].size[comp] = 2*local_length_j;
				data_for_S[field][dir].size[comp] = is_periodic[field][dir] ?
						 (size_t)(pow(nb_ranks_comm_i[dir],2) + 1)*local_length_j
					  : (size_t)(pow(nb_ranks_comm_i[dir]-1,2) + 1)*local_length_j ;
         } else if (dim == 3) {
            first_pass[field][dir].size[comp] =
                                    3*local_length_j*local_length_k;
            second_pass[field][dir].size[comp] =
	                                 2*local_length_j*local_length_k;
				data_for_S[field][dir].size[comp] = is_periodic[field][dir] ?
			    (size_t)(pow(nb_ranks_comm_i[dir],2) + 1)*local_length_j*local_length_k
			  : (size_t)(pow(nb_ranks_comm_i[dir]-1,2) + 1)*local_length_j*local_length_k ;
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
		data_for_S[field][dir].send =
							new double** [nb_comps[field]];
		data_for_S[field][dir].receive =
							new double** [nb_comps[field]];
      for (size_t comp = 0; comp < nb_comps[field]; comp++) {
         first_pass[field][dir].send[comp] =
                        new double* [1];
         first_pass[field][dir].receive[comp] =
                        new double* [nb_ranks_comm_i[dir]];
         second_pass[field][dir].send[comp] =
                        new double* [nb_ranks_comm_i[dir]];
         second_pass[field][dir].receive[comp] =
                        new double* [1];
			data_for_S[field][dir].send[comp] =
	                     new double* [1];
	      data_for_S[field][dir].receive[comp] =
	                     new double* [1];
			data_for_S[field][dir].send[comp][0] =
								new double[data_for_S[field][dir].size[comp]];
			data_for_S[field][dir].receive[comp][0] =
								new double[data_for_S[field][dir].size[comp]];
			first_pass[field][dir].send[comp][0] =
								new double[first_pass[field][dir].size[comp]];
			second_pass[field][dir].receive[comp][0] =
								new double[second_pass[field][dir].size[comp]];
         for (size_t i = 0; i < (size_t) nb_ranks_comm_i[dir]; i++) {
            first_pass[field][dir].receive[comp][i] =
                           new double[first_pass[field][dir].size[comp]];
            second_pass[field][dir].send[comp][i] =
                           new double[second_pass[field][dir].size[comp]];
         }
      }
   }
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: deallocate_mpi_variables ()
//---------------------------------------------------------------------------
{
   // Array declarations
   for (size_t field = 0; field < 2; field++) {
      for (size_t dir = 0; dir < dim; dir++) {
         for (size_t comp = 0; comp < nb_comps[field]; comp++) {
				delete [] data_for_S[field][dir].send[comp][0];
				delete [] data_for_S[field][dir].receive[comp][0];
				delete [] first_pass[field][dir].send[comp][0];
				delete [] second_pass[field][dir].receive[comp][0];
            for (size_t i = 0; i < (size_t) nb_ranks_comm_i[dir]; i++) {
               delete [] first_pass[field][dir].receive[comp][i];
               delete [] second_pass[field][dir].send[comp][i];
            }
            delete [] first_pass[field][dir].send[comp];
            delete [] first_pass[field][dir].receive[comp];
            delete [] second_pass[field][dir].send[comp];
            delete [] second_pass[field][dir].receive[comp];
				delete [] data_for_S[field][dir].send[comp];
            delete [] data_for_S[field][dir].receive[comp];
         }
         delete [] first_pass[field][dir].send;
         delete [] first_pass[field][dir].receive;
         delete [] second_pass[field][dir].send;
         delete [] second_pass[field][dir].receive;
			delete [] data_for_S[field][dir].send;
         delete [] data_for_S[field][dir].receive;
         delete [] first_pass[field][dir].size;
         delete [] second_pass[field][dir].size;
			delete [] data_for_S[field][dir].size;
      }
   }
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: free_DS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: free_DS_subcommunicators" ) ;
}




//---------------------------------------------------------------------------
void
DS_NavierStokes:: set_translation_vector()
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: set_translation_vector" ) ;

   MVQ_translation_vector.resize( dim );
   translation_direction = UF->primary_grid()->get_translation_direction() ;
   MVQ_translation_vector( translation_direction ) =
        UF->primary_grid()->get_translation_magnitude() ;
}




//---------------------------------------------------------------------------
void
DS_NavierStokes::build_links_translation()
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: build_links_translation" ) ;

   UF->create_transproj_interpolation() ;
   PF->create_transproj_interpolation() ;

}




//---------------------------------------------------------------------------
void
DS_NavierStokes::fields_projection()
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: fields_projection" ) ;

   FV_Mesh* pmesh = const_cast<FV_Mesh*>(UF->primary_grid()) ;

   pmesh->translation() ;

   UF->translation_projection( 0, 5, 0 ) ;
   UF->synchronize( 0 ) ;

   UF->translation_projection( 1, 5, 0 ) ;
   UF->synchronize( 1 ) ;

   UF->translation_projection( 2, 5, 0 ) ;
   UF->synchronize( 2 ) ;

   UF->translation_projection( 3, 5, 0 ) ;
	UF->synchronize( 3 ) ;

   UF->translation_projection( 4, 5 ) ;
	UF->synchronize( 4 ) ;

   PF->translation_projection( 0, 2, 0 ) ;
	PF->synchronize( 0 );

   PF->translation_projection( 1, 2 ) ;
	PF->synchronize( 1 );

}




//----------------------------------------------------------------------
double
DS_NavierStokes:: assemble_advection_Centered( size_t const& advecting_level,
                                                double const& coef,
                                                size_t const& advected_level,
                                                size_t const& i,
                                                size_t const& j,
                                                size_t const& k,
                                                size_t const& component )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: assemble_advection_Centered" );

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
DS_NavierStokes:: assemble_advection_Centered_FD(double const& coef,
													          size_t const& i,
                                                 size_t const& j,
                                                 size_t const& k,
                                                 size_t const& comp,
																 size_t const& level)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: assemble_advection_Centered_FD" );

	double dh = UF->primary_grid()->get_smallest_grid_size();
   double threshold = 0.0001*dh;
	FV_SHIFT_TRIPLET shift = UF->shift_staggeredToStaggered( comp );

	double dxC = UF->get_cell_size(i,comp,0) ;
	double dyC = UF->get_cell_size(j,comp,1) ;
	double dzC = (dim == 3) ? UF->get_cell_size(k,comp,2) : 1.;

	geomVector pt(3,0.);
	pt(0) = UF->get_DOF_coordinate( i, comp, 0 );
	pt(1) = UF->get_DOF_coordinate( j, comp, 1 );
	pt(2) = (dim == 3) ? UF->get_DOF_coordinate( k, comp, 2 ) : 0. ;

	size_t_vector face_vector(3,0);
	face_vector(0) = 1; face_vector(1) = 1; face_vector(2) = 0;

	doubleVector ui(3,0.);
	doubleVector dui(3,0.);
	doubleVector dC(3,1.);

	dC(0) = dxC; dC(1) = dyC; dC(2) = dzC;

	size_t p = UF->DOF_local_number(i,j,k,comp);

	size_t_array2D* intersect_vector = (is_solids) ?
						allrigidbodies->get_intersect_vector_on_grid(UF) : 0;

	bool nearRB = false;
	if (is_solids) {
		for (size_t dir = 0; dir < dim && !nearRB ; dir++) {
			nearRB = (intersect_vector->operator()(p,2*dir+0) == 1)
					|| (intersect_vector->operator()(p,2*dir+1) == 1);
		}
		size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);
		if (void_frac->operator()(p,0) != 0) return(0.);
	}

	for (size_t cpt = 0; cpt < dim; cpt++) {
		auto value = compute_first_derivative(comp,i,j,k,cpt,level);
		dui(cpt) = std::get<0>(value);
		dC(cpt) = std::get<1>(value);

		if (is_solids && nearRB) {
			if (cpt != comp) {
				size_t_vector i0(3);
				for (size_t l = 0; l < dim; l++) {
					size_t i0_temp;
					bool found = FV_Mesh::between( UF->get_DOF_coordinates_vector(cpt,l),
					pt(l) + threshold, i0_temp);
					if (found) i0(l) = i0_temp;
				}
				ui(cpt) = (dim == 2) ? allrigidbodies->
						Bilinear_interpolation(UF, cpt, &pt, i0, face_vector, {level})
												: allrigidbodies->
						Trilinear_interpolation(UF, cpt, &pt, i0, {level});
			} else {
				ui(cpt) = UF->DOF_value( i, j, k, comp, level );
			}
			// if (p == 10584) {
			// 	cout << cpt << "," << dui(cpt) << "," << dC(cpt) << "," << ui(cpt) << endl;
			// }
		} else {
			if (comp == 0) {
				if (cpt == 0) {
					ui(cpt) = UF->DOF_value( i, j, k, comp, level );
				} else if (cpt == 1) {
					double uToLe = UF->DOF_value(i+shift.i-1, j+shift.j, k, cpt, level );
			      double uToRi = UF->DOF_value(i+shift.i, j+shift.j, k, cpt, level );
					double uBoLe = UF->DOF_value(i+shift.i-1, j+shift.j-1, k, cpt, level );
			      double uBoRi = UF->DOF_value(i+shift.i, j+shift.j-1, k, cpt, level );
					if (UF->DOF_color( i, j, k, comp) == FV_BC_RIGHT) {
						ui(cpt) = 0.5 * (uToLe + uBoLe);
					} else if (UF->DOF_color( i, j, k, comp) == FV_BC_LEFT) {
						ui(cpt) = 0.5 * (uToRi + uBoRi);
					} else {
						ui(cpt) = 0.25 * (uToLe + uToRi + uBoLe + uBoRi);
					}
				} else if (cpt == 2) {
					double uFrLe = UF->DOF_value(i+shift.i-1, j, k+shift.k, cpt, level );
					double uFrRi = UF->DOF_value(i+shift.i, j, k+shift.k, cpt, level );
					double uBeLe = UF->DOF_value(i+shift.i-1, j, k+shift.k-1, cpt, level );
					double uBeRi = UF->DOF_value(i+shift.i, j, k+shift.k-1, cpt, level );
					if (UF->DOF_color( i, j, k, comp) == FV_BC_RIGHT) {
						ui(cpt) = 0.5 * (uFrLe + uBeLe);
					} else if (UF->DOF_color( i, j, k, comp) == FV_BC_LEFT) {
						ui(cpt) = 0.5 * (uFrRi + uBeRi);
					} else {
						ui(cpt) = 0.25 * (uFrLe + uFrRi + uBeLe + uBeRi);
					}
				}
			} else if (comp == 1) {
				if (cpt == 0) {
					double uToLe = UF->DOF_value(i+shift.i-1, j+shift.j, k, cpt, level );
			      double uToRi = UF->DOF_value(i+shift.i, j+shift.j, k, cpt, level );
					double uBoLe = UF->DOF_value(i+shift.i-1, j+shift.j-1, k, cpt, level );
			      double uBoRi = UF->DOF_value(i+shift.i, j+shift.j-1, k, cpt, level );
					if (UF->DOF_color( i, j, k, comp) == FV_BC_TOP) {
						ui(cpt) = 0.5 * (uBoLe + uBoRi);
					} else if (UF->DOF_color( i, j, k, comp) == FV_BC_BOTTOM) {
						ui(cpt) = 0.5 * (uToLe + uToRi);
					} else {
						ui(cpt) = 0.25 * (uToLe + uToRi + uBoLe + uBoRi);
					}
				} else if (cpt == 1) {
					ui(cpt) = UF->DOF_value( i, j, k, comp, level );
				} else if (cpt == 2) {
					double uFrBo = UF->DOF_value(i, j+shift.j-1, k+shift.k, cpt, level );
					double uFrTo = UF->DOF_value(i, j+shift.j, k+shift.k, cpt, level );
					double uBeBo = UF->DOF_value(i, j+shift.j-1, k+shift.k-1, cpt, level );
					double uBeTo = UF->DOF_value(i, j+shift.j, k+shift.k-1, cpt, level );
					if (UF->DOF_color( i, j, k, comp) == FV_BC_TOP) {
						ui(cpt) = 0.5 * (uFrBo + uBeBo);
					} else if (UF->DOF_color( i, j, k, comp) == FV_BC_BOTTOM) {
						ui(cpt) = 0.5 * (uFrTo + uBeTo);
					} else {
						ui(cpt) = 0.25 * (uFrBo + uFrTo + uBeBo + uBeTo);
					}
				}
			} else if (comp == 2) {
				if (cpt == 0) {
					double uFrLe = UF->DOF_value(i+shift.i-1, j, k+shift.k, cpt, level );
					double uFrRi = UF->DOF_value(i+shift.i, j, k+shift.k, cpt, level );
					double uBeLe = UF->DOF_value(i+shift.i-1, j, k+shift.k-1, cpt, level );
					double uBeRi = UF->DOF_value(i+shift.i, j, k+shift.k-1, cpt, level );
					if (UF->DOF_color( i, j, k, comp) == FV_BC_FRONT) {
						ui(cpt) = 0.5 * (uBeLe + uBeRi);
					} else if (UF->DOF_color( i, j, k, comp) == FV_BC_BEHIND) {
						ui(cpt) = 0.5 * (uFrLe + uFrRi);
					} else {
						ui(cpt) = 0.25 * (uFrLe + uFrRi + uBeLe + uBeRi);
					}
				} else if (cpt == 1) {
					double uFrBo = UF->DOF_value(i, j+shift.j-1, k+shift.k, cpt, level );
					double uFrTo = UF->DOF_value(i, j+shift.j, k+shift.k, cpt, level );
					double uBeBo = UF->DOF_value(i, j+shift.j-1, k+shift.k-1, cpt, level );
					double uBeTo = UF->DOF_value(i, j+shift.j, k+shift.k-1, cpt, level );
					if (UF->DOF_color( i, j, k, comp) == FV_BC_FRONT) {
						ui(cpt) = 0.5 * (uBeBo + uBeTo);
					} else if (UF->DOF_color( i, j, k, comp) == FV_BC_BEHIND) {
						ui(cpt) = 0.5 * (uFrBo + uFrTo);
					} else {
						ui(cpt) = 0.25 * (uFrBo + uFrTo + uBeBo + uBeTo);
					}
				} else if (cpt == 2) {
					ui(cpt) = UF->DOF_value( i, j, k, comp, level );
				}
			}
		}
	}

	double beta = 1.;

	beta = MAC::min(dC(0)/dxC,dC(1)/dyC);
	if (dim == 3) beta = MAC::min(beta, dC(2)/dzC);
	beta *= dxC * dyC * dzC;

	return ( coef*(ui(0)*dui(0)+ui(1)*dui(1)+ui(2)*dui(2))*beta);
}




//----------------------------------------------------------------------
double
DS_NavierStokes:: assemble_advection_Centered_CutCell(FV_TimeIterator const* t_it,
																		double const& coef,
																		size_t const& i,
                                                	   size_t const& j,
                                                		size_t const& k,
                                                		size_t const& comp,
																	 	size_t const& level)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: assemble_advection_Centered_CutCell" );

	double dh = UF->primary_grid()->get_smallest_grid_size();
   double threshold = 0.0001*dh;

	intVector color(6);
	color(0) = FV_BC_LEFT; color(1) = FV_BC_RIGHT;
	color(2) = FV_BC_BOTTOM; color(3) = FV_BC_TOP;
	if (dim == 3) {color(4) = FV_BC_BEHIND; color(5) = FV_BC_FRONT;}

	doubleArray2D* face_fraction = allrigidbodies->get_CC_face_fraction(UF);
	doubleArray3D* face_centroid = allrigidbodies->get_CC_face_centroid(UF);

	doubleVector flux(6,0.);

	// Area of the faces without even considering the RB
	doubleVector GridfaceArea(3,0.);
   double dxC = UF->get_cell_size(i,comp,0) ;
   double dyC = UF->get_cell_size(j,comp,1) ;
   double dzC = (dim == 3) ? UF->get_cell_size(k,comp,2) : 1.;
	GridfaceArea(0) = dyC * dzC ;
	GridfaceArea(1) = dxC * dzC ;
	GridfaceArea(2) = dxC * dyC ;

   double ValueC = UF->DOF_value( i, j, k, comp, level );

	size_t_vector face_vector(3,0);
	face_vector(0) = 1; face_vector(1) = 1; face_vector(2) = 0;

	size_t p = UF->DOF_local_number(i, j, k, comp);

	for (size_t dir = 0; dir < dim; dir++) {
		for (size_t side = 0; side < 2; side++) {
			if (comp == dir) {
				if ( UF->DOF_color( i, j, k, comp ) == color(2*dir+side) ) {
		         flux(2*dir+side) = ValueC * ValueC * GridfaceArea(dir);
				} else {
					if ((face_fraction->operator()(p,2*dir+side) <= GridfaceArea(dir) + threshold)
					 && (face_fraction->operator()(p,2*dir+side) != 0.) ) {
						geomVector pt(3);
						pt(0) = face_centroid->operator()(p,2*dir+side,0);
						pt(1) = face_centroid->operator()(p,2*dir+side,1);
						pt(2) = face_centroid->operator()(p,2*dir+side,2);

						size_t_vector i0(3);
						for (size_t l = 0; l < dim; l++) {
				         size_t i0_temp;
				         bool found = FV_Mesh::between( UF->get_DOF_coordinates_vector(comp,l),
				                                        pt(l) + threshold, i0_temp);
				         if (found) i0(l) = i0_temp;
						}
						double uface = (dim == 2) ? allrigidbodies->
							Bilinear_interpolation(UF, comp, &pt, i0, face_vector, {level})
														  : allrigidbodies->
						   Trilinear_interpolation(UF, comp, &pt, i0, {level});
						flux(2*dir+side) = uface * uface * face_fraction->operator()(p,2*dir+side);
					}
				}
			} else {
				if ((face_fraction->operator()(p,2*dir+side) <= GridfaceArea(dir) + threshold)
				 && (face_fraction->operator()(p,2*dir+side) != 0.) ) {
					geomVector pt(3);
					pt(0) = face_centroid->operator()(p,2*dir+side,0);
					pt(1) = face_centroid->operator()(p,2*dir+side,1);
					pt(2) = face_centroid->operator()(p,2*dir+side,2);

					size_t_vector i0(3);
					// Advector value
					for (size_t l = 0; l < dim; l++) {
						size_t i0_temp;
						bool found = FV_Mesh::between( UF->get_DOF_coordinates_vector(dir,l),
																 pt(l) + threshold, i0_temp);
						if (found) i0(l) = i0_temp;
					}
					double uAdvector = (dim == 2) ? allrigidbodies->
						Bilinear_interpolation(UF, dir, &pt, i0, face_vector, {level})
															: allrigidbodies->
						Trilinear_interpolation(UF, dir, &pt, i0, {level});
					// Advected value
					for (size_t l = 0; l < dim; l++) {
						size_t i0_temp;
						bool found = FV_Mesh::between( UF->get_DOF_coordinates_vector(comp,l),
																 pt(l) + threshold, i0_temp);
						if (found) i0(l) = i0_temp;
					}
					double uAdvected = (dim == 2) ? allrigidbodies->
						Bilinear_interpolation(UF, comp, &pt, i0, face_vector, {level})
										  					: allrigidbodies->
					  	Trilinear_interpolation(UF, comp, &pt, i0, {level});
					flux(2*dir+side) = uAdvected * uAdvector * face_fraction->operator()(p,2*dir+side);
				}
			}
		}
	}

	// RB contribution
	intVector* ownerID = allrigidbodies->get_CC_ownerID(UF);
	doubleArray2D* normalRB = allrigidbodies->get_CC_RB_normal(UF);
	doubleVector* RBarea = allrigidbodies->get_CC_RB_area(UF);
	double fsurf = 0.;
	// // Effect of volume change, based on Seo and Mittal 2011; and Arthur Ghigo CutCell
	doubleArray2D* CC_vol = allrigidbodies->get_CC_cell_volume(UF);
	double dFV = 0.;

	if ((ownerID->operator()(p) != -1) && (RBarea->operator()(p) != 0.)) {
		size_t local_parID = (size_t) ownerID->operator()(p);
		geomVector const* pgc = allrigidbodies->get_gravity_centre(local_parID);
      geomVector pt(pgc->operator()(0),
                    pgc->operator()(1),
                    pgc->operator()(2));
      geomVector rbVel = allrigidbodies->rigid_body_velocity(local_parID,pt);
		geomVector normVec(normalRB->operator()(p,0),
								 normalRB->operator()(p,1),
							    normalRB->operator()(p,2));
	   normVec *= 1./normVec.calcNorm();

		fsurf = rbVel(comp) * rbVel.operator,(normVec) * RBarea->operator()(p);

		dFV = (CC_vol->operator()(p,0) - CC_vol->operator()(p,1))
			 * (ValueC - rbVel(comp)) / t_it->time_step();
	}

   return ( coef * ((flux(3) - flux(2)) + (flux(1) - flux(0)) + (flux(5) - flux(4)) - fsurf + dFV) );
}




//----------------------------------------------------------------------
double
DS_NavierStokes:: assemble_advection_Upwind(
  size_t const& advecting_level, double const& coef, size_t const& advected_level,
  size_t const& i, size_t const& j, size_t const& k, size_t const& component )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: assemble_advection_Upwind" );

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
DS_NavierStokes:: assemble_advection_TVD(
  size_t const& advecting_level, double const& coef, size_t const& advected_level,
         size_t const& i, size_t const& j, size_t const& k, size_t const& component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: assemble_advection_TVD" );

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
