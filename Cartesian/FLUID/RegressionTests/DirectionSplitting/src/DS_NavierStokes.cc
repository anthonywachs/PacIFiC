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
   , AdvectionTimeAccuracy( fromDS.AdvectionTimeAccuracy_ )
   , rho( fromDS.rho_ )
	, b_restart ( fromDS.b_restart_ )
   , is_solids( fromDS.is_solids_ )
	, is_stressCal ( fromDS.is_stressCal_ )
	, ViscousStressOrder ( fromDS.ViscousStressOrder_ )
	, surface_cell_scale ( fromDS.surface_cell_scale_ )
	, is_surfacestressOUT ( fromDS.is_surfacestressOUT_ )
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

   // Timing routines
   if ( my_rank == is_master ) {
     SCT_insert_app("Matrix_Assembly&Initialization");
	  SCT_insert_app("Matrix_RE_Assembly&Initialization");
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

   if (is_solids) {
		// Build void frac and intersection variable
		allrigidbodies->build_solid_variables_on_fluid_grid(PF);
		allrigidbodies->build_solid_variables_on_fluid_grid(UF);
		// Compute void fraction for pressure and velocity field
		allrigidbodies->compute_void_fraction_on_grid(PF);
		allrigidbodies->compute_void_fraction_on_grid(UF);
		// Compute intersection with RB for pressure and velocity field
		allrigidbodies->compute_grid_intersection_with_rigidbody(PF);
		allrigidbodies->compute_grid_intersection_with_rigidbody(UF);

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
   // write_output_field(UF,1);
   // SCT_get_elapsed_time( "Writing CSV" );

	if (is_surfacestressOUT)
   	allrigidbodies->write_surface_discretization_for_all_RB();

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


   if ((is_par_motion) && (is_solids)) {
		if ( my_rank == is_master )
			SCT_set_start( "Matrix_RE_Assembly&Initialization" );
		// Solve equation of motion for all RB and update pos,vel
		allrigidbodies->solve_RB_equation_of_motion(t_it);
		allrigidbodies->compute_halo_zones_for_all_rigid_body();
		// Compute void fraction for pressure and velocity field
		allrigidbodies->compute_void_fraction_on_grid(PF);
		allrigidbodies->compute_void_fraction_on_grid(UF);
		// Compute intersection with RB for pressure and velocity field
		allrigidbodies->compute_grid_intersection_with_rigidbody(PF);
		allrigidbodies->compute_grid_intersection_with_rigidbody(UF);

		// Field initialization
		vector<size_t> vec{ 0, 1, 3};
		if (dim == 3) vec.push_back(4);
		initialize_grid_nodes_on_rigidbody(vec);

	   // Direction splitting
      // Assemble 1D tridiagonal matrices
      assemble_1D_matrices(PF,t_it);
      assemble_1D_matrices(UF,t_it);
		if ( my_rank == is_master )
			SCT_get_elapsed_time( "Matrix_RE_Assembly&Initialization" );
   }

	//  Projection_Translation
	if ( b_projection_translation )
		b_grid_has_been_translated_at_previous_time = false;

   // Perform matrix level operations before each time step
   GLOBAL_EQ->at_each_time_step( );

}

//---------------------------------------------------------------------------
void
DS_NavierStokes:: do_after_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: do_after_inner_iterations_stage" ) ;

   // Compute velocity change over the time step
   double velocity_time_change = GLOBAL_EQ->compute_DS_velocity_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master ) cout << "velocity change = " <<
     	MAC::doubleToString( ios::scientific, 5, velocity_time_change ) << endl;

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
double
DS_NavierStokes:: assemble_field_matrix ( FV_DiscreteField const* FF,
                                           FV_TimeIterator const* t_it,
                                           double const& gamma,
                                           size_t const& comp,
                                           size_t const& dir,
                                           size_t const& j,
                                           size_t const& k,
                                           size_t const& r_index )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: assemble_field_matrix" ) ;

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
         } else if (void_frac->operator()(p) != 0) {
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
DS_NavierStokes:: assemble_field_schur_matrix (size_t const& comp
                                              , size_t const& dir
                                              , double const& Aee_diagcoef
                                              , size_t const& field
                                              , size_t const& r_index )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: assemble_field_schur_matrix" ) ;
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
                     MPI_DOUBLE, (int)i, 0, DS_Comm_i[dir], &status ) ;

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
                                    MPI_DOUBLE, 0, 0, DS_Comm_i[dir] ) ;

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
DS_NavierStokes:: assemble_1D_matrices ( FV_DiscreteField const* FF
                                        , FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: assemble_1D_matrices" ) ;

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

         size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(field,dir,comp);

         for (size_t j = min_unknown_index(dir_j);
                    j <= max_unknown_index(dir_j);++j) {
            for (size_t k = local_min_k; k <= local_max_k; ++k) {
               size_t r_index = row_index->operator()(j,k);
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
  if ( my_rank == is_master ) SCT_set_start( "Pressure stress" );
  if (is_stressCal && (t_it->iteration_number() % stressCalFreq == 0)) {
     allrigidbodies->compute_pressure_force_and_torque_for_allRB();
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
DS_NavierStokes:: unpack_compute_ue_pack(size_t const& comp
                                        , size_t const& dir
                                        , size_t const& p
                                        , size_t const& field)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: unpack_compute_ue_pack" ) ;

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
                                            , size_t const& dir)
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
                                                            ,comp,dir,p);
         }
      }

   } else {
      // Send the packed data to master
      MPI_Send( first_pass[field][dir].send[comp][rank_in_i[dir]],
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
         MPI_Recv( second_pass[field][dir].send[comp][rank_in_i[dir]],
             (int) first_pass[field][dir].size[comp],
             MPI_DOUBLE, 0, 0, DS_Comm_i[dir], &status ) ;

         // Solve the system of equations in each proc
         for (size_t j = local_min_j; j <= local_max_j; j++) {
            for (size_t k = local_min_k; k <= local_max_k; k++) {
               size_t p = row_index->operator()(j,k);

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
DS_NavierStokes:: compute_adv_component ( size_t const& comp,
                                           size_t const& i,
                                           size_t const& j,
                                           size_t const& k)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: compute_adv_component" ) ;
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
DS_NavierStokes:: divergence_of_U( size_t const& i
                                  , size_t const& j
                                  , size_t const& k
                                  , size_t const& component
                                  , size_t const& level)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: divergence_of_U" );

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
DS_NavierStokes:: assemble_DS_un_at_rhs ( FV_TimeIterator const* t_it,
                                           double const& gamma)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: assemble_DS_un_at_rhs" ) ;

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
   size_t_vector* void_frac = (is_solids) ?
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
                  if (void_frac->operator()(p) != 0) {
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
                  if (void_frac->operator()(p) != 0) {
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
DS_NavierStokes:: Solve_i_in_jk ( FV_DiscreteField* FF
                                 , FV_TimeIterator const* t_it
                                 , size_t const& dir_i
                                 , size_t const& dir_j
                                 , size_t const& dir_k
                                 , double const& gamma)
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
              data_packing (FF,j,k,fe,comp,dir_i);
           }
        }
        solve_interface_unknowns ( FF, gamma, t_it, comp, dir_i);

     } else if (is_periodic[field][dir_i] == 0) {
        // Serial mode with non-periodic condition
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j){
           for (size_t k=local_min_k; k <= local_max_k; ++k) {
              size_t r_index = row_index->operator()(j,k);
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
DS_NavierStokes:: data_packing ( FV_DiscreteField const* FF
                                , size_t const& j
                                , size_t const& k
                                , double const& fe
                                , size_t const& comp
                                , size_t const& dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: data_packing" ) ;

   size_t field = (FF == PF) ? 0 : 1 ;

   LocalVector* VEC = GLOBAL_EQ->get_VEC(field) ;
   size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(field,dir,comp);

   double *packed_data = first_pass[field][dir].send[comp][rank_in_i[dir]];

   // Pack the data
   size_t vec_pos = row_index->operator()(j,k);

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
DS_NavierStokes::initialize_grid_nodes_on_rigidbody( vector<size_t> const& list )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_NavierStokes::initialize_grid_nodes_on_rigidbody" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  // Vector for solid presence
  size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(UF);

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
              if (void_frac->operator()(p) != 0) {
                 size_t par_id = void_frac->operator()(p) - 1;
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
   size_t_vector* void_frac = (is_solids) ?
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
DS_NavierStokes:: calculate_velocity_divergence ( size_t const& i,
                                                   size_t const& j,
                                                   size_t const& k,
                                                   size_t const& level,
                                                   FV_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_NavierStokes:: calculate_velocity_divergence" ) ;

   doubleVector grad(3,0);

   double beta = assemble_velocity_gradients(grad,i,j,k,level);

   double value = beta*(grad(0) + grad(1) + grad(2));

   return(value);
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
  size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(PF);

  size_t comp = 0;

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
           size_t p = PF->DOF_local_number(i,j,k,comp);
           node.bound_cell->set_item(p,0.);
           if (void_frac->operator()(p) != 0) {
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
  size_t_vector* void_frac = allrigidbodies->get_void_fraction_on_grid(PF);

  size_t comp = 0;

  for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
     for (size_t j = min_unknown_index(1);j <= max_unknown_index(1); ++j) {
        for (size_t k = min_unknown_index(2);k <= max_unknown_index(2); ++k) {
           size_t p = PF->DOF_local_number(i,j,k,comp);
           if ((void_frac->operator()(p) != 0)
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

  size_t_vector* void_frac = (is_solids) ?
  						allrigidbodies->get_void_fraction_on_grid(PF) : 0;
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
DS_NavierStokes:: NS_pressure_update ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_NavierStokes:: NS_pressure_update" ) ;

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
DS_NavierStokes:: NS_final_step ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: NS_final_step" ) ;

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




//----------------------------------------------------------------------
void
DS_NavierStokes::write_output_field(FV_DiscreteField const* FF)
//----------------------------------------------------------------------
{
  ofstream outputFile ;

  std::ostringstream os2;
  os2 << "./DS_results/intersection_data.csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());

  size_t field = (FF == UF) ? 1 : 0 ;

  size_t i,j,k;
  outputFile << "x,y,z,void_frac,left,lv"
             << ",right,rv"
             << ",bottom,bov"
             << ",top,tv"
             << ",behind,bev"
             << ",front,fv" << endl;

  size_t_vector min_index(dim,0);
  size_t_vector max_index(dim,0);

  size_t_vector* void_frac = (is_solids) ?
                    allrigidbodies->get_void_fraction_on_grid(FF) : 0;
  size_t_array2D* intersect_vector = (is_solids) ?
                    allrigidbodies->get_intersect_vector_on_grid(FF) : 0;
  doubleArray2D* intersect_distance = (is_solids) ?
                    allrigidbodies->get_intersect_distance_on_grid(FF) : 0;

  for (size_t comp=0;comp<nb_comps[field];comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_index(l) = 0;
        max_index(l) = FF->get_local_nb_dof( comp, l );
     }

     size_t local_min_k = 0;
     size_t local_max_k = 1;

     if (dim == 3) {
        local_min_k = min_index(2);
        local_max_k = max_index(2);
     }

     for (i=min_index(0);i<max_index(0);++i) {
        double xC = FF->get_DOF_coordinate( i, comp, 0 ) ;
        for (j=min_index(1);j<max_index(1);++j) {
           double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;
           for (k=local_min_k;k<local_max_k;++k) {
              double zC = 0.;
              if (dim == 3) zC = FF->get_DOF_coordinate( k, comp, 2 ) ;
              size_t p = FF->DOF_local_number(i,j,k,comp);

              outputFile << xC << "," << yC << "," << zC
              << "," << void_frac->operator()(p)
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

  doubleVector* divergence = GLOBAL_EQ->get_node_divergence(0);

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
           double vel_div = divergence->operator()(p);
           max_divu = MAC::max( MAC::abs(vel_div), max_divu );
           div_velocity += vel_div * ( dx * dy );
        } else {
           for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
              double dz = PF->get_cell_size( k, 0, 2 );
              size_t p = PF->DOF_local_number(i,j,k,0);
              double vel_div = divergence->operator()(p);
              max_divu = MAC::max( MAC::abs(vel_div), max_divu );
              div_velocity += vel_div * ( dx * dy * dz );
           }
        }
     }
  }

  div_velocity = macCOMM->sum( div_velocity ) ;
  max_divu = macCOMM->max( max_divu ) ;
  if ( my_rank == is_master )
    MAC::out() << "Norm L2 div(u) = "
               << MAC::doubleToString( ios::scientific, 12, div_velocity )
               << " Max div(u) = "
               << MAC::doubleToString( ios::scientific, 12, max_divu )
               << endl;
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

  size_t_vector* void_frac = (is_solids) ?
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
DS_NavierStokes::output_L2norm_velocity( size_t const& level )
//----------------------------------------------------------------------
{
  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  size_t_vector* void_frac = (is_solids) ?
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
DS_NavierStokes:: deallocate_mpi_variables ()
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
DS_NavierStokes:: synchronize_velocity_field ( size_t level )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: synchronize_velocity_field" ) ;

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
DS_NavierStokes:: synchronize_pressure_field ( size_t level )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_NavierStokes:: synchronize_pressure_field" ) ;

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
