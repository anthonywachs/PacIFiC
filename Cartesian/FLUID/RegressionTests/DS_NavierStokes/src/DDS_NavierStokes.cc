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
   , rho( 1. )
   , is_solids( false )
{
   MAC_LABEL( "DDS_NavierStokes:: DDS_NavierStokes" ) ;
   MAC_ASSERT( UF->discretization_type() == "staggered" ) ;
   MAC_ASSERT( PF->discretization_type() == "centered" ) ;
   MAC_ASSERT( UF->storage_depth() == 5 ) ;

   // Call of MAC_Communicator routine to set the rank of each proces and
   // the number of processes during execution
   pelCOMM = MAC_Exec::communicator();
   my_rank = pelCOMM->rank();
   nb_procs = pelCOMM->nb_ranks();
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

   // Read the presence of particles
   if ( exp->has_entry( "Particles" ) )
     is_solids = exp->bool_data( "Particles" ) ;

   if (is_solids) {
      Npart = exp->int_data( "NParticles" ) ;
      insertion_type = exp->string_data( "InsertionType" ) ;
      MAC_ASSERT( insertion_type == "file" ) ;
      solid_filename = exp->string_data( "Particle_FileName" ) ;
      loc_thres = exp->double_data( "Local_threshold" ) ;
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
   if ( AdvectionScheme != "Upwind" && AdvectionScheme != "TVD" )
   {
     string error_message="   - Upwind\n   - TVD";
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


   // Advection term time accuracy
   if ( exp->has_entry( "AdvectionTimeAccuracy" ) )
     AdvectionTimeAccuracy = exp->int_data( "AdvectionTimeAccuracy" );
   if ( AdvectionTimeAccuracy != 1 && AdvectionTimeAccuracy != 2 )
   {
     string error_message ="   - 1\n   - 2\n   ";
     MAC_Error::object()->raise_bad_data_value( exp,
	"AdvectionTimeAccuracy", error_message );
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

   // Build the matrix system
   MAC_ModuleExplorer* se = exp->create_subexplorer( 0,"DDS_NavierStokesSystem" ) ;
   GLOBAL_EQ = DDS_NavierStokesSystem::create( this, se, UF, PF ) ;
   se->destroy() ;

   // Timing routines
   if ( my_rank == is_master )
   {
     SCT_insert_app("Matrix_Assembly&Initialization");
     SCT_insert_app("Pressure predictor");
//     SCT_insert_app("Explicit Velocity step");
//     SCT_insert_app("Velocity x update");
//     SCT_insert_app("Velocity y update");
//     SCT_insert_app("Velocity z update");
     SCT_insert_app("Velocity update");
     SCT_insert_app("Penalty Step");
     SCT_insert_app("Pressure Update");
     SCT_insert_app("Writing CSV");
     SCT_get_elapsed_time("Objects_Creation");
   }
}




//---------------------------------------------------------------------------
DDS_NavierStokes:: ~DDS_NavierStokes( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: ~DDS_NavierStokes" ) ;

   free_DDS_subcommunicators() ;

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

   allocate_mpi_variables(PF,0);
   allocate_mpi_variables(UF,1);

   // Generate solid particles if required
   if (is_solids) {
      Solids_generation(0);
      Solids_generation(1);
      node_property_calculation(PF,0);
      node_property_calculation(UF,1);
//      nodes_temperature_initialization(0);
//      nodes_temperature_initialization(1);
//      nodes_temperature_initialization(2);
//      if (dim == 3) nodes_temperature_initialization(3);
   }
   

   // Initialize velocity vector at the matrix level
   GLOBAL_EQ->initialize_DS_velocity();
   GLOBAL_EQ->initialize_DS_pressure();

   // Direction splitting
   // Assemble 1D tridiagonal matrices
   assemble_1D_matrices(t_it);

   if ( my_rank == is_master ) SCT_get_elapsed_time( "Matrix_Assembly&Initialization" );

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
   // SCT_get_elapsed_time( "Writing CSV" );

   output_L2norm_velocity(0);
   output_L2norm_pressure(0);
   //error_with_analytical_solution();

   if ( my_rank == is_master )
   {
     double cputime = CT_get_elapsed_time();
     cout << endl << "Full problem" << endl;
     write_elapsed_time_smhd(cout,cputime,"Computation time");
     SCT_get_summary(cout,cputime);
   }

   deallocate_mpi_variables(0);
   deallocate_mpi_variables(1);
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

   double vel_divergence = get_velocity_divergence();

   double cfl = UF->compute_CFL( t_it, 0 );
   if ( my_rank == is_master )
      MAC::out() << "CFL: "<< cfl <<endl;

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

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: error_with_analytical_solution ( )
//---------------------------------------------------------------------------
{
   MAC_LABEL(
   "DDS_NavierStokes:: error_with_analytical_solution" ) ;

   // Parameters
   double x, y, z;
   size_t cpp; 
   double bodyterm=0., height;

   // Periodic pressure gradient
   if ( UF->primary_grid()->is_periodic_flow() ) {
      cpp = UF->primary_grid()->get_periodic_flow_direction() ;
      bodyterm = UF->primary_grid()->get_periodic_pressure_drop() / 
               ( UF->primary_grid()->get_main_domain_max_coordinate( cpp )
               - UF->primary_grid()->get_main_domain_min_coordinate( cpp ) ) ;
   }

   height = ( UF->primary_grid()->get_main_domain_max_coordinate(1) - UF->primary_grid()->get_main_domain_min_coordinate(1) ) /2. ;

   for (size_t comp=0;comp<nb_comps[1];++comp) {
      // Get nb of local dof
      size_t_vector local_dof_number( dim, 0 );
      for (size_t l=0;l<dim;++l)
         local_dof_number(l) = UF->get_local_nb_dof( comp, l ) ;

      // Compute error
      double computed_field = 0., analytical_solution = 0.;
      double error_L2 = 0.;
      for (size_t i=0;i<local_dof_number(0);++i) {
         x = UF->get_DOF_coordinate( i, comp, 0 ) ;
         for (size_t j=0;j<local_dof_number(1);++j) {
            y = UF->get_DOF_coordinate( j, comp, 1 ) ;

            if ( dim == 2 ) {
               size_t k = 0 ;
               computed_field = UF->DOF_value( i, j, k, comp, 0 ) ;

               if (comp == cpp) {
                  analytical_solution = (bodyterm/(2.*mu))*((y-height)*(y-height)-height*height);
               }

               if ( UF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) )
	          error_L2 += MAC::sqr( computed_field - analytical_solution ) * UF->get_cell_measure( i, j, k, comp ) ;
	    } else {
               for (size_t k=0;k<local_dof_number(2);++k) {
                  z = UF->get_DOF_coordinate( k, comp, 2 ) ;
                  computed_field = UF->DOF_value( i, j, k, comp, 0 ) ;

                  if (comp == cpp) {
                     analytical_solution = (bodyterm/(2.*mu))*(y*y-height*height);
                  }

                  if ( UF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) )
                     error_L2 += MAC::sqr( computed_field - analytical_solution ) * UF->get_cell_measure( i, j, k, comp ) ;
	       }
	    }
         }
      }

      error_L2 = pelCOMM->sum( error_L2 );
      error_L2 = MAC::sqrt(error_L2);

      if ( my_rank == 0 )
         cout << "L2 Error with analytical solution field " << UF->name() << ", component " << comp << " = " << std::fixed << std::setprecision(16) << error_L2 << endl;
   }

}

//---------------------------------------------------------------------------
size_t
DDS_NavierStokes:: return_row_index (
  FV_DiscreteField const* FF,
  size_t const& comp,
  size_t const& dir,
  size_t const& j,
  size_t const& k )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: return_row_index" ) ;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   size_t p=0;

   if (dim == 2) {
      if (dir == 0) {
         p = j - min_unknown_index(1);
      } else if (dir == 1) {
         p = j - min_unknown_index(0);
      }
   } else if (dim == 3) {
      if (dir == 0) {
         p = (j-min_unknown_index(1))+(1+max_unknown_index(1)-min_unknown_index(1))*(k-min_unknown_index(2));
      } else if (dir == 1) {
         p = (j-min_unknown_index(0))+(1+max_unknown_index(0)-min_unknown_index(0))*(k-min_unknown_index(2));
      } else if (dir == 2) {
         p = (j-min_unknown_index(0))+(1+max_unknown_index(0)-min_unknown_index(0))*(k-min_unknown_index(1));
      }
   }

   return(p);
}

//---------------------------------------------------------------------------
size_t
DDS_NavierStokes:: return_node_index (
  FV_DiscreteField const* FF,
  size_t const& comp,
  size_t const& i,
  size_t const& j,
  size_t const& k )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: return_node_index" ) ;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   size_t_vector i_length(dim,0);
   for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) - 1;
        max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) + 1;
        i_length(l) = 1 + max_unknown_index(l) - min_unknown_index(l);
   }

   size_t local_min_k = 0;
   if (dim == 3) local_min_k = min_unknown_index(2);

   size_t p = (j-min_unknown_index(1)) + i_length(1)*(i-min_unknown_index(0)) + i_length(0)*i_length(1)*(k-local_min_k);

   return(p);

}

//---------------------------------------------------------------------------
double
DDS_NavierStokes:: level_set_function (FV_DiscreteField const* FF, size_t const& m, size_t const& comp, double const& xC, double const& yC, double const& zC, size_t const& type, size_t const& field)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: level_set_solids" ) ;

  PartInput solid = GLOBAL_EQ->get_solid(field);

  double xp = solid.coord[comp]->item(m,0);
  double yp = solid.coord[comp]->item(m,1);
  double zp = solid.coord[comp]->item(m,2);
  double Rp = solid.size[comp]->item(m);

  doubleVector delta(3,0);

  delta(0) = xC-xp;
  delta(1) = yC-yp;
  delta(2) = 0;
  if (dim == 3) delta(2) = zC-zp;

  // Displacement correction in case of periodic boundary condition in any or all directions
  for (size_t dir=0;dir<dim;dir++) {
     if (is_periodic[field][dir]) {
        double isize = FF->primary_grid()->get_main_domain_max_coordinate(dir) - FF->primary_grid()->get_main_domain_min_coordinate(dir);
        delta(dir) = delta(dir) - round(delta(dir)/isize)*isize;
     }
  }

  double level_set = 0.;
  // Type 0 is for circular/spherical solids in 2D/3D system
  if (type == 0) {
     level_set = pow(pow(delta(0),2.)+pow(delta(1),2.)+pow(delta(2),2.),0.5)-Rp;
  }

  return(level_set);

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: Solids_generation (size_t const& field)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: Solids_generation" ) ;

  // Structure of particle input data
  PartInput solid = GLOBAL_EQ->get_solid(field);

  double xp,yp,zp,Rp,vx,vy,vz,Tp,off;

  for (size_t comp=0;comp<nb_comps[field];comp++) {
     ifstream inFile;
     std::ostringstream os2;
     os2 << "./InputFiles/" << solid_filename;
     std::string filename = os2.str();

     inFile.open(filename.c_str());
     string line;
     getline(inFile,line);
     for (size_t i=0;i<Npart;i++) {
        inFile >> xp >> yp >> zp >> Rp >> vx >> vy >> vz >> Tp >> off;
        solid.coord[comp]->set_item(i,0,xp);
        solid.coord[comp]->set_item(i,1,yp);
        solid.coord[comp]->set_item(i,2,zp);
        solid.size[comp]->set_item(i,Rp);
        solid.vel[comp]->set_item(i,0,vx);
        solid.vel[comp]->set_item(i,1,vy);
        solid.vel[comp]->set_item(i,2,vz);
        solid.temp[comp]->set_item(i,Tp);
        solid.inside[comp]->set_item(i,off);
     }
     inFile.close();
  }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: node_property_calculation (FV_DiscreteField const* FF, size_t const& field )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: node_property_calculation" ) ;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  PartInput solid = GLOBAL_EQ->get_solid(field);
  NodeProp node = GLOBAL_EQ->get_node_property(field);

  for (size_t comp=0;comp<nb_comps[field];comp++) {
     // Get local min and max indices; Calculation on the rows next to the proc as well
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) - 1 ;
        max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) + 1 ;
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_unknown_index(2);
        local_max_k = max_unknown_index(2);
     }

     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
        double xC = FF->get_DOF_coordinate( i, comp, 0 ) ;
        double dx = FF->get_cell_size(i,comp,0) ;
        for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
           double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;
           double dy = FF->get_cell_size(j,comp,1) ;
           for (size_t k=local_min_k;k<=local_max_k;++k) {
              double zC = 0.;
              double dC = min(dx,dy);
              if (dim == 3) {
                 zC = FF->get_DOF_coordinate( k, comp, 2 ) ;
                 double dz = FF->get_cell_size(k,comp,2) ;
                 dC = min(dC,dz);
              }
              size_t p = return_node_index(FF,comp,i,j,k);
              for (size_t m=0;m<Npart;m++) {
                 double level_set = level_set_function(FF,m,comp,xC,yC,zC,0,field);
                 level_set *= solid.inside[comp]->item(m);

                 // level_set is xb, if local critical time scale is 0.01 of the global time scale 
                 // then the node is considered inside the solid object
                 // (xb/dC)^2 = 0.01 --> (xb/xC) = 0.1
                 if (level_set <= pow(loc_thres,0.5)*dC) {
                 //if (level_set <= 1.E-1*dC) {
                    node.void_frac[comp]->set_item(p,1.);
                    node.parID[comp]->set_item(p,m);
                    break;
                 }
              }
           }
        }
     }
     // Level 0 is for the intersection matrix corresponding to fluid side
     if (field == 0) {
        assemble_intersection_matrix(PF,comp,0,field);
     } else if (field == 1) {
        assemble_intersection_matrix(UF,comp,0,field);
     }
     // Level 1 is for the intersection matrix corresponding to solids side
//     assemble_intersection_matrix(comp,1);

//     BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(0);
//     b_intersect[0].value[comp]->print_items(MAC::out(),0);
  }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_intersection_matrix ( FV_DiscreteField const* FF, size_t const& comp, size_t const& level, size_t const& field)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: assemble_intersection_matrix" ) ;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  size_t_vector ipos(3,0);
  size_t_array2D local_unknown_extents(dim,2,0);
  size_t_array2D node_neigh(dim,2,0);

  PartInput solid = GLOBAL_EQ->get_solid(field);
  NodeProp node = GLOBAL_EQ->get_node_property(field);
  BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(field,level);

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) - 1 ;
     max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) + 1 ;
     local_unknown_extents(l,0) = 0;
     local_unknown_extents(l,1) = (max_unknown_index(l)-min_unknown_index(l));
  }

  size_t local_min_k = 0;
  size_t local_max_k = 0;

  if (dim == 3) {
     local_min_k = min_unknown_index(2);
     local_max_k = max_unknown_index(2);
  }

  for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
     ipos(0) = i - min_unknown_index(0);
     for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
        ipos(1) = j - min_unknown_index(1);
        for (size_t k=local_min_k;k<=local_max_k;++k) {
           ipos(2) = k - local_min_k;
           size_t p = return_node_index(FF,comp,i,j,k);

           double center_void_frac = 0.;
           if (level == 0) {
              center_void_frac = 1.;
           } else if (level == 1) {
              center_void_frac = 0.;
           }

           if (node.void_frac[comp]->item(p) != center_void_frac) {
              node_neigh(0,0) = return_node_index(FF,comp,i-1,j,k);
              node_neigh(0,1) = return_node_index(FF,comp,i+1,j,k);
              node_neigh(1,0) = return_node_index(FF,comp,i,j-1,k);
              node_neigh(1,1) = return_node_index(FF,comp,i,j+1,k);
              if (dim == 3) {
                 node_neigh(2,0) = return_node_index(FF,comp,i,j,k-1);
                 node_neigh(2,1) = return_node_index(FF,comp,i,j,k+1);
              }

              for (size_t dir=0;dir<dim;dir++) {
                 size_t ii,jj,kk;
                 if (dir == 0) {
                    ii = i;jj = j; kk = k;
                 } else if (dir == 1) {
                    ii = j;jj = i; kk = k;
                 } else if (dir == 2) {
                    ii = k;jj = i; kk = j;
                 }
                 for (size_t off=0;off<2;off++) {
                    size_t left, right;
                    if (off == 0) {
                       left = ii-1; right = ii;
                    } else if (off == 1) {
                       left = ii; right = ii+1;
                    }

                    if ((node.void_frac[comp]->item(node_neigh(dir,off)) != node.void_frac[comp]->item(p)) && (ipos(dir) != local_unknown_extents(dir,off))) {
                       double xb = find_intersection(FF,left,right,jj,kk,comp,dir,off,field);
                       b_intersect[dir].offset[comp]->set_item(p,off,1);
                       b_intersect[dir].value[comp]->set_item(p,off,xb);
                       size_t par_id = node.parID[comp]->item(node_neigh(dir,off));
                       double field_value = solid.temp[comp]->item(par_id);
                       b_intersect[dir].field_var[comp]->set_item(p,off,field_value);
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
DDS_NavierStokes:: find_intersection ( FV_DiscreteField const* FF, size_t const& left, size_t const& right, size_t const& yconst, size_t const& zconst, size_t const& comp, size_t const& dir, size_t const& off, size_t const& field)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: find_intersection" ) ;

  PartInput solid = GLOBAL_EQ->get_solid(field);
  NodeProp node = GLOBAL_EQ->get_node_property(field);

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  size_t_vector side(2,0);

  side(0) = left;
  side(1) = right;

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) - 1;
     max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) + 1;
  }

  double funl=0., func=0., funr=0.;

  double xleft = FF->get_DOF_coordinate( side(0), comp, dir ) ;
  double xright = FF->get_DOF_coordinate( side(1), comp, dir ) ;

  double yvalue=0.,zvalue=0.;
  size_t p;

  if (dir == 0) {
     yvalue = FF->get_DOF_coordinate( yconst, comp, 1 ) ;
     if (dim == 3) zvalue = FF->get_DOF_coordinate( zconst, comp, 2 ) ;
     p = return_node_index(FF,comp,side(off),yconst,zconst);
  } else if (dir == 1) {
     yvalue = FF->get_DOF_coordinate( yconst, comp, 0 ) ;
     if (dim == 3) zvalue = FF->get_DOF_coordinate( zconst, comp, 2 ) ;
     p = return_node_index(FF,comp,yconst,side(off),zconst);
  } else if (dir == 2) {
     yvalue = FF->get_DOF_coordinate( yconst, comp, 0 ) ;
     if (dim == 3) zvalue = FF->get_DOF_coordinate( zconst, comp, 1 ) ;
     p = return_node_index(FF,comp,yconst,zconst,side(off));
  }

  size_t id = node.parID[comp]->item(p);

  double xcenter;

  if (dir == 0) {
     funl = level_set_function(FF,id,comp,xleft,yvalue,zvalue,0,field);
     funr = level_set_function(FF,id,comp,xright,yvalue,zvalue,0,field);
  } else if (dir == 1) {
     funl = level_set_function(FF,id,comp,yvalue,xleft,zvalue,0,field);
     funr = level_set_function(FF,id,comp,yvalue,xright,zvalue,0,field);
  } else if (dir == 2) {
     funl = level_set_function(FF,id,comp,yvalue,zvalue,xleft,0,field);
     funr = level_set_function(FF,id,comp,yvalue,zvalue,xright,0,field);
  }

  // In case both the points are on the same side of solid interface
  // This will occur when the point just outside the solid interface will be considered inside the solid
  // This condition enables the intersection with the interface using the point in fluid and the ACTUAL node in the solid 
  if (funl*funr > 0.) {
     double dx = FF->get_cell_size(side(off),comp,dir) ;
     if (off == 0) {
        xleft = xleft - dx;
     } else {
        xright = xright + dx;
     }
  }

  if (dir == 0) {
     funl = level_set_function(FF,id,comp,xleft,yvalue,zvalue,0,field);
     funr = level_set_function(FF,id,comp,xright,yvalue,zvalue,0,field);
  } else if (dir == 1) {
     funl = level_set_function(FF,id,comp,yvalue,xleft,zvalue,0,field);
     funr = level_set_function(FF,id,comp,yvalue,xright,zvalue,0,field);
  } else if (dir == 2) {
     funl = level_set_function(FF,id,comp,yvalue,zvalue,xleft,0,field);
     funr = level_set_function(FF,id,comp,yvalue,zvalue,xright,0,field);
  }

  // If the shifted point is also physically outside the solid then xb = dx
  if (funl*funr > 0.) {
     xcenter = FF->get_DOF_coordinate( side(off), comp, dir ) ;
  } else {
     // Bisection method algorithm
     while (MAC::abs(xright-xleft) > 1.E-15) {
        xcenter = (xleft+xright)/2.;
        if (dir == 0) {
           funl = level_set_function(FF,id,comp,xleft,yvalue,zvalue,0,field);
           func = level_set_function(FF,id,comp,xcenter,yvalue,zvalue,0,field);
        } else if (dir == 1) {
           funl = level_set_function(FF,id,comp,yvalue,xleft,zvalue,0,field);
           func = level_set_function(FF,id,comp,yvalue,xcenter,zvalue,0,field);
        } else if (dir == 2) {
           funl = level_set_function(FF,id,comp,yvalue,zvalue,xleft,0,field);
           func = level_set_function(FF,id,comp,yvalue,zvalue,xcenter,0,field);
        }

        if ((func == 1.E-16) || ((xcenter-xleft)/2. <= 1.E-16)) break;

        if (func*funl >= 1.E-16) {
           xleft = xcenter;
        } else {
           xright = xcenter;
        }
     }
  }

  if (off == 0) {
     xcenter = MAC::abs(xcenter - FF->get_DOF_coordinate( side(1), comp, dir ));
  } else if (off == 1) {
     xcenter = MAC::abs(xcenter - FF->get_DOF_coordinate( side(0), comp, dir ));
  }

  return (xcenter);
}


//---------------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_field_matrix (
  FV_DiscreteField const* FF,
  FV_TimeIterator const* t_it,
  double const& gamma,
  size_t const& comp,
  size_t const& dir, 
  size_t const& field,
  size_t const& j,
  size_t const& k,
  size_t const& r_index )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: assemble_field_matrix" ) ;

   // Parameters
   double dxr,dxl,dx,xR,xL,xC,right=0.,left=0.,center=0.;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   // Perform assembling
   size_t m, i;
   TDMatrix* A = GLOBAL_EQ-> get_A(field);

   double Aee_diagcoef=0.;

   for (m=0,i=min_unknown_index(dir);i<=max_unknown_index(dir);++i,++m) {
      xC= FF->get_DOF_coordinate( i, comp, dir) ;
      xR= FF->get_DOF_coordinate( i+1,comp, dir) ;
      xL= FF->get_DOF_coordinate( i-1, comp, dir) ;

      dx = FF->get_cell_size( i,comp, dir);

      dxr= xR - xC;
      dxl= xC - xL;

      size_t k_min, k_max;
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
         unsteady_term = rho*(FF->get_cell_size(i,comp,dir))/(t_it->time_step());
      }

      center = - (right+left);

      if (dim == 2) {
         k_min = 0; k_max = 0;
      } else {
         k_min = min_unknown_index(2); k_max = max_unknown_index(2);
      }

      bool r_bound = false;
      bool l_bound = false;
      // All the proc will have open right bound, except last proc for non periodic systems
      if ((is_periodic[field][dir] != 1) && (rank_in_i[dir] == nb_ranks_comm_i[dir]-1)) r_bound = true;
      // All the proc will have open left bound, except first proc for non periodic systems
      if ((is_periodic[field][dir] != 1) && (rank_in_i[dir] == 0)) l_bound = true;

      // Since, this function is used in all directions; 
      // ii, jj, and kk are used to convert the passed arguments corresponding to correct direction
      size_t ii=0,jj=0,kk=0;

      // Condition for handling the pressure neumann conditions at wall
      if (i==min_unknown_index(dir) && l_bound) {
         if (dir == 0) {
            ii = i-1; jj = min_unknown_index(1); kk = k_min;
         } else if (dir == 1) {
            ii = min_unknown_index(0); jj = i-1; kk = k_min;
         } else if (dir == 2) {
            ii = min_unknown_index(0); jj = min_unknown_index(1); kk = i-1;
         }
         if (FF->DOF_in_domain(ii,jj,kk,comp) && FF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
            value = center;
         } else {
            value = -right;
         }
      } else if (i==max_unknown_index(dir) && r_bound) {
         if (dir == 0) {
            ii = i+1; jj = max_unknown_index(1); kk = k_max;
         } else if (dir == 1) {
            ii = max_unknown_index(0); jj = i+1; kk = k_max;
         } else if (dir == 2) {
            ii = max_unknown_index(0); jj = max_unknown_index(1); kk = i+1;
         }
         if (FF->DOF_in_domain(ii,jj,kk,comp) && FF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
            value = center;
         } else {
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
         // For last index, Aee comes from this proc as it is interface unknown wrt this proc
         A[dir].ie[comp][r_index]->set_item(m-1,rank_in_i[dir],left);
         Aee_diagcoef = value;
         A[dir].ei[comp][r_index]->set_item(rank_in_i[dir],m-1,left);
      }

      // Set Aii_sub_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_periodic[field][dir] != 1)) {
         if (i > min_unknown_index(dir)) A[dir].ii_sub[comp][r_index]->set_item(m-1,left);
      } else {
         if (i<max_unknown_index(dir)) {
            if (i>min_unknown_index(dir)) {
               A[dir].ii_sub[comp][r_index]->set_item(m-1,left);
            }
         }
      }

      // Set Aii_super_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_periodic[field][dir] != 1)) {
         if (i < max_unknown_index(dir)) A[dir].ii_super[comp][r_index]->set_item(m,right);
      } else {
         if (i < max_unknown_index(dir)-1) {
            A[dir].ii_super[comp][r_index]->set_item(m,right);
         }
      }

      // Set Aii_main_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_periodic[field][dir] != 1)) {
         A[dir].ii_main[comp][r_index]->set_item(m,value);
      } else {
         if (i<max_unknown_index(dir)) {
            A[dir].ii_main[comp][r_index]->set_item(m,value);
         }
      }
   } // End of for loop

   GLOBAL_EQ->pre_thomas_treatment(comp,dir,A,r_index);

   return (Aee_diagcoef);
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_field_schur_matrix (struct TDMatrix *A, size_t const& comp, size_t const& dir, double const& Aee_diagcoef, size_t const& field, size_t const& r_index )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: assemble_field_schur_matrix" ) ;
   // Compute the product matrix for each proc

   if (nb_ranks_comm_i[dir]>1) {

      ProdMatrix* Ap = GLOBAL_EQ->get_Ap(field);

      GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir,field,r_index);

      LA_SeqMatrix* product_matrix = Ap[dir].ei_ii_ie[comp];
      LA_SeqMatrix* receive_matrix = product_matrix->create_copy(this,product_matrix);

      if ( rank_in_i[dir] == 0 ) {
         A[dir].ee[comp][r_index]->set_item(0,0,Aee_diagcoef);
   	 for (size_t i=1;i<nb_ranks_comm_i[dir];++i) {

            // Create the container to receive
            size_t nbrows = product_matrix->nb_rows();
            size_t nb_received_data = pow(nbrows,2)+1;
            double * received_data = new double [nb_received_data];

            // Receive the data
            static MPI_Status status ;
            MPI_Recv( received_data, nb_received_data, MPI_DOUBLE, i, 0,
                            DDS_Comm_i[dir], &status ) ;

            // Transfer the received data to the receive matrix
            for (int k=0;k<nbrows;k++) {
               for (int j=0;j<nbrows;j++) {
                  // Assemble the global product matrix by adding contributions from all the procs
                  receive_matrix->add_to_item(k,j,received_data[k*(nbrows)+j]);
               }
            }

  	    if (is_periodic[field][dir] == 0) {
               if (i<nb_ranks_comm_i[dir]-1) {
                  // Assemble the global Aee matrix 
                  // No periodic condition in x. So no fe contribution from last proc
                  A[dir].ee[comp][r_index]->set_item(i,i,received_data[nb_received_data-1]);
               }   
            } else {
               // Assemble the global Aee matrix
               // Periodic condition in x. So there is fe contribution from last proc
               A[dir].ee[comp][r_index]->set_item(i,i,received_data[nb_received_data-1]); 
            }
            delete [] received_data;
         }
      } else {
         // Create the packed data container
         size_t nbrows = product_matrix->nb_rows();
         size_t nb_send_data = pow(nbrows,2)+1;
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

         // Fill the last element of packed data with the diagonal coefficient Aee
         packed_data[nb_send_data-1] = Aee_diagcoef;

         // Send the data
         MPI_Send( packed_data, nb_send_data, MPI_DOUBLE, 0, 0, DDS_Comm_i[dir] ) ;

         delete [] packed_data;
      }

      // Assemble the schlur complement in the master proc

      if (rank_in_i[dir] == 0){
         TDMatrix* Schur = GLOBAL_EQ-> get_Schur(field);
         size_t nb_row = Schur[dir].ii_main[comp][r_index]->nb_rows();
         for (int p = 0; p < nb_row; p++) {
            Schur[dir].ii_main[comp][r_index]->set_item(p,A[dir].ee[comp][r_index]->item(p,p)-receive_matrix->item(p,p));
            if (p < nb_row-1) Schur[dir].ii_super[comp][r_index]->set_item(p,-receive_matrix->item(p,p+1));
            if (p > 0) Schur[dir].ii_sub[comp][r_index]->set_item(p-1,-receive_matrix->item(p,p-1));
            // In case of periodic and multi-processor, there will be a variant of Tridiagonal matrix instead of normal format
            if (is_periodic[field][dir] == 1) {
               Schur[dir].ie[comp][r_index]->set_item(p,0,-receive_matrix->item(p,nb_row));
               Schur[dir].ei[comp][r_index]->set_item(0,p,-receive_matrix->item(nb_row,p));
            }
         }
         // Pre-thomas treatment on Schur complement
         GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur,r_index);

         // In case of periodic and multi-processor, there will be a variant of Tridiagonal matrix instead of normal format
         // So, Schur complement of Schur complement is calculated
         if (is_periodic[field][dir] == 1) {
            Schur[dir].ee[comp][r_index]->set_item(0,0,A[dir].ee[comp][r_index]->item(nb_row,nb_row)-receive_matrix->item(nb_row,nb_row));

            ProdMatrix* SchurP = GLOBAL_EQ->get_SchurP(field);
            GLOBAL_EQ->compute_product_matrix_interior(Schur,SchurP,comp,0,dir,r_index);

            TDMatrix* DoubleSchur = GLOBAL_EQ-> get_DoubleSchur(field);
            size_t nb_row = DoubleSchur[dir].ii_main[comp][r_index]->nb_rows();
            DoubleSchur[dir].ii_main[comp][r_index]->set_item(0,Schur[dir].ee[comp][r_index]->item(0,0)-SchurP[dir].ei_ii_ie[comp]->item(0,0));
         }
      }
   } else if (is_periodic[field][dir] == 1) {
      // Condition for single processor in any direction with periodic boundary conditions
      ProdMatrix* Ap = GLOBAL_EQ->get_Ap(field);
      GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir, field,r_index);

      LA_SeqMatrix* product_matrix = Ap[dir].ei_ii_ie[comp];
      LA_SeqMatrix* receive_matrix = product_matrix->create_copy(this,product_matrix);

      A[dir].ee[comp][r_index]->set_item(0,0,Aee_diagcoef);

      TDMatrix* Schur = GLOBAL_EQ-> get_Schur(field);
      size_t nb_row = Schur[dir].ii_main[comp][r_index]->nb_rows();
      for (int p = 0; p < nb_row; p++) {
         Schur[dir].ii_main[comp][r_index]->set_item(p,A[dir].ee[comp][r_index]->item(p,p)-receive_matrix->item(p,p));
         if (p < nb_row-1) Schur[dir].ii_super[comp][r_index]->set_item(p,-receive_matrix->item(p,p+1));
         if (p > 0) Schur[dir].ii_sub[comp][r_index]->set_item(p-1,-receive_matrix->item(p,p-1));
      }
      GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur,r_index); 
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_1D_matrices ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: assemble_1D_matrices" ) ;

   double gamma = mu/2.0;

   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   // Assemble the matrices for pressure field(0) and velocity(1) field
   for (size_t field=0;field<2;field++) {
      for (size_t comp=0;comp<nb_comps[field];comp++) {
         // Get local min and max indices
         for (size_t l=0;l<dim;++l) {
            if (field == 0) {
               min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( comp, l ) ;
               max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( comp, l ) ;
            } else if (field == 1) {
               min_unknown_index(l) = UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
               max_unknown_index(l) = UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
            }
         }
         for (size_t dir=0;dir<dim;dir++) {
            size_t dir_j, dir_k;
            size_t local_min_k = 0;
            size_t local_max_k = 0;

            if (dir == 0) {
               dir_j = 1; dir_k = 2;
            } else if (dir == 1) {
               dir_j = 0; dir_k = 2;
            } else if (dir == 2) {
               dir_j = 0; dir_k = 1;
            }

            if (dim == 3) {
               local_min_k = min_unknown_index(dir_k);
               local_max_k = max_unknown_index(dir_k);
            }

            for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j) {
               for (size_t k=local_min_k; k <= local_max_k; ++k) {
                  size_t r_index;
                  double Aee_diagcoef;
                  if (field == 0) {
                     r_index = return_row_index (PF,comp,dir,j,k);
                     Aee_diagcoef = assemble_field_matrix (PF,t_it,gamma,comp,dir,0,j,k,r_index);
                  } else if (field == 1) {
                     r_index = return_row_index (UF,comp,dir,j,k);
                     Aee_diagcoef = assemble_field_matrix (UF,t_it,gamma,comp,dir,1,j,k,r_index);
                  }
                  TDMatrix* A = GLOBAL_EQ-> get_A(field);
                  assemble_field_schur_matrix (A,comp,dir,Aee_diagcoef,field,r_index);
               }
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

   size_t i, j, k;

  double value=0.;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  // First Equation

  // Get local min and max indices
  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_on_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_on_proc( 0, l ) ;
  }

  for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
     for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
        if (dim == 2) {
           k=0;
           // Set P*_n as sum of P_(n-1/2)+phi_(n-1/2)
           value = PF->DOF_value( i, j, k, 0, 0 ) + PF->DOF_value( i, j, k, 0, 1 );
           PF->set_DOF_value( i, j, k, 0, 1, value);
        } else {
           for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
              // Set P*_n as sum of P_(n-1/2)+phi_(n-1/2)
              value = PF->DOF_value( i, j, k, 0, 0 ) + PF->DOF_value( i, j, k, 0, 1 );
              PF->set_DOF_value( i, j, k, 0, 1, value);
           }
        }
     }
  }
  
  PF->set_neumann_DOF_values();
}

//---------------------------------------------------------------------------
double
DDS_NavierStokes:: compute_un_component ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k, size_t const& dir, size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: compute_un_component" ) ;

   double xhr,xhl,xright,xleft,yhr,yhl,yright,yleft;
   double zhr,zhl,zright,zleft, value=0.;

   if (dir == 0) {
      xhr= UF->get_DOF_coordinate( i+1,comp, 0 ) - UF->get_DOF_coordinate( i, comp, 0 ) ;
      xhl= UF->get_DOF_coordinate( i, comp, 0 ) - UF->get_DOF_coordinate( i-1, comp, 0 ) ;
      xright = UF->DOF_value( i+1, j, k, comp, level ) - UF->DOF_value( i, j, k, comp, level ) ;
      xleft = UF->DOF_value( i, j, k, comp, level ) - UF->DOF_value( i-1, j, k, comp, level ) ;

      //xvalue = xright/xhr - xleft/xhl;
      if (UF->DOF_in_domain( i-1, j, k, comp) && UF->DOF_in_domain( i+1, j, k, comp))
         value = xright/xhr - xleft/xhl;
      else if (UF->DOF_in_domain( i-1, j, k, comp))
         value = - xleft/xhl;
      else
         value = xright/xhr;
   } else if (dir == 1) {
      yhr= UF->get_DOF_coordinate( j+1,comp, 1 ) - UF->get_DOF_coordinate( j, comp, 1 ) ;
      yhl= UF->get_DOF_coordinate( j, comp, 1 ) - UF->get_DOF_coordinate( j-1, comp, 1 ) ;
      yright = UF->DOF_value( i, j+1, k, comp, level ) - UF->DOF_value( i, j, k, comp, level ) ;
      yleft = UF->DOF_value( i, j, k, comp, level ) - UF->DOF_value( i, j-1, k, comp, level ) ;

      //yvalue = yright/yhr - yleft/yhl;
      if (UF->DOF_in_domain(i, j-1, k, comp) && UF->DOF_in_domain(i, j+1, k, comp))
         value = yright/yhr - yleft/yhl;
      else if(UF->DOF_in_domain(i, j-1, k, comp))
         value = - yleft/yhl;
      else
         value = yright/yhr;
   } else if (dir == 2) {
      zhr= UF->get_DOF_coordinate( k+1,comp, 2 ) - UF->get_DOF_coordinate( k, comp, 2 ) ;
      zhl= UF->get_DOF_coordinate( k, comp, 2 ) - UF->get_DOF_coordinate( k-1, comp, 2 ) ;
      zright = UF->DOF_value( i, j, k+1, comp, level ) - UF->DOF_value( i, j, k, comp, level ) ;
      zleft = UF->DOF_value( i, j, k, comp, level ) - UF->DOF_value( i, j, k-1, comp, level ) ;

      //zvalue = zright/zhr - zleft/zhl;
      if (UF->DOF_in_domain(i, j, k-1, comp) && UF->DOF_in_domain(i, j, k+1, comp))
         value = zright/zhr - zleft/zhl;
      else if(UF->DOF_in_domain(i, j, k-1, comp))
         value = - zleft/zhl;
      else
         value = zright/zhr;
   }

   return(value);

}


//---------------------------------------------------------------------------
double
DDS_NavierStokes:: velocity_local_rhs ( size_t const& j, size_t const& k, double const& gamma, FV_TimeIterator const* t_it, size_t const& comp, size_t const& dir)
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

   size_t i,pos;
   int m;

   // Compute VEC_rhs_x = rhs in x
   double dC,hr=0,hl=0;
   double fe=0.;

   // Vector for fi
   LocalVector* VEC = GLOBAL_EQ->get_VEC(1);

   for (i=min_unknown_index(dir);i<=max_unknown_index(dir);++i) {
     double value=0.;
     pos = i - min_unknown_index(dir);

     // Get contribution of un
     hl= UF->get_DOF_coordinate(i,comp,dir) - UF->get_DOF_coordinate(i-1,comp,dir) ;
     hr= UF->get_DOF_coordinate(i+1,comp,dir) - UF->get_DOF_coordinate(i,comp,dir) ;

     dC = UF->get_cell_size(i,comp,dir);

     // x direction
     if (dir == 0) {
        value = compute_un_component(comp,i,j,k,dir,3);
     // y direction
     } else if (dir == 1) {
        if (dim == 2) {
           value = compute_un_component(comp,j,i,k,dir,1);
        } else if (dim == 3) {
           value = compute_un_component(comp,j,i,k,dir,4);
        }

     // z direction
     } else if (dir == 2) {
        value = compute_un_component(comp,j,k,i,dir,1);
     }

     double temp_val=0.;
     if (dir == 0) {
        temp_val = (UF->DOF_value(i,j,k,comp,0)*dC*rho)/(t_it->time_step()) - gamma*value;
     } else if (dir == 1) {
        temp_val = (UF->DOF_value(j,i,k,comp,3)*dC*rho)/(t_it->time_step()) - gamma*value;
     } else if (dir == 2) {
        temp_val = (UF->DOF_value(j,k,i,comp,4)*dC*rho)/(t_it->time_step()) - gamma*value;
     }

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

   // Since, this function is used in all directions; 
   // ii, jj, and kk are used to convert the passed arguments corresponding to correct direction
   size_t ii=0,jj=0,kk=0;

   // Effect of boundary conditions in case of non-periodic direction
   m = int(min_unknown_index(dir)) - 1;

   if (dir == 0) {
      ii = m; jj = j; kk = k;
   } else if (dir == 1) {
      ii = j; jj = m; kk = k;
   } else if (dir == 2) {
      ii = j; jj = k; kk = m;
   }

   if ( UF->DOF_in_domain(ii,jj,kk,comp))
      if ( UF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
         double ai = 1/(UF->get_DOF_coordinate(m+1,comp,dir) - UF->get_DOF_coordinate(m,comp,dir));
         double dirichlet_value = UF->DOF_value(ii,jj,kk,comp,1) ;
         VEC[dir].local_T[comp]->add_to_item( 0, + gamma * ai * dirichlet_value );
      }

   m = int(max_unknown_index(dir)) + 1;

   if (dir == 0) {
      ii = m; jj = j; kk = k;
   } else if (dir == 1) {
      ii = j; jj = m; kk = k;
   } else if (dir == 2) {
      ii = j; jj = k; kk = m;
   }

   if ( UF->DOF_in_domain(ii,jj,kk,comp))
      if ( UF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
      double ai = 1/(UF->get_DOF_coordinate(m,comp,dir) - UF->get_DOF_coordinate(m-1,comp,dir));
      double dirichlet_value = UF->DOF_value(ii,jj,kk,comp,1) ;
      VEC[dir].local_T[comp]->add_to_item( VEC[dir].local_T[comp]->nb_rows()-1 , + gamma * ai * dirichlet_value );
   }

   return fe;
}     

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: unpack_compute_ue_pack(size_t const& comp, size_t const& dir, size_t const& p, size_t const& field)
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
      VEC[dir].T[comp]->set_item(nb_ranks_comm_i[dir]-1,first_pass[field][dir].send[comp][rank_in_i[dir]][3*p]);
   VEC[dir].T[comp]->set_item(0,first_pass[field][dir].send[comp][rank_in_i[dir]][3*p+1]);
   VEC[dir].interface_T[comp]->set_item(0,first_pass[field][dir].send[comp][rank_in_i[dir]][3*p+2]);

   // Vec_temp might contain previous values

   for (size_t i=1;i<nb_ranks_comm_i[dir];i++) {
      if (i!=nb_ranks_comm_i[dir]-1) {
         VEC[dir].T[comp]->add_to_item(i-1,first_pass[field][dir].receive[comp][i][3*p]);
         VEC[dir].T[comp]->add_to_item(i,first_pass[field][dir].receive[comp][i][3*p+1]);
         VEC[dir].interface_T[comp]->set_item(i,first_pass[field][dir].receive[comp][i][3*p+2]);  // Assemble the interface rhs fe
      } else {
         if (is_periodic[field][dir] ==0) {
            VEC[dir].T[comp]->add_to_item(i-1,first_pass[field][dir].receive[comp][i][3*p]);
         } else{
            VEC[dir].T[comp]->add_to_item(i-1,first_pass[field][dir].receive[comp][i][3*p]);
            // If periodic in x, last proc has an interface unknown
            VEC[dir].T[comp]->add_to_item(i,first_pass[field][dir].receive[comp][i][3*p+1]);
            VEC[dir].interface_T[comp]->set_item(i,first_pass[field][dir].receive[comp][i][3*p+2]);
         }
      }
   }

   for (size_t i=0;i<nb_interface_unknowns;i++) {
      VEC[dir].interface_T[comp]->set_item(i,VEC[dir].interface_T[comp]->item(i)-VEC[dir].T[comp]->item(i)); // Get fe - Aei*xi to solve for ue
   }

   // Solve for ue (interface unknowns) in the master proc
   DS_interface_unknown_solver(VEC[dir].interface_T[comp],comp,dir,field,p);

   for (size_t i=1;i<nb_ranks_comm_i[dir];++i) {
      if (i != nb_ranks_comm_i[dir]-1) {
         second_pass[field][dir].send[comp][i][2*p+0]=VEC[dir].interface_T[comp]->item(i-1);
         second_pass[field][dir].send[comp][i][2*p+1]=VEC[dir].interface_T[comp]->item(i);
      } else {
         second_pass[field][dir].send[comp][i][2*p+0]=VEC[dir].interface_T[comp]->item(i-1);
         if (is_periodic[field][dir])
            second_pass[field][dir].send[comp][i][2*p+1]=VEC[dir].interface_T[comp]->item(i);
         else
            second_pass[field][dir].send[comp][i][2*p+1]=0;
      }
   }
}

//----------------------------------------------------------------------
void
DDS_NavierStokes::DS_interface_unknown_solver( LA_SeqVector* interface_rhs, size_t const& comp, size_t const& dir, size_t const& field, size_t const& r_index )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_interface_unknown_solver" ) ;

   TDMatrix* Schur = GLOBAL_EQ->get_Schur(field);

   // Condition for variant of Tridiagonal Schur complement in Perioidic direction with multi-processor 
   if ((is_periodic[field][dir] == 1) && (nb_ranks_comm_i[dir] != 1)) {
      LocalVector* Schur_VEC = GLOBAL_EQ->get_Schur_VEC(field);
      TDMatrix* DoubleSchur = GLOBAL_EQ-> get_DoubleSchur(field);

      // Transfer interface_rhs to Schur VEC (i.e. S_fi and S_fe)
      size_t nrows = Schur_VEC[dir].local_T[comp]->nb_rows();
      for (size_t i = 0; i < nrows; i++) {
          Schur_VEC[dir].local_T[comp]->set_item(i,interface_rhs->item(i));
      }
      Schur_VEC[dir].interface_T[comp]->set_item(0,interface_rhs->item(nrows));

      // Calculate Sei*(Sii)-1*S_fi
      compute_Aei_ui(Schur,Schur_VEC,comp,dir,r_index);

      // Calculate S_fe - Sei*(Sii)-1*S_fi
      Schur_VEC[dir].interface_T[comp]->set_item(0,Schur_VEC[dir].interface_T[comp]->item(0)-Schur_VEC[dir].T[comp]->item(0));

      // Calculate S_ue, using Schur complement of Schur complement
      GLOBAL_EQ->mod_thomas_algorithm(DoubleSchur, Schur_VEC[dir].interface_T[comp], comp, dir,r_index);

      // Calculate S_fi-Sie*S_ue
      Schur[dir].ie[comp][r_index]->multiply_vec_then_add(Schur_VEC[dir].interface_T[comp],Schur_VEC[dir].local_T[comp],-1.0,1.0);

      // Calculate S_ui
      GLOBAL_EQ->mod_thomas_algorithm(Schur, Schur_VEC[dir].local_T[comp], comp, dir,r_index);

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
DDS_NavierStokes:: unpack_ue(size_t const& comp, double * received_data, size_t const& dir, int const& p, size_t const& field)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: unpack_ue" ) ;
 
   LocalVector* VEC = GLOBAL_EQ->get_VEC(field);

   if (rank_in_i[dir] != nb_ranks_comm_i[dir]-1) {
      VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,received_data[2*p]);
      VEC[dir].interface_T[comp]->set_item(rank_in_i[dir],received_data[2*p+1]);
   } else {
      if (is_periodic[field][dir] ==0) {
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,received_data[2*p]);
      } else {
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,received_data[2*p]);
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir],received_data[2*p+1]);
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: solve_interface_unknowns ( FV_DiscreteField* FF, double const& gamma,  FV_TimeIterator const* t_it, size_t const& comp, size_t const& dir, size_t const& field)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: solve_interface_unknowns" ) ;

   size_t i,j,p;
   size_t k =0;

   //first_pass[field][dir_i].send[comp][rank_in_i[dir_i]], first_pass[field][dir_i].size[comp],

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   TDMatrix* A = GLOBAL_EQ-> get_A(field);
   LocalVector* VEC = GLOBAL_EQ->get_VEC(field);

   // Array declaration for sending data from master to all slaves
   size_t local_length_j=0, local_length_k=0;
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

   local_length_j = (local_max_j-local_min_j+1);
   local_length_k = (local_max_k-local_min_k+1);

   // Send and receive the data first pass
   if ( rank_in_i[dir] == 0 ) {
      if (nb_ranks_comm_i[dir] != 1) {      
         for (i=1;i<nb_ranks_comm_i[dir];++i) {
            // Receive the data
            static MPI_Status status;
            MPI_Recv( first_pass[field][dir].receive[comp][i], first_pass[field][dir].size[comp], MPI_DOUBLE, i, 0,
                              DDS_Comm_i[dir], &status ) ;
         }
      }

      // Solve system of interface unknowns for each y
      if (dim == 2) {
	 for (j=local_min_j;j<=local_max_j;j++) {

            p = j-local_min_j;
  
            unpack_compute_ue_pack(comp,dir,p,field);

            // Need to have the original rhs function assembled for corrosponding j,k pair
            double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir,field);

            // Setup RHS = fi - Aie*xe for solving ui
            A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_NavierStokes_solver(FF,j,k,min_unknown_index(dir),comp,dir,field,p);
         }
      } else {
         for (k=local_min_k;k<=local_max_k;k++) {
            for (j=local_min_j;j<=local_max_j;j++) {
	      
	       p = (j-local_min_j)+local_length_j*(k-local_min_k);

               unpack_compute_ue_pack(comp,dir,p,field);

     	       // Need to have the original rhs function assembled for corrosponding j,k pair
               double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir,field);

               // Setup RHS = fi - Aie*xe for solving ui
               A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

               // Solve ui and transfer solution into distributed vector
               GLOBAL_EQ->DS_NavierStokes_solver(FF,j,k,min_unknown_index(dir),comp,dir,field,p);
            }
         }
      }

   } else {
      // Send the packed data to master
      MPI_Send( first_pass[field][dir].send[comp][rank_in_i[dir]], first_pass[field][dir].size[comp], MPI_DOUBLE, 0, 0, DDS_Comm_i[dir] ) ;
   }

   // Send the data from master iff multi processor are used
   if (nb_ranks_comm_i[dir] != 1) {
      if ( rank_in_i[dir] == 0 ) {
         for (i=1;i<nb_ranks_comm_i[dir];++i) {
            MPI_Send( second_pass[field][dir].send[comp][i], second_pass[field][dir].size[comp], MPI_DOUBLE, i, 0, DDS_Comm_i[dir] ) ;
         }
      } else {
         // Receive the data
         static MPI_Status status ;
         MPI_Recv( second_pass[field][dir].send[comp][rank_in_i[dir]], first_pass[field][dir].size[comp], MPI_DOUBLE, 0, 0, DDS_Comm_i[dir], &status ) ;

         // Solve the system of equations in each proc

         if (dim == 2) {
            for (j = local_min_j;j<=local_max_j;j++) {
               p = j-local_min_j;

	       unpack_ue(comp,second_pass[field][dir].send[comp][rank_in_i[dir]],dir,p,field);

               // Need to have the original rhs function assembled for corrosponding j,k pair
               double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir,field);
 
               // Setup RHS = fi - Aie*xe for solving ui
               A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

               // Solve ui and transfer solution into distributed vector
               GLOBAL_EQ->DS_NavierStokes_solver(FF,j,k,min_unknown_index(dir),comp,dir,field,p);
            }
         } else {
            for (k = local_min_k;k<=local_max_k;k++) {
               for (j = local_min_j;j<=local_max_j;j++) {
                  p = (j-local_min_j)+local_length_j*(k-local_min_k);
   
                  unpack_ue(comp,second_pass[field][dir].send[comp][rank_in_i[dir]],dir,p,field);

                  // Need to have the original rhs function assembled for corrosponding j,k pair
                  double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir,field);

                  // Setup RHS = fi - Aie*xe for solving ui
                  A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

                  // Solve ui and transfer solution into distributed vector
                  GLOBAL_EQ->DS_NavierStokes_solver(FF,j,k,min_unknown_index(dir),comp,dir,field,p);
               }
            }
         }
      }
   }
}

//---------------------------------------------------------------------------
double
DDS_NavierStokes:: compute_p_component ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k)
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

   if (dim == 2) {
      if (comp == 0) {
         value = (PF->DOF_value( shift.i+i, j, 0, 0, 1 ) - PF->DOF_value( shift.i+i-1, j, 0, 0, 1 ))*dyC;
      } else {
         value = (PF->DOF_value( i, shift.j+j, 0, 0, 1 ) - PF->DOF_value( i, shift.j+j-1, 0, 0, 1 ))*dxC;
      }
   } else if (dim == 3) {
      double dzC = UF->get_cell_size( k, comp, 2 ) ;
      if (comp == 0) {
         value = (PF->DOF_value( shift.i+i, j, k, 0, 1 ) - PF->DOF_value( shift.i+i-1, j, k, 0, 1 ))*dyC*dzC;
      } else if(comp==1) {
         value = (PF->DOF_value( i, shift.j+j, k, 0, 1 ) - PF->DOF_value( i, shift.j+j-1, k, 0, 1 ))*dxC*dzC;
      } else {
         value = (PF->DOF_value( i, j, shift.k+k, 0, 1 ) - PF->DOF_value( i, j, shift.k+k-1, 0, 1 ))*dxC*dyC;
      }
   }
   return(value);
}

//---------------------------------------------------------------------------
double
DDS_NavierStokes:: compute_adv_component ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: compute_adv_component" ) ;
   double ugradu = 0., value = 0.;

   if ( AdvectionScheme == "TVD" )
      ugradu = assemble_advection_TVD(1,rho,1,i,j,k,comp);
   else
      ugradu = assemble_advection_Upwind(1,rho,1,i,j,k,comp);

   if ( AdvectionTimeAccuracy == 1 ) {
      value = ugradu;
   } else {
      value = 1.5*ugradu - 0.5*UF->DOF_value(i,j,k,comp,2);
      UF->set_DOF_value(i,j,k,comp,2,ugradu);
   }

   return(value);
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_DS_un_at_rhs (
        FV_TimeIterator const* t_it, double const& gamma)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: assemble_DS_un_at_rhs" ) ;
   size_t i, j, k;

   double dxC,xC,dyC,yC,dzC,zC;
   double pvalue =0., xvalue=0.,yvalue=0.,zvalue=0.,rhs=0.,bodyterm=0.,adv_value = 0.;
   int cpp=-1;

   // Periodic pressure gradient
   if ( UF->primary_grid()->is_periodic_flow() ) {
      cpp = UF->primary_grid()->get_periodic_flow_direction() ;
      bodyterm = UF->primary_grid()->get_periodic_pressure_drop() /
               ( UF->primary_grid()->get_main_domain_max_coordinate( cpp )
               - UF->primary_grid()->get_main_domain_min_coordinate( cpp ) ) ;
   }

   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   for (size_t comp=0;comp<nb_comps[1];comp++) {
      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) = UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) = UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
         // Compute VEC_rhs_x = rhs in x
         dxC = UF->get_cell_size( i, comp, 0 ) ;
         xC = UF->get_DOF_coordinate( i, comp, 0 ) ;
         for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
            dyC = UF->get_cell_size( j, comp, 1 ) ;
            yC = UF->get_DOF_coordinate( j, comp, 1 ) ;
            if (dim ==2 ) {
               k=0;
               // Dxx for un
               xvalue = compute_un_component(comp,i,j,k,0,3);
               // Dyy for un
               yvalue = compute_un_component(comp,i,j,k,1,1);
               // Pressure contribution
	       pvalue = compute_p_component(comp,i,j,k);
               // Advection contribution
	       adv_value = compute_adv_component(comp,i,j,k);
          
               rhs = gamma*(xvalue*dyC + yvalue*dxC) - pvalue - adv_value
                   + (UF->DOF_value( i, j, k, comp, 1 )*dxC*dyC*rho)/(t_it -> time_step());

               if ( cpp >= 0 && cpp==comp ) rhs += - bodyterm*dxC*dyC;  

               UF->set_DOF_value( i, j, k, comp, 0, rhs*(t_it -> time_step())/(dxC*dyC*rho));
            } else {
               for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
                  dzC = UF->get_cell_size( k, comp, 2 ) ;
                  zC = UF->get_DOF_coordinate( k, comp, 2 ) ;
                  // Dxx for un
                  xvalue = compute_un_component(comp,i,j,k,0,3);
                  // Dyy for un
                  yvalue = compute_un_component(comp,i,j,k,1,4);
                  // Dzz for un
                  zvalue = compute_un_component(comp,i,j,k,2,1);
                  // Pressure contribution
                  pvalue = compute_p_component(comp,i,j,k);
                  // Advection contribution
                  adv_value = compute_adv_component(comp,i,j,k);

                  rhs = gamma*(xvalue*dyC*dzC + yvalue*dxC*dzC + zvalue*dxC*dyC) - pvalue - adv_value + (UF->DOF_value( i, j, k, comp, 1 )*dxC*dyC*dzC*rho)/(t_it -> time_step());

                  if ( cpp >= 0 && cpp==comp ) rhs += - bodyterm*dxC*dyC*dzC;

		  UF->set_DOF_value( i, j, k, comp, 0, rhs*(t_it -> time_step())/(dxC*dyC*dzC*rho));
               }
            }
         }
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: Solve_i_in_jk ( FV_DiscreteField* FF, FV_TimeIterator const* t_it, size_t const& dir_i, size_t const& dir_j, size_t const& dir_k, double const& gamma, size_t const& field )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_NavierStokes:: Solve_i_in_jk" ) ;
  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  for (size_t comp=0;comp<nb_comps[field];comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_unknown_index(dir_k);
        local_max_k = max_unknown_index(dir_k);
     }

     LocalVector* VEC = GLOBAL_EQ->get_VEC(field) ;
     TDMatrix* A = GLOBAL_EQ->get_A(field);

     // Solve in i
     if ((nb_ranks_comm_i[dir_i]>1)||(is_periodic[field][dir_i] == 1)) {
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j) {
           for (size_t k=local_min_k; k <= local_max_k; ++k) {
              size_t r_index = return_row_index (FF,comp,dir_i,j,k);
              // Assemble fi and return fe for each proc locally        
              double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir_i,field);
              // Calculate Aei*ui in each proc locally
              compute_Aei_ui(A,VEC,comp,dir_i,r_index);
              // Pack Aei_ui and fe for sending it to master
              data_packing (FF,j,k,fe,comp,dir_i,field);
           }
        }
        solve_interface_unknowns ( FF, gamma, t_it, comp, dir_i,field );

     } else if (is_periodic[field][dir_i] == 0) {  // Serial mode with non-periodic condition
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j) {
           for (size_t k=local_min_k; k <= local_max_k; ++k) {
              size_t r_index = return_row_index (FF,comp,dir_i,j,k);
              double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir_i,field);
              GLOBAL_EQ->DS_NavierStokes_solver(FF,j,k,min_unknown_index(dir_i),comp,dir_i,field,r_index);
           }
        }
     }
  }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: data_packing ( FV_DiscreteField const* FF, size_t const& j, size_t const& k, double const& fe, size_t const& comp, size_t const& dir, size_t const& field)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: data_packing" ) ;
   LocalVector* VEC = GLOBAL_EQ->get_VEC(field) ;

   double *packed_data = first_pass[field][dir].send[comp][rank_in_i[dir]];

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   // Pack the data
   size_t vec_pos=0;
   if (dir == 0) {
      if (dim == 2) {
         vec_pos=j-min_unknown_index(1);
      } else {
         vec_pos=(j-min_unknown_index(1))+(max_unknown_index(1)-min_unknown_index(1)+1)*(k-min_unknown_index(2));
      }
   } else if (dir == 1) {
      if (dim == 2) {
         vec_pos=j-min_unknown_index(0);
      } else {
         vec_pos=(j-min_unknown_index(0))+(max_unknown_index(0)-min_unknown_index(0)+1)*(k-min_unknown_index(2));
      }
   } else if (dir == 2) {
      vec_pos=(j-min_unknown_index(0))+(max_unknown_index(0)-min_unknown_index(0)+1)*(k-min_unknown_index(1));
   }

   if (rank_in_i[dir] == 0) {
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_periodic[field][dir])
          packed_data[3*vec_pos+0] = VEC[dir].T[comp]->item(nb_ranks_comm_i[dir]-1);
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

   packed_data[3*vec_pos+2] = fe; // Send the fe values and 0 for last proc

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: compute_Aei_ui (struct TDMatrix* arr, struct LocalVector* VEC, size_t const& comp, size_t const& dir, size_t const& r_index)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: compute_Aei_ui" ) ;
   // create a replica of local rhs vector in local solution vector
   for (size_t i=0;i<VEC[dir].local_T[comp]->nb_rows();i++){
      VEC[dir].local_solution_T[comp]->set_item(i,VEC[dir].local_T[comp]->item(i));
   }

   // Solve for ui locally and put it in local solution vector
   GLOBAL_EQ->mod_thomas_algorithm(arr, VEC[dir].local_solution_T[comp], comp, dir,r_index);

   for (size_t i=0;i<VEC[dir].T[comp]->nb_rows();i++){
          VEC[dir].T[comp]->set_item(i,0);
   }

   // Calculate Aei*ui in each proc locally and put it in T vector
   arr[dir].ei[comp][r_index]->multiply_vec_then_add(VEC[dir].local_solution_T[comp],VEC[dir].T[comp]);

}

//---------------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_local_rhs ( size_t const& j, size_t const& k, double const& gamma, FV_TimeIterator const* t_it, size_t const& comp, size_t const& dir, size_t const& field )
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

   Solve_i_in_jk(UF,t_it,0,1,2,gamma,1);
   // Synchronize the distributed DS solution vector
   GLOBAL_EQ->synchronize_DS_solution_vec();
   // Tranfer back to field
   UF->update_free_DOFs_value( 3, GLOBAL_EQ->get_solution_DS_velocity() ) ;

   Solve_i_in_jk(UF,t_it,1,0,2,gamma,1);
   // Synchronize the distributed DS solution vector
   GLOBAL_EQ->synchronize_DS_solution_vec();
   // Tranfer back to field
   if (dim == 2) {
      UF->update_free_DOFs_value( 0 , GLOBAL_EQ->get_solution_DS_velocity() ) ;
   } else if (dim == 3) {
      UF->update_free_DOFs_value( 4 , GLOBAL_EQ->get_solution_DS_velocity() ) ;
   }


   if (dim == 3) {
      Solve_i_in_jk(UF,t_it,2,0,1,gamma,1);
      // Synchronize the distributed DS solution vector
      GLOBAL_EQ->synchronize_DS_solution_vec();
      // Tranfer back to field
      UF->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_DS_velocity() ) ;
   }

}

//---------------------------------------------------------------------------
double
DDS_NavierStokes:: pressure_local_rhs ( size_t const& j, size_t const& k, FV_TimeIterator const* t_it, size_t const& dir)
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

   size_t i,pos;
   FV_SHIFT_TRIPLET shift = PF->shift_staggeredToCentered() ;

   // Compute VEC_rhs_x = rhs in x
   double xhr,xright,yhr,yright,dx,zhr,zright;
   double fe=0.;
   double xvalue = 0.,yvalue=0.,zvalue=0.,value=0.;

   // Vector for fi
   LocalVector* VEC = GLOBAL_EQ->get_VEC(0);

   for (i=min_unknown_index(dir);i<=max_unknown_index(dir);++i) {
      dx = PF->get_cell_size( i, 0, dir );
      if (dir == 0) {
         // Dxx for un
         xhr= UF->get_DOF_coordinate( shift.i+i,0, 0 ) - UF->get_DOF_coordinate( shift.i+i-1, 0, 0 ) ;
         xright = UF->DOF_value( shift.i+i, j, k, 0, 0 ) - UF->DOF_value( shift.i+i-1, j, k, 0, 0 ) ;
         xvalue = xright/xhr;

         // Dyy for un
         yhr= UF->get_DOF_coordinate( shift.j+j,1, 1 ) - UF->get_DOF_coordinate( shift.j+j-1, 1, 1 ) ;
         yright = UF->DOF_value( i, shift.j+j, k, 1, 0 ) - UF->DOF_value( i, shift.j+j-1, k, 1, 0 ) ;
         yvalue = yright/yhr;

         if (dim == 3) {
            // Dzz for un
            zhr= UF->get_DOF_coordinate( shift.k+k,2, 2 ) - UF->get_DOF_coordinate( shift.k+k-1, 2, 2 ) ;
            zright = UF->DOF_value( i, j, shift.k+k, 2, 0 ) - UF->DOF_value( i, j, shift.k+k-1, 2, 0 ) ;
            zvalue = zright/zhr;
         }

         // Assemble the bodyterm
         if (dim == 2) {
            value = -(rho*(xvalue + yvalue)*dx)/(t_it -> time_step());
         } else {
            value = -(rho*(xvalue + yvalue + zvalue)*dx)/(t_it -> time_step());
         }
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
DDS_NavierStokes:: NS_pressure_update ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: NS_pressure_update" ) ;

  double gamma=mu/2.0;

  Solve_i_in_jk (PF,t_it,0,1,2,gamma,0);
  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec_P();
  // Tranfer back to field
  PF->update_free_DOFs_value( 1, GLOBAL_EQ->get_solution_DS_pressure() ) ;
  
  Solve_i_in_jk (PF,t_it,1,0,2,gamma,0);
  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec_P();
  // Tranfer back to field
  PF->update_free_DOFs_value( 1, GLOBAL_EQ->get_solution_DS_pressure() ) ;
 
  if (dim == 3) { 
     Solve_i_in_jk (PF,t_it,2,0,1,gamma,0);
     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec_P();
     // Tranfer back to field
     PF->update_free_DOFs_value( 1, GLOBAL_EQ->get_solution_DS_pressure() ) ;
  }
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: NS_final_step ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL(
   "DDS_NavierStokes:: NS_final_step" ) ;

   size_t i, j, k;

   double value=0.;

   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   double xhr,xright,xvalue1=0.,xvalue2=0.,xvalue=0.;
   double yhr,yright,yvalue1=0.,yvalue2=0.,yvalue=0.;
   FV_SHIFT_TRIPLET shift = PF->shift_staggeredToCentered() ;

   // Get local min and max indices 
   // When we are running in parallel, the unknowns in the overlapping region are not solved. So we need to include 
   // them here by calling get_min_index_unknown_on_proc() instead of get_min_index_unknown_handled_by_proc().

   for (size_t l=0;l<dim;++l) {
       min_unknown_index(l) = PF->get_min_index_unknown_on_proc( 0, l ) ;
       max_unknown_index(l) = PF->get_max_index_unknown_on_proc( 0, l ) ;
   }

   for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
      for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
         if (dim ==2 ) {
            k=0;
          
            xhr= UF->get_DOF_coordinate( shift.i+i,0, 0 ) - UF->get_DOF_coordinate( shift.i+i-1, 0, 0 ) ;
            // Divergence of un+1 (x component)
            xright = UF->DOF_value( shift.i+i, j, k, 0, 0 ) - UF->DOF_value( shift.i+i-1, j, k, 0, 0 ) ;
            xvalue1 = xright/xhr;
            // Divergence of un (x component)
            xright = UF->DOF_value( shift.i+i, j, k, 0, 1 ) - UF->DOF_value( shift.i+i-1, j, k, 0, 1 ) ;
            xvalue2 = xright/xhr;      
            xvalue = xvalue1+xvalue2;

            yhr= UF->get_DOF_coordinate( shift.j+j,1, 1 ) - UF->get_DOF_coordinate( shift.j+j-1, 1, 1 ) ;
            // Divergence of un+1 (y component)
            yright = UF->DOF_value( i, shift.j+j, k, 1, 0 ) - UF->DOF_value( i, shift.j+j-1, k, 1, 0 ) ;
            yvalue1 = yright/yhr;
            // Divergence of un (y component)
            yright = UF->DOF_value( i, shift.j+j, k, 1, 1 ) - UF->DOF_value( i, shift.j+j-1, k, 1, 1 ) ;
            yvalue2 = yright/yhr;
            yvalue = yvalue1+yvalue2;

            // Assemble the bodyterm
            value = PF->DOF_value( i, j, k, 0, 0 ) + PF->DOF_value( i, j, k, 0, 1 ) - 0.5*kai*mu*(xvalue + yvalue);

            PF->set_DOF_value( i, j, k, 0, 0, value);
        } else {
            for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
               xhr= UF->get_DOF_coordinate( shift.i+i,0, 0 ) - UF->get_DOF_coordinate( shift.i+i-1, 0, 0 ) ;
               // Divergence of un+1 (x component)
               xright = UF->DOF_value( shift.i+i, j, k, 0, 0 ) - UF->DOF_value( shift.i+i-1, j, k, 0, 0 ) ;
               xvalue1 = xright/xhr;
               // Divergence of un (x component)
               xright = UF->DOF_value( shift.i+i, j, k, 0, 1 ) - UF->DOF_value( shift.i+i-1, j, k, 0, 1 ) ;
               xvalue2 = xright/xhr;
               xvalue = xvalue1+xvalue2;

               yhr= UF->get_DOF_coordinate( shift.j+j,1, 1 ) - UF->get_DOF_coordinate( shift.j+j-1, 1, 1 ) ;
               // Divergence of un+1 (y component)            
               yright = UF->DOF_value( i, shift.j+j, k, 1, 0 ) - UF->DOF_value( i, shift.j+j-1, k, 1, 0 ) ;
               yvalue1 = yright/yhr;
               // Divergence of un (y component)
               yright = UF->DOF_value( i, shift.j+j, k, 1, 1 ) - UF->DOF_value( i, shift.j+j-1, k, 1, 1 ) ;
               yvalue2 = yright/yhr;
               yvalue = yvalue1+yvalue2;

               double zhr,zright,zvalue1=0.,zvalue2=0.,zvalue=0.;
               zhr= UF->get_DOF_coordinate( shift.k+k,2, 2 ) - UF->get_DOF_coordinate( shift.k+k-1, 2, 2 ) ;
               // Divergence of un+1 (z component)
               zright = UF->DOF_value( i, j, shift.k+k, 2, 0 ) - UF->DOF_value( i, j, shift.k+k-1, 2, 0 ) ;
               zvalue1 = zright/zhr;
               // Divergence of un (z component)
               zright = UF->DOF_value( i, j, shift.k+k, 2, 1 ) - UF->DOF_value( i, j, shift.k+k-1, 2, 1 ) ;
               zvalue2 = zright/zhr;            
               zvalue = zvalue1+zvalue2;

               // Assemble the bodyterm
               value = PF->DOF_value( i, j, k, 0, 0 ) + PF->DOF_value( i, j, k, 0, 1 ) - 0.5*kai*mu*(xvalue + yvalue+ zvalue);
               PF->set_DOF_value( i, j, k, 0, 0, value);
            }
         }
      }
   }

   // Propagate values to the boundaries depending on BC conditions
   PF->set_neumann_DOF_values();

   UF->copy_DOFs_value( 0, 1 );

}

//----------------------------------------------------------------------
void
DDS_NavierStokes::write_output_field(FV_DiscreteField const* FF, size_t const& field)
//----------------------------------------------------------------------
{
  ofstream outputFile ;

  std::ostringstream os2;
  os2 << "/home/goyal001/Documents/Computing/MAC-Test/DS_results/output_" << my_rank << ".csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());

  size_t i,j,k;
  outputFile << "x,y,z,par_ID,void_frac,left,lv,right,rv,bottom,bov,top,tv,behind,bev,front,fv" << endl;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  PartInput solid = GLOBAL_EQ->get_solid(field);
  NodeProp node = GLOBAL_EQ->get_node_property(field);
  BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(field,0);

  for (size_t comp=0;comp<nb_comps[field];comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) - 1;
        max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) + 1;
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_unknown_index(2);
        local_max_k = max_unknown_index(2);
     }

     for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
        double xC = FF->get_DOF_coordinate( i, comp, 0 ) ;
        for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
           double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;
           for (k=local_min_k;k<=local_max_k;++k) {
              double zC = 0.;
              if (dim == 3) zC = FF->get_DOF_coordinate( k, comp, 2 ) ;
              size_t p = return_node_index(FF,comp,i,j,k);
              size_t id = node.parID[comp]->item(p);
              size_t voidf = node.void_frac[comp]->item(p);

              outputFile << xC << "," << yC << "," << zC << "," << id << "," << voidf;
              for (size_t dir = 0; dir < dim; dir++) {
                  for (size_t off = 0; off < 2; off++) {
                      outputFile << "," << b_intersect[dir].offset[comp]->item(p,off) << "," << b_intersect[dir].value[comp]->item(p,off);
                  }
              }
              outputFile << endl;
           }
        }
     }
  }
  outputFile.close();
}

//----------------------------------------------------------------------
double
DDS_NavierStokes::get_velocity_divergence(void)
//----------------------------------------------------------------------
{

  size_t i,j,k;
  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  double dux, duy, dx, dy, dz, duz=0.;
  double div_velocity = 0.;
  double cell_div=0.,max_divu=0.;

  FV_SHIFT_TRIPLET shift = PF->shift_staggeredToCentered() ;
  for (size_t l=0;l<dim;++l)
    min_unknown_index(l) =
	PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
  for (size_t l=0;l<dim;++l)
    max_unknown_index(l) =
	PF->get_max_index_unknown_handled_by_proc( 0, l ) ;

  for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
  {
    dx = PF->get_cell_size( i, 0, 0 );
    for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
    {
      dy = PF->get_cell_size( j, 0, 1 );
      if( dim == 2 )
      {
        k=0;

        // Divergence of u (x component)
        dux = UF->DOF_value( shift.i+i, j, k, 0, 0 ) 
		- UF->DOF_value( shift.i+i-1, j, k, 0, 0 ) ;

        // Divergence of u (y component)
        duy = UF->DOF_value( i, shift.j+j, k, 1, 0 ) 
		- UF->DOF_value( i, shift.j+j-1, k, 1, 0 ) ;

        cell_div = dux * dy + duy * dx;		
	max_divu = MAC::max( MAC::abs(cell_div) / ( dx * dy ), max_divu );
        div_velocity += cell_div * cell_div / ( dx * dy );
      }
      else
      {
        for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	{
          dz = PF->get_cell_size( k, 0, 2 );
	  
	  // Divergence of u (x component)
          dux = UF->DOF_value( shift.i+i, j, k, 0, 0 ) 
	  	- UF->DOF_value( shift.i+i-1, j, k, 0, 0 ) ;

          // Divergence of u (y component)
          duy = UF->DOF_value( i, shift.j+j, k, 1, 0 ) 
	  	- UF->DOF_value( i, shift.j+j-1, k, 1, 0 ) ;
           
          // Divergence of u(z component)
          duz = UF->DOF_value( i, j, shift.k+k, 2, 0 ) 
	  	- UF->DOF_value( i, j, shift.k+k-1, 2, 0 ) ;
          

          cell_div = dux * dy * dz + duy * dx * dz + duz * dx * dy ;		
	  max_divu = MAC::max( MAC::abs(cell_div) / ( dx * dy * dz ), 
	  	max_divu );
          div_velocity += cell_div * cell_div / ( dx * dy * dz );
        }
      }
    }
  }

  FV_Mesh const* primary_mesh = UF->primary_grid() ;
/*  double domain_measure = dim == 2 ? 
  	primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
  	* primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 1 ):
	primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
	* ( primary_mesh->get_main_domain_max_coordinate(2)
		- primary_mesh->get_main_domain_min_coordinate(2) );*/
  
  div_velocity = pelCOMM->sum( div_velocity ) ;
//  div_velocity = MAC::sqrt( div_velocity / domain_measure );
  div_velocity = MAC::sqrt( div_velocity );
  max_divu = pelCOMM->max( max_divu ) ;
  if ( my_rank == is_master )
    MAC::out() << "Norm L2 div(u) = "<< MAC::doubleToString( ios::scientific, 12, div_velocity ) << " Max div(u) = " << MAC::doubleToString( ios::scientific, 12, max_divu ) << endl;

  return(max_divu);
}




//----------------------------------------------------------------------
void
DDS_NavierStokes::output_L2norm_pressure( size_t const& level )
//----------------------------------------------------------------------
{
  size_t i,j,k;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  double dx,dy;
  double L2normP = 0.;
  double cell_P=0.,max_P=0.;

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  }

  for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
     for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
        if (dim ==2 ) {
           k=0;
           dx = PF->get_cell_size( i,0, 0 );
           dy = PF->get_cell_size( j,0, 1 );
           cell_P = PF->DOF_value( i, j, k, 0, level );
           max_P = MAC::max( MAC::abs(cell_P), max_P );
           L2normP += cell_P * cell_P * dx * dy;
        } else {
           double dz=0.;
           for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {          
              dx = PF->get_cell_size( i,0, 0 );
              dy = PF->get_cell_size( j,0, 1 );
              dz = PF->get_cell_size( k,0, 2 );
              cell_P = PF->DOF_value( i, j, k, 0, level );
              max_P = MAC::max( MAC::abs(cell_P), max_P );
              L2normP += cell_P * cell_P * dx * dy * dz;
           }
        }
     }
  }

  FV_Mesh const* primary_mesh = UF->primary_grid() ;
/*  double domain_measure = dim == 2 ? 
                          primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
                        * primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 1 ):
                          primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
                      * ( primary_mesh->get_main_domain_max_coordinate(2)
                        - primary_mesh->get_main_domain_min_coordinate(2) );*/
  
  L2normP = pelCOMM->sum( L2normP ) ;
//  L2normP = MAC::sqrt( L2normP / domain_measure );
  L2normP = MAC::sqrt( L2normP );  
  max_P = pelCOMM->max( max_P ) ;
  if ( my_rank == is_master )
      MAC::out()<< "Norm L2 P = "<< MAC::doubleToString( ios::scientific, 12, L2normP ) << " Max P = " << MAC::doubleToString( ios::scientific, 12, max_P ) << endl;
      
}

//----------------------------------------------------------------------
void
DDS_NavierStokes::output_L2norm_velocity( size_t const& level )
//----------------------------------------------------------------------
{
  size_t i,j,k;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  double dx,dy;

  for (size_t comp=0;comp<nb_comps[1];comp++) {
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) = UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     double L2normU = 0.;
     double cell_U=0.,max_U=0.;

     for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
        for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
           if (dim ==2 ) {
              k=0;
              dx = UF->get_cell_size( i,comp, 0 );
              dy = UF->get_cell_size( j,comp, 1 );
              cell_U = UF->DOF_value( i, j, k, comp, level );
              max_U = MAC::max( MAC::abs(cell_U), max_U );
              L2normU += cell_U * cell_U * dx * dy;
           } else {
              double dz=0.;
              for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {          
                 dx = UF->get_cell_size( i,comp, 0 );
                 dy = UF->get_cell_size( j,comp, 1 );
                 dz = UF->get_cell_size( k,comp, 2 );
                 cell_U = UF->DOF_value( i, j, k, comp, level );
                 max_U = MAC::max( MAC::abs(cell_U), max_U );
                 L2normU += cell_U * cell_U * dx * dy * dz;
              }
           }
        }
     }

     FV_Mesh const* primary_mesh = UF->primary_grid() ;
/*     double domain_measure = dim == 2 ? 
                             primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
                           * primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 1 ):
                             primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
                         * ( primary_mesh->get_main_domain_max_coordinate(2)
                           - primary_mesh->get_main_domain_min_coordinate(2) );*/
    
     L2normU = pelCOMM->sum( L2normU ) ;
//     L2normP = MAC::sqrt( L2normP / domain_measure );
     L2normU = MAC::sqrt( L2normU );  
     max_U = pelCOMM->max( max_U ) ;
     if ( my_rank == is_master )
        MAC::out()<< "Component: "<<comp<< " Norm L2 U = "<< MAC::doubleToString( ios::scientific, 12, L2normU ) << " Max U = " << MAC::doubleToString( ios::scientific, 12, max_U ) << endl;
  }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: create_DDS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: create_DDS_subcommunicators" ) ;

   int color = 0, key = 0;
   int const* MPI_coordinates_world = UF->primary_grid()->get_MPI_coordinates() ;
   int const* MPI_number_of_coordinates = UF->primary_grid()->get_domain_decomposition() ;

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
      color = MPI_coordinates_world[1] + MPI_coordinates_world[2]*MPI_number_of_coordinates[1] ;
      key = MPI_coordinates_world[0];
      // Split by direction in x
      processor_splitting (color,key,0);

      // Assign color and key for splitting in y
      color = MPI_coordinates_world[2] + MPI_coordinates_world[0]*MPI_number_of_coordinates[2];
      key = MPI_coordinates_world[1];
      // Split by direction in y
      processor_splitting (color,key,1);

      // Assign color and key for splitting in y
      color = MPI_coordinates_world[0] + MPI_coordinates_world[1]*MPI_number_of_coordinates[0];;
      key = MPI_coordinates_world[2];

      // Split by direction in y
      processor_splitting (color,key,2);
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: processor_splitting ( int const& color, int const& key, size_t const& dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: processor_splitting" ) ;

   MPI_Comm_split(MPI_COMM_WORLD, color, key, &DDS_Comm_i[dir]);
   MPI_Comm_size( DDS_Comm_i[dir], &nb_ranks_comm_i[dir] ) ;
   MPI_Comm_rank( DDS_Comm_i[dir], &rank_in_i[dir] ) ;

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: allocate_mpi_variables (FV_DiscreteField const* FF, size_t const& field)
//---------------------------------------------------------------------------
{

   for (size_t dir = 0; dir < dim; dir++) {
      first_pass[field][dir].size = new int [nb_comps[field]];
      second_pass[field][dir].size = new int [nb_comps[field]];
      for (size_t comp = 0; comp < nb_comps[field]; comp++) {
         size_t local_min_j=0, local_max_j=0;
         size_t local_min_k=0, local_max_k=0;

         // Get local min and max indices
         size_t_vector min_unknown_index(dim,0);
         size_t_vector max_unknown_index(dim,0);
         for (size_t l=0;l<dim;++l) {
            min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
            max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
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
            first_pass[field][dir].size[comp] = 3*local_length_j*local_length_k;
            second_pass[field][dir].size[comp] = 2*local_length_j*local_length_k;
         }
      }
   }

   // Array declarations
   for (size_t dir = 0; dir < dim; dir++) {
      first_pass[field][dir].send = new double** [nb_comps[field]];
      first_pass[field][dir].receive = new double** [nb_comps[field]];
      second_pass[field][dir].send = new double** [nb_comps[field]];
      second_pass[field][dir].receive = new double** [nb_comps[field]];
      for (size_t comp = 0; comp < nb_comps[field]; comp++) {
         first_pass[field][dir].send[comp] = new double* [nb_ranks_comm_i[dir]];
         first_pass[field][dir].receive[comp] = new double* [nb_ranks_comm_i[dir]];
         second_pass[field][dir].send[comp] = new double* [nb_ranks_comm_i[dir]];
         second_pass[field][dir].receive[comp] = new double* [nb_ranks_comm_i[dir]];
         for (size_t i = 0; i < nb_ranks_comm_i[dir]; i++) {
            first_pass[field][dir].send[comp][i] = new double[first_pass[field][dir].size[comp]];
            first_pass[field][dir].receive[comp][i] = new double[first_pass[field][dir].size[comp]];
            second_pass[field][dir].send[comp][i] = new double[second_pass[field][dir].size[comp]];
            second_pass[field][dir].receive[comp][i] = new double[second_pass[field][dir].size[comp]];
         }
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: deallocate_mpi_variables (size_t const& field)
//---------------------------------------------------------------------------
{
   // Array declarations
   for (size_t dir = 0; dir < dim; dir++) {
      for (size_t comp = 0; comp < nb_comps[field]; comp++) {
         for (size_t i = 0; i < nb_ranks_comm_i[dir]; i++) {
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

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: free_DDS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: free_DDS_subcommunicators" ) ;
}

//----------------------------------------------------------------------
double 
DDS_NavierStokes:: assemble_advection_Upwind( 
  size_t const& advecting_level, double const& coef, size_t const& advected_level,
  size_t const& i, size_t const& j, size_t const& k, size_t const& component ) const
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
      flux = (fto - fbo) * dxC + (fri - fle) * dyC;
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
