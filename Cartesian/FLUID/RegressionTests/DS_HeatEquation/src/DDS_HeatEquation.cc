#include <DDS_HeatEquation.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DiscreteField.hh>
#include <DDS_HeatEquationSystem.hh>
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


DDS_HeatEquation const* DDS_HeatEquation::PROTOTYPE
                                                 = new DDS_HeatEquation() ;


//---------------------------------------------------------------------------
DDS_HeatEquation:: DDS_HeatEquation( void )
//--------------------------------------------------------------------------
   : FV_OneStepIteration( "DDS_HeatEquation" )
   , ComputingTime("Solver")
{
   MAC_LABEL( "DDS_HeatEquation:: DDS_HeatEquation" ) ;

}




//---------------------------------------------------------------------------
DDS_HeatEquation*
DDS_HeatEquation:: create_replica( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   DDS_HeatEquation* result =
                        new DDS_HeatEquation( a_owner, dom, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;

}

//---------------------------------------------------------------------------
DDS_HeatEquation:: DDS_HeatEquation( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp )
   , ComputingTime("Solver")
   , TF ( dom->discrete_field( "temperature" ) )
   , TF_ERROR( 0 )
   , TF_DS_ERROR( 0 )
   , GLOBAL_EQ( 0 )
   , peclet( 1. )
   , b_bodyterm( false )
   , is_firstorder( false )
{
   MAC_LABEL( "DDS_HeatEquation:: DDS_HeatEquation" ) ;

   // Call of MAC_Communicator routine to set the rank of each proces and
   // the number of processes during execution
   pelCOMM = MAC_Exec::communicator();
   my_rank = pelCOMM->rank();
   nb_procs = pelCOMM->nb_ranks();
   is_master = 0;
   is_iperiodic[0] = false;
   is_iperiodic[1] = false;
   is_iperiodic[2] = false;

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
   dim = TF->primary_grid()->nb_space_dimensions() ;
   nb_comps = TF->nb_components() ;

   if ( dim == 1 )
   {
     string error_message="Space dimension should either 2 or 3";
     MAC_Error::object()->raise_bad_data_value(exp,
	"nb_space_dimensions",
	error_message );
   }


   // Create the Direction Splitting subcommunicators
   create_DDS_subcommunicators();

   // Read Peclet number
   if ( exp->has_entry( "Peclet" ) )
   {
     peclet = exp->double_data( "Peclet" ) ;
     exp->test_data( "Peclet", "Peclet>0." ) ;
   }

   // Read with or without body term
   if ( exp->has_entry( "BodyTerm" ) )
     b_bodyterm = exp->bool_data( "BodyTerm" ) ;

   // Implement first order/second order in time
   if ( exp->has_entry( "FirstOrder" ) )
     is_firstorder = exp->bool_data( "FirstOrder" ) ;

   // Periodic boundary condition check
   periodic_comp = TF->primary_grid()->get_periodic_directions();
   is_iperiodic[0] = periodic_comp->operator()( 0 );
   is_iperiodic[1] = periodic_comp->operator()( 1 );
   if(dim >2) {
      is_iperiodic[2] = periodic_comp->operator()( 2 );
   }

   // Build the matrix system
   MAC_ModuleExplorer* se =
	exp->create_subexplorer( 0,
	"DDS_HeatEquationSystem" ) ;
   GLOBAL_EQ = DDS_HeatEquationSystem::create( this, se, TF ) ;
   se->destroy() ;


   // Duplicates field for error calculation in case of sinusoidal
   // solution with body term 3.pi^2.sin(pi.x).sin(pi.y).sin(pi.z)
   // on a [0:1]x[0:1]x[0:1] domain
   if ( b_bodyterm )
   {
     const_cast<FV_DomainAndFields*>(dom)->duplicate_field(
   	"temperature", "tf_error" ) ;

     TF_ERROR = dom->discrete_field( "tf_error" ) ;
     TF_ERROR->set_BC_values_modif_status( true ) ;
   }

   const_cast<FV_DomainAndFields*>(dom)->duplicate_field(
    "temperature", "tf_ds_error" ) ;

     TF_DS_ERROR = dom->discrete_field( "tf_ds_error" ) ;
     TF_DS_ERROR->set_BC_values_modif_status( true ) ;


   // Timing routines
   if ( my_rank == is_master )
   {
     SCT_insert_app("Matrix_Assembly&Initialization");
     SCT_insert_app("Matrix_Solution");
     SCT_insert_app("DS_Solution");
     SCT_insert_app("Solver first step");
     SCT_insert_app("Solver x solution");
     SCT_insert_app("Transfer x solution");
     SCT_insert_app("Solver y solution");
     SCT_insert_app("Transfer y solution");
     SCT_insert_app("Solver z solution");
     SCT_insert_app("Transfer z solution");
     SCT_get_elapsed_time("Objects_Creation");
   }
}




//---------------------------------------------------------------------------
DDS_HeatEquation:: ~DDS_HeatEquation( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: ~DDS_HeatEquation" ) ;

   free_DDS_subcommunicators() ;

}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: do_before_time_stepping( FV_TimeIterator const* t_it,
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: do_before_time_stepping" ) ;

   start_total_timer( "DDS_HeatEquation:: do_before_time_stepping" ) ;

   if ( my_rank == is_master ) SCT_set_start("Matrix_Assembly&Initialization");

   FV_OneStepIteration::do_before_time_stepping( t_it, basename ) ;

   // Direction splitting
   // Assemble 1D tridiagonal matrices and schur complement calculation
   assemble_temperature_and_schur(t_it);

   GLOBAL_EQ->initialize_temperature();

   if ( my_rank == is_master ) SCT_get_elapsed_time( "Matrix_Assembly&Initialization" );

   stop_total_timer() ;

}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: do_before_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: do_before_inner_iterations_stage" ) ;

   start_total_timer( "DDS_HeatEquation:: do_before_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_before_inner_iterations_stage( t_it ) ;

   // Perform matrix level operations before each time step
   GLOBAL_EQ->at_each_time_step( );

   stop_total_timer() ;

}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: do_one_inner_iteration( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: do_one_inner_iteration" ) ;
   MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "DDS_HeatEquation:: do_one_inner_iteration" ) ;
   start_solving_timer() ;

   if ( my_rank == is_master ) SCT_set_start("DS_Solution");
   
   // Solve heat equation using direction splitting 
   HeatEquation_DirectionSplittingSolver(t_it);
   
   if ( my_rank == is_master ) SCT_get_elapsed_time( "DS_Solution" );

   stop_solving_timer() ;
   stop_total_timer() ;

}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: do_after_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: do_after_inner_iterations_stage" ) ;

   start_total_timer( "DDS_HeatEquation:: do_after_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_after_inner_iterations_stage( t_it ) ;

   // Compute temperature change over the time step
   double temperature_time_change = GLOBAL_EQ->compute_temperature_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master )
     cout << "Temperature change = " <<
     	MAC::doubleToString( ios::scientific, 5, temperature_time_change )
	<< endl;

   stop_total_timer() ;

}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: do_after_time_stepping" ) ;

   //write_output_field();

   // Elapsed time by sub-problems
   if ( my_rank == is_master )
   {
     double cputime = CT_get_elapsed_time();
     cout << endl << "Full problem" << endl;
     write_elapsed_time_smhd(cout,cputime,"Computation time");
     SCT_get_summary(cout,cputime);
   }
   //DS_error_with_analytical_solution(TF,TF_DS_ERROR);
   GLOBAL_EQ->display_debug();
   output_l2norm();
}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: do_additional_savings" ) ;

   start_total_timer( "DDS_HeatEquation:: do_additional_savings" ) ;

   // if ( b_bodyterm )
   // {
   //   error_with_analytical_solution( TF, TF_ERROR );
   // }

   stop_total_timer() ;

}




//----------------------------------------------------------------------
void
DDS_HeatEquation::write_output_field()
//----------------------------------------------------------------------
{

  double value;
  ofstream outputFile ;

  outputFile.open("/home/path/Results/output_5.txt", std::ios_base::out | std::ios_base::trunc) ;
  size_t i,j,k;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  for(size_t comp=0;comp<nb_comps;comp++)
  {
    // Get local min and max indices
    for (size_t l=0;l<dim;++l)
      min_unknown_index(l) =
       TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
    for (size_t l=0;l<dim;++l)
      max_unknown_index(l) =
       TF->get_max_index_unknown_handled_by_proc( comp, l ) ;

    for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
    {
      for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
      {
        if(dim ==2 )
        {
          k=0;
          value = TF->DOF_value( i, j, k, comp, 0 );
             outputFile << std::fixed << std::setprecision(8) << value ;
             outputFile<<endl;
        }
        else{
          for (k=min_unknown_index(2);k<=max_unknown_index(2);++k){
            value = TF->DOF_value( i, j, k, comp, 0 );
             outputFile << std::fixed << std::setprecision(8) << value;
             outputFile<<endl;
          }
        }
      }
    }
  }
  outputFile.close();
}

//---------------------------------------------------------------------------
double
DDS_HeatEquation:: bodyterm_value ( double xC, double yC, double zC) 
//---------------------------------------------------------------------------
{
   double bodyterm = 0.;
   if (b_bodyterm) {
      if (dim == 2) {
         bodyterm = 2. * MAC::pi() * MAC::pi() * MAC::sin( MAC::pi() * xC )
                                               * MAC::sin( MAC::pi() * yC );
      } else if (dim == 3) {
         bodyterm = 3. * MAC::pi() * MAC::pi() * MAC::sin( MAC::pi() * xC )
                                               * MAC::sin( MAC::pi() * yC ) 
                                               * MAC::sin( MAC::pi() * zC );
      }
   }

   return (bodyterm);
}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: assemble_DS_un_at_rhs (
        FV_TimeIterator const* t_it, double const gamma)
//---------------------------------------------------------------------------
{
  double dxC, dyC, dzC, xC, yC, zC=0.;
  double xvalue=0.,yvalue=0.,zvalue=0.,rhs=0., bodyterm=0.;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     size_t i, j, k;
     for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
        // Compute VEC_rhs_x = rhs in x
        dxC = TF->get_cell_size( i, comp, 0 ) ;
        xC = TF->get_DOF_coordinate( i, comp, 0 ) ;
        for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
           dyC = TF->get_cell_size( j, comp, 1 ) ;
	   yC = TF->get_DOF_coordinate( j, comp, 1 ) ;
           if (dim ==2 ) {
              k = 0;
              // Dxx for un
              xvalue = compute_un_component(comp,i,j,k,0);
              // Dyy for un
              yvalue = compute_un_component(comp,i,j,k,1);
	      // Bodyterm for rhs
	      bodyterm = bodyterm_value(xC,yC,zC);

              rhs = gamma*(xvalue*dyC + yvalue*dxC) + (TF->DOF_value( i, j, k, comp, 1 )*dxC*dyC)/(t_it -> time_step());
              TF->set_DOF_value( i, j, k, comp, 0, rhs*(t_it -> time_step())/(dxC*dyC)+gamma*bodyterm*(t_it -> time_step()));

           } else {
              for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
                 dzC = TF->get_cell_size( k, comp, 2 ) ;
	         zC = TF->get_DOF_coordinate( k, comp, 2 ) ;
                 // Dxx for un
                 xvalue = compute_un_component(comp,i,j,k,0);
                 // Dyy for un
                 yvalue = compute_un_component(comp,i,j,k,1);
                 // Dzz for un
                 zvalue = compute_un_component(comp,i,j,k,2);
	         // Bodyterm for rhs
	         bodyterm = bodyterm_value(xC,yC,zC);

                 rhs = gamma*(xvalue*dyC*dzC + yvalue*dxC*dzC + zvalue*dxC*dyC)
                               + (TF->DOF_value( i, j, k, comp, 1 )*dxC*dyC*dzC)/(t_it -> time_step());
                 TF->set_DOF_value( i, j, k, comp, 0, rhs*(t_it -> time_step())/(dxC*dyC*dzC)+gamma*bodyterm*(t_it -> time_step()));
              }
           }
        }
     }
  }
}
  


//---------------------------------------------------------------------------
double
DDS_HeatEquation:: compute_un_component ( size_t const& comp, size_t i, size_t j, size_t k, size_t const dir)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatEquation:: compute_un_component" ) ;

   double xhr,xhl,xright,xleft,yhr,yhl,yright,yleft;
   double zhr,zhl,zright,zleft, value=0.;

   if (dir == 0) {
      xhr= TF->get_DOF_coordinate( i+1,comp, 0 ) - TF->get_DOF_coordinate( i, comp, 0 ) ;
      xhl= TF->get_DOF_coordinate( i, comp, 0 ) - TF->get_DOF_coordinate( i-1, comp, 0 ) ;
      xright = TF->DOF_value( i+1, j, k, comp, 1 ) - TF->DOF_value( i, j, k, comp, 1 ) ;
      xleft = TF->DOF_value( i, j, k, comp, 1 ) - TF->DOF_value( i-1, j, k, comp, 1 ) ;

      //xvalue = xright/xhr - xleft/xhl;
      if (TF->DOF_in_domain( i-1, j, k, comp) && TF->DOF_in_domain( i+1, j, k, comp))
         value = xright/xhr - xleft/xhl;
      else if (TF->DOF_in_domain( i-1, j, k, comp))
         value = - xleft/xhl;
      else
         value = xright/xhr;
   } else if (dir == 1) {
      yhr= TF->get_DOF_coordinate( j+1,comp, 1 ) - TF->get_DOF_coordinate( j, comp, 1 ) ;
      yhl= TF->get_DOF_coordinate( j, comp, 1 ) - TF->get_DOF_coordinate( j-1, comp, 1 ) ;
      yright = TF->DOF_value( i, j+1, k, comp, 1 ) - TF->DOF_value( i, j, k, comp, 1 ) ;
      yleft = TF->DOF_value( i, j, k, comp, 1 ) - TF->DOF_value( i, j-1, k, comp, 1 ) ;

      //yvalue = yright/yhr - yleft/yhl;
      if (TF->DOF_in_domain(i, j-1, k, comp) && TF->DOF_in_domain(i, j+1, k, comp))
         value = yright/yhr - yleft/yhl;
      else if(TF->DOF_in_domain(i, j-1, k, comp))
         value = - yleft/yhl;
      else
         value = yright/yhr;
   } else if (dir == 2) {
      zhr= TF->get_DOF_coordinate( k+1,comp, 2 ) - TF->get_DOF_coordinate( k, comp, 2 ) ;
      zhl= TF->get_DOF_coordinate( k, comp, 2 ) - TF->get_DOF_coordinate( k-1, comp, 2 ) ;
      zright = TF->DOF_value( i, j, k+1, comp, 1 ) - TF->DOF_value( i, j, k, comp, 1 ) ;
      zleft = TF->DOF_value( i, j, k, comp, 1 ) - TF->DOF_value( i, j, k-1, comp, 1 ) ;

      //zvalue = zright/zhr - zleft/zhl;
      if (TF->DOF_in_domain(i, j, k-1, comp) && TF->DOF_in_domain(i, j, k+1, comp))
         value = zright/zhr - zleft/zhl;
      else if(TF->DOF_in_domain(i, j, k-1, comp))
         value = - zleft/zhl;
      else
         value = zright/zhr;
   }

   return(value);
	   
}

//---------------------------------------------------------------------------
double
DDS_HeatEquation:: assemble_temperature_matrix (
  FV_DiscreteField const* FF,
  FV_TimeIterator const* t_it,
  double gamma,
  size_t const& comp,
  size_t const dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: assemble_temperature_matrix" ) ;

   if ( my_rank == is_master ) cout << "Temperature matrix in " << dir << endl;

   // Parameters
   double dxr,dxl,xR,xL,xC,right,left,center = 0. ;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
	min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
	max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   // Perform assembling
   size_t m, i;
   TDMatrix* A = GLOBAL_EQ-> get_A();

   double Aee_diagcoef=0.;

   for (m=0,i=min_unknown_index(dir);i<=max_unknown_index(dir);++i,++m) {
       xC = FF->get_DOF_coordinate( i, comp, dir ) ;
       xR = FF->get_DOF_coordinate( i+1, comp, dir ) ;
       xL = FF->get_DOF_coordinate( i-1, comp, dir ) ;

       dxr = xR - xC;
       dxl = xC - xL;

       right = -gamma/(dxr);
       left = -gamma/(dxl);
       center = - (right+left);

       // add unsteady term
       size_t k;
       double value;
       double unsteady_term = (FF->get_cell_size(i,comp,dir))/(t_it->time_step());


       if (dim == 2) {
	  k = 0;
       } else {
	  k = min_unknown_index(2);
       }

       value = center;

       value = value + unsteady_term;

       bool r_bound = false;
       bool l_bound = false;
       // All the proc will have open right bound, except last proc for non periodic systems
       if ((is_iperiodic[dir] != 1) && (rank_in_i[dir] == nb_ranks_comm_i[dir]-1)) r_bound = true;
       // All the proc will have open left bound, except first proc for non periodic systems
       if ((is_iperiodic[dir] != 1) && (rank_in_i[dir] == 0)) l_bound = true;

       // Set Aie, Aei and Ae 
       if ((!l_bound) && (i == min_unknown_index(dir))) {
          // Periodic boundary condition at minimum unknown index
          // First proc has non zero value in Aie,Aei for first & last index
	  if (rank_in_i[dir] == 0) {
             A[dir].ie[comp]->set_item(m,nb_ranks_comm_i[dir]-1,left);
      	     A[dir].ei[comp]->set_item(nb_ranks_comm_i[dir]-1,m,right);			
	  } else {
             A[dir].ie[comp]->set_item(m,rank_in_i[dir]-1,left);
	     A[dir].ei[comp]->set_item(rank_in_i[dir]-1,m,right);			
	  }
       }

       if ((!r_bound) && (i == max_unknown_index(dir))) {
          // Periodic boundary condition at maximum unknown index
          // For last index, Aee comes from this proc as it is interface unknown wrt this proc
          A[dir].ie[comp]->set_item(m-1,rank_in_i[dir],right);			
          Aee_diagcoef = value;
          A[dir].ei[comp]->set_item(rank_in_i[dir],m-1,left);
       }

       // Set Aii_sub_diagonal
       if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_iperiodic[dir] != 1)) {
          if (i > min_unknown_index(dir)) A[dir].ii_sub[comp]->set_item(m-1,left);
       } else {
          if (i<max_unknown_index(dir)) {
             if (i>min_unknown_index(dir)) {
                A[dir].ii_sub[comp]->set_item(m-1,left);
	     }
	  }
       }

       // Set Aii_super_diagonal
       if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_iperiodic[dir] != 1)) {
          if (i < max_unknown_index(dir)) A[dir].ii_super[comp]->set_item(m,right);
       } else {
          if (i < max_unknown_index(dir)-1) {
             A[dir].ii_super[comp]->set_item(m,right);
          }
       }

       // Set Aii_main_diagonal
       if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_iperiodic[dir] != 1)) {
          A[dir].ii_main[comp]->set_item(m,value);
       } else {
          if (i<max_unknown_index(dir)) {
             A[dir].ii_main[comp]->set_item(m,value);
          }
       }
   }  // End of for loop

   GLOBAL_EQ->pre_thomas_treatment(comp,dir,A);

   return(Aee_diagcoef);
}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: assemble_schur_matrix (struct TDMatrix *A, size_t const& comp, size_t const dir, double Aee_diagcoef )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: assemble_schur_matrix" ) ;
   // Compute the product matrix for each proc

   if ( my_rank == is_master ) cout << "Schur matrix in " << dir << endl;

   if (nb_ranks_comm_i[dir]>1) {

      ProdMatrix* Ap = GLOBAL_EQ->get_Ap();

      GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir);

      LA_SeqMatrix* product_matrix = Ap[dir].ei_ii_ie[comp];
      LA_SeqMatrix* receive_matrix = product_matrix->create_copy(this,product_matrix);

      if ( rank_in_i[dir] == 0 ) {
         A[dir].ee[comp]->set_item(0,0,Aee_diagcoef);
         for (size_t i=1;i<nb_ranks_comm_i[dir];++i) {

             // Create the container to receive
             size_t nbrows = product_matrix->nb_rows();
             size_t nb_received_data = nbrows*nbrows+1;
             double * received_data = new double [nb_received_data];

             // Receive the data
             static MPI_Status status ;
             MPI_Recv( received_data, nb_received_data, MPI_DOUBLE, i, 0,
                                       DDS_Comm_i[dir], &status ) ;

             // Transfer the received data to the receive matrix
             for (size_t k=0;k<nbrows;k++) {
                for (size_t j=0;j<nbrows;j++) {
                   // Assemble the global product matrix by adding contributions from all the procs
                   receive_matrix->add_to_item(k,j,received_data[k*(nbrows)+j]);
                }
             }

   	     if (is_iperiodic[dir] == 0) {
                if (i<nb_ranks_comm_i[dir]-1) {
                   // Assemble the global Aee matrix
                   // No periodic condition in x. So no fe contribution from last proc
                   A[dir].ee[comp]->set_item(i,i,received_data[nb_received_data-1]);
                }
             } else{
                // Assemble the global Aee matrix
                // Periodic condition in x. So there is fe contribution from last proc
                A[dir].ee[comp]->set_item(i,i,received_data[nb_received_data-1]);
             }
	     delete [] received_data;

         }
      } else {
         // Create the packed data container

         size_t nbrows = product_matrix->nb_rows();
         size_t nb_send_data = nbrows*nbrows+1;
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

      if (rank_in_i[dir] == 0) {
	 TDMatrix* Schur = GLOBAL_EQ-> get_Schur();
         size_t nb_row = Schur[dir].ii_main[comp]->nb_rows();
         for (int p = 0; p < nb_row; p++) {
            Schur[dir].ii_main[comp]->set_item(p,A[dir].ee[comp]->item(p,p)-receive_matrix->item(p,p));
            if (p < nb_row-1) Schur[dir].ii_super[comp]->set_item(p,-receive_matrix->item(p,p+1));
            if (p > 0) Schur[dir].ii_sub[comp]->set_item(p-1,-receive_matrix->item(p,p-1));
	    // In case of periodic and multi-processor, there will be a variant of Tridiagonal matrix instead of normal format
            if (is_iperiodic[dir] == 1) {
               Schur[dir].ie[comp]->set_item(p,0,-receive_matrix->item(p,nb_row)); 
               Schur[dir].ei[comp]->set_item(0,p,-receive_matrix->item(nb_row,p)); 
	    }
	 }

	 // Pre-thomas treatment on Schur complement
	 GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur);

	 // In case of periodic and multi-processor, there will be a variant of Tridiagonal matrix instead of normal format
	 // So, Schur complement of Schur complement is calculated
	 if (is_iperiodic[dir] == 1) {
            Schur[dir].ee[comp]->set_item(0,0,A[dir].ee[comp]->item(nb_row,nb_row)-receive_matrix->item(nb_row,nb_row));

	    ProdMatrix* SchurP = GLOBAL_EQ->get_SchurP();
            GLOBAL_EQ->compute_product_matrix_interior(Schur,SchurP,comp,0,dir);

            TDMatrix* DoubleSchur = GLOBAL_EQ-> get_DoubleSchur();
            size_t nb_row = DoubleSchur[dir].ii_main[comp]->nb_rows();
            DoubleSchur[dir].ii_main[comp]->set_item(0,Schur[dir].ee[comp]->item(0,0)-SchurP[dir].ei_ii_ie[comp]->item(0,0));
	 }

      }
   } else if (is_iperiodic[dir] == 1) {
      // Condition for single processor in any direction with periodic boundary conditions
      ProdMatrix* Ap = GLOBAL_EQ->get_Ap();
      GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir);

      LA_SeqMatrix* product_matrix = Ap[dir].ei_ii_ie[comp];
      LA_SeqMatrix* receive_matrix = product_matrix->create_copy(this,product_matrix);

      A[dir].ee[comp]->set_item(0,0,Aee_diagcoef);

      TDMatrix* Schur = GLOBAL_EQ-> get_Schur();
      size_t nb_row = Schur[dir].ii_main[comp]->nb_rows();
      for (int p = 0; p < nb_row; p++) {
         Schur[dir].ii_main[comp]->set_item(p,A[dir].ee[comp]->item(p,p)-receive_matrix->item(p,p));
         if (p < nb_row-1) Schur[dir].ii_super[comp]->set_item(p,-receive_matrix->item(p,p+1));
         if (p > 0) Schur[dir].ii_sub[comp]->set_item(p-1,-receive_matrix->item(p,p-1));
      }
      GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur);
   }
}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: assemble_temperature_and_schur ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: assemble_temperature_and_schur" ) ;

   double gamma;

   TDMatrix* A = GLOBAL_EQ-> get_A();

   // Assemble temperature matrix and schur complement for each component
   for (size_t comp=0;comp<nb_comps;comp++) {
        if (is_firstorder) {
                gamma = 1.0/peclet;
        } else {
                gamma = 1.0/2.0/peclet;
        }

	for (size_t i = 0; i<dim; i++) {
            double Aee_diagcoef = assemble_temperature_matrix (TF,t_it,gamma,comp,i);
	    assemble_schur_matrix(A,comp,i,Aee_diagcoef);
	}
   }
}

//---------------------------------------------------------------------------
double
DDS_HeatEquation:: assemble_local_rhs ( size_t const& j, size_t const& k, double gamma, FV_TimeIterator const* t_it, size_t const& comp, size_t const dir)
//---------------------------------------------------------------------------
{
   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc(comp,l) ;
     max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc(comp,l) ;
   }

   size_t i,pos;
   int m;

   // Compute VEC_rhs_x = rhs in x
   double dC,hr=0,hl=0,right=0,left=0,xC,yC,zC=0;
   double fe=0.;

   // Vector for fi
   LocalVector* VEC = GLOBAL_EQ->get_VEC();

   for (i=min_unknown_index(dir);i<=max_unknown_index(dir);++i) {
     double value=0.;
     pos = i - min_unknown_index(dir);

     // Get contribution of un
     hl= TF->get_DOF_coordinate(i,comp,dir) - TF->get_DOF_coordinate(i-1,comp,dir) ;
     hr= TF->get_DOF_coordinate(i+1,comp,dir) - TF->get_DOF_coordinate(i,comp,dir) ;

     dC = TF->get_cell_size(i,comp,dir) ;

     if (!is_firstorder) {
	// x direction
        if (dir == 0) {
           value = compute_un_component(comp,i,j,k,dir);
	// y direction
        } else if (dir == 1) {
           value = compute_un_component(comp,j,i,k,dir);
	// z direction
        } else if (dir == 2) {
           value = compute_un_component(comp,j,k,i,dir);
        }
     } else {
	if ((b_bodyterm) && (dir==0)) {
	   xC = TF->get_DOF_coordinate( i, comp, 0 ) ;
	   yC = TF->get_DOF_coordinate( j, comp, 1 ) ;
	   if (dim == 3) zC = TF->get_DOF_coordinate( k, comp, 2 ) ;
	   // Add bodyterm for first step of Crank_Nicolson scheme
           value = -bodyterm_value(xC,yC,zC)*dC/gamma;
	}
     }

     double temp_val=0.;
     if (dir == 0) {
        temp_val = (TF->DOF_value(i,j,k,comp,0)*dC)/(t_it->time_step()) - gamma*value;
     } else if (dir == 1) {
        temp_val = (TF->DOF_value(j,i,k,comp,0)*dC)/(t_it->time_step()) - gamma*value;
     } else if (dir == 2) {
        temp_val = (TF->DOF_value(j,k,i,comp,0)*dC)/(t_it->time_step()) - gamma*value;
     }

     if (is_iperiodic[dir] == 0) {
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

   if ( TF->DOF_in_domain(ii,jj,kk,comp))
      if ( TF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
         double ai = 1/(TF->get_DOF_coordinate(m+1,comp,dir) - TF->get_DOF_coordinate(m,comp,dir));
         double dirichlet_value = TF->DOF_value(ii,jj,kk,comp,1) ;
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

   if ( TF->DOF_in_domain(ii,jj,kk,comp))
      if ( TF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
      double ai = 1/(TF->get_DOF_coordinate(m,comp,dir) - TF->get_DOF_coordinate(m-1,comp,dir));
      double dirichlet_value = TF->DOF_value(ii,jj,kk,comp,1) ;
      VEC[dir].local_T[comp]->add_to_item( VEC[dir].local_T[comp]->nb_rows()-1 , + gamma * ai * dirichlet_value );
   }

   return fe;
}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: compute_Aei_ui (struct TDMatrix* arr, struct LocalVector* VEC, size_t const& comp, size_t const dir)
//---------------------------------------------------------------------------
{
   // create a replica of local rhs vector in local solution vector
   for (size_t i=0;i<VEC[dir].local_T[comp]->nb_rows();i++){
      VEC[dir].local_solution_T[comp]->set_item(i,VEC[dir].local_T[comp]->item(i));
   }

   // Solve for ui locally and put it in local solution vector
   GLOBAL_EQ->mod_thomas_algorithm(arr, VEC[dir].local_solution_T[comp], comp, dir);

   for (size_t i=0;i<VEC[dir].T[comp]->nb_rows();i++){
          VEC[dir].T[comp]->set_item(i,0);
   }

   // Calculate Aei*ui in each proc locally and put it in T vector
   arr[dir].ei[comp]->multiply_vec_then_add(VEC[dir].local_solution_T[comp],VEC[dir].T[comp]);

}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: data_packing ( size_t const& j, size_t const& k, double * packed_data, double fe, size_t const& comp, size_t const dir)
//---------------------------------------------------------------------------
{
   LocalVector* VEC = GLOBAL_EQ->get_VEC() ;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
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
      if(is_iperiodic[dir])
          packed_data[3*vec_pos+0] = VEC[dir].T[comp]->item(nb_ranks_comm_i[dir]-1);
      else
          packed_data[3*vec_pos+0] = 0;

      packed_data[3*vec_pos+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);

   } else if (rank_in_i[dir] == nb_ranks_comm_i[dir]-1) {
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_iperiodic[dir])
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
DDS_HeatEquation:: unpack_compute_ue_pack(size_t const& comp, double ** all_received_data, double * packed_data, double ** all_send_data, size_t const dir, size_t p)
//---------------------------------------------------------------------------
{
   LocalVector* VEC = GLOBAL_EQ->get_VEC() ;

   size_t nb_interface_unknowns = VEC[dir].T[comp]->nb_rows();
   for (size_t i=0;i<nb_interface_unknowns;i++) {
       VEC[dir].T[comp]->set_item(i,0);
       VEC[dir].interface_T[comp]->set_item(i,0);
   }

   // If periodic in x, first proc contributes to last interface unknown
   if (is_iperiodic[dir])
      VEC[dir].T[comp]->set_item(nb_ranks_comm_i[dir]-1,packed_data[3*p]);
   VEC[dir].T[comp]->set_item(0,packed_data[3*p+1]);
   VEC[dir].interface_T[comp]->set_item(0,packed_data[3*p+2]);

   // Vec_temp might contain previous values
   for (size_t i=1;i<nb_ranks_comm_i[dir];i++) {
      if (i!=nb_ranks_comm_i[dir]-1) {
         VEC[dir].T[comp]->add_to_item(i-1,all_received_data[i][3*p]);
         VEC[dir].T[comp]->add_to_item(i,all_received_data[i][3*p+1]);
         VEC[dir].interface_T[comp]->set_item(i,all_received_data[i][3*p+2]);  // Assemble the interface rhs fe
      } else {
         if (is_iperiodic[dir] ==0) {
            VEC[dir].T[comp]->add_to_item(i-1,all_received_data[i][3*p]);
         } else {
            VEC[dir].T[comp]->add_to_item(i-1,all_received_data[i][3*p]);
            // If periodic in x, last proc has an interface unknown
            VEC[dir].T[comp]->add_to_item(i,all_received_data[i][3*p+1]);
            VEC[dir].interface_T[comp]->set_item(i,all_received_data[i][3*p+2]);
         }
      }
   }

   for (size_t i=0;i<nb_interface_unknowns;i++) {
      VEC[dir].interface_T[comp]->set_item(i,VEC[dir].interface_T[comp]->item(i)-VEC[dir].T[comp]->item(i)); // Get fe - Aei*xi to solve for ue
   }

   // Solve for ue (interface unknowns) in the master proc
   DS_interface_unknown_solver(VEC[dir].interface_T[comp], comp, dir);

   // Pack the interface_rhs_x into the appropriate send_data
   for (size_t i=1;i<nb_ranks_comm_i[dir];++i) {
      if (i!=nb_ranks_comm_i[dir]-1) {
         all_send_data[i][2*p+0] = VEC[dir].interface_T[comp]->item(i-1);
         all_send_data[i][2*p+1] = VEC[dir].interface_T[comp]->item(i);
      } else {
         all_send_data[i][2*p+0] = VEC[dir].interface_T[comp]->item(i-1);
         if (is_iperiodic[dir])
            all_send_data[i][2*p+1] = VEC[dir].interface_T[comp]->item(i);
         else
            all_send_data[i][2*p+1] = 0;
      }
   }
   
}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: unpack_ue(size_t const& comp, double * received_data, size_t const dir, int p)
//---------------------------------------------------------------------------
{
   LocalVector* VEC = GLOBAL_EQ->get_VEC() ;

   if (rank_in_i[dir] != nb_ranks_comm_i[dir]-1) {
      VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,received_data[2*p]);
      VEC[dir].interface_T[comp]->set_item(rank_in_i[dir],received_data[2*p+1]);
   } else {
      if (is_iperiodic[dir] ==0) {
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,received_data[2*p]);
      } else {
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,received_data[2*p]);
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir],received_data[2*p+1]);
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: solve_interface_unknowns ( double * packed_data, size_t nb_received_data, double gamma,  FV_TimeIterator const* t_it, size_t const& comp, size_t const dir)
//---------------------------------------------------------------------------
{
   size_t i,j,p;
   size_t k = 0;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   TDMatrix* A = GLOBAL_EQ->get_A();
   LocalVector* VEC = GLOBAL_EQ->get_VEC() ;

   // Array declaration for sending data from master to all slaves
   double ** all_send_data = new double* [nb_ranks_comm_i[dir]];
   size_t nb_send_data=0, local_length_j=0, local_length_k=0;
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

   if (dim != 3) {
      nb_send_data = 2*local_length_j;
   } else if (dim == 3) {
      nb_send_data = 2*local_length_j*local_length_k;
   }

   for (p = 0; p < nb_ranks_comm_i[dir]; ++p) {
      all_send_data[p] = new double[nb_send_data];
   }


   // Send and receive the data first pass
   if ( rank_in_i[dir] == 0 ) {
      // Array declaration for receiving data from all slaves
      double ** all_received_data = new double* [nb_ranks_comm_i[dir]];
      for(p = 0; p < nb_ranks_comm_i[dir]; ++p) {
          all_received_data[p] = new double[nb_received_data];
      }
      // Receiving data from all the slave procs iff multi processors are used
      if (nb_ranks_comm_i[dir] != 1) {
         for (i=1;i<nb_ranks_comm_i[dir];++i) {
            static MPI_Status status;
            MPI_Recv( all_received_data[i], nb_received_data, MPI_DOUBLE, i, 0,
                                             DDS_Comm_i[dir], &status ) ;
         }
      }

      // Solve system of interface unknowns for each y
      if (dim == 2) {
         for (j=local_min_j;j<=local_max_j;j++) {

     	    p = j-local_min_j;

	    unpack_compute_ue_pack(comp,all_received_data,packed_data,all_send_data,dir,p); 

            // Need to have the original rhs function assembled for corrosponding j,k pair
            double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir);

            // Setup RHS = fi - Aie*ue for solving ui
            A[dir].ie[comp]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir),comp,dir);

         }
      } else {
         for (k=local_min_k;k<=local_max_k;k++) {
            for (j=local_min_j;j<=local_max_j;j++) {

   	       p = (j-local_min_j)+local_length_j*(k-local_min_k);

	       unpack_compute_ue_pack(comp,all_received_data,packed_data,all_send_data,dir,p); 

               // Need to have the original rhs function assembled for corrosponding j,k pair
               double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir);

               // Setup RHS = fi - Aie*ue for solving ui
               A[dir].ie[comp]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

               // Solve ui and transfer solution into distributed vector
               GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir),comp,dir);
            }
         }
      }

      for (p = 0; p < nb_ranks_comm_i[dir]; ++p) delete [] all_received_data[p];
      delete [] all_received_data;

   } else {
      // Send the packed data to master
      MPI_Send( packed_data, nb_received_data, MPI_DOUBLE, 0, 0, DDS_Comm_i[dir] ) ;
   }

   // Send the data from master iff multi processor are used
   if (nb_ranks_comm_i[dir] != 1) {
      if ( rank_in_i[dir] == 0 ) {
         for (i=1;i<nb_ranks_comm_i[dir];++i) {
            MPI_Send( all_send_data[i], nb_send_data, MPI_DOUBLE, i, 0, DDS_Comm_i[dir] ) ;
         }
      } else {
         // Create the container to receive the ue
         if (dim==2)
            nb_received_data = 2*local_length_j;
         else
            nb_received_data = 2*local_length_j*local_length_k;

         double * received_data = new double [nb_received_data];

         static MPI_Status status ;
         MPI_Recv( received_data, nb_received_data, MPI_DOUBLE, 0, 0,
                             DDS_Comm_i[dir], &status ) ;

         // Solve the system of equations in each proc

         if (dim == 2) {
            for (j = local_min_j;j<=local_max_j;j++) {
               p = j-local_min_j;

               unpack_ue(comp,received_data,dir,p);

               // Need to have the original rhs function assembled for corrosponding j,k pair
               double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir);

               // Setup RHS = fi - Aie*xe for solving ui
               A[dir].ie[comp]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

               // Solve ui and transfer solution into distributed vector
               GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir),comp,dir);
            }
         } else {
            for (k = local_min_k;k<=local_max_k;k++) {
               for (j = local_min_j;j<=local_max_j;j++) {
                  p = (j-local_min_j)+local_length_j*(k-local_min_k);

	          unpack_ue(comp,received_data,dir,p);

                  // Need to have the original rhs function assembled for corrosponding j,k pair
                  double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir);

                  // Setup RHS = fi - Aie*xe for solving ui
                  A[dir].ie[comp]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

                  // Solve ui and transfer solution into distributed vector
                  GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir),comp,dir);
               }
            }
         }
         delete [] received_data;
      }
   }

   for (p = 0; p < nb_ranks_comm_i[dir]; ++p) delete [] all_send_data[p];
   delete [] all_send_data;

}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: Solve_i_in_jk ( FV_TimeIterator const* t_it, double gamma, size_t const dir_i, size_t const dir_j, size_t const dir_k )
//---------------------------------------------------------------------------
{
  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     size_t nb_send_data=0;
     if (dim == 2) {
        nb_send_data = 3*(max_unknown_index(dir_j)-min_unknown_index(dir_j)+1);
     } else if (dim == 3) {
        nb_send_data = 3*(max_unknown_index(dir_j)-min_unknown_index(dir_j)+1)*(max_unknown_index(dir_k)-min_unknown_index(dir_k)+1);
        local_min_k = min_unknown_index(dir_k);
        local_max_k = max_unknown_index(dir_k);
     }

     double * packed_data = new double[nb_send_data];
     LocalVector* VEC = GLOBAL_EQ->get_VEC() ;
     TDMatrix* A = GLOBAL_EQ->get_A();

     // Solve in i
     if ((nb_ranks_comm_i[dir_i]>1)||(is_iperiodic[dir_i] == 1)) {
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j) {
	   for (size_t k=local_min_k; k <= local_max_k; ++k) {
              // Assemble fi and return fe for each proc locally        
              double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir_i);
              // Calculate Aei*ui in each proc locally
              compute_Aei_ui(A,VEC,comp,dir_i);
              // Pack Aei_ui and fe for sending it to master
              data_packing (j,k,packed_data,fe,comp,dir_i);
	   }
        }
        solve_interface_unknowns ( packed_data, nb_send_data, gamma, t_it,comp,dir_i );

     } else if (is_iperiodic[dir_i] == 0) {  // Serial mode with non-periodic condition
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j) {
           for (size_t k=local_min_k; k <= local_max_k; ++k) {
              double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir_i);
              GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir_i),comp,dir_i);
           }
        }
     }
     delete [] packed_data;
  }
}

//----------------------------------------------------------------------
void
DDS_HeatEquation::DS_interface_unknown_solver(LA_SeqVector* interface_rhs, size_t const& comp, size_t const dir)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: DS_interface_unknown_solver" ) ;

   TDMatrix* Schur = GLOBAL_EQ-> get_Schur();

   // Condition for variant of Tridiagonal Schur complement in Perioidic direction with multi-processor 
   if ((is_iperiodic[dir] == 1) && (nb_ranks_comm_i[dir] != 1)) {
      LocalVector* Schur_VEC = GLOBAL_EQ->get_Schur_VEC();
      TDMatrix* DoubleSchur = GLOBAL_EQ-> get_DoubleSchur();
      
      // Transfer interface_rhs to Schur VEC (i.e. S_fi and S_fe)
      size_t nrows = Schur_VEC[dir].local_T[comp]->nb_rows();
      for (size_t i = 0; i < nrows; i++) {
          Schur_VEC[dir].local_T[comp]->set_item(i,interface_rhs->item(i));
      }
      Schur_VEC[dir].interface_T[comp]->set_item(0,interface_rhs->item(nrows));

      // Calculate Sei*(Sii)-1*S_fi
      compute_Aei_ui(Schur,Schur_VEC,comp,dir);

      // Calculate S_fe - Sei*(Sii)-1*S_fi
      Schur_VEC[dir].interface_T[comp]->set_item(0,Schur_VEC[dir].interface_T[comp]->item(0)-Schur_VEC[dir].T[comp]->item(0));

      // Calculate S_ue, using Schur complement of Schur complement
      GLOBAL_EQ->mod_thomas_algorithm(DoubleSchur, Schur_VEC[dir].interface_T[comp], comp, dir);

      // Calculate S_fi-Sie*S_ue
      Schur[dir].ie[comp]->multiply_vec_then_add(Schur_VEC[dir].interface_T[comp],Schur_VEC[dir].local_T[comp],-1.0,1.0);

      // Calculate S_ui
      GLOBAL_EQ->mod_thomas_algorithm(Schur, Schur_VEC[dir].local_T[comp], comp, dir);

      // Transfer back the solution to interface_rhs
      for (size_t i = 0; i < nrows; i++) {
          interface_rhs->set_item(i,Schur_VEC[dir].local_T[comp]->item(i));
      }
      interface_rhs->set_item(nrows,Schur_VEC[dir].interface_T[comp]->item(0));
   } else {
      GLOBAL_EQ->mod_thomas_algorithm(Schur, interface_rhs, comp, dir);
   }

}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: HeatEquation_DirectionSplittingSolver ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatEquation:: HeatEquation_DirectionSplittingSolver" ) ;

  double gamma=1.0/peclet;

  TF->copy_DOFs_value( 0, 1 );
  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  if (!is_firstorder) {
     // First Equation
     if ( my_rank == is_master ) SCT_set_start("Solver first step");
     assemble_DS_un_at_rhs (t_it,gamma);
     if ( my_rank == is_master ) SCT_get_elapsed_time("Solver first step");
     // Update gamma based for invidual direction
     gamma = 1.0/2.0/peclet;
  }

  if ( my_rank == is_master ) SCT_set_start("Solver x solution");
  // Solve x-direction(i.e. 0) in y(i.e. 1) and z(i.e. 2)
  Solve_i_in_jk (t_it,gamma,0,1,2);
  if ( my_rank == is_master ) SCT_get_elapsed_time("Solver x solution");
  if ( my_rank == is_master ) SCT_set_start("Transfer x solution");
  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec();
  // Tranfer back to field
  TF->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_DS_temperature() ) ;
  if ( my_rank == is_master ) SCT_get_elapsed_time("Transfer x solution");


  if ( my_rank == is_master ) SCT_set_start("Solver y solution");
  // Solve y-direction(i.e. 1) in x(i.e. 0) and z(i.e. 2)
  Solve_i_in_jk (t_it,gamma,1,0,2);
  if ( my_rank == is_master ) SCT_get_elapsed_time("Solver y solution");
  if ( my_rank == is_master ) SCT_set_start("Transfer y solution");
  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec();
  // Tranfer back to field
  TF->update_free_DOFs_value( 0 , GLOBAL_EQ->get_solution_DS_temperature() ) ;
  if ( my_rank == is_master ) SCT_get_elapsed_time("Transfer y solution");

  if (dim == 3) {
     if ( my_rank == is_master ) SCT_set_start("Solver z solution");
     // Solve z-direction(i.e. 2) in x(i.e. 0) and y(i.e. 1)
     Solve_i_in_jk (t_it,gamma,2,0,1);
     if ( my_rank == is_master ) SCT_get_elapsed_time("Solver z solution");
     if ( my_rank == is_master ) SCT_set_start("Transfer z solution");
     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec();
     // Tranfer back to field
     TF->update_free_DOFs_value( 0 , GLOBAL_EQ->get_solution_DS_temperature() ) ;
     if ( my_rank == is_master ) SCT_get_elapsed_time("Transfer z solution");
  }
}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: output_l2norm ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatEquation:: output_l2norm" ) ;

   // Parameters
   size_t local_number,i,j,k;

   LA_Vector const* DS_Solution = GLOBAL_EQ->get_solution_DS_temperature();

   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   for (size_t comp=0;comp<nb_comps;comp++) {
      double computed_DS_field=0., computed_DS_L2 = 0.;

      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
         for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {

            if (dim == 2) {
               k=0;
               local_number = TF->DOF_local_number( i, j, k, comp );
               computed_DS_field = DS_Solution->item(local_number);

               if ( TF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                  computed_DS_L2 += computed_DS_field*computed_DS_field
                               * TF->get_cell_measure( i, j, k, comp ) ;
               }
            } else {
               for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
                  local_number = TF->DOF_local_number( i, j, k, comp );
                  computed_DS_field = TF->DOF_value( i, j, k, comp, 0 ) ;

                  if ( TF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                     computed_DS_L2 += computed_DS_field*computed_DS_field
                                  * TF->get_cell_measure( i, j, k, comp ) ;
                  }
               }
            }
         }
      }

      computed_DS_L2 = pelCOMM->sum( computed_DS_L2 ) ;
      computed_DS_L2 = MAC::sqrt(computed_DS_L2);

      if (my_rank == is_master) {
         cout << "L2 Norm for component "<<comp<<" with DS = " << std::fixed << std::setprecision(10) << computed_DS_L2<<endl;
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: DS_error_with_analytical_solution ( FV_DiscreteField const* FF,FV_DiscreteField* FF_ERROR )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatEquation:: DS_error_with_analytical_solution" ) ;

   // Parameters
   size_t i,j,k;

   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   for (size_t comp=0;comp<nb_comps;comp++) {

      double x,y,z,analytical_solution=0., computed_DS_field=0.,error_L2 = 0. ,denom = 0., norm_error = 0.;
      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
	 min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
         x = FF->get_DOF_coordinate( i, comp, 0 ) ;
         for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
            y = FF->get_DOF_coordinate( j, comp, 1 ) ;
            if (dim == 2) {
               k=0;

               analytical_solution = MAC::sin( MAC::pi() * x )* MAC::sin( MAC::pi() * y ) ;
               computed_DS_field = FF->DOF_value( i, j, k, comp, 0 ) ;

               //Get error between computed solutions with Direction splitting
               FF_ERROR->set_DOF_value( i, j, k, comp, 0, MAC::abs( computed_DS_field - analytical_solution) ) ;

               if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                  error_L2 += MAC::sqr( computed_DS_field - analytical_solution)
                                         * FF->get_cell_measure( i, j, k, comp ) ;
	       }
               denom += MAC::sqr( analytical_solution) * FF->get_cell_measure( i, j, k, comp ) ;
            } else {
               for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
                  z = FF->get_DOF_coordinate( k, comp, 2 ) ;
                  analytical_solution = MAC::sin( MAC::pi() * x )* MAC::sin( MAC::pi() * y ) * MAC::sin( MAC::pi() * z ) ;
                  computed_DS_field = FF->DOF_value( i, j, k, comp, 0 ) ;
 
       		  //Get error between computed solutions with Direction splitting
	    	  FF_ERROR->set_DOF_value( i, j, k, comp, 0, MAC::abs( computed_DS_field- analytical_solution ) ) ;

                  if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                     error_L2 += MAC::sqr( computed_DS_field - analytical_solution)
                                            * FF->get_cell_measure( i, j, k, comp ) ;
		  }
                  denom += MAC::sqr( analytical_solution) * FF->get_cell_measure( i, j, k, comp ) ;
               }
            }
         }
      }

      error_L2 = pelCOMM->sum( error_L2 ) ;
      error_L2 = MAC::sqrt(error_L2);
      denom = pelCOMM->sum( denom ) ;
      denom = MAC::sqrt(denom);
      norm_error = error_L2/denom;

      if (my_rank == is_master) {
         cout << "L2 Error with analytical solution for "<<comp<<" DS (Normalized) = " << norm_error<<endl;
         cout << "L2 Error with analytical solution for "<<comp<<" DS = " << error_L2<<endl;
      }
   }
}




//---------------------------------------------------------------------------
void
DDS_HeatEquation:: create_DDS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: create_DDS_subcommunicators" ) ;

   int color = 0, key = 0;
   //int const* number_of_subdomains_per_direction = TF->primary_grid()->get_domain_decomposition() ;
   int const* MPI_coordinates_world = TF->primary_grid()->get_MPI_coordinates() ;
   int const* MPI_number_of_coordinates = TF->primary_grid()->get_domain_decomposition() ;

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
DDS_HeatEquation:: processor_splitting ( int color, int key, size_t const dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: processor_splitting" ) ;

   MPI_Comm_split(MPI_COMM_WORLD, color, key, &DDS_Comm_i[dir]);
   MPI_Comm_size( DDS_Comm_i[dir], &nb_ranks_comm_i[dir] ) ;
   MPI_Comm_rank( DDS_Comm_i[dir], &rank_in_i[dir] ) ;

}

//---------------------------------------------------------------------------
void
DDS_HeatEquation:: free_DDS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: free_DDS_subcommunicators" ) ;


}
