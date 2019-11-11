#include <DDS_NavierStokesSystem.hh>
#include <LA_Matrix.hh>
#include <LA_Vector.hh>
#include <LA_Scatter.hh>
#include <LA_SeqVector.hh>
#include <LA_SeqMatrix.hh>
#include <LA_Solver.hh>
#include <LA_MatrixIterator.hh>
#include <LA_CRSmatrix.hh>
#include <intVector.hh>
#include <MAC.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Timer.hh>
#include <MAC_Vector.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <FV_DiscreteField.hh>
#include <FV_SystemNumbering.hh>
#include <FV_Mesh.hh>
#include <iostream>
#include <math.h>
// Additions
#include <stdio.h>
#include <stdlib.h>


//----------------------------------------------------------------------
DDS_NavierStokesSystem*
DDS_NavierStokesSystem:: create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_UF,
  FV_DiscreteField* mac_PF )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   MAC_CHECK_PRE( mac_UF != 0 ) ;

   DDS_NavierStokesSystem* result =
         new DDS_NavierStokesSystem( a_owner, exp, mac_UF, mac_PF ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;

   return( result ) ;

}




//----------------------------------------------------------------------
DDS_NavierStokesSystem:: DDS_NavierStokesSystem(
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_UF,
  FV_DiscreteField* mac_PF )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , UF( mac_UF )
   , PF( mac_PF )
   , MAT_velocityUnsteadyPlusDiffusion_1D_y( 0 )
   , U_is_xperiodic( false )
   , U_is_yperiodic( false )
   , U_is_zperiodic( false )
   , P_is_xperiodic( false )
   , P_is_yperiodic( false )
   , P_is_zperiodic( false )
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DDS_NavierStokesSystem" ) ;

   int const* MPI_coordinates_world = UF->primary_grid()->get_MPI_coordinates() ;
   int const* MPI_max_coordinates_world = UF->primary_grid()->get_domain_decomposition() ;

   proc_pos_in_x = MPI_coordinates_world[0];
   proc_pos_in_y = MPI_coordinates_world[1];
   proc_pos_in_z = MPI_coordinates_world[2];
   nb_procs_in_x = MPI_max_coordinates_world[0];
   nb_procs_in_y = MPI_max_coordinates_world[1];
   nb_procs_in_z = MPI_max_coordinates_world[2];

   dim = UF->primary_grid()->nb_space_dimensions() ;
   nb_comps = UF->nb_components() ;

   // Periodic boundary condition check for velocity
   U_periodic_comp = UF->primary_grid()->get_periodic_directions();
   U_is_xperiodic = U_periodic_comp->operator()( 0 );
   U_is_yperiodic = U_periodic_comp->operator()( 1 );
   if(dim >2)
      U_is_zperiodic = U_periodic_comp->operator()( 2 ); 

   // Periodic boundary condition check for pressure
   P_periodic_comp = PF->primary_grid()->get_periodic_directions();
   P_is_xperiodic = P_periodic_comp->operator()( 0 );
   P_is_yperiodic = P_periodic_comp->operator()( 1 );
   if(dim >2)
      P_is_zperiodic = P_periodic_comp->operator()( 2 ); 

   // Build the matrices & vectors
   build_system(exp) ;

   re_initialize() ;

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem:: build_system( MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: build_system" ) ;

   // velocity unsteady
   MAT_A_velocityUnsteady = LA_Matrix::make( this,
         exp->create_subexplorer( this,
	 "MAT_A_velocityUnsteady" ) ) ;
   VEC_rhs_A_velocityUnsteady =
     	MAT_A_velocityUnsteady->create_vector(this) ;

   // velocity Laplacian
   MAT_D_velocityUnsteadyPlusDiffusion = LA_Matrix::make( this,
         exp->create_subexplorer( this,
	 "MAT_D_velocityDiffusion"  ) ) ;
   VEC_rhs_D_velocityDiffusionPlusBodyTerm =
     	MAT_D_velocityUnsteadyPlusDiffusion->create_vector(this) ;

   // Unknowns vectors
   VEC_UF = MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;
   VEC_UF_previoustime =
     MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;
   VEC_UF_timechange =
     MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;

   VEC_DS_UF = MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;
   VEC_DS_UF_previoustime = MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;
   VEC_DS_UF_timechange = MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;

   VEC_DS_PF = MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;

   // Local vector

   UF_LOC = LA_SeqVector::create( this, 0 ) ;
   UF_DS_LOC = LA_SeqVector::create( this, 0 ) ;
   PF_DS_LOC = LA_SeqVector::create( this, 0 ) ;

   // Solvers
   SOLVER_velocity = LA_Solver::make( this,
               exp->create_subexplorer( this,
	       "SOLVER_velocity" ) ) ;
   SOLVER_velocity->set_initial_guess_nonzero( true );

   // velocity numbering
   UF_NUM = FV_SystemNumbering::create( this, UF ) ;
   PF_NUM = FV_SystemNumbering::create( this, PF ) ;

   // Direction splitting matrices
   MAT_velocityUnsteadyPlusDiffusion_1D_y =
   	LA_SeqMatrix::make( this, exp->create_subexplorer( this,
	 "MAT_1DLAP_generic" ) );

  Aii_x_main_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
  Aii_x_super_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
  Aii_x_sub_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
  Aii_x_mod_super_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

  U_vec_xu = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
  U_vec_xv = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ; 

	Aie_x = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
	Aei_x = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
	Aei_Aii_Aie_product_x = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;

	product_result_x = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
	Aii_Aie_product_x = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
	VEC_rhs_velocity_1D_x = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
	VEC_local_temp_x = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
	VEC_local_solution_temp_x = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
	VEC_temp_x = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
	VEC_interface_temp_x = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

	Aee_x = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
	schlur_complement_x = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
	schlur_complement_x_ref = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;

   Aii_y_main_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   Aii_y_super_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   Aii_y_sub_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   Aii_y_mod_super_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

   U_vec_yu = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   U_vec_yv = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ; 

   Aie_y = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
   Aei_y = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
   Aei_Aii_Aie_product_y = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;

   product_result_y = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   Aii_Aie_product_y = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   VEC_rhs_velocity_1D_y = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   VEC_local_temp_y = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   VEC_local_solution_temp_y = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   VEC_temp_y = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   VEC_interface_temp_y = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

   Aee_y = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
   schlur_complement_y = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
   schlur_complement_y_ref = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;

   if(dim == 3){
      Aii_z = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
      Aii_z_main_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      Aii_z_super_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      Aii_z_sub_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      Aii_z_mod_super_diagonal = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

      U_vec_zu = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      U_vec_zv = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ; 

      Aii_z_ref = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;

      Aie_z = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
      Aei_z = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
      Aei_Aii_Aie_product_z = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;

      product_result_z = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      Aii_Aie_product_z = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      VEC_rhs_velocity_1D_z = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      VEC_local_temp_z = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      VEC_local_solution_temp_z = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      VEC_temp_z = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      VEC_interface_temp_z = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

      Aee_z = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
      schlur_complement_z = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
      schlur_complement_z_ref = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
   }

   Aii_x_main_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
             ->create_vector( this ) ;
   Aii_x_super_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
       ->create_vector( this ) ;
   Aii_x_sub_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
       ->create_vector( this ) ;
   Aii_x_mod_super_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
       ->create_vector( this ) ;

  P_vec_xu = MAT_velocityUnsteadyPlusDiffusion_1D_y
             ->create_vector( this ) ;
  P_vec_xv = MAT_velocityUnsteadyPlusDiffusion_1D_y
             ->create_vector( this ) ;

  Aie_x_P =
     MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
     MAT_velocityUnsteadyPlusDiffusion_1D_y );

  Aei_x_P =
     MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
     MAT_velocityUnsteadyPlusDiffusion_1D_y );

 Aei_Aii_Aie_product_x_P =
     MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
     MAT_velocityUnsteadyPlusDiffusion_1D_y );
  product_result_x_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
     ->create_vector( this ) ;
 Aii_Aie_product_x_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
    ->create_vector( this ) ;

if(proc_pos_in_x == 0){
   Aee_x_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
   MAT_velocityUnsteadyPlusDiffusion_1D_y );
   schlur_complement_x_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
   MAT_velocityUnsteadyPlusDiffusion_1D_y );
   schlur_complement_x_ref_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
   MAT_velocityUnsteadyPlusDiffusion_1D_y );
 }

VEC_local_temp_x_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
VEC_local_solution_temp_x_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
VEC_temp_x_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
VEC_interface_temp_x_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;

Aii_y_main_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
             ->create_vector( this ) ;
Aii_y_super_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
    ->create_vector( this ) ;
Aii_y_sub_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
    ->create_vector( this ) ;
Aii_y_mod_super_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
       ->create_vector( this ) ;

P_vec_yu = MAT_velocityUnsteadyPlusDiffusion_1D_y
           ->create_vector( this ) ;
P_vec_yv = MAT_velocityUnsteadyPlusDiffusion_1D_y
           ->create_vector( this ) ;

Aie_y_P =
    MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
    MAT_velocityUnsteadyPlusDiffusion_1D_y );

Aei_y_P =
    MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
    MAT_velocityUnsteadyPlusDiffusion_1D_y );

 Aei_Aii_Aie_product_y_P =
    MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
    MAT_velocityUnsteadyPlusDiffusion_1D_y );
product_result_y_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
    ->create_vector( this ) ;
 Aii_Aie_product_y_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
       ->create_vector( this ) ;

if(proc_pos_in_y == 0){
    Aee_y_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
    MAT_velocityUnsteadyPlusDiffusion_1D_y );
    schlur_complement_y_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
    MAT_velocityUnsteadyPlusDiffusion_1D_y );
    schlur_complement_y_ref_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
    MAT_velocityUnsteadyPlusDiffusion_1D_y );
 }

VEC_local_temp_y_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
VEC_local_solution_temp_y_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
VEC_temp_y_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
VEC_interface_temp_y_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;

if(dim>2)
{
     Aii_z_main_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
         ->create_vector( this ) ;
     Aii_z_super_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
         ->create_vector( this ) ;
     Aii_z_sub_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
         ->create_vector( this ) ;
     Aii_z_mod_super_diagonal_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
            ->create_vector( this ) ;

     P_vec_zu = MAT_velocityUnsteadyPlusDiffusion_1D_y
               ->create_vector( this ) ;
     P_vec_zv = MAT_velocityUnsteadyPlusDiffusion_1D_y
               ->create_vector( this ) ;

     Aie_z_P =
        MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
        MAT_velocityUnsteadyPlusDiffusion_1D_y );

     Aei_z_P =
        MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
        MAT_velocityUnsteadyPlusDiffusion_1D_y );

      Aei_Aii_Aie_product_z_P =
        MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
        MAT_velocityUnsteadyPlusDiffusion_1D_y );
     product_result_z_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
        ->create_vector( this ) ;
      Aii_Aie_product_z_P = MAT_velocityUnsteadyPlusDiffusion_1D_y
          ->create_vector( this ) ;

     if(proc_pos_in_z == 0){
         Aee_z_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
         MAT_velocityUnsteadyPlusDiffusion_1D_y );
         schlur_complement_z_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
         MAT_velocityUnsteadyPlusDiffusion_1D_y );
         schlur_complement_z_ref_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
         MAT_velocityUnsteadyPlusDiffusion_1D_y );
     }

     VEC_local_temp_z_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
     VEC_local_solution_temp_z_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
     VEC_temp_z_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
     VEC_interface_temp_z_P = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
 }

  for (size_t comp=0;comp<nb_comps;++comp)
  {
    Aii_x_main_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
         ->create_vector( this ) ;
    Aii_x_super_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
         ->create_vector( this ) ;
    Aii_x_sub_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
         ->create_vector( this ) ;
    Aii_x_mod_super_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
         ->create_vector( this ) ;

    U_vec_xu[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
               ->create_vector( this ) ;
    U_vec_xv[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
               ->create_vector( this ) ;

    Aie_x[comp] =
       MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
       MAT_velocityUnsteadyPlusDiffusion_1D_y );

    Aei_x[comp] =
       MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
       MAT_velocityUnsteadyPlusDiffusion_1D_y );

 	  Aei_Aii_Aie_product_x[comp] =
       MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
       MAT_velocityUnsteadyPlusDiffusion_1D_y );
    product_result_x[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
       ->create_vector( this ) ;
 	  Aii_Aie_product_x[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
 		  ->create_vector( this ) ;

	  if(proc_pos_in_x == 0){
		   Aee_x[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
		   MAT_velocityUnsteadyPlusDiffusion_1D_y );
       schlur_complement_x[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
       MAT_velocityUnsteadyPlusDiffusion_1D_y );
       schlur_complement_x_ref[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
       MAT_velocityUnsteadyPlusDiffusion_1D_y );
    }


    // Direction splitting vectors
    VEC_rhs_velocity_1D_x[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
  	  ->create_vector( this ) ;

    VEC_local_temp_x[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
    VEC_local_solution_temp_x[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
    VEC_temp_x[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
    VEC_interface_temp_x[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;

    Aii_y_main_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
           ->create_vector( this ) ;
    Aii_y_super_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
        ->create_vector( this ) ;
    Aii_y_sub_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
        ->create_vector( this ) ;
    Aii_y_mod_super_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
           ->create_vector( this ) ;

    U_vec_yu[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
               ->create_vector( this ) ;
    U_vec_yv[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
               ->create_vector( this ) ;

    Aie_y[comp] =
        MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
        MAT_velocityUnsteadyPlusDiffusion_1D_y );

    Aei_y[comp] =
        MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
        MAT_velocityUnsteadyPlusDiffusion_1D_y );

     Aei_Aii_Aie_product_y[comp] =
        MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
        MAT_velocityUnsteadyPlusDiffusion_1D_y );
    product_result_y[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
        ->create_vector( this ) ;
     Aii_Aie_product_y[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
           ->create_vector( this ) ;
     VEC_rhs_velocity_1D_y[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
        ->create_vector( this ) ;

    if(proc_pos_in_y == 0){
        Aee_y[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
        MAT_velocityUnsteadyPlusDiffusion_1D_y );
        schlur_complement_y[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
        MAT_velocityUnsteadyPlusDiffusion_1D_y );
        schlur_complement_y_ref[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
        MAT_velocityUnsteadyPlusDiffusion_1D_y );
     }

    VEC_local_temp_y[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
    VEC_local_solution_temp_y[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
    VEC_temp_y[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
    VEC_interface_temp_y[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;

    if(dim>2)
    {
       VEC_rhs_velocity_1D_z[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
       ->create_vector( this ) ;

       Aii_z[comp] =
       MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
       MAT_velocityUnsteadyPlusDiffusion_1D_y );

       Aii_z_main_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
           ->create_vector( this ) ;
       Aii_z_super_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
           ->create_vector( this ) ;
       Aii_z_sub_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
           ->create_vector( this ) ;
       Aii_z_mod_super_diagonal[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
              ->create_vector( this ) ;

       U_vec_zu[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
                 ->create_vector( this ) ;
       U_vec_zv[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
                 ->create_vector( this ) ;

       Aii_z_ref[comp] =
          MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
          MAT_velocityUnsteadyPlusDiffusion_1D_y );

       Aie_z[comp] =
          MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
          MAT_velocityUnsteadyPlusDiffusion_1D_y );

       Aei_z[comp] =
          MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
          MAT_velocityUnsteadyPlusDiffusion_1D_y );

        Aei_Aii_Aie_product_z[comp] =
          MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
          MAT_velocityUnsteadyPlusDiffusion_1D_y );
       product_result_z[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
          ->create_vector( this ) ;
        Aii_Aie_product_z[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y
            ->create_vector( this ) ;

       if(proc_pos_in_z == 0){
           Aee_z[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
           MAT_velocityUnsteadyPlusDiffusion_1D_y );
           schlur_complement_z[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
           MAT_velocityUnsteadyPlusDiffusion_1D_y );
           schlur_complement_z_ref[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_copy( this,
           MAT_velocityUnsteadyPlusDiffusion_1D_y );
       }

       VEC_local_temp_z[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
       VEC_local_solution_temp_z[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
       VEC_temp_z[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
       VEC_interface_temp_z[comp] = MAT_velocityUnsteadyPlusDiffusion_1D_y->create_vector( this ) ;
    }

  }



}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem:: re_initialize( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: re_initialize" ) ;

   size_t UF_glob = UF->nb_global_unknowns() ;
   size_t UF_loc = UF->nb_local_unknowns() ;

   size_t pf_glob = PF->nb_global_unknowns() ;
   size_t pf_loc = PF->nb_local_unknowns() ;

   // velocity unsteady
   MAT_A_velocityUnsteady->re_initialize( UF_glob, UF_glob ) ;
   VEC_rhs_A_velocityUnsteady->re_initialize( UF_glob ) ;

   // velocity Laplacian
   MAT_D_velocityUnsteadyPlusDiffusion->re_initialize( UF_glob, UF_glob ) ;
   VEC_rhs_D_velocityDiffusionPlusBodyTerm->re_initialize( UF_glob ) ;

   // Unknowns vectors
   VEC_UF->re_initialize( UF_glob ) ;
   VEC_UF_previoustime->re_initialize( UF_glob ) ;
   VEC_UF_timechange->re_initialize( UF_glob ) ;

   VEC_DS_UF->re_initialize( UF_glob ) ;
   VEC_DS_UF_previoustime->re_initialize( UF_glob ) ;
   VEC_DS_UF_timechange->re_initialize( UF_glob ) ;

   VEC_DS_PF->re_initialize( pf_glob ) ;

   // Local vector
   UF_LOC->re_initialize( UF_loc ) ;
   UF_DS_LOC->re_initialize( UF_loc ) ;

   PF_DS_LOC->re_initialize( pf_loc ) ;

   // velocity numbering
   UF_NUM->define_scatter( VEC_UF ) ;
   PF_NUM->define_scatter( VEC_DS_PF ) ;

   // Initialize Direction splitting matrices & vectors for pressure 

   size_t_vector nb_unknowns_handled_by_proc( dim, 0 );
   for (size_t l=0;l<dim;++l)
     nb_unknowns_handled_by_proc( l ) =
      1 + PF->get_max_index_unknown_handled_by_proc( 0, l )
   - PF->get_min_index_unknown_handled_by_proc( 0, l ) ;

  if(P_is_xperiodic!=1 || nb_procs_in_x == 1)
  {
     if(proc_pos_in_x == nb_procs_in_x-1)
     {
        Aii_x_main_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 0 ));
        Aii_x_super_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 0 )-1);
        Aii_x_sub_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 0 )-1);
        Aii_x_mod_super_diagonal_P
        ->re_initialize(
        nb_unknowns_handled_by_proc( 0 )-1);
        
        P_vec_xu->re_initialize(
        nb_unknowns_handled_by_proc( 0 ));
        P_vec_xv->re_initialize(
        nb_unknowns_handled_by_proc( 0 ));

        Aie_x_P->re_initialize(
        nb_unknowns_handled_by_proc( 0 ),
        nb_procs_in_x-1);
        Aei_x_P->re_initialize(
        nb_procs_in_x-1,
        nb_unknowns_handled_by_proc( 0 ) );
        product_result_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 ) );

        VEC_local_temp_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )) ;
        VEC_local_solution_temp_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )) ;
      }
     else{

        Aii_x_main_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 0 )-1);
        Aii_x_super_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 0 )-2);
        Aii_x_sub_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 0 )-2);
        Aii_x_mod_super_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 0 )-2);

        Aie_x_P->re_initialize(
        nb_unknowns_handled_by_proc( 0 )-1,
        nb_procs_in_x-1 );
        Aei_x_P->re_initialize(
        nb_procs_in_x-1,
        nb_unknowns_handled_by_proc( 0 )-1 );
        product_result_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1 );
        VEC_local_temp_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
        VEC_local_solution_temp_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
      }
      Aii_Aie_product_x_P->re_initialize(
       nb_procs_in_x-1);

     Aei_Aii_Aie_product_x_P->re_initialize(
        nb_procs_in_x-1,
        nb_procs_in_x-1 );

     VEC_interface_temp_x_P->re_initialize( nb_procs_in_x-1 ) ;
     VEC_temp_x_P->re_initialize( nb_procs_in_x-1 ) ;

     if(proc_pos_in_x == 0){
       Aee_x_P->re_initialize(
          nb_procs_in_x-1,
          nb_procs_in_x-1 );
       schlur_complement_x_P->re_initialize(
            nb_procs_in_x-1,
            nb_procs_in_x-1 );
       schlur_complement_x_ref_P->re_initialize(
            nb_procs_in_x-1,
            nb_procs_in_x-1 );
      }

  }
  else{
    Aii_x_main_diagonal_P->re_initialize(
    nb_unknowns_handled_by_proc( 0 )-1);
    Aii_x_super_diagonal_P->re_initialize(
    nb_unknowns_handled_by_proc( 0 )-2);
    Aii_x_sub_diagonal_P->re_initialize(
    nb_unknowns_handled_by_proc( 0 )-2);
    Aii_x_mod_super_diagonal_P->re_initialize(
    nb_unknowns_handled_by_proc( 0 )-2);

    Aie_x_P->re_initialize(
    nb_unknowns_handled_by_proc( 0 )-1,
    nb_procs_in_x );
    Aei_x_P->re_initialize(
    nb_procs_in_x,
    nb_unknowns_handled_by_proc( 0 )-1 );
    product_result_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1 );
    VEC_local_temp_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
    VEC_local_solution_temp_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
      
    Aii_Aie_product_x_P->re_initialize(
       nb_procs_in_x);

    Aei_Aii_Aie_product_x_P->re_initialize(
        nb_procs_in_x,
        nb_procs_in_x );

    VEC_interface_temp_x_P->re_initialize( nb_procs_in_x ) ;
    VEC_temp_x_P->re_initialize( nb_procs_in_x ) ;

    if(proc_pos_in_x == 0){
       Aee_x_P->re_initialize(
          nb_procs_in_x,
          nb_procs_in_x );
       schlur_complement_x_P->re_initialize(
            nb_procs_in_x,
            nb_procs_in_x );
       schlur_complement_x_ref_P->re_initialize(
            nb_procs_in_x,
            nb_procs_in_x );
    }
  }
  if(P_is_yperiodic!=1 || nb_procs_in_y == 1)
  {
      if(proc_pos_in_y == nb_procs_in_y-1){
        Aii_y_main_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 1 ));
        Aii_y_super_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 1 )-1);
        Aii_y_sub_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 1 )-1);
        Aii_y_mod_super_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 1 )-1);

        P_vec_yu->re_initialize(
        nb_unknowns_handled_by_proc( 1 ));
        P_vec_yv->re_initialize(
        nb_unknowns_handled_by_proc( 1 ));

        Aie_y_P->re_initialize(
        nb_unknowns_handled_by_proc( 1 ),
        nb_procs_in_y-1);
        Aei_y_P->re_initialize(
        nb_procs_in_y-1,
        nb_unknowns_handled_by_proc( 1 ) );
        product_result_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 ) );

        VEC_local_temp_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )) ;
        VEC_local_solution_temp_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )) ;
     }
     else{
        Aii_y_main_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 1 )-1);
        Aii_y_super_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 1 )-2);
        Aii_y_sub_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 1 )-2);
        Aii_y_mod_super_diagonal_P->re_initialize(
        nb_unknowns_handled_by_proc( 1 )-2);

        Aie_y_P->re_initialize(
        nb_unknowns_handled_by_proc( 1 )-1,
        nb_procs_in_y-1 );
        Aei_y_P->re_initialize(
        nb_procs_in_y-1,
        nb_unknowns_handled_by_proc( 1 )-1 );
        product_result_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1 );
        VEC_local_temp_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
        VEC_local_solution_temp_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
     }

      Aii_Aie_product_y_P->re_initialize(
         nb_procs_in_y-1);

      Aei_Aii_Aie_product_y_P->re_initialize(
           nb_procs_in_y-1,
           nb_procs_in_y-1 );


      if(proc_pos_in_y == 0){
         Aee_y_P->re_initialize(
              nb_procs_in_y-1,
              nb_procs_in_y-1 );
         schlur_complement_y_P->re_initialize(
              nb_procs_in_y-1,
              nb_procs_in_y-1 );
         schlur_complement_y_ref_P->re_initialize(
              nb_procs_in_y-1,
              nb_procs_in_y-1 );
      }

     VEC_interface_temp_y_P->re_initialize( nb_procs_in_y-1 ) ;
     VEC_temp_y_P->re_initialize( nb_procs_in_y-1 ) ;
   }
   else{
      Aii_y_main_diagonal_P->re_initialize(
      nb_unknowns_handled_by_proc( 1 )-1);
      Aii_y_super_diagonal_P->re_initialize(
      nb_unknowns_handled_by_proc( 1 )-2);
      Aii_y_sub_diagonal_P->re_initialize(
      nb_unknowns_handled_by_proc( 1 )-2);
      Aii_y_mod_super_diagonal_P->re_initialize(
      nb_unknowns_handled_by_proc( 1 )-2);

      Aie_y_P->re_initialize(
      nb_unknowns_handled_by_proc( 1 )-1,
      nb_procs_in_y );
      Aei_y_P->re_initialize(
      nb_procs_in_y,
      nb_unknowns_handled_by_proc( 1 )-1 );
      product_result_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1 );
      VEC_local_temp_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
      VEC_local_solution_temp_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;

      Aii_Aie_product_y_P->re_initialize(
           nb_procs_in_y);

      Aei_Aii_Aie_product_y_P->re_initialize(
             nb_procs_in_y,
             nb_procs_in_y );


      if(proc_pos_in_y == 0){
         Aee_y_P->re_initialize(
              nb_procs_in_y,
              nb_procs_in_y );
         schlur_complement_y_P->re_initialize(
              nb_procs_in_y,
              nb_procs_in_y );
         schlur_complement_y_ref_P->re_initialize(
              nb_procs_in_y,
              nb_procs_in_y );
      }

      VEC_interface_temp_y_P->re_initialize( nb_procs_in_y ) ;
      VEC_temp_y_P->re_initialize( nb_procs_in_y) ;
   }
   if(dim>2)
   {
      if(P_is_zperiodic !=1 || nb_procs_in_z == 1)
      {

        if(proc_pos_in_z == nb_procs_in_z-1){

           Aii_z_main_diagonal_P->re_initialize(
           nb_unknowns_handled_by_proc( 2 ));
           Aii_z_super_diagonal_P->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-1);
           Aii_z_sub_diagonal_P->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-1);
           Aii_z_mod_super_diagonal_P->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-1);

           P_vec_zu->re_initialize(
           nb_unknowns_handled_by_proc( 2 ));
           P_vec_zv->re_initialize(
           nb_unknowns_handled_by_proc( 2 ));

           Aie_z_P->re_initialize(
           nb_unknowns_handled_by_proc( 2 ),
           nb_procs_in_z-1);
           Aei_z_P->re_initialize(
           nb_procs_in_z-1,
           nb_unknowns_handled_by_proc( 2 ) );
           product_result_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 ) );

           VEC_local_temp_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )) ;
           VEC_local_solution_temp_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )) ;

        }
        else{

           Aii_z_main_diagonal_P->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-1);
           Aii_z_super_diagonal_P->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-2);
           Aii_z_sub_diagonal_P->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-2);
           Aii_z_mod_super_diagonal_P->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-2);

           Aie_z_P->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-1,
           nb_procs_in_z-1 );
           Aei_z_P->re_initialize(
           nb_procs_in_z-1,
           nb_unknowns_handled_by_proc( 2 )-1 );
           product_result_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1 );
           VEC_local_temp_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
           VEC_local_solution_temp_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
        }

        if(proc_pos_in_z == 0){
            Aee_z_P->re_initialize(
                 nb_procs_in_z-1,
                 nb_procs_in_z-1 );
            schlur_complement_z_P->re_initialize(
                 nb_procs_in_z-1,
                 nb_procs_in_z-1 );
            schlur_complement_z_ref_P->re_initialize(
                 nb_procs_in_z-1,
                 nb_procs_in_z-1 );
        }

        Aii_Aie_product_z_P->re_initialize(
           nb_procs_in_z-1);

        Aei_Aii_Aie_product_z_P->re_initialize(
           nb_procs_in_z-1,
           nb_procs_in_z-1 );

        VEC_interface_temp_z_P->re_initialize( nb_procs_in_z-1 ) ;
        VEC_temp_z_P->re_initialize( nb_procs_in_z-1 ) ;
    }
    else{
         Aii_z_main_diagonal_P->re_initialize(
         nb_unknowns_handled_by_proc( 2 )-1);
         Aii_z_super_diagonal_P->re_initialize(
         nb_unknowns_handled_by_proc( 2 )-2);
         Aii_z_sub_diagonal_P->re_initialize(
         nb_unknowns_handled_by_proc( 2 )-2);
         Aii_z_mod_super_diagonal_P->re_initialize(
         nb_unknowns_handled_by_proc( 2 )-2);

         Aie_z_P->re_initialize(
         nb_unknowns_handled_by_proc( 2 )-1,
         nb_procs_in_z);
         Aei_z_P->re_initialize(
         nb_procs_in_z,
         nb_unknowns_handled_by_proc( 2 )-1 );
         product_result_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1 );
         VEC_local_temp_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
         VEC_local_solution_temp_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
         
         if(proc_pos_in_z == 0){
              Aee_z_P->re_initialize(
                   nb_procs_in_z,
                   nb_procs_in_z );
              schlur_complement_z_P->re_initialize(
                   nb_procs_in_z,
                   nb_procs_in_z );
              schlur_complement_z_ref_P->re_initialize(
                   nb_procs_in_z,
                   nb_procs_in_z );
          }

         Aii_Aie_product_z_P->re_initialize(
             nb_procs_in_z);

         Aei_Aii_Aie_product_z_P->re_initialize(
             nb_procs_in_z,
             nb_procs_in_z );

         VEC_interface_temp_z_P->re_initialize( nb_procs_in_z ) ;
         VEC_temp_z_P->re_initialize( nb_procs_in_z ) ;
    }
  }
   // Initialize Direction splitting matrices & vectors for velocity

	for (size_t comp=0;comp<nb_comps;++comp)
	{

   size_t_vector nb_unknowns_handled_by_proc( dim, 0 );
   for (size_t l=0;l<dim;++l)
     nb_unknowns_handled_by_proc( l ) =
      1 + UF->get_max_index_unknown_handled_by_proc( comp, l )
   - UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   if(U_is_xperiodic!=1 || nb_procs_in_x == 1){
      if(proc_pos_in_x == nb_procs_in_x-1)
       {

          Aii_x_main_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 ));
          Aii_x_super_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 )-1);
          Aii_x_sub_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 )-1);
          Aii_x_mod_super_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 )-1);

          U_vec_xu[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 ));
          U_vec_xv[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 ));

          Aie_x[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 ),
          nb_procs_in_x-1);
          Aei_x[comp]->re_initialize(
          nb_procs_in_x-1,
          nb_unknowns_handled_by_proc( 0 ) );
          product_result_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 ) );

          VEC_local_temp_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )) ;
          VEC_local_solution_temp_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )) ;
        }
       else{

          Aii_x_main_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 )-1);
          Aii_x_super_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 )-2);
          Aii_x_sub_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 )-2);
          Aii_x_mod_super_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 )-2);

          Aie_x[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 0 )-1,
          nb_procs_in_x-1 );
          Aei_x[comp]->re_initialize(
          nb_procs_in_x-1,
          nb_unknowns_handled_by_proc( 0 )-1 );
          product_result_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1 );
          VEC_local_temp_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
          VEC_local_solution_temp_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
        }
        Aii_Aie_product_x[comp]->re_initialize(
         nb_procs_in_x-1);

       Aei_Aii_Aie_product_x[comp]->re_initialize(
          nb_procs_in_x-1,
          nb_procs_in_x-1 );

       VEC_interface_temp_x[comp]->re_initialize( nb_procs_in_x-1 ) ;
       VEC_temp_x[comp]->re_initialize( nb_procs_in_x-1 ) ;

       if(proc_pos_in_x == 0){
         Aee_x[comp]->re_initialize(
            nb_procs_in_x-1,
            nb_procs_in_x-1 );
           schlur_complement_x[comp]->re_initialize(
                nb_procs_in_x-1,
                nb_procs_in_x-1 );
           schlur_complement_x_ref[comp]->re_initialize(
                nb_procs_in_x-1,
                nb_procs_in_x-1 );
        }

      VEC_rhs_velocity_1D_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 ) );
   }
   else{
      Aii_x_main_diagonal[comp]->re_initialize(
      nb_unknowns_handled_by_proc( 0 )-1);
      Aii_x_super_diagonal[comp]->re_initialize(
      nb_unknowns_handled_by_proc( 0 )-2);
      Aii_x_sub_diagonal[comp]->re_initialize(
      nb_unknowns_handled_by_proc( 0 )-2);
      Aii_x_mod_super_diagonal[comp]->re_initialize(
      nb_unknowns_handled_by_proc( 0 )-2);

      Aie_x[comp]->re_initialize(
      nb_unknowns_handled_by_proc( 0 )-1,
      nb_procs_in_x );
      Aei_x[comp]->re_initialize(
      nb_procs_in_x,
      nb_unknowns_handled_by_proc( 0 )-1 );
      product_result_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1 );
      VEC_local_temp_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
      VEC_local_solution_temp_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
        
      Aii_Aie_product_x[comp]->re_initialize(
         nb_procs_in_x);

      Aei_Aii_Aie_product_x[comp]->re_initialize(
          nb_procs_in_x,
          nb_procs_in_x );

      VEC_interface_temp_x[comp]->re_initialize( nb_procs_in_x ) ;
      VEC_temp_x[comp]->re_initialize( nb_procs_in_x ) ;

      if(proc_pos_in_x == 0){
         Aee_x[comp]->re_initialize(
            nb_procs_in_x,
            nb_procs_in_x );
         schlur_complement_x[comp]->re_initialize(
              nb_procs_in_x,
              nb_procs_in_x );
         schlur_complement_x_ref[comp]->re_initialize(
              nb_procs_in_x,
              nb_procs_in_x );
      }

      VEC_rhs_velocity_1D_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 ) );
   }
       
   if(U_is_yperiodic!=1 || nb_procs_in_y == 1)
   {
      if(proc_pos_in_y == nb_procs_in_y-1){
          Aii_y_main_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 ));
          Aii_y_super_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-1);
          Aii_y_sub_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-1);
          Aii_y_mod_super_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-1);

          U_vec_yu[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 ));
          U_vec_yv[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 ));

          Aie_y[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 ),
          nb_procs_in_y-1);
          Aei_y[comp]->re_initialize(
          nb_procs_in_y-1,
          nb_unknowns_handled_by_proc( 1 ) );
          product_result_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 ) );

          VEC_local_temp_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )) ;
          VEC_local_solution_temp_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )) ;
       }
       else{
          Aii_y_main_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-1);
          Aii_y_super_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-2);
          Aii_y_sub_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-2);
          Aii_y_mod_super_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-2);

          Aie_y[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-1,
          nb_procs_in_y-1 );
          Aei_y[comp]->re_initialize(
          nb_procs_in_y-1,
          nb_unknowns_handled_by_proc( 1 )-1 );
          product_result_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1 );
          VEC_local_temp_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
          VEC_local_solution_temp_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
       }

       MAT_velocityUnsteadyPlusDiffusion_1D_y->re_initialize(
          nb_unknowns_handled_by_proc( 1 ),
          nb_unknowns_handled_by_proc( 1 ) );

        Aii_Aie_product_y[comp]->re_initialize(
           nb_procs_in_y-1);

        Aei_Aii_Aie_product_y[comp]->re_initialize(
             nb_procs_in_y-1,
             nb_procs_in_y-1 );


        if(proc_pos_in_y == 0){
           Aee_y[comp]->re_initialize(
                nb_procs_in_y-1,
                nb_procs_in_y-1 );
           schlur_complement_y[comp]->re_initialize(
                nb_procs_in_y-1,
                nb_procs_in_y-1 );
           schlur_complement_y_ref[comp]->re_initialize(
                nb_procs_in_y-1,
                nb_procs_in_y-1 );
        }

       VEC_rhs_velocity_1D_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 ) );

       VEC_interface_temp_y[comp]->re_initialize( nb_procs_in_y-1 ) ;
       VEC_temp_y[comp]->re_initialize( nb_procs_in_y-1 ) ;
   }
   else
   {  
        Aii_y_main_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-1);
          Aii_y_super_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-2);
          Aii_y_sub_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-2);
          Aii_y_mod_super_diagonal[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-2);

          Aie_y[comp]->re_initialize(
          nb_unknowns_handled_by_proc( 1 )-1,
          nb_procs_in_y);
          Aei_y[comp]->re_initialize(
          nb_procs_in_y,
          nb_unknowns_handled_by_proc( 1 )-1 );
          product_result_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1 );
          VEC_local_temp_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
          VEC_local_solution_temp_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
       
       MAT_velocityUnsteadyPlusDiffusion_1D_y->re_initialize(
          nb_unknowns_handled_by_proc( 1 ),
          nb_unknowns_handled_by_proc( 1 ) );

        Aii_Aie_product_y[comp]->re_initialize(
           nb_procs_in_y);

        Aei_Aii_Aie_product_y[comp]->re_initialize(
             nb_procs_in_y,
             nb_procs_in_y );


        if(proc_pos_in_y == 0){
           Aee_y[comp]->re_initialize(
                nb_procs_in_y,
                nb_procs_in_y );
           schlur_complement_y[comp]->re_initialize(
                nb_procs_in_y,
                nb_procs_in_y );
           schlur_complement_y_ref[comp]->re_initialize(
                nb_procs_in_y,
                nb_procs_in_y );
        }

       VEC_rhs_velocity_1D_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 ) );

       VEC_interface_temp_y[comp]->re_initialize( nb_procs_in_y) ;
       VEC_temp_y[comp]->re_initialize( nb_procs_in_y) ;
   }
       

   if(dim>2)
   {
      if(U_is_zperiodic!=1 || nb_procs_in_z == 1)
      {
        if(proc_pos_in_z == nb_procs_in_z-1){
             Aii_z[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 ),
             nb_unknowns_handled_by_proc( 2 ));

             Aii_z_main_diagonal[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 ));
             Aii_z_super_diagonal[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 )-1);
             Aii_z_sub_diagonal[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 )-1);
             Aii_z_mod_super_diagonal[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 )-1);

             U_vec_zu[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 ));
             U_vec_zv[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 ));

             Aii_z_ref[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 ),
             nb_unknowns_handled_by_proc( 2 ));

             Aie_z[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 ),
             nb_procs_in_z-1);
             Aei_z[comp]->re_initialize(
             nb_procs_in_z-1,
             nb_unknowns_handled_by_proc( 2 ) );
             product_result_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 ) );

             VEC_local_temp_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )) ;
             VEC_local_solution_temp_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )) ;

          }
          else{
             Aii_z[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 )-1,
             nb_unknowns_handled_by_proc( 2 )-1 );

             Aii_z_main_diagonal[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 )-1);
             Aii_z_super_diagonal[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 )-2);
             Aii_z_sub_diagonal[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 )-2);
             Aii_z_mod_super_diagonal[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 )-2);

             Aii_z_ref[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 )-1,
             nb_unknowns_handled_by_proc( 2 )-1 );

             Aie_z[comp]->re_initialize(
             nb_unknowns_handled_by_proc( 2 )-1,
             nb_procs_in_z-1 );
             Aei_z[comp]->re_initialize(
             nb_procs_in_z-1,
             nb_unknowns_handled_by_proc( 2 )-1 );
             product_result_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1 );
             VEC_local_temp_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
             VEC_local_solution_temp_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
          }

          if(proc_pos_in_z == 0){
              Aee_z[comp]->re_initialize(
                   nb_procs_in_z-1,
                   nb_procs_in_z-1 );
              schlur_complement_z[comp]->re_initialize(
                   nb_procs_in_z-1,
                   nb_procs_in_z-1 );
              schlur_complement_z_ref[comp]->re_initialize(
                   nb_procs_in_z-1,
                   nb_procs_in_z-1 );
          }

          Aii_Aie_product_z[comp]->re_initialize(
             nb_procs_in_z-1);

          Aei_Aii_Aie_product_z[comp]->re_initialize(
             nb_procs_in_z-1,
             nb_procs_in_z-1 );

          VEC_interface_temp_z[comp]->re_initialize( nb_procs_in_z-1 ) ;
          VEC_temp_z[comp]->re_initialize( nb_procs_in_z-1 ) ;

          VEC_rhs_velocity_1D_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 ) );
      }
      else{
           Aii_z[comp]->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-1,
           nb_unknowns_handled_by_proc( 2 )-1 );

           Aii_z_main_diagonal[comp]->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-1);
           Aii_z_super_diagonal[comp]->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-2);
           Aii_z_sub_diagonal[comp]->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-2);
           Aii_z_mod_super_diagonal[comp]->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-2);

           Aii_z_ref[comp]->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-1,
           nb_unknowns_handled_by_proc( 2 )-1 );

           Aie_z[comp]->re_initialize(
           nb_unknowns_handled_by_proc( 2 )-1,
           nb_procs_in_z);
           Aei_z[comp]->re_initialize(
           nb_procs_in_z,
           nb_unknowns_handled_by_proc( 2 )-1 );
           product_result_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1 );
           VEC_local_temp_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
           VEC_local_solution_temp_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;

           if(proc_pos_in_z == 0){
              Aee_z[comp]->re_initialize(
                   nb_procs_in_z,
                   nb_procs_in_z );
              schlur_complement_z[comp]->re_initialize(
                   nb_procs_in_z,
                   nb_procs_in_z );
              schlur_complement_z_ref[comp]->re_initialize(
                   nb_procs_in_z,
                   nb_procs_in_z );
           }

           Aii_Aie_product_z[comp]->re_initialize(
             nb_procs_in_z);

           Aei_Aii_Aie_product_z[comp]->re_initialize(
             nb_procs_in_z,
             nb_procs_in_z );

           VEC_interface_temp_z[comp]->re_initialize( nb_procs_in_z ) ;
           VEC_temp_z[comp]->re_initialize( nb_procs_in_z ) ;

           VEC_rhs_velocity_1D_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 ) );
      }
          
    }

  } // Closing of for loop over components

}




//----------------------------------------------------------------------
DDS_NavierStokesSystem:: ~DDS_NavierStokesSystem( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: ~DDS_NavierStokesSystem" ) ;
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::initialize_velocity( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: initialize_velocity" ) ;

   UF->extract_unknown_DOFs_value( 0, UF_LOC ) ;
   UF_NUM->scatter()->set( UF_LOC, VEC_UF ) ;
         
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::initialize_DS_velocity( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: initialize_velocity" ) ;

   UF->extract_unknown_DOFs_value( 0, UF_DS_LOC ) ;
   UF_NUM->scatter()->set( UF_DS_LOC, VEC_DS_UF ) ;
         
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::initialize_DS_pressure( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: initialize_pressure" ) ;
   
   PF->extract_unknown_DOFs_value( 0, PF_DS_LOC ) ;
   PF_NUM->scatter()->set( PF_DS_LOC, VEC_DS_PF ) ;
         
}




//----------------------------------------------------------------------
LA_SeqVector const*
DDS_NavierStokesSystem:: get_solution_velocity( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_solution_velocity" ) ;

   UF_NUM->scatter()->get( VEC_UF, UF_LOC ) ;

   LA_SeqVector const* result = UF_LOC ;

   return( result ) ;

}




//----------------------------------------------------------------------
LA_SeqVector const*
DDS_NavierStokesSystem:: get_solution_DS_velocity( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_solution_DS_velocity" ) ;

   UF_NUM->scatter()->get( VEC_DS_UF, UF_DS_LOC ) ;

   LA_SeqVector const* DS_result = UF_DS_LOC ;

   return( DS_result ) ;

}




//----------------------------------------------------------------------
LA_SeqVector const*
DDS_NavierStokesSystem:: get_solution_DS_velocity_P( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_solution_DS_velocity_P" ) ;

   PF_NUM->scatter()->get( VEC_DS_PF, PF_DS_LOC ) ;

   LA_SeqVector const* DS_result = PF_DS_LOC ;

   return( DS_result ) ;

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem:: finalize_constant_matrices( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: finalize_constant_matrices" ) ;

   bool same_pattern = false ;
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();

   if ( macCOMM->rank() == 0 )
     cout << "Unsteady + laplacian matrix" << endl;

   MAT_D_velocityUnsteadyPlusDiffusion->add_Mat(
      	MAT_A_velocityUnsteady, 1., same_pattern ) ;

   SOLVER_velocity->set_matrix( MAT_D_velocityUnsteadyPlusDiffusion );

}




//----------------------------------------------------------------------
bool
DDS_NavierStokesSystem:: NavierStokes_solver( void  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: NavierStokes_solver" ) ;

   // Synchronize matrices & vectors
   VEC_UF->synchronize() ;

   // Compute velocity unsteady rhs
   MAT_A_velocityUnsteady->multiply_vec_then_add( VEC_UF,
   	VEC_rhs_A_velocityUnsteady ) ;

   // Add velocity BC rhs
   VEC_rhs_A_velocityUnsteady->sum(
   	VEC_rhs_D_velocityDiffusionPlusBodyTerm ) ;

   // Solve unsteady laplacian problem
   SOLVER_velocity->solve( VEC_rhs_A_velocityUnsteady,
   	VEC_UF );
   MAC_ASSERT( SOLVER_velocity->solution_is_achieved() ) ;

   return( true ) ;

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::at_each_time_step( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: at_each_time_step" ) ;

   // Store velocity at previous time
   VEC_UF->synchronize() ;
   VEC_UF_previoustime->set( VEC_UF ) ;
   VEC_DS_UF->synchronize() ;
   VEC_DS_UF_previoustime->set( VEC_DS_UF ) ;
   VEC_DS_PF->synchronize() ;

}




//----------------------------------------------------------------------
double
DDS_NavierStokesSystem:: compute_velocity_change( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: compute_velocity_change" ) ;

   VEC_UF->synchronize() ;
   VEC_UF_timechange->set( VEC_UF ) ;
   VEC_UF_timechange->sum( VEC_UF_previoustime, -1.0 ) ;

   double norm_UF = VEC_UF->two_norm() ;
   double time_change = VEC_UF_timechange->two_norm() ;
   if ( norm_UF > 1e-4 ) time_change /= norm_UF;

   return ( time_change ) ;

}




//----------------------------------------------------------------------
double
DDS_NavierStokesSystem:: compute_directionsplitting_velocity_change( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: "
   	"compute_directionsplitting_velocity_change" ) ;

   VEC_DS_UF->synchronize() ;
   VEC_DS_UF_timechange->set( VEC_DS_UF ) ;
   VEC_DS_UF_timechange->sum( VEC_DS_UF_previoustime, -1.0 ) ;

   double norm_UF = VEC_DS_UF->two_norm() ;
   double time_change = VEC_DS_UF_timechange->two_norm() ;
   // if(proc_pos_in_x == 0 && proc_pos_in_y == 0)
   //    MAC::out()<<"Time change "<<time_change<<endl;
   if ( norm_UF > 1e-4 ) time_change /= norm_UF;

   return ( time_change ) ;

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem:: assemble_velocity_unsteady_matrix(
	double const& coef )
//----------------------------------------------------------------------
{
   MAC_LABEL(
   "REG_ProjectionNavierStokesSystem:: assemble_velocity_unsteady_matrix" ) ;

   UF->assemble_mass_matrix( coef, MAT_A_velocityUnsteady );
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem:: assemble_velocity_diffusion_matrix_rhs(
	double const& coef_lap )
//----------------------------------------------------------------------
{
   MAC_LABEL(
   "DDS_NavierStokesSystem:: assemble_velocity_diffusion_matrix_rhs" ) ;

   UF->assemble_constantcoef_laplacian_matrix( coef_lap,
	MAT_D_velocityUnsteadyPlusDiffusion,
	VEC_rhs_D_velocityDiffusionPlusBodyTerm );

}




//----------------------------------------------------------------------
LA_Vector*
DDS_NavierStokesSystem::get_diffrhs_plus_bodyterm_vector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_diffrhs_plus_bodyterm_vector" ) ;

   return ( VEC_rhs_D_velocityDiffusionPlusBodyTerm ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_temp_x( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_x" ) ;

   return ( VEC_local_temp_x[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_solution_temp_x( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_x" ) ;

   return ( VEC_local_solution_temp_x[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_interface_temp_x( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_x" ) ;

   return ( VEC_interface_temp_x[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_temp_x( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_x" ) ;

   return ( VEC_temp_x[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_temp_y( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_y" ) ;

   return ( VEC_local_temp_y[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_solution_temp_y( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_y" ) ;

   return ( VEC_local_solution_temp_y[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_interface_temp_y( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_y" ) ;

   return ( VEC_interface_temp_y[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_temp_y( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_y" ) ;

   return ( VEC_temp_y[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_temp_z( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_z" ) ;

   return ( VEC_local_temp_z[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_solution_temp_z( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_z" ) ;

   return ( VEC_local_solution_temp_z[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_interface_temp_z( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_z" ) ;

   return ( VEC_interface_temp_z[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_temp_z( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_z" ) ;

   return ( VEC_temp_z[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_U_vec_xu( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_U_vec_xu" ) ;

   return ( U_vec_xu[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_U_vec_xv( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_U_vec_xv" ) ;

   return ( U_vec_xv[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_U_vec_yu( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_U_vec_yu" ) ;

   return ( U_vec_yu[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_U_vec_yv( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_U_vec_yv" ) ;

   return ( U_vec_yv[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_U_vec_zu( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_U_vec_zu" ) ;

   return ( U_vec_zu[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_U_vec_zv( size_t const& c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_U_vec_zv" ) ;

   return ( U_vec_zv[c] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_P_vec_xu( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_P_vec_xu" ) ;

   return ( P_vec_xu ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_P_vec_xv( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_P_vec_xv" ) ;

   return ( P_vec_xv ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_P_vec_yu( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_P_vec_yu" ) ;

   return ( P_vec_yu ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_P_vec_yv( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_P_vec_yv" ) ;

   return ( P_vec_yv ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_P_vec_zu( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_P_vec_zu" ) ;

   return ( P_vec_zu ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_P_vec_zv( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_P_vec_zv" ) ;

   return ( P_vec_zv ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_temp_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_x_P" ) ;

   return ( VEC_local_temp_x_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_solution_temp_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_x_P" ) ;

   return ( VEC_local_solution_temp_x_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_interface_temp_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_x_P" ) ;

   return ( VEC_interface_temp_x_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_temp_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_x_P" ) ;

   return ( VEC_temp_x_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_temp_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_y_P" ) ;

   return ( VEC_local_temp_y_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_solution_temp_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_y_P" ) ;

   return ( VEC_local_solution_temp_y_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_interface_temp_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_y_P" ) ;

   return ( VEC_interface_temp_y_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_temp_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_temp_y_P" ) ;

   return ( VEC_temp_y_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_temp_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_z_P" ) ;

   return ( VEC_local_temp_z_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_solution_temp_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_z_P" ) ;

   return ( VEC_local_solution_temp_z_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_interface_temp_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_z_P" ) ;

   return ( VEC_interface_temp_z_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_temp_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_z_P" ) ;

   return ( VEC_temp_z_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_rhs_DS_velocity_x( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_rhs_DS_velocity_x" ) ;

   return ( VEC_rhs_velocity_1D_x[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_rhs_DS_velocity_y( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_rhs_DS_velocity_y" ) ;

   return ( VEC_rhs_velocity_1D_y[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_rhs_DS_velocity_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_rhs_DS_velocity_z" ) ;

   return ( VEC_rhs_velocity_1D_z[comp] ) ;

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_Aii_x_ref( size_t const& comp)
//----------------------------------------------------------------------
{

   size_t nrows = Aii_x_main_diagonal[comp] -> nb_rows() ;

   double temp = Aii_x_main_diagonal[comp]->item(0);
   Aii_x_mod_super_diagonal[comp]-> set_item(0,Aii_x_super_diagonal[comp] ->item(0)/temp);

   //  // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=Aii_x_sub_diagonal[comp] ->item(m-1);
     b=Aii_x_main_diagonal[comp] ->item(m);
     prevc=Aii_x_mod_super_diagonal[comp] ->item(m-1);

     if(m<nrows-1){
         c=Aii_x_super_diagonal[comp] ->item(m);
         Aii_x_mod_super_diagonal[comp] ->set_item(m,c/(b - a*prevc));
     }
   }
   if(nb_procs_in_x == 1 && U_is_xperiodic == 1)
      DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], U_vec_xu[comp]);   
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_schlur_x_ref( size_t const& comp)
//----------------------------------------------------------------------
{
   size_t nrows = schlur_complement_x[comp] -> nb_rows() ;
   double temp = schlur_complement_x[comp]->item(0,0);
   schlur_complement_x_ref[comp]-> set_item(0,0,schlur_complement_x[comp] ->item(0,0));
   schlur_complement_x_ref[comp]-> set_item(0,1,schlur_complement_x[comp] ->item(0,1)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=schlur_complement_x[comp] ->item(m,m-1);
     b=schlur_complement_x[comp] ->item(m,m);
     c=schlur_complement_x[comp] ->item(m,m+1);
     prevc=schlur_complement_x_ref[comp] ->item(m-1,m);

     if(m<nrows-1)
         schlur_complement_x_ref[comp] ->set_item(m,m+1,c/(b - a*prevc));
      schlur_complement_x_ref[comp] ->set_item(m,m,b);
      schlur_complement_x_ref[comp] ->set_item(m,m-1,a);
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_Aii_y_ref( size_t const& comp)
//----------------------------------------------------------------------
{
   size_t nrows = Aii_y_main_diagonal[comp] -> nb_rows() ;
   double temp = Aii_y_main_diagonal[comp]->item(0);
   Aii_y_mod_super_diagonal[comp]-> set_item(0,Aii_y_super_diagonal[comp] ->item(0)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=Aii_y_sub_diagonal[comp] ->item(m-1);
     b=Aii_y_main_diagonal[comp] ->item(m);

     prevc=Aii_y_mod_super_diagonal[comp] ->item(m-1);

     if(m<nrows-1){
         c=Aii_y_super_diagonal[comp] ->item(m);
         Aii_y_mod_super_diagonal[comp] ->set_item(m,c/(b - a*prevc));
     }

   }
   if(nb_procs_in_y == 1 && U_is_yperiodic == 1)
      DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], U_vec_yu[comp]);   
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_schlur_y_ref( size_t const& comp)
//----------------------------------------------------------------------
{
   size_t nrows = schlur_complement_y[comp] -> nb_rows() ;
   double temp = schlur_complement_y[comp]->item(0,0);
   schlur_complement_y_ref[comp]-> set_item(0,0,schlur_complement_y[comp] ->item(0,0));
   schlur_complement_y_ref[comp]-> set_item(0,1,schlur_complement_y[comp] ->item(0,1)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=schlur_complement_y[comp] ->item(m,m-1);
     b=schlur_complement_y[comp] ->item(m,m);
     c=schlur_complement_y[comp] ->item(m,m+1);
     prevc=schlur_complement_y_ref[comp] ->item(m-1,m);

     if(m<nrows-1)
         schlur_complement_y_ref[comp] ->set_item(m,m+1,c/(b - a*prevc));
      schlur_complement_y_ref[comp] ->set_item(m,m,b);
      schlur_complement_y_ref[comp] ->set_item(m,m-1,a);
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_Aii_z_ref( size_t const& comp)
//----------------------------------------------------------------------
{
   size_t nrows = Aii_z_main_diagonal[comp] -> nb_rows() ;
   double temp = Aii_z_main_diagonal[comp]->item(0);
   Aii_z_mod_super_diagonal[comp]-> set_item(0,Aii_z_super_diagonal[comp] ->item(0)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=Aii_z_sub_diagonal[comp] ->item(m-1);
     b=Aii_z_main_diagonal[comp] ->item(m);

     prevc=Aii_z_mod_super_diagonal[comp] ->item(m-1);

     if(m<nrows-1){
         c=Aii_z_super_diagonal[comp] ->item(m);
         Aii_z_mod_super_diagonal[comp] ->set_item(m,c/(b - a*prevc));
     }
   }
   if(nb_procs_in_z == 1 && U_is_zperiodic == 1)
      DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], U_vec_zu[comp]);   
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_schlur_z_ref( size_t const& comp)
//----------------------------------------------------------------------
{
   size_t nrows = schlur_complement_z[comp] -> nb_rows() ;
   double temp = schlur_complement_z[comp]->item(0,0);
   schlur_complement_z_ref[comp]-> set_item(0,0,schlur_complement_z[comp] ->item(0,0));
   schlur_complement_z_ref[comp]-> set_item(0,1,schlur_complement_z[comp] ->item(0,1)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=schlur_complement_z[comp] ->item(m,m-1);
     b=schlur_complement_z[comp] ->item(m,m);
     c=schlur_complement_z[comp] ->item(m,m+1);
     prevc=schlur_complement_z_ref[comp] ->item(m-1,m);

     if(m<nrows-1)
         schlur_complement_z_ref[comp] ->set_item(m,m+1,c/(b - a*prevc));
      schlur_complement_z_ref[comp] ->set_item(m,m,b);
      schlur_complement_z_ref[comp] ->set_item(m,m-1,a);
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_Aii_x_ref_P( void )
//----------------------------------------------------------------------
{

   size_t nrows = Aii_x_main_diagonal_P -> nb_rows() ;

   double temp = Aii_x_main_diagonal_P->item(0);
   Aii_x_mod_super_diagonal_P-> set_item(0,Aii_x_super_diagonal_P ->item(0)/temp);

   //  // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=Aii_x_sub_diagonal_P ->item(m-1);
     b=Aii_x_main_diagonal_P ->item(m);
     prevc=Aii_x_mod_super_diagonal_P ->item(m-1);

     if(m<nrows-1){
         c=Aii_x_super_diagonal_P ->item(m);
         Aii_x_mod_super_diagonal_P ->set_item(m,c/(b - a*prevc));
     }
   }
   if(nb_procs_in_x == 1 && P_is_xperiodic == 1)
      DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, P_vec_xu);   

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_schlur_x_ref_P( void )
//----------------------------------------------------------------------
{
   size_t nrows = schlur_complement_x_P -> nb_rows() ;
   double temp = schlur_complement_x_P->item(0,0);
   schlur_complement_x_ref_P-> set_item(0,0,schlur_complement_x_P ->item(0,0));
   schlur_complement_x_ref_P-> set_item(0,1,schlur_complement_x_P ->item(0,1)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=schlur_complement_x_P ->item(m,m-1);
     b=schlur_complement_x_P ->item(m,m);
     c=schlur_complement_x_P ->item(m,m+1);
     prevc=schlur_complement_x_ref_P ->item(m-1,m);

     if(m<nrows-1)
         schlur_complement_x_ref_P ->set_item(m,m+1,c/(b - a*prevc));
      schlur_complement_x_ref_P ->set_item(m,m,b);
      schlur_complement_x_ref_P ->set_item(m,m-1,a);
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_Aii_y_ref_P( void )
//----------------------------------------------------------------------
{
   size_t nrows = Aii_y_main_diagonal_P -> nb_rows() ;
   double temp = Aii_y_main_diagonal_P->item(0);
   Aii_y_mod_super_diagonal_P-> set_item(0,Aii_y_super_diagonal_P ->item(0)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=Aii_y_sub_diagonal_P ->item(m-1);
     b=Aii_y_main_diagonal_P ->item(m);

     prevc=Aii_y_mod_super_diagonal_P ->item(m-1);

     if(m<nrows-1){
         c=Aii_y_super_diagonal_P ->item(m);
         Aii_y_mod_super_diagonal_P ->set_item(m,c/(b - a*prevc));
     }

   }
   if(nb_procs_in_y == 1 && P_is_yperiodic == 1)
      DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, P_vec_yu);   
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_schlur_y_ref_P( void )
//----------------------------------------------------------------------
{
   size_t nrows = schlur_complement_y_P -> nb_rows() ;
   double temp = schlur_complement_y_P->item(0,0);
   schlur_complement_y_ref_P-> set_item(0,0,schlur_complement_y_P ->item(0,0));
   schlur_complement_y_ref_P-> set_item(0,1,schlur_complement_y_P ->item(0,1)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=schlur_complement_y_P ->item(m,m-1);
     b=schlur_complement_y_P ->item(m,m);
     c=schlur_complement_y_P ->item(m,m+1);
     prevc=schlur_complement_y_ref_P ->item(m-1,m);

     if(m<nrows-1)
         schlur_complement_y_ref_P ->set_item(m,m+1,c/(b - a*prevc));
      schlur_complement_y_ref_P ->set_item(m,m,b);
      schlur_complement_y_ref_P ->set_item(m,m-1,a);
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_Aii_z_ref_P( void )
//----------------------------------------------------------------------
{
   size_t nrows = Aii_z_main_diagonal_P -> nb_rows() ;
   double temp = Aii_z_main_diagonal_P->item(0);
   Aii_z_mod_super_diagonal_P-> set_item(0,Aii_z_super_diagonal_P ->item(0)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=Aii_z_sub_diagonal_P ->item(m-1);
     b=Aii_z_main_diagonal_P ->item(m);

     prevc=Aii_z_mod_super_diagonal_P ->item(m-1);

     if(m<nrows-1){
         c=Aii_z_super_diagonal_P ->item(m);
         Aii_z_mod_super_diagonal_P ->set_item(m,c/(b - a*prevc));
     }
   }
   if(nb_procs_in_z == 1 && P_is_zperiodic == 1)
      DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, P_vec_zu);   

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_schlur_z_ref_P( void )
//----------------------------------------------------------------------
{
   size_t nrows = schlur_complement_z_P -> nb_rows() ;
   double temp = schlur_complement_z_P->item(0,0);
   schlur_complement_z_ref_P-> set_item(0,0,schlur_complement_z_P ->item(0,0));
   schlur_complement_z_ref_P-> set_item(0,1,schlur_complement_z_P ->item(0,1)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=schlur_complement_z_P ->item(m,m-1);
     b=schlur_complement_z_P ->item(m,m);
     c=schlur_complement_z_P ->item(m,m+1);
     prevc=schlur_complement_z_ref_P ->item(m-1,m);

     if(m<nrows-1)
         schlur_complement_z_ref_P ->set_item(m,m+1,c/(b - a*prevc));
      schlur_complement_z_ref_P ->set_item(m,m,b);
      schlur_complement_z_ref_P ->set_item(m,m-1,a);
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::thomas_algorithm( LA_SeqMatrix* mat_A, LA_SeqVector* rhs)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: thomas_algorithm" ) ;

   size_t nrows = mat_A -> nb_rows() ;
   double temp = mat_A->item(0,0);
   rhs-> set_item(0,rhs ->item(0)/temp);

    // Perform Forward Elimination
   size_t m;
   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,d,prevd,prevc;
     a=mat_A ->item(m,m-1);
     b=mat_A ->item(m,m);
     d=rhs->item(m);
     prevc=mat_A ->item(m-1,m);
     prevd=rhs->item(m-1);

     rhs -> set_item(m, (d-a*prevd)/(b-a*prevc));
   }

   //Perform backward substitution
   if(nrows>1){
      rhs->set_item(nrows-1,rhs->item(nrows-1));
      for (m = nrows-2; m< nrows-1;m--) {
         double c,nextd;
         c=mat_A ->item(m,m+1);
         nextd=rhs->item(m+1);
         rhs->add_to_item(m,-c*nextd);
      }
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::mod_thomas_algorithm( LA_SeqVector* x,LA_SeqVector* y,LA_SeqVector* z,LA_SeqVector* rhs)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: thomas_algorithm" ) ;

   size_t nrows = y -> nb_rows() ;
   double temp = y->item(0);
   rhs-> set_item(0,rhs ->item(0)/temp);

    // Perform Forward Elimination
   size_t m;
   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,d,prevd,prevc;
     a=x ->item(m-1);
     b=y ->item(m);
     d=rhs->item(m);
     prevc=z ->item(m-1);
     prevd=rhs->item(m-1);

     rhs -> set_item(m, (d-a*prevd)/(b-a*prevc));
   }

   //Perform backward substitution
   if(nrows>1){
      rhs->set_item(nrows-1,rhs->item(nrows-1));
      for (m = nrows-2; m< nrows-1;m--) {
         double c,nextd;
         c=z ->item(m);
         nextd=rhs->item(m+1);
         rhs->add_to_item(m,-c*nextd);
      }
   }
}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_main_diag_in_x( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x" ) ;

   return ( Aii_x_main_diagonal[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_super_diag_in_x( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x" ) ;

   return ( Aii_x_super_diagonal[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_sub_diag_in_x( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x" ) ;

   return ( Aii_x_sub_diagonal[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aie_in_x( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_x" ) ;

   return ( Aie_x[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aei_in_x( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_x" ) ;

   return ( Aei_x[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aee_matrix_in_x( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aee_matrix" ) ;

   return ( Aee_x[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_schlur_complement_in_x( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_x" ) ;

   return ( schlur_complement_x[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aei_Aii_Aie_product_in_x( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_x" ) ;

   return ( Aei_Aii_Aie_product_x[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_main_diag_in_y( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x" ) ;

   return ( Aii_y_main_diagonal[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_super_diag_in_y( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x" ) ;

   return ( Aii_y_super_diagonal[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_sub_diag_in_y( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x" ) ;

   return ( Aii_y_sub_diagonal[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aie_in_y( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_y" ) ;

   return ( Aie_y[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aei_in_y( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_y" ) ;

   return ( Aei_y[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aee_matrix_in_y( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aee_matrix_in_y" ) ;

   return ( Aee_y[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_schlur_complement_in_y( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_y" ) ;

   return ( schlur_complement_y[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aei_Aii_Aie_product_in_y( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_y" ) ;

   return ( Aei_Aii_Aie_product_y[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aii_in_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_z" ) ;

   return ( Aii_z[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_main_diag_in_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_z" ) ;

   return ( Aii_z_main_diagonal[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_super_diag_in_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x" ) ;

   return ( Aii_z_super_diagonal[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_sub_diag_in_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x" ) ;

   return ( Aii_z_sub_diagonal[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aie_in_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_z" ) ;

   return ( Aie_z[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aei_in_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_z" ) ;

   return ( Aei_z[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aee_matrix_in_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aee_matrix_in_z" ) ;

   return ( Aee_z[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_schlur_complement_in_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_z" ) ;

   return ( schlur_complement_z[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aei_Aii_Aie_product_in_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_z" ) ;

   return ( Aei_Aii_Aie_product_z[comp] ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_main_diag_in_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x_P" ) ;

   return ( Aii_x_main_diagonal_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_super_diag_in_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x_P" ) ;

   return ( Aii_x_super_diagonal_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_sub_diag_in_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_x_P" ) ;

   return ( Aii_x_sub_diagonal_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aie_in_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_x_P" ) ;

   return ( Aie_x_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aei_in_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_x_P" ) ;

   return ( Aei_x_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aee_matrix_in_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aee_matrix_P" ) ;

   return ( Aee_x_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_schlur_complement_in_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_x_P" ) ;

   return ( schlur_complement_x_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aei_Aii_Aie_product_in_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_x_P" ) ;

   return ( Aei_Aii_Aie_product_x_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_main_diag_in_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_y_P" ) ;

   return ( Aii_y_main_diagonal_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_super_diag_in_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_y_P" ) ;

   return ( Aii_y_super_diagonal_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_sub_diag_in_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_y_P" ) ;

   return ( Aii_y_sub_diagonal_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aie_in_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_y_P" ) ;

   return ( Aie_y_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aei_in_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_y_P" ) ;

   return ( Aei_y_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aee_matrix_in_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aee_matrix_P" ) ;

   return ( Aee_y_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_schlur_complement_in_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_y_P" ) ;

   return ( schlur_complement_y_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aei_Aii_Aie_product_in_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_y_P" ) ;

   return ( Aei_Aii_Aie_product_y_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_main_diag_in_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_z_P" ) ;

   return ( Aii_z_main_diagonal_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_super_diag_in_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_z_P" ) ;

   return ( Aii_z_super_diagonal_P ) ;

}




//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_sub_diag_in_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_in_z_P" ) ;

   return ( Aii_z_sub_diagonal_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aie_in_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_z_P" ) ;

   return ( Aie_z_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aei_in_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_z_P" ) ;

   return ( Aei_z_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aee_matrix_in_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aee_matrix_P" ) ;

   return ( Aee_z_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_schlur_complement_in_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_z_P" ) ;

   return ( schlur_complement_z_P ) ;

}




//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aei_Aii_Aie_product_in_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_z_P" ) ;

   return ( Aei_Aii_Aie_product_z_P ) ;

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_x_local_unknown_solver( LA_SeqVector* rhs, size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_x_local_unknown_solver" ) ;

   // Solve the DS splitting problem in
   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], rhs);

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_x_interface_unknown_solver( LA_SeqVector* rhs, size_t const& comp  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_x_interface_unknown_solver" ) ;

   // Solve the DS splitting problem in
   DDS_NavierStokesSystem::thomas_algorithm( schlur_complement_x_ref[comp], rhs);

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_y_local_unknown_solver( LA_SeqVector* rhs, size_t const& comp  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_y_local_unknown_solver" ) ;

   // Solve the DS splitting problem in
   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], rhs);
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_y_interface_unknown_solver( LA_SeqVector* rhs, size_t const& comp  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_y_interface_unknown_solver" ) ;

   // Solve the DS splitting problem in

   DDS_NavierStokesSystem::thomas_algorithm( schlur_complement_y_ref[comp], rhs);
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_z_local_unknown_solver( LA_SeqVector* rhs, size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_z_local_unknown_solver" ) ;

   // Solve the DS splitting problem in
	 DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], rhs);
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_z_interface_unknown_solver( LA_SeqVector* rhs, size_t const& comp)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_z_interface_unknown_solver" ) ;

   // Solve the DS splitting problem in z

   DDS_NavierStokesSystem::thomas_algorithm( schlur_complement_z_ref[comp], rhs);
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_x_local_unknown_solver_P( LA_SeqVector* rhs )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_x_local_unknown_solver_P" ) ;

   // Solve the DS splitting problem in x
   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, rhs);

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_x_interface_unknown_solver_P( LA_SeqVector* rhs  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_x_interface_unknown_solver_P" ) ;

   // Solve the DS splitting problem in x
   DDS_NavierStokesSystem::thomas_algorithm( schlur_complement_x_ref_P, rhs);

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_y_local_unknown_solver_P( LA_SeqVector* rhs )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_y_local_unknown_solver_P" ) ;

   // Solve the DS splitting problem in x
   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, rhs);

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_y_interface_unknown_solver_P( LA_SeqVector* rhs  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_y_interface_unknown_solver_P" ) ;

   // Solve the DS splitting problem in x
   DDS_NavierStokesSystem::thomas_algorithm( schlur_complement_y_ref_P, rhs);

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_z_local_unknown_solver_P( LA_SeqVector* rhs )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_z_local_unknown_solver_P" ) ;

   // Solve the DS splitting problem in x
   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, rhs);

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_z_interface_unknown_solver_P( LA_SeqVector* rhs  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_z_interface_unknown_solver_P" ) ;

   // Solve the DS splitting problem in x
   DDS_NavierStokesSystem::thomas_algorithm( schlur_complement_z_ref_P, rhs);

}




//----------------------------------------------------------------------
double
DDS_NavierStokesSystem::compute_vector_transpose_product(
  LA_SeqVector* a, LA_SeqVector* b)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquationSystem:: DS_HeatEquation_x_solver" ) ;
   size_t i;
   double result = 0.;
   for(i=0;i<a->nb_rows();i++)
      result += a->item(i)*b->item(i);

   return result;
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_x_solver(
	size_t const& j, size_t const& k, size_t const& min_i, LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_x_solver" ) ;

   // Solve the DS splitting problem in

   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], rhs);

   // Transfer in the distributed vector
   size_t nb_local_unk_x = rhs->nb_rows();
   size_t m, i, global_number_in_distributed_vector;

   for (m=0;m<nb_local_unk_x;++m)
   {
     i = min_i + m ;
     global_number_in_distributed_vector =
     	UF->DOF_global_number( i, j, k, comp );
     VEC_DS_UF->set_item( global_number_in_distributed_vector,
     	rhs->item( m ) );
   }

   // Put the interface unknowns in distributed vector
   if(nb_procs_in_x > 1){
      if(U_is_xperiodic == 1){
          i = min_i + nb_local_unk_x ;
          global_number_in_distributed_vector =
          UF->DOF_global_number( i, j, k, comp );
          VEC_DS_UF->set_item( global_number_in_distributed_vector,
             interface_rhs->item( proc_pos_in_x ) );
      }
      else if(proc_pos_in_x != nb_procs_in_x-1){
          i = min_i + nb_local_unk_x ;
          global_number_in_distributed_vector =
          UF->DOF_global_number( i, j, k, comp );
          VEC_DS_UF->set_item( global_number_in_distributed_vector,
             interface_rhs->item( proc_pos_in_x ) );
      }
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_x_solver_periodic(
  size_t const& j, size_t const& k, size_t const& min_i, LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_x_solver" ) ;

   // Solve the DS splitting problem in

   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], rhs);

   double scale = (compute_vector_transpose_product(U_vec_xv[comp],rhs))/(1+compute_vector_transpose_product(U_vec_xv[comp],U_vec_xu[comp]) );

   rhs->sum(U_vec_xu[comp],-1.0*scale);

   // Transfer in the distributed vector
   size_t nb_local_unk_x = rhs->nb_rows();
   size_t m, i, global_number_in_distributed_vector;

   for (m=0;m<nb_local_unk_x;++m)
   {
     i = min_i + m ;
     global_number_in_distributed_vector =
      UF->DOF_global_number( i, j, k, comp );
     VEC_DS_UF->set_item( global_number_in_distributed_vector,
      rhs->item( m ) );
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_y_solver(
   size_t const& i, size_t const& k, size_t const& min_j ,LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_y_solver" ) ;

   // Solve the DS splitting problem in
   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], rhs);

   // Transfer in the distributed vector
   size_t nb_local_unk_y = rhs->nb_rows();
   size_t m, j, global_number_in_distributed_vector;
   for (m=0;m<nb_local_unk_y;++m)
   {
     j = min_j + m ;
     global_number_in_distributed_vector =
      UF->DOF_global_number( i, j, k, comp );
     VEC_DS_UF->set_item( global_number_in_distributed_vector,
      rhs->item( m ) );
   }

   // Put the interface unknowns in distributed vector
   if(nb_procs_in_y > 1){
      if(U_is_yperiodic == 1){
          j = min_j + nb_local_unk_y ;
          global_number_in_distributed_vector =
          UF->DOF_global_number( i, j, k, comp);
          VEC_DS_UF->set_item( global_number_in_distributed_vector,
             interface_rhs->item( proc_pos_in_y ) );
      }
      else if(proc_pos_in_y != nb_procs_in_y-1){
          j = min_j + nb_local_unk_y ;
          global_number_in_distributed_vector =
          UF->DOF_global_number( i, j, k, comp);
          VEC_DS_UF->set_item( global_number_in_distributed_vector,
             interface_rhs->item( proc_pos_in_y ) );
      }
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_y_solver_periodic(
   size_t const& i, size_t const& k, size_t const& min_j ,LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_y_solver" ) ;

   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], rhs);

   double scale = (compute_vector_transpose_product(U_vec_yv[comp],rhs))/(1+compute_vector_transpose_product(U_vec_yv[comp],U_vec_yu[comp]) );

   rhs->sum(U_vec_yu[comp],-1.0*scale);

   // Transfer in the distributed vector
   size_t nb_local_unk_y = rhs->nb_rows();
   size_t m, j, global_number_in_distributed_vector;
   for (m=0;m<nb_local_unk_y;++m)
   {
     j = min_j + m ;
     global_number_in_distributed_vector =
      UF->DOF_global_number( i, j, k, comp );
     VEC_DS_UF->set_item( global_number_in_distributed_vector,
      rhs->item( m ) );
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_z_solver(
	size_t const& i, size_t const& j, size_t const& min_k ,LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_z_solver" ) ;

   // Solve the DS splitting problem in

	 DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], rhs);

   // Transfer in the distributed vector
   size_t nb_local_unk_z = rhs->nb_rows();
   size_t m, k, global_number_in_distributed_vector;
   for (m=0;m<nb_local_unk_z;++m)
   {
     k = min_k + m ;
     global_number_in_distributed_vector =
     	UF->DOF_global_number( i, j, k, comp );
     VEC_DS_UF->set_item( global_number_in_distributed_vector,
     	rhs->item( m ) );
   }

    // Put the interface unknowns in distributed vector
   if(nb_procs_in_z > 1){
      if(U_is_zperiodic == 1){
          k = min_k + nb_local_unk_z ;
          global_number_in_distributed_vector =
          UF->DOF_global_number( i, j, k, comp );
          VEC_DS_UF->set_item( global_number_in_distributed_vector,
             interface_rhs->item( proc_pos_in_z ) );
      }
      else if(proc_pos_in_z != nb_procs_in_z-1){
          k = min_k + nb_local_unk_z ;
          global_number_in_distributed_vector =
          UF->DOF_global_number( i, j, k, comp );
          VEC_DS_UF->set_item( global_number_in_distributed_vector,
             interface_rhs->item( proc_pos_in_z ) );
      } 
   }
   
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_z_solver_periodic(
  size_t const& i, size_t const& j, size_t const& min_k ,LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_z_solver" ) ;

   // Solve the DS splitting problem in

   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], rhs);

   double scale = (compute_vector_transpose_product(U_vec_zv[comp],rhs))/(1+compute_vector_transpose_product(U_vec_zv[comp],U_vec_zu[comp]) );

   rhs->sum(U_vec_zu[comp],-1.0*scale);

   // Transfer in the distributed vector
   size_t nb_local_unk_z = rhs->nb_rows();
   size_t m, k, global_number_in_distributed_vector;
   for (m=0;m<nb_local_unk_z;++m)
   {
     k = min_k + m ;
     global_number_in_distributed_vector =
      UF->DOF_global_number( i, j, k, comp );
     VEC_DS_UF->set_item( global_number_in_distributed_vector,
      rhs->item( m ) );
   }
   
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_x_solver_P(
  size_t const& j, size_t const& k, size_t const& min_i, LA_SeqVector* rhs, LA_SeqVector* interface_rhs )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_x_solver_P" ) ;

   // Solve the DS splitting problem in

   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, rhs);
   
   // Transfer in the distributed vector for pressure
   size_t nb_local_unk_x = rhs->nb_rows();
   size_t m, i, global_number_in_distributed_vector;

   for (m=0;m<nb_local_unk_x;++m)
   {
     i = min_i + m ;
     global_number_in_distributed_vector =
      PF->DOF_global_number( i, j, k, 0 );
     VEC_DS_PF->set_item( global_number_in_distributed_vector,
      rhs->item( m ) );
   }

   // Put the interface unknowns in distributed vector for pressure
   if(nb_procs_in_x > 1){
    if(P_is_xperiodic == 1){
      i = min_i + nb_local_unk_x ;
      global_number_in_distributed_vector =
      PF->DOF_global_number( i, j, k, 0 );
      VEC_DS_PF->set_item( global_number_in_distributed_vector,
         interface_rhs->item( proc_pos_in_x ) );
    }
    else if(proc_pos_in_x != nb_procs_in_x-1){
      i = min_i + nb_local_unk_x ;
      global_number_in_distributed_vector =
      PF->DOF_global_number( i, j, k, 0 );
      VEC_DS_PF->set_item( global_number_in_distributed_vector,
         interface_rhs->item( proc_pos_in_x ) );
    } 
   }
   
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_x_solver_P_periodic(
  size_t const& j, size_t const& k, size_t const& min_i, LA_SeqVector* rhs, LA_SeqVector* interface_rhs )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_x_solver_P" ) ;

   // Solve the DS splitting problem in

   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, rhs);

   double scale = (compute_vector_transpose_product(P_vec_xv,rhs))/(1+compute_vector_transpose_product(P_vec_xv,P_vec_xu) );

   rhs->sum(P_vec_xu,-1.0*scale);

   // Transfer in the distributed vector for pressure
   size_t nb_local_unk_x = rhs->nb_rows();
   size_t m, i, global_number_in_distributed_vector;

   for (m=0;m<nb_local_unk_x;++m)
   {
     i = min_i + m ;
     global_number_in_distributed_vector =
      PF->DOF_global_number( i, j, k, 0 );
     VEC_DS_PF->set_item( global_number_in_distributed_vector,
      rhs->item( m ) );
   }
   
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_y_solver_P(
   size_t const& i, size_t const& k, size_t const& min_j ,LA_SeqVector* rhs, LA_SeqVector* interface_rhs )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_y_solver_P" ) ;

   // Solve the DS splitting problem in
   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, rhs);

   /*if(i >= 8){
      rhs->print_items(MAC::out(),0); 
   }*/
   // Transfer in the distributed vector for pressure
   size_t nb_local_unk_y = rhs->nb_rows();
   size_t m, j, global_number_in_distributed_vector;
   for (m=0;m<nb_local_unk_y;++m)
   {
     j = min_j + m ;
     global_number_in_distributed_vector =
      PF->DOF_global_number( i, j, k, 0 );
     VEC_DS_PF->set_item( global_number_in_distributed_vector,
      rhs->item( m ) );
   }

   // Put the interface unknowns in distributed vector for pressure
   if(nb_procs_in_y > 1){
    if(P_is_yperiodic == 1){
      j = min_j + nb_local_unk_y ;
      global_number_in_distributed_vector =
      PF->DOF_global_number( i, j, k, 0);
      VEC_DS_PF->set_item( global_number_in_distributed_vector,
         interface_rhs->item( proc_pos_in_y ) );
    }
    else if(proc_pos_in_y != nb_procs_in_y-1){
      j = min_j + nb_local_unk_y ;
      global_number_in_distributed_vector =
      PF->DOF_global_number( i, j, k, 0);
      VEC_DS_PF->set_item( global_number_in_distributed_vector,
         interface_rhs->item( proc_pos_in_y ) );
    }
   }

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_y_solver_P_periodic(
   size_t const& i, size_t const& k, size_t const& min_j ,LA_SeqVector* rhs, LA_SeqVector* interface_rhs )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_y_solver_P" ) ;

   // Solve the DS splitting problem in
   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, rhs);

   double scale = (compute_vector_transpose_product(P_vec_yv,rhs))/(1+compute_vector_transpose_product(P_vec_yv,P_vec_yu) );

   rhs->sum(P_vec_yu,-1.0*scale);

   // Transfer in the distributed vector for pressure
   size_t nb_local_unk_y = rhs->nb_rows();
   size_t m, j, global_number_in_distributed_vector;
   for (m=0;m<nb_local_unk_y;++m)
   {
     j = min_j + m ;
     global_number_in_distributed_vector =
      PF->DOF_global_number( i, j, k, 0 );
     VEC_DS_PF->set_item( global_number_in_distributed_vector,
      rhs->item( m ) );
   }

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_z_solver_P(
  size_t const& i, size_t const& j, size_t const& min_k ,LA_SeqVector* rhs, LA_SeqVector* interface_rhs )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_z_solver_P" ) ;

   // Solve the DS splitting problem in

   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, rhs);

   // Transfer in the distributed vector for pressure
   size_t nb_local_unk_z = rhs->nb_rows();
   size_t m, k, global_number_in_distributed_vector;
   for (m=0;m<nb_local_unk_z;++m)
   {
     k = min_k + m ;
     global_number_in_distributed_vector =
      PF->DOF_global_number( i, j, k, 0 );
     VEC_DS_PF->set_item( global_number_in_distributed_vector,
      rhs->item( m ) );
   }

    // Put the interface unknowns in distributed vector for pressure
   if(nb_procs_in_z > 1){
      if(P_is_zperiodic == 1){
          k = min_k + nb_local_unk_z ;
          global_number_in_distributed_vector =
          PF->DOF_global_number( i, j, k, 0 );
          VEC_DS_PF->set_item( global_number_in_distributed_vector,
             interface_rhs->item( proc_pos_in_z ) );
      }
      else if(proc_pos_in_z != nb_procs_in_z-1){
          k = min_k + nb_local_unk_z ;
          global_number_in_distributed_vector =
          PF->DOF_global_number( i, j, k, 0 );
          VEC_DS_PF->set_item( global_number_in_distributed_vector,
             interface_rhs->item( proc_pos_in_z ) );
      }
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_z_solver_P_periodic(
  size_t const& i, size_t const& j, size_t const& min_k ,LA_SeqVector* rhs, LA_SeqVector* interface_rhs )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_z_solver_P" ) ;

   // Solve the DS splitting problem in

   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, rhs);

   double scale = (compute_vector_transpose_product(P_vec_zv,rhs))/(1+compute_vector_transpose_product(P_vec_zv,P_vec_zu) );

   rhs->sum(P_vec_zu,-1.0*scale);

   // Transfer in the distributed vector for pressure
   size_t nb_local_unk_z = rhs->nb_rows();
   size_t m, k, global_number_in_distributed_vector;
   for (m=0;m<nb_local_unk_z;++m)
   {
     k = min_k + m ;
     global_number_in_distributed_vector =
      PF->DOF_global_number( i, j, k, 0 );
     VEC_DS_PF->set_item( global_number_in_distributed_vector,
      rhs->item( m ) );
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::synchronize_DS_solution_vec( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: synchronize_DS_solution_vec" ) ;

   VEC_DS_UF->synchronize();
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::synchronize_solution_vec( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: synchronize_solution_vec" ) ;

   VEC_UF->synchronize();
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::synchronize_DS_solution_vec_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: synchronize_DS_solution_vec" ) ;

   VEC_DS_PF->synchronize();
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_product_matrix_x( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: compute_product_matrix_x" ) ;

   int i;
   if(U_is_xperiodic != 1)
   {
       if(proc_pos_in_x == nb_procs_in_x-1)
       {

          // Get appropriate column of Aie
          Aie_x[comp] -> extract_col( proc_pos_in_x-1, product_result_x[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], product_result_x[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_x[comp]->multiply_vec_then_add(product_result_x[comp],Aii_Aie_product_x[comp]);
          for(i=0;i<nb_procs_in_x-1;i++){
             Aei_Aii_Aie_product_x[comp]->set_item(i,proc_pos_in_x-1,Aii_Aie_product_x[comp]->item(i));
          }
       }

       else if(proc_pos_in_x == 0){

          // Get appropriate column of Aie
          Aie_x[comp] -> extract_col( proc_pos_in_x, product_result_x[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], product_result_x[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_x[comp]->multiply_vec_then_add(product_result_x[comp],Aii_Aie_product_x[comp]);
          for(i=0;i<nb_procs_in_x-1;i++){
          Aei_Aii_Aie_product_x[comp]->set_item(i,proc_pos_in_x,Aii_Aie_product_x[comp]->item(i));
          }
       }
       else{

          // Get appropriate column of Aie
          Aie_x[comp] -> extract_col( proc_pos_in_x-1, product_result_x[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], product_result_x[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_x[comp]->multiply_vec_then_add(product_result_x[comp],Aii_Aie_product_x[comp]);
          for(i=0;i<nb_procs_in_x-1;i++){
             Aei_Aii_Aie_product_x[comp]->set_item(i,proc_pos_in_x-1,Aii_Aie_product_x[comp]->item(i));
          }

          // Get appropriate column of Aie
          Aie_x[comp] -> extract_col( proc_pos_in_x, product_result_x[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], product_result_x[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_x[comp]->multiply_vec_then_add(product_result_x[comp],Aii_Aie_product_x[comp]);
          for(i=0;i<nb_procs_in_x-1;i++){
             Aei_Aii_Aie_product_x[comp]->set_item(i,proc_pos_in_x,Aii_Aie_product_x[comp]->item(i));
          }
       } 
   }
   else{
    
    if(proc_pos_in_x == 0){

      // Get appropriate column of Aie
      Aie_x[comp] -> extract_col( proc_pos_in_x, product_result_x[comp] );

      // Get inv(Aii)*Aie for for appropriate column of Aie
      DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], product_result_x[comp]);

      // Get product of Aei*inv(Aii)*Aie for appropriate column
      Aei_x[comp]->multiply_vec_then_add(product_result_x[comp],Aii_Aie_product_x[comp]);
      for(i=0;i<nb_procs_in_x;i++){
      Aei_Aii_Aie_product_x[comp]->set_item(i,proc_pos_in_x,Aii_Aie_product_x[comp]->item(i));
      }

      // Get appropriate column of Aie
      Aie_x[comp] -> extract_col( nb_procs_in_x-1, product_result_x[comp] );

      // Get inv(Aii)*Aie for for appropriate column of Aie
      DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], product_result_x[comp]);

      // Get product of Aei*inv(Aii)*Aie for appropriate column
      Aei_x[comp]->multiply_vec_then_add(product_result_x[comp],Aii_Aie_product_x[comp]);
      for(i=0;i<nb_procs_in_x;i++){
      Aei_Aii_Aie_product_x[comp]->set_item(i,nb_procs_in_x-1,Aii_Aie_product_x[comp]->item(i));
      }
    }
    else{

      // Get appropriate column of Aie
      Aie_x[comp] -> extract_col( proc_pos_in_x-1, product_result_x[comp] );

      // Get inv(Aii)*Aie for for appropriate column of Aie
      DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], product_result_x[comp]);

      // Get product of Aei*inv(Aii)*Aie for appropriate column
      Aei_x[comp]->multiply_vec_then_add(product_result_x[comp],Aii_Aie_product_x[comp]);
      for(i=0;i<nb_procs_in_x;i++){
         Aei_Aii_Aie_product_x[comp]->set_item(i,proc_pos_in_x-1,Aii_Aie_product_x[comp]->item(i));
      }

      // Get appropriate column of Aie
      Aie_x[comp] -> extract_col( proc_pos_in_x, product_result_x[comp] );

      // Get inv(Aii)*Aie for for appropriate column of Aie
      DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal[comp],Aii_x_main_diagonal[comp],Aii_x_mod_super_diagonal[comp], product_result_x[comp]);

      // Get product of Aei*inv(Aii)*Aie for appropriate column
      Aei_x[comp]->multiply_vec_then_add(product_result_x[comp],Aii_Aie_product_x[comp]);
      for(i=0;i<nb_procs_in_x;i++){
         Aei_Aii_Aie_product_x[comp]->set_item(i,proc_pos_in_x,Aii_Aie_product_x[comp]->item(i));
      }
   }
   }
   

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_product_matrix_y( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: compute_product_matrix_y" ) ;

   int i;
   if(U_is_yperiodic != 1){
      if(proc_pos_in_y == nb_procs_in_y-1){

          // Get appropriate column of Aie
          Aie_y[comp] -> extract_col( proc_pos_in_y-1, product_result_y[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], product_result_y[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_y[comp]->multiply_vec_then_add(product_result_y[comp],Aii_Aie_product_y[comp]);
          for(i=0;i<nb_procs_in_y-1;i++){
             Aei_Aii_Aie_product_y[comp]->set_item(i,proc_pos_in_y-1,Aii_Aie_product_y[comp]->item(i));
          }
       }
       else if(proc_pos_in_y==0){

          // Get appropriate column of Aie
          Aie_y[comp] -> extract_col( proc_pos_in_y, product_result_y[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], product_result_y[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_y[comp]->multiply_vec_then_add(product_result_y[comp],Aii_Aie_product_y[comp]);
          for(i=0;i<nb_procs_in_y-1;i++){
             Aei_Aii_Aie_product_y[comp]->set_item(i,proc_pos_in_y,Aii_Aie_product_y[comp]->item(i));
          }
       }
       else{

          // Get appropriate column of Aie
          Aie_y[comp] -> extract_col( proc_pos_in_y-1, product_result_y[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], product_result_y[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_y[comp]->multiply_vec_then_add(product_result_y[comp],Aii_Aie_product_y[comp]);
          for(i=0;i<nb_procs_in_y-1;i++){
             Aei_Aii_Aie_product_y[comp]->set_item(i,proc_pos_in_y-1,Aii_Aie_product_y[comp]->item(i));
          }

          // Get appropriate column of Aie
          Aie_y[comp] -> extract_col( proc_pos_in_y, product_result_y[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], product_result_y[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_y[comp]->multiply_vec_then_add(product_result_y[comp],Aii_Aie_product_y[comp]);
          for(i=0;i<nb_procs_in_y-1;i++){
             Aei_Aii_Aie_product_y[comp]->set_item(i,proc_pos_in_y,Aii_Aie_product_y[comp]->item(i));
          }
        }
   }
   else{
       if(proc_pos_in_y==0){

          // Get appropriate column of Aie
          Aie_y[comp] -> extract_col( proc_pos_in_y, product_result_y[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], product_result_y[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_y[comp]->multiply_vec_then_add(product_result_y[comp],Aii_Aie_product_y[comp]);
          for(i=0;i<nb_procs_in_y;i++){
             Aei_Aii_Aie_product_y[comp]->set_item(i,proc_pos_in_y,Aii_Aie_product_y[comp]->item(i));
          }

          // Get appropriate column of Aie
          Aie_y[comp] -> extract_col( nb_procs_in_y-1, product_result_y[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], product_result_y[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_y[comp]->multiply_vec_then_add(product_result_y[comp],Aii_Aie_product_y[comp]);
          for(i=0;i<nb_procs_in_y;i++){
             Aei_Aii_Aie_product_y[comp]->set_item(i,nb_procs_in_y-1,Aii_Aie_product_y[comp]->item(i));
          }
       }
       else{

          // Get appropriate column of Aie
          Aie_y[comp] -> extract_col( proc_pos_in_y-1, product_result_y[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], product_result_y[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_y[comp]->multiply_vec_then_add(product_result_y[comp],Aii_Aie_product_y[comp]);
          for(i=0;i<nb_procs_in_y;i++){
             Aei_Aii_Aie_product_y[comp]->set_item(i,proc_pos_in_y-1,Aii_Aie_product_y[comp]->item(i));
          }

          // Get appropriate column of Aie
          Aie_y[comp] -> extract_col( proc_pos_in_y, product_result_y[comp] );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal[comp],Aii_y_main_diagonal[comp],Aii_y_mod_super_diagonal[comp], product_result_y[comp]);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_y[comp]->multiply_vec_then_add(product_result_y[comp],Aii_Aie_product_y[comp]);
          for(i=0;i<nb_procs_in_y;i++){
             Aei_Aii_Aie_product_y[comp]->set_item(i,proc_pos_in_y,Aii_Aie_product_y[comp]->item(i));
          }
        }
   }       

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_product_matrix_z( size_t const& comp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: compute_product_matrix_z" ) ;

   int i;

   if(U_is_zperiodic !=1){
      if(proc_pos_in_z == nb_procs_in_z-1){

        // Get appropriate column of Aie
        Aie_z[comp] -> extract_col( proc_pos_in_z-1, product_result_z[comp] );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], product_result_z[comp]);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z[comp]->multiply_vec_then_add(product_result_z[comp],Aii_Aie_product_z[comp]);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z[comp]->set_item(i,proc_pos_in_z-1,Aii_Aie_product_z[comp]->item(i));
        }
     }
     else if(proc_pos_in_z==0){

        // Get appropriate column of Aie
        Aie_z[comp] -> extract_col( proc_pos_in_z, product_result_z[comp] );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], product_result_z[comp]);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z[comp]->multiply_vec_then_add(product_result_z[comp],Aii_Aie_product_z[comp]);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z[comp]->set_item(i,proc_pos_in_z,Aii_Aie_product_z[comp]->item(i));
        }
     }
     else{

        // Get appropriate column of Aie
        Aie_z[comp] -> extract_col( proc_pos_in_z-1, product_result_z[comp] );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], product_result_z[comp]);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z[comp]->multiply_vec_then_add(product_result_z[comp],Aii_Aie_product_z[comp]);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z[comp]->set_item(i,proc_pos_in_z-1,Aii_Aie_product_z[comp]->item(i));
        }

        // Get appropriate column of Aie
        Aie_z[comp] -> extract_col( proc_pos_in_z, product_result_z[comp] );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], product_result_z[comp]);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z[comp]->multiply_vec_then_add(product_result_z[comp],Aii_Aie_product_z[comp]);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z[comp]->set_item(i,proc_pos_in_z,Aii_Aie_product_z[comp]->item(i));
        }
      }
   }
   else{
     if(proc_pos_in_z==0){

        // Get appropriate column of Aie
        Aie_z[comp] -> extract_col( proc_pos_in_z, product_result_z[comp] );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], product_result_z[comp]);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z[comp]->multiply_vec_then_add(product_result_z[comp],Aii_Aie_product_z[comp]);
        for(i=0;i<nb_procs_in_z;i++){
           Aei_Aii_Aie_product_z[comp]->set_item(i,proc_pos_in_z,Aii_Aie_product_z[comp]->item(i));
        }

        // Get appropriate column of Aie
        Aie_z[comp] -> extract_col( nb_procs_in_z-1, product_result_z[comp] );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], product_result_z[comp]);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z[comp]->multiply_vec_then_add(product_result_z[comp],Aii_Aie_product_z[comp]);
        for(i=0;i<nb_procs_in_z;i++){
           Aei_Aii_Aie_product_z[comp]->set_item(i,nb_procs_in_z-1,Aii_Aie_product_z[comp]->item(i));
        }
     }
     else{

        // Get appropriate column of Aie
        Aie_z[comp] -> extract_col( proc_pos_in_z-1, product_result_z[comp] );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], product_result_z[comp]);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z[comp]->multiply_vec_then_add(product_result_z[comp],Aii_Aie_product_z[comp]);
        for(i=0;i<nb_procs_in_z;i++){
           Aei_Aii_Aie_product_z[comp]->set_item(i,proc_pos_in_z-1,Aii_Aie_product_z[comp]->item(i));
        }

        // Get appropriate column of Aie
        Aie_z[comp] -> extract_col( proc_pos_in_z, product_result_z[comp] );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal[comp],Aii_z_main_diagonal[comp],Aii_z_mod_super_diagonal[comp], product_result_z[comp]);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z[comp]->multiply_vec_then_add(product_result_z[comp],Aii_Aie_product_z[comp]);
        for(i=0;i<nb_procs_in_z;i++){
           Aei_Aii_Aie_product_z[comp]->set_item(i,proc_pos_in_z,Aii_Aie_product_z[comp]->item(i));
        }
      }
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_product_matrix_x_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: compute_product_matrix_x_P" ) ;

   int i;
   if(P_is_xperiodic != 1)
   {
       if(proc_pos_in_x == nb_procs_in_x-1){

          // Get appropriate column of Aie
          Aie_x_P -> extract_col( proc_pos_in_x-1, product_result_x_P );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, product_result_x_P);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_x_P->multiply_vec_then_add(product_result_x_P,Aii_Aie_product_x_P);
          for(i=0;i<nb_procs_in_x-1;i++){
             Aei_Aii_Aie_product_x_P->set_item(i,proc_pos_in_x-1,Aii_Aie_product_x_P->item(i));
          }
       }

       else if(proc_pos_in_x == 0){

          // Get appropriate column of Aie
          Aie_x_P -> extract_col( proc_pos_in_x, product_result_x_P );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, product_result_x_P);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_x_P->multiply_vec_then_add(product_result_x_P,Aii_Aie_product_x_P);
          for(i=0;i<nb_procs_in_x-1;i++){
          Aei_Aii_Aie_product_x_P->set_item(i,proc_pos_in_x,Aii_Aie_product_x_P->item(i));
        }
       }
       else{

          // Get appropriate column of Aie
          Aie_x_P -> extract_col( proc_pos_in_x-1, product_result_x_P );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, product_result_x_P);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_x_P->multiply_vec_then_add(product_result_x_P,Aii_Aie_product_x_P);
        for(i=0;i<nb_procs_in_x-1;i++){
           Aei_Aii_Aie_product_x_P->set_item(i,proc_pos_in_x-1,Aii_Aie_product_x_P->item(i));
        }

          // Get appropriate column of Aie
        Aie_x_P -> extract_col( proc_pos_in_x, product_result_x_P );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, product_result_x_P);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_x_P->multiply_vec_then_add(product_result_x_P,Aii_Aie_product_x_P);
        for(i=0;i<nb_procs_in_x-1;i++){
           Aei_Aii_Aie_product_x_P->set_item(i,proc_pos_in_x,Aii_Aie_product_x_P->item(i));
        }
       }
    }
    else
    {
      if(proc_pos_in_x == 0){

          // Get appropriate column of Aie
          Aie_x_P -> extract_col( proc_pos_in_x, product_result_x_P );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, product_result_x_P);

          // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_x_P->multiply_vec_then_add(product_result_x_P,Aii_Aie_product_x_P);
          for(i=0;i<nb_procs_in_x-1;i++){
            Aei_Aii_Aie_product_x_P->set_item(i,proc_pos_in_x,Aii_Aie_product_x_P->item(i));
          }

          // Get appropriate column of Aie
          Aie_x_P -> extract_col( nb_procs_in_x-1, product_result_x_P );

            // Get inv(Aii)*Aie for for appropriate column of Aie
            DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, product_result_x_P);

            // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_x_P->multiply_vec_then_add(product_result_x_P,Aii_Aie_product_x_P);
          for(i=0;i<nb_procs_in_x-1;i++){
             Aei_Aii_Aie_product_x_P->set_item(i,nb_procs_in_x-1,Aii_Aie_product_x_P->item(i));
          }
       }
       else{

          // Get appropriate column of Aie
          Aie_x_P -> extract_col( proc_pos_in_x-1, product_result_x_P );

          // Get inv(Aii)*Aie for for appropriate column of Aie
          DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, product_result_x_P);

            // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_x_P->multiply_vec_then_add(product_result_x_P,Aii_Aie_product_x_P);
          for(i=0;i<nb_procs_in_x-1;i++){
             Aei_Aii_Aie_product_x_P->set_item(i,proc_pos_in_x-1,Aii_Aie_product_x_P->item(i));
          }

          // Get appropriate column of Aie
          Aie_x_P -> extract_col( proc_pos_in_x, product_result_x_P );

            // Get inv(Aii)*Aie for for appropriate column of Aie
            DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_x_sub_diagonal_P,Aii_x_main_diagonal_P,Aii_x_mod_super_diagonal_P, product_result_x_P);

            // Get product of Aei*inv(Aii)*Aie for appropriate column
          Aei_x_P->multiply_vec_then_add(product_result_x_P,Aii_Aie_product_x_P);
          for(i=0;i<nb_procs_in_x-1;i++){
             Aei_Aii_Aie_product_x_P->set_item(i,proc_pos_in_x,Aii_Aie_product_x_P->item(i));
          }
       }
    }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_product_matrix_y_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: compute_product_matrix_y_P" ) ;

   int i;

   if(P_is_yperiodic != 1)
   {
    if(proc_pos_in_y == nb_procs_in_y-1){
        // Get appropriate column of Aie
        Aie_y_P -> extract_col( proc_pos_in_y-1, product_result_y_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, product_result_y_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_y_P->multiply_vec_then_add(product_result_y_P,Aii_Aie_product_y_P);
        for(i=0;i<nb_procs_in_y-1;i++){
           Aei_Aii_Aie_product_y_P->set_item(i,proc_pos_in_y-1,Aii_Aie_product_y_P->item(i));
        }
     }
     if(proc_pos_in_y==0){

        // Get appropriate column of Aie
        Aie_y_P -> extract_col( proc_pos_in_y, product_result_y_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, product_result_y_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_y_P->multiply_vec_then_add(product_result_y_P,Aii_Aie_product_y_P);
        for(i=0;i<nb_procs_in_y-1;i++){
           Aei_Aii_Aie_product_y_P->set_item(i,proc_pos_in_y,Aii_Aie_product_y_P->item(i));
        }
     }
     else{

        // Get appropriate column of Aie
        Aie_y_P -> extract_col( proc_pos_in_y-1, product_result_y_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, product_result_y_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_y_P->multiply_vec_then_add(product_result_y_P,Aii_Aie_product_y_P);
        for(i=0;i<nb_procs_in_y-1;i++){
           Aei_Aii_Aie_product_y_P->set_item(i,proc_pos_in_y-1,Aii_Aie_product_y_P->item(i));
        }

        // Get appropriate column of Aie
        Aie_y_P -> extract_col( proc_pos_in_y, product_result_y_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, product_result_y_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_y_P->multiply_vec_then_add(product_result_y_P,Aii_Aie_product_y_P);
        for(i=0;i<nb_procs_in_y-1;i++){
           Aei_Aii_Aie_product_y_P->set_item(i,proc_pos_in_y,Aii_Aie_product_y_P->item(i));
        }
      }
   }
   else
   {
      if(proc_pos_in_y==0){

        // Get appropriate column of Aie
        Aie_y_P -> extract_col( proc_pos_in_y, product_result_y_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, product_result_y_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_y_P->multiply_vec_then_add(product_result_y_P,Aii_Aie_product_y_P);
        for(i=0;i<nb_procs_in_y-1;i++){
           Aei_Aii_Aie_product_y_P->set_item(i,proc_pos_in_y,Aii_Aie_product_y_P->item(i));
        }

        // Get appropriate column of Aie
        Aie_y_P -> extract_col( nb_procs_in_y-1, product_result_y_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, product_result_y_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_y_P->multiply_vec_then_add(product_result_y_P,Aii_Aie_product_y_P);
        for(i=0;i<nb_procs_in_y-1;i++){
           Aei_Aii_Aie_product_y_P->set_item(i,nb_procs_in_y-1,Aii_Aie_product_y_P->item(i));
        }
     }
     else{

        // Get appropriate column of Aie
        Aie_y_P -> extract_col( proc_pos_in_y-1, product_result_y_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, product_result_y_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_y_P->multiply_vec_then_add(product_result_y_P,Aii_Aie_product_y_P);
        for(i=0;i<nb_procs_in_y-1;i++){
           Aei_Aii_Aie_product_y_P->set_item(i,proc_pos_in_y-1,Aii_Aie_product_y_P->item(i));
        }

        // Get appropriate column of Aie
        Aie_y_P -> extract_col( proc_pos_in_y, product_result_y_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_y_sub_diagonal_P,Aii_y_main_diagonal_P,Aii_y_mod_super_diagonal_P, product_result_y_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_y_P->multiply_vec_then_add(product_result_y_P,Aii_Aie_product_y_P);
        for(i=0;i<nb_procs_in_y-1;i++){
           Aei_Aii_Aie_product_y_P->set_item(i,proc_pos_in_y,Aii_Aie_product_y_P->item(i));
        }
      }
   }

}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_product_matrix_z_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: compute_product_matrix_z_P" ) ;

   int i;

   if(P_is_zperiodic != 1){
    if(proc_pos_in_z == nb_procs_in_z-1){

        // Get appropriate column of Aie
        Aie_z_P -> extract_col( proc_pos_in_z-1, product_result_z_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, product_result_z_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z_P->multiply_vec_then_add(product_result_z_P,Aii_Aie_product_z_P);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z_P->set_item(i,proc_pos_in_z-1,Aii_Aie_product_z_P->item(i));
        }
     }
     else if(proc_pos_in_z==0){

        // Get appropriate column of Aie
        Aie_z_P -> extract_col( proc_pos_in_z, product_result_z_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, product_result_z_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z_P->multiply_vec_then_add(product_result_z_P,Aii_Aie_product_z_P);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z_P->set_item(i,proc_pos_in_z,Aii_Aie_product_z_P->item(i));
        }
     }
     else{

        // Get appropriate column of Aie
        Aie_z_P -> extract_col( proc_pos_in_z-1, product_result_z_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, product_result_z_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z_P->multiply_vec_then_add(product_result_z_P,Aii_Aie_product_z_P);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z_P->set_item(i,proc_pos_in_z-1,Aii_Aie_product_z_P->item(i));
        }

        // Get appropriate column of Aie
        Aie_z_P -> extract_col( proc_pos_in_z, product_result_z_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, product_result_z_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z_P->multiply_vec_then_add(product_result_z_P,Aii_Aie_product_z_P);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z_P->set_item(i,proc_pos_in_z,Aii_Aie_product_z_P->item(i));
        }
      }
   }
   else{
      if(proc_pos_in_z==0){

        // Get appropriate column of Aie
        Aie_z_P -> extract_col( proc_pos_in_z, product_result_z_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, product_result_z_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z_P->multiply_vec_then_add(product_result_z_P,Aii_Aie_product_z_P);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z_P->set_item(i,proc_pos_in_z,Aii_Aie_product_z_P->item(i));
        }

        // Get appropriate column of Aie
        Aie_z_P -> extract_col( nb_procs_in_z-1, product_result_z_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, product_result_z_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z_P->multiply_vec_then_add(product_result_z_P,Aii_Aie_product_z_P);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z_P->set_item(i,nb_procs_in_z-1,Aii_Aie_product_z_P->item(i));
        }
     }
     else{

        // Get appropriate column of Aie
        Aie_z_P -> extract_col( proc_pos_in_z-1, product_result_z_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, product_result_z_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z_P->multiply_vec_then_add(product_result_z_P,Aii_Aie_product_z_P);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z_P->set_item(i,proc_pos_in_z-1,Aii_Aie_product_z_P->item(i));
        }

        // Get appropriate column of Aie
        Aie_z_P -> extract_col( proc_pos_in_z, product_result_z_P );

        // Get inv(Aii)*Aie for for appropriate column of Aie
        DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_z_sub_diagonal_P,Aii_z_main_diagonal_P,Aii_z_mod_super_diagonal_P, product_result_z_P);

        // Get product of Aei*inv(Aii)*Aie for appropriate column
        Aei_z_P->multiply_vec_then_add(product_result_z_P,Aii_Aie_product_z_P);
        for(i=0;i<nb_procs_in_z-1;i++){
           Aei_Aii_Aie_product_z_P->set_item(i,proc_pos_in_z,Aii_Aie_product_z_P->item(i));
        }
      }
   }
}




//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::display_debug(void)
//----------------------------------------------------------------------
{
  //VEC_DS_PF->print_items(MAC::out(),0);
  // if(proc_pos_in_x == 0)
  // Aii_z_main_diagonal[0]->print_items(MAC::out(),0);

}
