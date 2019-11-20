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
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DDS_NavierStokesSystem" ) ;

   int const* MPI_coordinates_world = UF->primary_grid()->get_MPI_coordinates() ;
   int const* MPI_max_coordinates_world = UF->primary_grid()->get_domain_decomposition() ;

   is_periodic[0][0] = false;
   is_periodic[0][1] = false;
   is_periodic[0][2] = false;
   is_periodic[1][0] = false;
   is_periodic[1][1] = false;
   is_periodic[1][2] = false;

   proc_pos_in_i[0] = MPI_coordinates_world[0];
   proc_pos_in_i[1] = MPI_coordinates_world[1];
   proc_pos_in_i[2] = MPI_coordinates_world[2];
   nb_procs_in_i[0] = MPI_max_coordinates_world[0];
   nb_procs_in_i[1] = MPI_max_coordinates_world[1];
   nb_procs_in_i[2] = MPI_max_coordinates_world[2];

   dim = UF->primary_grid()->nb_space_dimensions() ;
   nb_comps = UF->nb_components() ;

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

if(proc_pos_in_i[0] == 0){
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

if(proc_pos_in_i[1] == 0){
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

     if(proc_pos_in_i[2] == 0){
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

	  if(proc_pos_in_i[0] == 0){
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

    if(proc_pos_in_i[1] == 0){
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

       if(proc_pos_in_i[2] == 0){
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

  if(is_periodic[0][0]==0 || nb_procs_in_i[0] == 1)
  {
     if(proc_pos_in_i[0] == nb_procs_in_i[0]-1)
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
        nb_procs_in_i[0]-1);
        Aei_x_P->re_initialize(
        nb_procs_in_i[0]-1,
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
        nb_procs_in_i[0]-1 );
        Aei_x_P->re_initialize(
        nb_procs_in_i[0]-1,
        nb_unknowns_handled_by_proc( 0 )-1 );
        product_result_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1 );
        VEC_local_temp_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
        VEC_local_solution_temp_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
      }
      Aii_Aie_product_x_P->re_initialize(
       nb_procs_in_i[0]-1);

     Aei_Aii_Aie_product_x_P->re_initialize(
        nb_procs_in_i[0]-1,
        nb_procs_in_i[0]-1 );

     VEC_interface_temp_x_P->re_initialize( nb_procs_in_i[0]-1 ) ;
     VEC_temp_x_P->re_initialize( nb_procs_in_i[0]-1 ) ;

     if(proc_pos_in_i[0] == 0){
       Aee_x_P->re_initialize(
          nb_procs_in_i[0]-1,
          nb_procs_in_i[0]-1 );
       schlur_complement_x_P->re_initialize(
            nb_procs_in_i[0]-1,
            nb_procs_in_i[0]-1 );
       schlur_complement_x_ref_P->re_initialize(
            nb_procs_in_i[0]-1,
            nb_procs_in_i[0]-1 );
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
    nb_procs_in_i[0] );
    Aei_x_P->re_initialize(
    nb_procs_in_i[0],
    nb_unknowns_handled_by_proc( 0 )-1 );
    product_result_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1 );
    VEC_local_temp_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
    VEC_local_solution_temp_x_P->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
      
    Aii_Aie_product_x_P->re_initialize(
       nb_procs_in_i[0]);

    Aei_Aii_Aie_product_x_P->re_initialize(
        nb_procs_in_i[0],
        nb_procs_in_i[0] );

    VEC_interface_temp_x_P->re_initialize( nb_procs_in_i[0] ) ;
    VEC_temp_x_P->re_initialize( nb_procs_in_i[0] ) ;

    if(proc_pos_in_i[0] == 0){
       Aee_x_P->re_initialize(
          nb_procs_in_i[0],
          nb_procs_in_i[0] );
       schlur_complement_x_P->re_initialize(
            nb_procs_in_i[0],
            nb_procs_in_i[0] );
       schlur_complement_x_ref_P->re_initialize(
            nb_procs_in_i[0],
            nb_procs_in_i[0] );
    }
  }
  if(is_periodic[0][1]!=1 || nb_procs_in_i[1] == 1)
  {
      if(proc_pos_in_i[1] == nb_procs_in_i[1]-1){
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
        nb_procs_in_i[1]-1);
        Aei_y_P->re_initialize(
        nb_procs_in_i[1]-1,
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
        nb_procs_in_i[1]-1 );
        Aei_y_P->re_initialize(
        nb_procs_in_i[1]-1,
        nb_unknowns_handled_by_proc( 1 )-1 );
        product_result_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1 );
        VEC_local_temp_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
        VEC_local_solution_temp_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
     }

      Aii_Aie_product_y_P->re_initialize(
         nb_procs_in_i[1]-1);

      Aei_Aii_Aie_product_y_P->re_initialize(
           nb_procs_in_i[1]-1,
           nb_procs_in_i[1]-1 );


      if(proc_pos_in_i[1] == 0){
         Aee_y_P->re_initialize(
              nb_procs_in_i[1]-1,
              nb_procs_in_i[1]-1 );
         schlur_complement_y_P->re_initialize(
              nb_procs_in_i[1]-1,
              nb_procs_in_i[1]-1 );
         schlur_complement_y_ref_P->re_initialize(
              nb_procs_in_i[1]-1,
              nb_procs_in_i[1]-1 );
      }

     VEC_interface_temp_y_P->re_initialize( nb_procs_in_i[1]-1 ) ;
     VEC_temp_y_P->re_initialize( nb_procs_in_i[1]-1 ) ;
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
      nb_procs_in_i[1] );
      Aei_y_P->re_initialize(
      nb_procs_in_i[1],
      nb_unknowns_handled_by_proc( 1 )-1 );
      product_result_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1 );
      VEC_local_temp_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
      VEC_local_solution_temp_y_P->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;

      Aii_Aie_product_y_P->re_initialize(
           nb_procs_in_i[1]);

      Aei_Aii_Aie_product_y_P->re_initialize(
             nb_procs_in_i[1],
             nb_procs_in_i[1] );


      if(proc_pos_in_i[1] == 0){
         Aee_y_P->re_initialize(
              nb_procs_in_i[1],
              nb_procs_in_i[1] );
         schlur_complement_y_P->re_initialize(
              nb_procs_in_i[1],
              nb_procs_in_i[1] );
         schlur_complement_y_ref_P->re_initialize(
              nb_procs_in_i[1],
              nb_procs_in_i[1] );
      }

      VEC_interface_temp_y_P->re_initialize( nb_procs_in_i[1] ) ;
      VEC_temp_y_P->re_initialize( nb_procs_in_i[1]) ;
   }
   if(dim>2)
   {
      if(is_periodic[0][2] !=1 || nb_procs_in_i[2] == 1)
      {

        if(proc_pos_in_i[2] == nb_procs_in_i[2]-1){

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
           nb_procs_in_i[2]-1);
           Aei_z_P->re_initialize(
           nb_procs_in_i[2]-1,
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
           nb_procs_in_i[2]-1 );
           Aei_z_P->re_initialize(
           nb_procs_in_i[2]-1,
           nb_unknowns_handled_by_proc( 2 )-1 );
           product_result_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1 );
           VEC_local_temp_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
           VEC_local_solution_temp_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
        }

        if(proc_pos_in_i[2] == 0){
            Aee_z_P->re_initialize(
                 nb_procs_in_i[2]-1,
                 nb_procs_in_i[2]-1 );
            schlur_complement_z_P->re_initialize(
                 nb_procs_in_i[2]-1,
                 nb_procs_in_i[2]-1 );
            schlur_complement_z_ref_P->re_initialize(
                 nb_procs_in_i[2]-1,
                 nb_procs_in_i[2]-1 );
        }

        Aii_Aie_product_z_P->re_initialize(
           nb_procs_in_i[2]-1);

        Aei_Aii_Aie_product_z_P->re_initialize(
           nb_procs_in_i[2]-1,
           nb_procs_in_i[2]-1 );

        VEC_interface_temp_z_P->re_initialize( nb_procs_in_i[2]-1 ) ;
        VEC_temp_z_P->re_initialize( nb_procs_in_i[2]-1 ) ;
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
         nb_procs_in_i[2]);
         Aei_z_P->re_initialize(
         nb_procs_in_i[2],
         nb_unknowns_handled_by_proc( 2 )-1 );
         product_result_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1 );
         VEC_local_temp_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
         VEC_local_solution_temp_z_P->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
         
         if(proc_pos_in_i[2] == 0){
              Aee_z_P->re_initialize(
                   nb_procs_in_i[2],
                   nb_procs_in_i[2] );
              schlur_complement_z_P->re_initialize(
                   nb_procs_in_i[2],
                   nb_procs_in_i[2] );
              schlur_complement_z_ref_P->re_initialize(
                   nb_procs_in_i[2],
                   nb_procs_in_i[2] );
          }

         Aii_Aie_product_z_P->re_initialize(
             nb_procs_in_i[2]);

         Aei_Aii_Aie_product_z_P->re_initialize(
             nb_procs_in_i[2],
             nb_procs_in_i[2] );

         VEC_interface_temp_z_P->re_initialize( nb_procs_in_i[2] ) ;
         VEC_temp_z_P->re_initialize( nb_procs_in_i[2] ) ;
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
   if(is_periodic[1][0]!=1 || nb_procs_in_i[0] == 1){
      if(proc_pos_in_i[0] == nb_procs_in_i[0]-1)
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
          nb_procs_in_i[0]-1);
          Aei_x[comp]->re_initialize(
          nb_procs_in_i[0]-1,
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
          nb_procs_in_i[0]-1 );
          Aei_x[comp]->re_initialize(
          nb_procs_in_i[0]-1,
          nb_unknowns_handled_by_proc( 0 )-1 );
          product_result_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1 );
          VEC_local_temp_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
          VEC_local_solution_temp_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
        }
        Aii_Aie_product_x[comp]->re_initialize(
         nb_procs_in_i[0]-1);

       Aei_Aii_Aie_product_x[comp]->re_initialize(
          nb_procs_in_i[0]-1,
          nb_procs_in_i[0]-1 );

       VEC_interface_temp_x[comp]->re_initialize( nb_procs_in_i[0]-1 ) ;
       VEC_temp_x[comp]->re_initialize( nb_procs_in_i[0]-1 ) ;

       if(proc_pos_in_i[0] == 0){
         Aee_x[comp]->re_initialize(
            nb_procs_in_i[0]-1,
            nb_procs_in_i[0]-1 );
           schlur_complement_x[comp]->re_initialize(
                nb_procs_in_i[0]-1,
                nb_procs_in_i[0]-1 );
           schlur_complement_x_ref[comp]->re_initialize(
                nb_procs_in_i[0]-1,
                nb_procs_in_i[0]-1 );
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
      nb_procs_in_i[0] );
      Aei_x[comp]->re_initialize(
      nb_procs_in_i[0],
      nb_unknowns_handled_by_proc( 0 )-1 );
      product_result_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1 );
      VEC_local_temp_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
      VEC_local_solution_temp_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 )-1) ;
        
      Aii_Aie_product_x[comp]->re_initialize(
         nb_procs_in_i[0]);

      Aei_Aii_Aie_product_x[comp]->re_initialize(
          nb_procs_in_i[0],
          nb_procs_in_i[0] );

      VEC_interface_temp_x[comp]->re_initialize( nb_procs_in_i[0] ) ;
      VEC_temp_x[comp]->re_initialize( nb_procs_in_i[0] ) ;

      if(proc_pos_in_i[0] == 0){
         Aee_x[comp]->re_initialize(
            nb_procs_in_i[0],
            nb_procs_in_i[0] );
         schlur_complement_x[comp]->re_initialize(
              nb_procs_in_i[0],
              nb_procs_in_i[0] );
         schlur_complement_x_ref[comp]->re_initialize(
              nb_procs_in_i[0],
              nb_procs_in_i[0] );
      }

      VEC_rhs_velocity_1D_x[comp]->re_initialize( nb_unknowns_handled_by_proc( 0 ) );
   }
       
   if(is_periodic[1][1]!=1 || nb_procs_in_i[1] == 1)
   {
      if(proc_pos_in_i[1] == nb_procs_in_i[1]-1){
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
          nb_procs_in_i[1]-1);
          Aei_y[comp]->re_initialize(
          nb_procs_in_i[1]-1,
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
          nb_procs_in_i[1]-1 );
          Aei_y[comp]->re_initialize(
          nb_procs_in_i[1]-1,
          nb_unknowns_handled_by_proc( 1 )-1 );
          product_result_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1 );
          VEC_local_temp_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
          VEC_local_solution_temp_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
       }

       MAT_velocityUnsteadyPlusDiffusion_1D_y->re_initialize(
          nb_unknowns_handled_by_proc( 1 ),
          nb_unknowns_handled_by_proc( 1 ) );

        Aii_Aie_product_y[comp]->re_initialize(
           nb_procs_in_i[1]-1);

        Aei_Aii_Aie_product_y[comp]->re_initialize(
             nb_procs_in_i[1]-1,
             nb_procs_in_i[1]-1 );


        if(proc_pos_in_i[1] == 0){
           Aee_y[comp]->re_initialize(
                nb_procs_in_i[1]-1,
                nb_procs_in_i[1]-1 );
           schlur_complement_y[comp]->re_initialize(
                nb_procs_in_i[1]-1,
                nb_procs_in_i[1]-1 );
           schlur_complement_y_ref[comp]->re_initialize(
                nb_procs_in_i[1]-1,
                nb_procs_in_i[1]-1 );
        }

       VEC_rhs_velocity_1D_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 ) );

       VEC_interface_temp_y[comp]->re_initialize( nb_procs_in_i[1]-1 ) ;
       VEC_temp_y[comp]->re_initialize( nb_procs_in_i[1]-1 ) ;
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
          nb_procs_in_i[1]);
          Aei_y[comp]->re_initialize(
          nb_procs_in_i[1],
          nb_unknowns_handled_by_proc( 1 )-1 );
          product_result_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1 );
          VEC_local_temp_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
          VEC_local_solution_temp_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 )-1) ;
       
       MAT_velocityUnsteadyPlusDiffusion_1D_y->re_initialize(
          nb_unknowns_handled_by_proc( 1 ),
          nb_unknowns_handled_by_proc( 1 ) );

        Aii_Aie_product_y[comp]->re_initialize(
           nb_procs_in_i[1]);

        Aei_Aii_Aie_product_y[comp]->re_initialize(
             nb_procs_in_i[1],
             nb_procs_in_i[1] );


        if(proc_pos_in_i[1] == 0){
           Aee_y[comp]->re_initialize(
                nb_procs_in_i[1],
                nb_procs_in_i[1] );
           schlur_complement_y[comp]->re_initialize(
                nb_procs_in_i[1],
                nb_procs_in_i[1] );
           schlur_complement_y_ref[comp]->re_initialize(
                nb_procs_in_i[1],
                nb_procs_in_i[1] );
        }

       VEC_rhs_velocity_1D_y[comp]->re_initialize( nb_unknowns_handled_by_proc( 1 ) );

       VEC_interface_temp_y[comp]->re_initialize( nb_procs_in_i[1]) ;
       VEC_temp_y[comp]->re_initialize( nb_procs_in_i[1]) ;
   }
       

   if(dim>2)
   {
      if(is_periodic[1][2]!=1 || nb_procs_in_i[2] == 1)
      {
        if(proc_pos_in_i[2] == nb_procs_in_i[2]-1){
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
             nb_procs_in_i[2]-1);
             Aei_z[comp]->re_initialize(
             nb_procs_in_i[2]-1,
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
             nb_procs_in_i[2]-1 );
             Aei_z[comp]->re_initialize(
             nb_procs_in_i[2]-1,
             nb_unknowns_handled_by_proc( 2 )-1 );
             product_result_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1 );
             VEC_local_temp_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
             VEC_local_solution_temp_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
          }

          if(proc_pos_in_i[2] == 0){
              Aee_z[comp]->re_initialize(
                   nb_procs_in_i[2]-1,
                   nb_procs_in_i[2]-1 );
              schlur_complement_z[comp]->re_initialize(
                   nb_procs_in_i[2]-1,
                   nb_procs_in_i[2]-1 );
              schlur_complement_z_ref[comp]->re_initialize(
                   nb_procs_in_i[2]-1,
                   nb_procs_in_i[2]-1 );
          }

          Aii_Aie_product_z[comp]->re_initialize(
             nb_procs_in_i[2]-1);

          Aei_Aii_Aie_product_z[comp]->re_initialize(
             nb_procs_in_i[2]-1,
             nb_procs_in_i[2]-1 );

          VEC_interface_temp_z[comp]->re_initialize( nb_procs_in_i[2]-1 ) ;
          VEC_temp_z[comp]->re_initialize( nb_procs_in_i[2]-1 ) ;

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
           nb_procs_in_i[2]);
           Aei_z[comp]->re_initialize(
           nb_procs_in_i[2],
           nb_unknowns_handled_by_proc( 2 )-1 );
           product_result_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1 );
           VEC_local_temp_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;
           VEC_local_solution_temp_z[comp]->re_initialize( nb_unknowns_handled_by_proc( 2 )-1) ;

           if(proc_pos_in_i[2] == 0){
              Aee_z[comp]->re_initialize(
                   nb_procs_in_i[2],
                   nb_procs_in_i[2] );
              schlur_complement_z[comp]->re_initialize(
                   nb_procs_in_i[2],
                   nb_procs_in_i[2] );
              schlur_complement_z_ref[comp]->re_initialize(
                   nb_procs_in_i[2],
                   nb_procs_in_i[2] );
           }

           Aii_Aie_product_z[comp]->re_initialize(
             nb_procs_in_i[2]);

           Aei_Aii_Aie_product_z[comp]->re_initialize(
             nb_procs_in_i[2],
             nb_procs_in_i[2] );

           VEC_interface_temp_z[comp]->re_initialize( nb_procs_in_i[2] ) ;
           VEC_temp_z[comp]->re_initialize( nb_procs_in_i[2] ) ;

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
   // if(proc_pos_in_i[0] == 0 && proc_pos_in_i[1] == 0)
   //    MAC::out()<<"Time change "<<time_change<<endl;
   if ( norm_UF > 1e-4 ) time_change /= norm_UF;

   return ( time_change ) ;

}

//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_temp( size_t const& c, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_x_P" ) ;
         return ( VEC_local_temp_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_y_P" ) ;
         return ( VEC_local_temp_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_z_P" ) ;
         return ( VEC_local_temp_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_x" ) ;
         return ( VEC_local_temp_x[c] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_y" ) ;
         return ( VEC_local_temp_y[c] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_temp_z" ) ;
         return ( VEC_local_temp_z[c] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_local_solution_temp( size_t const& c, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_solution_temp_x_P" ) ;
         return ( VEC_local_solution_temp_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_solution_temp_y_P" ) ;
         return ( VEC_local_solution_temp_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_solution_temp_z_P" ) ;
         return ( VEC_local_solution_temp_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_solution_temp_x" ) ;
         return ( VEC_local_solution_temp_x[c] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_solution_temp_y" ) ;
         return ( VEC_local_solution_temp_y[c] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_local_solution_temp_z" ) ;
         return ( VEC_local_solution_temp_z[c] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_interface_temp( size_t const& c, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_x_P" ) ;
         return ( VEC_interface_temp_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_y_P" ) ;
         return ( VEC_interface_temp_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_z_P" ) ;
         return ( VEC_interface_temp_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_x" ) ;
         return ( VEC_interface_temp_x[c] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_y" ) ;
         return ( VEC_interface_temp_y[c] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_z" ) ;
         return ( VEC_interface_temp_z[c] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_temp( size_t const& c, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) { 
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_temp_x_P" ) ;
         return ( VEC_temp_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_temp_y_P" ) ;
         return ( VEC_temp_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_temp_z_P" ) ;
         return ( VEC_temp_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_x" ) ;
         return ( VEC_temp_x[c] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_y" ) ;
         return ( VEC_temp_y[c] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_interface_temp_z" ) ;
         return ( VEC_temp_z[c] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_product_result( size_t const& c, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_product_result_x_P" ) ;
         return ( product_result_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_product_result_y_P" ) ;
         return ( product_result_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_product_result_z_P" ) ;
         return ( product_result_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_product_result_x" ) ;
         return ( product_result_x[c] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_product_result_y" ) ;
         return ( product_result_y[c] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_product_result_z" ) ;
         return ( product_result_z[c] ) ;
      }
   }
}

//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_Aii_ref( size_t const& comp, size_t const dir, size_t const field)
//----------------------------------------------------------------------
{
   LA_SeqVector* Aii_main_diagonal = get_aii_main_diag(comp,dir,field);
   LA_SeqVector* Aii_super_diagonal = get_aii_super_diag(comp,dir,field);
   LA_SeqVector* Aii_sub_diagonal = get_aii_sub_diag(comp,dir,field);
   LA_SeqVector* Aii_mod_super_diagonal = get_aii_mod_super_diag(comp,dir,field);

   size_t nrows = Aii_main_diagonal->nb_rows() ;
   double temp = Aii_main_diagonal->item(0);
   Aii_mod_super_diagonal->set_item(0,Aii_super_diagonal->item(0)/temp);

   //  // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m) {
      double a,b,c,prevc;
      a=Aii_sub_diagonal->item(m-1);
      b=Aii_main_diagonal->item(m);
      prevc=Aii_mod_super_diagonal->item(m-1);

      if (m<nrows-1){
         c=Aii_super_diagonal->item(m);
         Aii_mod_super_diagonal->set_item(m,c/(b - a*prevc));
      }
   }
}

//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_schlur_ref( size_t const& comp, size_t const dir, size_t const field)
//----------------------------------------------------------------------
{
   LA_SeqMatrix* schlur_complement = get_schlur_complement(comp,dir,field);	
   LA_SeqMatrix* schlur_complement_ref = get_schlur_complement_ref(comp,dir,field);	
   size_t nrows = schlur_complement->nb_rows() ;
   double temp = schlur_complement->item(0,0);
   schlur_complement_ref-> set_item(0,0,schlur_complement->item(0,0));
   schlur_complement_ref-> set_item(0,1,schlur_complement->item(0,1)/temp);

    // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a=schlur_complement->item(m,m-1);
     b=schlur_complement->item(m,m);
     c=schlur_complement->item(m,m+1);
     prevc=schlur_complement_ref->item(m-1,m);

     if(m<nrows-1)
         schlur_complement_ref->set_item(m,m+1,c/(b - a*prevc));
      schlur_complement_ref->set_item(m,m,b);
      schlur_complement_ref->set_item(m,m-1,a);
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

   if(nrows>1){
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
DDS_NavierStokesSystem::get_aii_main_diag( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_main_diag_in_P_in_x" ) ;
         return ( Aii_x_main_diagonal_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_main_diag_in_P_in_y" ) ;
         return ( Aii_y_main_diagonal_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_main_diag_in_P_in_z" ) ;
         return ( Aii_z_main_diagonal_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_main_diag_in_x" ) ;
         return ( Aii_x_main_diagonal[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_main_diag_in_y" ) ;
         return ( Aii_y_main_diagonal[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_main_diag_in_z" ) ;
         return ( Aii_z_main_diagonal[comp] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_super_diag( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_super_diag_in_x_in_P" ) ;
         return ( Aii_x_super_diagonal_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_super_diag_in_y_in_P" ) ;
         return ( Aii_y_super_diagonal_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_super_diag_in_z_in_P" ) ;
         return ( Aii_z_super_diagonal_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_super_diag_in_x" ) ;
         return ( Aii_x_super_diagonal[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_super_diag_in_y" ) ;
         return ( Aii_y_super_diagonal[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_super_diag_in_z" ) ;
         return ( Aii_z_super_diagonal[comp] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_mod_super_diag( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_super_diag_in_x_in_P" ) ;
         return ( Aii_x_mod_super_diagonal_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_super_diag_in_y_in_P" ) ;
         return ( Aii_y_mod_super_diagonal_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_super_diag_in_z_in_P" ) ;
         return ( Aii_z_mod_super_diagonal_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_mod_super_diag_in_x" ) ;
         return ( Aii_x_mod_super_diagonal[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_mod_super_diag_in_y" ) ;
         return ( Aii_y_mod_super_diagonal[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_mod_super_diag_in_z" ) ;
         return ( Aii_z_mod_super_diagonal[comp] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_aii_sub_diag( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_sub_diag_in_x_in_P" ) ;
         return ( Aii_x_sub_diagonal_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_sub_diag_in_y_in_P" ) ;
         return ( Aii_y_sub_diagonal_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_sub_diag_in_z_in_P" ) ;
         return ( Aii_z_sub_diagonal_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_sub_diag_in_x" ) ;
         return ( Aii_x_sub_diagonal[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_sub_diag_in_y" ) ;
         return ( Aii_y_sub_diagonal[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aii_sub_diag_in_z" ) ;
         return ( Aii_z_sub_diagonal[comp] ) ;
      }
   }

}

//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aie( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_x_in_P" ) ;
         return ( Aie_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_y_in_P" ) ;
         return ( Aie_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_z_in_P" ) ;
         return ( Aie_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_x" ) ;
         return ( Aie_x[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_y" ) ;
         return ( Aie_y[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aie_in_z" ) ;
         return ( Aie_z[comp] ) ;
      }
   } 
}

//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_aei( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_x_in_P" ) ;
         return ( Aei_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_y_in_P" ) ;
         return ( Aei_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_z_in_P" ) ;
         return ( Aei_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_x" ) ;
         return ( Aei_x[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_y" ) ;
         return ( Aei_y[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aei_in_z" ) ;
         return ( Aei_z[comp] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aee_matrix( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aee_in_x_in_P" ) ;
         return ( Aee_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aee_in_y_in_P" ) ;
         return ( Aee_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aee_in_z_in_P" ) ;
         return ( Aee_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aee_in_x" ) ;
         return ( Aee_x[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aee_in_y" ) ;
         return ( Aee_y[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_aee_in_z" ) ;
         return ( Aee_z[comp] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_schlur_complement( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_x_in_P" ) ;
         return ( schlur_complement_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_y_in_P" ) ;
         return ( schlur_complement_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_z_in_P" ) ;
         return ( schlur_complement_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_x" ) ;
         return ( schlur_complement_x[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_y" ) ;
         return ( schlur_complement_y[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_in_z" ) ;
         return ( schlur_complement_z[comp] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_schlur_complement_ref( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_x_ref_P" ) ;
         return ( schlur_complement_x_ref_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_y_ref_P" ) ;
         return ( schlur_complement_y_ref_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_z_ref_P" ) ;
         return ( schlur_complement_z_ref_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_x_ref" ) ;
         return ( schlur_complement_x_ref[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_y_ref" ) ;
         return ( schlur_complement_y_ref[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_schlur_complement_z_ref" ) ;
         return ( schlur_complement_z_ref[comp] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqMatrix*
DDS_NavierStokesSystem::get_Aei_Aii_Aie_product( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_x_in_P" ) ;
         return ( Aei_Aii_Aie_product_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_y_in_P" ) ;
         return ( Aei_Aii_Aie_product_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_z_in_P" ) ;
         return ( Aei_Aii_Aie_product_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_x" ) ;
         return ( Aei_Aii_Aie_product_x[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_y" ) ;
         return ( Aei_Aii_Aie_product_y[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aei_Aii_Aie_product_in_z" ) ;
         return ( Aei_Aii_Aie_product_z[comp] ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_SeqVector*
DDS_NavierStokesSystem::get_Aii_Aie_product( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   if (field == 0) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aii_Aie_product_in_x_in_P" ) ;
         return ( Aii_Aie_product_x_P ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aii_Aie_product_in_y_in_P" ) ;
         return ( Aii_Aie_product_y_P ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aii_Aie_product_in_z_in_P" ) ;
         return ( Aii_Aie_product_z_P ) ;
      }
   } else if (field == 1) {
      if (dir == 0) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aii_Aie_product_in_x" ) ;
         return ( Aii_Aie_product_x[comp] ) ;
      } else if (dir == 1) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aii_Aie_product_in_y" ) ;
         return ( Aii_Aie_product_y[comp] ) ;
      } else if (dir == 2) {
         MAC_LABEL( "DDS_NavierStokesSystem:: get_Aii_Aie_product_in_z" ) ;
         return ( Aii_Aie_product_z[comp] ) ;
      }
   }
}

//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_local_unknown_solver( LA_SeqVector* rhs, size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_local_unknown_solver" ) ;

   LA_SeqVector* Aii_main_diagonal = get_aii_main_diag(comp,dir,field);
   LA_SeqVector* Aii_sub_diagonal = get_aii_sub_diag(comp,dir,field);
   LA_SeqVector* Aii_mod_super_diagonal = get_aii_mod_super_diag(comp,dir,field);

   // Solve the DS splitting problem in
   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_sub_diagonal,Aii_main_diagonal,Aii_mod_super_diagonal, rhs);

}

//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_interface_unknown_solver( LA_SeqVector* rhs, size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_x_interface_unknown_solver" ) ;

   LA_SeqMatrix* schlur_complement_ref = get_schlur_complement_ref(comp,dir,field);

   // Solve the DS splitting problem in
   DDS_NavierStokesSystem::thomas_algorithm( schlur_complement_ref, rhs);

}

//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::DS_NavierStokes_solver(FV_DiscreteField* FF
	,size_t const& j, size_t const& k, size_t const& min_i, LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: DS_NavierStokes_solver" ) ;

   LA_SeqVector* Aii_main_diagonal = get_aii_main_diag(comp,dir,field);
   LA_SeqVector* Aii_sub_diagonal = get_aii_sub_diag(comp,dir,field);
   LA_SeqVector* Aii_mod_super_diagonal = get_aii_mod_super_diag(comp,dir,field);
   // Solve the DS splitting problem in

   DDS_NavierStokesSystem::mod_thomas_algorithm( Aii_sub_diagonal,Aii_main_diagonal,Aii_mod_super_diagonal, rhs);

   // Transfer in the distributed vector
   size_t nb_local_unk = rhs->nb_rows();
   // Since, this function is used in all directions;
   // ii, jj, and kk are used to convert the passed arguments corresponding to correct direction
   size_t ii=0,jj=0,kk=0;
   
   size_t m, i, global_number_in_distributed_vector;

   for (m=0;m<nb_local_unk;++m) {
      i = min_i + m ;

      if (dir == 0) {
         ii = i; jj = j; kk = k;
      } else if (dir == 1) {
         ii = j; jj = i; kk = k;
      } else if (dir == 2) {
         ii = j; jj = k; kk = i;
      }      

      global_number_in_distributed_vector = FF->DOF_global_number( ii, jj, kk, comp );
      if (field == 0) {
         VEC_DS_PF->set_item( global_number_in_distributed_vector, rhs->item( m ) );
      } else if (field == 1) { 
         VEC_DS_UF->set_item( global_number_in_distributed_vector, rhs->item( m ) );
      }
   }

   // Put the interface unknowns in distributed vector
   i = min_i + nb_local_unk ;
   if (dir == 0) {
      ii = i; jj = j; kk = k;
   } else if (dir == 1) {
      ii = j; jj = i; kk = k;
   } else if (dir == 2) {
      ii = j; jj = k; kk = i;
   }

   if ((is_periodic[field][dir] == 1)) {
      global_number_in_distributed_vector = FF->DOF_global_number( ii, jj, kk, comp );
      if (field == 0) {
         VEC_DS_PF->set_item( global_number_in_distributed_vector,interface_rhs->item( proc_pos_in_i[dir] ) );
      } else if (field == 1) {
         VEC_DS_UF->set_item( global_number_in_distributed_vector, interface_rhs->item( proc_pos_in_i[dir] ) );
      }
   } else if ((is_periodic[field][dir] == 0) && (proc_pos_in_i[dir] != nb_procs_in_i[dir]-1)) {
      global_number_in_distributed_vector = FF->DOF_global_number( ii, jj, kk, comp );
      if (field == 0) {
         VEC_DS_PF->set_item( global_number_in_distributed_vector,interface_rhs->item( proc_pos_in_i[dir] ) );
      } else if (field == 1) {
         VEC_DS_UF->set_item( global_number_in_distributed_vector, interface_rhs->item( proc_pos_in_i[dir] ) );
      }
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
DDS_NavierStokesSystem::compute_product_matrix_interior(size_t const& comp, size_t const column, size_t const dir, size_t const field)
//----------------------------------------------------------------------
{
   LA_SeqVector* Aii_main_diagonal = get_aii_main_diag(comp,dir,field);
   LA_SeqVector* Aii_mod_super_diagonal = get_aii_mod_super_diag(comp,dir,field);
   LA_SeqVector* Aii_sub_diagonal = get_aii_sub_diag(comp,dir,field);
   LA_SeqVector* product_result = get_product_result(comp,dir,field);
   LA_SeqVector* Aii_Aie = get_Aii_Aie_product(comp,dir,field);
   LA_SeqMatrix* Aei_Aii_Aie = get_Aei_Aii_Aie_product(comp,dir,field);

   LA_SeqMatrix* Aie = get_aie(comp,dir,field);
   LA_SeqMatrix* Aei = get_aei(comp,dir,field);


  // Get appropriate column of Aie
  Aie->extract_col(column, product_result);

  // Get inv(Aii)*Aie for for appropriate column of Aie
  mod_thomas_algorithm(Aii_sub_diagonal,Aii_main_diagonal,Aii_mod_super_diagonal, product_result);

  // Get product of Aei*inv(Aii)*Aie for appropriate column
  Aei->multiply_vec_then_add(product_result,Aii_Aie);

  size_t nb_procs;

  nb_procs = nb_procs_in_i[dir];

  size_t int_unknown = Aii_Aie->nb_rows();

  for (size_t i = 0; i < int_unknown; i++){
      Aei_Aii_Aie->set_item(i,column,Aii_Aie->item(i));
  }

}

//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::compute_product_matrix( size_t const& comp, size_t const dir, size_t const field )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokesSystem:: compute_product_matrix" ) ;

   size_t proc_pos, nb_procs;

   proc_pos = proc_pos_in_i[dir];
   nb_procs = nb_procs_in_i[dir];

   if (proc_pos == nb_procs - 1){
      // Condition for serial processor and multi processor
      if (proc_pos == 0) {
         compute_product_matrix_interior(comp,proc_pos,dir,field);
      } else {
         compute_product_matrix_interior(comp,proc_pos-1,dir,field);
         if (is_periodic[field][dir] == 1) compute_product_matrix_interior(comp,proc_pos,dir,field);
      }
   }else if(proc_pos == 0){
      compute_product_matrix_interior(comp,proc_pos,dir,field);
      if (is_periodic[field][dir] == 1) compute_product_matrix_interior(comp,nb_procs-1,dir,field);
   }else{
      compute_product_matrix_interior(comp,proc_pos-1,dir,field);
      compute_product_matrix_interior(comp,proc_pos,dir,field);
   }
}

//----------------------------------------------------------------------
void
DDS_NavierStokesSystem::display_debug(void)
//----------------------------------------------------------------------
{
  // VEC_DS_UF->print_items(MAC::out(),0);
  // if(proc_pos_in_i[0] == 0)
  //   schlur_complement_x[1]->print_items(MAC::out(),0);
  // Aii_x_super_diagonal_P->print_items(MAC::out(),0);

}
