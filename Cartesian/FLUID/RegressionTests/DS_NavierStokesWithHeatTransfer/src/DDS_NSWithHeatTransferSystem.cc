#include <DDS_NSWithHeatTransferSystem.hh>
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
DDS_NSWithHeatTransferSystem*
DDS_NSWithHeatTransferSystem:: create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_UF,
        FV_DiscreteField* mac_PF,
        struct NavierStokes2System const& transfer )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   MAC_CHECK_PRE( mac_UF != 0 ) ;

   DDS_NSWithHeatTransferSystem* result =
         new DDS_NSWithHeatTransferSystem( a_owner, exp, mac_UF, mac_PF, transfer ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;

   return( result ) ;

}




//----------------------------------------------------------------------
DDS_NSWithHeatTransferSystem:: DDS_NSWithHeatTransferSystem(
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_UF,
        FV_DiscreteField* mac_PF,
        struct NavierStokes2System const& fromNS )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , UF( mac_UF )
   , PF( mac_PF )
   , MAT_velocityUnsteadyPlusDiffusion_1D( 0 )
   , is_solids ( fromNS.is_solids_ )
   , is_stressCal (fromNS.is_stressCal_ )
   , Npart (fromNS.Npart_ )
   , level_set_type (fromNS.level_set_type_ )
   , Nmax (fromNS.Npoints_ )
   , ar (fromNS.ar_ ) 
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: DDS_NSWithHeatTransferSystem" ) ;

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
   nb_comps[0] = PF->nb_components() ;
   nb_comps[1] = UF->nb_components() ;

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
DDS_NSWithHeatTransferSystem:: build_system( MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: build_system" ) ;

   // velocity Laplacian
   MAT_D_velocityUnsteadyPlusDiffusion = LA_Matrix::make( this,exp->create_subexplorer( this,"MAT_D_velocityDiffusion"  ) ) ;
   VEC_rhs_D_velocityDiffusionPlusBodyTerm = MAT_D_velocityUnsteadyPlusDiffusion->create_vector(this) ;

   // Unknowns vectors
   VEC_DS_UF = MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;
   VEC_DS_UF_previoustime = MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;
   VEC_DS_UF_timechange = MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;

   VEC_DS_PF = MAT_D_velocityUnsteadyPlusDiffusion->create_vector( this ) ;

   // Local vector
   UF_DS_LOC = LA_SeqVector::create( this, 0 ) ;
   PF_DS_LOC = LA_SeqVector::create( this, 0 ) ;

   // velocity numbering
   UF_NUM = FV_SystemNumbering::create( this, UF ) ;
   PF_NUM = FV_SystemNumbering::create( this, PF ) ;

   // Direction splitting matrices
   MAT_velocityUnsteadyPlusDiffusion_1D = LA_SeqMatrix::make( this, exp->create_subexplorer( this,"MAT_1DLAP_generic" ) );

   // Structure for the particle surface discretization
   surface.coordinate = (LA_SeqMatrix*) malloc(sizeof(LA_SeqMatrix)) ;
   surface.area = (LA_SeqVector*) malloc(sizeof(LA_SeqVector)) ;
   surface.normal = (LA_SeqMatrix*) malloc(sizeof(LA_SeqMatrix)) ;

   for (size_t field = 0; field < 2; field++) {

      // Structure for the particle input data
      solid[field].coord = (LA_SeqMatrix**) malloc(nb_comps[field] * sizeof(LA_SeqMatrix*)) ;
      solid[field].size = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
      solid[field].thetap = (LA_SeqMatrix**) malloc(nb_comps[field] * sizeof(LA_SeqMatrix*)) ;
      solid[field].vel = (LA_SeqMatrix**) malloc(nb_comps[field] * sizeof(LA_SeqMatrix*)) ;
      solid[field].ang_vel = (LA_SeqMatrix**) malloc(nb_comps[field] * sizeof(LA_SeqMatrix*)) ;
      solid[field].temp = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
      solid[field].inside = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;

      // Vector to store the presence/absence of particle on the field variable
      node[field].void_frac = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
      node[field].parID = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
      node[field].bound_cell = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
      

      for (size_t dir = 0; dir < dim; dir++) {
         // Spacial discretization matrices
         A[field][dir].ii_main = (LA_SeqVector***) malloc(nb_comps[field] * sizeof(LA_SeqVector**)) ;
         A[field][dir].ii_super = (LA_SeqVector***) malloc(nb_comps[field] * sizeof(LA_SeqVector**)) ;
         A[field][dir].ii_sub = (LA_SeqVector***) malloc(nb_comps[field] * sizeof(LA_SeqVector**)) ;
         A[field][dir].ie = (LA_SeqMatrix***) malloc(nb_comps[field] * sizeof(LA_SeqMatrix**)) ;
         A[field][dir].ei = (LA_SeqMatrix***) malloc(nb_comps[field] * sizeof(LA_SeqMatrix**)) ;
         A[field][dir].ee = (LA_SeqMatrix***) malloc(nb_comps[field] * sizeof(LA_SeqMatrix**)) ;

         // Product matrices of spacial discretization
         Ap[field][dir].ei_ii_ie = (LA_SeqMatrix**) malloc(nb_comps[field] * sizeof(LA_SeqMatrix*)) ;
         Ap[field][dir].result = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
         Ap[field][dir].ii_ie = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;

         Ap_proc0[field][dir].ei_ii_ie = (LA_SeqMatrix**) malloc(nb_comps[field] * sizeof(LA_SeqMatrix*)) ;

         // VEC to store local/interface solution and RHS
         VEC[field][dir].local_T = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
         VEC[field][dir].local_solution_T = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
         VEC[field][dir].T = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
         VEC[field][dir].interface_T = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;

         // Schur complement matrices
         Schur[field][dir].ii_main = (LA_SeqVector***) malloc(nb_comps[field] * sizeof(LA_SeqVector**)) ;
         Schur[field][dir].ii_super = (LA_SeqVector***) malloc(nb_comps[field] * sizeof(LA_SeqVector**)) ;
         Schur[field][dir].ii_sub = (LA_SeqVector***) malloc(nb_comps[field] * sizeof(LA_SeqVector**)) ;
         Schur[field][dir].ie = (LA_SeqMatrix***) malloc(nb_comps[field] * sizeof(LA_SeqMatrix**)) ;
         Schur[field][dir].ei = (LA_SeqMatrix***) malloc(nb_comps[field] * sizeof(LA_SeqMatrix**)) ;
         Schur[field][dir].ee = (LA_SeqMatrix***) malloc(nb_comps[field] * sizeof(LA_SeqMatrix**)) ;

         // Matrix for Schur complement of Schur complement
         DoubleSchur[field][dir].ii_main = (LA_SeqVector***) malloc(nb_comps[field] * sizeof(LA_SeqVector**)) ;

         // Product of Schur complement matrices
         SchurP[field][dir].ei_ii_ie = (LA_SeqMatrix**) malloc(nb_comps[field] * sizeof(LA_SeqMatrix*)) ;
         SchurP[field][dir].result = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
         SchurP[field][dir].ii_ie = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;

         // VEC to store local/interface solution and RHS for Schur complement
         Schur_VEC[field][dir].local_T = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
         Schur_VEC[field][dir].local_solution_T = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;
         Schur_VEC[field][dir].T = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ; 
         Schur_VEC[field][dir].interface_T = (LA_SeqVector**) malloc(nb_comps[field] * sizeof(LA_SeqVector*)) ;

         for (size_t j=0;j<2;j++) {
            b_intersect[field][j][dir].offset = (LA_SeqMatrix**) malloc(nb_comps[field] * sizeof(LA_SeqMatrix*)) ;
            b_intersect[field][j][dir].value = (LA_SeqMatrix**) malloc(nb_comps[field] * sizeof(LA_SeqMatrix*)) ;
            b_intersect[field][j][dir].field_var = (LA_SeqMatrix**) malloc(nb_comps[field] * sizeof(LA_SeqMatrix*)) ;
         }

         for (size_t comp = 0; comp < nb_comps[field]; comp++) {
            size_t_vector nb_unknowns_handled_by_proc( dim, 0 );
            size_t nb_index=0;
            for (size_t l=0;l<dim;++l) {
               if (field == 0) {
                  nb_unknowns_handled_by_proc( l ) =
                                  1 + PF->get_max_index_unknown_handled_by_proc( comp, l )
                                    - PF->get_min_index_unknown_handled_by_proc( comp, l ) ;
               } else if (field == 1) {
                  nb_unknowns_handled_by_proc( l ) =
                                  1 + UF->get_max_index_unknown_handled_by_proc( comp, l )
                                    - UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
               }

            }
            if (dir == 0) {
               if (dim == 2) {
                  nb_index = nb_unknowns_handled_by_proc(1);
               } else if (dim == 3) {
                  nb_index = nb_unknowns_handled_by_proc(1)*nb_unknowns_handled_by_proc(2);
               }
            } else if (dir == 1) {
               if (dim == 2) {
                  nb_index = nb_unknowns_handled_by_proc(0);
               } else if (dim == 3) {
                  nb_index = nb_unknowns_handled_by_proc(0)*nb_unknowns_handled_by_proc(2);
               }
            } else if (dir == 2) {
               nb_index = nb_unknowns_handled_by_proc(0)*nb_unknowns_handled_by_proc(1);
            }

            A[field][dir].ii_main[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
            A[field][dir].ii_super[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
            A[field][dir].ii_sub[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
            A[field][dir].ie[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;
            A[field][dir].ei[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;
            A[field][dir].ee[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;

            Schur[field][dir].ii_main[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
            Schur[field][dir].ii_super[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
            Schur[field][dir].ii_sub[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
            Schur[field][dir].ie[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;
            Schur[field][dir].ei[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;
            Schur[field][dir].ee[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;

            DoubleSchur[field][dir].ii_main[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;

            for (size_t index = 0; index < nb_index; index++) {
               A[field][dir].ii_main[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               A[field][dir].ii_super[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               A[field][dir].ii_sub[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               A[field][dir].ie[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
               A[field][dir].ei[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );

               if (proc_pos_in_i[dir] == 0) {
                  A[field][dir].ee[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
                  Schur[field][dir].ii_main[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
                  Schur[field][dir].ii_super[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
                  Schur[field][dir].ii_sub[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
                  Schur[field][dir].ie[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
                  Schur[field][dir].ei[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
                  Schur[field][dir].ee[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );

                  DoubleSchur[field][dir].ii_main[comp][index] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               }
            }
         }
      }
   }

   surface.coordinate = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
   surface.area = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
   surface.normal = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );

   for (size_t field = 0; field < 2; field++) {
      for (size_t dir = 0; dir < dim; dir++) {
         for (size_t comp=0;comp<nb_comps[field];++comp) {
            Ap[field][dir].ei_ii_ie[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
            Ap[field][dir].result[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            Ap[field][dir].ii_ie[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;

            Ap_proc0[field][dir].ei_ii_ie[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );

            VEC[field][dir].local_T[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            VEC[field][dir].local_solution_T[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            VEC[field][dir].T[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            VEC[field][dir].interface_T[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;

            if (dir == 0) {
               solid[field].coord[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
               solid[field].size[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               solid[field].thetap[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
               solid[field].temp[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               solid[field].vel[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
               solid[field].ang_vel[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
               solid[field].inside[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               node[field].void_frac[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               node[field].parID[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               node[field].bound_cell[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            }

            for (size_t j=0;j<2;j++) {
               b_intersect[field][j][dir].offset[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
               b_intersect[field][j][dir].value[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
               b_intersect[field][j][dir].field_var[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
            }

            if (proc_pos_in_i[dir] == 0) {
               SchurP[field][dir].ei_ii_ie[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_copy( this,MAT_velocityUnsteadyPlusDiffusion_1D );
               SchurP[field][dir].result[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               SchurP[field][dir].ii_ie[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;

               Schur_VEC[field][dir].local_T[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               Schur_VEC[field][dir].local_solution_T[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               Schur_VEC[field][dir].T[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               Schur_VEC[field][dir].interface_T[comp] = MAT_velocityUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            }         
         }
      }
   }
}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem:: re_initialize( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: re_initialize" ) ;

   size_t UF_glob = UF->nb_global_unknowns() ;
   size_t UF_loc = UF->nb_local_unknowns() ;

   size_t pf_glob = PF->nb_global_unknowns() ;
   size_t pf_loc = PF->nb_local_unknowns() ;

   // velocity Laplacian
   MAT_D_velocityUnsteadyPlusDiffusion->re_initialize( UF_glob, UF_glob ) ;
   VEC_rhs_D_velocityDiffusionPlusBodyTerm->re_initialize( UF_glob ) ;

   // Unknowns vectors
   VEC_DS_UF->re_initialize( UF_glob ) ;
   VEC_DS_UF_previoustime->re_initialize( UF_glob ) ;
   VEC_DS_UF_timechange->re_initialize( UF_glob ) ;

   VEC_DS_PF->re_initialize( pf_glob ) ;

   // Local vector
   UF_DS_LOC->re_initialize( UF_loc ) ;

   PF_DS_LOC->re_initialize( pf_loc ) ;

   // velocity numbering
   UF_NUM->define_scatter( VEC_DS_UF ) ;
   PF_NUM->define_scatter( VEC_DS_PF ) ;

   // Initialize Direction splitting matrices & vectors for pressure 
   size_t nb_procs, proc_pos;

   if (is_solids && is_stressCal) {
      if (dim == 3) {
         if (level_set_type == "Sphere") {
            surface.coordinate->re_initialize(2*(size_t)Nmax,3);
            surface.area->re_initialize(2*(size_t)Nmax);
            surface.normal->re_initialize(2*(size_t)Nmax,3);
	 } else if (level_set_type == "Cube") {
            surface.coordinate->re_initialize(6*(size_t)pow(Nmax,2),3);
            surface.area->re_initialize(6*(size_t)pow(Nmax,2));
            surface.normal->re_initialize(6*(size_t)pow(Nmax,2),3);
	 } else if (level_set_type == "Cylinder") {
            double Npm1 = round(pow(MAC::sqrt(Nmax) - MAC::sqrt(MAC::pi()/ar),2.));
            double dh = 1. - MAC::sqrt(Npm1/Nmax);
            double Nr = round(2./dh);
	    size_t Ncyl = (size_t) (2*Nmax + Nr*(Nmax - Npm1));
            surface.coordinate->re_initialize(Ncyl,3);
            surface.area->re_initialize(Ncyl);
            surface.normal->re_initialize(Ncyl,3);
	 }
      } else {
	 if (level_set_type == "Sphere") {
            surface.coordinate->re_initialize((size_t)Nmax,3);
            surface.area->re_initialize((size_t)Nmax);
            surface.normal->re_initialize((size_t)Nmax,3);
	 } else if (level_set_type == "Cube") {
            surface.coordinate->re_initialize(4*(size_t)Nmax,3);
            surface.area->re_initialize(4*(size_t)Nmax);
            surface.normal->re_initialize(4*(size_t)Nmax,3);
	 }
      }
   }


   for (size_t field = 0; field < 2; field++) {
      for (size_t comp = 0; comp < nb_comps[field]; comp++) {
         size_t_vector nb_unknowns_handled_by_proc( dim, 0 );
         size_t_vector nb_unknowns_on_proc( dim, 0 );

         size_t nb_total_unknown = 1;

         for (size_t l = 0;l < dim; l++) { 
            if (field == 0) {
               nb_unknowns_handled_by_proc( l ) = 1 + PF->get_max_index_unknown_handled_by_proc( comp, l )
                                                    - PF->get_min_index_unknown_handled_by_proc( comp, l ) ;
               nb_unknowns_on_proc( l ) = 1 + PF->get_max_index_unknown_on_proc( comp, l )
                                            - PF->get_min_index_unknown_on_proc( comp, l ) ;
            } else if (field == 1) {
               nb_unknowns_handled_by_proc( l ) = 1 + UF->get_max_index_unknown_handled_by_proc( comp, l )
                                                    - UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
               nb_unknowns_on_proc( l ) = 1 + UF->get_max_index_unknown_on_proc( comp, l )
                                            - UF->get_min_index_unknown_on_proc( comp, l ) ;
            }
         }

         for (size_t l = 0;l < dim; l++) {
            nb_total_unknown *= (2+nb_unknowns_on_proc(l));
            size_t nb_index=0;
            if (l == 0) {
               if (dim == 2) {
                  nb_index = nb_unknowns_handled_by_proc(1);
               } else if (dim == 3) {
                  nb_index = nb_unknowns_handled_by_proc(1)*nb_unknowns_handled_by_proc(2);
               }
            } else if (l == 1) {
               if (dim == 2) {
                  nb_index = nb_unknowns_handled_by_proc(0);
               } else if (dim == 3) {
                  nb_index = nb_unknowns_handled_by_proc(0)*nb_unknowns_handled_by_proc(2);
               }
            } else if (l == 2) {
               nb_index = nb_unknowns_handled_by_proc(0)*nb_unknowns_handled_by_proc(1);
            }
             
            nb_procs = nb_procs_in_i[l];
            proc_pos = proc_pos_in_i[l];

            if ((l == 0) && (is_solids)) {
               // Presence of solid and only once
               solid[field].coord[comp]->re_initialize(Npart,3);
               solid[field].size[comp]->re_initialize(Npart);
               solid[field].thetap[comp]->re_initialize(Npart,3);
               solid[field].vel[comp]->re_initialize(Npart,3);
               solid[field].ang_vel[comp]->re_initialize(Npart,3);
               solid[field].temp[comp]->re_initialize(Npart);
               solid[field].inside[comp]->re_initialize(Npart);
            }

            if (is_periodic[field][l] != 1) {
               if (proc_pos == nb_procs-1) {
                  // Non-periodic and last processor
                  for (size_t index = 0; index < nb_index; index++) {
                     A[field][l].ii_main[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l ));
                     A[field][l].ii_super[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1);
                     A[field][l].ii_sub[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1);
                     A[field][l].ie[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l ),nb_procs-1);
                     A[field][l].ei[comp][index]->re_initialize(nb_procs-1,nb_unknowns_handled_by_proc( l ) );
                  }

                  Ap[field][l].result[comp]->re_initialize( nb_unknowns_handled_by_proc( l ) );
                  VEC[field][l].local_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )) ;
                  VEC[field][l].local_solution_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )) ;

               } else {
                  // Non-periodic for processor expect last
                  for (size_t index = 0; index < nb_index; index++) {
                     A[field][l].ii_main[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1);
                     A[field][l].ii_super[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-2);
                     A[field][l].ii_sub[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-2);
                     A[field][l].ie[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1,nb_procs-1 );
                     A[field][l].ei[comp][index]->re_initialize(nb_procs-1,nb_unknowns_handled_by_proc( l )-1 );
                  }

                  Ap[field][l].result[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1 );
                  VEC[field][l].local_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1) ;
                  VEC[field][l].local_solution_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1) ;
               }

               if (l == 1) MAT_velocityUnsteadyPlusDiffusion_1D->re_initialize(nb_unknowns_handled_by_proc( l ),nb_unknowns_handled_by_proc( l ) );
               Ap[field][l].ii_ie[comp]->re_initialize(nb_procs-1);
               Ap[field][l].ei_ii_ie[comp]->re_initialize(nb_procs-1,nb_procs-1 );
               Ap_proc0[field][l].ei_ii_ie[comp]->re_initialize(nb_procs-1,nb_procs-1 );
               VEC[field][l].interface_T[comp]->re_initialize( nb_procs-1 ) ;
               VEC[field][l].T[comp]->re_initialize( nb_procs-1 ) ;

               if (proc_pos == 0) {
                  // Master processor
                  for (size_t index = 0; index < nb_index; index++) {
                     A[field][l].ee[comp][index]->re_initialize(nb_procs-1,nb_procs-1 );
                     if (nb_procs != 1) {
                        Schur[field][l].ii_main[comp][index]->re_initialize(nb_procs-1);
                        Schur[field][l].ii_super[comp][index]->re_initialize(nb_procs-2);
                        Schur[field][l].ii_sub[comp][index]->re_initialize(nb_procs-2);
                     }
                  }
               }
   
            } else {
               // Periodic domain
               for (size_t index = 0; index < nb_index; index++) {
                  A[field][l].ii_main[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1);
                  A[field][l].ii_super[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-2);
                  A[field][l].ii_sub[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-2);
                  A[field][l].ie[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1,nb_procs );
                  A[field][l].ei[comp][index]->re_initialize(nb_procs,nb_unknowns_handled_by_proc( l )-1 );
               }

               Ap[field][l].result[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1 );
               Ap[field][l].ii_ie[comp]->re_initialize(nb_procs);
               Ap[field][l].ei_ii_ie[comp]->re_initialize(nb_procs,nb_procs );
               Ap_proc0[field][l].ei_ii_ie[comp]->re_initialize(nb_procs,nb_procs );

               VEC[field][l].local_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1) ;
               VEC[field][l].local_solution_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1) ;
               VEC[field][l].interface_T[comp]->re_initialize( nb_procs ) ;
               VEC[field][l].T[comp]->re_initialize( nb_procs ) ;

               if (proc_pos == 0) {
                  // Master processor
                  for (size_t index = 0; index < nb_index; index++) {
                     A[field][l].ee[comp][index]->re_initialize(nb_procs,nb_procs );
                  }
                  if (nb_procs != 1) {
                     // Mutli processor with periodic domain
                     // Condition where schur complement won't be a standard tridiagonal matrix but a variation
                     for (size_t index = 0; index < nb_index; index++) {
                        Schur[field][l].ii_main[comp][index]->re_initialize(nb_procs-1);
                        Schur[field][l].ii_super[comp][index]->re_initialize(nb_procs-2);
                        Schur[field][l].ii_sub[comp][index]->re_initialize(nb_procs-2);
                        Schur[field][l].ie[comp][index]->re_initialize(nb_procs-1,1);
                        Schur[field][l].ei[comp][index]->re_initialize(1,nb_procs-1);
                        Schur[field][l].ee[comp][index]->re_initialize(1,1);
                        DoubleSchur[field][l].ii_main[comp][index]->re_initialize(1);
                     }
 
                     SchurP[field][l].result[comp]->re_initialize(nb_procs-1);
                     SchurP[field][l].ii_ie[comp]->re_initialize(1);
                     SchurP[field][l].ei_ii_ie[comp]->re_initialize(1,1);

                     Schur_VEC[field][l].local_T[comp]->re_initialize(nb_procs-1) ;
                     Schur_VEC[field][l].local_solution_T[comp]->re_initialize(nb_procs-1) ;
                     Schur_VEC[field][l].interface_T[comp]->re_initialize(1) ;
                     Schur_VEC[field][l].T[comp]->re_initialize(1) ;
                  } else {
                     // Serial mode with periodic domain
                     for (size_t index = 0; index < nb_index; index++) {
                        Schur[field][l].ii_main[comp][index]->re_initialize(nb_procs);
                        Schur[field][l].ii_super[comp][index]->re_initialize(nb_procs-1);
                        Schur[field][l].ii_sub[comp][index]->re_initialize(nb_procs-1);
                     }
                  }
               }
            }
         }
         if (is_solids) {
            node[field].void_frac[comp]->re_initialize( nb_total_unknown ) ;
            node[field].parID[comp]->re_initialize( nb_total_unknown ) ;
            node[field].bound_cell[comp]->re_initialize( nb_total_unknown ) ;
            for (size_t i=0;i<dim;i++) {
               for (size_t j=0;j<2;j++) {
                  b_intersect[field][j][i].offset[comp]->re_initialize( nb_total_unknown,2 ) ;      // Column0 for left and Column1 for right
                  b_intersect[field][j][i].value[comp]->re_initialize( nb_total_unknown,2 ) ;      // Column0 for left and Column1 for right
                  b_intersect[field][j][i].field_var[comp]->re_initialize( nb_total_unknown,2 ) ;      // Column0 for left and Column1 for right
               }
            }
         }
      }
   }
}

//----------------------------------------------------------------------
DDS_NSWithHeatTransferSystem:: ~DDS_NSWithHeatTransferSystem( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: ~DDS_NSWithHeatTransferSystem" ) ;
}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::initialize_DS_velocity( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: initialize_DS_velocity" ) ;

   UF->extract_unknown_DOFs_value( 0, UF_DS_LOC ) ;
   UF_NUM->scatter()->set( UF_DS_LOC, VEC_DS_UF ) ;
         
}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::initialize_DS_pressure( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: initialize_DS_pressure" ) ;
   
   PF->extract_unknown_DOFs_value( 0, PF_DS_LOC ) ;
   PF_NUM->scatter()->set( PF_DS_LOC, VEC_DS_PF ) ;
         
}

//----------------------------------------------------------------------
LA_SeqVector const*
DDS_NSWithHeatTransferSystem:: get_solution_DS_velocity( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_solution_DS_velocity" ) ;

   UF_NUM->scatter()->get( VEC_DS_UF, UF_DS_LOC ) ;

   LA_SeqVector const* DS_result = UF_DS_LOC ;

   return( DS_result ) ;

}

//----------------------------------------------------------------------
LA_SeqVector const*
DDS_NSWithHeatTransferSystem:: get_solution_DS_pressure( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_solution_DS_pressure" ) ;

   PF_NUM->scatter()->get( VEC_DS_PF, PF_DS_LOC ) ;

   LA_SeqVector const* DS_result = PF_DS_LOC ;

   return( DS_result ) ;

}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::at_each_time_step( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: at_each_time_step" ) ;

   // Store velocity at previous time
   VEC_DS_UF->synchronize() ;
   VEC_DS_UF_previoustime->set( VEC_DS_UF ) ;
   VEC_DS_PF->synchronize() ;

}

//----------------------------------------------------------------------
double
DDS_NSWithHeatTransferSystem:: compute_DS_velocity_change( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: compute_DS_velocity_change" ) ;

   VEC_DS_UF->synchronize() ;
   VEC_DS_UF_timechange->set( VEC_DS_UF ) ;
   VEC_DS_UF_timechange->sum( VEC_DS_UF_previoustime, -1.0 ) ;

   double norm_UF = VEC_DS_UF->two_norm() ;
   double time_change = VEC_DS_UF_timechange->two_norm() ;
   if ( norm_UF > 1e-4 ) time_change /= norm_UF;

   return ( time_change ) ;
}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::pre_thomas_treatment( size_t const& comp, size_t const& dir, struct TDMatrix *arr, size_t const& r_index)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: pre_thomas_treatment" ) ;

   size_t nrows = arr[dir].ii_main[comp][r_index]->nb_rows() ;

   double temp = arr[dir].ii_main[comp][r_index]->item(0);
   if (nrows > 1) arr[dir].ii_super[comp][r_index]->set_item(0,arr[dir].ii_super[comp][r_index]->item(0)/temp);

   //  // Perform Forward Elimination
   size_t m;

   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,c,prevc;
     a = arr[dir].ii_sub[comp][r_index]->item(m-1);
     b = arr[dir].ii_main[comp][r_index]->item(m);
     prevc = arr[dir].ii_super[comp][r_index]->item(m-1);

     if(m<nrows-1){
         c = arr[dir].ii_super[comp][r_index]->item(m);
         arr[dir].ii_super[comp][r_index]->set_item(m,c/(b - a*prevc));
     }
   }
}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::mod_thomas_algorithm(TDMatrix *arr, LA_SeqVector* rhs, size_t const& comp, size_t const& dir, size_t const& r_index)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquationSystem:: mod_thomas_algorithm" ) ;

   size_t nrows = arr[dir].ii_main[comp][r_index] -> nb_rows() ;
   double temp = arr[dir].ii_main[comp][r_index]->item(0);
   rhs-> set_item(0,rhs ->item(0)/temp);

    // Perform Forward Elimination
   size_t m;
   /// Showing problem for last row elimination
   for (m=1;m<nrows;++m)
   {
     double a,b,d,prevd,prevc;
     a=arr[dir].ii_sub[comp][r_index]->item(m-1);
     b=arr[dir].ii_main[comp][r_index]->item(m);
     d=rhs->item(m);
     prevc=arr[dir].ii_super[comp][r_index]->item(m-1);
     prevd=rhs->item(m-1);

     rhs -> set_item(m, (d-a*prevd)/(b-a*prevc));
   }

   //Perform backward substitution
   if(nrows>1){
      rhs->set_item(nrows-1,rhs->item(nrows-1));
      for (m = nrows-2; m< nrows-1;m--) {
         double c,nextd;
         c=arr[dir].ii_super[comp][r_index]->item(m);
         nextd=rhs->item(m+1);
         rhs->add_to_item(m,-c*nextd);
      }
   }
}

//----------------------------------------------------------------------
TDMatrix*
DDS_NSWithHeatTransferSystem::get_A(size_t const& field)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_A" ) ;
   return (A[field]) ;
}

//----------------------------------------------------------------------
PartInput
DDS_NSWithHeatTransferSystem::get_solid(size_t const& field)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_solid" ) ;
   return (solid[field]) ;
}

//----------------------------------------------------------------------
SurfaceDiscretize
DDS_NSWithHeatTransferSystem::get_surface()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_surface" ) ;
   return (surface) ;
}

//----------------------------------------------------------------------
BoundaryBisec*
DDS_NSWithHeatTransferSystem::get_b_intersect(size_t const& field, size_t const& level)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_b_intersect" ) ;
   return (b_intersect[field][level]) ;
}

//----------------------------------------------------------------------
NodeProp
DDS_NSWithHeatTransferSystem::get_node_property(size_t const& field)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquationSystem:: get_node_property" ) ;
   return (node[field]) ;
}

//----------------------------------------------------------------------
TDMatrix*
DDS_NSWithHeatTransferSystem::get_Schur(size_t const& field)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_Schur" ) ;
   return (Schur[field]) ;
}

//----------------------------------------------------------------------
TDMatrix*
DDS_NSWithHeatTransferSystem::get_DoubleSchur(size_t const& field)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_DoubleSchur" ) ;
   return (DoubleSchur[field]) ;
}

//----------------------------------------------------------------------
ProdMatrix*
DDS_NSWithHeatTransferSystem::get_Ap(size_t const& field)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_Ap" ) ;
   return (Ap[field]) ;
}

//----------------------------------------------------------------------
ProdMatrix*
DDS_NSWithHeatTransferSystem::get_Ap_proc0(size_t const& field)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_Ap" ) ;
   return (Ap_proc0[field]) ;
}

//----------------------------------------------------------------------
ProdMatrix*
DDS_NSWithHeatTransferSystem::get_SchurP(size_t const& field)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_SchurP" ) ;
   return (SchurP[field]) ;
}

//----------------------------------------------------------------------
LocalVector*
DDS_NSWithHeatTransferSystem::get_Schur_VEC(size_t const& field)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_Schur_VEC" ) ;
   return (Schur_VEC[field]) ;
}

//----------------------------------------------------------------------
LocalVector*
DDS_NSWithHeatTransferSystem::get_VEC(size_t const& field)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: get_VEC" ) ;
   return (VEC[field]) ;
}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::update_global_P_vector(size_t const& i, size_t const& j, size_t const& k, double const& value)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: update_global_P_vector" ) ;

   size_t global_number_in_distributed_vector = PF->DOF_global_number(i,j,k,0);

   VEC_DS_PF->set_item( global_number_in_distributed_vector, value);

}
//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::DS_NavierStokes_solver(FV_DiscreteField* FF
	,size_t const& j, size_t const& k, size_t const& min_i, size_t const& comp, size_t const& dir, size_t const& field, size_t const& r_index )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: DS_NavierStokes_solver" ) ;

   LocalVector* rhs = get_VEC(field);
   TDMatrix* arr = get_A(field);

   size_t nb_procs, proc_pos;

   proc_pos = proc_pos_in_i[dir];
   nb_procs = nb_procs_in_i[dir];

   // Solve the DS splitting problem in
   DDS_NSWithHeatTransferSystem::mod_thomas_algorithm( arr, rhs[dir].local_T[comp], comp, dir, r_index);

   // Transfer in the distributed vector
   size_t nb_local_unk = rhs[dir].local_T[comp]->nb_rows();
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
         VEC_DS_PF->set_item( global_number_in_distributed_vector, rhs[dir].local_T[comp]->item( m ) );
      } else if (field == 1) { 
         VEC_DS_UF->set_item( global_number_in_distributed_vector, rhs[dir].local_T[comp]->item( m ) );
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
         VEC_DS_PF->set_item( global_number_in_distributed_vector,rhs[dir].interface_T[comp]->item( proc_pos ) );
      } else if (field == 1) {
         VEC_DS_UF->set_item( global_number_in_distributed_vector, rhs[dir].interface_T[comp]->item( proc_pos ) );
      }
   } else if ((is_periodic[field][dir] == 0) && (proc_pos != nb_procs-1)) {
      global_number_in_distributed_vector = FF->DOF_global_number( ii, jj, kk, comp );
      if (field == 0) {
         VEC_DS_PF->set_item( global_number_in_distributed_vector,rhs[dir].interface_T[comp]->item( proc_pos ) );
      } else if (field == 1) {
         VEC_DS_UF->set_item( global_number_in_distributed_vector, rhs[dir].interface_T[comp]->item( proc_pos ) );
      }
   }
}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::synchronize_DS_solution_vec( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: synchronize_DS_solution_vec" ) ;

   VEC_DS_UF->synchronize();
}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::synchronize_DS_solution_vec_P( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: synchronize_DS_solution_vec" ) ;

   VEC_DS_PF->synchronize();
}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::compute_product_matrix_interior(struct TDMatrix *arr, struct ProdMatrix *prr, size_t const& comp, size_t const& column,size_t const& dir,size_t const& r_index)
//----------------------------------------------------------------------
{

  MAC_LABEL( "DDS_NSWithHeatTransferSystem:: compute_product_matrix_interior" ) ;

  // Get appropriate column of Aie
  arr[dir].ie[comp][r_index]->extract_col(column, prr[dir].result[comp]);

  // Get inv(Aii)*Aie for for appropriate column of Aie
  mod_thomas_algorithm(arr, prr[dir].result[comp], comp, dir,r_index);

  // Get product of Aei*inv(Aii)*Aie for appropriate column
  arr[dir].ei[comp][r_index]->multiply_vec_then_add(prr[dir].result[comp],prr[dir].ii_ie[comp]);

  size_t int_unknown = prr[dir].ii_ie[comp]->nb_rows();

  for (size_t i = 0; i < int_unknown; i++){
      prr[dir].ei_ii_ie[comp]->set_item(i,column,prr[dir].ii_ie[comp]->item(i));
  }

}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::compute_product_matrix(struct TDMatrix *arr, struct ProdMatrix *prr, size_t const& comp, size_t const& dir, size_t const& field, size_t const& r_index )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransferSystem:: compute_product_matrix" ) ;

   size_t proc_pos, nb_procs;

   proc_pos = proc_pos_in_i[dir];
   nb_procs = nb_procs_in_i[dir];

   if (proc_pos == nb_procs - 1){
      // Condition for serial processor and multi processor
      if (proc_pos == 0) {
         compute_product_matrix_interior(arr,prr,comp,proc_pos,dir,r_index);
      } else {
         compute_product_matrix_interior(arr,prr,comp,proc_pos-1,dir,r_index);
         if (is_periodic[field][dir] == 1) compute_product_matrix_interior(arr,prr,comp,proc_pos,dir,r_index);
      }
   }else if(proc_pos == 0){
      compute_product_matrix_interior(arr,prr,comp,proc_pos,dir,r_index);
      if (is_periodic[field][dir] == 1) compute_product_matrix_interior(arr,prr,comp,nb_procs-1,dir,r_index);
   }else{
      compute_product_matrix_interior(arr,prr,comp,proc_pos-1,dir,r_index);
      compute_product_matrix_interior(arr,prr,comp,proc_pos,dir,r_index);
   }
}

//----------------------------------------------------------------------
void
DDS_NSWithHeatTransferSystem::display_debug(void)
//----------------------------------------------------------------------
{
  // VEC_DS_UF->print_items(MAC::out(),0);

}
