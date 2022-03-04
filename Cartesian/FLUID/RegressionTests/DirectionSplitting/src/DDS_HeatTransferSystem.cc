#include <DDS_HeatTransferSystem.hh>
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
DDS_HeatTransferSystem*
DDS_HeatTransferSystem:: create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_tf,
        struct HeatTransfer2System const& transfer )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   MAC_CHECK_PRE( mac_tf != 0 ) ;

   DDS_HeatTransferSystem* result =
         new DDS_HeatTransferSystem( a_owner, exp, mac_tf, transfer ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;

   return( result ) ;

}




//----------------------------------------------------------------------
DDS_HeatTransferSystem:: DDS_HeatTransferSystem(
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_tf,
        struct HeatTransfer2System const& fromHE )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , TF( mac_tf )
   , MAT_TemperatureUnsteadyPlusDiffusion_1D( 0 )
   , is_solids ( fromHE.is_solids_ )
   , is_stressCal (fromHE.is_stressCal_ )
   , Npart (fromHE.Npart_ )
   , level_set_type (fromHE.level_set_type_ )
   , Nmax (fromHE.Npoints_ )
   , ar (fromHE.ar_ ) 
{
   MAC_LABEL( "DDS_HeatTransferSystem:: DDS_HeatTransferSystem" ) ;

   int const* MPI_coordinates_world = TF->primary_grid()->get_MPI_coordinates() ;
   int const* MPI_max_coordinates_world = TF->primary_grid()->get_domain_decomposition() ;

   is_iperiodic[0] = false;
   is_iperiodic[1] = false;
   is_iperiodic[2] = false;

   proc_pos_in_i[0] = MPI_coordinates_world[0];
   proc_pos_in_i[1] = MPI_coordinates_world[1];
   proc_pos_in_i[2] = MPI_coordinates_world[2];
   nb_procs_in_i[0] = MPI_max_coordinates_world[0];
   nb_procs_in_i[1] = MPI_max_coordinates_world[1];
   nb_procs_in_i[2] = MPI_max_coordinates_world[2];

   dim = TF->primary_grid()->nb_space_dimensions() ;
   nb_comps = TF->nb_components() ;

   periodic_comp = TF->primary_grid()->get_periodic_directions();
   is_iperiodic[0] = periodic_comp->operator()( 0 );
   is_iperiodic[1] = periodic_comp->operator()( 1 );
   if(dim >2)
      is_iperiodic[2] = periodic_comp->operator()( 2 );

   // Build the matrices & vectors
   build_system(exp) ;
   re_initialize() ;

}




//----------------------------------------------------------------------
void
DDS_HeatTransferSystem:: build_system( MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: build_system" ) ;

   // Temperature Laplacian
   MAT_D_TemperatureUnsteadyPlusDiffusion = LA_Matrix::make( this,
         exp->create_subexplorer( this,
	 "MAT_D_TemperatureDiffusion"  ) ) ;

   VEC_DS_TF = MAT_D_TemperatureUnsteadyPlusDiffusion->create_vector( this ) ;
   VEC_DS_TF_previoustime = MAT_D_TemperatureUnsteadyPlusDiffusion->create_vector( this ) ;
   VEC_DS_TF_timechange = MAT_D_TemperatureUnsteadyPlusDiffusion->create_vector( this ) ;

   // Local vector
   TF_DS_LOC = LA_SeqVector::create( this, 0 ) ;

   // Temperature numbering
   TF_NUM = FV_SystemNumbering::create( this, TF ) ;

   // Direction splitting matrices
   MAT_TemperatureUnsteadyPlusDiffusion_1D =
   	LA_SeqMatrix::make( this, exp->create_subexplorer( this,"MAT_1DLAP_generic" ) );

   // Structure for the particle surface discretization
   surface.coordinate = (LA_SeqMatrix*) malloc(sizeof(LA_SeqMatrix)) ;
   surface.area = (LA_SeqVector*) malloc(sizeof(LA_SeqVector)) ;
   surface.normal = (LA_SeqMatrix*) malloc(sizeof(LA_SeqMatrix)) ;

   // Structure for the particle input data
   solid.coord = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
   solid.size = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   solid.thetap = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
   solid.vel = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
   solid.ang_vel = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
   solid.temp = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   solid.inside = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

   // Vector to store the presence/absence of particle on the field variable
   node.void_frac = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
   node.parID = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

   for (size_t dir = 0; dir < dim; dir++) {
      // Spacial discretization matrices
      A[dir].ii_main = (LA_SeqVector***) malloc(nb_comps * sizeof(LA_SeqVector**)) ;
      A[dir].ii_super = (LA_SeqVector***) malloc(nb_comps * sizeof(LA_SeqVector**)) ;
      A[dir].ii_sub = (LA_SeqVector***) malloc(nb_comps * sizeof(LA_SeqVector**)) ;
      A[dir].ie = (LA_SeqMatrix***) malloc(nb_comps * sizeof(LA_SeqMatrix**)) ;
      A[dir].ei = (LA_SeqMatrix***) malloc(nb_comps * sizeof(LA_SeqMatrix**)) ;
      A[dir].ee = (LA_SeqMatrix***) malloc(nb_comps * sizeof(LA_SeqMatrix**)) ;

      // Schur complement matrices
      Schur[dir].ii_main = (LA_SeqVector***) malloc(nb_comps * sizeof(LA_SeqVector**)) ;
      Schur[dir].ii_super = (LA_SeqVector***) malloc(nb_comps * sizeof(LA_SeqVector**)) ;
      Schur[dir].ii_sub = (LA_SeqVector***) malloc(nb_comps * sizeof(LA_SeqVector**)) ;
      Schur[dir].ie = (LA_SeqMatrix***) malloc(nb_comps * sizeof(LA_SeqMatrix**)) ;
      Schur[dir].ei = (LA_SeqMatrix***) malloc(nb_comps * sizeof(LA_SeqMatrix**)) ;
      Schur[dir].ee = (LA_SeqMatrix***) malloc(nb_comps * sizeof(LA_SeqMatrix**)) ;

      // Matrix for Schur complement of Schur complement
      DoubleSchur[dir].ii_main = (LA_SeqVector***) malloc(nb_comps * sizeof(LA_SeqVector**)) ;


      for (size_t comp = 0; comp < nb_comps; comp++) {
         size_t_vector nb_unknowns_handled_by_proc( dim, 0 );
         size_t nb_index=0;
         for (size_t l=0;l<dim;++l) {
            nb_unknowns_handled_by_proc( l ) =
                            1 + TF->get_max_index_unknown_handled_by_proc( comp, l )
                              - TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
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

         A[dir].ii_main[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
         A[dir].ii_super[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
         A[dir].ii_sub[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
         A[dir].ie[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;
         A[dir].ei[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;
         A[dir].ee[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;
         
         Schur[dir].ii_main[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
         Schur[dir].ii_super[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
         Schur[dir].ii_sub[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;
         Schur[dir].ie[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;
         Schur[dir].ei[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;
         Schur[dir].ee[comp] = (LA_SeqMatrix**) malloc(nb_index * sizeof(LA_SeqMatrix*)) ;

         DoubleSchur[dir].ii_main[comp] = (LA_SeqVector**) malloc(nb_index * sizeof(LA_SeqVector*)) ;

         for (size_t index = 0; index < nb_index; index++) {
            A[dir].ii_main[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            A[dir].ii_super[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            A[dir].ii_sub[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            A[dir].ie[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
            A[dir].ei[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );

            if (proc_pos_in_i[dir] == 0) {
               A[dir].ee[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
               Schur[dir].ii_main[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               Schur[dir].ii_super[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               Schur[dir].ii_sub[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
               Schur[dir].ie[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
               Schur[dir].ei[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
               Schur[dir].ee[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );

               DoubleSchur[dir].ii_main[comp][index] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            }
         }
      }
         

      // Product matrices of spacial discretization
      Ap[dir].ei_ii_ie = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
      Ap[dir].result = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      Ap[dir].ii_ie = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

      Ap_proc0[dir].ei_ii_ie = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;

      // VEC to store local/interface solution and RHS 
      VEC[dir].local_T = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      VEC[dir].local_solution_T = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      VEC[dir].T = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      VEC[dir].interface_T = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

      // Product of Schur complement matrices
      SchurP[dir].ei_ii_ie = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
      SchurP[dir].result = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      SchurP[dir].ii_ie = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

      // VEC to store local/interface solution and RHS for Schur complement
      Schur_VEC[dir].local_T = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      Schur_VEC[dir].local_solution_T = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      Schur_VEC[dir].T = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;
      Schur_VEC[dir].interface_T = (LA_SeqVector**) malloc(nb_comps * sizeof(LA_SeqVector*)) ;

      for (size_t j=0;j<2;j++) {
         b_intersect[j][dir].offset = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
         b_intersect[j][dir].value = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
         b_intersect[j][dir].field = (LA_SeqMatrix**) malloc(nb_comps * sizeof(LA_SeqMatrix*)) ;
      }
   }

   surface.coordinate = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
   surface.area = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
   surface.normal = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );

   for (size_t dir=0;dir<dim;++dir) {
      for (size_t comp=0;comp<nb_comps;++comp) {
         Ap[dir].ei_ii_ie[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
         Ap[dir].result[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
         Ap[dir].ii_ie[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;

         Ap_proc0[dir].ei_ii_ie[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );

         VEC[dir].local_T[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
         VEC[dir].local_solution_T[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
         VEC[dir].T[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
         VEC[dir].interface_T[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;

         if (dir == 0) {
            solid.coord[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
            solid.size[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            solid.thetap[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
            solid.temp[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            solid.vel[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
            solid.ang_vel[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
            solid.inside[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            node.void_frac[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            node.parID[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
         }

         for (size_t j=0;j<2;j++) {
            b_intersect[j][dir].offset[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
            b_intersect[j][dir].value[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
            b_intersect[j][dir].field[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
         }

         if (proc_pos_in_i[dir] == 0) {
            SchurP[dir].ei_ii_ie[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_copy( this,MAT_TemperatureUnsteadyPlusDiffusion_1D );
            SchurP[dir].result[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            SchurP[dir].ii_ie[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;

            Schur_VEC[dir].local_T[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            Schur_VEC[dir].local_solution_T[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            Schur_VEC[dir].T[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
            Schur_VEC[dir].interface_T[comp] = MAT_TemperatureUnsteadyPlusDiffusion_1D->create_vector( this ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
DDS_HeatTransferSystem:: re_initialize( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: re_initialize" ) ;

   size_t tf_glob = TF->nb_global_unknowns() ;
   size_t tf_loc = TF->nb_local_unknowns() ;

   // Temperature Laplacian
   MAT_D_TemperatureUnsteadyPlusDiffusion->re_initialize( tf_glob, tf_glob ) ;

   // Unknowns vectors
   VEC_DS_TF->re_initialize( tf_glob ) ;
   VEC_DS_TF_previoustime->re_initialize( tf_glob ) ;
   VEC_DS_TF_timechange->re_initialize( tf_glob ) ;

   // Local vector
   TF_DS_LOC->re_initialize( tf_loc ) ;

   // Temperature numbering
   TF_NUM->define_scatter( VEC_DS_TF ) ;

   // Direction splitting matrices & vectors
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

   for (size_t comp=0;comp<nb_comps;++comp) {
      size_t_vector nb_unknowns_handled_by_proc( dim, 0 );
      size_t_vector nb_dof_on_proc( dim, 0 );

      size_t nb_total_unknown = 1;

      for (size_t l=0;l<dim;++l) {
         nb_unknowns_handled_by_proc( l ) =
                            1 + TF->get_max_index_unknown_handled_by_proc( comp, l )
                              - TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
	 nb_dof_on_proc( l ) = TF->get_local_nb_dof( comp, l ) ;
      }

    
      for (size_t l=0;l<dim;++l) {
         nb_total_unknown *= 1+nb_dof_on_proc(l);

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
            solid.coord[comp]->re_initialize(Npart,3);
            solid.size[comp]->re_initialize(Npart);
            solid.thetap[comp]->re_initialize(Npart,9);
            solid.vel[comp]->re_initialize(Npart,3);
            solid.ang_vel[comp]->re_initialize(Npart,3);
            solid.temp[comp]->re_initialize(Npart);
            solid.inside[comp]->re_initialize(Npart);
         }

         if (is_iperiodic[l] != 1) {
            if (proc_pos == nb_procs-1) {
	       // Non-periodic and last processor
               for (size_t index = 0; index < nb_index; index++) {
                  A[l].ii_main[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l ));
                  A[l].ii_super[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1);
                  A[l].ii_sub[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1);
                  A[l].ie[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l ),nb_procs-1);
                  A[l].ei[comp][index]->re_initialize(nb_procs-1,nb_unknowns_handled_by_proc( l ) );
               }

               Ap[l].result[comp]->re_initialize( nb_unknowns_handled_by_proc( l ) );
               VEC[l].local_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )) ;
               VEC[l].local_solution_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )) ;

            } else {
	       // Non-periodic for processor expect last
               for (size_t index = 0; index < nb_index; index++) {
                  A[l].ii_main[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1);
                  A[l].ii_super[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-2);
                  A[l].ii_sub[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-2);
                  A[l].ie[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1,nb_procs-1 );
                  A[l].ei[comp][index]->re_initialize(nb_procs-1,nb_unknowns_handled_by_proc( l )-1 );
               }

               Ap[l].result[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1 );
               VEC[l].local_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1) ;
               VEC[l].local_solution_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1) ;
            }

            if (l == 1) MAT_TemperatureUnsteadyPlusDiffusion_1D->re_initialize(nb_unknowns_handled_by_proc( l ),nb_unknowns_handled_by_proc( l ) );
            Ap[l].ii_ie[comp]->re_initialize(nb_procs-1);
            Ap[l].ei_ii_ie[comp]->re_initialize(nb_procs-1,nb_procs-1 );
            Ap_proc0[l].ei_ii_ie[comp]->re_initialize(nb_procs-1,nb_procs-1 );
            VEC[l].interface_T[comp]->re_initialize( nb_procs-1 ) ;
            VEC[l].T[comp]->re_initialize( nb_procs-1 ) ;

            if (proc_pos == 0) {
	       // Master processor
               for (size_t index = 0; index < nb_index; index++) {
                  A[l].ee[comp][index]->re_initialize(nb_procs-1,nb_procs-1 );
                  if (nb_procs != 1) {
                     Schur[l].ii_main[comp][index]->re_initialize(nb_procs-1);
                     Schur[l].ii_super[comp][index]->re_initialize(nb_procs-2);
                     Schur[l].ii_sub[comp][index]->re_initialize(nb_procs-2);
                  }
               }
            }

         } else {             
            // Periodic domain 
            for (size_t index = 0; index < nb_index; index++) {
               A[l].ii_main[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1);
               A[l].ii_super[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-2);
               A[l].ii_sub[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-2);
               A[l].ie[comp][index]->re_initialize(nb_unknowns_handled_by_proc( l )-1,nb_procs );
               A[l].ei[comp][index]->re_initialize(nb_procs,nb_unknowns_handled_by_proc( l )-1 );
            }

	    Ap[l].result[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1 );
            Ap[l].ii_ie[comp]->re_initialize(nb_procs);
            Ap[l].ei_ii_ie[comp]->re_initialize(nb_procs,nb_procs );

            Ap_proc0[l].ei_ii_ie[comp]->re_initialize(nb_procs,nb_procs );

	    VEC[l].local_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1) ;
            VEC[l].local_solution_T[comp]->re_initialize( nb_unknowns_handled_by_proc( l )-1) ;
            VEC[l].interface_T[comp]->re_initialize( nb_procs ) ;
            VEC[l].T[comp]->re_initialize( nb_procs ) ;

            if (proc_pos == 0) {
	       // Master processor
               for (size_t index = 0; index < nb_index; index++) {
                  A[l].ee[comp][index]->re_initialize(nb_procs,nb_procs );
               }
	       if (nb_procs != 1) {
		  // Mutli processor with periodic domain
		  // Condition where schur complement won't be a standard tridiagonal matrix but a variation
                  for (size_t index = 0; index < nb_index; index++) {
                     Schur[l].ii_main[comp][index]->re_initialize(nb_procs-1);
                     Schur[l].ii_super[comp][index]->re_initialize(nb_procs-2);
                     Schur[l].ii_sub[comp][index]->re_initialize(nb_procs-2);
                     Schur[l].ie[comp][index]->re_initialize(nb_procs-1,1);
                     Schur[l].ei[comp][index]->re_initialize(1,nb_procs-1);
                     Schur[l].ee[comp][index]->re_initialize(1,1);
		     DoubleSchur[l].ii_main[comp][index]->re_initialize(1);
                  }

		  SchurP[l].result[comp]->re_initialize(nb_procs-1);
                  SchurP[l].ii_ie[comp]->re_initialize(1);
                  SchurP[l].ei_ii_ie[comp]->re_initialize(1,1);

      		  Schur_VEC[l].local_T[comp]->re_initialize(nb_procs-1) ;
                  Schur_VEC[l].local_solution_T[comp]->re_initialize(nb_procs-1) ;
                  Schur_VEC[l].interface_T[comp]->re_initialize(1) ;
                  Schur_VEC[l].T[comp]->re_initialize(1) ;

	       } else {
                  for (size_t index = 0; index < nb_index; index++) {
                     // Serial mode with periodic domain
                     Schur[l].ii_main[comp][index]->re_initialize(nb_procs);
                     Schur[l].ii_super[comp][index]->re_initialize(nb_procs-1);
                     Schur[l].ii_sub[comp][index]->re_initialize(nb_procs-1);
                  }
	       }
            }
         }
      }
      if (is_solids) {
         node.void_frac[comp]->re_initialize( nb_total_unknown ) ;
         node.parID[comp]->re_initialize( nb_total_unknown ) ;
         for (size_t i=0;i<dim;i++) {
            for (size_t j=0;j<2;j++) {
               b_intersect[j][i].offset[comp]->re_initialize( nb_total_unknown,2 ) ;      // Column0 for left and Column1 for right
               b_intersect[j][i].value[comp]->re_initialize( nb_total_unknown,2 ) ;      // Column0 for left and Column1 for right
               b_intersect[j][i].field[comp]->re_initialize( nb_total_unknown,2 ) ;      // Column0 for left and Column1 for right
            }
         }
      }
   }
}

//----------------------------------------------------------------------
DDS_HeatTransferSystem:: ~DDS_HeatTransferSystem( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: ~DDS_HeatTransferSystem" ) ;
}

//----------------------------------------------------------------------
void
DDS_HeatTransferSystem::initialize_temperature( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: initialize_temperature" ) ;

   TF->extract_unknown_DOFs_value( 0, TF_DS_LOC ) ;
   TF_NUM->scatter()->set( TF_DS_LOC, VEC_DS_TF ) ;

}

//----------------------------------------------------------------------
LA_SeqVector const*
DDS_HeatTransferSystem:: get_solution_DS_temperature( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_solution_DS_temperature" ) ;

   TF_NUM->scatter()->get( VEC_DS_TF, TF_DS_LOC ) ;

   LA_SeqVector const* DS_result = TF_DS_LOC ;

   return( DS_result ) ;

}

//----------------------------------------------------------------------
void
DDS_HeatTransferSystem::at_each_time_step( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: at_each_time_step" ) ;

   // Store temperature at previous time
   VEC_DS_TF->synchronize() ;
   VEC_DS_TF_previoustime->set( VEC_DS_TF ) ;

}

//----------------------------------------------------------------------
double
DDS_HeatTransferSystem:: compute_temperature_change( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: compute_temperature_change" ) ;

   VEC_DS_TF->synchronize() ;
   VEC_DS_TF_timechange->set( VEC_DS_TF ) ;
   VEC_DS_TF_timechange->sum( VEC_DS_TF_previoustime, -1.0 ) ;

   double norm_TF = VEC_DS_TF->two_norm() ;
   double time_change = VEC_DS_TF_timechange->two_norm() ;
   if ( norm_TF > 1e-4 ) time_change /= norm_TF;

   return ( time_change ) ;

}

//----------------------------------------------------------------------
void
DDS_HeatTransferSystem::pre_thomas_treatment( size_t const& comp, size_t const& dir, struct TDMatrix *arr, size_t const& r_index)
//----------------------------------------------------------------------
{
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
DDS_HeatTransferSystem::mod_thomas_algorithm(TDMatrix *arr, LA_SeqVector* rhs, size_t const& comp, size_t const& dir, size_t const& r_index)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: mod_thomas_algorithm" ) ;

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
DDS_HeatTransferSystem::get_DoubleSchur()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_Schur" ) ;
   return (DoubleSchur) ;
}

//----------------------------------------------------------------------
TDMatrix*
DDS_HeatTransferSystem::get_Schur()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_Schur" ) ;
   return (Schur) ;
}

//----------------------------------------------------------------------
TDMatrix*
DDS_HeatTransferSystem::get_A()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_A" ) ;
   return (A) ;
}

//----------------------------------------------------------------------
ProdMatrix*
DDS_HeatTransferSystem::get_Ap()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_Ap" ) ;
   return (Ap) ;
}

//----------------------------------------------------------------------
ProdMatrix*
DDS_HeatTransferSystem::get_Ap_proc0()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_Ap" ) ;
   return (Ap_proc0) ;
}


//----------------------------------------------------------------------
PartInput
DDS_HeatTransferSystem::get_solid()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_solid" ) ;
   return (solid) ;
}

//----------------------------------------------------------------------
SurfaceDiscretize
DDS_HeatTransferSystem::get_surface()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_surface" ) ;
   return (surface) ;
}

//----------------------------------------------------------------------
ProdMatrix*
DDS_HeatTransferSystem::get_SchurP()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_SchurP" ) ;
   return (SchurP) ;
}

//----------------------------------------------------------------------
LocalVector*
DDS_HeatTransferSystem::get_VEC()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_VEC" ) ;
   return (VEC) ;
}

//----------------------------------------------------------------------
BoundaryBisec*
DDS_HeatTransferSystem::get_b_intersect(size_t const& level)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_b_intersect" ) ;
   return (b_intersect[level]) ;
}

//----------------------------------------------------------------------
NodeProp
DDS_HeatTransferSystem::get_node_property()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_voidfrac_VEC" ) ;
   return (node) ;
}

//----------------------------------------------------------------------
LocalVector*
DDS_HeatTransferSystem::get_Schur_VEC()
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: get_VEC" ) ;
   return (Schur_VEC) ;
}

//----------------------------------------------------------------------
void
DDS_HeatTransferSystem::DS_HeatEquation_solver(
        size_t const& j, size_t const& k, size_t const& min_i, size_t const& comp, size_t const& dir, size_t const& r_index)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: DS_HeatEquation_solver" ) ;

   LocalVector* rhs = get_VEC();
   TDMatrix* arr = get_A();

   size_t nb_procs, proc_pos;

   proc_pos = proc_pos_in_i[dir];
   nb_procs = nb_procs_in_i[dir];

   // Solve the DS splitting problem in
   DDS_HeatTransferSystem::mod_thomas_algorithm(arr, rhs[dir].local_T[comp], comp, dir,r_index);

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

      global_number_in_distributed_vector = TF->DOF_global_number( ii, jj, kk, comp );
      VEC_DS_TF->set_item( global_number_in_distributed_vector,rhs[dir].local_T[comp]->item( m ) );
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

   if ((is_iperiodic[dir] == 1)) {
      global_number_in_distributed_vector =
      TF->DOF_global_number( ii, jj, kk, comp );
      VEC_DS_TF->set_item( global_number_in_distributed_vector,
          rhs[dir].interface_T[comp]->item( proc_pos ) );
   } else if ((is_iperiodic[dir] == 0) && (proc_pos != nb_procs-1)) {
      global_number_in_distributed_vector =
      TF->DOF_global_number( ii, jj, kk, comp );
      VEC_DS_TF->set_item( global_number_in_distributed_vector,
          rhs[dir].interface_T[comp]->item( proc_pos ) );
   }
}

//----------------------------------------------------------------------
void
DDS_HeatTransferSystem::synchronize_DS_solution_vec( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: synchronize_DS_solution_vec" ) ;

   VEC_DS_TF->synchronize();
}

//----------------------------------------------------------------------
void
DDS_HeatTransferSystem::compute_product_matrix(struct TDMatrix *arr, struct ProdMatrix *prr, size_t const& comp, size_t const& dir, size_t const& r_index )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransferSystem:: compute_product_matrix" ) ;

   size_t proc_pos, nb_procs;

   proc_pos = proc_pos_in_i[dir];
   nb_procs = nb_procs_in_i[dir];

   if (proc_pos == nb_procs - 1){
      // Condition for serial processor and multi processor
      if (proc_pos == 0) {
         compute_product_matrix_interior(arr,prr,comp,proc_pos,dir,r_index);
      } else {
         compute_product_matrix_interior(arr,prr,comp,proc_pos-1,dir,r_index);
         if (is_iperiodic[dir] == 1) compute_product_matrix_interior(arr,prr,comp,proc_pos,dir,r_index);
      }
   }else if(proc_pos == 0){
      compute_product_matrix_interior(arr,prr,comp,proc_pos,dir,r_index);
      if (is_iperiodic[dir] == 1) compute_product_matrix_interior(arr,prr,comp,nb_procs-1,dir,r_index);
   }else{
      compute_product_matrix_interior(arr,prr,comp,proc_pos-1,dir,r_index);
      compute_product_matrix_interior(arr,prr,comp,proc_pos,dir,r_index);
   }

}

//----------------------------------------------------------------------
void
DDS_HeatTransferSystem::compute_product_matrix_interior(struct TDMatrix *arr,struct ProdMatrix *prr, size_t const& comp, size_t const& column,size_t const& dir,size_t const& r_index)
//----------------------------------------------------------------------
{

  // Get appropriate column of Aie
  arr[dir].ie[comp][r_index] -> extract_col(column, prr[dir].result[comp]);

  // Get inv(Aii)*Aie for for appropriate column of Aie
  DDS_HeatTransferSystem::mod_thomas_algorithm(arr, prr[dir].result[comp], comp, dir,r_index);

  // Get product of Aei*inv(Aii)*Aie for appropriate column
  arr[dir].ei[comp][r_index]->multiply_vec_then_add(prr[dir].result[comp],prr[dir].ii_ie[comp]);

  size_t int_unknown = prr[dir].ii_ie[comp]->nb_rows();

  for (size_t i = 0; i < int_unknown; i++){
      prr[dir].ei_ii_ie[comp]->set_item(i,column,prr[dir].ii_ie[comp]->item(i));
  }

}

//----------------------------------------------------------------------
void
DDS_HeatTransferSystem::display_debug(void)
//----------------------------------------------------------------------
{
   //Aii_y_main_diagonal[0]->print_items(MAC::out(),0);
}
