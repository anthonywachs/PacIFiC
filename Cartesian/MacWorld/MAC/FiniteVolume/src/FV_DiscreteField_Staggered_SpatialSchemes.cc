#include <FV_DiscreteField_Staggered.hh>
#include <FV_Mesh.hh>
#include <FV_BoundaryCondition.hh>
#include <FV_DomainBuilder.hh>
#include <FV_TimeIterator.hh>
#include <FV.hh>
#include <MAC_Communicator.hh>
#include <MAC.hh>
#include <MAC_Exec.hh>
#include <MAC_Communicator.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_DoubleArray2D.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Error.hh>
#include <MAC_String.hh>
#include <MAC_Variable.hh>
#include <LA_Vector.hh>
#include <LA_Matrix.hh>
#include <doubleArray2D.hh>
#include <intArray3D.hh>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

using std::cout ; 
using std::endl ;
using std::string ; 
using std::ostringstream ;
using std::ofstream ;
using std::ios ;


//----------------------------------------------------------------------
double
FV_DiscreteField_Staggered:: compute_CFL( FV_TimeIterator const* t_it,
	size_t level ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Staggered:: compute_CFL" );   
   MAC_CHECK_PRE( level < STO_DEPTH ) ;
   
   double local_cfl = 0., CFL = 0., dl, velo = 0., dt = t_it->time_step() ;
   size_t component;
  
   size_t_vector max_index(DIM,0);

   // X component
   component = 0;
   
   // Get local max indices
   for (size_t l=0;l<DIM;++l) 
     max_index(l) = (*local_dof_number)[component](l);

   for (size_t i=0;i<max_index(0);++i)
   {          
     dl = (*local_cell_size)[component][component](i) ;      
     for (size_t j=0;j<max_index(1);++j)
     {
       if ( DIM == 2 )
       {
         size_t k = 0 ;
	 velo = fabs( (*VALUES)[level][component]( i, j, k ) );
         local_cfl = velo * dt / dl;
         CFL = local_cfl > CFL ? local_cfl : CFL; 
       }
       else
       {
         for (size_t k=0;k<max_index(2);++k)
	 {
	   velo = fabs( (*VALUES)[level][component]( i, j, k ) );
           local_cfl = velo * dt / dl;
           CFL = local_cfl > CFL ? local_cfl : CFL;
	 }
       }
     }
   }
    
   // Y component
   component = 1;
   
   // Get local max indices
   for (size_t l=0;l<DIM;++l) 
     max_index(l) = (*local_dof_number)[component](l);

   for (size_t j=0;j<max_index(1);++j)
   {          
     dl = (*local_cell_size)[component][component](j) ;      
     for (size_t i=0;i<max_index(0);++i)
     {
       if ( DIM == 2 )
       {
         size_t k = 0 ;
	 velo = fabs( (*VALUES)[level][component]( i, j, k ) );
         local_cfl = velo * dt / dl;
         CFL = local_cfl > CFL ? local_cfl : CFL; 
       }
       else
       {
         for (size_t k=0;k<max_index(2);++k)
	 {
	   velo = fabs( (*VALUES)[level][component]( i, j, k ) );
           local_cfl = velo * dt / dl;
           CFL = local_cfl > CFL ? local_cfl : CFL;
	 }
       }
     }
   }

   if ( DIM == 3 )
   {
     // Z component
     component = 2;
     
     // Get local max indices
     for (size_t l=0;l<DIM;++l) 
       max_index(l) = (*local_dof_number)[component](l);

     for (size_t k=0;k<max_index(2);++k)
     {          
       dl = (*local_cell_size)[component][component](k) ;
       for (size_t i=0;i<max_index(0);++i)
       {
         for (size_t j=0;j<max_index(1);++j)
	 {
	   velo = fabs( (*VALUES)[level][component]( i, j, k ) );
           local_cfl = velo * dt / dl;
           CFL = local_cfl > CFL ? local_cfl : CFL;
	 }
       }
     }
   }
   
   double collective_CFL = MAC_Exec::communicator()->max(CFL);

   return collective_CFL;

}




// //----------------------------------------------------------------------
// void 
// FV_DiscreteField_Staggered:: assemble_advection_Upwind( 
// 	FV_DiscreteField const* AdvectingField,
// 	size_t advecting_level, double const& coef, size_t advected_level,
// 	LA_Vector *VEC_rhs ) const
// //----------------------------------------------------------------------
// {
//    MAC_LABEL( "FV_DiscreteField_Staggered:: assemble_advection_Upwind" );   
//    MAC_CHECK_PRE( advected_level < STO_DEPTH ) ;
//    MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
//    MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ;       
// 
//    // Parameters
//    size_t_vector min_unknown_index(DIM,0);
//    size_t_vector max_unknown_index(DIM,0);
//    size_t center_pos_in_matrix = 0, component ;
//    double dxC, dyC, dzC;
//    double AdvectedValueC = 0., AdvectedValueRi = 0., AdvectedValueLe = 0.,
//    	AdvectedValueTo = 0., AdvectedValueBo = 0., AdvectedValueFr = 0., 
// 	AdvectedValueBe = 0,
//    	AdvectorValueC = 0., AdvectorValueRi = 0., AdvectorValueLe = 0.,
//    	AdvectorValueTo = 0., AdvectorValueBo = 0., AdvectorValueFr = 0., 
// 	AdvectorValueBe = 0, AdvectorValueToLe = 0., AdvectorValueToRi = 0., 
// 	AdvectorValueBoLe = 0., AdvectorValueBoRi = 0., AdvectorValueFrLe = 0., 
// 	AdvectorValueFrRi = 0., AdvectorValueBeLe = 0., AdvectorValueBeRi = 0.,
// 	AdvectorValueFrTo = 0., AdvectorValueFrBo = 0., AdvectorValueBeTo = 0., 
// 	AdvectorValueBeBo = 0.,
// 	ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
// 	fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
//    FV_SHIFT_TRIPLET shift ;
// 
//    // Nullify vector
//    VEC_rhs->nullify();
// 
//    // Assemble vector
//    for( component=0; component<NB_COMPS; ++component )
//    {
//      // Get local min and max indices
//      for( size_t l=0; l<DIM; ++l )
//      {
//        min_unknown_index(l) = 
//        	(*min_index_unknown_handled_by_proc)[component](l);
//        max_unknown_index(l) =
//         (*max_index_unknown_handled_by_proc)[component](l);
//      }
//      shift = shift_staggeredToStaggered( component ) ;
// 
//      // Perform assembling
//      for( size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
//      {
//        dxC = (*local_cell_size)[component][0](i) ;      
//        for( size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
//        {
// 	 dyC = (*local_cell_size)[component][1](j) ; 
//  
//          if( DIM == 2 )
// 	 {
// 	   size_t k = 0 ;
// 	   center_pos_in_matrix = DOF_global_number( i, j, k, component );
//            AdvectedValueC = DOF_value( i, j, k, component, 
// 	   	advected_level );
// 	   AdvectorValueC = AdvectingField->DOF_value( i, j, k, component, 
// 	   	advecting_level );
// 	   
// 	   // The First Component (u)
// 	   if ( component == 0 )
// 	   {
// 	     // Right (U_X)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	       fri = AdvectorValueC * AdvectedValueC;
// 	     else
// 	     {
//                AdvectedValueRi = DOF_value(
//                    i+1, j, k, component, advected_level );
// 	       AdvectorValueRi = AdvectingField->DOF_value( 
// 	           i+1, j, k, component, advecting_level );
// 	       ur = 0.5 * ( AdvectorValueC + AdvectorValueRi );
// 	       if ( ur > 0. ) fri = ur * AdvectedValueC;
// 	       else fri = ur * AdvectedValueRi;
// 	     }
// 	         
// 	     // Left (U_X)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	       fle = AdvectorValueC * AdvectedValueC;
// 	     else
// 	     {
//                AdvectedValueLe = DOF_value( 
// 	       	     i-1, j, k, component, advected_level );
// 	       AdvectorValueLe = AdvectingField->DOF_value( 
// 	           i-1, j, k, component, advecting_level );
// 	       ul = 0.5 * ( AdvectorValueC + AdvectorValueLe );
// 	       if ( ul > 0. ) fle = ul * AdvectedValueLe;
// 	       else fle = ul * AdvectedValueC;
// 	     }
// 	     
// 	     // Top (U_Y)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	       fto = 0.;
// 	     else
// 	     {
//                AdvectedValueTo = DOF_value(
// 	       	     i, j+1, k, component, advected_level );
// 	       AdvectorValueToLe = AdvectingField->DOF_value(
// 	           i+shift.i-1, j+shift.j, k, 1, advecting_level );
// 	       AdvectorValueToRi = AdvectingField->DOF_value(
// 	           i+shift.i, j+shift.j, k, 1, advecting_level );
// 	       vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
// 	       if ( vt > 0. ) fto = vt * AdvectedValueC;
// 	       else fto = vt * AdvectedValueTo;
// 	     }
// 	 
// 	     // Bottom (U_Y)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	       fbo = 0.;
// 	     else
// 	     {
// 	       AdvectedValueBo = DOF_value(
// 	       	     i, j-1, k, component, advected_level );
// 	       AdvectorValueBoLe = AdvectingField->DOF_value(
// 	           i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
// 	       AdvectorValueBoRi = AdvectingField->DOF_value(
// 	           i+shift.i, j+shift.j-1, k, 1, advecting_level );
// 	       vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
// 	       if ( vb > 0. ) fbo = vb * AdvectedValueBo;
// 	       else fbo = vb * AdvectedValueC;
// 	     }
// 
// 	   }
// 	   // The second Component (v)
// 	   else
// 	   {
// 	     // Right (V_X)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	       fri = 0.;
// 	     else
// 	     {
// 	       AdvectedValueRi = DOF_value(
// 	       	     i+1, j, k, component, advected_level );
// 	       AdvectorValueToRi = AdvectingField->DOF_value(
// 	           i+shift.i, j+shift.j, k, 0, advecting_level );
// 	       AdvectorValueBoRi = AdvectingField->DOF_value(
// 	           i+shift.i, j+shift.j-1, k, 0, advecting_level );
// 	       ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
// 	       if ( ur > 0. ) fri = ur * AdvectedValueC;
// 	       else fri = ur * AdvectedValueRi;
// 	     }
// 	         
// 	     // Left (V_X)
// 	     if ( DOF_color(i, j, k, component ) == FV_BC_LEFT )
// 	       fle = 0.;
// 	     else
// 	     {
// 	       AdvectedValueLe = DOF_value(
// 	       	     i-1, j, k, component, advected_level );
// 	       AdvectorValueToLe = AdvectingField->DOF_value(
// 	           i+shift.i-1, j+shift.j, k, 0, advecting_level );
// 	       AdvectorValueBoLe = AdvectingField->DOF_value(
// 	           i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
// 	       ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
// 	       if ( ul > 0. ) fle = ul * AdvectedValueLe;
// 	       else fle = ul * AdvectedValueC;
// 	     }
// 	 
// 	     // Top (V_Y)
// 	     if ( DOF_color(i, j, k, component ) == FV_BC_TOP )
// 	       fto = AdvectorValueC * AdvectedValueC;
// 	     else
// 	     {
// 	       AdvectedValueTo = DOF_value(
// 	       	     i, j+1, k, component, advected_level );
// 	       AdvectorValueTo = AdvectingField->DOF_value(
// 	           i, j+1, k, component, advecting_level );
// 	       vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
// 	       if ( vt > 0. ) fto = vt * AdvectedValueC;
// 	       else fto = vt * AdvectedValueTo;
// 	     }
// 	 
// 	     // Bottom (V_Y)
// 	     if ( DOF_color(i, j, k, component ) == FV_BC_BOTTOM )
// 	       fbo = AdvectorValueC * AdvectedValueC;
// 	     else
// 	     {
// 	       AdvectedValueBo = DOF_value(
// 	       	     i, j-1, k, component, advected_level );
// 	       AdvectorValueBo = AdvectingField->DOF_value(
// 	           i, j-1, k, component, advecting_level );
// 	       vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
// 	       if ( vb > 0. ) fbo = vb * AdvectedValueBo;
// 	       else fbo = vb * AdvectedValueC;
// 	     }
// 	   }
// 
//            flux = (fto - fbo) * dxC + (fri - fle) * dyC;
//            VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
// 	 }
// 	 else // DIM = 3
// 	 {
// 	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
// 	   {
// 	     dzC = (*local_cell_size)[component][2](k) ; 	     
// 	     center_pos_in_matrix = DOF_global_number( i, j, k, component );
//              AdvectedValueC = DOF_value( i, j, k, component, 
// 	     	advected_level );
//              AdvectorValueC = AdvectingField->DOF_value( i, j, k, component, 
// 	     	advecting_level );
// 
// 	     // The First Component (u)
// 	     if ( component == 0 )
// 	     {
// 	       // Right (U_X)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	         fri = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {
//                  AdvectedValueRi = DOF_value(
// 		       i+1, j, k, component, advected_level );
// 		 AdvectorValueRi = AdvectingField->DOF_value(
// 		     i+1, j, k, component, advecting_level );
//                  ur = 0.5 * ( AdvectorValueRi + AdvectorValueC );
//                  if ( ur > 0. ) fri = ur * AdvectedValueC;
//                  else fri = ur * AdvectedValueRi;
// 	       }
// 	         
// 	       // Left (U_X)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	         fle = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {
//                  AdvectedValueLe = DOF_value(
// 		       i-1, j, k, component, advected_level );
// 		 AdvectorValueLe = AdvectingField->DOF_value(
// 		     i-1, j, k, component, advecting_level );
//                  ul = 0.5 * ( AdvectorValueLe + AdvectorValueC );
//                  if ( ul > 0. ) fle = ul * AdvectedValueLe;
//                  else fle = ul * AdvectedValueC;
// 	       }
// 	 
// 	       // Top (U_Y)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	         fto = 0.;
// 	       else
// 	       {
// 	         AdvectedValueTo = DOF_value(
// 		       i, j+1, k, component, advected_level );
// 	         AdvectorValueToLe = AdvectingField->DOF_value( 
// 		     i+shift.i-1, j+shift.j, k, 1, advecting_level );
// 	         AdvectorValueToRi = AdvectingField->DOF_value( 
// 		     i+shift.i, j+shift.j, k, 1, advecting_level );
// 	         vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
// 	         if ( vt > 0. ) fto = vt * AdvectedValueC;
// 	         else fto = vt * AdvectedValueTo;
// 	       }
// 	 
// 	       // Bottom (U_Y)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	         fbo = 0.;
// 	       else
// 	       {
// 	         AdvectedValueBo = DOF_value(
// 		       i, j-1, k, component, advected_level );
// 		 AdvectorValueBoLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
// 	         AdvectorValueBoRi = AdvectingField->DOF_value(
// 		     i+shift.i, j+shift.j-1, k, 1, advecting_level );
// 	         vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
// 	         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
// 	         else fbo = vb * AdvectedValueC;
// 	       }
//       
// 	       // Front (U_Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
// 	         ffr = 0.;
// 	       else
// 	       {
//                  AdvectedValueFr = DOF_value(
// 		       i, j, k+1, component, advected_level );
// 		 AdvectorValueFrLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j, k+shift.k, 2, advecting_level );
// 	         AdvectorValueFrRi = AdvectingField->DOF_value(
// 		     i+shift.i, j, k+shift.k, 2, advecting_level );
// 	         wf = 0.5 * ( AdvectorValueFrLe + AdvectorValueFrRi );
// 	         if ( wf > 0. ) ffr = wf * AdvectedValueC;
// 	         else ffr = wf * AdvectedValueFr;
// 	       }
// 	   
// 	       // Behind (U_Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND)
// 	         fbe = 0.;
// 	       else
// 	       {
// 	         AdvectedValueBe = DOF_value(
// 		       i, j, k-1, component, advected_level );
// 		 AdvectorValueBeLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j, k+shift.k-1, 2, advecting_level );
// 	         AdvectorValueBeRi = AdvectingField->DOF_value(
// 		     i+shift.i, j, k+shift.k-1, 2, advecting_level );
// 	         wb = 0.5 * ( AdvectorValueBeLe + AdvectorValueBeRi );
// 	         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
// 	         else fbe = wb * AdvectedValueC;
// 	       }
// 	     }
// 	     // The Second Component (v)
// 	     else if ( component == 1 )
// 	     {
// 	       // Right (V_X)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	         fri = 0.;
// 	       else
// 	       {	     
// 	         AdvectedValueRi = DOF_value(
// 		       i+1, j, k, component, advected_level );
// 	         AdvectorValueToRi = AdvectingField->DOF_value(
// 		     i+shift.i, j+shift.j, k, 0, advecting_level );
// 	         AdvectorValueBoRi = AdvectingField->DOF_value(
// 		     i+shift.i, j+shift.j-1, k, 0, advecting_level );
// 	         ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
// 	         if ( ur > 0. ) fri = ur * AdvectedValueC;
// 	         else fri = ur * AdvectedValueRi;
// 	       }
// 	         
// 	       // Left (V_X)
// 	       if ( DOF_color(i, j, k, component ) == FV_BC_LEFT )
// 	         fle = 0.;
// 	       else
// 	       {	     
// 	         AdvectedValueLe = DOF_value(
// 		       i-1, j, k, component, advected_level );
// 	         AdvectorValueToLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j+shift.j, k, 0, advecting_level );
// 	         AdvectorValueBoLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
// 	         ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
// 	         if ( ul > 0. ) fle = ul * AdvectedValueLe;
// 	         else fle = ul * AdvectedValueC;
// 	       }
// 	 
// 	       // Top (V_Y)
// 	       if ( DOF_color(i, j, k, component ) == FV_BC_TOP )
// 	         fto = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {
// 	         AdvectedValueTo = DOF_value(
// 		       i, j+1, k, component, advected_level );
// 		 AdvectorValueTo = AdvectingField->DOF_value(
// 		     i, j+1, k, component, advecting_level );
// 	         vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
// 	         if ( vt > 0. ) fto = vt * AdvectedValueC;
// 	         else fto = vt * AdvectedValueTo;
// 	       }
// 	 
// 	       // Bottom (V_Y)
// 	       if ( DOF_color(i, j, k, component ) == FV_BC_BOTTOM )
// 	         fbo = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {
// 	         AdvectedValueBo = DOF_value(
// 		       i, j-1, k, component, advected_level );
// 		 AdvectorValueBo = AdvectingField->DOF_value(
// 		     i, j-1, k, component, advecting_level );
// 	         vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
// 	         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
// 	         else fbo = vb * AdvectedValueC;
// 	       }
//       
// 	       // Front (V_Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
// 	         ffr = 0.;
// 	       else
// 	       {
// 	         AdvectedValueFr = DOF_value(
// 		       i, j, k+1, component, advected_level );
// 		 AdvectorValueFrTo = AdvectingField->DOF_value(
// 		     i, j+shift.j, k+shift.k, 2, advecting_level );
// 	         AdvectorValueFrBo = AdvectingField->DOF_value(
// 		     i, j+shift.j-1, k+shift.k, 2, advecting_level );
// 	         wf = 0.5 * ( AdvectorValueFrTo + AdvectorValueFrBo );
// 	         if ( wf > 0. ) ffr = wf * AdvectedValueC;
// 	         else ffr = wf * AdvectedValueFr;
// 	       }
// 	   
// 	       // Behind (V_Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
// 	         fbe = 0.;
// 	       else
// 	       {	     
// 	         AdvectedValueBe = DOF_value(
// 		       i, j, k-1, component, advected_level );
// 		 AdvectorValueBeTo = AdvectingField->DOF_value(
// 		     i, j+shift.j, k+shift.k-1, 2, advecting_level );
// 	         AdvectorValueBeBo = AdvectingField->DOF_value(
// 		     i, j+shift.j-1, k+shift.k-1, 2, advecting_level );
// 	         wb = 0.5 * ( AdvectorValueBeTo + AdvectorValueBeBo );
// 	         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
// 	         else fbe = wb * AdvectedValueC;
// 	       }
// 	     }
// 	     // The Third Component (w)
// 	     else
// 	     {
// 	       // Right (W_X)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	         fri = 0.;
// 	       else
// 	       {	     
// 	         AdvectedValueRi = DOF_value(
// 		       i+1, j, k, component, advected_level );
// 	         AdvectorValueFrRi = AdvectingField->DOF_value(
//   		     i+shift.i, j, k+shift.k, 0, advecting_level );
// 	         AdvectorValueBeRi = AdvectingField->DOF_value(
// 		     i+shift.i, j, k+shift.k-1, 0, advecting_level );
// 	         ur = 0.5 * ( AdvectorValueFrRi + AdvectorValueBeRi );
// 	         if ( ur > 0. ) fri = ur * AdvectedValueC;
// 	         else fri = ur * AdvectedValueRi;
// 	       }
// 	         
// 	       // Left (W_X)
// 	       if ( DOF_color(i, j, k, component ) == FV_BC_LEFT )
// 	         fle = 0.;
// 	       else
// 	       {	     
//                  AdvectedValueLe = DOF_value(
// 		       i-1, j, k, component, advected_level );
// 		 AdvectorValueFrLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j, k+shift.k, 0, advecting_level );
// 	         AdvectorValueBeLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j, k+shift.k-1, 0, advecting_level );
// 	         ul = 0.5 * ( AdvectorValueFrLe + AdvectorValueBeLe );
// 	         if ( ul > 0. ) fle = ul * AdvectedValueLe;
// 	         else fle = ul * AdvectedValueC;
// 	       }
// 
// 	       // Top (W_Y)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	         fto = 0.;
// 	       else
// 	       {
// 	         AdvectedValueTo = DOF_value(
// 		       i, j+1, k, component, advected_level );
// 	         AdvectorValueFrTo = AdvectingField->DOF_value(
// 		 		i, j+shift.j, k+shift.k, 1, advecting_level );
// 	         AdvectorValueBeTo = AdvectingField->DOF_value(
// 		 		i, j+shift.j, k+shift.k-1, 1, advecting_level );
// 	         vt = 0.5 * ( AdvectorValueFrTo + AdvectorValueBeTo );
// 	         if ( vt > 0. ) fto = vt * AdvectedValueC;
// 	         else fto = vt * AdvectedValueTo;
// 	       }
// 	 
// 	       // Bottom (W_Y)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	         fbo = 0.;
// 	       else
// 	       {	     
// 	         AdvectedValueBo = DOF_value(
// 		       i, j-1, k, component, advected_level );
// 	         AdvectorValueFrBo = AdvectingField->DOF_value(
// 		     i, j+shift.j-1, k+shift.k, 1, advecting_level );
// 	         AdvectorValueBeBo = AdvectingField->DOF_value(
// 		     i, j+shift.j-1, k+shift.k-1, 1, advecting_level );
// 	         vb = 0.5 * ( AdvectorValueFrBo + AdvectorValueBeBo );
// 	         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
// 	         else fbo = vb * AdvectedValueC;
// 	       }
//       
// 	       // Front (W_Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
// 	         ffr = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {
// 	         AdvectedValueFr = DOF_value(
// 		       i, j, k+1, component, advected_level );
// 	         AdvectorValueFr = AdvectingField->DOF_value(
// 		     i, j, k+1, component, advecting_level );
// 	         wf = 0.5 * ( AdvectorValueFr + AdvectorValueC );
// 	         if ( wf > 0. ) ffr = wf * AdvectedValueC;
// 	         else ffr = wf * AdvectedValueFr;
// 	       }
// 	   
// 	       // Behind (W_Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
// 	         fbe = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {
// 	         AdvectedValueBe = DOF_value(
// 		       i, j, k-1, component, advected_level );
// 	         AdvectorValueBe = AdvectingField->DOF_value(
// 		     i, j, k-1, component, advecting_level );
// 	         wb = 0.5 * ( AdvectorValueBe + AdvectorValueC );
// 	         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
// 	         else fbe = wb * AdvectedValueC;
// 	       }
// 	     }
// 	     
//              flux = (fto - fbo) * dxC * dzC
// 	     		+ (fri - fle) * dyC * dzC
// 			+ (ffr - fbe) * dxC * dyC;
//              VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
// 	   }
// 	 }
//        }
//      }
//    }
// 
//    // Synchronize vector for parallel usage
//    VEC_rhs->synchronize() ;   
// 
// }	
//----------------------------------------------------------------------
void 
FV_DiscreteField_Staggered:: assemble_advection_Upwind( 
	FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t advected_level,
	LA_Vector *VEC_rhs ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Staggered:: assemble_advection_Upwind" );   
   MAC_CHECK_PRE( advected_level < STO_DEPTH ) ;
   MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
   MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ;       

   // Parameters
   size_t_vector min_unknown_index(DIM,0);
   size_t_vector max_unknown_index(DIM,0);
   size_t center_pos_in_matrix = 0, component ;
   double dxC = 0., dyC = 0., dzC = 0.;
   double AdvectedValueC = 0., AdvectedValueRi = 0., AdvectedValueLe = 0.,
   	AdvectedValueTo = 0., AdvectedValueBo = 0., AdvectedValueFr = 0., 
	AdvectedValueBe = 0,
   	AdvectorValueC = 0., AdvectorValueRi = 0., AdvectorValueLe = 0.,
   	AdvectorValueTo = 0., AdvectorValueBo = 0., AdvectorValueFr = 0., 
	AdvectorValueBe = 0, AdvectorValueToLe = 0., AdvectorValueToRi = 0., 
	AdvectorValueBoLe = 0., AdvectorValueBoRi = 0., AdvectorValueFrLe = 0., 
	AdvectorValueFrRi = 0., AdvectorValueBeLe = 0., AdvectorValueBeRi = 0.,
	AdvectorValueFrTo = 0., AdvectorValueFrBo = 0., AdvectorValueBeTo = 0., 
	AdvectorValueBeBo = 0.,
	ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
	fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
   FV_SHIFT_TRIPLET shift ;

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

   // Nullify vector
   VEC_rhs->nullify();

   // Assemble vector
   for( component=0; component<NB_COMPS; ++component )
   {
     // Get local min and max indices
     for( size_t l=0; l<DIM; ++l )
     {
       min_unknown_index(l) = 
       	(*min_index_unknown_handled_by_proc)[component](l);
       max_unknown_index(l) =
        (*max_index_unknown_handled_by_proc)[component](l);
     }
     shift = shift_staggeredToStaggered( component ) ;

     // Perform assembling
     for( size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {
       dxC = (*local_cell_size)[component][0](i) ;      
       for( size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
	 dyC = (*local_cell_size)[component][1](j) ; 
 
         if( DIM == 2 )
	 {
	   size_t k = 0 ;
	   center_pos_in_matrix = DOF_global_number( i, j, k, component );
           AdvectedValueC = DOF_value( i, j, k, component, 
	   	advected_level );
	   AdvectorValueC = AdvectingField->DOF_value( i, j, k, component, 
	   	advecting_level );
	   
	   // The First Component (u)
	   if ( component == 0 )
	   {
	     // Right (U_X)
	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
	       fri = AdvectorValueC * AdvectedValueC;
	     else
	     {
               AdvectedValueRi = DOF_value(
                   i+1, j, k, component, advected_level );
	       AdvectorValueRi = AdvectingField->DOF_value( 
	           i+1, j, k, component, advecting_level );
	       ur = 0.5 * ( AdvectorValueC + AdvectorValueRi );
	       if ( ur > 0. ) fri = ur * AdvectedValueC;
	       else fri = ur * AdvectedValueRi;
	     }
	         
	     // Left (U_X)
	     if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
	       fle = AdvectorValueC * AdvectedValueC;
	     else
	     {
               AdvectedValueLe = DOF_value( 
	       	     i-1, j, k, component, advected_level );
	       AdvectorValueLe = AdvectingField->DOF_value( 
	           i-1, j, k, component, advecting_level );
	       ul = 0.5 * ( AdvectorValueC + AdvectorValueLe );
	       if ( ul > 0. ) fle = ul * AdvectedValueLe;
	       else fle = ul * AdvectedValueC;
	     }
	     
	     // Top (U_Y)
	     AdvectedValueTo = DOF_value(
	       	     i, j+1, k, component, advected_level );
	     AdvectorValueToLe = AdvectingField->DOF_value(
	           i+shift.i-1, j+shift.j, k, 1, advecting_level );
	     AdvectorValueToRi = AdvectingField->DOF_value(
	           i+shift.i, j+shift.j, k, 1, advecting_level );
	     vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
	     if ( vt > 0. ) fto = vt * AdvectedValueC;
	     else fto = vt * AdvectedValueTo;
	 
	     // Bottom (U_Y)
	     AdvectedValueBo = DOF_value(
	       	     i, j-1, k, component, advected_level );
	     AdvectorValueBoLe = AdvectingField->DOF_value(
	           i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
	     AdvectorValueBoRi = AdvectingField->DOF_value(
	           i+shift.i, j+shift.j-1, k, 1, advecting_level );
	     vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
	     if ( vb > 0. ) fbo = vb * AdvectedValueBo;
	     else fbo = vb * AdvectedValueC;
	   }
	   
	   // The second Component (v)
	   else
	   {
	     // Right (V_X)
	     AdvectedValueRi = DOF_value(
	       	     i+1, j, k, component, advected_level );
	     AdvectorValueToRi = AdvectingField->DOF_value(
	           i+shift.i, j+shift.j, k, 0, advecting_level );
	     AdvectorValueBoRi = AdvectingField->DOF_value(
	           i+shift.i, j+shift.j-1, k, 0, advecting_level );
	     ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
	     if ( ur > 0. ) fri = ur * AdvectedValueC;
	     else fri = ur * AdvectedValueRi;
	         
	     // Left (V_X)
	     AdvectedValueLe = DOF_value(
	       	     i-1, j, k, component, advected_level );
	     AdvectorValueToLe = AdvectingField->DOF_value(
	           i+shift.i-1, j+shift.j, k, 0, advecting_level );
	     AdvectorValueBoLe = AdvectingField->DOF_value(
	           i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
	     ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
	     if ( ul > 0. ) fle = ul * AdvectedValueLe;
	     else fle = ul * AdvectedValueC;
	 
	     // Top (V_Y)
	     if ( DOF_color(i, j, k, component ) == FV_BC_TOP )
	       fto = AdvectorValueC * AdvectedValueC;
	     else
	     {
	       AdvectedValueTo = DOF_value(
	       	     i, j+1, k, component, advected_level );
	       AdvectorValueTo = AdvectingField->DOF_value(
	           i, j+1, k, component, advecting_level );
	       vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
	       if ( vt > 0. ) fto = vt * AdvectedValueC;
	       else fto = vt * AdvectedValueTo;
	     }
	 
	     // Bottom (V_Y)
	     if ( DOF_color(i, j, k, component ) == FV_BC_BOTTOM )
	       fbo = AdvectorValueC * AdvectedValueC;
	     else
	     {
	       AdvectedValueBo = DOF_value(
	       	     i, j-1, k, component, advected_level );
	       AdvectorValueBo = AdvectingField->DOF_value(
	           i, j-1, k, component, advecting_level );
	       vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
	       if ( vb > 0. ) fbo = vb * AdvectedValueBo;
	       else fbo = vb * AdvectedValueC;
	     }
	   }

           flux = (fto - fbo) * dxC + (fri - fle) * dyC;
           VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
	 }
	 else // DIM = 3
	 {
	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	   {
	     dzC = (*local_cell_size)[component][2](k) ; 	     
	     center_pos_in_matrix = DOF_global_number( i, j, k, component );
             AdvectedValueC = DOF_value( i, j, k, component, 
	     	advected_level );
             AdvectorValueC = AdvectingField->DOF_value( i, j, k, component, 
	     	advecting_level );

	     // The First Component (u)
	     if ( component == 0 )
	     {
	       // Right (U_X)
	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
	         fri = AdvectorValueC * AdvectedValueC;
	       else
	       {
                 AdvectedValueRi = DOF_value(
		       i+1, j, k, component, advected_level );
		 AdvectorValueRi = AdvectingField->DOF_value(
		     i+1, j, k, component, advecting_level );
                 ur = 0.5 * ( AdvectorValueRi + AdvectorValueC );
                 if ( ur > 0. ) fri = ur * AdvectedValueC;
                 else fri = ur * AdvectedValueRi;
	       }
	         
	       // Left (U_X)
	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
	         fle = AdvectorValueC * AdvectedValueC;
	       else
	       {
                 AdvectedValueLe = DOF_value(
		       i-1, j, k, component, advected_level );
		 AdvectorValueLe = AdvectingField->DOF_value(
		     i-1, j, k, component, advecting_level );
                 ul = 0.5 * ( AdvectorValueLe + AdvectorValueC );
                 if ( ul > 0. ) fle = ul * AdvectedValueLe;
                 else fle = ul * AdvectedValueC;
	       }
	 
	       // Top (U_Y)
	       AdvectedValueTo = DOF_value(
		       i, j+1, k, component, advected_level );
	       AdvectorValueToLe = AdvectingField->DOF_value( 
		     i+shift.i-1, j+shift.j, k, 1, advecting_level );
	       AdvectorValueToRi = AdvectingField->DOF_value( 
		     i+shift.i, j+shift.j, k, 1, advecting_level );
	       vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
	       if ( vt > 0. ) fto = vt * AdvectedValueC;
	       else fto = vt * AdvectedValueTo;
	 
	       // Bottom (U_Y)
	       AdvectedValueBo = DOF_value(
		       i, j-1, k, component, advected_level );
	       AdvectorValueBoLe = AdvectingField->DOF_value(
		     i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
	       AdvectorValueBoRi = AdvectingField->DOF_value(
		     i+shift.i, j+shift.j-1, k, 1, advecting_level );
	       vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
	       if ( vb > 0. ) fbo = vb * AdvectedValueBo;
	       else fbo = vb * AdvectedValueC;
      
	       // Front (U_Z)
	       AdvectedValueFr = DOF_value(
		       i, j, k+1, component, advected_level );
	       AdvectorValueFrLe = AdvectingField->DOF_value(
		     i+shift.i-1, j, k+shift.k, 2, advecting_level );
	       AdvectorValueFrRi = AdvectingField->DOF_value(
		     i+shift.i, j, k+shift.k, 2, advecting_level );
	       wf = 0.5 * ( AdvectorValueFrLe + AdvectorValueFrRi );
	       if ( wf > 0. ) ffr = wf * AdvectedValueC;
	       else ffr = wf * AdvectedValueFr;
	   
	       // Behind (U_Z)
	       AdvectedValueBe = DOF_value(
		       i, j, k-1, component, advected_level );
	       AdvectorValueBeLe = AdvectingField->DOF_value(
		     i+shift.i-1, j, k+shift.k-1, 2, advecting_level );
	       AdvectorValueBeRi = AdvectingField->DOF_value(
		     i+shift.i, j, k+shift.k-1, 2, advecting_level );
	       wb = 0.5 * ( AdvectorValueBeLe + AdvectorValueBeRi );
	       if ( wb > 0. ) fbe = wb * AdvectedValueBe;
	       else fbe = wb * AdvectedValueC;
	     }
	     
	     // The Second Component (v)
	     else if ( component == 1 )
	     {
	       // Right (V_X)
	       AdvectedValueRi = DOF_value(
		       i+1, j, k, component, advected_level );
	       AdvectorValueToRi = AdvectingField->DOF_value(
		     i+shift.i, j+shift.j, k, 0, advecting_level );
	       AdvectorValueBoRi = AdvectingField->DOF_value(
		     i+shift.i, j+shift.j-1, k, 0, advecting_level );
	       ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
	       if ( ur > 0. ) fri = ur * AdvectedValueC;
	       else fri = ur * AdvectedValueRi;
	         
	       // Left (V_X)
	       AdvectedValueLe = DOF_value(
		       i-1, j, k, component, advected_level );
	       AdvectorValueToLe = AdvectingField->DOF_value(
		     i+shift.i-1, j+shift.j, k, 0, advecting_level );
	       AdvectorValueBoLe = AdvectingField->DOF_value(
		     i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
	       ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
	       if ( ul > 0. ) fle = ul * AdvectedValueLe;
	       else fle = ul * AdvectedValueC;
	 
	       // Top (V_Y)
	       if ( DOF_color(i, j, k, component ) == FV_BC_TOP )
	         fto = AdvectorValueC * AdvectedValueC;
	       else
	       {
	         AdvectedValueTo = DOF_value(
		       i, j+1, k, component, advected_level );
		 AdvectorValueTo = AdvectingField->DOF_value(
		     i, j+1, k, component, advecting_level );
	         vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
	         if ( vt > 0. ) fto = vt * AdvectedValueC;
	         else fto = vt * AdvectedValueTo;
	       }
	 
	       // Bottom (V_Y)
	       if ( DOF_color(i, j, k, component ) == FV_BC_BOTTOM )
	         fbo = AdvectorValueC * AdvectedValueC;
	       else
	       {
	         AdvectedValueBo = DOF_value(
		       i, j-1, k, component, advected_level );
		 AdvectorValueBo = AdvectingField->DOF_value(
		     i, j-1, k, component, advecting_level );
	         vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
	         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
	         else fbo = vb * AdvectedValueC;
	       }
      
	       // Front (V_Z)
	       AdvectedValueFr = DOF_value(
		       i, j, k+1, component, advected_level );
	       AdvectorValueFrTo = AdvectingField->DOF_value(
		     i, j+shift.j, k+shift.k, 2, advecting_level );
	       AdvectorValueFrBo = AdvectingField->DOF_value(
		     i, j+shift.j-1, k+shift.k, 2, advecting_level );
	       wf = 0.5 * ( AdvectorValueFrTo + AdvectorValueFrBo );
	       if ( wf > 0. ) ffr = wf * AdvectedValueC;
	       else ffr = wf * AdvectedValueFr;
	   
	       // Behind (V_Z)
	       AdvectedValueBe = DOF_value(
		       i, j, k-1, component, advected_level );
	       AdvectorValueBeTo = AdvectingField->DOF_value(
		     i, j+shift.j, k+shift.k-1, 2, advecting_level );
	       AdvectorValueBeBo = AdvectingField->DOF_value(
		     i, j+shift.j-1, k+shift.k-1, 2, advecting_level );
	       wb = 0.5 * ( AdvectorValueBeTo + AdvectorValueBeBo );
	       if ( wb > 0. ) fbe = wb * AdvectedValueBe;
	       else fbe = wb * AdvectedValueC;
	     }
	     
	     // The Third Component (w)
	     else
	     {
	       // Right (W_X)
	       AdvectedValueRi = DOF_value(
		       i+1, j, k, component, advected_level );
	       AdvectorValueFrRi = AdvectingField->DOF_value(
  		     i+shift.i, j, k+shift.k, 0, advecting_level );
	       AdvectorValueBeRi = AdvectingField->DOF_value(
		     i+shift.i, j, k+shift.k-1, 0, advecting_level );
	       ur = 0.5 * ( AdvectorValueFrRi + AdvectorValueBeRi );
	       if ( ur > 0. ) fri = ur * AdvectedValueC;
	       else fri = ur * AdvectedValueRi;
	         
	       // Left (W_X)
	       AdvectedValueLe = DOF_value(
		       i-1, j, k, component, advected_level );
	       AdvectorValueFrLe = AdvectingField->DOF_value(
		     i+shift.i-1, j, k+shift.k, 0, advecting_level );
	       AdvectorValueBeLe = AdvectingField->DOF_value(
		     i+shift.i-1, j, k+shift.k-1, 0, advecting_level );
	       ul = 0.5 * ( AdvectorValueFrLe + AdvectorValueBeLe );
	       if ( ul > 0. ) fle = ul * AdvectedValueLe;
	       else fle = ul * AdvectedValueC;

	       // Top (W_Y)
	       AdvectedValueTo = DOF_value(
		       i, j+1, k, component, advected_level );
	       AdvectorValueFrTo = AdvectingField->DOF_value(
		 		i, j+shift.j, k+shift.k, 1, advecting_level );
	       AdvectorValueBeTo = AdvectingField->DOF_value(
		 		i, j+shift.j, k+shift.k-1, 1, advecting_level );
	       vt = 0.5 * ( AdvectorValueFrTo + AdvectorValueBeTo );
	       if ( vt > 0. ) fto = vt * AdvectedValueC;
	       else fto = vt * AdvectedValueTo;
	 
	       // Bottom (W_Y)
	       AdvectedValueBo = DOF_value(
		       i, j-1, k, component, advected_level );
	       AdvectorValueFrBo = AdvectingField->DOF_value(
		     i, j+shift.j-1, k+shift.k, 1, advecting_level );
	       AdvectorValueBeBo = AdvectingField->DOF_value(
		     i, j+shift.j-1, k+shift.k-1, 1, advecting_level );
	       vb = 0.5 * ( AdvectorValueFrBo + AdvectorValueBeBo );
	       if ( vb > 0. ) fbo = vb * AdvectedValueBo;
	       else fbo = vb * AdvectedValueC;
      
	       // Front (W_Z)
	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
	         ffr = AdvectorValueC * AdvectedValueC;
	       else
	       {
	         AdvectedValueFr = DOF_value(
		       i, j, k+1, component, advected_level );
	         AdvectorValueFr = AdvectingField->DOF_value(
		     i, j, k+1, component, advecting_level );
	         wf = 0.5 * ( AdvectorValueFr + AdvectorValueC );
	         if ( wf > 0. ) ffr = wf * AdvectedValueC;
	         else ffr = wf * AdvectedValueFr;
	       }
	   
	       // Behind (W_Z)
	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
	         fbe = AdvectorValueC * AdvectedValueC;
	       else
	       {
	         AdvectedValueBe = DOF_value(
		       i, j, k-1, component, advected_level );
	         AdvectorValueBe = AdvectingField->DOF_value(
		     i, j, k-1, component, advecting_level );
	         wb = 0.5 * ( AdvectorValueBe + AdvectorValueC );
	         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
	         else fbe = wb * AdvectedValueC;
	       }
	     }
	     
             flux = (fto - fbo) * dxC * dzC
	     		+ (fri - fle) * dyC * dzC
			+ (ffr - fbe) * dxC * dyC;
             VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
	   }
	 }
       }
     }
   }

   // Synchronize vector for parallel usage
   VEC_rhs->synchronize() ;   

}



// //----------------------------------------------------------------------
// void 
// FV_DiscreteField_Staggered:: assemble_advection_TVD( 
// 	FV_DiscreteField const* AdvectingField,
// 	size_t advecting_level, double const& coef, size_t advected_level,
// 	LA_Vector *VEC_rhs ) const
// //----------------------------------------------------------------------
// {
//    MAC_LABEL( "FV_DiscreteField_Staggered:: assemble_advection_TVD" );   
//    MAC_CHECK_PRE( advected_level < STO_DEPTH ) ;
//    MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
//    MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ; 
//    
//    // Parameters
//    size_t_vector min_unknown_index(DIM,0);
//    size_t_vector max_unknown_index(DIM,0);
//    size_t center_pos_in_matrix = 0, component = 0 ;
//    double xC = 0., yC = 0., zC = 0., xr, xR, xl, xL, yt, yT, yb, yB,
// 	zf, zF, zb, zB;
//    double dxC, dyC, dzC, dxr = 0., dxl, dxCr, dxCl, dxRr, dxR, dxLl,
//    	dyt = 0., dyb, dyCt, dyCb, dyTt, dyT, dyBb,
// 	dzf = 0., dzb, dzCf, dzCb, dzFf, dzF, dzBb;
// 
//    double AdvectedValueC = 0., AdvectedValueRi = 0., AdvectedValueLe = 0.,
//    	AdvectedValueTo = 0., AdvectedValueBo = 0., AdvectedValueFr = 0., 
// 	AdvectedValueBe = 0, AdvectedValueLeLe=0., AdvectedValueRiRi=0., 
// 	AdvectedValueBoBo=0., AdvectedValueToTo=0., AdvectedValueBeBe=0.,
// 	AdvectedValueFrFr=0.,  
//    	AdvectorValueC = 0., AdvectorValueRi = 0., AdvectorValueLe = 0.,
//    	AdvectorValueTo = 0., AdvectorValueBo = 0., AdvectorValueFr = 0., 
// 	AdvectorValueBe = 0, AdvectorValueToLe = 0., AdvectorValueToRi = 0., 
// 	AdvectorValueBoLe = 0., AdvectorValueBoRi = 0., AdvectorValueFrLe = 0., 
// 	AdvectorValueFrRi = 0., AdvectorValueBeLe = 0., AdvectorValueBeRi = 0.,
// 	AdvectorValueFrTo = 0., AdvectorValueFrBo = 0., AdvectorValueBeTo = 0., 
// 	AdvectorValueBeBo = 0.;
//    double ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
// 	fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
//    double cRip12, cLip12, cRim12, cLim12,
//    	thetaC, thetaRi, thetaLe, thetaTo, thetaBo, thetaFr, thetaBe;
//    FV_SHIFT_TRIPLET shift ;
// 
//    // Nullify vector
//    VEC_rhs->nullify();
// 
//    // Assemble vector
//    for (component=0;component<NB_COMPS;++component)
//    {
//      // Get local min and max indices
//      for( size_t l=0; l<DIM; ++l )
//      {
//        min_unknown_index(l) = 
//        	(*min_index_unknown_handled_by_proc)[component](l);
//        max_unknown_index(l) =
//         (*max_index_unknown_handled_by_proc)[component](l);
//      }
//      shift = shift_staggeredToStaggered( component ) ;
// 	     	
//      // Perform assembling
//      for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
//      {          
//        xC = get_DOF_coordinate( i, component, 0 );
//        dxC = get_cell_size( i, component, 0 ) ;      
//        for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
//        {
// 	 yC = get_DOF_coordinate( j, component, 1 );
// 	 dyC = get_cell_size( j, component, 1 ) ; 
//  
//          if ( DIM == 2 )
// 	 {
// 	   size_t k = 0 ;
// 	   center_pos_in_matrix = DOF_global_number( i, j, k, component );
//            AdvectorValueC = AdvectingField->DOF_value( i, j, k, component, 
// 	   	advecting_level );
//            AdvectedValueC = DOF_value( i, j, k, component, advected_level );
// 	   
// 	   // The First component (u)
// 	   if ( component == 0 )
// 	   {	       	     
// 	     // Right and Left
// 	     // --------------
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	     {
// 	       AdvectorValueRi = AdvectorValueC;
// 	       AdvectedValueRi = AdvectedValueC;
// 	     }
// 	     else
// 	     {
// 	       AdvectorValueRi = AdvectingField->DOF_value( i+1, j, k, 
// 	       		component, advecting_level );	     
// 	       AdvectedValueRi = DOF_value(
// 	       		i+1, j, k, component, advected_level );	     
// 	     }
// 	     
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	     {
// 	       AdvectorValueLe = AdvectorValueC;
// 	       AdvectedValueLe = AdvectedValueC;
// 	     }
// 	     else
// 	     {
// 	       AdvectorValueLe = AdvectingField->DOF_value( i-1, j, k, 
// 	       		component, advecting_level );
// 	       AdvectedValueLe = DOF_value(
// 	       		i-1, j, k, component, advected_level );
// 	     }
// 
// 	     thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
// 	   	( AdvectedValueC - AdvectedValueLe ) /
// 		( AdvectedValueRi - AdvectedValueC ) : 1.e20;
// 	     
// 	     // Right (X)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	       fri = AdvectorValueC * AdvectedValueC;
// 	     else
// 	     {	     
//                ur = 0.5 * ( AdvectorValueRi + AdvectorValueC );
// 	       if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT )
// 	       {
//                  if ( ur > 0. ) fri = ur * AdvectedValueC;
//                  else fri = ur * AdvectedValueRi;
// 	       }
// 	       else
// 	       {
// 	         xr = get_DOF_coordinate( i+shift.i, 1, 0 );
// 	         xR = get_DOF_coordinate( i+1, component, 0 );
// 	         dxCr = xr - xC;
// 	         dxr  = xR - xC;
// 		 cLip12 = AdvectedValueC + ( dxCr / dxr )
// 		 	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	   		* ( AdvectedValueRi - AdvectedValueC );
// 			
// 	         dxRr = xR - xr;
// 	         dxR = get_cell_size( i+1, component, 0 );
//                  AdvectedValueRiRi = DOF_value( 
// 	     	       i+2, j, k, component, advected_level );
// 	     
// 	         thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
// 		 	1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
// 			( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
// 	         cRip12 = AdvectedValueRi - ( dxRr / dxR ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaRi)
// 			* ( AdvectedValueRiRi - AdvectedValueRi );
// 	         fri = 0.5 * ( ur * ( cRip12 + cLip12 )
// 		 		- fabs(ur) * ( cRip12 - cLip12 ) );	 
// 	       }
//              }
// 
// 	     // Left (X)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	       fle = AdvectorValueC * AdvectedValueC;
// 	     else
// 	     {	     
// 	       ul = 0.5 * ( AdvectorValueLe + AdvectorValueC );
// 	       if ( DOF_color(i-1, j, k, component ) == FV_BC_LEFT )
// 	       {
// 	         if ( ul > 0. ) fle = ul * AdvectedValueLe;
// 	         else fle = ul * AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         xl = get_DOF_coordinate( i+shift.i-1, 1, 0 );
// 	         xL = get_DOF_coordinate( i-1, component, 0 );
// 	         dxl  = xC - xL;
// 	         dxLl = xl - xL;
// 		 
//                  AdvectedValueLeLe = DOF_value( 
// 	     	       i-2, j, k, component, advected_level );
// 	         
// 		 thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
// 		 	( AdvectedValueLe - AdvectedValueLeLe )/
// 			( AdvectedValueC - AdvectedValueLe ) : 1.e20;
// 	         cLim12 = AdvectedValueLe + ( dxLl / dxl ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaLe)
// 			* ( AdvectedValueC - AdvectedValueLe );
// 	         if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
// 		   cRim12 = AdvectedValueC;
// 	         else
// 		 {
// 	           xR = get_DOF_coordinate( i+1, component, 0 );
// 		   dxr  = xR - xC;
// 		   dxCl = xC - xl;
// 		   
// 	           cRim12 = AdvectedValueC - ( dxCl / dxr ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaC)
// 		 	* ( AdvectedValueRi - AdvectedValueC );
// 		 }
// 	         fle = 0.5 * ( ul * ( cRim12 + cLim12 )
// 	       		- fabs(ul) * ( cRim12 - cLim12 ) );
// 	       }
//              }
// 	 
// 	     // Top and Bottom
// 	     // --------------
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	     {
// 	       AdvectorValueTo = AdvectorValueC;
// 	       AdvectedValueTo = AdvectedValueC;
// 	     }
// 	     else
// 	     {
// 	       AdvectorValueTo = AdvectingField->DOF_value( i, j+1, k, 
// 	       		component, advecting_level );	     
//                AdvectedValueTo = DOF_value( 
// 	     	     i, j+1, k, component, advected_level );
// 	     }
// 
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	     {
// 	       AdvectorValueBo = AdvectorValueC;
// 	       AdvectedValueBo = AdvectedValueC;
// 	     }
// 	     else
//              {
// 	       AdvectorValueBo = AdvectingField->DOF_value( i, j-1, k, 
// 	       	component, advecting_level );
//                AdvectedValueBo = DOF_value( 
// 	     	     i, j-1, k, component, advected_level );
// 	     }
// 
// 	     thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
// 	 	( AdvectedValueC - AdvectedValueBo ) /
// 		( AdvectedValueTo - AdvectedValueC ) : 1.e20;
// 
// 	     // Top (Y)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	       fto = 0.;
// 	     else
// 	     {
// 	       AdvectorValueToLe = AdvectingField->DOF_value( i+shift.i-1, 
// 	       		j+shift.j, k, 1, advecting_level );
// 	       AdvectorValueToRi = AdvectingField->DOF_value( i+shift.i, 
// 	       		j+shift.j, k, 1, advecting_level );
// 	       vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
// 	       if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP
// 	       	|| DOF_color( i, j+1, k, component ) == FV_BC_TOP_LEFT
// 	       	|| DOF_color( i, j+1, k, component ) == FV_BC_TOP_RIGHT )
// 	       {
// 	         if ( vt > 0. ) fto = vt * AdvectedValueC;
// 	         else fto = vt * AdvectedValueTo;
// 	       }
// 	       else
// 	       {
// 	         yt = get_DOF_coordinate( j+shift.j, 1, 1 );
// 		 yT = get_DOF_coordinate( j+1, component, 1 );
// 		 dyCt = yt - yC;
// 		 dyt  = yT - yC;
// 		 
// 		 cLip12 = AdvectedValueC + ( dyCt / dyt ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	 		* ( AdvectedValueTo - AdvectedValueC );
// 	         dyTt = yT - yt;
// 	         dyT = get_cell_size( j+1, component, 1 );
// 	         
//                  AdvectedValueToTo = DOF_value( 
// 	     	       i, j+2, k, component, advected_level );
// 	   
// 	         thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
// 		 	1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
// 			( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
// 		 cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaTo)
// 	       		* ( AdvectedValueToTo - AdvectedValueTo );
// 		 
// 		 fto = 0.5 * ( vt * ( cRip12 + cLip12 )
// 	     		- fabs(vt) * ( cRip12 - cLip12 ) );
// 	       }
// 	     }
// 
// 	     // Bottom (Y)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	       fbo = 0.;
// 	     else
// 	     {
// 	       AdvectorValueBoLe = AdvectingField->DOF_value(
// 	           i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
// 	       AdvectorValueBoRi = AdvectingField->DOF_value( 
// 	           i+shift.i, j+shift.j-1, k, 1, advecting_level );
// 	       vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
// 	       if ( DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
// 	       	|| DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM_LEFT
// 	       	|| DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM_RIGHT )
// 	       {
// 	         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
// 	         else fbo = vb * AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         yb = get_DOF_coordinate( j+shift.j-1, 1, 1 );
// 	         yB = get_DOF_coordinate( j-1, component, 1 );
// 	         dyb  = yC - yB;
// 		 if ( DOF_color( i, j, k, component ) == FV_BC_TOP)
// 		   cRim12 = AdvectedValueC;
// 		 else
// 		 {
// 		   yT = get_DOF_coordinate( j+1, component, 1 );
// 		   dyt  = yT - yC;
// 		   dyCb = yC - yb;
// 	           cRim12 = AdvectedValueC - ( dyCb / dyt ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	     		* ( AdvectedValueTo - AdvectedValueC );
// 	         }
// 		 dyBb = yb - yB;
//                  AdvectedValueBoBo = DOF_value( 
// 	     	       i, j-2, k, component, advected_level );
// 	   
// 	         thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
// 	   		( AdvectedValueBo - AdvectedValueBoBo ) /
// 			( AdvectedValueC - AdvectedValueBo ) : 1.e20;
// 	         cLim12 = AdvectedValueBo + ( dyBb / dyb ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaBo)
// 	       		* ( AdvectedValueC - AdvectedValueBo );
// 	         fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
// 	     		- fabs(vb) * ( cRim12 - cLim12 ) );
// 	       }
// 	     }
// 	   }
// 	   // The second component (v)
// 	   else
// 	   {
// 	     // Right and Left
// 	     // --------------
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	     {
// 	       AdvectorValueRi = AdvectorValueC;
// 	       AdvectedValueRi = AdvectedValueC;
// 	     }
// 	     else
// 	     {
// 	       AdvectorValueRi = AdvectingField->DOF_value( i+1, j, k, 
// 	       		component, advecting_level );	     
//                AdvectedValueRi = DOF_value( 
// 	     	     i+1, j, k, component, advected_level );
// 	     }
// 
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	     {
// 	       AdvectorValueLe = AdvectorValueC;
// 	       AdvectedValueLe = AdvectedValueC;
// 	     }
// 	     else
//              {
// 	       AdvectorValueLe = AdvectingField->DOF_value( i-1, j, k, 
// 	       		component, advecting_level );
//                AdvectedValueLe = DOF_value( 
// 	     	     i-1, j, k, component, advected_level );
// 	     }
// 
// 	     thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
// 	     	( AdvectedValueC - AdvectedValueLe ) /
// 		( AdvectedValueRi - AdvectedValueC ) : 1.e20;
// 
// 	     // Right (X)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	       fri = 0.;
// 	     else
// 	     {	     
// 	       AdvectorValueToRi = AdvectingField->DOF_value( 
// 	           i+shift.i, j+shift.j, k, 0, advecting_level );
// 	       AdvectorValueBoRi = AdvectingField->DOF_value( 
// 	           i+shift.i, j+shift.j-1, k, 0, advecting_level );
// 	       ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
// 	       if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
// 	       	|| DOF_color( i+1, j, k, component ) == FV_BC_BOTTOM_RIGHT
// 	       	|| DOF_color( i+1, j, k, component ) == FV_BC_TOP_RIGHT )
// 	       {
// 	         if ( ur > 0. ) fri = ur * AdvectedValueC;
// 	         else fri = ur * AdvectedValueRi;
// 	       }
// 	       else
// 	       {
// 	         xr = get_DOF_coordinate( i+shift.i, 0, 0 );
// 		 xR = get_DOF_coordinate( i+1, component, 0 );
// 		 dxCr = xr - xC;
// 		 dxr  = xR - xC;
// 		 
// 		 cLip12 = AdvectedValueC + ( dxCr / dxr ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	 	     	* ( AdvectedValueRi - AdvectedValueC );
// 		 
// 		 dxRr = xR - xr;
// 	         dxR = get_cell_size( i+1, component, 0 );
//                  AdvectedValueRiRi = DOF_value( 
// 	     	       i+2, j, k, component, advected_level );
// 	   
// 	         thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
// 		 	1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
// 			( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
// 	         cRip12 = AdvectedValueRi - ( dxRr / dxR ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaRi)
// 	       		* ( AdvectedValueRiRi - AdvectedValueRi );
// 	         fri = 0.5 * ( ur * ( cRip12 + cLip12 )
// 	     		- fabs(ur) * ( cRip12 - cLip12 ) );
// 	       }
// 	     }
// 	         
// 	     // Left (X)
// 	     if ( DOF_color(i, j, k, component ) == FV_BC_LEFT )
// 	       fle = 0.;
// 	     else
// 	     {
// 	       AdvectorValueToLe = AdvectingField->DOF_value(
// 	           i+shift.i-1, j+shift.j, k, 0, advecting_level );
// 	       AdvectorValueBoLe = AdvectingField->DOF_value(
// 	           i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
// 	       ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
// 	       if ( DOF_color( i-1, j, k, component ) == FV_BC_LEFT
// 	       	|| DOF_color( i-1, j, k, component ) == FV_BC_BOTTOM_LEFT
// 	       	|| DOF_color( i-1, j, k, component ) == FV_BC_TOP_LEFT)
// 	       {
// 	         if ( ul > 0. ) fle = ul * AdvectedValueLe;
// 	         else fle = ul * AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         xl = get_DOF_coordinate( i+shift.i-1, 0, 0 );
// 		 xL = get_DOF_coordinate( i-1, component, 0 );
// 		 dxl  = xC - xL;
// 		 if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
// 		   cRim12 = AdvectedValueC;
// 		 else
// 		 {
// 		   xR = get_DOF_coordinate( i+1, component, 0 );
// 		   dxr  = xR - xC;
// 		   dxCl = xC - xl;
// 		   cRim12 = AdvectedValueC
// 		       - ( dxCl / dxr ) * FV_DiscreteField::SuperBee_phi(thetaC)
// 	     	       * ( AdvectedValueRi - AdvectedValueC );
// 		 }
// 	         dxLl = xl - xL;
//                  AdvectedValueLeLe = DOF_value( 
// 	     	       i-2, j, k, component, advected_level );
// 	   
// 	         thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
// 	   		( AdvectedValueLe - AdvectedValueLeLe ) /
// 			( AdvectedValueC - AdvectedValueLe ) : 1.e20;
// 	         cLim12 = AdvectedValueLe + ( dxLl / dxl ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaLe)
// 	       		* ( AdvectedValueC - AdvectedValueLe );
// 	         fle = 0.5 * ( ul * ( cRim12 + cLim12 )
// 	     		- fabs(ul) * ( cRim12 - cLim12 ) );
// 	       }
// 	     }
// 	 
// 	     // Top and Bottom
// 	     // --------------
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	     {
// 	       AdvectorValueTo = AdvectorValueC;
// 	       AdvectedValueTo = AdvectedValueC;
// 	     }
// 	     else
// 	     {
// 	       AdvectorValueTo = AdvectingField->DOF_value( i, j+1, k, 
// 	       		component, advecting_level );	     
//                AdvectedValueTo = DOF_value( 
// 	     	     i, j+1, k, component, advected_level );
// 	     }
// 
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	     {
// 	       AdvectorValueBo = AdvectorValueC;
// 	       AdvectedValueBo = AdvectedValueC;
// 	     }
// 	     else
// 	     {
//                AdvectorValueBo = AdvectingField->DOF_value( i, j-1, k, 
// 	       		component, advecting_level );
//                AdvectedValueBo = DOF_value( 
// 	     	       i, j-1, k, component, advected_level );
// 	     }
// 
// 	     thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
// 	   	( AdvectedValueC - AdvectedValueBo ) /
// 		( AdvectedValueTo - AdvectedValueC ) : 1.e20;
// 	     
// 	     // Top (Y)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	       fto = AdvectorValueC * AdvectedValueC;
// 	     else
// 	     {	     
//                vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
// 	       if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP )
// 	       {
// 	         if ( vt > 0. ) fto = vt * AdvectedValueC;
// 	         else fto = vt * AdvectedValueTo;
// 	       }
// 	       else
// 	       {
// 	         yt = get_DOF_coordinate( j+shift.j, 0, 1 );
// 	         yT = get_DOF_coordinate( j+1, component, 1 );
// 	         dyCt = yt - yC;
// 	         dyt  = yT - yC;
// 		 cLip12 = AdvectedValueC + ( dyCt / dyt ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	   		* ( AdvectedValueTo - AdvectedValueC );
// 
// 	         dyTt = yT - yt;
// 	         dyT = get_cell_size( j+1, component, 1 );
//                  AdvectedValueToTo = DOF_value( 
// 	     	       i, j+2, k, component, advected_level );
// 	     
// 	         thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
// 		 	1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
// 		     	( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
// 	         cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaTo)
// 			* ( AdvectedValueToTo - AdvectedValueTo );
// 	         fto = 0.5 * ( vt * ( cRip12 + cLip12 )
// 	       		- fabs(vt) * ( cRip12 - cLip12 ) );
// 	       }	   
//              }
// 
// 	     // Bottom (Y)
// 	     if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	       fbo = AdvectorValueC * AdvectedValueC;
// 	     else
// 	     {
// 	       vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
// 	       if ( DOF_color(i,j-1,k,component) == FV_BC_BOTTOM )
// 	       {
// 	         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
// 	         else fbo = vb * AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         yb = get_DOF_coordinate( j+shift.j-1, 0, 1 );
// 	         yB = get_DOF_coordinate( j-1, component, 1 );
// 	         dyb  = yC - yB;
// 	       
// 	         dyBb = yb - yB;
//                  AdvectedValueBoBo = DOF_value( 
// 	     	       i, j-2, k, component, advected_level );
// 	     
// 	         thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
// 		 	( AdvectedValueBo - AdvectedValueBoBo ) /
// 			( AdvectedValueC - AdvectedValueBo ) : 1.e20;
// 	         cLim12 = AdvectedValueBo + ( dyBb / dyb ) 
// 		 	* FV_DiscreteField::SuperBee_phi(thetaBo)
// 			* ( AdvectedValueC - AdvectedValueBo );
// 	   
// 	         if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 		   cRim12 = AdvectedValueC;
// 	         else
// 		 {
// 		   yT = get_DOF_coordinate( j+1, component, 1 );
// 		   dyt  = yT - yC;
// 		   dyCb = yC - yb;
// 	           cRim12 = AdvectedValueC - ( dyCb / dyt )
// 		   	* FV_DiscreteField::SuperBee_phi(thetaC)
// 		 	* ( AdvectedValueTo - AdvectedValueC );
// 	         }
// 		 fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
// 	       		- fabs(vb) * ( cRim12 - cLim12 ) );
// 	       }
//              }
// 	   }
// 
//            flux = ( fto - fbo ) * dxC + ( fri - fle ) * dyC;
//            VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
// 	 }
// 	 else // DIM = 3
// 	 {
// 	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
// 	   {
// 	     zC = get_DOF_coordinate( k, component, 2 ) ;
// 	     dzC = get_cell_size( k, component, 2 ) ;	     
// 	     center_pos_in_matrix = DOF_global_number( i, j, k, component );
//              AdvectorValueC = AdvectingField->DOF_value( i, j, k, component, 
// 	     	advecting_level );
//              AdvectedValueC = DOF_value( 
// 	           i, j, k, component, advected_level );
// 
// 	     // The First component (u)
// 	     if ( component == 0 )
// 	     {
// 	       // Right and Left
// 	       // --------------
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	       {
// 	         AdvectorValueRi = AdvectorValueC;
// 	         AdvectedValueRi = AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         AdvectorValueRi = AdvectingField->DOF_value( i+1, j, k, 
// 		 	component, advecting_level );
//                  AdvectedValueRi = DOF_value( 
// 	     	       	i+1, j, k, component, advected_level );
// 	       }
// 
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	       {
// 	         AdvectorValueLe = AdvectorValueC;
// 	         AdvectedValueLe = AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         AdvectorValueLe = AdvectingField->DOF_value( i-1, j, k, 
// 		 	component, advecting_level );
//                  AdvectedValueLe = DOF_value( 
// 	     	       i-1, j, k, component, advected_level );
// 	      }
// 
// 	       thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
// 	       		( AdvectedValueC - AdvectedValueLe ) /
// 			( AdvectedValueRi - AdvectedValueC ) : 1.e20;
// 	     
// 	       // Right (X)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	         fri = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {	     
// 	         ur = 0.5 * ( AdvectorValueRi + AdvectorValueC );
// 	         if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT )
// 		 {
//                    if ( ur > 0. ) fri = ur * AdvectedValueC;
//                    else fri = ur * AdvectedValueRi;
// 		 }
// 	         else
// 	         {		 
// 		   xr = get_DOF_coordinate( i+shift.i, 1, 0 );
// 	           xR = get_DOF_coordinate( i+1, component, 0 );
// 	           dxCr = xr - xC;
// 	           dxr  = xR - xC;
// 		   cLip12 = AdvectedValueC + ( dxCr / dxr ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	   		* ( AdvectedValueRi - AdvectedValueC );
// 			
// 	           dxRr = xR - xr;
// 	           dxR = get_cell_size( i+1, component, 0 );
//                    AdvectedValueRiRi = DOF_value( 
// 	     	         i+2, j, k, component, advected_level );
// 	     
// 	           thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
// 		   	1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
// 			( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20 ;
// 	           cRip12 = AdvectedValueRi - ( dxRr / dxR ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaRi)
// 			* ( AdvectedValueRiRi - AdvectedValueRi );
// 	           fri = 0.5 * ( ur * ( cRip12 + cLip12 )
// 	       		- fabs(ur) * ( cRip12 - cLip12 ) );
// 		}
//                }
// 
// 	       // Left (X)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	         fle = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {	     
// 	         ul = 0.5 * ( AdvectorValueLe + AdvectorValueC );
// 	         if ( DOF_color( i-1, j, k, component ) == FV_BC_LEFT )
// 		 {
// 	           if ( ul > 0. ) fle = ul * AdvectedValueLe;
// 	           else fle = ul * AdvectedValueC;
// 		 }
// 	         else
// 	         {
// 		   xl = get_DOF_coordinate( i+shift.i-1, 1, 0 );
// 		   xL = get_DOF_coordinate( i-1, component, 0 );
// 		   dxl  = xC - xL;		   
// 		   dxLl = xl - xL;
//                    AdvectedValueLeLe = DOF_value( 
// 	     	       i-2, j, k, component, advected_level );
// 
// 	           thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
// 		 	( AdvectedValueLe - AdvectedValueLeLe ) /
// 			( AdvectedValueC - AdvectedValueLe ) : 1.e20;
// 	           cLim12 = AdvectedValueLe + ( dxLl / dxl ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaLe)
// 			* ( AdvectedValueC - AdvectedValueLe );
// 	   
// 	           if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
// 		     cRim12 = AdvectedValueC;
// 	           else
// 		   {
// 	             xR = get_DOF_coordinate( i+1, component, 0 );
// 		     dxr  = xR - xC;
// 		     dxCl = xC - xl;
// 		     cRim12 = AdvectedValueC - ( dxCl / dxr )
// 			* FV_DiscreteField::SuperBee_phi(thetaC)
// 			* ( AdvectedValueRi - AdvectedValueC );
// 	           }
// 	           fle = 0.5 * ( ul * ( cRim12 + cLim12 )
// 	       		- fabs(ul) * ( cRim12 - cLim12 ) );
// 	         }
//                }
// 
// 	 
// 	       // Top and Bottom
// 	       // --------------
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	       {
// 	         AdvectorValueTo = AdvectorValueC;
// 	         AdvectedValueTo = AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         AdvectorValueTo = AdvectingField->DOF_value( i, j+1, k, 
// 		 	component, advecting_level );	     
//                  AdvectedValueTo = DOF_value( 
// 	     	       	i, j+1, k, component, advected_level );
//                }
// 	       
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	       {
// 	         AdvectorValueBo = AdvectorValueC;
// 	         AdvectedValueBo = AdvectedValueC;
// 	       }
// 	       else
//                {
// 	         AdvectorValueBo = AdvectingField->DOF_value( i, j-1, k, 
// 		 	component, advecting_level );
//                  AdvectedValueBo = DOF_value( 
// 	     	       	i, j-1, k, component, advected_level );
//                }
// 	       
// 	       thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
// 			( AdvectedValueC - AdvectedValueBo ) /
// 			( AdvectedValueTo - AdvectedValueC ) : 1.e20;
// 
// 	       // Top (Y)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	         fto = 0.;
// 	       else
// 	       {
//                  AdvectorValueToLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j+shift.j, k, 1, advecting_level );
// 	         AdvectorValueToRi = AdvectingField->DOF_value( 
// 		     i+shift.i, j+shift.j, k, 1, advecting_level );
// 	         vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
// 	         if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP
// 	       	  	|| DOF_color( i, j+1, k, component ) == FV_BC_TOP_LEFT
// 	       	  	|| DOF_color( i, j+1, k, component ) 
// 				== FV_BC_TOP_RIGHT )
// 	         {
// 	           if ( vt > 0. ) fto = vt * AdvectedValueC;
// 	           else fto = vt * AdvectedValueTo;
// 	         }
// 	         else
// 	         {
// 	           yt = get_DOF_coordinate( j+shift.j, 1, 1 );
// 		   yT = get_DOF_coordinate( j+1, component, 1 );
// 		   dyCt = yt - yC;
// 		   dyt  = yT - yC;
// 		   cLip12 = AdvectedValueC + ( dyCt / dyt ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	 		* ( AdvectedValueTo - AdvectedValueC );
// 	           dyTt = yT - yt;
// 	           dyT = get_cell_size( j+1, component, 1 );
//                    AdvectedValueToTo = DOF_value( 
// 	     	         i, j+2, k, component, advected_level );
// 	   
// 	           thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
// 		   	1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
// 			( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
// 	           cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaTo)
// 	       		* ( AdvectedValueToTo - AdvectedValueTo );
// 		   fto = 0.5 * ( vt * ( cRip12 + cLip12 )
// 	     		- fabs(vt) * ( cRip12 - cLip12 ) );
// 	         }
// 	       }
// 
// 	       // Bottom (Y)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	         fbo = 0.;
// 	       else
// 	       {	     
// 	         AdvectorValueBoLe = AdvectingField->DOF_value( 
// 		     i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
// 	         AdvectorValueBoRi = AdvectingField->DOF_value( 
// 		     i+shift.i, j+shift.j-1, k,	1, advecting_level );
// 	         vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
// 	         if ( DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
// 	       	  	|| DOF_color( i, j-1, k, component ) 
// 				== FV_BC_BOTTOM_LEFT
// 	          	|| DOF_color( i, j-1, k, component ) 
// 				== FV_BC_BOTTOM_RIGHT )
// 	         {
// 	           if ( vb > 0. ) fbo = vb * AdvectedValueBo;
// 	           else fbo = vb * AdvectedValueC;
// 	         }
// 	         else
// 	         {
// 	           yb = get_DOF_coordinate( j+shift.j-1, 1, 1 );
// 	           yB = get_DOF_coordinate( j-1, component, 1 );
// 	           dyb  = yC - yB;
// 		   if ( DOF_color( i, j, k, component ) == FV_BC_TOP ) 
// 		     cRim12 = AdvectedValueC;
// 		   else
// 		   {
// 		     yT = get_DOF_coordinate( j+1, component, 1 );
// 		     dyt  = yT - yC;
// 		     dyCb = yC - yb;
// 		     cRim12 = AdvectedValueC - ( dyCb / dyt )
// 		     	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	     		* ( AdvectedValueTo - AdvectedValueC );
// 	           }
// 		   dyBb = yb - yB;
//                    AdvectedValueBoBo = DOF_value( 
// 	     	         i, j-2, k, component, advected_level );
// 	   
// 	           thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
// 	   		( AdvectedValueBo - AdvectedValueBoBo ) /
// 			( AdvectedValueC - AdvectedValueBo ) : 1.e20;
// 	           cLim12 = AdvectedValueBo + ( dyBb / dyb )
// 		   	* FV_DiscreteField::SuperBee_phi(thetaBo)
// 	       		* ( AdvectedValueC - AdvectedValueBo );
// 	           fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
// 	     		- fabs(vb) * ( cRim12 - cLim12 ) );
// 	         }
// 	       }
// 
// 	       // Front and Behind
// 	       // ----------------
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
// 	       {
// 	         AdvectorValueFr = AdvectorValueC;
// 	         AdvectedValueFr = AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         AdvectorValueFr = AdvectingField->DOF_value( i, j, k+1, 
// 		 	component, advecting_level );	     
//                  AdvectedValueFr = DOF_value( 
// 	     	       	i, j, k+1, component, advected_level );
//                }
// 	       
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
// 	       {
// 	         AdvectorValueBe = AdvectorValueC;
// 	         AdvectedValueBe = AdvectedValueC;
// 	       }
// 	       else
//                {
// 	         AdvectorValueBe = AdvectingField->DOF_value( i, j, k-1, 
// 		 	component, advecting_level );
//                  AdvectedValueBe = DOF_value( 
// 	     	       	i, j, k-1, component, advected_level );
//                }
// 	       
// 	       thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ? 
// 	       		( AdvectedValueC - AdvectedValueBe ) /
// 			( AdvectedValueFr - AdvectedValueC ) : 1.e20;
// 	       
// 	       // Front (Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
// 	         ffr = 0.;
// 	       else
// 	       {
// 	         AdvectorValueFrLe = AdvectingField->DOF_value( 
// 		     i+shift.i-1, j, k+shift.k, 2, advecting_level );
// 	         AdvectorValueFrRi = AdvectingField->DOF_value( 
// 		     i+shift.i, j, k+shift.k, 2, advecting_level );
// 	         wf = 0.5 * (AdvectorValueFrLe + AdvectorValueFrRi);
// 	         if ( DOF_color( i, j, k+1, component ) == FV_BC_FRONT
// 	       		|| DOF_color( i, j, k+1, component ) 
// 				== FV_BC_FRONT_LEFT 
// 	       		|| DOF_color( i, j, k+1, component ) 
// 				== FV_BC_FRONT_RIGHT )
// 		 {
// 		   if ( wf > 0. ) ffr = wf * AdvectedValueC;
// 	           else ffr = wf * AdvectedValueFr;
// 		 }
// 		 else
// 		 {
// 		   zf = get_DOF_coordinate( k+shift.k, 2, 2 );
// 		   zF = get_DOF_coordinate( k+1, component, 2 );
// 		   dzCf = zf - zC;
// 		   dzf  = zF - zC;
// 		   cLip12 = AdvectedValueC + ( dzCf / dzf )
// 		   	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	 		* ( AdvectedValueFr - AdvectedValueC );
// 		   dzFf = zF - zf;
// 		   dzF = get_cell_size( k+1, component, 2 );
//                    AdvectedValueFrFr = DOF_value( 
// 	     	         i, j, k+2, component, advected_level );
// 
// 		   thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 
// 		   	1.e-20 ? ( AdvectedValueFr - AdvectedValueC ) /
// 			( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
// 	           cRip12 = AdvectedValueFr - ( dzFf / dzF )
// 		 	* FV_DiscreteField::SuperBee_phi(thetaFr)
// 	       		* ( AdvectedValueFrFr - AdvectedValueFr );
// 		   ffr = 0.5 * ( wf * ( cRip12 + cLip12 )
// 	     		- fabs(wf) * ( cRip12 - cLip12 ) );
// 		 }
// 	       }
// 
// 	       // Behind (Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
// 	         fbe = 0.;
// 	       else
// 	       {	     
// 	         AdvectorValueBeLe = AdvectingField->DOF_value( 
// 		     i+shift.i-1, j, k+shift.k-1, 2, advecting_level );
// 	         AdvectorValueBeRi = AdvectingField->DOF_value(
// 		     i+shift.i, j, k+shift.k-1, 2, advecting_level );
// 	         wb = 0.5 * ( AdvectorValueBeLe + AdvectorValueBeRi );
// 		 if ( DOF_color( i, j, k-1, component ) == FV_BC_BEHIND
// 		 	|| DOF_color( i, j, k-1, component ) 
// 				== FV_BC_BEHIND_LEFT
// 		 	|| DOF_color( i, j, k-1, component ) 
// 				== FV_BC_BEHIND_RIGHT )
// 		 {
// 	           if ( wb > 0. ) fbe = wb * AdvectedValueBe;
// 	           else fbe = wb * AdvectedValueC;
// 		 }
// 		 else
// 		 {
// 		   zb = get_DOF_coordinate( k+shift.k-1, 2, 2 );
// 		   zB = get_DOF_coordinate( k-1, component, 2 );
// 		   dzb  = zC - zB;
// 		   if ( DOF_color( i, j, k, component ) == FV_BC_FRONT ) 
// 		     cRim12 = AdvectedValueC;
// 		   else
// 		   {
// 		     zF = get_DOF_coordinate( k+1, component, 2 );
// 		     dzf  = zF - zC;
// 		     dzCb = zC - zb;
// 		     cRim12 = AdvectedValueC - ( dzCb / dzf )
// 			* FV_DiscreteField::SuperBee_phi(thetaC)
// 			* ( AdvectedValueFr - AdvectedValueC );
// 		   }
// 		   dzBb = zb - zB;
//                    AdvectedValueBeBe = DOF_value( 
// 	     	         i, j, k-2, component, advected_level );
// 	   
// 	           thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
// 	   		( AdvectedValueBe - AdvectedValueBeBe ) /
// 			( AdvectedValueC - AdvectedValueBe ) : 1.e20;
// 	           cLim12 = AdvectedValueBe + ( dzBb / dzb )
// 		 	* FV_DiscreteField::SuperBee_phi(thetaBe)
// 	       		* ( AdvectedValueC - AdvectedValueBe );
// 		   fbe = 0.5 * ( wb * ( cRim12 + cLim12 )
// 	     		- fabs(wb) * ( cRim12 - cLim12 ) );
// 		 }
// 	       }
// 	     }
// 	     // The Second component (v)
// 	     else if ( component == 1 )
// 	     {
// 	       // Right and Left
// 	       // --------------
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	       {
// 	         AdvectorValueRi = AdvectorValueC;
// 	         AdvectedValueRi = AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         AdvectorValueRi = AdvectingField->DOF_value( i+1, j, k, 
// 		 	component, advecting_level );	     
//                  AdvectedValueRi = DOF_value( 
// 	     	      	i+1, j, k, component, advected_level );
//                }
// 	       
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	       {
// 	         AdvectorValueLe = AdvectorValueC;
// 	         AdvectedValueLe = AdvectedValueC;
// 	       }
// 	       else
// 	       {
//                  AdvectorValueLe = AdvectingField->DOF_value( i-1, j, k, 
// 		 	component, advecting_level );
//                  AdvectedValueLe = DOF_value( 
// 	     	       	i-1, j, k, component, advected_level );
// 	       }
// 
// 	       thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
// 	       		( AdvectedValueC - AdvectedValueLe ) /
// 			( AdvectedValueRi - AdvectedValueC ) : 1.e20;
// 
// 	       // Right (X)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	         fri = 0.;
// 	       else
// 	       {
// 	         AdvectorValueToRi = AdvectingField->DOF_value( 
// 		     i+shift.i, j+shift.j, k, 0, advecting_level );
// 	         AdvectorValueBoRi = AdvectingField->DOF_value(
// 		     i+shift.i, j+shift.j-1, k,	0, advecting_level );
// 	         ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
// 	         if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
// 	       		|| DOF_color( i+1, j, k, component ) 
// 				== FV_BC_BOTTOM_RIGHT
// 	       		|| DOF_color( i+1, j, k, component ) 
// 				== FV_BC_TOP_RIGHT )
// 		 {
// 		   if ( ur > 0. ) fri = ur * AdvectedValueC;
// 		   else fri = ur * AdvectedValueRi;
// 		 }
// 		 else
// 		 {
// 		   xr = get_DOF_coordinate( i+shift.i, 0, 0 );
// 		   xR = get_DOF_coordinate( i+1, component, 0 );
// 		   dxCr = xr - xC;
// 		   dxr  = xR - xC;
// 		   cLip12 = AdvectedValueC + ( dxCr / dxr ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	 		* ( AdvectedValueRi - AdvectedValueC );
// 		   
// 		   dxRr = xR - xr;
// 	           dxR = get_cell_size( i+1, component, 0 );
//                    AdvectedValueRiRi = DOF_value( 
// 	     	         i+2, j, k, component, advected_level );
// 	   
// 	           thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
// 		   	1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
// 			( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
// 	           cRip12 = AdvectedValueRi - ( dxRr / dxR )
// 			* FV_DiscreteField::SuperBee_phi(thetaRi)
// 			* ( AdvectedValueRiRi - AdvectedValueRi );
// 		   fri = 0.5 * ( ur * ( cRip12 + cLip12 )
// 	     		- fabs(ur) * ( cRip12 - cLip12 ) );
// 		 }
// 	       }
// 	       
// 	       // Left (X)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	         fle = 0.;
// 	       else
// 	       {
// 	         AdvectorValueToLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j+shift.j, k, 0, advecting_level );
// 	         AdvectorValueBoLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
// 	         ul = 0.5 * (AdvectorValueToLe + AdvectorValueBoLe);
// 	         if ( DOF_color( i-1, j, k, component ) == FV_BC_LEFT
// 	       		|| DOF_color( i-1, j, k, component ) 
// 				== FV_BC_BOTTOM_LEFT
// 	       		|| DOF_color( i-1, j, k, component ) 
// 				== FV_BC_TOP_LEFT )
// 		 {
// 	           if ( ul > 0. ) fle = ul * AdvectedValueLe;
// 	           else fle = ul * AdvectedValueC;
// 		 }
// 		 else
// 		 {
// 		   xl = get_DOF_coordinate( i+shift.i-1, 0, 0 );
// 		   xL = get_DOF_coordinate( i-1, component, 0 );
// 		   dxl  = xC - xL;
// 		   if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
// 		     cRim12 = AdvectedValueC;
// 		   else
// 		   {
// 		     xR = get_DOF_coordinate( i+1, component, 0 );
// 		     dxr  = xR - xC;
// 		     dxCl = xC - xl;
// 		     cRim12 = AdvectedValueC - ( dxCl / dxr )
// 		     	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	     		* ( AdvectedValueRi - AdvectedValueC );
// 		   }
// 		   
// 		   dxLl = xl - xL;
//                    AdvectedValueLeLe = DOF_value( 
// 	     	         i-2, j, k, component, advected_level );
// 
// 	           thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
// 	   		( AdvectedValueLe - AdvectedValueLeLe ) /
// 			( AdvectedValueC - AdvectedValueLe ) : 1.e20;
// 	           cLim12 = AdvectedValueLe + ( dxLl / dxl )
// 			* FV_DiscreteField::SuperBee_phi(thetaLe)
// 			* ( AdvectedValueC - AdvectedValueLe );
// 		   fle = 0.5 * ( ul * ( cRim12 + cLim12 )
// 	     		- fabs(ul) * ( cRim12 - cLim12 ) );
// 		 }
// 	       }
// 
// 
// 	       // Top and Bottom
// 	       // --------------
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	       {
// 	         AdvectorValueTo = AdvectorValueC;
// 	         AdvectedValueTo = AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         AdvectorValueTo = AdvectingField->DOF_value( i, j+1, k, 
// 		 	component, advecting_level );	     
//                  AdvectedValueTo = DOF_value( 
// 	     	       	i, j+1, k, component, advected_level );
// 	       }
// 
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	       {
// 	         AdvectorValueBo = AdvectorValueC;
// 	         AdvectedValueBo = AdvectedValueC;
// 	       }
// 	       else
// 	       {
//                  AdvectorValueBo = AdvectingField->DOF_value( i, j-1, k, 
// 		 	component, advecting_level );
//                  AdvectedValueBo = DOF_value( 
// 	     	       	i, j-1, k, component, advected_level );
//                }
// 	       thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
// 	       		( AdvectedValueC - AdvectedValueBo ) /
// 			( AdvectedValueTo - AdvectedValueC ) : 1.e20;
// 	     
// 	       // Top (Y)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	         fto = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {	     
// 	         vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
// 	         if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP )
// 		 {
// 	           if ( vt > 0. ) fto = vt * AdvectedValueC;
// 	           else fto = vt * AdvectedValueTo;
// 		 }
// 	         else
// 	         {
// 		   yt = get_DOF_coordinate( j+shift.j, 0, 1 );
// 	           yT = get_DOF_coordinate( j+1, component, 1 );
// 	           dyCt = yt - yC;
// 	           dyt  = yT - yC;
// 		   cLip12 = AdvectedValueC + ( dyCt / dyt )
// 		   	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	   		* ( AdvectedValueTo - AdvectedValueC );
// 			
// 	           dyTt = yT - yt;
// 	           dyT = get_cell_size( j+1, component, 1 );
//                    AdvectedValueToTo = DOF_value( 
// 	     	         i, j+2, k, component, advected_level );
// 	     
// 	           thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
// 		   	1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
// 		       ( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
// 	           cRip12 = AdvectedValueTo - ( dyTt / dyT )
// 			* FV_DiscreteField::SuperBee_phi(thetaTo)
// 			* ( AdvectedValueToTo - AdvectedValueTo );
// 	           fto = 0.5 * ( vt * ( cRip12 + cLip12 )
// 	       		- fabs(vt) * ( cRip12 - cLip12 ) );
// 		 }
//                }
// 	 
// 	       // Bottom (Y)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	         fbo = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {
// 	         vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
// 		 if ( DOF_color(i, j-1, k, component ) == FV_BC_BOTTOM )
// 		 {
// 		   if ( vb > 0. ) fbo = vb * AdvectedValueBo;
// 	           else fbo = vb * AdvectedValueC;
// 		 }
// 	         else
// 	         {
// 	           yb = get_DOF_coordinate( j+shift.j-1, 0, 1 );
// 	           yB = get_DOF_coordinate( j-1, component, 1 );
// 	           dyb  = yC - yB;
// 		   if ( DOF_color( i, j, k, component ) == FV_BC_TOP ) 
// 		     cRim12 = AdvectedValueC;
// 	           else
// 		   {
// 	             yT = get_DOF_coordinate( j+1, component, 1 );
// 		     dyt  = yT - yC;
// 		     dyCb = yC - yb;
// 	             cRim12 = AdvectedValueC - ( dyCb / dyt )
// 			* FV_DiscreteField::SuperBee_phi(thetaC)
// 			* ( AdvectedValueTo - AdvectedValueC );
// 	           }
// 		   dyBb = yb - yB;
//                    AdvectedValueBoBo = DOF_value( 
// 	     	         i, j-2, k, component, advected_level );
// 	     
// 	           thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
// 		 	( AdvectedValueBo - AdvectedValueBoBo ) /
// 			( AdvectedValueC - AdvectedValueBo ) : 1.e20;
// 	           cLim12 = AdvectedValueBo + ( dyBb / dyb )
// 			* FV_DiscreteField::SuperBee_phi(thetaBo)
// 			* ( AdvectedValueC - AdvectedValueBo );
// 	           fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
// 	       		- fabs(vb) * ( cRim12 - cLim12 ) );
// 	         }
//                }
// 
// 
// 	       // Front and Behind
// 	       // ----------------
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
// 	       {
// 	         AdvectorValueFr = AdvectorValueC;
// 	         AdvectedValueFr = AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         AdvectorValueFr = AdvectingField->DOF_value( i, j, k+1, 
// 		 	component, advecting_level );	     
//                  AdvectedValueFr = DOF_value( 
// 	     	       	i, j, k+1, component, advected_level );
// 	       }
// 
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
// 	       {
// 	         AdvectorValueBe = AdvectorValueC;
// 	         AdvectedValueBe = AdvectedValueC;
// 	       }
// 	       else
// 	       {
//                  AdvectorValueBe = AdvectingField->DOF_value( i, j, k-1, 
// 		 	component, advecting_level );
//                  AdvectedValueBe = DOF_value( 
// 	     	       	i, j, k-1, component, advected_level );
// 	       }
// 
// 	       thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ? 
// 	       		( AdvectedValueC - AdvectedValueBe ) /
// 			( AdvectedValueFr - AdvectedValueC ) : 1.e20;
// 	       
// 	       // Front (Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
// 	         ffr = 0.;
// 	       else
// 	       {
//                  AdvectorValueFrBo = AdvectingField->DOF_value( 
// 		     i, j+shift.j-1, k+shift.k, 2, advecting_level );
// 	         AdvectorValueFrTo = AdvectingField->DOF_value( 
// 		     i, j+shift.j, k+shift.k, 2, advecting_level );
//                  wf = 0.5 * ( AdvectorValueFrBo + AdvectorValueFrTo );
// 	         if ( DOF_color( i, j, k+1, component ) == FV_BC_FRONT
// 	       		|| DOF_color( i, j, k+1, component ) 
// 				== FV_BC_FRONT_BOTTOM
// 	       		|| DOF_color( i, j, k+1, component ) 
// 				== FV_BC_FRONT_TOP )
// 		 {
// 		   if ( wf > 0. ) ffr = wf * AdvectedValueC;
// 	           else ffr = wf * AdvectedValueFr;
// 		 }
// 		 else
// 		 {
// 		   zf = get_DOF_coordinate( k+shift.k, 2, 2 );
// 		   zF = get_DOF_coordinate( k+1, component, 2 );
// 		   dzCf = zf - zC;
// 		   dzf  = zF - zC;
// 		   cLip12 = AdvectedValueC + ( dzCf / dzf )
// 			* FV_DiscreteField::SuperBee_phi(thetaC)
// 			* ( AdvectedValueFr - AdvectedValueC );
// 		   dzFf = zF - zf;
// 		   dzF = get_cell_size( k+1, component, 2 );
//                    AdvectedValueFrFr = DOF_value( 
// 	     	         i, j, k+2, component, advected_level );
// 
// 		   thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 
// 		   	1.e-20 ? ( AdvectedValueFr - AdvectedValueC ) /
// 			( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
// 	           cRip12 = AdvectedValueFr - ( dzFf / dzF )
// 			* FV_DiscreteField::SuperBee_phi(thetaFr)
// 			* ( AdvectedValueFrFr - AdvectedValueFr );
// 		   ffr = 0.5 * ( wf * ( cRip12 + cLip12 )
// 	     		- fabs(wf) * ( cRip12 - cLip12 ) );
// 		 }
// 	       }
// 
// 	       // Behind (Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
// 	         fbe = 0.;
// 	       else
// 	       {	     
// 	         AdvectorValueBeBo = AdvectingField->DOF_value( 
// 		     i, j+shift.j-1, k+shift.k-1, 2, advecting_level );
// 	         AdvectorValueBeTo = AdvectingField->DOF_value(
// 		     i, j+shift.j, k+shift.k-1, 2, advecting_level );
// 	         wb = 0.5 * ( AdvectorValueBeBo + AdvectorValueBeTo );
// 		 if ( DOF_color( i, j, k-1, component ) == FV_BC_BEHIND
// 	       		|| DOF_color( i, j, k-1, component ) 
// 				== FV_BC_BEHIND_BOTTOM
// 	       		|| DOF_color( i, j, k-1, component ) 
// 				== FV_BC_BEHIND_TOP )
// 		 {
// 	           if ( wb > 0. ) fbe = wb * AdvectedValueBe;
// 	           else fbe = wb * AdvectedValueC;
// 		 }
// 		 else
// 		 {
// 		   zb = get_DOF_coordinate( k+shift.k-1, 2, 2 );
// 		   zB = get_DOF_coordinate( k-1, component, 2 );
// 		   dzb  = zC - zB;
// 		   if ( DOF_color( i, j, k, component ) == FV_BC_FRONT ) 
// 		     cRim12 = AdvectedValueC;
// 		   else
// 		   {
// 		     zF = get_DOF_coordinate( k+1, component, 2 );
// 		     dzf  = zF - zC;
// 		     dzCb = zC - zb;
// 		     cRim12 = AdvectedValueC - ( dzCb / dzf )
// 			* FV_DiscreteField::SuperBee_phi(thetaC)
// 			* ( AdvectedValueFr - AdvectedValueC );
// 		   }
// 		   dzBb = zb - zB;
//                    AdvectedValueBeBe = DOF_value( 
// 	     	         i, j, k-2, component, advected_level );
// 	   
// 	           thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
// 	   		( AdvectedValueBe - AdvectedValueBeBe ) /
// 			( AdvectedValueC - AdvectedValueBe ) : 1.e20;
// 	           cLim12 = AdvectedValueBe + ( dzBb / dzb )
// 			* FV_DiscreteField::SuperBee_phi(thetaBe)
// 			* ( AdvectedValueC - AdvectedValueBe );
// 		   fbe = 0.5 * ( wb * ( cRim12 + cLim12 )
// 	     		- fabs(wb) * ( cRim12 - cLim12 ) );
// 		 }
// 	       }
// 	     }
// 	     // The Third component (w)
// 	     else
// 	     {
// 	       // Right and Left
// 	       // --------------
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	       {
// 	         AdvectorValueRi = AdvectorValueC;
// 	         AdvectedValueRi = AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         AdvectorValueRi = AdvectingField->DOF_value( i+1, j, k, 
// 		 	component, advecting_level );
//                  AdvectedValueRi = DOF_value( 
// 	     	       	i+1, j, k, component, advected_level );
// 	       }
// 
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	       {
// 	         AdvectorValueLe = AdvectorValueC;
// 	         AdvectedValueLe = AdvectedValueC;
// 	       }
// 	       else
// 	       {
//                  AdvectorValueLe = AdvectingField->DOF_value( i-1, j, k, 
// 		 	component, advecting_level );
//                  AdvectedValueLe = DOF_value( 
// 	     	       	i-1, j, k, component, advected_level );
// 	       }
// 
// 	       thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
// 	       		( AdvectedValueC - AdvectedValueLe ) /
// 			( AdvectedValueRi - AdvectedValueC ) : 1.e20;
// 
// 	       // Right (X)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	         fri = 0.;
// 	       else
// 	       {
// 	         AdvectorValueFrRi = AdvectingField->DOF_value(
// 		     i+shift.i, j, k+shift.k, 0, advecting_level );
// 	         AdvectorValueBeRi = AdvectingField->DOF_value( 
// 		     i+shift.i, j, k+shift.k-1,	0, advecting_level );
// 	         ur = 0.5 * (AdvectorValueFrRi + AdvectorValueBeRi);
// 	         if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
// 	       		|| DOF_color( i+1, j, k, component ) 
// 				== FV_BC_BEHIND_RIGHT
// 	       		|| DOF_color( i+1, j, k, component ) 
// 				== FV_BC_FRONT_RIGHT )
// 		 {
// 		   if ( ur > 0. ) fri = ur * AdvectedValueC;
// 		   else fri = ur * AdvectedValueRi;
// 		 }
// 		 else
// 		 {
// 		   xr = get_DOF_coordinate( i+shift.i, 0, 0 );
// 		   xR = get_DOF_coordinate( i+1, component, 0 );
// 		   dxCr = xr - xC;
// 		   dxr  = xR - xC;
// 		   cLip12 = AdvectedValueC + ( dxCr / dxr ) 
// 		       * FV_DiscreteField::SuperBee_phi(thetaC)
// 	 	       * ( AdvectedValueRi - AdvectedValueC );
// 		   
// 		   dxRr = xR - xr;
// 	           dxR = get_cell_size( i+1, component, 0 );
//                    AdvectedValueRiRi = DOF_value( 
// 	     	         i+2, j, k, component, advected_level );
// 	   
// 	           thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
// 		   	1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
// 			( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
// 	           cRip12 = AdvectedValueRi - ( dxRr / dxR )
// 			* FV_DiscreteField::SuperBee_phi(thetaRi)
// 			* ( AdvectedValueRiRi - AdvectedValueRi );
// 		   fri = 0.5 * ( ur * ( cRip12 + cLip12 )
// 	     		- fabs(ur) * ( cRip12 - cLip12 ) );
// 		 }
// 	       }
// 	       
// 	       // Left (X)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	         fle = 0.;
// 	       else
// 	       {
// 	         AdvectorValueFrLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j, k+shift.k, 0, advecting_level );
// 	         AdvectorValueBeLe = AdvectingField->DOF_value(
// 		     i+shift.i-1, j, k+shift.k-1, 0, advecting_level );
// 	         ul = 0.5 * ( AdvectorValueFrLe + AdvectorValueBeLe );
// 	         if ( DOF_color( i-1, j, k, component ) == FV_BC_LEFT
// 	       		|| DOF_color( i-1, j, k, component ) 
// 				== FV_BC_BEHIND_LEFT
// 	       		|| DOF_color( i-1, j, k, component ) 
// 				== FV_BC_FRONT_LEFT )
// 		 {
// 	           if ( ul > 0. ) fle = ul * AdvectedValueLe;
// 	           else fle = ul * AdvectedValueC;
// 		 }
// 		 else
// 		 {
// 		   xl = get_DOF_coordinate( i+shift.i-1, 0, 0 );
// 		   xL = get_DOF_coordinate( i-1, component, 0 );
// 		   dxl  = xC - xL;
// 		   if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
// 		     cRim12 = AdvectedValueC;
// 		   else
// 		   {
// 		     xR = get_DOF_coordinate( i+1, component, 0 );
// 		     dxr  = xR - xC;
// 		     dxCl = xC - xl;
// 		     cRim12 = AdvectedValueC - ( dxCl / dxr )
// 			* FV_DiscreteField::SuperBee_phi(thetaC)
// 			* ( AdvectedValueRi - AdvectedValueC );
// 		   }
// 		   
// 		   dxLl = xl - xL;
//                    AdvectedValueLeLe = DOF_value( 
// 	     	         i-2, j, k, component, advected_level );
// 
// 	           thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
// 	   		( AdvectedValueLe - AdvectedValueLeLe ) /
// 			( AdvectedValueC - AdvectedValueLe ) : 1.e20;
// 	           cLim12 = AdvectedValueLe + ( dxLl / dxl )
// 			* FV_DiscreteField::SuperBee_phi(thetaLe)
// 			* ( AdvectedValueC - AdvectedValueLe );
// 		   fle = 0.5 * ( ul * ( cRim12 + cLim12 )
// 	     		- fabs(ul) * ( cRim12 - cLim12 ) );
// 		 }
// 	       }
// 
// 
// 	       // Top and Bottom
// 	       // --------------
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	       {
// 	         AdvectorValueTo = AdvectorValueC;
// 	         AdvectedValueTo = AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         AdvectorValueTo = AdvectingField->DOF_value( i, j+1, k, 
// 		 	component, advecting_level );
//                  AdvectedValueTo = DOF_value( 
// 	     	       	i, j+1, k, component, advected_level );
// 	       }
// 
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	       {
// 	         AdvectorValueBo = AdvectorValueC;
// 	         AdvectedValueBo = AdvectedValueC;
// 	       }
// 	       else
// 	       {
//                  AdvectorValueBo = AdvectingField->DOF_value( i, j-1, k, 
// 		 	component, advecting_level );
//                  AdvectedValueBo = DOF_value( 
// 	     	       	i, j-1, k, component, advected_level );
// 	       }
// 
// 	       thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
// 	       		( AdvectedValueC - AdvectedValueBo ) /
// 			( AdvectedValueTo - AdvectedValueC ) : 1.e20;
// 
// 	       // Top (Y)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	         fto = 0.;
// 	       else
// 	       {
//                  AdvectorValueBeTo = AdvectingField->DOF_value(
// 		     i, j+shift.j, k+shift.k-1,	1, advecting_level );
//                  AdvectorValueFrTo = AdvectingField->DOF_value( 
// 		     i, j+shift.j, k+shift.k, 1, advecting_level );
// 	         vt = 0.5 * ( AdvectorValueBeTo + AdvectorValueFrTo );
// 	         if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP
// 	       		|| DOF_color( i, j+1, k, component ) 
// 				== FV_BC_BEHIND_TOP
// 	       		|| DOF_color( i, j+1, k, component ) 
// 				== FV_BC_FRONT_TOP )
// 	         {
// 	           if ( vt > 0. ) fto = vt * AdvectedValueC;
// 	           else fto = vt * AdvectedValueTo;
// 	         }
// 	         else
// 	         {
// 	           yt = get_DOF_coordinate( j+shift.j, 1, 1 );
// 		   yT = get_DOF_coordinate( j+1, component, 1 );
// 		   dyCt = yt - yC;
// 		   dyt  = yT - yC;
// 		   cLip12 = AdvectedValueC + ( dyCt / dyt ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	 		* ( AdvectedValueTo - AdvectedValueC );
// 	           dyTt = yT - yt;
// 	           dyT = get_cell_size( j+1, component, 1 );
//                    AdvectedValueToTo = DOF_value( 
// 	     	         i, j+2, k, component, advected_level );
// 	   
// 	           thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
// 		   	1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
// 			( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
// 	           cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaTo)
// 	       		* ( AdvectedValueToTo - AdvectedValueTo );
// 		   fto = 0.5 * ( vt * ( cRip12 + cLip12 )
// 	     		- fabs(vt) * ( cRip12 - cLip12 ) );
// 	         }
// 	       }
// 
// 	       // Bottom (Y)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	         fbo = 0.;
// 	       else
// 	       {	     
// 	         AdvectorValueBeBo = AdvectingField->DOF_value( 
// 		     i, j+shift.j-1, k+shift.k-1, 1, advecting_level );
// 	         AdvectorValueFrBo = AdvectingField->DOF_value( 
// 		     i, j+shift.j-1, k+shift.k, 1, advecting_level );
// 	         vb = 0.5 * ( AdvectorValueBeBo + AdvectorValueFrBo );
// 	         if ( DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
// 	       		|| DOF_color( i, j-1, k, component ) 
// 				== FV_BC_BEHIND_BOTTOM
// 	       		|| DOF_color( i, j-1, k, component ) 
// 				== FV_BC_FRONT_BOTTOM )
// 	         {
// 	           if ( vb > 0. ) fbo = vb * AdvectedValueBo;
// 	           else fbo = vb * AdvectedValueC;
// 	         }
// 	         else
// 	         {
// 	           yb = get_DOF_coordinate( j+shift.j-1, 1, 1 );
// 	           yB = get_DOF_coordinate( j-1, component, 1 );
// 	           dyb  = yC - yB;
// 		   if ( DOF_color( i, j, k, component ) == FV_BC_TOP ) 
// 		     cRim12 = AdvectedValueC;
// 		   else
// 		   {
// 		     yT = get_DOF_coordinate( j+1, component, 1 );
// 		     dyt  = yT - yC;
// 		     dyCb = yC - yb;
// 		     cRim12 = AdvectedValueC - ( dyCb / dyt )
// 			* FV_DiscreteField::SuperBee_phi(thetaC)
// 			* ( AdvectedValueTo - AdvectedValueC );
// 	           }
// 		   dyBb = yb - yB;
//                    AdvectedValueBoBo = DOF_value( 
// 	     	         i, j-2, k, component, advected_level );
// 	   
// 	           thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
// 	   		( AdvectedValueBo - AdvectedValueBoBo ) /
// 			( AdvectedValueC - AdvectedValueBo ) : 1.e20;
// 	           cLim12 = AdvectedValueBo + ( dyBb / dyb )
// 		   	* FV_DiscreteField::SuperBee_phi(thetaBo)
// 	       		* ( AdvectedValueC - AdvectedValueBo );
// 	           fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
// 	     		- fabs(vb) * ( cRim12 - cLim12 ) );
// 	         }
// 	       }
// 
// 
// 	       // Front and Behind
// 	       // ----------------
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
// 	       {
// 	         AdvectorValueFr = AdvectorValueC;
// 	         AdvectedValueFr = AdvectedValueC;
// 	       }
// 	       else
// 	       {
// 	         AdvectorValueFr = AdvectingField->DOF_value( i, j, k+1, 
// 		 	component, advecting_level );	     
//                  AdvectedValueFr = DOF_value( 
// 	     	       	i, j, k+1, component, advected_level );
// 	       }
// 
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
// 	       {
// 	         AdvectorValueBe = AdvectorValueC;
// 	         AdvectedValueBe = AdvectedValueC;
// 	       }
// 	       else
// 	       {
//                  AdvectorValueBe = AdvectingField->DOF_value( i, j, k-1, 
// 		 	component, advecting_level );
//                  AdvectedValueBe = DOF_value( 
// 	     	       	i, j, k-1, component, advected_level );
// 	       }
// 
// 	       thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ? 
// 	       		( AdvectedValueC - AdvectedValueBe ) /
// 			( AdvectedValueFr - AdvectedValueC ) : 1.e20;
// 	       
// 	       // Front (Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
// 	         ffr = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {	     
// 	         wf = 0.5 * ( AdvectorValueFr + AdvectorValueC );
// 	         if ( DOF_color( i, j, k+1, component ) == FV_BC_FRONT )
// 		 {
// 	           if ( wf > 0. ) ffr = wf * AdvectedValueC;
// 	           else ffr = wf * AdvectedValueFr;
// 		 }
// 	         else
// 	         {
// 		   zf = get_DOF_coordinate( k+shift.k, 0, 2 );
// 	           zF = get_DOF_coordinate( k+1, component, 2 );
// 	           dzCf = zf - zC;
// 	           dzf  = zF - zC;
// 		   cLip12 = AdvectedValueC + ( dzCf / dzf ) 
// 		   	* FV_DiscreteField::SuperBee_phi(thetaC)
// 	   		* ( AdvectedValueFr - AdvectedValueC );	   
// 
// 	           dzFf = zF - zf;
// 	           dzF = get_cell_size( k+1, component, 2 );
//                    AdvectedValueFrFr = DOF_value( 
// 	     	         i, j, k+2, component, advected_level );
// 
// 	           thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 
// 		   	1.e-20 ? ( AdvectedValueFr - AdvectedValueC ) /
// 			( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
// 	           cRip12 = AdvectedValueFr - ( dzFf / dzF )
// 			* FV_DiscreteField::SuperBee_phi(thetaFr)
// 			* ( AdvectedValueFrFr - AdvectedValueFr );
// 	           ffr = 0.5 * ( wf * ( cRip12 + cLip12 )
// 	       		- fabs(wf) * ( cRip12 - cLip12 ) );
// 		 } 
//                }
// 
// 	       // Behind (Z)
// 	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
// 	         fbe = AdvectorValueC * AdvectedValueC;
// 	       else
// 	       {	     
// 	         wb = 0.5 * (AdvectorValueBe + AdvectorValueC);
// 	         if ( DOF_color(i, j, k-1, component ) == FV_BC_BEHIND )
// 		 {
// 	           if (wb > 0.) fbe = wb * AdvectedValueBe;
// 	           else fbe = wb * AdvectedValueC;
// 		 }
// 	         else
// 	         {
// 	           zb = get_DOF_coordinate( k+shift.k-1, 0, 2 );
// 	           zB = get_DOF_coordinate( k-1, component, 2 );
// 	           dzb  = zC - zB;
// 		   if ( DOF_color( i, j, k, component ) == FV_BC_FRONT ) 
// 		     cRim12 = AdvectedValueC;
// 	           else
// 		   {
// 	             zF = get_DOF_coordinate( k+1, component, 2 );
// 		     dzf  = zF - zC;
// 		     dzCb = zC - zb;
// 	             cRim12 = AdvectedValueC - ( dzCb / dzf )
// 			* FV_DiscreteField::SuperBee_phi(thetaC)
// 			* ( AdvectedValueFr - AdvectedValueC );
// 		   }
// 		   dzBb = zb - zB;
//                    AdvectedValueBeBe = DOF_value( 
// 	     	         i, j, k-2, component, advected_level );
// 	     
// 	           thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
// 		 	( AdvectedValueBe - AdvectedValueBeBe ) /
// 			( AdvectedValueC - AdvectedValueBe ) : 1.e20;
// 	           cLim12 = AdvectedValueBe + ( dzBb / dzb )
// 		   	* FV_DiscreteField::SuperBee_phi(thetaBe)
// 		 	* ( AdvectedValueC - AdvectedValueBe );
// 	           fbe = 0.5 * ( wb * ( cRim12 + cLim12 )
// 	       		- fabs(wb) * ( cRim12 - cLim12 ) );
// 	         }
//                }
// 	     }
// 	     
//              flux = ( fto - fbo ) * dxC * dzC
//                   + ( fri - fle ) * dyC * dzC
//                   + ( ffr - fbe ) * dxC * dyC;
//              VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
// 	   }
// 	 }
//        }
//      }
//    }
// 
//    // Synchronize vector for parallel usage
//    VEC_rhs->synchronize() ;  
//    
// }
//----------------------------------------------------------------------
void 
FV_DiscreteField_Staggered:: assemble_advection_TVD( 
	FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t advected_level,
	LA_Vector *VEC_rhs ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Staggered:: assemble_advection_TVD" );   
   MAC_CHECK_PRE( advected_level < STO_DEPTH ) ;
   MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
   MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ; 
   
   // Parameters
   size_t_vector min_unknown_index(DIM,0);
   size_t_vector max_unknown_index(DIM,0);
   size_t center_pos_in_matrix = 0, component = 0 ;
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
	AdvectedValueFrFr=0.,  
   	AdvectorValueC = 0., AdvectorValueRi = 0., AdvectorValueLe = 0.,
   	AdvectorValueTo = 0., AdvectorValueBo = 0., AdvectorValueFr = 0., 
	AdvectorValueBe = 0, AdvectorValueToLe = 0., AdvectorValueToRi = 0., 
	AdvectorValueBoLe = 0., AdvectorValueBoRi = 0., AdvectorValueFrLe = 0., 
	AdvectorValueFrRi = 0., AdvectorValueBeLe = 0., AdvectorValueBeRi = 0.,
	AdvectorValueFrTo = 0., AdvectorValueFrBo = 0., AdvectorValueBeTo = 0., 
	AdvectorValueBeBo = 0.;
   double ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
	fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
   double cRip12 = 0., cLip12 = 0., cRim12 = 0., cLim12 = 0., thetaC = 0., 
   	thetaRi = 0., thetaLe = 0., thetaTo = 0., thetaBo = 0., thetaFr = 0., 
	thetaBe = 0.;
   FV_SHIFT_TRIPLET shift ;

   // Nullify vector
   VEC_rhs->nullify();

   // Assemble vector
   for (component=0;component<NB_COMPS;++component)
   {
     // Get local min and max indices
     for( size_t l=0; l<DIM; ++l )
     {
       min_unknown_index(l) = 
       	(*min_index_unknown_handled_by_proc)[component](l);
       max_unknown_index(l) =
        (*max_index_unknown_handled_by_proc)[component](l);
     }
     shift = shift_staggeredToStaggered( component ) ;
	     	
     // Perform assembling
     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {          
       xC = get_DOF_coordinate( i, component, 0 );
       dxC = get_cell_size( i, component, 0 ) ;      
       for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
	 yC = get_DOF_coordinate( j, component, 1 );
	 dyC = get_cell_size( j, component, 1 ) ; 
 
         if ( DIM == 2 )
	 {
	   size_t k = 0 ;
	   center_pos_in_matrix = DOF_global_number( i, j, k, component );
           AdvectorValueC = AdvectingField->DOF_value( i, j, k, component, 
	   	advecting_level );
           AdvectedValueC = DOF_value( i, j, k, component, advected_level );
	   
	   // The First component (u)
	   if ( component == 0 )
	   {	       	     
	     // Right and Left
	     // --------------
	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
	     {
	       AdvectorValueRi = AdvectorValueC;
	       AdvectedValueRi = AdvectedValueC;
	     }
	     else
	     {
	       AdvectorValueRi = AdvectingField->DOF_value( i+1, j, k, 
	       		component, advecting_level );	     
	       AdvectedValueRi = DOF_value(
	       		i+1, j, k, component, advected_level );	     
	     }
	     
	     if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
	     {
	       AdvectorValueLe = AdvectorValueC;
	       AdvectedValueLe = AdvectedValueC;
	     }
	     else
	     {
	       AdvectorValueLe = AdvectingField->DOF_value( i-1, j, k, 
	       		component, advecting_level );
	       AdvectedValueLe = DOF_value(
	       		i-1, j, k, component, advected_level );
	     }

	     thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
	   	( AdvectedValueC - AdvectedValueLe ) /
		( AdvectedValueRi - AdvectedValueC ) : 1.e20;
	     
	     // Right (X)
	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
	       fri = AdvectorValueC * AdvectedValueC;
	     else
	     {	     
               ur = 0.5 * ( AdvectorValueRi + AdvectorValueC );
	       if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT )
	       {
                 if ( ur > 0. ) fri = ur * AdvectedValueC;
                 else fri = ur * AdvectedValueRi;
	       }
	       else
	       {
	         xr = get_DOF_coordinate( i+shift.i, 1, 0 );
	         xR = get_DOF_coordinate( i+1, component, 0 );
	         dxCr = xr - xC;
	         dxr  = xR - xC;
		 cLip12 = AdvectedValueC + ( dxCr / dxr )
		 	* FV_DiscreteField::SuperBee_phi(thetaC)
	   		* ( AdvectedValueRi - AdvectedValueC );
			
	         dxRr = xR - xr;
	         dxR = get_cell_size( i+1, component, 0 );
                 AdvectedValueRiRi = DOF_value( 
	     	       i+2, j, k, component, advected_level );
	     
	         thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
		 	1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
			( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
	         cRip12 = AdvectedValueRi - ( dxRr / dxR ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaRi)
			* ( AdvectedValueRiRi - AdvectedValueRi );
	         fri = 0.5 * ( ur * ( cRip12 + cLip12 )
		 		- fabs(ur) * ( cRip12 - cLip12 ) );	 
	       }
             }

	     // Left (X)
	     if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
	       fle = AdvectorValueC * AdvectedValueC;
	     else
	     {	     
	       ul = 0.5 * ( AdvectorValueLe + AdvectorValueC );
	       if ( DOF_color(i-1, j, k, component ) == FV_BC_LEFT )
	       {
	         if ( ul > 0. ) fle = ul * AdvectedValueLe;
	         else fle = ul * AdvectedValueC;
	       }
	       else
	       {
	         xl = get_DOF_coordinate( i+shift.i-1, 1, 0 );
	         xL = get_DOF_coordinate( i-1, component, 0 );
	         dxl  = xC - xL;
	         dxLl = xl - xL;
		 
                 AdvectedValueLeLe = DOF_value( 
	     	       i-2, j, k, component, advected_level );
	         
		 thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
		 	( AdvectedValueLe - AdvectedValueLeLe )/
			( AdvectedValueC - AdvectedValueLe ) : 1.e20;
	         cLim12 = AdvectedValueLe + ( dxLl / dxl ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaLe)
			* ( AdvectedValueC - AdvectedValueLe );
	         if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
		   cRim12 = AdvectedValueC;
	         else
		 {
	           xR = get_DOF_coordinate( i+1, component, 0 );
		   dxr  = xR - xC;
		   dxCl = xC - xl;
		   
	           cRim12 = AdvectedValueC - ( dxCl / dxr ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaC)
		 	* ( AdvectedValueRi - AdvectedValueC );
		 }
	         fle = 0.5 * ( ul * ( cRim12 + cLim12 )
	       		- fabs(ul) * ( cRim12 - cLim12 ) );
	       }
             }
	 
	     // Top and Bottom
	     // --------------
	     if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
	     {
	       AdvectorValueTo = AdvectorValueC;
	       AdvectedValueTo = AdvectedValueC;
	     }
	     else
	     {
	       AdvectorValueTo = AdvectingField->DOF_value( i, j+1, k, 
	       		component, advecting_level );	     
               AdvectedValueTo = DOF_value( 
	     	     i, j+1, k, component, advected_level );
	     }

	     if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
	     {
	       AdvectorValueBo = AdvectorValueC;
	       AdvectedValueBo = AdvectedValueC;
	     }
	     else
             {
	       AdvectorValueBo = AdvectingField->DOF_value( i, j-1, k, 
	       		component, advecting_level );
               AdvectedValueBo = DOF_value( 
			i, j-1, k, component, advected_level );
	     }

	     thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
	 	( AdvectedValueC - AdvectedValueBo ) /
		( AdvectedValueTo - AdvectedValueC ) : 1.e20;

	     // Top (Y)
	     AdvectorValueToLe = AdvectingField->DOF_value( i+shift.i-1, 
	       		j+shift.j, k, 1, advecting_level );
	     AdvectorValueToRi = AdvectingField->DOF_value( i+shift.i, 
	       		j+shift.j, k, 1, advecting_level );
	     vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
	     if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP
	       	|| DOF_color( i, j+1, k, component ) == FV_BC_TOP_LEFT
	       	|| DOF_color( i, j+1, k, component ) == FV_BC_TOP_RIGHT )
	     {
	       if ( vt > 0. ) fto = vt * AdvectedValueC;
	       else fto = vt * AdvectedValueTo;
	     }
	     else
	     {
	       yt = get_DOF_coordinate( j+shift.j, 1, 1 );
	       yT = get_DOF_coordinate( j+1, component, 1 );
	       dyCt = yt - yC;
	       dyt  = yT - yC;
		 
	       cLip12 = AdvectedValueC + ( dyCt / dyt ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaC)
	 		* ( AdvectedValueTo - AdvectedValueC );
	       dyTt = yT - yt;
	       dyT = get_cell_size( j+1, component, 1 );
	         
               AdvectedValueToTo = DOF_value( 
	     	       i, j+2, k, component, advected_level );
	   
	       thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
		 	1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
			( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
	       cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaTo)
	       		* ( AdvectedValueToTo - AdvectedValueTo );
		 
	       fto = 0.5 * ( vt * ( cRip12 + cLip12 )
	     		- fabs(vt) * ( cRip12 - cLip12 ) );
	     }

	     // Bottom (Y)
	     AdvectorValueBoLe = AdvectingField->DOF_value(
	           i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
	     AdvectorValueBoRi = AdvectingField->DOF_value( 
	           i+shift.i, j+shift.j-1, k, 1, advecting_level );
	     vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
	     if ( DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
	       	|| DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM_LEFT
	       	|| DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM_RIGHT )
	     {
	       if ( vb > 0. ) fbo = vb * AdvectedValueBo;
	       else fbo = vb * AdvectedValueC;
	     }
	     else
	     {
	       yb = get_DOF_coordinate( j+shift.j-1, 1, 1 );
	       yB = get_DOF_coordinate( j-1, component, 1 );
	       dyb  = yC - yB;
	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP)
		 cRim12 = AdvectedValueC;
	       else
	       {
		 yT = get_DOF_coordinate( j+1, component, 1 );
		 dyt  = yT - yC;
		 dyCb = yC - yb;
	         cRim12 = AdvectedValueC - ( dyCb / dyt ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaC)
	     		* ( AdvectedValueTo - AdvectedValueC );
	       }
	       dyBb = yb - yB;
               AdvectedValueBoBo = DOF_value( 
	     	       i, j-2, k, component, advected_level );
	   
	       thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
	   		( AdvectedValueBo - AdvectedValueBoBo ) /
			( AdvectedValueC - AdvectedValueBo ) : 1.e20;
	       cLim12 = AdvectedValueBo + ( dyBb / dyb ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaBo)
	       		* ( AdvectedValueC - AdvectedValueBo );
	       fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
	     		- fabs(vb) * ( cRim12 - cLim12 ) );
	     }
	   }
	   
	   // The second component (v)
	   else
	   {
	     // Right and Left
	     // --------------
	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
	     {
	       AdvectorValueRi = AdvectorValueC;
	       AdvectedValueRi = AdvectedValueC;
	     }
	     else
	     {
	       AdvectorValueRi = AdvectingField->DOF_value( i+1, j, k, 
	       		component, advecting_level );	     
               AdvectedValueRi = DOF_value( 
	     	     i+1, j, k, component, advected_level );
	     }

	     if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
	     {
	       AdvectorValueLe = AdvectorValueC;
	       AdvectedValueLe = AdvectedValueC;
	     }
	     else
             {
	       AdvectorValueLe = AdvectingField->DOF_value( i-1, j, k, 
	       		component, advecting_level );
               AdvectedValueLe = DOF_value( 
	     	     i-1, j, k, component, advected_level );
	     }

	     thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
	     	( AdvectedValueC - AdvectedValueLe ) /
		( AdvectedValueRi - AdvectedValueC ) : 1.e20;

	     // Right (X)
	     AdvectorValueToRi = AdvectingField->DOF_value( 
	           i+shift.i, j+shift.j, k, 0, advecting_level );
	     AdvectorValueBoRi = AdvectingField->DOF_value( 
	           i+shift.i, j+shift.j-1, k, 0, advecting_level );
	     ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
	     if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
	       	|| DOF_color( i+1, j, k, component ) == FV_BC_BOTTOM_RIGHT
	       	|| DOF_color( i+1, j, k, component ) == FV_BC_TOP_RIGHT )
	     {
	       if ( ur > 0. ) fri = ur * AdvectedValueC;
	       else fri = ur * AdvectedValueRi;
	     }
	     else
	     {
	       xr = get_DOF_coordinate( i+shift.i, 0, 0 );
	       xR = get_DOF_coordinate( i+1, component, 0 );
	       dxCr = xr - xC;
	       dxr  = xR - xC;
		 
	       cLip12 = AdvectedValueC + ( dxCr / dxr ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaC)
	 	     	* ( AdvectedValueRi - AdvectedValueC );
		 
	       dxRr = xR - xr;
	       dxR = get_cell_size( i+1, component, 0 );
	       AdvectedValueRiRi = DOF_value( 
	     	       i+2, j, k, component, advected_level );
	   
	       thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
		 	1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
			( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
	       cRip12 = AdvectedValueRi - ( dxRr / dxR ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaRi)
	       		* ( AdvectedValueRiRi - AdvectedValueRi );
	       fri = 0.5 * ( ur * ( cRip12 + cLip12 )
	     		- fabs(ur) * ( cRip12 - cLip12 ) );
	     }
	         
	     // Left (X)
	     AdvectorValueToLe = AdvectingField->DOF_value(
	           i+shift.i-1, j+shift.j, k, 0, advecting_level );
	     AdvectorValueBoLe = AdvectingField->DOF_value(
	           i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
	     ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
	     if ( DOF_color( i-1, j, k, component ) == FV_BC_LEFT
	       	|| DOF_color( i-1, j, k, component ) == FV_BC_BOTTOM_LEFT
	       	|| DOF_color( i-1, j, k, component ) == FV_BC_TOP_LEFT)
	     {
	       if ( ul > 0. ) fle = ul * AdvectedValueLe;
	       else fle = ul * AdvectedValueC;
	     }
	     else
	     {
	       xl = get_DOF_coordinate( i+shift.i-1, 0, 0 );
	       xL = get_DOF_coordinate( i-1, component, 0 );
	       dxl  = xC - xL;
	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
		 cRim12 = AdvectedValueC;
	       else
	       {
		 xR = get_DOF_coordinate( i+1, component, 0 );
		 dxr  = xR - xC;
		 dxCl = xC - xl;
		 cRim12 = AdvectedValueC
			- ( dxCl / dxr ) 
			* FV_DiscreteField::SuperBee_phi(thetaC)
			* ( AdvectedValueRi - AdvectedValueC );
	       }
	       dxLl = xl - xL;
               AdvectedValueLeLe = DOF_value( 
	     	       i-2, j, k, component, advected_level );
	   
	       thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
	   		( AdvectedValueLe - AdvectedValueLeLe ) /
			( AdvectedValueC - AdvectedValueLe ) : 1.e20;
	       cLim12 = AdvectedValueLe + ( dxLl / dxl ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaLe)
	       		* ( AdvectedValueC - AdvectedValueLe );
	       fle = 0.5 * ( ul * ( cRim12 + cLim12 )
	     		- fabs(ul) * ( cRim12 - cLim12 ) );
	     }
	 
	     // Top and Bottom
	     // --------------
	     if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
	     {
	       AdvectorValueTo = AdvectorValueC;
	       AdvectedValueTo = AdvectedValueC;
	     }
	     else
	     {
	       AdvectorValueTo = AdvectingField->DOF_value( i, j+1, k, 
	       		component, advecting_level );	     
               AdvectedValueTo = DOF_value( 
	     	     i, j+1, k, component, advected_level );
	     }

	     if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
	     {
	       AdvectorValueBo = AdvectorValueC;
	       AdvectedValueBo = AdvectedValueC;
	     }
	     else
	     {
               AdvectorValueBo = AdvectingField->DOF_value( i, j-1, k, 
	       		component, advecting_level );
               AdvectedValueBo = DOF_value( 
	     	       i, j-1, k, component, advected_level );
	     }

	     thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
	   	( AdvectedValueC - AdvectedValueBo ) /
		( AdvectedValueTo - AdvectedValueC ) : 1.e20;
	     
	     // Top (Y)
	     if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
	       fto = AdvectorValueC * AdvectedValueC;
	     else
	     {	     
               vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
	       if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP )
	       {
	         if ( vt > 0. ) fto = vt * AdvectedValueC;
	         else fto = vt * AdvectedValueTo;
	       }
	       else
	       {
	         yt = get_DOF_coordinate( j+shift.j, 0, 1 );
	         yT = get_DOF_coordinate( j+1, component, 1 );
	         dyCt = yt - yC;
	         dyt  = yT - yC;
		 cLip12 = AdvectedValueC + ( dyCt / dyt ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaC)
	   		* ( AdvectedValueTo - AdvectedValueC );

	         dyTt = yT - yt;
	         dyT = get_cell_size( j+1, component, 1 );
                 AdvectedValueToTo = DOF_value( 
	     	       i, j+2, k, component, advected_level );
	     
	         thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
		 	1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
		     	( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
	         cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaTo)
			* ( AdvectedValueToTo - AdvectedValueTo );
	         fto = 0.5 * ( vt * ( cRip12 + cLip12 )
	       		- fabs(vt) * ( cRip12 - cLip12 ) );
	       }	   
             }

	     // Bottom (Y)
	     if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
	       fbo = AdvectorValueC * AdvectedValueC;
	     else
	     {
	       vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
	       if ( DOF_color(i,j-1,k,component) == FV_BC_BOTTOM )
	       {
	         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
	         else fbo = vb * AdvectedValueC;
	       }
	       else
	       {
	         yb = get_DOF_coordinate( j+shift.j-1, 0, 1 );
	         yB = get_DOF_coordinate( j-1, component, 1 );
	         dyb  = yC - yB;
	       
	         dyBb = yb - yB;
                 AdvectedValueBoBo = DOF_value( 
	     	       i, j-2, k, component, advected_level );
	     
	         thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
		 	( AdvectedValueBo - AdvectedValueBoBo ) /
			( AdvectedValueC - AdvectedValueBo ) : 1.e20;
	         cLim12 = AdvectedValueBo + ( dyBb / dyb ) 
		 	* FV_DiscreteField::SuperBee_phi(thetaBo)
			* ( AdvectedValueC - AdvectedValueBo );
	   
	         if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
		   cRim12 = AdvectedValueC;
	         else
		 {
		   yT = get_DOF_coordinate( j+1, component, 1 );
		   dyt  = yT - yC;
		   dyCb = yC - yb;
	           cRim12 = AdvectedValueC - ( dyCb / dyt )
		   	* FV_DiscreteField::SuperBee_phi(thetaC)
		 	* ( AdvectedValueTo - AdvectedValueC );
	         }
		 fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
	       		- fabs(vb) * ( cRim12 - cLim12 ) );
	       }
             }
	   }

           flux = ( fto - fbo ) * dxC + ( fri - fle ) * dyC;
           VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
	 }
	 else // DIM = 3
	 {
	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	   {
	     zC = get_DOF_coordinate( k, component, 2 ) ;
	     dzC = get_cell_size( k, component, 2 ) ;	     
	     center_pos_in_matrix = DOF_global_number( i, j, k, component );
             AdvectorValueC = AdvectingField->DOF_value( i, j, k, component, 
	     	advecting_level );
             AdvectedValueC = DOF_value( 
	           i, j, k, component, advected_level );

	     // The First component (u)
	     if ( component == 0 )
	     {
	       // Right and Left
	       // --------------
	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
	       {
	         AdvectorValueRi = AdvectorValueC;
	         AdvectedValueRi = AdvectedValueC;
	       }
	       else
	       {
	         AdvectorValueRi = AdvectingField->DOF_value( i+1, j, k, 
		 	component, advecting_level );
                 AdvectedValueRi = DOF_value( 
	     	       	i+1, j, k, component, advected_level );
	       }

	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
	       {
	         AdvectorValueLe = AdvectorValueC;
	         AdvectedValueLe = AdvectedValueC;
	       }
	       else
	       {
	         AdvectorValueLe = AdvectingField->DOF_value( i-1, j, k, 
		 	component, advecting_level );
                 AdvectedValueLe = DOF_value( 
	     	       i-1, j, k, component, advected_level );
	      }

	       thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
	       		( AdvectedValueC - AdvectedValueLe ) /
			( AdvectedValueRi - AdvectedValueC ) : 1.e20;
	     
	       // Right (X)
	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
	         fri = AdvectorValueC * AdvectedValueC;
	       else
	       {	     
	         ur = 0.5 * ( AdvectorValueRi + AdvectorValueC );
	         if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT )
		 {
                   if ( ur > 0. ) fri = ur * AdvectedValueC;
                   else fri = ur * AdvectedValueRi;
		 }
	         else
	         {		 
		   xr = get_DOF_coordinate( i+shift.i, 1, 0 );
	           xR = get_DOF_coordinate( i+1, component, 0 );
	           dxCr = xr - xC;
	           dxr  = xR - xC;
		   cLip12 = AdvectedValueC + ( dxCr / dxr ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaC)
	   		* ( AdvectedValueRi - AdvectedValueC );
			
	           dxRr = xR - xr;
	           dxR = get_cell_size( i+1, component, 0 );
                   AdvectedValueRiRi = DOF_value( 
	     	         i+2, j, k, component, advected_level );
	     
	           thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
		   	1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
			( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20 ;
	           cRip12 = AdvectedValueRi - ( dxRr / dxR ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaRi)
			* ( AdvectedValueRiRi - AdvectedValueRi );
	           fri = 0.5 * ( ur * ( cRip12 + cLip12 )
	       		- fabs(ur) * ( cRip12 - cLip12 ) );
		}
               }

	       // Left (X)
	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
	         fle = AdvectorValueC * AdvectedValueC;
	       else
	       {	     
	         ul = 0.5 * ( AdvectorValueLe + AdvectorValueC );
	         if ( DOF_color( i-1, j, k, component ) == FV_BC_LEFT )
		 {
	           if ( ul > 0. ) fle = ul * AdvectedValueLe;
	           else fle = ul * AdvectedValueC;
		 }
	         else
	         {
		   xl = get_DOF_coordinate( i+shift.i-1, 1, 0 );
		   xL = get_DOF_coordinate( i-1, component, 0 );
		   dxl  = xC - xL;		   
		   dxLl = xl - xL;
                   AdvectedValueLeLe = DOF_value( 
	     	       i-2, j, k, component, advected_level );

	           thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
		 	( AdvectedValueLe - AdvectedValueLeLe ) /
			( AdvectedValueC - AdvectedValueLe ) : 1.e20;
	           cLim12 = AdvectedValueLe + ( dxLl / dxl ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaLe)
			* ( AdvectedValueC - AdvectedValueLe );
	   
	           if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
		     cRim12 = AdvectedValueC;
	           else
		   {
	             xR = get_DOF_coordinate( i+1, component, 0 );
		     dxr  = xR - xC;
		     dxCl = xC - xl;
		     cRim12 = AdvectedValueC - ( dxCl / dxr )
			* FV_DiscreteField::SuperBee_phi(thetaC)
			* ( AdvectedValueRi - AdvectedValueC );
	           }
	           fle = 0.5 * ( ul * ( cRim12 + cLim12 )
	       		- fabs(ul) * ( cRim12 - cLim12 ) );
	         }
               }

	 
	       // Top and Bottom
	       // --------------
	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
	       {
	         AdvectorValueTo = AdvectorValueC;
	         AdvectedValueTo = AdvectedValueC;
	       }
	       else
	       {
	         AdvectorValueTo = AdvectingField->DOF_value( i, j+1, k, 
		 	component, advecting_level );	     
                 AdvectedValueTo = DOF_value( 
	     	       	i, j+1, k, component, advected_level );
               }
	       
	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
	       {
	         AdvectorValueBo = AdvectorValueC;
	         AdvectedValueBo = AdvectedValueC;
	       }
	       else
               {
	         AdvectorValueBo = AdvectingField->DOF_value( i, j-1, k, 
		 	component, advecting_level );
                 AdvectedValueBo = DOF_value( 
	     	       	i, j-1, k, component, advected_level );
               }
	       
	       thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
			( AdvectedValueC - AdvectedValueBo ) /
			( AdvectedValueTo - AdvectedValueC ) : 1.e20;

	       // Top (Y)
	       AdvectorValueToLe = AdvectingField->DOF_value(
		     i+shift.i-1, j+shift.j, k, 1, advecting_level );
	       AdvectorValueToRi = AdvectingField->DOF_value( 
		     i+shift.i, j+shift.j, k, 1, advecting_level );
	       vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
	       if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP
	       	  	|| DOF_color( i, j+1, k, component ) == FV_BC_TOP_LEFT
	       	  	|| DOF_color( i, j+1, k, component ) 
				== FV_BC_TOP_RIGHT )
	       {
	         if ( vt > 0. ) fto = vt * AdvectedValueC;
	         else fto = vt * AdvectedValueTo;
	       }
	       else
	       {
	         yt = get_DOF_coordinate( j+shift.j, 1, 1 );
		 yT = get_DOF_coordinate( j+1, component, 1 );
		 dyCt = yt - yC;
		 dyt  = yT - yC;
		 cLip12 = AdvectedValueC + ( dyCt / dyt ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaC)
	 		* ( AdvectedValueTo - AdvectedValueC );
	         dyTt = yT - yt;
	         dyT = get_cell_size( j+1, component, 1 );
                 AdvectedValueToTo = DOF_value( 
	     	         i, j+2, k, component, advected_level );
	   
	         thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
		   	1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
			( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
	         cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaTo)
	       		* ( AdvectedValueToTo - AdvectedValueTo );
		 fto = 0.5 * ( vt * ( cRip12 + cLip12 )
	     		- fabs(vt) * ( cRip12 - cLip12 ) );
	       }

	       // Bottom (Y)
	       AdvectorValueBoLe = AdvectingField->DOF_value( 
		     i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
	       AdvectorValueBoRi = AdvectingField->DOF_value( 
		     i+shift.i, j+shift.j-1, k,	1, advecting_level );
	       vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
	       if ( DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
	       	  	|| DOF_color( i, j-1, k, component ) 
				== FV_BC_BOTTOM_LEFT
	          	|| DOF_color( i, j-1, k, component ) 
				== FV_BC_BOTTOM_RIGHT )
	       {
	         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
	         else fbo = vb * AdvectedValueC;
	       }
	       else
	       {
	         yb = get_DOF_coordinate( j+shift.j-1, 1, 1 );
	         yB = get_DOF_coordinate( j-1, component, 1 );
	         dyb  = yC - yB;
		 if ( DOF_color( i, j, k, component ) == FV_BC_TOP ) 
		   cRim12 = AdvectedValueC;
		 else
		 {
		   yT = get_DOF_coordinate( j+1, component, 1 );
		   dyt  = yT - yC;
		   dyCb = yC - yb;
		   cRim12 = AdvectedValueC - ( dyCb / dyt )
		     	* FV_DiscreteField::SuperBee_phi(thetaC)
	     		* ( AdvectedValueTo - AdvectedValueC );
	         }
		 dyBb = yb - yB;
                 AdvectedValueBoBo = DOF_value( 
	     	         i, j-2, k, component, advected_level );
	   
	         thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
	   		( AdvectedValueBo - AdvectedValueBoBo ) /
			( AdvectedValueC - AdvectedValueBo ) : 1.e20;
	         cLim12 = AdvectedValueBo + ( dyBb / dyb )
		   	* FV_DiscreteField::SuperBee_phi(thetaBo)
	       		* ( AdvectedValueC - AdvectedValueBo );
	         fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
	     		- fabs(vb) * ( cRim12 - cLim12 ) );
	       }

	       // Front and Behind
	       // ----------------
	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
	       {
	         AdvectorValueFr = AdvectorValueC;
	         AdvectedValueFr = AdvectedValueC;
	       }
	       else
	       {
	         AdvectorValueFr = AdvectingField->DOF_value( i, j, k+1, 
		 	component, advecting_level );	     
                 AdvectedValueFr = DOF_value( 
	     	       	i, j, k+1, component, advected_level );
               }
	       
	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
	       {
	         AdvectorValueBe = AdvectorValueC;
	         AdvectedValueBe = AdvectedValueC;
	       }
	       else
               {
	         AdvectorValueBe = AdvectingField->DOF_value( i, j, k-1, 
		 	component, advecting_level );
                 AdvectedValueBe = DOF_value( 
	     	       	i, j, k-1, component, advected_level );
               }
	       
	       thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ? 
	       		( AdvectedValueC - AdvectedValueBe ) /
			( AdvectedValueFr - AdvectedValueC ) : 1.e20;
	       
	       // Front (Z)
	       AdvectorValueFrLe = AdvectingField->DOF_value( 
		     i+shift.i-1, j, k+shift.k, 2, advecting_level );
	       AdvectorValueFrRi = AdvectingField->DOF_value( 
		     i+shift.i, j, k+shift.k, 2, advecting_level );
	       wf = 0.5 * (AdvectorValueFrLe + AdvectorValueFrRi);
	       if ( DOF_color( i, j, k+1, component ) == FV_BC_FRONT
	       		|| DOF_color( i, j, k+1, component ) 
				== FV_BC_FRONT_LEFT 
	       		|| DOF_color( i, j, k+1, component ) 
				== FV_BC_FRONT_RIGHT )
	       {
		 if ( wf > 0. ) ffr = wf * AdvectedValueC;
	         else ffr = wf * AdvectedValueFr;
	       }
	       else
	       {
		 zf = get_DOF_coordinate( k+shift.k, 2, 2 );
		 zF = get_DOF_coordinate( k+1, component, 2 );
		 dzCf = zf - zC;
		 dzf  = zF - zC;
		 cLip12 = AdvectedValueC + ( dzCf / dzf )
		   	* FV_DiscreteField::SuperBee_phi(thetaC)
	 		* ( AdvectedValueFr - AdvectedValueC );
		 dzFf = zF - zf;
		 dzF = get_cell_size( k+1, component, 2 );
                 AdvectedValueFrFr = DOF_value( 
	     	         i, j, k+2, component, advected_level );

		 thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 
		   	1.e-20 ? ( AdvectedValueFr - AdvectedValueC ) /
			( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
	         cRip12 = AdvectedValueFr - ( dzFf / dzF )
		 	* FV_DiscreteField::SuperBee_phi(thetaFr)
	       		* ( AdvectedValueFrFr - AdvectedValueFr );
		 ffr = 0.5 * ( wf * ( cRip12 + cLip12 )
	     		- fabs(wf) * ( cRip12 - cLip12 ) );
	       }

	       // Behind (Z)
	       AdvectorValueBeLe = AdvectingField->DOF_value( 
		     i+shift.i-1, j, k+shift.k-1, 2, advecting_level );
	       AdvectorValueBeRi = AdvectingField->DOF_value(
		     i+shift.i, j, k+shift.k-1, 2, advecting_level );
	       wb = 0.5 * ( AdvectorValueBeLe + AdvectorValueBeRi );
	       if ( DOF_color( i, j, k-1, component ) == FV_BC_BEHIND
		 	|| DOF_color( i, j, k-1, component ) 
				== FV_BC_BEHIND_LEFT
		 	|| DOF_color( i, j, k-1, component ) 
				== FV_BC_BEHIND_RIGHT )
	       {
	         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
	         else fbe = wb * AdvectedValueC;
	       }
	       else
	       {
		 zb = get_DOF_coordinate( k+shift.k-1, 2, 2 );
		 zB = get_DOF_coordinate( k-1, component, 2 );
		 dzb  = zC - zB;
		 if ( DOF_color( i, j, k, component ) == FV_BC_FRONT ) 
		   cRim12 = AdvectedValueC;
		 else
		 {
		   zF = get_DOF_coordinate( k+1, component, 2 );
		   dzf  = zF - zC;
		   dzCb = zC - zb;
		   cRim12 = AdvectedValueC - ( dzCb / dzf )
			* FV_DiscreteField::SuperBee_phi(thetaC)
			* ( AdvectedValueFr - AdvectedValueC );
		 }
		 dzBb = zb - zB;
                 AdvectedValueBeBe = DOF_value( 
			i, j, k-2, component, advected_level );
	   
	         thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
	   		( AdvectedValueBe - AdvectedValueBeBe ) /
			( AdvectedValueC - AdvectedValueBe ) : 1.e20;
	         cLim12 = AdvectedValueBe + ( dzBb / dzb )
		 	* FV_DiscreteField::SuperBee_phi(thetaBe)
	       		* ( AdvectedValueC - AdvectedValueBe );
		 fbe = 0.5 * ( wb * ( cRim12 + cLim12 )
	     		- fabs(wb) * ( cRim12 - cLim12 ) );
	       }
	     }
	     
	     // The Second component (v)
	     else if ( component == 1 )
	     {
	       // Right and Left
	       // --------------
	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
	       {
	         AdvectorValueRi = AdvectorValueC;
	         AdvectedValueRi = AdvectedValueC;
	       }
	       else
	       {
	         AdvectorValueRi = AdvectingField->DOF_value( i+1, j, k, 
		 	component, advecting_level );	     
                 AdvectedValueRi = DOF_value( 
	     	      	i+1, j, k, component, advected_level );
               }
	       
	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
	       {
	         AdvectorValueLe = AdvectorValueC;
	         AdvectedValueLe = AdvectedValueC;
	       }
	       else
	       {
                 AdvectorValueLe = AdvectingField->DOF_value( i-1, j, k, 
		 	component, advecting_level );
                 AdvectedValueLe = DOF_value( 
	     	       	i-1, j, k, component, advected_level );
	       }

	       thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
	       		( AdvectedValueC - AdvectedValueLe ) /
			( AdvectedValueRi - AdvectedValueC ) : 1.e20;

	       // Right (X)
	       AdvectorValueToRi = AdvectingField->DOF_value( 
		     i+shift.i, j+shift.j, k, 0, advecting_level );
	       AdvectorValueBoRi = AdvectingField->DOF_value(
		     i+shift.i, j+shift.j-1, k,	0, advecting_level );
	       ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
	       if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
	       		|| DOF_color( i+1, j, k, component ) 
				== FV_BC_BOTTOM_RIGHT
	       		|| DOF_color( i+1, j, k, component ) 
				== FV_BC_TOP_RIGHT )
	       {
		 if ( ur > 0. ) fri = ur * AdvectedValueC;
		 else fri = ur * AdvectedValueRi;
	       }
	       else
	       {
		 xr = get_DOF_coordinate( i+shift.i, 0, 0 );
		 xR = get_DOF_coordinate( i+1, component, 0 );
		 dxCr = xr - xC;
		 dxr  = xR - xC;
		 cLip12 = AdvectedValueC + ( dxCr / dxr ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaC)
	 		* ( AdvectedValueRi - AdvectedValueC );
		   
		 dxRr = xR - xr;
	         dxR = get_cell_size( i+1, component, 0 );
                 AdvectedValueRiRi = DOF_value( 
			i+2, j, k, component, advected_level );
	   
	         thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
		   	1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
			( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
	         cRip12 = AdvectedValueRi - ( dxRr / dxR )
			* FV_DiscreteField::SuperBee_phi(thetaRi)
			* ( AdvectedValueRiRi - AdvectedValueRi );
		 fri = 0.5 * ( ur * ( cRip12 + cLip12 )
	     		- fabs(ur) * ( cRip12 - cLip12 ) );
	       }
	       
	       // Left (X)
	       AdvectorValueToLe = AdvectingField->DOF_value(
		     i+shift.i-1, j+shift.j, k, 0, advecting_level );
	       AdvectorValueBoLe = AdvectingField->DOF_value(
		     i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
	       ul = 0.5 * (AdvectorValueToLe + AdvectorValueBoLe);
	       if ( DOF_color( i-1, j, k, component ) == FV_BC_LEFT
	       		|| DOF_color( i-1, j, k, component ) 
				== FV_BC_BOTTOM_LEFT
	       		|| DOF_color( i-1, j, k, component ) 
				== FV_BC_TOP_LEFT )
	       {
	         if ( ul > 0. ) fle = ul * AdvectedValueLe;
	         else fle = ul * AdvectedValueC;
	       }
	       else
	       {
		 xl = get_DOF_coordinate( i+shift.i-1, 0, 0 );
		 xL = get_DOF_coordinate( i-1, component, 0 );
		 dxl  = xC - xL;
		 if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
		   cRim12 = AdvectedValueC;
		 else
		 {
		   xR = get_DOF_coordinate( i+1, component, 0 );
		   dxr  = xR - xC;
		   dxCl = xC - xl;
		   cRim12 = AdvectedValueC - ( dxCl / dxr )
			* FV_DiscreteField::SuperBee_phi(thetaC)
	     		* ( AdvectedValueRi - AdvectedValueC );
		 }
		   
		 dxLl = xl - xL;
                 AdvectedValueLeLe = DOF_value( 
			i-2, j, k, component, advected_level );

	         thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
	   		( AdvectedValueLe - AdvectedValueLeLe ) /
			( AdvectedValueC - AdvectedValueLe ) : 1.e20;
	         cLim12 = AdvectedValueLe + ( dxLl / dxl )
			* FV_DiscreteField::SuperBee_phi(thetaLe)
			* ( AdvectedValueC - AdvectedValueLe );
		 fle = 0.5 * ( ul * ( cRim12 + cLim12 )
	     		- fabs(ul) * ( cRim12 - cLim12 ) );
	       }


	       // Top and Bottom
	       // --------------
	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
	       {
	         AdvectorValueTo = AdvectorValueC;
	         AdvectedValueTo = AdvectedValueC;
	       }
	       else
	       {
	         AdvectorValueTo = AdvectingField->DOF_value( i, j+1, k, 
		 	component, advecting_level );	     
                 AdvectedValueTo = DOF_value( 
	     	       	i, j+1, k, component, advected_level );
	       }

	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
	       {
	         AdvectorValueBo = AdvectorValueC;
	         AdvectedValueBo = AdvectedValueC;
	       }
	       else
	       {
                 AdvectorValueBo = AdvectingField->DOF_value( i, j-1, k, 
		 	component, advecting_level );
                 AdvectedValueBo = DOF_value( 
	     	       	i, j-1, k, component, advected_level );
               }
	       thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
	       		( AdvectedValueC - AdvectedValueBo ) /
			( AdvectedValueTo - AdvectedValueC ) : 1.e20;
	     
	       // Top (Y)
	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
	         fto = AdvectorValueC * AdvectedValueC;
	       else
	       {	     
	         vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
	         if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP )
		 {
	           if ( vt > 0. ) fto = vt * AdvectedValueC;
	           else fto = vt * AdvectedValueTo;
		 }
	         else
	         {
		   yt = get_DOF_coordinate( j+shift.j, 0, 1 );
	           yT = get_DOF_coordinate( j+1, component, 1 );
	           dyCt = yt - yC;
	           dyt  = yT - yC;
		   cLip12 = AdvectedValueC + ( dyCt / dyt )
		   	* FV_DiscreteField::SuperBee_phi(thetaC)
	   		* ( AdvectedValueTo - AdvectedValueC );
			
	           dyTt = yT - yt;
	           dyT = get_cell_size( j+1, component, 1 );
                   AdvectedValueToTo = DOF_value( 
	     	         i, j+2, k, component, advected_level );
	     
	           thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
		   	1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
		       ( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
	           cRip12 = AdvectedValueTo - ( dyTt / dyT )
			* FV_DiscreteField::SuperBee_phi(thetaTo)
			* ( AdvectedValueToTo - AdvectedValueTo );
	           fto = 0.5 * ( vt * ( cRip12 + cLip12 )
	       		- fabs(vt) * ( cRip12 - cLip12 ) );
		 }
               }
	 
	       // Bottom (Y)
	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
	         fbo = AdvectorValueC * AdvectedValueC;
	       else
	       {
	         vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
		 if ( DOF_color(i, j-1, k, component ) == FV_BC_BOTTOM )
		 {
		   if ( vb > 0. ) fbo = vb * AdvectedValueBo;
	           else fbo = vb * AdvectedValueC;
		 }
	         else
	         {
	           yb = get_DOF_coordinate( j+shift.j-1, 0, 1 );
	           yB = get_DOF_coordinate( j-1, component, 1 );
	           dyb  = yC - yB;
		   if ( DOF_color( i, j, k, component ) == FV_BC_TOP ) 
		     cRim12 = AdvectedValueC;
	           else
		   {
	             yT = get_DOF_coordinate( j+1, component, 1 );
		     dyt  = yT - yC;
		     dyCb = yC - yb;
	             cRim12 = AdvectedValueC - ( dyCb / dyt )
			* FV_DiscreteField::SuperBee_phi(thetaC)
			* ( AdvectedValueTo - AdvectedValueC );
	           }
		   dyBb = yb - yB;
                   AdvectedValueBoBo = DOF_value( 
	     	         i, j-2, k, component, advected_level );
	     
	           thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
		 	( AdvectedValueBo - AdvectedValueBoBo ) /
			( AdvectedValueC - AdvectedValueBo ) : 1.e20;
	           cLim12 = AdvectedValueBo + ( dyBb / dyb )
			* FV_DiscreteField::SuperBee_phi(thetaBo)
			* ( AdvectedValueC - AdvectedValueBo );
	           fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
	       		- fabs(vb) * ( cRim12 - cLim12 ) );
	         }
               }


	       // Front and Behind
	       // ----------------
	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
	       {
	         AdvectorValueFr = AdvectorValueC;
	         AdvectedValueFr = AdvectedValueC;
	       }
	       else
	       {
	         AdvectorValueFr = AdvectingField->DOF_value( i, j, k+1, 
		 	component, advecting_level );	     
                 AdvectedValueFr = DOF_value( 
	     	       	i, j, k+1, component, advected_level );
	       }

	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
	       {
	         AdvectorValueBe = AdvectorValueC;
	         AdvectedValueBe = AdvectedValueC;
	       }
	       else
	       {
                 AdvectorValueBe = AdvectingField->DOF_value( i, j, k-1, 
		 	component, advecting_level );
                 AdvectedValueBe = DOF_value( 
	     	       	i, j, k-1, component, advected_level );
	       }

	       thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ? 
	       		( AdvectedValueC - AdvectedValueBe ) /
			( AdvectedValueFr - AdvectedValueC ) : 1.e20;
	       
	       // Front (Z)
	       AdvectorValueFrBo = AdvectingField->DOF_value( 
		     i, j+shift.j-1, k+shift.k, 2, advecting_level );
	       AdvectorValueFrTo = AdvectingField->DOF_value( 
		     i, j+shift.j, k+shift.k, 2, advecting_level );
               wf = 0.5 * ( AdvectorValueFrBo + AdvectorValueFrTo );
	       if ( DOF_color( i, j, k+1, component ) == FV_BC_FRONT
	       		|| DOF_color( i, j, k+1, component ) 
				== FV_BC_FRONT_BOTTOM
	       		|| DOF_color( i, j, k+1, component ) 
				== FV_BC_FRONT_TOP )
	       {
		 if ( wf > 0. ) ffr = wf * AdvectedValueC;
	         else ffr = wf * AdvectedValueFr;
	       }
	       else
	       {
		 zf = get_DOF_coordinate( k+shift.k, 2, 2 );
		 zF = get_DOF_coordinate( k+1, component, 2 );
		 dzCf = zf - zC;
		 dzf  = zF - zC;
		 cLip12 = AdvectedValueC + ( dzCf / dzf )
			* FV_DiscreteField::SuperBee_phi(thetaC)
			* ( AdvectedValueFr - AdvectedValueC );
		 dzFf = zF - zf;
		 dzF = get_cell_size( k+1, component, 2 );
                 AdvectedValueFrFr = DOF_value( 
			i, j, k+2, component, advected_level );

		 thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 
		   	1.e-20 ? ( AdvectedValueFr - AdvectedValueC ) /
			( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
	         cRip12 = AdvectedValueFr - ( dzFf / dzF )
			* FV_DiscreteField::SuperBee_phi(thetaFr)
			* ( AdvectedValueFrFr - AdvectedValueFr );
		 ffr = 0.5 * ( wf * ( cRip12 + cLip12 )
	     		- fabs(wf) * ( cRip12 - cLip12 ) );
	       }

	       // Behind (Z)
	       AdvectorValueBeBo = AdvectingField->DOF_value( 
		     i, j+shift.j-1, k+shift.k-1, 2, advecting_level );
	       AdvectorValueBeTo = AdvectingField->DOF_value(
		     i, j+shift.j, k+shift.k-1, 2, advecting_level );
	       wb = 0.5 * ( AdvectorValueBeBo + AdvectorValueBeTo );
	       if ( DOF_color( i, j, k-1, component ) == FV_BC_BEHIND
	       		|| DOF_color( i, j, k-1, component ) 
				== FV_BC_BEHIND_BOTTOM
	       		|| DOF_color( i, j, k-1, component ) 
				== FV_BC_BEHIND_TOP )
	       {
	         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
	         else fbe = wb * AdvectedValueC;
	       }
	       else
	       {
		 zb = get_DOF_coordinate( k+shift.k-1, 2, 2 );
		 zB = get_DOF_coordinate( k-1, component, 2 );
		 dzb  = zC - zB;
		 if ( DOF_color( i, j, k, component ) == FV_BC_FRONT ) 
		   cRim12 = AdvectedValueC;
		 else
		 {
		   zF = get_DOF_coordinate( k+1, component, 2 );
		   dzf  = zF - zC;
		   dzCb = zC - zb;
		   cRim12 = AdvectedValueC - ( dzCb / dzf )
			* FV_DiscreteField::SuperBee_phi(thetaC)
			* ( AdvectedValueFr - AdvectedValueC );
		 }
		 dzBb = zb - zB;
                 AdvectedValueBeBe = DOF_value( 
	     	         i, j, k-2, component, advected_level );
	   
	         thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
	   		( AdvectedValueBe - AdvectedValueBeBe ) /
			( AdvectedValueC - AdvectedValueBe ) : 1.e20;
	         cLim12 = AdvectedValueBe + ( dzBb / dzb )
			* FV_DiscreteField::SuperBee_phi(thetaBe)
			* ( AdvectedValueC - AdvectedValueBe );
		 fbe = 0.5 * ( wb * ( cRim12 + cLim12 )
	     		- fabs(wb) * ( cRim12 - cLim12 ) );
	       }
	     }
	     
	     // The Third component (w)
	     else
	     {
	       // Right and Left
	       // --------------
	       if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
	       {
	         AdvectorValueRi = AdvectorValueC;
	         AdvectedValueRi = AdvectedValueC;
	       }
	       else
	       {
	         AdvectorValueRi = AdvectingField->DOF_value( i+1, j, k, 
		 	component, advecting_level );
                 AdvectedValueRi = DOF_value( 
	     	       	i+1, j, k, component, advected_level );
	       }

	       if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
	       {
	         AdvectorValueLe = AdvectorValueC;
	         AdvectedValueLe = AdvectedValueC;
	       }
	       else
	       {
                 AdvectorValueLe = AdvectingField->DOF_value( i-1, j, k, 
		 	component, advecting_level );
                 AdvectedValueLe = DOF_value( 
	     	       	i-1, j, k, component, advected_level );
	       }

	       thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
	       		( AdvectedValueC - AdvectedValueLe ) /
			( AdvectedValueRi - AdvectedValueC ) : 1.e20;

	       // Right (X)
	       AdvectorValueFrRi = AdvectingField->DOF_value(
		     i+shift.i, j, k+shift.k, 0, advecting_level );
	       AdvectorValueBeRi = AdvectingField->DOF_value( 
		     i+shift.i, j, k+shift.k-1,	0, advecting_level );
	       ur = 0.5 * (AdvectorValueFrRi + AdvectorValueBeRi);
	       if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
	       		|| DOF_color( i+1, j, k, component ) 
				== FV_BC_BEHIND_RIGHT
	       		|| DOF_color( i+1, j, k, component ) 
				== FV_BC_FRONT_RIGHT )
	       {
		 if ( ur > 0. ) fri = ur * AdvectedValueC;
		 else fri = ur * AdvectedValueRi;
	       }
	       else
	       {
		 xr = get_DOF_coordinate( i+shift.i, 0, 0 );
		 xR = get_DOF_coordinate( i+1, component, 0 );
		 dxCr = xr - xC;
		 dxr  = xR - xC;
		 cLip12 = AdvectedValueC + ( dxCr / dxr ) 
			* FV_DiscreteField::SuperBee_phi(thetaC)
			* ( AdvectedValueRi - AdvectedValueC );
		   
		 dxRr = xR - xr;
	         dxR = get_cell_size( i+1, component, 0 );
                 AdvectedValueRiRi = DOF_value( 
			i+2, j, k, component, advected_level );
	   
	         thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
		   	1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
			( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
	         cRip12 = AdvectedValueRi - ( dxRr / dxR )
			* FV_DiscreteField::SuperBee_phi(thetaRi)
			* ( AdvectedValueRiRi - AdvectedValueRi );
		 fri = 0.5 * ( ur * ( cRip12 + cLip12 )
	     		- fabs(ur) * ( cRip12 - cLip12 ) );
	       }
	       
	       // Left (X)
	       AdvectorValueFrLe = AdvectingField->DOF_value(
		     i+shift.i-1, j, k+shift.k, 0, advecting_level );
	       AdvectorValueBeLe = AdvectingField->DOF_value(
		     i+shift.i-1, j, k+shift.k-1, 0, advecting_level );
	       ul = 0.5 * ( AdvectorValueFrLe + AdvectorValueBeLe );
	       if ( DOF_color( i-1, j, k, component ) == FV_BC_LEFT
	       		|| DOF_color( i-1, j, k, component ) 
				== FV_BC_BEHIND_LEFT
	       		|| DOF_color( i-1, j, k, component ) 
				== FV_BC_FRONT_LEFT )
	       {
	         if ( ul > 0. ) fle = ul * AdvectedValueLe;
	         else fle = ul * AdvectedValueC;
	       }
	       else
	       {
		 xl = get_DOF_coordinate( i+shift.i-1, 0, 0 );
		 xL = get_DOF_coordinate( i-1, component, 0 );
		 dxl  = xC - xL;
		 if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
		   cRim12 = AdvectedValueC;
		 else
		 {
		   xR = get_DOF_coordinate( i+1, component, 0 );
		   dxr  = xR - xC;
		   dxCl = xC - xl;
		   cRim12 = AdvectedValueC - ( dxCl / dxr )
			* FV_DiscreteField::SuperBee_phi(thetaC)
			* ( AdvectedValueRi - AdvectedValueC );
		 }
		   
		 dxLl = xl - xL;
                 AdvectedValueLeLe = DOF_value( 
	     	         i-2, j, k, component, advected_level );

	         thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
	   		( AdvectedValueLe - AdvectedValueLeLe ) /
			( AdvectedValueC - AdvectedValueLe ) : 1.e20;
	         cLim12 = AdvectedValueLe + ( dxLl / dxl )
			* FV_DiscreteField::SuperBee_phi(thetaLe)
			* ( AdvectedValueC - AdvectedValueLe );
		 fle = 0.5 * ( ul * ( cRim12 + cLim12 )
	     		- fabs(ul) * ( cRim12 - cLim12 ) );
	       }

	       // Top and Bottom
	       // --------------
	       if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
	       {
	         AdvectorValueTo = AdvectorValueC;
	         AdvectedValueTo = AdvectedValueC;
	       }
	       else
	       {
	         AdvectorValueTo = AdvectingField->DOF_value( i, j+1, k, 
		 	component, advecting_level );
                 AdvectedValueTo = DOF_value( 
	     	       	i, j+1, k, component, advected_level );
	       }

	       if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
	       {
	         AdvectorValueBo = AdvectorValueC;
	         AdvectedValueBo = AdvectedValueC;
	       }
	       else
	       {
                 AdvectorValueBo = AdvectingField->DOF_value( i, j-1, k, 
		 	component, advecting_level );
                 AdvectedValueBo = DOF_value( 
	     	       	i, j-1, k, component, advected_level );
	       }

	       thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
	       		( AdvectedValueC - AdvectedValueBo ) /
			( AdvectedValueTo - AdvectedValueC ) : 1.e20;

	       // Top (Y)
	       AdvectorValueBeTo = AdvectingField->DOF_value(
		     i, j+shift.j, k+shift.k-1,	1, advecting_level );
               AdvectorValueFrTo = AdvectingField->DOF_value( 
		     i, j+shift.j, k+shift.k, 1, advecting_level );
	       vt = 0.5 * ( AdvectorValueBeTo + AdvectorValueFrTo );
	       if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP
	       		|| DOF_color( i, j+1, k, component ) 
				== FV_BC_BEHIND_TOP
	       		|| DOF_color( i, j+1, k, component ) 
				== FV_BC_FRONT_TOP )
	       {
	         if ( vt > 0. ) fto = vt * AdvectedValueC;
	         else fto = vt * AdvectedValueTo;
	       }
	       else
	       {
	         yt = get_DOF_coordinate( j+shift.j, 1, 1 );
		 yT = get_DOF_coordinate( j+1, component, 1 );
		 dyCt = yt - yC;
		 dyt  = yT - yC;
		 cLip12 = AdvectedValueC + ( dyCt / dyt ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaC)
	 		* ( AdvectedValueTo - AdvectedValueC );
	         dyTt = yT - yt;
	         dyT = get_cell_size( j+1, component, 1 );
                 AdvectedValueToTo = DOF_value( 
	     	         i, j+2, k, component, advected_level );
	   
	         thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
		   	1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
			( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
	         cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaTo)
	       		* ( AdvectedValueToTo - AdvectedValueTo );
		 fto = 0.5 * ( vt * ( cRip12 + cLip12 )
	     		- fabs(vt) * ( cRip12 - cLip12 ) );
	       }

	       // Bottom (Y)
	       AdvectorValueBeBo = AdvectingField->DOF_value( 
		     i, j+shift.j-1, k+shift.k-1, 1, advecting_level );
	       AdvectorValueFrBo = AdvectingField->DOF_value( 
		     i, j+shift.j-1, k+shift.k, 1, advecting_level );
	       vb = 0.5 * ( AdvectorValueBeBo + AdvectorValueFrBo );
	       if ( DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
	       		|| DOF_color( i, j-1, k, component ) 
				== FV_BC_BEHIND_BOTTOM
	       		|| DOF_color( i, j-1, k, component ) 
				== FV_BC_FRONT_BOTTOM )
	       {
	         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
	         else fbo = vb * AdvectedValueC;
	       }
	       else
	       {
	         yb = get_DOF_coordinate( j+shift.j-1, 1, 1 );
	         yB = get_DOF_coordinate( j-1, component, 1 );
	         dyb  = yC - yB;
		 if ( DOF_color( i, j, k, component ) == FV_BC_TOP ) 
		   cRim12 = AdvectedValueC;
		 else
		 {
		   yT = get_DOF_coordinate( j+1, component, 1 );
		   dyt  = yT - yC;
		   dyCb = yC - yb;
		   cRim12 = AdvectedValueC - ( dyCb / dyt )
			* FV_DiscreteField::SuperBee_phi(thetaC)
			* ( AdvectedValueTo - AdvectedValueC );
	         }
		 dyBb = yb - yB;
                 AdvectedValueBoBo = DOF_value( 
	     	         i, j-2, k, component, advected_level );
	   
	         thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
	   		( AdvectedValueBo - AdvectedValueBoBo ) /
			( AdvectedValueC - AdvectedValueBo ) : 1.e20;
	         cLim12 = AdvectedValueBo + ( dyBb / dyb )
		   	* FV_DiscreteField::SuperBee_phi(thetaBo)
	       		* ( AdvectedValueC - AdvectedValueBo );
	         fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
	     		- fabs(vb) * ( cRim12 - cLim12 ) );
	       }


	       // Front and Behind
	       // ----------------
	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
	       {
	         AdvectorValueFr = AdvectorValueC;
	         AdvectedValueFr = AdvectedValueC;
	       }
	       else
	       {
	         AdvectorValueFr = AdvectingField->DOF_value( i, j, k+1, 
		 	component, advecting_level );	     
                 AdvectedValueFr = DOF_value( 
	     	       	i, j, k+1, component, advected_level );
	       }

	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
	       {
	         AdvectorValueBe = AdvectorValueC;
	         AdvectedValueBe = AdvectedValueC;
	       }
	       else
	       {
                 AdvectorValueBe = AdvectingField->DOF_value( i, j, k-1, 
		 	component, advecting_level );
                 AdvectedValueBe = DOF_value( 
	     	       	i, j, k-1, component, advected_level );
	       }

	       thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ? 
	       		( AdvectedValueC - AdvectedValueBe ) /
			( AdvectedValueFr - AdvectedValueC ) : 1.e20;
	       
	       // Front (Z)
	       if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
	         ffr = AdvectorValueC * AdvectedValueC;
	       else
	       {	     
	         wf = 0.5 * ( AdvectorValueFr + AdvectorValueC );
	         if ( DOF_color( i, j, k+1, component ) == FV_BC_FRONT )
		 {
	           if ( wf > 0. ) ffr = wf * AdvectedValueC;
	           else ffr = wf * AdvectedValueFr;
		 }
	         else
	         {
		   zf = get_DOF_coordinate( k+shift.k, 0, 2 );
	           zF = get_DOF_coordinate( k+1, component, 2 );
	           dzCf = zf - zC;
	           dzf  = zF - zC;
		   cLip12 = AdvectedValueC + ( dzCf / dzf ) 
		   	* FV_DiscreteField::SuperBee_phi(thetaC)
	   		* ( AdvectedValueFr - AdvectedValueC );	   

	           dzFf = zF - zf;
	           dzF = get_cell_size( k+1, component, 2 );
                   AdvectedValueFrFr = DOF_value( 
	     	         i, j, k+2, component, advected_level );

	           thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 
		   	1.e-20 ? ( AdvectedValueFr - AdvectedValueC ) /
			( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
	           cRip12 = AdvectedValueFr - ( dzFf / dzF )
			* FV_DiscreteField::SuperBee_phi(thetaFr)
			* ( AdvectedValueFrFr - AdvectedValueFr );
	           ffr = 0.5 * ( wf * ( cRip12 + cLip12 )
	       		- fabs(wf) * ( cRip12 - cLip12 ) );
		 } 
               }

	       // Behind (Z)
	       if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
	         fbe = AdvectorValueC * AdvectedValueC;
	       else
	       {	     
	         wb = 0.5 * (AdvectorValueBe + AdvectorValueC);
	         if ( DOF_color(i, j, k-1, component ) == FV_BC_BEHIND )
		 {
	           if (wb > 0.) fbe = wb * AdvectedValueBe;
	           else fbe = wb * AdvectedValueC;
		 }
	         else
	         {
	           zb = get_DOF_coordinate( k+shift.k-1, 0, 2 );
	           zB = get_DOF_coordinate( k-1, component, 2 );
	           dzb  = zC - zB;
		   if ( DOF_color( i, j, k, component ) == FV_BC_FRONT ) 
		     cRim12 = AdvectedValueC;
	           else
		   {
	             zF = get_DOF_coordinate( k+1, component, 2 );
		     dzf  = zF - zC;
		     dzCb = zC - zb;
	             cRim12 = AdvectedValueC - ( dzCb / dzf )
			* FV_DiscreteField::SuperBee_phi(thetaC)
			* ( AdvectedValueFr - AdvectedValueC );
		   }
		   dzBb = zb - zB;
                   AdvectedValueBeBe = DOF_value( 
	     	         i, j, k-2, component, advected_level );
	     
	           thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
		 	( AdvectedValueBe - AdvectedValueBeBe ) /
			( AdvectedValueC - AdvectedValueBe ) : 1.e20;
	           cLim12 = AdvectedValueBe + ( dzBb / dzb )
		   	* FV_DiscreteField::SuperBee_phi(thetaBe)
		 	* ( AdvectedValueC - AdvectedValueBe );
	           fbe = 0.5 * ( wb * ( cRim12 + cLim12 )
	       		- fabs(wb) * ( cRim12 - cLim12 ) );
	         }
               }
	     }
	     
             flux = ( fto - fbo ) * dxC * dzC
                  + ( fri - fle ) * dyC * dzC
                  + ( ffr - fbe ) * dxC * dyC;
             VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
	   }
	 }
       }
     }
   }

   // Synchronize vector for parallel usage
   VEC_rhs->synchronize() ;  
   
}   



//----------------------------------------------------------------------
void
FV_DiscreteField_Staggered:: assemble_tauGradv_tensor_divergence_matrix( 
      	FV_DiscreteField const* DD,
	double const& coef, LA_Matrix* MAT ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "FV_DiscreteField_Staggered:: assemble_tauGradv_tensor_divergence_matrix" );

   MAC_ASSERT( DD->discretization_type() == "tensor" ) ;
   
   // Parameters
   size_t_vector min_unknown_index(DIM,0);
   size_t_vector max_unknown_index(DIM,0);
   size_t center_pos_in_matrix = 0, pos_in_matrix = 0, comp ;
   double dxC, dyC, dzC ;
   double arx, alx, aty, aby, afz, abz ;   
   
   FV_SHIFT_TRIPLET shift ;
   
   for (comp=0;comp<NB_COMPS;++comp)
   {
     // Get local min and max indices
     for (size_t l=0;l<DIM;++l) 
       min_unknown_index(l) = 
       	get_min_index_unknown_handled_by_proc( comp, l ) ;
     for (size_t l=0;l<DIM;++l) 
       max_unknown_index(l) = 
       	get_max_index_unknown_handled_by_proc( comp, l ) ;

     // We can use this shift function because 
     shift = shift_staggeredToStaggered( comp ) ;
          
     // Perform assembling
     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {          
       dxC = get_cell_size( i, comp, 0 ) ;      
       for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
	 dyC = get_cell_size( j, comp, 1 ) ; 
 
         if ( DIM == 2 )
	 {
	   size_t k = 0 ;
	   center_pos_in_matrix = DOF_global_number( i, j, k, comp );
	   
	   // The First Component (u)
	   if ( comp == 0 )
	   {
	     // Right (X) - dxx (comp = 0)
	     if ( DOF_color( i, j, k, comp) != FV_BC_RIGHT )
	     {
	       arx = coef * dyC ;
	       pos_in_matrix = DD->DOF_global_number(
	     			i+shift.i, j, k, comp );
	       MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, arx );
	     }
	     	     
	     // Left (X) - dxx (comp = 0)
	     if ( DOF_color( i, j, k, comp) != FV_BC_LEFT )
	     {
	       alx = - coef * dyC ;
	       pos_in_matrix = DD->DOF_global_number(
	     			i+shift.i-1, j, k, comp );
	       MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, alx );
	     }
	     	   
	     // Top (Y) - dxy (comp = 2)
	     if ( DOF_color( i, j, k, comp) != FV_BC_TOP )
	     {
	       aty = coef * dxC ;
	       pos_in_matrix = DD->DOF_global_number(
	     			i, j+shift.j, k, 2 );
	       MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, aty );
	     }	   
	   
	     // Bottom (Y) - dxy (comp = 2)
	     if ( DOF_color( i, j, k, comp) != FV_BC_BOTTOM )
	     {
	       aby = - coef * dxC ;
	       pos_in_matrix = DD->DOF_global_number(
	     			i, j+shift.j-1, k, 2 );
	       MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, aby );
	     }	   
	   }
	   // The second Component (v)
	   else
	   {
	     // Right (X) - dxy (comp = 2)
	     if ( DOF_color( i, j, k, comp) != FV_BC_RIGHT )
	     {
	       arx = coef * dyC ;
	       pos_in_matrix = DD->DOF_global_number(
	     			i+shift.i, j, k, 2 );
	       MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, arx );
	     }	   
	     
	     // Left (X) - dxy (comp = 2)
	     if ( DOF_color( i, j, k, comp) != FV_BC_LEFT )
	     {
	       alx = - coef * dyC ;
	       pos_in_matrix = DD->DOF_global_number(
	     			i+shift.i-1, j, k, 2 );
	       MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, alx );
	     }	   
	   
	     // Top (Y) - dyy (comp = 1)
	     if ( DOF_color( i, j, k, comp) != FV_BC_TOP )
	     {
	       aty = coef * dxC ;
	       pos_in_matrix = DD->DOF_global_number(
	     			i, j+shift.j, k, comp );
	       MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, aty );
	     }	   
	   
	     // Bottom (Y) - dyy (comp = 1)
	     if ( DOF_color( i, j, k, comp) != FV_BC_BOTTOM )
	     {
	       aby = - coef * dxC ;
	       pos_in_matrix = DD->DOF_global_number(
	     			i, j+shift.j-1, k, comp );
	       MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, aby );
	     }	   
	   }
	 }
	 else
	 {
	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	   {
	     dzC = get_cell_size( k, comp, 2 ) ;	     
	     center_pos_in_matrix = DOF_global_number( i, j, k, comp );

	     // The First Component (u)
	     if ( comp == 0 )
	     {
	       // Right (X) - dxx (0)
	       if ( DOF_color( i, j, k, comp) != FV_BC_RIGHT )
	       {
	         arx = coef * dyC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i+shift.i, j, k, comp );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, arx );
	       }
	       // Left (X) - dxx (0)
	       if ( DOF_color( i, j, k, comp) != FV_BC_LEFT )
	       {
	         alx = - coef * dyC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i+shift.i-1, j, k, comp );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, alx );
	       }
	   
	       // Top (Y) - dxy (comp = 3)
	       if ( DOF_color( i, j, k, comp) != FV_BC_TOP )
	       {
	         aty = coef * dxC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j+shift.j, k, 3 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, aty );
	       }
	   
	       // Bottom (Y) - dxy (comp = 3)
	       if ( DOF_color( i, j, k, comp) != FV_BC_BOTTOM )
	       {
	         aby = - coef * dxC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j+shift.j-1, k, 3 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, aby );
	       }

	       // Front (Z) - dxz (comp = 4)
	       if ( DOF_color( i, j, k, comp) != FV_BC_FRONT )
	       {
	         afz = coef * dxC * dyC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j, k+shift.k, 4 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, afz );
	       }
				     
	       // Behind (Z) - dxz (comp = 4)
	       if ( DOF_color( i, j, k, comp) != FV_BC_BEHIND )
	       {
	         abz = - coef * dxC * dyC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j, k+shift.k-1, 4 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, abz );
	       }
	     }
	     // The Second Component (v)
	     else if ( comp == 1 )
	     {
	       // Right (X) - dxy (comp = 3)
	       if ( DOF_color( i, j, k, comp) != FV_BC_RIGHT )
	       {
	         arx = coef * dyC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i+shift.i, j, k, 3 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, arx );
	       }
	 
	       // Left (X) - dxy (comp = 3)
	       if ( DOF_color( i, j, k, comp) != FV_BC_LEFT )
	       {
	         alx = - coef * dyC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i+shift.i-1, j, k, 3 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, alx );
	       }
	   
	       // Top (Y) - dyy (comp = 1)
	       if ( DOF_color( i, j, k, comp) != FV_BC_TOP )
	       {
	         aty = coef * dxC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j+shift.j, k, comp );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, aty );
	       }
	   
	       // Bottom (Y) - dyy (comp = 1)
	       if ( DOF_color( i, j, k, comp) != FV_BC_BOTTOM )
	       {
	         aby = - coef * dxC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j+shift.j-1, k, comp );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, aby );
	       }

	       // Front (Z) - dyz (comp = 5)
	       if ( DOF_color( i, j, k, comp) != FV_BC_FRONT )
	       {
	         afz = coef * dxC * dyC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j, k+shift.k, 5 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, afz );
	       }
				     
	       // Behind (Z) - dyz (comp = 5)
	       if ( DOF_color( i, j, k, comp) != FV_BC_BEHIND )
	       {
	         abz = - coef * dxC * dyC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j, k+shift.k-1, 5 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, abz );
	       }
	     }
	     // The Third Component (w)
	     else
	     {
	       // Right (X) - dxz (comp = 4)
	       if ( DOF_color( i, j, k, comp) != FV_BC_RIGHT )
	       {
	         arx = coef * dyC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i+shift.i, j, k, 4 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, arx );
	       }
	 
	       // Left (X) - dxz (comp = 4)
	       if ( DOF_color( i, j, k, comp) != FV_BC_LEFT )
	       {
	         alx = - coef * dyC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i+shift.i-1, j, k, 4 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, alx );
	       }
	   
	       // Top (Y) - dyz (comp = 5)
	       if ( DOF_color( i, j, k, comp) != FV_BC_TOP )
	       {
	         aty = coef * dxC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j+shift.j, k, 5 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, aty );
	       }
	   
	       // Bottom (Y) - dyz (comp = 5)
	       if ( DOF_color( i, j, k, comp) != FV_BC_BOTTOM )
	       {
	         aby = - coef * dxC * dzC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j+shift.j-1, k, 5 );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, aby );
	       }

	       // Front (Z) - dzz (comp = 2)
	       if ( DOF_color( i, j, k, comp) != FV_BC_FRONT )
	       {
	         afz = coef * dxC * dyC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j, k+shift.k, comp );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, afz );
	       }
				     
	       // Behind (Z) - dzz (comp = 2)
	       if ( DOF_color( i, j, k, comp) != FV_BC_BEHIND )
	       {
	         abz = - coef * dxC * dyC ;
	         pos_in_matrix = DD->DOF_global_number(
	     			i, j, k+shift.k-1, comp );
	         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, abz );
	       }
	     }
	   }
	 }
       }
     }
   }      

   // Synchronize matrix for parallel usage
   MAT->synchronize() ;

}
