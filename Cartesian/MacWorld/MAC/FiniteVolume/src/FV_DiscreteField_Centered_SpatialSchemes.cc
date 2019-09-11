#include <FV_DiscreteField_Centered.hh>
#include <FV_Mesh.hh>
#include <FV_BoundaryCondition.hh>
#include <FV_DomainBuilder.hh>
#include <FV.hh>
#include <MAC_Exec.hh>
#include <MAC_Communicator.hh>
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
#include <cmath>

using std::cout ; 
using std::endl ;
using std::string ; 
using std::ostringstream ;


//----------------------------------------------------------------------
void
FV_DiscreteField_Centered:: assemble_pDivv_matrix( FV_DiscreteField const* UU,
	double const& coef, LA_Matrix *MAT, LA_Vector *VEC_rhs ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: assemble_pDivv_matrix" );   
   MAC_ASSERT( coef == 1. || coef == -1. ) ;
   MAC_ASSERT( UU->discretization_type() == "staggered" ) ;   

   // Parameters
   size_t_vector min_unknown_index(DIM,0);
   size_t_vector max_unknown_index(DIM,0);
   double dxC, dyC, dzC ;
   double arx, alx, aty, aby, afz, abz ;   
   size_t center_pos_in_matrix = 0 ; 
   
   // Get local min and max indices
   for( size_t l=0; l<DIM; ++l )
   {
     min_unknown_index(l) = 
       	get_min_index_unknown_handled_by_proc( 0, l ) ;
     max_unknown_index(l) = 
       	get_max_index_unknown_handled_by_proc( 0, l ) ;
   }
   FV_SHIFT_TRIPLET shift = shift_staggeredToCentered() ;

   // Perform assembling
   for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
   {          
     dxC = get_cell_size( i, 0, 0 ) ;      
     for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
     {
       dyC = get_cell_size( j, 0, 1 ) ; 
 
       if ( DIM == 2 )
       {
	 size_t k = 0 ;
	 center_pos_in_matrix = DOF_global_number( i, j, k, 0 );

	 // Right (X)
	 arx = coef * dyC ;
	 one_DOF_pDivv( UU, MAT, VEC_rhs, i+shift.i, j, k, 0, 
	 	center_pos_in_matrix, arx ) ;
	 
	 // Left (X)
	 alx = - coef * dyC ;
	 one_DOF_pDivv( UU, MAT, VEC_rhs, i+shift.i-1, j, k, 0, 
	 	center_pos_in_matrix, alx );
	   
	 // Top (Y)
	 aty = coef * dxC ;
	 one_DOF_pDivv( UU, MAT, VEC_rhs, i, j+shift.j, k, 1, 
	 	center_pos_in_matrix, aty ) ;
	   
	 // Bottom (Y)
	 aby = - coef * dxC ;
	 one_DOF_pDivv( UU, MAT, VEC_rhs, i, j+shift.j-1, k, 1, 
	 	center_pos_in_matrix, aby );
       }
       else
       {
         for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	 {
	   dzC = UU->get_cell_size( k, 0, 2 ) ;	     
	   center_pos_in_matrix = DOF_global_number( i, j, k, 0 );
	   
	   // Right (X)
	   arx = coef * dyC * dzC ;
	   one_DOF_pDivv( UU, MAT, VEC_rhs, i+shift.i, j, k, 0, 
	   	center_pos_in_matrix, arx ) ;
	 
	   // Left (X)
	   alx = - coef * dyC * dzC ;
	   one_DOF_pDivv( UU, MAT, VEC_rhs, i+shift.i-1, j, k, 0, 
	   	center_pos_in_matrix, alx );
	   
	   // Top (Y)
	   aty = coef * dxC * dzC ;
	   one_DOF_pDivv( UU, MAT, VEC_rhs, i, j+shift.j, k, 1, 
	   	center_pos_in_matrix, aty ) ;
	   
	   // Bottom (Y)
	   aby = - coef * dxC * dzC ;
	   one_DOF_pDivv( UU, MAT, VEC_rhs, i, j+shift.j-1, k, 1, 
	   	center_pos_in_matrix, aby );

	   // Front (Z)
	   afz = coef * dxC * dyC ;
	   one_DOF_pDivv( UU, MAT, VEC_rhs, i, j, k+shift.k, 2, 
	   	center_pos_in_matrix, afz );
				     
	   // Behind (Z)
	   abz = - coef * dxC * dyC ;
	   one_DOF_pDivv( UU, MAT, VEC_rhs, i, j, k+shift.k-1, 2, 
	   	center_pos_in_matrix, abz );
	 }
       }
     } 
   } 

   // Synchronize matrix and vector for parallel usage
   MAT->synchronize() ;
   VEC_rhs->synchronize() ;

}




//---------------------------------------------------------------------------
void
FV_DiscreteField_Centered:: one_DOF_pDivv( FV_DiscreteField const* UU,
	LA_Matrix *MAT, LA_Vector *VEC_rhs, 
	int i, int j, int k, size_t component, 
	size_t center_pos_in_matrix, double ai ) const
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_DiscreteField_Centered:: one_DOF_pDivv" ) ;

   if ( UU->DOF_is_unknown( i, j, k, component ) )
   {
     size_t pos_in_matrix = UU->DOF_global_number( i, j, k, component );      
     MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, ai );
   }
   else if ( UU->DOF_has_imposed_Dirichlet_value( i, j, k, component ) )
   {
     double dirichlet_value = UU->DOF_value( i, j, k, component, 0 ) ;
     VEC_rhs->add_to_item( center_pos_in_matrix, - ai * dirichlet_value );
   }
   else
   {
     ostringstream mesg ;
     mesg << "Assembling problem in MAC_NavierStokes:: one_DOF_pDivv "
     		<< endl;
     MAC_Error::object()->raise_plain( mesg.str() ) ; 
   }     
   
}




// //----------------------------------------------------------------------
// void 
// FV_DiscreteField_Centered:: assemble_advection_Upwind( 
// 	FV_DiscreteField const* AdvectingField,
// 	size_t advecting_level, double const& coef, size_t advected_level,
// 	LA_Vector *VEC_rhs ) const
// //----------------------------------------------------------------------
// {
//    MAC_LABEL( "FV_DiscreteField_Centered:: assemble_advection_Upwind" );   
//    MAC_CHECK_PRE( advected_level < STO_DEPTH ) ;
//    MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
//    MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ;       
//    MAC_ASSERT( NB_COMPS == 1 ) ;
//    
//    // Parameters
//    size_t_vector min_unknown_index(DIM,0);
//    size_t_vector max_unknown_index(DIM,0); 
//    size_t center_pos_in_matrix = 0, component = 0 ;     
//    double dxC, dyC, dzC ;
//    double AdvectedvalueC = 0., AdvectedvalueRi = 0., AdvectedvalueLe = 0., 
//    	AdvectedvalueTo = 0., AdvectedvalueBo = 0.,
//    	AdvectedvalueFr = 0., AdvectedvalueBe = 0,
// 	ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
// 	fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
//    FV_SHIFT_TRIPLET shift ;
//    
//    // Nullify vector
//    VEC_rhs->nullify();
// 
//    // Get local min and max indices
//    for (size_t l=0;l<DIM;++l) 
//    {
//      min_unknown_index(l) = 
//        	(*min_index_unknown_handled_by_proc)[component](l);
//      max_unknown_index(l) =
//         (*max_index_unknown_handled_by_proc)[component](l);
//    }
//    
//    shift = shift_staggeredToCentered() ;
// 
//    // Perform assembling
//    for( size_t i=min_unknown_index(0); i<=max_unknown_index(0); ++i )
//    {          
//      dxC = (*local_cell_size)[component][0](i) ;
//      for( size_t j=min_unknown_index(1); j<=max_unknown_index(1); ++j )
//      {
//        dyC = (*local_cell_size)[component][1](j) ; 
// 
//        if ( DIM == 2 )
//        {
//          size_t k = 0;
// 	 center_pos_in_matrix = DOF_global_number( i, j, k, component );
// 	 AdvectedvalueC = DOF_value( i, j, k, component, advected_level );
// 	   
// 	 // Right (X)
// 	 if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	   fri = 0.;
// 	 else
// 	 {
// 	   AdvectedvalueRi = DOF_value( i+1, j, k, component, advected_level );
// 	   ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, 
// 	   	advecting_level );
// 	   if ( ur > 0. ) fri = ur * AdvectedvalueC;
// 	   else fri = ur * AdvectedvalueRi;
// 	 }
// 	 
// 	 // Left (X)
// 	 if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	   fle = 0.;
// 	 else
// 	 {
// 	   AdvectedvalueLe = DOF_value( i-1, j, k, component, advected_level );
// 	   ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, 
// 	   	advecting_level );
// 	   if ( ul > 0. ) fle = ul * AdvectedvalueLe;
// 	   else fle = ul * AdvectedvalueC;
// 	 }
// 	 
// 	 // Top (Y)
// 	 if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	   fto = 0.;
// 	 else
// 	 {
// 	   AdvectedvalueTo = DOF_value( i, j+1, k, component, advected_level );
// 	   vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, 
// 	   	advecting_level );
// 	   if ( vt > 0. ) fto = vt * AdvectedvalueC;
// 	   else fto = vt * AdvectedvalueTo;
// 	 }
// 	 
// 	 // Bottom (Y)
// 	 if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	   fbo = 0.;
// 	 else
// 	 {
// 	   AdvectedvalueBo = DOF_value( i, j-1, k, component, advected_level );
// 	   vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, 
// 	   	advecting_level );
// 	   if ( vb > 0. ) fbo = vb * AdvectedvalueBo;
// 	   else fbo = vb * AdvectedvalueC;
// 	 }
// 	 
//          flux = (fto - fbo) * dxC + (fri - fle) * dyC;	 
// 	 VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
//        }
//        else
//        {
//          for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
// 	 {
// 	   dzC = (*local_cell_size)[component][2](k) ; 
// 	   center_pos_in_matrix = DOF_global_number( i, j, k, component );
// 	   AdvectedvalueC = DOF_value( i, j, k, component, advected_level );
// 	   
// 	   // Right (X)
// 	   if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT )
// 	     fri = 0.;
// 	   else
// 	   {
// 	     AdvectedvalueRi = DOF_value( i+1, j, k, component, 
// 	     	advected_level );
// 	     ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, 
// 	     	advecting_level );
// 	     if ( ur > 0. ) fri = ur * AdvectedvalueC;
// 	     else fri = ur * AdvectedvalueRi;
// 	   }
// 	 
// 	   // Left (X)
// 	   if ( DOF_color( i, j, k, component ) == FV_BC_LEFT )
// 	     fle = 0.;
// 	   else
// 	   {
// 	     AdvectedvalueLe = DOF_value( i-1, j, k, component, 
// 	     	advected_level );
// 	     ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, 
// 	     	advecting_level );
// 	     if ( ul > 0. ) fle = ul * AdvectedvalueLe;
// 	     else fle = ul * AdvectedvalueC;
// 	   }
// 	 
// 	   // Top (Y)
// 	   if ( DOF_color( i, j, k, component ) == FV_BC_TOP )
// 	     fto = 0.;
// 	   else
// 	   {
// 	     AdvectedvalueTo = DOF_value( i, j+1, k, component, 
// 	     	advected_level );
// 	     vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, 
// 	     	advecting_level );
// 	     if ( vt > 0. ) fto = vt * AdvectedvalueC;
// 	     else fto = vt * AdvectedvalueTo;
// 	   }
// 	 
// 	   // Bottom (Y)
// 	   if ( DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
// 	     fbo = 0.;
// 	   else
// 	   {
// 	     AdvectedvalueBo = DOF_value( i, j-1, k, component, 
// 	     	advected_level );
// 	     vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, 
// 	     	advecting_level );
// 	     if ( vb > 0. ) fbo = vb * AdvectedvalueBo;
// 	     else fbo = vb * AdvectedvalueC;
// 	   }
// 
// 	   // Front (Z)
// 	   if ( DOF_color( i, j, k, component ) == FV_BC_FRONT )
// 	     ffr = 0.;
// 	   else
// 	   {
// 	     AdvectedvalueFr = DOF_value( i, j, k+1, component, 
// 	     	advected_level );
// 	     wf = AdvectingField->DOF_value( i, j, k+shift.k, 2, 
// 	     	advecting_level );
// 	     if ( wf > 0. ) ffr = wf * AdvectedvalueC;
// 	     else ffr = wf * AdvectedvalueFr;
// 	   }
// 	 
// 	   // Behind (Z)
// 	   if ( DOF_color( i, j, k, component ) == FV_BC_BEHIND )
// 	     fbe = 0.;
// 	   else
// 	   {
// 	     AdvectedvalueBe = DOF_value( i, j, k-1, component, 
// 	     	advected_level );
// 	     wb = AdvectingField->DOF_value( i, j, k+shift.k-1, 2, 
// 	     	advecting_level );
// 	     if ( wb > 0. ) fbe = wb * AdvectedvalueBe;
// 	     else fbe = wb * AdvectedvalueC;
// 	   }
// 	 
//            flux = (fto - fbo) * dxC * dzC
// 	   	+ (fri - fle) * dyC * dzC
// 		+ (ffr - fbe) * dxC * dyC;
// 	   VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
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
FV_DiscreteField_Centered:: assemble_advection_Upwind( 
	FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t advected_level,
	LA_Vector *VEC_rhs ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: assemble_advection_Upwind" );   
   MAC_CHECK_PRE( advected_level < STO_DEPTH ) ;
   MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
   MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ;       
   MAC_ASSERT( NB_COMPS == 1 ) ;
   
   // Parameters
   size_t_vector min_unknown_index(DIM,0);
   size_t_vector max_unknown_index(DIM,0); 
   size_t center_pos_in_matrix = 0, component = 0 ;     
   double dxC = 0., dyC = 0., dzC = 0.;
   double AdvectedvalueC = 0., AdvectedvalueRi = 0., AdvectedvalueLe = 0., 
   	AdvectedvalueTo = 0., AdvectedvalueBo = 0.,
   	AdvectedvalueFr = 0., AdvectedvalueBe = 0,
	ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
	fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
   FV_SHIFT_TRIPLET shift ;

   // Comment: cell centered unknowns always have a defined value at +1/-1
   // indices in all 3 directions. Whether one of the +1/-1 DOF values is on a
   // boundary or not, and whether that boundary has a Dirichlet or Neumann 
   // condition is irrelevant, this +1/-1 DOF always has the right value. 
   // For Neumann, this is guaranted by 
   // FV_BoundaryCondition:: set_free_DOF_values in 
   // FV_DiscreteField:: update_free_DOFs_value or 
   // FV_DiscreteField:: add_to_free_DOFs_value 
   
   // Nullify vector
   VEC_rhs->nullify();

   // Get local min and max indices
   for (size_t l=0;l<DIM;++l) 
   {
     min_unknown_index(l) = 
       	(*min_index_unknown_handled_by_proc)[component](l);
     max_unknown_index(l) =
        (*max_index_unknown_handled_by_proc)[component](l);
   }
   
   shift = shift_staggeredToCentered() ;

   // Perform assembling
   for( size_t i=min_unknown_index(0); i<=max_unknown_index(0); ++i )
   {          
     dxC = (*local_cell_size)[component][0](i) ;
     for( size_t j=min_unknown_index(1); j<=max_unknown_index(1); ++j )
     {
       dyC = (*local_cell_size)[component][1](j) ; 

       if ( DIM == 2 )
       {
         size_t k = 0;
	 center_pos_in_matrix = DOF_global_number( i, j, k, component );
	 AdvectedvalueC = DOF_value( i, j, k, component, advected_level );
	   
	 // Right (X)
	 AdvectedvalueRi = DOF_value( i+1, j, k, component, advected_level );
	 ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, 
	   	advecting_level );
	 if ( ur > 0. ) fri = ur * AdvectedvalueC;
	 else fri = ur * AdvectedvalueRi;
	 
	 // Left (X)
	 AdvectedvalueLe = DOF_value( i-1, j, k, component, advected_level );
	 ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, 
	   	advecting_level );
	 if ( ul > 0. ) fle = ul * AdvectedvalueLe;
	 else fle = ul * AdvectedvalueC;
	 
	 // Top (Y)
	 AdvectedvalueTo = DOF_value( i, j+1, k, component, advected_level );
	 vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, 
	   	advecting_level );
	 if ( vt > 0. ) fto = vt * AdvectedvalueC;
	 else fto = vt * AdvectedvalueTo;
	 
	 // Bottom (Y)
	 AdvectedvalueBo = DOF_value( i, j-1, k, component, advected_level );
	 vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, 
	   	advecting_level );
	 if ( vb > 0. ) fbo = vb * AdvectedvalueBo;
	 else fbo = vb * AdvectedvalueC;
	 
         flux = (fto - fbo) * dxC + (fri - fle) * dyC;	 
	 VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
       }
       else
       {
         for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	 {
	   dzC = (*local_cell_size)[component][2](k) ; 
	   center_pos_in_matrix = DOF_global_number( i, j, k, component );
	   AdvectedvalueC = DOF_value( i, j, k, component, advected_level );
	   
	   // Right (X)
	   AdvectedvalueRi = DOF_value( i+1, j, k, component, 
	     	advected_level );
	   ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, 
	     	advecting_level );
	   if ( ur > 0. ) fri = ur * AdvectedvalueC;
	   else fri = ur * AdvectedvalueRi;
	 
	   // Left (X)
	   AdvectedvalueLe = DOF_value( i-1, j, k, component, 
	     	advected_level );
	   ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, 
	     	advecting_level );
	   if ( ul > 0. ) fle = ul * AdvectedvalueLe;
	   else fle = ul * AdvectedvalueC;

	   // Top (Y)
	   AdvectedvalueTo = DOF_value( i, j+1, k, component, 
	     	advected_level );
	   vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, 
	     	advecting_level );
	   if ( vt > 0. ) fto = vt * AdvectedvalueC;
	   else fto = vt * AdvectedvalueTo;
	 
	   // Bottom (Y)
	   AdvectedvalueBo = DOF_value( i, j-1, k, component, 
	     	advected_level );
	   vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, 
	     	advecting_level );
	   if ( vb > 0. ) fbo = vb * AdvectedvalueBo;
	   else fbo = vb * AdvectedvalueC;
	   
	   // Front (Z)
	   AdvectedvalueFr = DOF_value( i, j, k+1, component, 
	     	advected_level );
	   wf = AdvectingField->DOF_value( i, j, k+shift.k, 2, 
	     	advecting_level );
	   if ( wf > 0. ) ffr = wf * AdvectedvalueC;
	   else ffr = wf * AdvectedvalueFr;
	 
	   // Behind (Z)
	   AdvectedvalueBe = DOF_value( i, j, k-1, component, 
	     	advected_level );
	   wb = AdvectingField->DOF_value( i, j, k+shift.k-1, 2, 
	     	advecting_level );
	   if ( wb > 0. ) fbe = wb * AdvectedvalueBe;
	   else fbe = wb * AdvectedvalueC;
	 
           flux = (fto - fbo) * dxC * dzC
	   	+ (fri - fle) * dyC * dzC
		+ (ffr - fbe) * dxC * dyC;
	   VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;
	 }
       } 
     }      	       	   
   }   

   // Synchronize vector for parallel usage
   VEC_rhs->synchronize() ;
    
}	




//----------------------------------------------------------------------
void 
FV_DiscreteField_Centered:: assemble_advection_TVD( 
	FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t advected_level,
	LA_Vector *VEC_rhs ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Centered:: assemble_advection_TVD" );   
   MAC_CHECK_PRE( advected_level < STO_DEPTH ) ;
   MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
   MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ; 
   
   // Parameters
   size_t_vector min_unknown_index(DIM,0);
   size_t_vector max_unknown_index(DIM,0); 
   size_t center_pos_in_matrix = 0, component = 0 ;  
   double xC = 0., yC = 0., zC = 0., 
   	xr = 0., xR = 0., xl = 0., xL = 0., yt = 0., yT = 0., yb = 0., yB = 0.,
	zf = 0., zF = 0., zb = 0., zB = 0.;
   double dxC = 0., dyC = 0., dzC = 0., 
   	dxr = 0., dxl = 0., dxCr = 0., dxCl = 0., dxRr = 0., dxR = 0., 
	dxLl = 0., dyt = 0., dyb = 0., dyCt = 0., dyCb = 0., dyTt = 0., 
	dyT = 0., dyBb = 0., dzf = 0., dzb = 0., dzCf = 0., dzCb = 0., 
	dzFf = 0., dzF = 0., dzBb = 0.;

   double AdvectedvalueC = 0., AdvectedvalueRi = 0., AdvectedvalueRiRi = 0., 
   	AdvectedvalueLe = 0.,  AdvectedvalueLeLe = 0., AdvectedvalueTo = 0., 
	AdvectedvalueToTo = 0., AdvectedvalueBo = 0., AdvectedvalueBoBo = 0.,
   	AdvectedvalueFr = 0., AdvectedvalueFrFr = 0., AdvectedvalueBe = 0., 
	AdvectedvalueBeBe = 0.;
   double ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
	fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
   double cRip12 = 0., cLip12 = 0., cRim12 = 0., cLim12 = 0.,
   	thetaC = 0., thetaRi = 0., thetaLe = 0., thetaTo = 0., thetaBo = 0., 
	thetaFr = 0., thetaBe = 0.;
   FV_SHIFT_TRIPLET shift ;
      
   // Nullify vector
   VEC_rhs->nullify();

   // Get local min and max indices
   for (size_t l=0;l<DIM;++l) 
   {
     min_unknown_index(l) = 
       	(*min_index_unknown_handled_by_proc)[component](l);
     max_unknown_index(l) =
        (*max_index_unknown_handled_by_proc)[component](l);
   }
   
   shift = shift_staggeredToCentered() ;
	
   // Perform assembling
   for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
   {          
     xC = get_DOF_coordinate( i, component, 0 );
     dxC = (*local_cell_size)[component][0](i) ;
     for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
     {
       yC = get_DOF_coordinate( j, component, 1 );
       dyC = (*local_cell_size)[component][1](j) ; 

       if ( DIM == 2 )
       {
         size_t k = 0 ;
	 center_pos_in_matrix = DOF_global_number( i, j, k, component );
	 AdvectedvalueC = DOF_value( i, j, k, component, advected_level );
	 
	 // Right and Left
	 // --------------
	 AdvectedvalueRi = DOF_value( i+1, j, k, component, advected_level );
	 AdvectedvalueLe = DOF_value( i-1, j, k, component, advected_level );
	 	 
	 thetaC = fabs( AdvectedvalueRi - AdvectedvalueC ) > 1.e-20  ? 
	 	( AdvectedvalueC - AdvectedvalueLe ) 
		/ ( AdvectedvalueRi - AdvectedvalueC ) : 1e20 ;

	 // Right (X)
	 ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );

	 if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT )
	 {
	   if ( ur > 0. ) fri = ur * AdvectedvalueC;
	   else fri = ur * AdvectedvalueRi;
	 }
	 else
	 {
	   xr = AdvectingField->get_DOF_coordinate( i+shift.i, 0, 0 );
	   xR = get_DOF_coordinate( i+1, component, 0 );
	   dxCr = xr - xC;
	   dxr  = xR - xC;
	   dxRr = xR - xr;
	   dxR = get_cell_size( i+1, component, 0 );
	   AdvectedvalueRiRi = DOF_value( i+2, j, k, component, 
	   	advected_level );
	     
	   thetaRi = fabs( AdvectedvalueRiRi - AdvectedvalueRi) > 1.e-20 ? 
		( AdvectedvalueRi - AdvectedvalueC) 
		/ ( AdvectedvalueRiRi - AdvectedvalueRi) : 1e20 ;
	   cRip12 = AdvectedvalueRi
		- ( dxRr / dxR ) * FV_DiscreteField::SuperBee_phi(thetaRi)
		* ( AdvectedvalueRiRi - AdvectedvalueRi );
           cLip12 = AdvectedvalueC + ( dxCr / dxr ) 
	   	* FV_DiscreteField::SuperBee_phi( thetaC )
		* ( AdvectedvalueRi - AdvectedvalueC );

	   fri = 0.5 * ( ur * ( cRip12 + cLip12 )
		- fabs(ur) * ( cRip12 - cLip12 ) ) ;
	 }
	 
	 // Left (X)
	 ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, 
	 	advecting_level );
	   
	 if ( DOF_color( i-1, j, k, component ) == FV_BC_LEFT )
	 {
	   if ( ul > 0. ) fle = ul * AdvectedvalueLe;
	   else fle = ul * AdvectedvalueC;
	 }
	 else
	 {
	   xl = AdvectingField->get_DOF_coordinate( i+shift.i-1, 0, 0 );
	   xL = get_DOF_coordinate( i-1, component, 0 );
	   dxCl = xC - xl;
	   dxl  = xC - xL;
	   dxLl = xl - xL;
	   AdvectedvalueLeLe = DOF_value( i-2, j, k, component, 
	   	advected_level );
	     
	   thetaLe = fabs( AdvectedvalueC - AdvectedvalueLe ) > 1.e-20 ?
			( AdvectedvalueLe - AdvectedvalueLeLe ) 
			/ ( AdvectedvalueC - AdvectedvalueLe ) : 1e20 ;
	   cLim12 = AdvectedvalueLe
		+ ( dxLl / dxl ) * FV_DiscreteField::SuperBee_phi( thetaLe )
		* ( AdvectedvalueC - AdvectedvalueLe ) ;
	   if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
	     cRim12 = AdvectedvalueC;
	   else
	   {
	     xR = get_DOF_coordinate( i+1, component, 0 );
	     dxr  = xR - xC;
	     cRim12 = AdvectedvalueC - ( dxCl / dxr ) 
	     	* FV_DiscreteField::SuperBee_phi( thetaC ) 
		* ( AdvectedvalueRi - AdvectedvalueC ) ;
	   }

	   fle = 0.5 * ( ul * ( cRim12 + cLim12 )
			- fabs(ul) * ( cRim12 - cLim12 ) ) ;
	 }
	 
	 // Top and Bottom
	 // --------------
	 AdvectedvalueTo = DOF_value( i, j+1, k, component, advected_level );
	 AdvectedvalueBo = DOF_value( i, j-1, k, component, advected_level );

	 thetaC = fabs( AdvectedvalueTo - AdvectedvalueC ) > 1.e-20 ? 
	   	( AdvectedvalueC - AdvectedvalueBo ) 
		/ ( AdvectedvalueTo - AdvectedvalueC) : 1e20 ;

	 // Top (Y)
	 vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );
	   
	 if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP )
	 {
	   if ( vt > 0. ) fto = vt * AdvectedvalueC;
	   else fto = vt * AdvectedvalueTo;
	 }
	 else
	 {
	   yt = AdvectingField->get_DOF_coordinate( j+shift.j, 1, 1 );
	   yT = get_DOF_coordinate( j+1, component, 1 );
	   dyCt = yt - yC;
	   dyt  = yT - yC;
	   dyTt = yT - yt;
	   dyT = get_cell_size( j+1, component, 1 );
           AdvectedvalueToTo = DOF_value( i, j+2, k, component, 
	   	advected_level );
	     
	   thetaTo = fabs( AdvectedvalueToTo - AdvectedvalueTo ) > 1.e-20 ? 
		( AdvectedvalueTo - AdvectedvalueC ) 
		/ ( AdvectedvalueToTo - AdvectedvalueTo) : 1e20 ;
	   cRip12 = AdvectedvalueTo
		- ( dyTt / dyT ) * FV_DiscreteField::SuperBee_phi( thetaTo )
		* ( AdvectedvalueToTo - AdvectedvalueTo );   
           cLip12 = AdvectedvalueC + ( dyCt / dyt ) 
	   	* FV_DiscreteField::SuperBee_phi( thetaC )
		* ( AdvectedvalueTo - AdvectedvalueC );
   
	   fto = 0.5 * ( vt * ( cRip12 + cLip12 )
	       		- fabs(vt) * ( cRip12 - cLip12 ) ) ;
	 }

	 // Bottom (Y)
	 vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, 
	 	advecting_level );
	   
	 if ( DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM )
	 {
	   if ( vb > 0. ) fbo = vb * AdvectedvalueBo;
	   else fbo = vb * AdvectedvalueC;
	 }
	 else
	 {
	   yb = AdvectingField->get_DOF_coordinate( j+shift.j-1, 1, 1 );
	   yB = get_DOF_coordinate( j-1, component, 1 );
	   dyCb = yC - yb;
	   dyb  = yC - yB;
	   dyBb = yb - yB;
	   AdvectedvalueBoBo = DOF_value( i, j-2, k, component, 
	   	advected_level );
	     
	   thetaBo = fabs( AdvectedvalueC - AdvectedvalueBo ) > 1.e-20 ?
		( AdvectedvalueBo - AdvectedvalueBoBo ) 
		/ ( AdvectedvalueC - AdvectedvalueBo ) : 1e20 ;
	   cLim12 = AdvectedvalueBo
		+ ( dyBb / dyb ) * FV_DiscreteField::SuperBee_phi( thetaBo )
		* ( AdvectedvalueC - AdvectedvalueBo );
	   if ( DOF_color( i, j, k, component ) == FV_BC_TOP ) 
	     cRim12 = AdvectedvalueC;
	   else
	   {
	     yT = get_DOF_coordinate( j+1, component, 1 );
	     dyt  = yT - yC;
	     cRim12 = AdvectedvalueC - ( dyCb / dyt ) 
	     	* FV_DiscreteField::SuperBee_phi( thetaC ) 
		* ( AdvectedvalueTo - AdvectedvalueC );
	   }
	     
	   fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
	       		- fabs(vb) * ( cRim12 - cLim12 ) ) ;
	 }

         flux = ( fto - fbo ) * dxC + ( fri - fle ) * dyC;	 
	 VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;	 
       }
       else
       {
         for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	 {
	   zC = get_DOF_coordinate( k, component, 2 );
	   dzC = (*local_cell_size)[component][2](k) ; 
	   center_pos_in_matrix = DOF_global_number( i, j, k, component );
	   AdvectedvalueC = DOF_value( i, j, k, component, advected_level );
	 
	   // Right and Left
	   // --------------
	   AdvectedvalueRi = DOF_value( i+1, j, k, component, advected_level );
	   AdvectedvalueLe = DOF_value( i-1, j, k, component, advected_level );
	 	 
	   thetaC = fabs( AdvectedvalueRi - AdvectedvalueC) > 1.e-20 ? 
	 	( AdvectedvalueC - AdvectedvalueLe ) 
		/ ( AdvectedvalueRi - AdvectedvalueC) : 1e20 ;

	   // Right (X)
	   ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, 
	   	advecting_level );

	   if ( DOF_color( i+1, j, k, component ) == FV_BC_RIGHT )
	   {
	     if ( ur > 0. ) fri = ur * AdvectedvalueC;
	     else fri = ur * AdvectedvalueRi;
	   }
	   else
	   {
	     xr = AdvectingField->get_DOF_coordinate( i+shift.i, 0, 0 );
	     xR = get_DOF_coordinate( i+1, component, 0 );
	     dxCr = xr - xC;
	     dxr  = xR - xC;
	     dxRr = xR - xr;
	     dxR = get_cell_size( i+1, component, 0 );
	     AdvectedvalueRiRi = DOF_value( i+2, j, k, component, 
	     	advected_level );
	     
	     thetaRi = fabs( AdvectedvalueRiRi - AdvectedvalueRi ) > 1.e-20 ? 
		( AdvectedvalueRi - AdvectedvalueC ) \
		/ ( AdvectedvalueRiRi - AdvectedvalueRi ) : 1e20 ;
	     cRip12 = AdvectedvalueRi
		- ( dxRr / dxR ) * FV_DiscreteField::SuperBee_phi( thetaRi )
		* ( AdvectedvalueRiRi - AdvectedvalueRi );
             cLip12 = AdvectedvalueC + ( dxCr / dxr ) 
	     	* FV_DiscreteField::SuperBee_phi( thetaC )
		* ( AdvectedvalueRi - AdvectedvalueC );

	     fri = 0.5 * ( ur * ( cRip12 + cLip12 )
	   		- fabs(ur) * ( cRip12 - cLip12 ) ) ;
	   }
	 
	   // Left (X)
	   ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, 
	   	advecting_level );
	   
	   if ( DOF_color( i-1, j, k, component ) == FV_BC_LEFT )
	   {
	     if ( ul > 0. ) fle = ul * AdvectedvalueLe;
	     else fle = ul * AdvectedvalueC;
	   }
	   else
	   {
	     xl = AdvectingField->get_DOF_coordinate( i+shift.i-1, 0, 0 );
	     xL = get_DOF_coordinate( i-1, component, 0 );
	     dxCl = xC - xl;
	     dxl  = xC - xL;
	     dxLl = xl - xL;
	     AdvectedvalueLeLe = DOF_value( i-2, j, k, component, 
	     	advected_level );
	     
	     thetaLe = fabs( AdvectedvalueC - AdvectedvalueLe) > 1.e-20 ?
		( AdvectedvalueLe - AdvectedvalueLeLe ) 
		/ ( AdvectedvalueC - AdvectedvalueLe) : 1e20 ;
	     cLim12 = AdvectedvalueLe
		+ ( dxLl / dxl ) * FV_DiscreteField::SuperBee_phi( thetaLe )
		* ( AdvectedvalueC - AdvectedvalueLe ) ;
	     if ( DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
	       cRim12 = AdvectedvalueC;
	     else
	     {
	       xR = get_DOF_coordinate( i+1, advected_level, 0 );
	       dxr  = xR - xC;
	       cRim12 = AdvectedvalueC - ( dxCl / dxr ) 
			* FV_DiscreteField::SuperBee_phi( thetaC )
			* ( AdvectedvalueRi - AdvectedvalueC );
	     }

	     fle = 0.5 * ( ul * ( cRim12 + cLim12 )
	       		- fabs(ul) * ( cRim12 - cLim12 ) ) ;
	   }
	 
	   // Top and Bottom
	   // --------------
	   AdvectedvalueTo = DOF_value( i, j+1, k, component, advected_level );
	   AdvectedvalueBo = DOF_value( i, j-1, k, component, advected_level );

	   thetaC = fabs( AdvectedvalueTo - AdvectedvalueC ) > 1.e-20 ? 
	   	( AdvectedvalueC - AdvectedvalueBo ) 
		/ ( AdvectedvalueTo - AdvectedvalueC ) : 1e20 ;

	   // Top (Y)
	   vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, 
	   	advecting_level );
	   
	   if ( DOF_color( i, j+1, k, component ) == FV_BC_TOP )
	   {
	     if ( vt > 0. ) fto = vt * AdvectedvalueC;
	     else fto = vt * AdvectedvalueTo;
	   }
	   else
	   {
	     yt = AdvectingField->get_DOF_coordinate( j+shift.j, 1, 1 );
	     yT = get_DOF_coordinate( j+1, component, 1 );
	     dyCt = yt - yC;
	     dyt  = yT - yC;
	     dyTt = yT - yt;
	     dyT = get_cell_size( j+1, component, 1 );
             AdvectedvalueToTo = DOF_value( i, j+2, k, component, 
	     	advected_level );
	     
	     thetaTo = fabs( AdvectedvalueToTo - AdvectedvalueTo) > 1.e-20 ? 
		( AdvectedvalueTo - AdvectedvalueC ) 
		/ ( AdvectedvalueToTo - AdvectedvalueTo ) : 1e20 ;
	     cRip12 = AdvectedvalueTo
		- ( dyTt / dyT ) * FV_DiscreteField::SuperBee_phi( thetaTo )
		* ( AdvectedvalueToTo - AdvectedvalueTo ) ;   
             cLip12 = AdvectedvalueC + ( dyCt / dyt ) 
		* FV_DiscreteField::SuperBee_phi(thetaC)
		* ( AdvectedvalueTo - AdvectedvalueC ) ;
   
	     fto = 0.5 * ( vt * ( cRip12 + cLip12 )
	       		- fabs(vt) * ( cRip12 - cLip12 ) ) ;
	   }

	   // Bottom (Y)
	   vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, 
	   	advecting_level );
	   
	   if ( DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM )
	   {
	     if ( vb > 0. ) fbo = vb * AdvectedvalueBo;
	     else fbo = vb * AdvectedvalueC;
	   }
	   else
	   {
	     yb = AdvectingField->get_DOF_coordinate( j+shift.j-1, 1, 1 );
	     yB = get_DOF_coordinate( j-1, component, 1 );
	     dyCb = yC - yb;
	     dyb  = yC - yB;
	     dyBb = yb - yB;
	     AdvectedvalueBoBo = DOF_value( i, j-2, k, component, 
	     	advected_level );
	     
	     thetaBo = fabs( AdvectedvalueC - AdvectedvalueBo) > 1.e-20 ?
		( AdvectedvalueBo - AdvectedvalueBoBo )
		/ ( AdvectedvalueC - AdvectedvalueBo ) : 1e20 ;
	     cLim12 = AdvectedvalueBo
		+ ( dyBb / dyb ) * FV_DiscreteField::SuperBee_phi( thetaBo )
		* ( AdvectedvalueC - AdvectedvalueBo ) ;
	     if ( DOF_color( i, j, k, component ) == FV_BC_TOP ) 
	       cRim12 = AdvectedvalueC;
	     else
	     {
	       yT = get_DOF_coordinate( j+1, component, 1 );
	       dyt  = yT - yC;
	       cRim12 = AdvectedvalueC - ( dyCb / dyt ) 
			* FV_DiscreteField::SuperBee_phi( thetaC )
			* ( AdvectedvalueTo - AdvectedvalueC ) ;
	     }
	     
	     fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
	       		- fabs(vb) * ( cRim12 - cLim12 ) ) ;
	   }

	 
	   // Front and Behind
	   // ----------------
	   AdvectedvalueFr = DOF_value( i, j, k+1, component, advected_level );
	   AdvectedvalueBe = DOF_value( i, j, k-1, component, advected_level );

	   thetaC = fabs( AdvectedvalueFr - AdvectedvalueC) > 1.e-20 ? 
	       	( AdvectedvalueC - AdvectedvalueBe )
		/ ( AdvectedvalueFr - AdvectedvalueC ) : 1e20 ;

	   // Front (Z)
	   wf = AdvectingField->DOF_value( i, j, k+shift.k, 2, 
	   	advecting_level );
	   
	   if ( DOF_color( i, j, k+1, component ) == FV_BC_FRONT )
	   {
	     if ( wf > 0. ) ffr = wf * AdvectedvalueC;
	     else ffr = wf * AdvectedvalueFr;
	   }
	   else
	   {
	     zf = AdvectingField->get_DOF_coordinate( k+shift.k, 2, 2 );
	     zF = get_DOF_coordinate( k+1, component, 2 );
	     dzCf = zf - zC;
	     dzf  = zF - zC;
	     dzFf = zF - zf;
	     dzF = get_cell_size( k+1, component, 2 );
             AdvectedvalueFrFr = DOF_value( i, j, k+2, component, 
	     	advected_level );
	     
	     thetaFr = fabs( AdvectedvalueFrFr - AdvectedvalueFr) > 1.e-20 ? 
		( AdvectedvalueFr - AdvectedvalueC )
		/ ( AdvectedvalueFrFr - AdvectedvalueFr ) : 1e20 ;
	     cRip12 = AdvectedvalueFr
		- ( dzFf / dzF ) * FV_DiscreteField::SuperBee_phi( thetaFr )
		* ( AdvectedvalueFrFr - AdvectedvalueFr ) ;   
             cLip12 = AdvectedvalueC + ( dzCf / dzf ) 
		* FV_DiscreteField::SuperBee_phi( thetaC )
		* ( AdvectedvalueFr - AdvectedvalueC ) ;
   
	     ffr = 0.5 * ( wf * ( cRip12 + cLip12 )
	       		- fabs(wf) * ( cRip12 - cLip12 ) ) ;
	   }

	   // Behind (Z)
	   wb = AdvectingField->DOF_value( i, j, k+shift.k-1, 2, 
	   	advecting_level );
	   
	   if ( DOF_color( i, j, k-1, component ) == FV_BC_BEHIND )
	   {
	     if ( wb > 0. ) fbe = wb * AdvectedvalueBe;
	     else fbe = wb * AdvectedvalueC;
	   }
	   else
	   {
	     zb = AdvectingField->get_DOF_coordinate( k+shift.k-1, 2, 2 );
	     zB = get_DOF_coordinate( k-1, component, 2 );
	     dzCb = zC - zb;
	     dzb  = zC - zB;
	     dzBb = zb - zB;
	     AdvectedvalueBeBe = DOF_value( i, j, k-2, component, 
	     	advected_level );
	     
	     thetaBe = fabs( AdvectedvalueC - AdvectedvalueBe) > 1.e-20 ?
		( AdvectedvalueBe - AdvectedvalueBeBe )
		/ ( AdvectedvalueC - AdvectedvalueBe ) : 1e20 ;
	     cLim12 = AdvectedvalueBe
		+ ( dzBb / dzb ) * FV_DiscreteField::SuperBee_phi( thetaBe )
		* ( AdvectedvalueC - AdvectedvalueBe ) ;
	     if ( DOF_color( i, j, k, component ) == FV_BC_FRONT ) 
	       cRim12 = AdvectedvalueC;
	     else
	     {
	       zF = get_DOF_coordinate( k+1, component, 2 );
	       dzf  = zF - zC;
	       cRim12 = AdvectedvalueC - ( dzCb / dzf ) 
	       		* FV_DiscreteField::SuperBee_phi( thetaC )
			* ( AdvectedvalueFr - AdvectedvalueC ) ;
	     }
	     
	     fbe = 0.5 * ( wb * ( cRim12 + cLim12 )
	       		- fabs(wb) * ( cRim12 - cLim12 ) ) ;
	   }

           flux = ( fto - fbo ) * dxC * dzC
	   	+ ( fri - fle ) * dyC * dzC
		+ ( ffr - fbe ) * dxC * dyC;
	   VEC_rhs->set_item( center_pos_in_matrix, coef * flux ) ;	 
	 }
       } 
     }      	       	   
   }
   
   // Synchronize vector for parallel usage
   VEC_rhs->synchronize() ;
         
}
