#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <FV_BoundaryCondition.hh>
#include <FV_DomainBuilder.hh>
#include <FV_TimeIterator.hh>
#include <FV.hh>
#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Data.hh>
#include <MAC_DoubleArray3D.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Int.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_ObjectWriter.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_String.hh>
#include <MAC_Vector.hh>
#include <MAC_Root.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Variable.hh>
#include <stringVector.hh>
#include <LA_SeqVector.hh>
#include <LA_Vector.hh>
#include <LA_Matrix.hh>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

using std::cout ; 
using std::endl ;
using std::string ; 
using std::ostringstream ;


//----------------------------------------------------------------------------
void FV_DiscreteField::assemble_constantcoef_laplacian_matrix( 
	double const& coef_lap,
	LA_Matrix *MAT, LA_Vector *VEC_rhs,
	bool const& rescale ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: assemble_constantcoef_laplacian_matrix" ) ;
  
   // Parameters
   size_t_vector min_unknown_index(DIM,0);
   size_t_vector max_unknown_index(DIM,0);   
   double xC, xR, xL, yC, yT, yB, dxr, dxl, dyt, dyb, dxC, dyC,
   	zC, zF, zB, dzf, dzb, dzC ;
   double ac, arx, alx, aty, aby, afz, abz ;
   size_t center_pos_in_matrix = 0 ;
   bool center_unknown_on_BC = false ;
   
   for( size_t component=0; component<NB_COMPS; ++component )
   {
     // Get local min and max indices
     for (size_t l=0;l<DIM;++l)
     { 
       min_unknown_index(l) =
           get_min_index_unknown_handled_by_proc( component, l );
       max_unknown_index(l) = 
           get_max_index_unknown_handled_by_proc( component, l );
     }

     // Perform assembling
     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {          
       xC = get_DOF_coordinate( i, component, 0 ) ;
       xR = get_DOF_coordinate_Assembling( i+1, component, 0 ) ;
       xL = get_DOF_coordinate_Assembling( i-1, component, 0 ) ;
       dxr = xR - xC ; 
       dxl = xC - xL ; 
       dxC = get_cell_size( i, component, 0 ) ;      
       for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
         yC = get_DOF_coordinate( j, component, 1 ) ; 
 	 yT = get_DOF_coordinate_Assembling( j+1, component, 1 ) ;
 	 yB = get_DOF_coordinate_Assembling( j-1, component, 1 ) ;
	 dyt = yT - yC ;
	 dyb = yC - yB ;
	 dyC = get_cell_size( j, component, 1 ) ; 
 
         if ( DIM == 2 )
	 {
	   size_t k = 0 ;
	   center_pos_in_matrix = DOF_global_number( i, j, k, component );
	   
	   if ( center_pos_in_matrix || !rescale )
	   {
	     center_unknown_on_BC = DOF_on_BC( i, j, k, component ) ;
	   
	     // Right (X)
	     arx = coef_lap * dyC / dxr ;
	     arx = one_DOF_laplacian( MAT, VEC_rhs, rescale, 
	   		i+1, j, k, component, 
	   		center_unknown_on_BC,
	   		center_pos_in_matrix, arx ) ;
	 
	     // Left (X)
	     alx = coef_lap * dyC / dxl ;
	     alx = one_DOF_laplacian( MAT, VEC_rhs, rescale, 
	   		i-1, j, k, component, 
	   		center_unknown_on_BC, 
	   		center_pos_in_matrix, alx ) ;
	   
	     // Top (Y)
	     aty = coef_lap * dxC / dyt ;
	     aty = one_DOF_laplacian( MAT, VEC_rhs, rescale, 
	   		i, j+1, k, component, 
	   		center_unknown_on_BC, 
	   		center_pos_in_matrix, aty ) ;
	   
	     // Bottom (Y)
	     aby = coef_lap * dxC / dyb ;
	     aby = one_DOF_laplacian( MAT, VEC_rhs, rescale, 
	   		i, j-1, k, component, 
	   		center_unknown_on_BC, 
	   		center_pos_in_matrix, aby ) ;
	   
	     // Center
	     ac = - arx - alx - aty - aby ;
	   }
	   else ac = 1. ;
	   
	   MAT->add_to_item( center_pos_in_matrix, center_pos_in_matrix, ac );
	 }
	 else // if dim==3
	 {
	   for( size_t k=min_unknown_index(2); k<=max_unknown_index(2); ++k )
	   {
             zC = get_DOF_coordinate( k, component, 2 ) ; 
 	     zF = get_DOF_coordinate_Assembling( k+1, component, 2 ) ;
 	     zB = get_DOF_coordinate_Assembling( k-1, component, 2 ) ;
	     dzf = zF - zC ;
	     dzb = zC - zB ;
	     dzC = get_cell_size( k, component, 2 ) ;
	     
	     center_pos_in_matrix = DOF_global_number( i, j, k, component );

	     if ( center_pos_in_matrix || !rescale )
	     {
	       center_unknown_on_BC = DOF_on_BC( i, j, k, component ) ;
	     
	       // Right (X)
	       arx = coef_lap * dyC * dzC / dxr ;
	       arx = one_DOF_laplacian( MAT, VEC_rhs, rescale, 
	     		i+1, j, k, component, 
	     		center_unknown_on_BC, 
	     		center_pos_in_matrix, arx ) ;
	 
	       // Left (X)
	       alx = coef_lap * dyC * dzC / dxl ;
	       alx = one_DOF_laplacian( MAT, VEC_rhs, rescale, 
	     		i-1, j, k, component, 
	     		center_unknown_on_BC, 
	     		center_pos_in_matrix, alx ) ;
	   
	       // Top (Y)
	       aty = coef_lap * dxC * dzC / dyt ;
	       aty = one_DOF_laplacian( MAT, VEC_rhs, rescale, 
	     		i, j+1, k, component, 
	     		center_unknown_on_BC, 
	     		center_pos_in_matrix, aty ) ;
	   
	       // Bottom (B)
	       aby = coef_lap * dxC * dzC / dyb ;
	       aby = one_DOF_laplacian( MAT, VEC_rhs, rescale, 
	     		i, j-1, k, component, 
	     		center_unknown_on_BC, 
	     		center_pos_in_matrix, aby ) ;
	   
	       // Front (Z)
	       afz = coef_lap * dxC * dyC / dzf ;
	       afz = one_DOF_laplacian( MAT, VEC_rhs, rescale, 
	     		i, j, k+1, component, 
	     		center_unknown_on_BC, 
	     		center_pos_in_matrix, afz ) ;
				     
	       // Behind (Z)
	       abz = coef_lap * dxC * dyC / dzb ;
	       abz = one_DOF_laplacian( MAT, VEC_rhs, rescale, 
	     		i, j, k-1, component, 
	     		center_unknown_on_BC, 
	     		center_pos_in_matrix, abz ) ;
				     
	       // Center
	       ac = - arx - alx - aty - aby - afz - abz ;
	     }
	     else ac = 1. ;
	     	       
	     MAT->add_to_item( center_pos_in_matrix, center_pos_in_matrix, ac );
	   }
	 }
       } 
     }      	       	   
   }  

   // Synchronize matrix and vector for parallel usage
   MAT->synchronize() ;
   VEC_rhs->synchronize() ;
   
}




//---------------------------------------------------------------------------
double
FV_DiscreteField:: one_DOF_laplacian( LA_Matrix *MAT, LA_Vector *VEC_rhs, 
	bool const& rescale,
	int i, int j, int k,
	size_t component, bool center_unknown_on_BC,
	size_t center_pos_in_matrix, double ai ) const
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_DiscreteField:: one_DOF_laplacian" ) ;

   bool assemble = true ;

   // When the centered unknown is located on a boundary, it implies that the BC
   // is of the Neumann type, hence if the DOF is outside the domain, this means
   // that the flux in the direction (centered DOF -> this DOF) is zero
   if ( center_unknown_on_BC )
     if ( !DOF_in_domain( i, j, k, component ) ) 
     {
       ai = 0. ;
       assemble = false ;
     }
   
   if ( assemble )  
   {
     if ( DOF_is_unknown( i, j, k, component ) )
     {
       size_t pos_in_matrix = DOF_global_number( i, j, k, component );  
       
       if ( pos_in_matrix || !rescale )
         MAT->add_to_item( center_pos_in_matrix, pos_in_matrix, ai );
     }
     else if ( DOF_has_imposed_Dirichlet_value( i, j, k, component ) )
     {
       double dirichlet_value = DOF_value( i, j, k, component, 0 ) ;
       VEC_rhs->add_to_item( center_pos_in_matrix, - ai * dirichlet_value );
     }
     else ai = 0. ;
   }
	 
   return ( ai ) ;
   
}  




//---------------------------------------------------------------------------
void
FV_DiscreteField:: assemble_mass_matrix( double const& coef, 
	LA_Matrix *MAT, double const& power_index ) const
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_DiscreteField:: assemble_mass_matrix" ) ;
   
   // Parameters
   double dxC, dyC, dzC;
   size_t center_pos_in_matrix = 0;
   size_t_vector min_unknown_index(DIM,0);
   size_t_vector max_unknown_index(DIM,0);
   
   for( size_t component=0; component<NB_COMPS; ++component )
   {
     // Get local min and max indices
     for( size_t l=0; l<DIM; ++l )
     {
       min_unknown_index(l) = 
          get_min_index_unknown_handled_by_proc( component, l );
       max_unknown_index(l) = 
          get_max_index_unknown_handled_by_proc( component, l );
     }

     // Perform assembling
     for( size_t i=min_unknown_index(0); i<=max_unknown_index(0); ++i )
     {
       dxC = get_cell_size( i, component, 0 );
       for( size_t j=min_unknown_index(1); j<=max_unknown_index(1); ++j )
       {
         dyC = get_cell_size( j, component, 1 ); 

         if ( DIM == 2 )
         {
           size_t k = 0 ;
           center_pos_in_matrix = DOF_global_number( i, j, k, component );
           MAT->set_item( center_pos_in_matrix, center_pos_in_matrix, 
                  coef * pow( dxC * dyC, power_index ) );
         }
         else
         {
           for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
           {
             dzC = get_cell_size( k, component, 2 );
             center_pos_in_matrix = DOF_global_number( i, j, k, component );
             MAT->set_item( center_pos_in_matrix, center_pos_in_matrix, 
                    coef * pow( dxC * dyC * dzC, power_index ) );
           }
         }
       }
     }
   }

   // Synchronize matrix for parallel usage
   MAT->synchronize() ; 
       
}




//---------------------------------------------------------------------------
void
FV_DiscreteField:: assemble_mass_vector( double const& coef, 
	LA_Vector *VEC, double const& power_index ) const
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_DiscreteField:: assemble_mass_vector" ) ;
   
   // Parameters
   double dxC, dyC, dzC;
   size_t center_pos_in_matrix = 0;
   size_t_vector min_unknown_index(DIM,0);
   size_t_vector max_unknown_index(DIM,0);

   for( size_t component=0; component<NB_COMPS; ++component )
   {
     // Get local min and max indices
     for( size_t l=0; l<DIM; ++l )
     {
       min_unknown_index(l) = 
          get_min_index_unknown_handled_by_proc( component, l );
       max_unknown_index(l) = 
          get_max_index_unknown_handled_by_proc( component, l );
     }

     // Perform assembling
     for( size_t i=min_unknown_index(0); i<=max_unknown_index(0); ++i )
     {
       dxC = get_cell_size( i, component, 0 );
       for( size_t j=min_unknown_index(1); j<=max_unknown_index(1); ++j )
       {
         dyC = get_cell_size( j, component, 1 ); 

         if ( DIM == 2 )
         {
           size_t k = 0 ;
           center_pos_in_matrix = DOF_global_number( i, j, k, component );
           VEC->set_item( center_pos_in_matrix, 
	   	coef * pow( dxC * dyC, power_index ) );
         }
         else
         {
           for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
           {
             dzC = get_cell_size( k, component, 2 );
             center_pos_in_matrix = DOF_global_number( i, j, k, component );
             VEC->set_item( center_pos_in_matrix, 
	     	coef * pow( dxC * dyC * dzC, power_index ) );
           }
         }
       }
     }
   }

   // Synchronize vector for parallel usage
   VEC->synchronize() ; 
       
}






   
//----------------------------------------------------------------------
void
FV_DiscreteField:: assemble_pDivv_matrix( FV_DiscreteField const* UU,
	double const& coef, LA_Matrix *MAT, LA_Vector *VEC_rhs ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: assemble_pDivv_matrix" );   

   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the "
   	<< "FV_DiscreteField_Centered type; "
	<< "method \"assemble_pDivv_matrix\" is "
	<< "implemented for a field of type "
	<< "FV_DiscreteField_Centered only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;     

}




//---------------------------------------------------------------------------
double 
FV_DiscreteField::SuperBee_phi( double const& theta )
//--------------------------------------------------------------------------- 
{
   MAC_LABEL( "FV_DiscreteField:: SuperBee_phi" ) ;

   double min1,min2,max1,max2;

   if (2. * theta < 1.) min1 = 2. * theta;
   else min1 = 1.;
   if (theta < 2.) min2 = theta;
   else min2 = 2.;
   
   if (min1 > min2) max1 = min1;
   else max1 = min2;
   if (max1 > 0.) max2 = max1;
   else max2 = 0.;

   return max2;
   
}




//----------------------------------------------------------------------
double
FV_DiscreteField:: compute_CFL( FV_TimeIterator const* t_it,
	size_t level ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: compute_CFL" );   

   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the "
   	<< "FV_DiscreteField_Staggered type; "
	<< "method \"compute_CFL\" is "
	<< "implemented for a field of type "
	<< "FV_DiscreteField_Staggered only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;     

}




//----------------------------------------------------------------------
void 
FV_DiscreteField:: assemble_advection_Upwind( 
	FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t advected_level,
	LA_Vector *VEC_rhs ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: assemble_advection_Upwind" );   

   ostringstream mesg ;
   mesg << "Field " << FNAME << " is of the "
   	<< FDISCRETIZATION << " type; "
	<< "method \"assemble_advection_Upwind\" is not "
	<< "implemented for fields of type "
	<< FDISCRETIZATION << " !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;     

}	




//----------------------------------------------------------------------
void 
FV_DiscreteField:: assemble_advection_TVD( 
	FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t advected_level,
	LA_Vector *VEC_rhs ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: assemble_advection_TVD" );   

   ostringstream mesg ;
   mesg << "Field " << FNAME << " is of the "
   	<< FDISCRETIZATION << " type; "
	<< "method \"assemble_advection_TVD\" is not "
	<< "implemented for fields of type "
	<< FDISCRETIZATION << " !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;     

}	




//----------------------------------------------------------------------
void
FV_DiscreteField:: assemble_tauGradv_tensor_divergence_matrix( 
      	FV_DiscreteField const* DD,
	double const& coef, LA_Matrix *MAT ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField:: assemble_tauGradv_tensor_divergence_matrix" );

   ostringstream mesg ;
   mesg << "Field " << FNAME << " is not of the "
   	<< "FV_DiscreteField_Staggered type; "
	<< "method \"assemble_tauGradv_tensor_divergence_matrix\" is "
	<< "implemented for a field of type "
	<< "FV_DiscreteField_Staggered only !!" << endl;
   MAC_Error::object()->raise_plain( mesg.str() ) ;     

}
