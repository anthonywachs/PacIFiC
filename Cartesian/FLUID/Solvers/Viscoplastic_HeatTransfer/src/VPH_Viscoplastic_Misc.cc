#include <VPH_Viscoplastic.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DomainBuilder.hh>
#include <FV_DiscreteField.hh>
#include <LA_Matrix.hh>
#include <LA_Vector.hh>
#include <LA_SeqVector.hh>
#include <MAC_Communicator.hh>
#include <MAC_DoubleVector.hh>
#include <math.h>


//---------------------------------------------------------------------------
void
VPH_Viscoplastic:: compute_strain_rate_tensor_D ( LA_Vector* VEC )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "VPH_Viscoplastic:: compute_strain_rate_tensor_D" ) ;

   // Parameters
   size_t_vector min_unknown_index( dim, 0 );
   size_t_vector max_unknown_index( dim, 0 );
   size_t nb_comps = DD->nb_components(), comp ;
   double value = 0. ;
   
   FV_SHIFT_TRIPLET shift ;

   for (comp=0;comp<nb_comps;++comp)
   {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) 
       min_unknown_index(l) = 
       	DD->get_min_index_unknown_handled_by_proc( comp, l ) ;
     for (size_t l=0;l<dim;++l) 
       max_unknown_index(l) = 
       	DD->get_max_index_unknown_handled_by_proc( comp, l ) ;

     shift = DD->shift_staggeredToTensor( comp ) ;
          
     // Perform assembling
     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {          
       for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
         if ( dim == 2 )
	 {
	   size_t k = 0 ;
	   
	   // The First Component (Dxx)
	   if ( comp == 0 )
	   {
	     value = ( UU->DOF_value( i+shift.i, j, k, comp, 0 )
	     		- UU->DOF_value( i+shift.i-1, j, k, comp, 0 ) )
			/ DD->get_cell_size( i, comp, 0 ) ;
	   }
	   
	   // The Second Component (Dyy)
	   else if ( comp == 1 )
	   {
	     value = ( UU->DOF_value( i, j+shift.j, k, comp, 0 )
	     		- UU->DOF_value( i, j+shift.j-1, k, comp, 0 ) )
			/ DD->get_cell_size( j, comp, 1 ) ;
	   }
	   
	   // The Third Component (Dxy)
	   else
	   {
	       value = 0.5 * ( ( UU->DOF_value( i, j+shift.j, k, 0, 0 )
	     		- UU->DOF_value( i, j+shift.j-1, k, 0, 0 ) )
			/ DD->get_cell_size( j, comp, 1 )
		+ ( UU->DOF_value( i+shift.i, j, k, 1, 0 )
	     		- UU->DOF_value( i+shift.i-1, j, k, 1, 0 ) )
			/ DD->get_cell_size( i, comp, 0 ) ) ;
	   }
   
	   // Set the value in the distributed vector
	   VEC->set_item( DD->DOF_global_number( i, j, k, comp ), value ) ;	
	 }
	 else
	 {
	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	   {
	     // The First Component (Dxx)
	     if ( comp == 0 )
	     {
	       value = ( UU->DOF_value( i+shift.i, j, k, comp, 0 )
	     		- UU->DOF_value( i+shift.i-1, j, k, comp, 0 ) )
			/ DD->get_cell_size( i, comp, 0 ) ;
	     }
	     
	     // The Second Component (Dyy)
	     else if ( comp == 1 )
	     {
	       value = ( UU->DOF_value( i, j+shift.j, k, comp, 0 )
	     		- UU->DOF_value( i, j+shift.j-1, k, comp, 0 ) )
			/ DD->get_cell_size( j, comp, 1 ) ;
	     }
	     
	     // The Third Component (Dzz)
	     else if ( comp == 2 )
	     {
	       value = ( UU->DOF_value( i, j, k+shift.k, comp, 0 )
	     		- UU->DOF_value( i, j, k+shift.k-1, comp, 0 ) )
			/ DD->get_cell_size( k, comp, 2 ) ;
	     }
	     
	     // The Fourth Component (Dxy)
	     else if ( comp == 3 )
	     {
	       value = 0.5 * ( ( UU->DOF_value( i, j+shift.j, k, 0, 0 )
	     		- UU->DOF_value( i, j+shift.j-1, k, 0, 0 ) )
			/ DD->get_cell_size( j, comp, 1 )
		+ ( UU->DOF_value( i+shift.i, j, k, 1, 0 )
	     		- UU->DOF_value( i+shift.i-1, j, k, 1, 0 ) )
			/ DD->get_cell_size( i, comp, 0 ) ) ;
	     }
	     
	     // The Fifth Component (Dxz)
	     else if ( comp == 4 )
	     {
	       value = 0.5 * ( ( UU->DOF_value( i, j, k+shift.k, 0, 0 )
	     		- UU->DOF_value( i, j, k+shift.k-1, 0, 0 ) )
			/ DD->get_cell_size( k, comp, 2 )
		+ ( UU->DOF_value( i+shift.i, j, k, 2, 0 )
	     		- UU->DOF_value( i+shift.i-1, j, k, 2, 0 ) )
			/ DD->get_cell_size( i, comp, 0 ) ) ;
	     }
	     
	     // The sixth Component (Dyz)
	     else
	     {
	       value = 0.5 * ( ( UU->DOF_value( i, j, k+shift.k, 1, 0 )
	     		- UU->DOF_value( i, j, k+shift.k-1, 1, 0 ) )
			/ DD->get_cell_size( k, comp, 2 )
		+ ( UU->DOF_value( i, j+shift.j, k, 2, 0 )
	     		- UU->DOF_value( i, j+shift.j-1, k, 2, 0 ) )
			/ DD->get_cell_size( j, comp, 1 ) ) ;
	     }

	     // Set the value in the distributed vector
	     VEC->set_item( DD->DOF_global_number( i, j, k, comp ), value ) ;
	   }
	 }
       }
     }
   }

   // Synchronize vector for parallel usage
   VEC->synchronize() ;
   
   // Copy back in field values
   DD->update_free_DOFs_value( D_level, 
   	GLOBAL_EQ->get_solution_strainrate_tensor() ) ;     

}




//---------------------------------------------------------------------------
void
VPH_Viscoplastic:: update_strain_rate_tensor_d ( LA_Vector* VEC, 
	double const& alpha, double const& beta, 
	size_t const& first_field_level )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "VPH_Viscoplastic:: update_strain_rate_tensor_d" ) ;

   // Parameters
   size_t_vector min_unknown_index( dim, 0 );
   size_t_vector max_unknown_index( dim, 0 );
   size_t nb_comps = DD->nb_components() ;
   double value = 0., norm = 0. ;
   FV_SHIFT_TRIPLET shift ;
   
   for (size_t comp=0;comp<nb_comps;++comp)
   {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) 
       min_unknown_index(l) = 
       	DD->get_min_index_unknown_handled_by_proc( comp, l ) ;
     for (size_t l=0;l<dim;++l) 
       max_unknown_index(l) = 
       	DD->get_max_index_unknown_handled_by_proc( comp, l ) ;

     shift = DD->shift_tensorToTensor( comp ) ;
     
     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {          
       for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
         if ( dim == 2 )
	 {
	   size_t k = 0 ;
	   
           // Compute the Euclidian norm of first_field+beta.D(u)
	   norm = Compute_euclidian_norm( i, j, comp, shift, beta, 
	   	first_field_level ) ;

           // Update each component of d
           if ( norm > yield_stress )
             value = ( 1. - yield_stress / norm ) / alpha
	     		* ( DD->DOF_value( i, j, k, comp, first_field_level )
				+ beta
				* DD->DOF_value( i, j, k, comp, D_level ) ) ;
 	   else
	     value = 0. ;
	   
	   // Set the value in the distributed vector
	   VEC->set_item( DD->DOF_global_number( i, j, k, comp ), value ) ;
	 }
	 else
	 {
	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	   {
             // Compute the Euclidian norm of first_field+beta.D(u)
             norm = Compute_euclidian_norm( i, j, k, comp, shift, beta,
	     	first_field_level ) ;

             // Update each component of d
             if ( norm > yield_stress )
               value = ( 1. - yield_stress / norm ) / alpha
	     		* ( DD->DOF_value( i, j, k, comp, first_field_level )
				+ beta
				* DD->DOF_value( i, j, k, comp, D_level ) ) ;
             else
	       value = 0. ;	   

             // Set the value in the distributed vector
	     VEC->set_item( DD->DOF_global_number( i, j, k, comp ), value ) ;
	   }
	 }
       }
     }
   }

   // Synchronize vector for parallel usage
   VEC->synchronize() ;
   
   // Copy back in field values
   DD->update_free_DOFs_value( d_level, 
   	GLOBAL_EQ->get_solution_strainrate_tensor() ) ;      

}




//---------------------------------------------------------------------------
double
VPH_Viscoplastic:: Compute_euclidian_norm (
		size_t i, size_t j, size_t comp,
		FV_SHIFT_TRIPLET shift, double const& beta, 
		size_t const& first_field_level )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "VPH_Viscoplastic:: Compute_euclidian_norm" ) ;

   size_t nb_comps = DD->nb_components() ;
   doubleVector LambdaPlusLagParamDotDD( nb_comps ); 
   double value_first_field = 0., value_Du = 0.;
   double norm = 0.;
   
   for (size_t localComp=0;localComp<nb_comps;++localComp)
   {		
     value_first_field = DD->interpolateOneCompOnAnotherComp(
		i, j, localComp, comp, first_field_level, shift ) ;
     value_Du = DD->interpolateOneCompOnAnotherComp(
		i, j, localComp, comp, D_level, shift ) ;
     LambdaPlusLagParamDotDD( localComp )
 	   	= value_first_field + beta * value_Du ;     		
   }

   norm = 0.5 * ( pow(LambdaPlusLagParamDotDD(0),2.)
     		+ pow(LambdaPlusLagParamDotDD(1),2.) )
		+ pow(LambdaPlusLagParamDotDD(2),2.) ;
   norm = sqrt( norm );

   return( norm );   

}




//---------------------------------------------------------------------------
double
VPH_Viscoplastic:: Compute_euclidian_norm (
		size_t i, size_t j, size_t k, size_t comp,
		FV_SHIFT_TRIPLET shift, double const& beta, 
		size_t const& first_field_level )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "VPH_Viscoplastic:: Compute_euclidian_norm" ) ;

   size_t nb_comps = DD->nb_components();
   doubleVector LambdaPlusLagParamDotDD(nb_comps);
   double value_first_field = 0., value_Du = 0.;   
   double norm = 0.;
   
   for (size_t localComp=0;localComp<nb_comps;++localComp)
   {
     value_first_field = DD->interpolateOneCompOnAnotherComp(
		i, j, k, localComp, comp, first_field_level, shift ) ;
     value_Du = DD->interpolateOneCompOnAnotherComp(
		i, j, k, localComp, comp, D_level, shift ) ;
     LambdaPlusLagParamDotDD( localComp )
 	   	= value_first_field + beta * value_Du ; 		
   }

   norm = 0.5 * ( pow(LambdaPlusLagParamDotDD(0),2.)
     		+ pow(LambdaPlusLagParamDotDD(1),2.)
		+ pow(LambdaPlusLagParamDotDD(2),2.) )
     		+ pow(LambdaPlusLagParamDotDD(3),2.)
		+ pow(LambdaPlusLagParamDotDD(4),2.)
		+ pow(LambdaPlusLagParamDotDD(5),2.);
   norm = sqrt( norm );

   return( norm );

}




//---------------------------------------------------------------------------
void
VPH_Viscoplastic:: update_Lagrange_multiplier ( LA_Vector* VEC )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_ViscoplasticALG2:: update_Lagrange_multiplier" ) ;

   // Parameters
   size_t_vector min_unknown_index( dim, 0 );
   size_t_vector max_unknown_index( dim, 0 );
   size_t nb_comps = DD->nb_components(), comp ;
   double value = 0. ;
   
   for (comp=0;comp<nb_comps;++comp)
   {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) 
       min_unknown_index(l) = 
       	DD->get_min_index_unknown_handled_by_proc( comp, l ) ;
     for (size_t l=0;l<dim;++l) 
       max_unknown_index(l) = 
       	DD->get_max_index_unknown_handled_by_proc( comp, l ) ;

     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {          
       for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
         if ( dim == 2 )
	 {
	   size_t k = 0 ;
	   value = DD->DOF_value( i, j, k, comp, lambda_level )
	   	+ ALG2_aug_param
			* ( DD->DOF_value( i, j, k, comp, D_level )
				- DD->DOF_value( i, j, k, comp, d_level ) ) ;
           
	   // Set the value in the distributed vector
           VEC->set_item( DD->DOF_global_number( i, j, k, comp ), value ) ;
	 }
	 else
	 {
	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	   {
	     value = DD->DOF_value( i, j, k, comp, lambda_level )
	   	+ ALG2_aug_param
			* ( DD->DOF_value( i, j, k, comp, D_level )
			- DD->DOF_value( i, j, k, comp, d_level ) ) ;

	     // Set the value in the distributed vector
             VEC->set_item( DD->DOF_global_number( i, j, k, comp ), value ) ;
	   }
	 }
       }
     }
   }

   // Synchronize vector for parallel usage
   VEC->synchronize() ;
   
   // Copy back in field values
   DD->update_free_DOFs_value( lambda_level, 
   	GLOBAL_EQ->get_solution_lambda_tensor() ) ;   
	   
}




//---------------------------------------------------------------------------
void
VPH_Viscoplastic:: update_tau ( LA_Vector* VEC )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "VPH_Viscoplastic:: update_tau" ) ;

   // Parameters
   size_t_vector min_unknown_index( dim, 0 );
   size_t_vector max_unknown_index( dim, 0 );
   size_t nb_comps = DD->nb_components(), comp ;
   double value = 0. ;

   for (comp=0;comp<nb_comps;++comp)
   {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) 
       min_unknown_index(l) = 
       	DD->get_min_index_unknown_handled_by_proc( comp, l ) ;
     for (size_t l=0;l<dim;++l) 
       max_unknown_index(l) = 
       	DD->get_max_index_unknown_handled_by_proc( comp, l ) ;

     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {          
       for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
         if ( dim == 2 )
	 {
	   size_t k = 0 ;
	   value = DD->DOF_value( i, j, k, comp, tau_hat_level )
	   	+ FISTA_1overL_param
			* ( DD->DOF_value( i, j, k, comp, D_level )
				- DD->DOF_value( i, j, k, comp, d_level ) ) ;

           // Set the value in the distributed vector
           VEC->set_item( DD->DOF_global_number( i, j, k, comp ), value ) ;
	 }
	 else
	 {
	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	   {
	     value = DD->DOF_value( i, j, k, comp, tau_hat_level )
	   	+ FISTA_1overL_param
			* ( DD->DOF_value( i, j, k, comp, D_level )
			- DD->DOF_value( i, j, k, comp, d_level ) ) ;

             // Set the value in the distributed vector
             VEC->set_item( DD->DOF_global_number( i, j, k, comp ), value ) ;
	   }
	 }
       }
     }
   }

   // Synchronize vector for parallel usage
   VEC->synchronize() ;
   
   // Copy back in field values
   DD->update_free_DOFs_value( tau_level, 
   	GLOBAL_EQ->get_solution_stress_tensor() ) ;   
	   
}




//---------------------------------------------------------------------------
void
VPH_Viscoplastic:: update_tau_hat ( LA_Vector* VEC, 
	double const& stepsize )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "VPH_Viscoplastic:: update_tau_hat" ) ;

   // Parameters
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   size_t nb_comps = DD->nb_components(), comp ;
   double value = 0. ;
   
   for (comp=0;comp<nb_comps;++comp)
   {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) 
       min_unknown_index(l) = 
       	DD->get_min_index_unknown_handled_by_proc( comp, l ) ;
     for (size_t l=0;l<dim;++l) 
       max_unknown_index(l) = 
       	DD->get_max_index_unknown_handled_by_proc( comp, l ) ;

     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {          
       for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
         if ( dim == 2 )
	 {
	   size_t k = 0 ;
	   value = DD->DOF_value( i, j, k, comp, tau_level );
	   value += stepsize * ( value
			- DD->DOF_value( i, j, k, comp, tau_old_level ) ) ;

           // Set the value in the distributed vector
           VEC->set_item( DD->DOF_global_number( i, j, k, comp ), value ) ;
	 }
	 else
	 {
	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	   {
	     value = DD->DOF_value( i, j, k, comp, tau_level );
	     value += stepsize * ( value
			- DD->DOF_value( i, j, k, comp, tau_old_level ) ) ;

             // Set the value in the distributed vector
             VEC->set_item( DD->DOF_global_number( i, j, k, comp ), value ) ;
	   }
	 }
       }
     }
   }

   // Synchronize vector for parallel usage
   VEC->synchronize() ;
   
   // Copy back in field values
   DD->update_free_DOFs_value( tau_hat_level, 
   	GLOBAL_EQ->get_solution_stress_tensor() ) ;  

}



//---------------------------------------------------------------------------
double
VPH_Viscoplastic:: compute_norm_Dminusd ( )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "VPH_Viscoplastic:: compute_norm_Dminusd" ) ;

   size_t_vector min_unknown_index( dim, 0 );
   size_t_vector max_unknown_index( dim, 0 );
   size_t nb_comps = DD->nb_components() ;
   double norm = 0., norm_collective = 0.;

   for (size_t comp=0;comp<nb_comps;++comp)
   {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) 
       min_unknown_index(l) = 
       	DD->get_min_index_unknown_handled_by_proc( comp, l ) ;
     for (size_t l=0;l<dim;++l) 
       max_unknown_index(l) = 
       	DD->get_max_index_unknown_handled_by_proc( comp, l ) ;

     norm = 0. ;

     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {          
       for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
         if ( dim == 2 )
	 {
	   size_t k = 0 ;
	   norm += pow( DD->DOF_value( i, j, k, comp, D_level )
	   		- DD->DOF_value( i, j, k, comp, d_level ), 2.) ;
	 }
	 else
	 {
	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	     norm += pow( DD->DOF_value( i, j, k, comp, D_level )
			- DD->DOF_value( i, j, k, comp, d_level ), 2.) ;
	 }
       }
     }

     norm_collective += sqrt( macCOMM->sum( norm ) );
     
   }

//   norm_collective = pelCOMM->sum( norm );
//   norm_collective += sqrt(norm_collective);

   return( norm_collective );   

}




//---------------------------------------------------------------------------
void
VPH_Viscoplastic:: assemble_buoyancy_rhs ( 
	LA_Vector* VEC_rhs )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "VPH_Viscoplastic:: assemble_buoyancy_rhs" ) ;
   
  // Parameters
  size_t nb_comps = UU->nb_components() ;
  double dxC, dyC, dzC, meanTemp ;
  size_t center_pos_in_matrix = 0 ;
  size_t_vector min_unknown_index( dim, 0 );
  size_t_vector max_unknown_index( dim, 0 );   
  FV_SHIFT_TRIPLET shift ;

  doubleVector gravityTerm( nb_comps ) ;
  doubleVector const& gg = gravity_vector->to_double_vector();  
  for ( size_t i = 0; i < nb_comps; ++i )
    gravityTerm( i ) = gg( i ) * density * boussinesq_thermal_expansion_coef ;

  for (size_t comp=0;comp<nb_comps;++comp)
  {
    // Get local min and max indices
    for (size_t l=0;l<dim;++l) 
      min_unknown_index(l) = 
       	UU->get_min_index_unknown_handled_by_proc( comp, l ) ;
    for (size_t l=0;l<dim;++l) 
      max_unknown_index(l) = 
       	UU->get_max_index_unknown_handled_by_proc( comp, l ) ;
	
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
    shift = UU->shift_staggeredToStaggered( comp ) ;

    // Perform assembling
    for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
    {          
      dxC = UU->get_cell_size( i, comp, 0 ) ;      
      for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
      {
	dyC = UU->get_cell_size( j, comp, 1 ) ; 
 
        if ( dim == 2 )
        {
	  size_t k = 0 ;
	  center_pos_in_matrix = UU->DOF_global_number( i, j, k, comp );
	   
	  // Interpolate temperature 
          if ( comp == 0 )
	  {
	    if ( UU->DOF_color( i, j, k, comp ) == FV_BC_RIGHT )
	      meanTemp = TT->DOF_value( i+shift.i, j, k, 0, 0) ;
	    else if ( UU->DOF_color( i, j, k, comp ) == FV_BC_LEFT )
	      meanTemp = TT->DOF_value( i+shift.i-1, j, k, 0, 0) ;
	    else 
	      meanTemp = 0.5 * ( 
	      		TT->DOF_value( i+shift.i, j, k, 0, 0) +
	       		TT->DOF_value( i+shift.i-1, j, k, 0, 0) ) ;
          }
	  else  // comp 1
	  {
	    if ( UU->DOF_color( i, j, k, comp ) == FV_BC_TOP )
	       meanTemp = TT->DOF_value( i, j+shift.j, k, 0, 0) ;
	     else if ( UU->DOF_color( i, j, k, comp ) == FV_BC_BOTTOM )
	       meanTemp = TT->DOF_value( i, j+shift.j-1, k, 0, 0) ;
	     else
	       meanTemp = 0.5 * ( 
	       		TT->DOF_value( i, j+shift.j, k, 0, 0) +
	       		TT->DOF_value( i, j+shift.j-1, k, 0, 0) ) ;
	  }

	  VEC_rhs->set_item( center_pos_in_matrix, 
		gravityTerm( comp ) 
		* ( meanTemp - boussinesq_reference_temperature ) 
		* dxC * dyC ) ;
        }
	else // dim = 3
	{
	  for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	  {
	    dzC = UU->get_cell_size( k, comp, 2 ) ;	     
	    center_pos_in_matrix = UU->DOF_global_number( i, j, k, comp );
	     
	    // Interpolate temperature 
            if ( comp == 0 )
	    {
	      if ( UU->DOF_color( i, j, k, comp ) == FV_BC_RIGHT )
	        meanTemp = TT->DOF_value( i+shift.i, j, k, 0, 0) ;
	      else if ( UU->DOF_color( i, j, k, comp ) == FV_BC_LEFT )
	        meanTemp = TT->DOF_value( i+shift.i-1, j, k, 0, 0) ;
	      else 
	        meanTemp = 0.5 * ( 
			TT->DOF_value( i+shift.i, j, k, 0, 0) +
	       		TT->DOF_value( i+shift.i-1, j, k, 0, 0) ) ;
	    }
	    
	    else if ( comp == 1 )
	    {
	      if ( UU->DOF_color( i, j, k, comp ) == FV_BC_TOP )
	        meanTemp = TT->DOF_value( i, j+shift.j, k, 0, 0) ;
	      else if ( UU->DOF_color( i, j, k, comp ) == FV_BC_BOTTOM )
	        meanTemp = TT->DOF_value( i, j+shift.j-1, k, 0, 0) ;
	      else
	        meanTemp = 0.5 * ( 
			TT->DOF_value( i, j+shift.j, k, 0, 0) +
	       		TT->DOF_value( i, j+shift.j-1, k, 0, 0) ) ;
	    }
	    
	    else
            {
	      if ( UU->DOF_color( i, j, k, comp ) == FV_BC_FRONT )
	        meanTemp = TT->DOF_value( i, j, k+shift.k, 0, 0) ;
	      else if ( UU->DOF_color( i, j, k, comp ) == FV_BC_BEHIND )
	        meanTemp = TT->DOF_value( i, j, k+shift.k-1, 0, 0) ;
	      else
	        meanTemp = 0.5 * ( 
			TT->DOF_value( i, j, k+shift.k, 0, 0) +
	       		TT->DOF_value( i, j, k+shift.k-1, 0, 0) ) ;
	    }

            VEC_rhs->set_item( center_pos_in_matrix, 
		gravityTerm( comp ) 
		* ( meanTemp - boussinesq_reference_temperature ) 
		* dxC * dyC * dzC ) ;
	  }
	}
      } 
    }      	       	   
  }

  // Synchronize vector for parallel usage
  VEC_rhs->synchronize();   
      
}
