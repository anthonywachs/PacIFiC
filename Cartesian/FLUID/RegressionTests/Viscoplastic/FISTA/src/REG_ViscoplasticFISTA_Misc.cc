#include <REG_ViscoplasticFISTA.hh>
#include <LA_Matrix.hh>
#include <LA_Vector.hh>
#include <LA_SeqVector.hh>
#include <MAC_Communicator.hh>
#include <math.h>


//---------------------------------------------------------------------------
void
REG_ViscoplasticFISTA:: compute_strain_rate_tensor_D ( LA_Vector* VEC )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_ViscoplasticFISTA:: compute_strain_rate_tensor_D" ) ;

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
REG_ViscoplasticFISTA:: update_strain_rate_tensor_d ( LA_Vector* VEC, 
	double const& alpha, double const& beta, 
	size_t const& first_field_level )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_ViscoplasticFISTA:: update_strain_rate_tensor_d" ) ;

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
REG_ViscoplasticFISTA:: Compute_euclidian_norm (
		size_t i, size_t j, size_t comp,
		FV_SHIFT_TRIPLET shift, double const& beta, 
		size_t const& first_field_level )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_ViscoplasticFISTA:: Compute_euclidian_norm" ) ;

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
REG_ViscoplasticFISTA:: Compute_euclidian_norm (
		size_t i, size_t j, size_t k, size_t comp,
		FV_SHIFT_TRIPLET shift, double const& beta, 
		size_t const& first_field_level )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_ViscoplasticFISTA:: Compute_euclidian_norm" ) ;

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
REG_ViscoplasticFISTA:: update_tau ( LA_Vector* VEC )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_ViscoplasticFISTA:: update_tau" ) ;

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
REG_ViscoplasticFISTA:: update_tau_hat ( LA_Vector* VEC, 
	double const& stepsize )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_ViscoplasticFISTA:: update_tau_hat" ) ;

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
REG_ViscoplasticFISTA:: compute_norm_Dminusd ( )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_ViscoplasticFISTA:: compute_norm_Dminusd" ) ;

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
