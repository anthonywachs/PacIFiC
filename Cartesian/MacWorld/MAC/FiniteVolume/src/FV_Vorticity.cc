#include <FV_Vorticity.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DiscreteField.hh>
#include <FV_DiscreteField_Centered.hh>
#include <FV_SystemNumbering.hh>
#include <FV_Mesh.hh>
#include <FV.hh>
#include <FV_TimeIterator.hh>
#include <MAC.hh>
#include <MAC_Root.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Vector.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <LA_Vector.hh>
#include <LA_Scatter.hh>
#include <LA_SeqVector.hh>
#include <LA_DistVector.hh>
#include <LA_Matrix.hh>
#include <size_t_vector.hh>
#include <math.h>

using std::cout ; 
using std::endl ;
using std::ostringstream ; 

FV_Vorticity const* FV_Vorticity::PROTOTYPE = 
	new FV_Vorticity() ;

struct FV_Vorticity_ERROR
{
   static void n1( std::string const& name ) ;   
} ;


//---------------------------------------------------------------------------
FV_Vorticity:: FV_Vorticity( void )
//--------------------------------------------------------------------------
   : FV_OneStepIteration( "FV_Vorticity" ) 
{
   MAC_LABEL( "FV_Vorticity:: FV_Vorticity" ) ;
}




//---------------------------------------------------------------------------
FV_Vorticity*
FV_Vorticity:: create_replica( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_Vorticity:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FV_Vorticity* result = 
                        new FV_Vorticity( a_owner, dom, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;   
}




//---------------------------------------------------------------------------
FV_Vorticity:: FV_Vorticity( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp ) 
   , b_at_each_TS( false )   
   , UU( 0 ) 
   , OM( 0 )
   , velocity_level( 0 ) 
   , velocity_name( "velocity" )
   , OM_NUM( 0 )
   , VEC_OM( 0 )
   , OM_LOC( 0 )
   , MAT_OM( 0 )
{
   MAC_LABEL( "FV_Vorticity:: FV_Vorticity" ) ;

   // MPI
   macCOMM = MAC_Exec::communicator();	   		
   my_rank = macCOMM->rank();
   nb_procs = macCOMM->nb_ranks();
   is_master = 0; 

   // Velocity field
   if ( exp->has_entry( "velocity" ) )
     velocity_name = exp->string_data( "velocity" ) ;         

   if ( !dom->has_discrete_field( velocity_name ) )
     FV_Vorticity_ERROR::n1( velocity_name ) ;

   UU = dom->discrete_field( velocity_name ) ;
   
   // Vorticity field        
   OM = FV_DiscreteField::create( this, dom->primary_grid(),
	"vorticity", "vorticity", UU->nb_components() == 2 ? 1 : 3, 1 ) ;
   OM->build_BCs( NULL, NULL ) ;
   OM->build_field_numbering() ;	 	 	 	 
   OM->out_endOfBuilding( FV::out(), 6, my_rank ) ;
   FV_DomainAndFields* ddom = const_cast<FV_DomainAndFields*>(dom);
   ddom->append_field( OM ) ;
   ddom->append_post_processing_field( OM, "at_vertices", "Vorticity" ) ;

   // Global vector for storage
   size_t om_glob = OM->nb_global_unknowns() ;
   size_t om_loc = OM->nb_local_unknowns() ; 
   
   // For periodic flows, only the PETSC scatter works, while the LA scatter
   // crashes. This is related to the fact that LA scatter does not allow to
   // have twice the same unknown global number on the same process. 
   // This is the case for periodic flows in FV Fiscrete Field. 
   // Conversely, the PETSC scatter does not care. 
   // The only way to safely create a PETSC implementation of VEC_OM is to use
   // of PETSC implementation of a dummy matrix first. It seems rather tricky to
   // create a PETSC matrix or vector on the fly, i.e. without using a 
   // MAC_ModuleExplorer. Hence, when the flow is periodic, the following 
   // module is required in the data file:
   //   MODULE Implementation
   //      concrete_name = "PETSc_MPIAIJ"
   //   END MODULE Implementation  
   // If this module is not specified, the code stops with an error message.
   if ( exp->has_module ( "Implementation" ) )
   {   
     MAC_ModuleExplorer* ee = exp->create_subexplorer(	0, "Implementation" ) ;
     string imp_type = ee->string_data( "concrete_name" ) ;
     if ( imp_type != "PETSc_MPIAIJ" )
     {
       string error_message="   - PETSc_MPIAIJ\n";
       MAC_Error::object()->raise_bad_data_value( exp, 
	"concrete_name", error_message );
     }       
     MAT_OM = LA_Matrix::make( this, ee ) ; 
     VEC_OM = MAT_OM->create_vector( this ) ;
     VEC_OM->re_initialize( om_glob ) ;    
     ee->destroy() ; ee=0 ;
   }
   else 
   {
     if ( UU->primary_grid()->is_periodic_domain() )
       MAC_Error::object()->raise_missing_module( exp, "MODULE Implementation"
	"\n       concrete_name = \"PETSc_MPIAIJ\""
	"\n    END MODULE Implementation" ) ;
     else  
       VEC_OM = LA_DistVector::create( this, om_glob, om_loc, 
     		LA::FromGlobalSize );
   }
   OM_LOC = LA_SeqVector::create( this, 0 ) ;
   OM_LOC->re_initialize( om_loc ) ;   
   OM_NUM = FV_SystemNumbering::create( this, OM ) ;
   OM_NUM->define_scatter( VEC_OM ) ;
   
   if (exp->has_entry( "compute_at_each_TS" ))
   {
     b_at_each_TS=exp->bool_data( "compute_at_each_TS" );
     if( my_rank == is_master ) cout << "b_at_each_TS "<<b_at_each_TS<<endl;
   } 
       
}




//---------------------------------------------------------------------------
FV_Vorticity:: ~FV_Vorticity( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_Vorticity:: ~FV_Vorticity" ) ;

}




//---------------------------------------------------------------------------
void
FV_Vorticity:: do_one_inner_iteration( FV_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_Vorticity:: do_one_inner_iteration" ) ;
   MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
   
   if( b_at_each_TS )
     prepare_vorticity_computation();
}




//---------------------------------------------------------------------------
void
FV_Vorticity:: do_before_time_stepping( FV_TimeIterator const* t_it, 
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_Vorticity:: do_before_time_stepping" ) ;
   
   start_total_timer( "FV_Vorticity:: do_before_time_stepping" ) ;

   stop_total_timer() ;  
}




//---------------------------------------------------------------------------
void
FV_Vorticity:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_Vorticity:: do_after_time_stepping" ) ;  
}




//---------------------------------------------------------------------------
void
FV_Vorticity:: do_before_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_Vorticity:: do_before_inner_iterations_stage" ) ;  
}




//---------------------------------------------------------------------------
void
FV_Vorticity:: do_after_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_Vorticity:: do_after_inner_iterations_stage" ) ;  
}
  



//---------------------------------------------------------------------------
void
FV_Vorticity:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_Vorticity:: do_additional_savings" ) ;

   start_total_timer( "FV_Vorticity:: do_additional_savings" ) ;
   
   if( !b_at_each_TS )
     prepare_vorticity_computation();

   stop_total_timer() ;    
}




//---------------------------------------------------------------------------
void
FV_Vorticity:: prepare_vorticity_computation( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_Vorticity:: prepare_vorticity_computation" ) ;

   if ( my_rank == is_master ) 
     FV::out() << "         Compute vorticity" << endl;

   // Set global vector to zero
   VEC_OM->nullify() ;
   
   // Check for mesh translation
   if ( OM->primary_grid()->is_translation_active() )
     OM->check_field_primary_meshes_coincide( true, FV::out(), 9 ) ;

   // Compute vorticity and store it in global vector
   if ( UU->nb_components() == 2 ) compute_vorticity_2D() ; 
   else compute_vorticity_3D() ; 
   VEC_OM->synchronize() ;
   
   // Copy back to local vector
   OM_NUM->scatter()->get( VEC_OM, OM_LOC ) ; 
   OM->update_free_DOFs_value( 0, OM_LOC ) ;
}




//---------------------------------------------------------------------------
void
FV_Vorticity:: compute_vorticity_2D( void )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( 
   "FV_Vorticity:: compute_vorticity_2D" ) ;

   // Parameters
   size_t dim = 2 ;
   double xR, xL, yT, yB, vorticity = 0. ;
   size_t pos_in_vector = 0 ;
   FV_SHIFT_TRIPLET shift = OM->shift_vorticityToStaggered( 0 ) ;
   
   // Get local min and max indices   
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) min_unknown_index(l) = 
       	OM->get_min_index_unknown_handled_by_proc( 0, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) max_unknown_index(l) = 
       	OM->get_max_index_unknown_handled_by_proc( 0, l ) ;
	
   // Compute vorticity and store into global vector
   for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
   {          
     xR = UU->get_DOF_coordinate_Assembling( i+shift.i, 1, 0 ) ;
     xL = UU->get_DOF_coordinate_Assembling( i+shift.i-1, 1, 0 ) ;     
     for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
     {
       pos_in_vector = OM->DOF_global_number( i, j, 0, 0 );
       yT = UU->get_DOF_coordinate_Assembling( j+shift.j, 0, 1 ) ;
       yB = UU->get_DOF_coordinate_Assembling( j+shift.j-1, 0, 1 ) ;
       
       vorticity = 
       	( UU->DOF_value( i+shift.i, j, 0, 1, velocity_level ) 
       	- UU->DOF_value( i+shift.i-1, j, 0, 1, velocity_level ) ) 
	/ ( xR - xL )
	- ( UU->DOF_value( i, j+shift.j, 0, 0, velocity_level ) 
       	- UU->DOF_value( i, j+shift.j-1, 0, 0, velocity_level ) ) 
	/ ( yT - yB ) ;
        
       VEC_OM->set_item( pos_in_vector, vorticity ) ;
     }      	       	   
   } 
}




//---------------------------------------------------------------------------
void
FV_Vorticity:: compute_vorticity_3D( void )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_Vorticity:: compute_vorticity_3D" ) ;

   // Parameters
   size_t dim = 3 ;
   double xR, xL, yT, yB, zF, zB, vorticity = 0. ;
   size_t pos_in_vector = 0 ;
   
   // Compute vorticity and store into global vector   
   for (size_t comp=0;comp<3;++comp)
   {
     // Get local min and max indices
     size_t_vector min_unknown_index(dim,0);
     for (size_t l=0;l<dim;++l) 
       min_unknown_index(l) = 
       	OM->get_min_index_unknown_handled_by_proc( comp, l ) ;
     size_t_vector max_unknown_index(dim,0);
     for (size_t l=0;l<dim;++l) 
       max_unknown_index(l) = 
       	OM->get_max_index_unknown_handled_by_proc( comp, l ) ;

     FV_SHIFT_TRIPLET shift = OM->shift_vorticityToStaggered( comp ) ;
	
     // Perform computing & storing
     switch(comp)
     {
       case 0:
         for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
           for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
           {
             yT = UU->get_DOF_coordinate_Assembling( j+shift.j, comp, 1 ) ;
             yB = UU->get_DOF_coordinate_Assembling( j+shift.j-1, comp, 1 ) ;
	     for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	     {
 	       zF = UU->get_DOF_coordinate_Assembling( k+shift.k, comp, 2 ) ;
 	       zB = UU->get_DOF_coordinate_Assembling( k+shift.k-1, comp, 2 ) ;
	       
	       vorticity = ( 
	       	UU->DOF_value( i, j+shift.j, k, 2, velocity_level )
       		- UU->DOF_value( i, j+shift.j-1, k, 2, velocity_level ) ) 
		/ ( yT - yB )
		- ( UU->DOF_value( i, j, k+shift.k, 1, velocity_level ) 
       		- UU->DOF_value( i, j, k+shift.k-1, 1, velocity_level ) ) 
		/ ( zF - zB ) ;
               
	       pos_in_vector = OM->DOF_global_number( i, j, k, comp );
	       VEC_OM->set_item( pos_in_vector, vorticity ) ;
	     }
	   }
         break ;	 
	       
       case 1:
         for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
         {          
           xR = UU->get_DOF_coordinate_Assembling( i+shift.i, comp, 0 ) ;
           xL = UU->get_DOF_coordinate_Assembling( i+shift.i-1, comp, 0 ) ;    
           for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
	     for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	     {
 	       zF = UU->get_DOF_coordinate_Assembling( k+shift.k, comp, 2 ) ;
 	       zB = UU->get_DOF_coordinate_Assembling( k+shift.k-1, comp, 2 ) ;
	       
	       vorticity = ( 
	       	UU->DOF_value( i, j, k+shift.k, 0, velocity_level )
       		- UU->DOF_value( i, j, k+shift.k-1, 0, velocity_level ) ) 
		/ ( zF - zB )
		- ( UU->DOF_value( i+shift.i, j, k, 2, velocity_level ) 
       		- UU->DOF_value( i+shift.i-1, j, k, 2, velocity_level ) ) 
		/ ( xR - xL ) ;
               
	       pos_in_vector = OM->DOF_global_number( i, j, k, comp );
	       VEC_OM->set_item( pos_in_vector, vorticity ) ;
	     }
	 }
         break ;
	       
       case 2:
         for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
         {          
           xR = UU->get_DOF_coordinate_Assembling( i+shift.i, comp, 0 ) ;
           xL = UU->get_DOF_coordinate_Assembling( i+shift.i-1, comp, 0 ) ;     
           for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
           {
             yT = UU->get_DOF_coordinate_Assembling( j+shift.j, comp, 1 ) ;
             yB = UU->get_DOF_coordinate_Assembling( j+shift.j-1, comp, 1 ) ;
	     for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	     {	       
	       vorticity = ( 
	       	UU->DOF_value( i+shift.i, j, k, 1, velocity_level )
       		- UU->DOF_value( i+shift.i-1, j, k, 1, velocity_level ) ) 
		/ ( xR - xL )
		- ( UU->DOF_value( i, j+shift.j, k, 0, velocity_level ) 
       		- UU->DOF_value( i, j+shift.j-1, k, 0, velocity_level ) ) 
		/ ( yT - yB ) ;
               
	       pos_in_vector = OM->DOF_global_number( i, j, k, comp );
	       VEC_OM->set_item( pos_in_vector, vorticity ) ;
	     }
	   }
	 }
         break ;	       	       
     }     
   }      
}




//----------------------------------------------------------------------------
void
FV_Vorticity:: do_more_post_processing(
        FV_DomainAndFields* dom,
        MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_Vorticity:: do_more_post_processing" ) ;

   if ( my_rank == is_master ) 
     FV::out() << "         Nothing to do in FV_Vorticity::"
     	"do_more_post_processing " << endl << endl;
}




//internal--------------------------------------------------------------
void
FV_Vorticity_ERROR:: n1( std::string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Field \"" << name << "\" is missing in MODULE interior_fields "
   	<< "to build FV_Vorticity" << endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}
