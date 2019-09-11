#include <REG_HeatTransfer.hh>
#include <MAC.hh>
#include <MAC_Exec.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Timer.hh>
#include <MAC_Vector.hh>
#include <MAC_Communicator.hh>
#include <MAC_ListIdentity.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Variable.hh>
#include <FV_DiscreteField.hh>
#include <FV_TimeIterator.hh>
#include <FV_DomainAndFields.hh>
#include <FV_Mesh.hh>
#include <LA_Vector.hh>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>



//----------------------------------------------------------------------
REG_HeatTransfer*
REG_HeatTransfer:: create( MAC_Object* a_owner,
          MAC_ModuleExplorer const* exp,
          struct NavierStokes2Temperature const& transfert )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatTransfer:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;   

   REG_HeatTransfer* result = 
         new REG_HeatTransfer( a_owner, exp, transfert ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;	  
         
   return( result ) ;
   
}




//----------------------------------------------------------------------
REG_HeatTransfer:: REG_HeatTransfer( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	struct NavierStokes2Temperature const& fromNS )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , TT ( fromNS.dom_->discrete_field( "temperature" ) )
   , UU ( fromNS.UU_ )
   , GLOBAL_EQ( 0 )
   , density( fromNS.density_ )
   , heat_capacity( exp->double_data( "Heat_capacity") )
   , thermal_conductivity( exp->double_data( "Thermal_conductivity") ) 
   , imposed_CFL( fromNS.imposed_CFL_ )
   , AdvectionScheme( "TVD" )
   , resultsDirectory( fromNS.resultsDirectory_ )
   , DiffusionTimeAccuracy( 1 )
   , AdvectionTimeAccuracy( 1 )
   , levelAdvectingVelocity( fromNS.levelAdvectingVelocity_ )
   , b_restart( fromNS.b_restart_ )
   , b_sourceTerm( false )
   , sourceTerm_CTX( 0 )
   , sourceTerm_formula( 0 )
   , sourceTerm_coords( 0 )    
{
   MAC_LABEL( "REG_HeatTransfer:: REG_HeatTransfer" ) ;
   MAC_ASSERT( TT->discretization_type() == "centered" ) ;

   // Call of MAC_Communicator routine to set the rank of each proces and
   // the number of processes during execution of REG_ProjNSWithHeatTransfer
   macCOMM = MAC_Exec::communicator();
   my_rank = macCOMM->rank();
   nb_ranks = macCOMM->nb_ranks();
   is_master = 0;

   
   // Advection scheme
   if ( exp->has_entry( "AdvectionScheme" ) )
     AdvectionScheme = exp->string_data( "AdvectionScheme" );
   if ( AdvectionScheme != "Upwind" && AdvectionScheme != "TVD" )
   {
     string error_message="   - Upwind\n   - TVD";
     MAC_Error::object()->raise_bad_data_value( exp,
        "AdvectionScheme", error_message );
   }
   if ( AdvectionScheme == "TVD"
   	&& UU->primary_grid()->get_security_bandwidth() < 2 )
   {
     string error_message="   >= 2 with TVD scheme";
     MAC_Error::object()->raise_bad_data_value( exp,
        "security_bandwidth", error_message );
   }   


   // Diffusion term time accuracy
   if ( exp->has_entry( "DiffusionTimeAccuracy" ) )
     DiffusionTimeAccuracy = exp->int_data( "DiffusionTimeAccuracy" );
   if ( DiffusionTimeAccuracy != 1 && DiffusionTimeAccuracy != 2 )
   {
     string error_message ="   - 1\n   - 2\n   ";
     MAC_Error::object()->raise_bad_data_value( exp,
         "ViscousTimeAccuracy", error_message );
   }


   // Advection term time accuracy
   if ( exp->has_entry( "AdvectionTimeAccuracy" ) )
     AdvectionTimeAccuracy = exp->int_data( "AdvectionTimeAccuracy" );
   if ( AdvectionTimeAccuracy != 1 && AdvectionTimeAccuracy != 2 )
   {
     string error_message ="   - 1\n   - 2\n   ";
     MAC_Error::object()->raise_bad_data_value( exp,
	"AdvectionTimeAccuracy", error_message );
   }      


   // Source term
   if ( exp->has_entry( "SourceTerm" ) )
   {
     b_sourceTerm = true;
     MAC_ContextSimple* c = MAC_ContextSimple::create( this ) ;
     sourceTerm_coords = MAC_DoubleVector::create( c, doubleVector(0) ) ;
     c->extend( MAC_Variable::object( "DV_X" ), sourceTerm_coords ) ;
     sourceTerm_CTX = c ; 
     
     MAC_DataWithContext const* DEFAULT_FORMULA = exp->abstract_data( 
   	sourceTerm_CTX, "SourceTerm", sourceTerm_CTX ) ;
     if ( DEFAULT_FORMULA->data_type() != MAC_Data::DoubleVector )
     {
        MAC_Error::object()->raise_bad_data_type( exp, "SourceTerm", 
		MAC_Data::DoubleVector ) ;
     }
     doubleVector dof_coordinates( TT->primary_grid()->nb_space_dimensions(), 
     	0. );
     sourceTerm_coords->set( dof_coordinates ) ;
     doubleVector vv = DEFAULT_FORMULA->to_double_vector( sourceTerm_CTX );
     if ( vv.size() != 1 ) 
       MAC_Error::object()->raise_data_error( exp, "SourceTerm", 
		"should have a single value" ) ;
     sourceTerm_formula = DEFAULT_FORMULA;
   }
   

   // Heat transfer problem parameters
   if ( my_rank == is_master )
   {
     MAC::out() << endl << "*** Heat transfer solver" << endl
     	<< endl;
     MAC::out() << "   Diffusion term time accuracy = " <<
     	DiffusionTimeAccuracy << endl;
     MAC::out() << "   Advection term time accuracy = " <<
     	AdvectionTimeAccuracy << endl;
     MAC::out() << "   Advection scheme = " << AdvectionScheme << endl;
     MAC::out() << "   Source term = " << ( b_sourceTerm ? "yes" : "no " )
	<< endl;
     MAC::out() << endl << endl;     
   }


   // Build the matrix system
   MAC_ModuleExplorer* se =
     exp->create_subexplorer( 0, "REG_HeatTransferSystem" ) ;
   GLOBAL_EQ = REG_HeatTransferSystem::create( this, se, TT, UU, 
   	DiffusionTimeAccuracy, AdvectionTimeAccuracy ) ;
   se->destroy() ;   
         
}




//---------------------------------------------------------------------------
void
REG_HeatTransfer:: do_one_inner_iteration( 
	FV_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatTransfer:: do_one_inner_iteration" ) ;

   // Compute CFL
   double computed_CFL = UU->compute_CFL( t_it, 0 );
   size_t n_advection_subtimesteps = unsigned( computed_CFL / imposed_CFL ) + 1;
   
   if ( n_advection_subtimesteps == 1 )
   {
     if ( my_rank == is_master )
     {
       MAC::out() << "-----------------------------------------" << 
       		"-------------" << endl;
       MAC::out() << "Temperature Advection-diffusion problem" << endl;
       MAC::out() << "-----------------------------------------" << 
       		"-------------" << endl;
       MAC::out() << "CFL : imposed = " << imposed_CFL << " computed = " 
       		<< computed_CFL << "  Nb of sub time steps = " << 
		n_advection_subtimesteps << endl;
     } 
     
     // Compute temperature advection rhs
     GLOBAL_EQ->assemble_temperature_advection( AdvectionScheme, 0, 
     	- density * heat_capacity, 0 ) ;

     // Compute the advection-diffusion rhs
     GLOBAL_EQ->compute_temperatureAdvectionDiffusion_rhs( 
	b_restart, t_it->iteration_number(), true ) ;	
     
     // Solve the diffusion system with advection in rhs 
     GLOBAL_EQ->TemperatureDiffusion_solver() ;    
   }
   else
   {   
     // Sub-problem: Advection problem
     if ( my_rank == is_master )
     {       
       MAC::out() << "-----------------------------------------" << 
       		"-------------" << endl;
       MAC::out() << "Sub-problem 1: Advection problem" << endl;
       MAC::out() << "-----------------------------------------" << 
       		"-------------" << endl;
       MAC::out() << "CFL : imposed = " << imposed_CFL << " computed = " 
       		<< computed_CFL << "  Nb of sub time steps = " << 
		n_advection_subtimesteps << endl;
       if ( AdvectionTimeAccuracy == 2 ) 
         MAC::out() << "Degenerates as order 1 scheme, affects the whole"
	 	" splitting algorithm" << endl;	  
     }
        	
     for (size_t i=0; i < n_advection_subtimesteps; ++i)
     {
       // Copy back temperature solution in field for i != 0 as for i = 0, the
       // field already has the right values of the temperature
       if ( i ) TT->update_free_DOFs_value( 0, 
       		GLOBAL_EQ->get_solution_temperature() ) ;
       
       // Compute temperature advection rhs
       GLOBAL_EQ->assemble_temperature_advection( AdvectionScheme, 0, 
       	- density * heat_capacity / double(n_advection_subtimesteps), 0 ) ;
	 
       // Solve advection problem
       GLOBAL_EQ->TemperatureAdvection_solver() ;
       
       // Store temperature advection rhs at t^n-1 in case of 2nd order scheme
       // Then if at next time the scheme is again 2nd (no sub-time stepping)
       // the explicit Adams-Bashforth term is properly computed
       if ( AdvectionTimeAccuracy == 2 && i == 0 )
         GLOBAL_EQ->store_ugradT_Nm2( n_advection_subtimesteps );
     }

      
     // Sub-problem: Diffusion problem     
     // Compute the diffusion rhs (note the false boolean as the last parameter)
     GLOBAL_EQ->compute_temperatureAdvectionDiffusion_rhs( 
        b_restart, t_it->iteration_number(), false ) ;
	
     // Solve the diffusion problem	
     if ( my_rank == is_master )
     {       
       MAC::out() << "-------------------------------" << endl;
       MAC::out() << "Sub-problem 2: Diffusion problem" << endl;
       MAC::out() << "-------------------------------" << endl; 
     }
     GLOBAL_EQ->TemperatureDiffusion_solver(); 
   }

   // Copy back temperature solution in field
   TT->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_temperature() ) ;  	
   
}




//---------------------------------------------------------------------------
void
REG_HeatTransfer:: do_before_time_stepping( 
	FV_TimeIterator const* t_it, 
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatTransfer:: do_before_time_stepping" ) ;
   
   // Temperature unsteady matrix
   if ( my_rank == is_master ) 
     MAC::out() << "            Temperature unsteady matrix" << endl;    
   GLOBAL_EQ->assemble_temperature_unsteady_matrix( 
   	density * heat_capacity / t_it->time_step() );

   // Temperature diffusion matrix and rhs
   // Note: we assemble here the total diffusion matrix, in case of 2nd order 
   // Crank-Nicholson scheme for the diffusion term, half thermal conductivity 
   // for the temperature operator is taken care of at the matrix level in
   // REG_HeatTransferSystem:: finalize_constant_matrices       
   if ( my_rank == is_master ) 
     MAC::out() << "            Temperature diffusion matrix & rhs" << endl;
   GLOBAL_EQ->assemble_temperature_diffusion_matrix_rhs( 
   	- thermal_conductivity ); 

   // Temperature source term rhs
   assemble_sourceTerm( GLOBAL_EQ->get_sourceTerm_vector() );
	
   // Synchronize and finalize matrices 
   GLOBAL_EQ->finalize_constant_matrices() ;

   // Initialize temperature
   GLOBAL_EQ->initialize_temperature() ;

   // Do additional reload
   do_additional_reload( basename ) ;	
	     
}




//---------------------------------------------------------------------------
void
REG_HeatTransfer:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatTransfer:: do_after_time_stepping" ) ;  
     
}




//---------------------------------------------------------------------------
void
REG_HeatTransfer:: do_before_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( 
   	"REG_HeatTransfer:: do_before_inner_iterations_stage" ) ;

   // Perform matrix level operations before each time step 
   GLOBAL_EQ->at_each_time_step( );  
      
}




//---------------------------------------------------------------------------
void
REG_HeatTransfer:: do_after_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatTransfer:: do_after_inner_iterations_stage" ) ;
   
   // Compute temperature change over the time step  
   double temperature_time_change = GLOBAL_EQ->compute_temperature_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master )
     MAC::out() << "         Temperature change = " << 
     	MAC::doubleToString( ios::scientific, 14, temperature_time_change ) 
	<< endl;
   
}
  



//---------------------------------------------------------------------------
void
REG_HeatTransfer:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "REG_HeatTransfer:: do_additional_savings" ) ;
  
}




//---------------------------------------------------------------------------
void
REG_HeatTransfer:: do_additional_save_for_restart( 
	FV_TimeIterator const* t_it,
      	size_t const& restartCycleNumber, std::string const& basename )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "REG_HeatTransfer:: do_additional_save_for_restart" ) ;
    
}



//----------------------------------------------------------------------
REG_HeatTransfer:: ~REG_HeatTransfer( void )
//----------------------------------------------------------------------
{}




//---------------------------------------------------------------------------
void 
REG_HeatTransfer:: do_additional_reload( string const& basename )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_HeatTransfer:: do_additional_reload" ) ;

}




//----------------------------------------------------------------------
void
REG_HeatTransfer::add_storable_objects( MAC_ListIdentity* list )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatTransfer:: add_storable_objects" ) ; 

   GLOBAL_EQ->add_storable_objects( list ) ;
   	
}




//----------------------------------------------------------------------
void
REG_HeatTransfer::assemble_sourceTerm( LA_Vector* VEC_rhs )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatTransfer:: assemble_sourceTerm" ) ; 

   if ( my_rank == is_master ) cout << "Temperature body term rhs "
   	<< endl;

   // Parameters
   double dxC, dyC, dzC, bodyterm = 0. ;
   size_t center_pos_in_matrix = 0 ;
   size_t dim = TT->primary_grid()->nb_space_dimensions();
   size_t k = 0 ;
   doubleVector dof_coordinates( dim, 0. );   
   
   // Get local min and max indices
   size_t_vector min_unknown_index( dim, 0 );
   for (size_t l=0;l<dim;++l) 
     min_unknown_index(l) = 
       	TT->get_min_index_unknown_handled_by_proc( 0, l ) ;
   size_t_vector max_unknown_index( dim, 0 );
   for (size_t l=0;l<dim;++l) 
     max_unknown_index(l) = 
       	TT->get_max_index_unknown_handled_by_proc( 0, l ) ;
	
   // Perform assembling
   for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
   {          
     dxC = TT->get_cell_size( i, 0, 0 ) ;      
     dof_coordinates(0) = TT->get_DOF_coordinate( i, 0, 0 ) ;       
     for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
     {
       dyC = TT->get_cell_size( j, 0, 1 ) ; 
       dof_coordinates(1) = TT->get_DOF_coordinate( j, 0, 1 ) ;
	 
       if ( dim == 2 )
       {
	 sourceTerm_coords->set( dof_coordinates ) ;
	 doubleVector const& val = sourceTerm_formula->to_double_vector() ;
	 bodyterm = val( 0 ) ;
	 center_pos_in_matrix = TT->DOF_global_number( i, j, k, 0 );
	 VEC_rhs->set_item( center_pos_in_matrix, bodyterm * dxC * dyC );
       }
       else
       {
	 for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	 {
	   dzC = TT->get_cell_size( k, 0, 2 ) ;
	   dof_coordinates(2) = TT->get_DOF_coordinate( k, 0, 2 ) ;
	   sourceTerm_coords->set( dof_coordinates ) ;
	   doubleVector const& val = sourceTerm_formula->to_double_vector() ;
	   bodyterm = val( 0 ) ;	     
	   center_pos_in_matrix = TT->DOF_global_number( i, j, k, 0 );
	   VEC_rhs->set_item( center_pos_in_matrix, 
	   	bodyterm * dxC * dyC * dzC );
	 }
       } 
     }      	       	   
   }   

  // Synchronize vector for parallel usage
  VEC_rhs->synchronize();    
   	
}
