#include <VPH_HeatTransfer.hh>
#include <MAC.hh>
#include <MAC_Exec.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Timer.hh>
#include <MAC_Vector.hh>
#include <MAC_Communicator.hh>
#include <MAC_ListIdentity.hh>
#include <FV_DiscreteField.hh>
#include <FV_TimeIterator.hh>
#include <FV_DomainAndFields.hh>
#include <FV_Mesh.hh>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>



//----------------------------------------------------------------------
VPH_HeatTransfer*
VPH_HeatTransfer:: create( MAC_Object* a_owner,
          MAC_ModuleExplorer const* exp,
          struct NavierStokes2Temperature const& transfert )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransfer:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;   

   VPH_HeatTransfer* result = 
         new VPH_HeatTransfer( a_owner, exp, transfert ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;	  
         
   return( result ) ;
   
}




//----------------------------------------------------------------------
VPH_HeatTransfer:: VPH_HeatTransfer( MAC_Object* a_owner,
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
{
   MAC_LABEL( "VPH_HeatTransfer:: VPH_HeatTransfer" ) ;
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
     MAC::out() << endl << endl;     
   }


   // Build the matrix system
   MAC_ModuleExplorer* se =
     exp->create_subexplorer( 0, "VPH_HeatTransferSystem" ) ;
   GLOBAL_EQ = VPH_HeatTransferSystem::create( this, se, TT, UU, 
   	DiffusionTimeAccuracy, AdvectionTimeAccuracy ) ;
   se->destroy() ;   
         
}




//---------------------------------------------------------------------------
void
VPH_HeatTransfer:: do_one_inner_iteration( 
	FV_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransfer:: do_one_inner_iteration" ) ;

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
VPH_HeatTransfer:: do_before_time_stepping( 
	FV_TimeIterator const* t_it, 
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransfer:: do_before_time_stepping" ) ;
   
   // Temperature unsteady matrix
   if ( my_rank == is_master ) 
     MAC::out() << "            Temperature unsteady matrix" << endl;    
   GLOBAL_EQ->assemble_temperature_unsteady_matrix( 
   	density * heat_capacity / t_it->time_step() );

   // Temperature diffusion matrix and rhs
   // Note: we assemble here the total diffusion matrix, in case of 2nd order 
   // Crank-Nicholson scheme for the diffusion term, half thermal conductivity 
   // for the temperature operator is taken care of at the matrix level in
   // VPH_HeatTransferSystem:: finalize_constant_matrices       
   if ( my_rank == is_master ) 
     MAC::out() << "            Temperature diffusion matrix & rhs" << endl;
   GLOBAL_EQ->assemble_temperature_diffusion_matrix_rhs( 
   	- thermal_conductivity ); 
	
   // Synchronize and finalize matrices 
   GLOBAL_EQ->finalize_constant_matrices() ;

   // Initialize temperature
   GLOBAL_EQ->initialize_temperature() ;

   // Do additional reload
   do_additional_reload( basename ) ;	
	     
}




//---------------------------------------------------------------------------
void
VPH_HeatTransfer:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransfer:: do_after_time_stepping" ) ;  
     
}




//---------------------------------------------------------------------------
void
VPH_HeatTransfer:: do_before_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( 
   	"VPH_HeatTransfer:: do_before_inner_iterations_stage" ) ;

   // Perform matrix level operations before each time step 
   GLOBAL_EQ->at_each_time_step( );  
      
}




//---------------------------------------------------------------------------
void
VPH_HeatTransfer:: do_after_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransfer:: do_after_inner_iterations_stage" ) ;
   
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
VPH_HeatTransfer:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "VPH_HeatTransfer:: do_additional_savings" ) ;
  
}




//---------------------------------------------------------------------------
void
VPH_HeatTransfer:: do_additional_save_for_restart( 
	FV_TimeIterator const* t_it,
      	size_t const& restartCycleNumber, std::string const& basename )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "VPH_HeatTransfer:: do_additional_save_for_restart" ) ;
    
}



//----------------------------------------------------------------------
VPH_HeatTransfer:: ~VPH_HeatTransfer( void )
//----------------------------------------------------------------------
{}




//---------------------------------------------------------------------------
void 
VPH_HeatTransfer:: do_additional_reload( string const& basename )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "VPH_HeatTransfer:: do_additional_reload" ) ;

}




//----------------------------------------------------------------------
void
VPH_HeatTransfer::add_storable_objects( MAC_ListIdentity* list )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransfer:: add_storable_objects" ) ; 

   GLOBAL_EQ->add_storable_objects( list ) ;
   	
}
