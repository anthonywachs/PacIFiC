#include <MAC_ApplicationRestorer.hh>

#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Communicator.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_Root.hh>
#include <MAC.hh>
#include <MAC_assertions.hh>

using std::string ;

MAC_ApplicationRestorer const* 
MAC_ApplicationRestorer:: PROTOTYPE = new MAC_ApplicationRestorer() ;



//-------------------------------------------------------------------------
MAC_ApplicationRestorer:: MAC_ApplicationRestorer( void )
//-------------------------------------------------------------------------
   : MAC_Application( "MAC_ApplicationRestorer" )
   , READER( 0 )
   , APPLI( 0 )
{
}



   
//---------------------------------------------------------------------------
MAC_ApplicationRestorer*
MAC_ApplicationRestorer:: create_replica( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ApplicationRestorer:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   MAC_ApplicationRestorer* result = new MAC_ApplicationRestorer( a_owner,
	exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
MAC_ApplicationRestorer:: MAC_ApplicationRestorer( 
	MAC_Object* a_owner, MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : MAC_Application( a_owner, exp, true, true )
   , READER( 0 )
   , APPLI(0)
{
   MAC_ModuleExplorer* ee = 0 ;

   ee = exp->create_subexplorer( 0, "MAC_ObjectReader" ) ;
   READER = MAC_ObjectReader::create( this, ee ) ;
   ee->destroy() ;

   MAC_Module* mod = create_modified_data_deck_module( exp ) ;   
   ee = MAC_ModuleExplorer::create( 0, mod ) ;
   double initial_time = READER->get_initial_time() ; 
   APPLI = MAC_Application::make( this, ee, initial_time ) ;
   ee->destroy() ;

   size_t cycle_number = 0 ;
   if ( exp->has_entry( "cycle_number" ) )
     cycle_number = exp->int_data( "cycle_number" ) ;

   if ( MAC_Exec::communicator()->rank() == 0 ) READER->print( MAC::out(), 0 );   
   READER->seek_cycle( cycle_number ) ;
   APPLI->restore_registered_objects( READER ) ;
   READER->close_cycle() ;
   if ( MAC_Exec::communicator()->rank() == 0 ) MAC::out() << std::endl;
   
   MAC_CHECK_POST( APPLI!=0 ) ;
}




//---------------------------------------------------------------------------
MAC_ApplicationRestorer:: ~MAC_ApplicationRestorer( void )
//---------------------------------------------------------------------------
{
}




//-------------------------------------------------------------------------
void
MAC_ApplicationRestorer:: run( std::string const& inputRestartFileName ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ApplicationRestorer:: run" ) ;

   APPLI->run( READER->input_file_name() ) ;
}




//---------------------------------------------------------------------------
MAC_Module*
MAC_ApplicationRestorer:: create_modified_data_deck_module( 
	MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ApplicationRestorer:: create_modified_data_deck_module" ) ;

   MAC_Module* m = 0 ;

   MAC_Module* header = READER->header_module() ;
   if( header->has_module( "MAC_Application" ) )
      m = header->module( "MAC_Application" ) ;
   else
      MAC_Error::object()->raise_plain( "invalid restart file" ) ; 

   MAC_ASSERT( m != 0 ) ;
   MAC_Module* result = m->create_clone( this ) ;

   result->remove_module( "MAC_ObjectWriter" ) ;
   
   if ( result->has_entry( "graphics_output_times" ) )
     result->remove_entry( "graphics_output_times" );
   if ( result->has_entry( "number_graphics_output_times" ) )
     result->remove_entry( "number_graphics_output_times" );     
   if ( result->has_entry( "state_saving_times" ) )
     result->remove_entry( "state_saving_times" );     
   if ( result->has_entry( "number_state_saving_times" ) )
     result->remove_entry( "number_state_saving_times" );     
            
   change_owner( MAC_Root::object(), result ) ;

   if ( exp->has_entry( "appendum_file_name" ) )
   {
      string appendum_file = exp->string_data( "appendum_file_name" ) ;
      MAC_Module* appendum =
         MAC_Module::create( 0, "ROOT", appendum_file,
                             MAC_Exec::execution_context() ) ;
      MAC_Module* o = appendum->module( "MAC_Application" ) ;
      result->merge_module( o ) ;
      appendum->destroy() ; appendum = 0 ;
   }
   else if ( exp->has_module( "MAC_Application" ) )
   {
     MAC_ModuleExplorer* se =
     	exp->create_subexplorer( 0, "MAC_Application" ) ;
     MAC_Module* o = se->create_clone_of_attached_module( 0 );
     result->merge_module( o ) ;
     o->destroy() ; o = 0 ;
     se->destroy() ; se = 0 ;
   }
   else
     MAC_Error::object()->raise_missing_keyword( exp,
	"appendum_file_name\" or \"MODULE MAC_Application" );
   
   MAC_CHECK( result != 0 ) ;
   MAC_CHECK( result->owner() == MAC_Root::object() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
std::string
MAC_ApplicationRestorer:: input_restart_file_name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ApplicationRestorer:: input_restart_file_name" ) ;

   return( READER->input_file_name() );
}
