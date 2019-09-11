#include <MAC_Application.hh>

#include <MAC.hh>
#include <MAC_Error.hh>
#include <MAC_ListIdentity.hh>
#include <MAC_ListIterator.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_ObjectWriter.hh>
#include <MAC_Root.hh>
#include <MAC_System.hh>
#include <MAC_assertions.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;

struct MAC_Application_ERROR
{
   static void n0( stringVector const& args ) ;
   static void n1( std::string const& class_name ) ;
   static void n2( std::string const& method ) ;   
} ;


//------------------------------------------------------------------------
bool MAC_Application::B_FOLLOW = false ;
bool MAC_Application::B_RELOAD = false ;
//------------------------------------------------------------------------


//----------------------------------------------------------------------
MAC_Application*
MAC_Application:: make( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: make" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   string name = exp->string_data( "concrete_name" ) ;
   MAC_Application const* proto =
      static_cast<MAC_Application const*>(
                                    plugins_map()->item( name ) ) ;
   MAC_ASSERT( proto->is_a_prototype() ) ;
      
   MAC_Application* result = proto->create_replica( a_owner, exp,
   	initial_time ) ;
   result->initialize_objects_storage( exp ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}





// Create application from args NOT SUPPORTED for now
// //----------------------------------------------------------------------
// MAC_Application*
// MAC_Application:: make( MAC_Object* a_owner,
//                         stringVector& args )
// //----------------------------------------------------------------------
// {
//    MAC_LABEL( "MAC_Application:: make" ) ;
//    MAC_CHECK_PRE( args.size() != 0 ) ;
// 
//    if( args( 0 ) != "-A" ) MAC_Application_ERROR::n0( args ) ;
//    args.remove_at( 0 ) ;
//    string name = args( 0 ) ;
//    args.remove_at( 0 ) ;
//    MAC_Application const* proto =
//       static_cast<MAC_Application const*>(
//                                     plugins_map()->item( name ) ) ;
//    MAC_ASSERT( proto->is_a_prototype() ) ;
//       
//    MAC_Application* result = proto->create_replica_from_args( a_owner, args ) ;
//    if( args.size() != 0 ) result->print_usage_then_exit() ;
//    
//    result->initialize_objects_storage( 0 ) ;
//    
//    MAC_CHECK_POST( result != 0 ) ;
//    MAC_CHECK_POST( result->owner() == a_owner ) ;
//    return( result ) ;
// }





//----------------------------------------------------------------------
MAC_Application:: MAC_Application( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	bool const& is_follow_,
	bool const& is_reload_ )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , IS_PROTO( false )
   , SAVER( 0 )
   , persistent_objects( 0 )
   , persistent_objects_it( 0 )
{
   persistent_objects = MAC_ListIdentity::create( this ) ;
   persistent_objects_it = MAC_ListIterator::create( this,
                                                     persistent_objects ) ;

   // As soon as one implemented MAC application set B_FOLLOW or B_RELOAD
   // to true, it cannot be unchanged by subsequent implementations of other
   // MAC applications.
   // Only MAC_ApplicationFieldReloader and MAC_ApplicationFollower are allowed 
   // to modify B_FOLLOW or B_RELOAD
   if ( MAC_Application::B_FOLLOW == false )
     MAC_Application::B_FOLLOW = is_follow_;
   if ( MAC_Application::B_RELOAD == false )
     MAC_Application::B_RELOAD = is_reload_;     

   MAC_CHECK_POST( owner() == a_owner ) ;
   MAC_CHECK_POST( !is_a_prototype() ) ;
}





//----------------------------------------------------------------------
MAC_Application:: MAC_Application( std::string const& name )
//----------------------------------------------------------------------
   : MAC_Object( plugins_map() )
   , IS_PROTO( true )
   , SAVER( 0 )
   , persistent_objects( 0 )
   , persistent_objects_it( 0 )
{
   MAC_LABEL( "MAC_Application:: MAC_Application" ) ;
   
   plugins_map()->register_item( name, this ) ;

   // This constructor is used for plugins only and leaves the static data
   // B_FOLLOW and B_RELOAD unchanged to their default value FALSE
   
   MAC_CHECK_POST( is_under_ownership_of( plugins_map() ) ) ;
   MAC_CHECK_POST( is_a_prototype() ) ;
}





//----------------------------------------------------------------------
MAC_Application:: ~MAC_Application( void )
//----------------------------------------------------------------------
{
}





//----------------------------------------------------------------------
void
MAC_Application:: register_storable_objects( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: register_storable_objects" ) ;

   add_storable_objects( persistent_objects ) ;
}




//----------------------------------------------------------------------
void
MAC_Application:: add_storable_objects( MAC_ListIdentity* list )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: add_storable_objects" ) ;
}




//----------------------------------------------------------------------
void
MAC_Application:: write_storable_objects( double const& time ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: write_storable_objects" ) ;

   if( SAVER != 0 )
   {
      SAVER->start_cycle( time ) ;

      MAC::out() << std::endl << "---Saving (cycle:"
                 << SAVER->cycle_number()
                 << ")---" << std::endl;		 

      SAVER->start_new_object( "MAC_Application" ) ;

      for( persistent_objects_it->start();
	   persistent_objects_it->is_valid() ;
	   persistent_objects_it->go_next() )
      {
	 persistent_objects_it->item()->save_state( SAVER ) ;
      }

      SAVER->finalize_object() ;

      SAVER->terminate_cycle() ;
   }
}




//----------------------------------------------------------------------
void
MAC_Application:: restore_registered_objects( MAC_ObjectReader* ret ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: restore_registered_objects" ) ;

   ret->start_object_retrieval( "MAC_Application" ) ;

   for( persistent_objects_it->start();
        persistent_objects_it->is_valid() ;
        persistent_objects_it->go_next() )
   {
      persistent_objects_it->item()->restore_state( ret ) ;
   }

   ret->end_object_retrieval() ;
}




//----------------------------------------------------------------------
bool
MAC_Application:: is_follow( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: is_follow" ) ;

   return( B_FOLLOW );
}




//----------------------------------------------------------------------
bool
MAC_Application:: is_reload( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: is_reload" ) ;

   return( B_RELOAD );
}




//----------------------------------------------------------------------
size_t
MAC_Application:: restart_cycle_number( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: restart_cycle_number" ) ;
   MAC_CHECK_PRE( SAVER != 0 ) ;
   
   if ( SAVER == 0 ) MAC_Application_ERROR::n2( "restart_cycle_number" );

   return( SAVER->cycle_number() );
}




//----------------------------------------------------------------------
std::string
MAC_Application:: output_restart_file_name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: output_restart_file_name" ) ;
   MAC_CHECK_PRE( SAVER != 0 ) ;

   if ( SAVER == 0 ) MAC_Application_ERROR::n2( "output_restart_file_name" );

   return( SAVER->get_current_file_name() );
}




//----------------------------------------------------------------------
std::string
MAC_Application:: input_restart_file_name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: output_restart_file_name" ) ;

   MAC_Application_ERROR::n2( "input_restart_file_name" );

   return( MAC::undefined_string );
}




//----------------------------------------------------------------------
void
MAC_Application:: swap_restart_file_names( std::string const& inputfilename ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: swap_restart_file_names" ) ;
   
   if ( SAVER != 0 )    
     SAVER->swap_restart_file_names( inputfilename );  
}




//----------------------------------------------------------------------
bool
MAC_Application:: force_write_restart_files_at_last_time( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: force_write_restart_files_at_last_time" ) ;

   return( SAVER ? SAVER->force_write_at_last_time() : false );
}




//----------------------------------------------------------------------
std::string const&
MAC_Application:: object_writer_module_name( void ) const
//----------------------------------------------------------------------
{
   static std::string const result = "MAC_ObjectWriter" ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
MAC_Application:: create_replica_PRE( MAC_Object* a_owner,
                                      MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_a_prototype() ) ;
   MAC_ASSERT( exp != 0 ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Application:: create_replica_POST( MAC_Application const* result,
                                       MAC_Object* a_owner,
                                       MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   MAC_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}




// Create application from args NOT SUPPORTED for now
// //----------------------------------------------------------------------
// bool
// MAC_Application:: create_replica_from_args_POST( MAC_Application const* result,
//                                                  MAC_Object* a_owner,
//                                                  stringVector& args ) const
// //----------------------------------------------------------------------
// {
//    MAC_ASSERT( result != 0 ) ;
//    MAC_ASSERT( result->owner() == a_owner ) ;
//    MAC_ASSERT( !result->is_a_prototype() ) ;
//    return( true ) ;
// }




//----------------------------------------------------------------------
bool
MAC_Application:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   return( true ) ;
}




// Create application from args NOT SUPPORTED for now
// //----------------------------------------------------------------------
// MAC_Application*
// MAC_Application:: create_replica_from_args( MAC_Object* a_owner,
//                                             stringVector& args ) const
// //----------------------------------------------------------------------
// {
//    MAC_LABEL( "MAC_Application:: create_replica_from_args" ) ;
//    
//    MAC_Application* result = 0 ;
// 
//    MAC_Application_ERROR::n1( type_name() ) ;
// 
//    MAC_CHECK_POST( create_replica_from_args_POST( result, a_owner, args ) ) ;
//    return( result ) ;
// }




//----------------------------------------------------------------------
bool
MAC_Application:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}




//----------------------------------------------------------------------
void
MAC_Application:: notify_error_in_arguments( void ) const
//----------------------------------------------------------------------
{
   print_usage_then_exit( -1 ) ;
}




//----------------------------------------------------------------------
void
MAC_Application:: print_usage( void ) const
//----------------------------------------------------------------------
{
}




//----------------------------------------------------------------------
void
MAC_Application:: print_options( void ) const
//----------------------------------------------------------------------
{
}




//----------------------------------------------------------------------
void
MAC_Application:: print_operands( void ) const
//----------------------------------------------------------------------
{
}




//----------------------------------------------------------------------
void
MAC_Application:: print_exit_status( void ) const
//----------------------------------------------------------------------
{
}




//----------------------------------------------------------------------
std::string
MAC_Application:: usage_title( std::string const& name ) const
//----------------------------------------------------------------------
{
   std::string result = "USAGE\n" ;
   result += "     <exe> -A " + name ;
   result += " [-h] " ;
   return( result ) ;
}




//----------------------------------------------------------------------
std::string
MAC_Application:: options_title( void ) const
//----------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "OPTIONS" << endl ;
   mesg << "     -h   Display help and exit" << endl << endl ;
   return mesg.str() ;
}




//----------------------------------------------------------------------
std::string
MAC_Application:: operands_title( void ) const
//----------------------------------------------------------------------
{
   std::string result = "OPERANDS\n" ;
   return( result ) ;
}




//----------------------------------------------------------------------
std::string
MAC_Application:: exit_status_title( void ) const
//----------------------------------------------------------------------
{
   std::string result = "EXIT STATUS\n" ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC_Application:: print_usage_then_exit( int exit_status ) const
//----------------------------------------------------------------------
{
   print_usage() ;
   MAC::out() << endl <<"OPTIONS" << endl ;
   MAC::out() << "     -h   Display help and exit" << endl << endl ;
   print_options() ;
   print_operands() ;
   print_exit_status() ;
   MAC::out() << endl ;
   MAC_System::exit( exit_status ) ;
}




//----------------------------------------------------------------------
void
MAC_Application:: initialize_objects_storage( MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   // exp->has_module( object_writer_module_name() ) guarantees that
   // only 1 application has a SAVER
   // For instance, MAC_ApplicationRestorer will not implemented a SAVER
   // but the application it implements will implement it
   if( ( exp != 0 ) && exp->has_module( object_writer_module_name() ) )
   {
      MAC_ModuleExplorer const* e =
         exp->create_subexplorer( 0, object_writer_module_name() ) ;
      SAVER = MAC_ObjectWriter::create( this, e, exp ) ;   
      e->destroy() ;
   }
   register_storable_objects() ;
}




//----------------------------------------------------------------------
MAC_ObjectRegister*
MAC_Application:: plugins_map( void )
//----------------------------------------------------------------------
{
   static MAC_ObjectRegister* result =
            MAC_ObjectRegister::create( MAC_Root::object(),
                                        "MAC_Application descendant" ) ;
   return( result ) ;
}




//internal---------------------------------------------------------------
void
MAC_Application_ERROR:: n0( stringVector const& args )
//internal---------------------------------------------------------------
{
   std::ostringstream os ;
   os << "invalid command-line arguments : " << endl << "  ";
   for( size_t i=0 ; i<args.size() ; ++i ) os << " " << args(i) ;
   MAC_Error::object()->raise_plain( os.str() ) ;
}




//internal---------------------------------------------------------------
void
MAC_Application_ERROR:: n1( std::string const& class_name )
//internal---------------------------------------------------------------
{
   std::ostringstream os ;
   os << "Class " << class_name << " derived from MAC_Application :" << endl ;
   os << "   instantiation from the command_line is impossible" << endl ;
   os << "   (\"create_replica_from_args\" has not been overridden)." ;
   MAC_Error::object()->raise_plain( os.str() ) ;
}




//internal---------------------------------------------------------------
void
MAC_Application_ERROR:: n2( std::string const& method )
//internal---------------------------------------------------------------
{
   std::ostringstream os ;
   os << "This instance of MAC_Application did not implement the " << endl ;
   os << " MAC_ObjectWriter SAVER" << endl ;
   os << " Method MAC_Application::" << method << " should not be called";
   os << " by this instance of MAC_Application. This is a design problem";
   MAC_Error::object()->raise_plain( os.str() ) ;
}
