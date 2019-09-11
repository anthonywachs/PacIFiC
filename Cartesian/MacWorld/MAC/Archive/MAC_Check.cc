#include <MAC_Check.hh>

#include <MAC.hh>
#include <MAC_Exec.hh>
#include <MAC_Error.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ModulePattern.hh>
#include <MAC_System.hh>
#include <MAC_assertions.hh>

#include <fstream>
#include <iostream>

MAC_Check const* MAC_Check::PROTOTYPE = new MAC_Check( "check" ) ;

/*
  Check validity of several data files from a pattern model.
  
  Two modes are available :
   - in interactive mode, application waits for a command line on standard
     input stream of the form : <filename> ;
   - else files to be checked are given on application command line.

   COMMAND LINE :
   Syntaxes recognized in command line mode are :
   check [-s] -i pattern.mac                        interactive mode
   check [-s] pattern.mac file1.mac ... filen.mac   non interactive mode

   DATA DECK :
   Data files recognized are :
   
   MODULE MAC_Application
      concrete_name = "check"
      pattern_filename = "pattern.mac"
      checked_files = < "file1.mac" file2.mac" ... > // Optionnal
      output_file = "expected.err"                   // Optionnal
   END MODULE MAC_Application
   interactive mode is actived when checked_files is not provided.
   FRAMEWORK INSTANTIATION
   Class wishing inherit from check must implement following methods :
    `do_check(a_owner, exp)' : method that process checking ;
                               
     Plug in methods : create_replica, create_replica_from_args, constructors.
*/

//----------------------------------------------------------------------
MAC_Check:: MAC_Check( std::string const& a_name )
//----------------------------------------------------------------------
   : MAC_Application( a_name )
   , PATTERN()
   , MY_ARGS( 0 )
   , FILE()
   , SILENT( false )
   , NAME( a_name )
{
}

//----------------------------------------------------------------------
MAC_Check:: MAC_Check( MAC_Object* a_owner,
                       std::string const& a_name,
                       MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : MAC_Application( a_owner, exp )
   , PATTERN()
   , MY_ARGS( 0 )
   , FILE( "" )
   , SILENT( false )
   , INTERACTIVE( false )
   , NAME( a_name )
{
   MAC_LABEL( "MAC_Check:: MAC_Check" ) ;
   PATTERN = exp->string_data( "pattern_filename" ) ;
   if(exp->has_entry("checked_files"))
   {
      MY_ARGS = exp->stringVector_data( "checked_files" ) ;
   }
   else
   {
      INTERACTIVE = true ;
   }
   if( exp->has_entry( "silent_mode" ) )
   {
      SILENT = exp->bool_data( "silent_mode" ) ;
   }
   
   if( exp->has_entry("output_file") )
   {
      FILE=exp->string_data( "output_file" ) ;
   }
}

//----------------------------------------------------------------------
MAC_Check:: MAC_Check( MAC_Object* a_owner,
                       std::string const& a_name,
                       stringVector& args )
//----------------------------------------------------------------------
   : MAC_Application( a_owner, 0 )
   , PATTERN()
   , MY_ARGS( 0 )
   , FILE( "" )
   , SILENT( false )
   , INTERACTIVE( false )
   , NAME( a_name )
{
   MAC_LABEL( "MAC_Check:: MAC_Check" ) ;

   for( size_t i=0 ; i<args.size() ; i++ )
   {
      if(args(i)=="-s" )
      {
         SILENT = true ;
         args.remove_at(i) ;
         break ;
      }
      else if(args(i)=="-i" )
      {
         INTERACTIVE = true ;
         args.remove_at(i) ;
         break ;
      }
   }
   
   if( !INTERACTIVE && args.size() < 2 )
   {
      MAC_Error::object()->raise_plain(
         "usage : " + NAME + " pattern.mac data1.mac [datai.mac...]" ) ;
   }
   else if( INTERACTIVE && args.size() != 1 )
   {
      MAC_Error::object()->raise_plain(
         "usage : " + NAME + " -i pattern.mac " ) ;
   }
   PATTERN = args(0) ;
   args.remove_at(0) ;
   MY_ARGS = args ;
   args.re_initialize(0) ;
}

//----------------------------------------------------------------------
MAC_Check:: ~MAC_Check( void )
//----------------------------------------------------------------------
{
   if( this == PROTOTYPE )
   {
      PROTOTYPE = 0 ;  
   }
}

//----------------------------------------------------------------------
MAC_Check* 
MAC_Check:: create_replica( MAC_Object* a_owner,
                            MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Check:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   MAC_Check* result = new MAC_Check( a_owner, NAME, exp ) ;

   MAC_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return result ;
}

//----------------------------------------------------------------------
MAC_Check* 
MAC_Check:: create_replica_from_args( MAC_Object* a_owner,
                                      stringVector& args ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Check:: create_replica_from_args" ) ;

   MAC_Check* result = new MAC_Check( a_owner, NAME, args ) ;

   MAC_CHECK_POST( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return result ;
}

//----------------------------------------------------------------------
void
MAC_Check:: run( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Check:: run" ) ;
   if( MAC_ModulePattern::build_pattern() )
   {
      MAC_ModulePattern::save_pattern() ;
   }
   
   MAC_ModulePattern::open_pattern_base( PATTERN ) ;
   
   if( !SILENT )
      MAC_Exec::out() << "****** " << NAME << " output" << std::endl << std::endl ;

   if( !INTERACTIVE )
   {
      if( !SILENT )
         MAC_Exec::out() << NAME << " : interactive mode is off"
                         << std::endl << std::endl ;
      for( size_t i=0 ; i<MY_ARGS.size() ; i++ )
      {
         process( MY_ARGS(i) ) ;         
      }
   }
   else 
   {
      bool end = false ;
      MAC::out() << NAME << " : interactive mode is on" << std::endl ;
      MAC::out() << " Enter name of file to check (QUIT to quit)" << std::endl ;
      while( !end ) 
      {
         std::string a_file("QUIT") ;
      
         std::cin >> a_file ;
         if( a_file == "QUIT" ) 
         {
            end = true ;
         }
         else
         {
            process( a_file ) ;
            MAC::out() << "That's all folks" << std::endl ;
            MAC::out().flush() ;
         }
      }
   }

   MAC_ModulePattern::close_pattern_base() ;
}

//----------------------------------------------------------------------
void
MAC_Check:: process( std::string const& file_to_parse )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Check:: process" ) ;
   
   if( !SILENT )
      MAC_Exec::out() << "*** Checking " << file_to_parse << std::endl ;
      
   stringVector args(1) ;
   args(0) = file_to_parse ;
   MAC_ModuleExplorer* validity = 0 ;
   MAC_Module const* mod = MAC_Exec::create_module_with_data( this, args ) ;
   if( mod!=0 ) // reading successful
   {
      MAC_ModuleExplorer* exp = MAC_ModuleExplorer::create( 0, mod ) ;
      validity = do_check( 0, exp ) ;
      exp->destroy() ; exp=0 ;
   }
   else
   {
      std::string const& parsed =
         MAC_Module::current_parsed_module_path_name() ;
         
      if( !parsed.empty() )
      {
         MAC_Module* res = MAC_Module::create( this, "ParsingError" ) ;
         MAC_Error::notify( res, "Parse error", parsed ) ;
         validity = MAC_ModuleExplorer::create(0, res ) ;
      }
   }
   if( !validity->is_empty() )
   {
      if( FILE.empty() )
      {
         MAC_Exec::out() << std::endl ;
         MAC_Error::object()->display_data_checking( validity ) ;
      }
      else
      {
         std::ofstream out( FILE.c_str() ) ;
         if( !out )
            MAC_Error::object()->raise_plain( "Unable to open "+FILE ) ;
            
         validity->print( out, 0 ) ;
         out.close() ;         
      }
      MAC_Exec::set_exit_code( 1 ) ;
   }
   else
   {
      if(FILE.empty())
      {
         if( !SILENT )
         MAC_Exec::out() << std::endl << "Pattern checking is successful"
                         << std::endl ;
      } else 
      {
         MAC_System::erase( FILE ) ;
      }
      
   }
   
   validity->destroy() ;
      
   if( mod!=0 ) destroy_possession( mod ) ;
   mod = 0 ;
}

//----------------------------------------------------------------------
MAC_ModuleExplorer*
MAC_Check:: do_check( MAC_Object* a_owner,
                      MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Check:: do_check" ) ;
   MAC_CHECK_PRE( do_check_PRE( a_owner, exp ) ) ;

   MAC_Module* a_mod = exp->create_clone_of_attached_module(0) ;
   
   MAC_ModuleExplorer* ctrl =
      MAC_ModuleExplorer::create( 0,
                                  a_mod,
                                  MAC_ModuleExplorer::verify ) ;
   
   MAC_ModuleExplorer* result = ctrl->validity(a_owner) ;

   ctrl->destroy() ;
   a_mod->destroy() ;
   
   MAC_CHECK_POST( do_check_POST( a_owner, result ) ) ;
   return result ;
}


//----------------------------------------------------------------------
bool
MAC_Check:: do_check_PRE( MAC_Object* a_owner,
                          MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( exp!=0 ) ;
   return true ;
}

//----------------------------------------------------------------------
bool
MAC_Check:: do_check_POST( MAC_Object* a_owner,
                            MAC_ModuleExplorer const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result!=0 ) ;
   MAC_ASSERT( result->owner()==a_owner ) ;
   return true ;
}

