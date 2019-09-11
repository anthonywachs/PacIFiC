// LEVEL macro is undefined in MAC_assertions.hh
int const a_compilation_level = LEVEL ;

#include <MAC.hh>
#include <MAC_Exec.hh>

#include <MAC_assertions.hh>
#include <MAC_Application.hh>
#include <MAC_Bool.hh>
#include <MAC_Communicator.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_Error.hh>
#include <MAC_Exceptions.hh>
#include <MAC_ExternalAPI.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_ModulePattern.hh>
#include <MAC_Object.hh>
#include <MAC_Root.hh>
#include <MAC_String.hh>
#include <MAC_System.hh>
#include <MAC_Timer.hh>
#include <MAC_Variable.hh>

#include <stringVector.hh>

#include <fstream>
#include <iostream>
#include <sstream>
using std::endl ;

//------------------------------------------------------------------------
bool MAC_Exec::VERBOSE = false ;
bool MAC_Exec::ONLY_HELP = false ;
bool MAC_Exec::SIGNAL_HANDLING = true ;
MAC_Timer* MAC_Exec::APPLI_TIMER = 0 ;
std::string MAC_Exec::EXE_FILE = "" ;
int MAC_Exec::EXIT_STATUS = 0 ;
MAC_ModuleExplorer::PatternStatus MAC_Exec::PATTERN_STATUS =
                                              MAC_ModuleExplorer::ignore ;
MAC_ContextSimple* MAC_Exec::EXEC_CONTEXT = 0 ;
std::istream* MAC_Exec::MAC_IN = 0 ;
std::ostream* MAC_Exec::MAC_OUT = 0 ;
std::ostream* MAC_Exec::MAC_ERR = 0 ;
bool MAC_Exec::EXTERNAL_API = true ;
//------------------------------------------------------------------------


//------------------------------------------------------------------------
size_t
MAC_Exec:: compilation_level( void )
//-------------------------------------------------------------------------
{
   return( (size_t) a_compilation_level ) ;
}




//------------------------------------------------------------------------
std::string const&
MAC_Exec:: date_of_compilation( void )
//-------------------------------------------------------------------------
{
   static std::string const result = __DATE__ ;
   return( result ) ;
}




//------------------------------------------------------------------------
std::string const&
MAC_Exec:: name_of_exe( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec:: name_of_exe" ) ;
   return( EXE_FILE );
}




//------------------------------------------------------------------------
std::string
MAC_Exec:: dynamic_check_list( void )
//-------------------------------------------------------------------------
{
   std::string result ;
   if( MAC_Assertion::is_handling_check( MAC_Assertion::Precondition ) )
   {
      result += "preconditions " ;
   }
   if( MAC_Assertion::is_handling_check( MAC_Assertion::Postcondition ) )
   {
      result += "postconditions " ;
   }
   if( MAC_Assertion::is_handling_check( MAC_Assertion::Invariant ) )
   {
      result += "invariants " ;
   }
   if( MAC_Assertion::is_handling_check( MAC_Assertion::Check ) )
   {
      result += "checks " ;
   }
   if( MAC_Assertion::is_handling_check( MAC_Assertion::Objects ) )
   {
      result += "objects " ;
   }
   
   return( result ) ;
}




//------------------------------------------------------------------------
int
MAC_Exec::exit_code( void )
//------------------------------------------------------------------------
{
   return( EXIT_STATUS ) ;
}




//------------------------------------------------------------------------
void
MAC_Exec::set_exit_code( int status )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec::set_exit_code" ) ;
   EXIT_STATUS = status ;
   MAC_CHECK_POST( exit_code() == status ) ;
}




//------------------------------------------------------------------------
MAC_Module const* 
MAC_Exec:: create_module_with_data( MAC_Object* a_owner,
                                    stringVector& args )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec:: create_module_with_data" ) ;
   MAC_Module* result = 0 ;
   try 
   { 
      MAC_CHECK_PRE( a_owner!=0 ) ; 

      if( !ONLY_HELP )
      {
         std::string dataFile = "" ;
         MAC_Module* main_module = 0 ;

         if( args.size()!=0 && args(0)=="-MAC_Application" )
         {
//             // la structure de modules est écrite en ligne....
//             MAC_ASSERT( args.size() > 1 ) ;
//             args.remove_at(0) ;
//             std::string api = args(0) ;
//             args.remove_at(0) ;
//             main_module = create_module_from_args( 0, api, args ) ;
//             dataFile = "" ;
//          Create module from args NOT SUPPORTED for now
	    MAC_Error::object()->raise_plain(
                     "Create module from args NOT SUPPORTED for now" ) ;
         }
         else if( args.size() == 1 )
         {
            dataFile = args(0) ;
            args.remove_at(0) ;
         }
      
         // Display MAC header.
         if( VERBOSE ) print( dataFile, MAC_Exec::out() ) ;
   
         if( ( main_module==0 ) && ( !dataFile.empty() ) )
         {
            // data-deck reading and storage in memory
            main_module = MAC_Module::create( 0, "MAIN", dataFile,
                                              execution_context() ) ;
         }

         if( main_module!=0 ) 
         {
            MAC_ModuleIterator* it = main_module->create_module_iterator( 0 ) ;
            it->start() ;
            if( !it->is_valid() ) 
            {
               MAC_Error::object()->raise_plain(
                     "Empty file " + dataFile +"\nNo root module found" ) ;
            }
            
            result = it->item() ; it->go_next() ;
            if( it->is_valid() ) 
            {
               MAC_Error::object()->raise_plain(
                     "No only one root module found in "+dataFile ) ;
            }
            it->destroy() ;
            main_module->remove_module( result->name() ) ;
            main_module->change_owner( a_owner, result ) ;
            main_module->destroy() ; main_module = 0 ;
         }
      }

      MAC_CHECK_POST( IMPLIES( result!=0, result->owner() == a_owner ) ) ;
   }
   catch( MAC_Exceptions::Error )
   {
      result = 0 ;
   }
   return( result ) ;  
}




//------------------------------------------------------------------------
MAC_Application* 
MAC_Exec:: create_application( MAC_Object* a_owner,
                               MAC_Module const* appli_mod,
                               stringVector& args )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec:: create_application" ) ;
   
   MAC_Application* result = 0 ;

   if( !ONLY_HELP && appli_mod!=0 )
   {
      // creation of an instance of a concrete subclass of MAC_Application
      try
      {
         MAC_ModuleExplorer* exp =
            MAC_ModuleExplorer::create( a_owner,
                                        appli_mod,
                                        PATTERN_STATUS ) ;
         if( PATTERN_STATUS==MAC_ModuleExplorer::verify )
         {
            MAC_ModuleExplorer* validity = exp->validity(0) ;
            if( !validity->is_empty() )
            {
               MAC_Error::object()->display_data_checking(
                  validity ) ;
               MAC_Error::object()->raise_plain( "Pattern check aborts" ) ;
            }
            validity->destroy() ; validity = 0 ;
         }
         result = MAC_Application::make( a_owner, exp ) ;
      }
      catch( MAC_Exceptions::Error )
      {
         result = 0 ;
      }
   }
// Create application from args NOT SUPPORTED for now
//    else if( !ONLY_HELP && args.size()!=0 )
//    {
//       try
//       {
//          result = MAC_Application::make( a_owner, args ) ;
//       }
//       catch( MAC_Exceptions::Error )
//       {
//          result = 0 ;
//       }
//    }

   try
   {
      MAC_CHECK_POST( IMPLIES( result!=0, result->owner() == a_owner ) ) ;
   }
   catch( MAC_Exceptions::Error )
   {
   }
   return( result ) ;
}




//-------------------------------------------------------------------------
MAC_Communicator const*
MAC_Exec:: communicator( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec:: communicator" ) ;
   static MAC_Communicator const* result = 0 ;
   if( result==0 )
   {
      std::string com_name = "MAC_SequentialCommunicator" ;
      MAC_Variable const* var = MAC_Variable::object( "BS_with_MPI" ) ;
      if( execution_context()->has_variable( var ) )
      {
         MAC_ASSERT( execution_context()->value( var )->to_bool() ) ;
         com_name = "EXT_MPIcommunicator" ;
      }
      result = MAC_Communicator::object( com_name ) ;
   }
   MAC_CHECK_POST( result != 0 ) ;
   
   // This post-condition verifies that no static objet needs communicator
   // before MPI one has been created.
   MAC_CHECK_POST( IMPLIES(
	execution_context()->has_variable(MAC_Variable::object("BS_with_MPI") ),
	result->name()=="EXT_MPIcommunicator" ) ) ;
   MAC_CHECK_POST( result->is_under_ownership_of( MAC_Root::object() ) ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
MAC_Context const*
MAC_Exec:: execution_context( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec:: execution_context" ) ;
   
   static bool first = true ;
   if( first )
   {
      EXEC_CONTEXT = MAC_ContextSimple::create( MAC_Root::object() ) ;
      first = false ;
   }
   MAC_Context const* result = EXEC_CONTEXT ;
   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->owner()==MAC_Root::object() ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
void
MAC_Exec:: add_variable_to_execution_context( MAC_Variable const* a_variable,
                                              MAC_Data* a_value )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec:: add_variable_to_execution_context" ) ;
   MAC_CHECK_PRE( a_variable!=0 ) ;
   MAC_CHECK_PRE( a_value!=0 && a_value->owner()==0 ) ;
   MAC_CHECK_PRE( !MAC_Exec::execution_context()->has_variable( a_variable ) ) ;

   // Creation of the context if needed :
   MAC_Context const* ct = MAC_Exec::execution_context() ;
   MAC_ASSERT( ct != 0 ) ;

   a_value->set_owner( EXEC_CONTEXT ) ;
   EXEC_CONTEXT->extend( a_variable, a_value ) ;

   MAC_CHECK_POST( MAC_Exec::execution_context()->has_variable( a_variable ) ) ;
   MAC_CHECK_POST( MAC_Exec::execution_context()->value( a_variable )==a_value ) ;  
}




//------------------------------------------------------------------------
std::ostream&
MAC_Exec::out( void )
//------------------------------------------------------------------------
{
   return( MAC_OUT==0 ? std::cout : *MAC_OUT ) ;
}




//------------------------------------------------------------------------
std::ostream&
MAC_Exec::err( void )
//------------------------------------------------------------------------
{
   return( MAC_ERR==0 ? std::cerr : *MAC_ERR ) ;
}




//------------------------------------------------------------------------
std::istream&
MAC_Exec::in( void )
//------------------------------------------------------------------------
{
   return( MAC_IN==0 ? std::cin : *MAC_IN ) ;
}




//------------------------------------------------------------------------
int
MAC_Exec:: initialize( int argc, char* argv[], stringVector& args )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec:: initialize" ) ;
   MAC_CHECK( args.size()==0 ) ;

   int result = 0 ;

   try
   {
      ONLY_HELP = false ;

      for( size_t i=0 ; i<(size_t)argc ; i++ )
      {
         std::string arg =  argv[i] ;
         
         if( arg == "-no_external_API" )
            EXTERNAL_API = false ;
      }
      if ( EXTERNAL_API )  MAC_ExternalAPI::initialize_all_APIs( argc, argv ) ;

      args.re_initialize( argc ) ;
      for( size_t i=0 ; i<(size_t)argc ; i++ )
      {
         args(i) = argv[i] ;
      }

      // Parse command line arguments specific to MAC.
      parse_arguments( args ) ;

      if( SIGNAL_HANDLING ) MAC_System::exception_trapping() ;
      
      if( VERBOSE )
      {
         APPLI_TIMER = MAC_Timer::create( MAC_Root::object() ) ;
         APPLI_TIMER->start() ;
      }
      if( communicator()->nb_ranks()>1 )
      {
         std::cout << "Process #" << communicator()->rank()
                   << " : launched on " << MAC_System::host_name()
                   << " (" << MAC_System::process_id() << ")." <<  std::endl ;
      }
   }
   catch( MAC_Exceptions::Error )
   {
      result = 1 ;
   }

   return( result ) ;
}




//------------------------------------------------------------------------
void
MAC_Exec:: run_application( MAC_Application* appli  )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec::run_application" ) ;

   if( ONLY_HELP )
   {
      MAC_Exec::print_usage( MAC_Exec::err() ) ;
   }
   else
   {
      try
      {
         if( appli!=0 )
         {
            // program core exection : the tasks specific to the concrete 
            // application are performed
            appli->run() ;
         }
      }
      catch( MAC_Exceptions::Error )
      {
      }
   }
}




//------------------------------------------------------------------------
void
MAC_Exec:: terminate( void )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec::terminate" ) ;

   try
   {
      size_t const rank = communicator()->rank() ;
      size_t const nb_ranks = communicator()->nb_ranks() ;
   
      if( APPLI_TIMER !=0 )
      {
         APPLI_TIMER->stop() ;
         MAC_Exec::out() << "*** Elapsed time in second: " << std::endl ;
         MAC_Exec::out() << "user " << APPLI_TIMER->time() ;
         MAC_Exec::out() << std::endl << std::endl ;
      }
      if( PATTERN_STATUS == MAC_ModuleExplorer::build &&
          MAC_ModulePattern::build_pattern() )
      {
         MAC_ModulePattern::save_pattern() ;
      }
      // termination of all objects belonging to a ownership tree whose
      // root node is not the NULL object
      MAC_Root::cleanup() ;

      if(EXTERNAL_API) MAC_ExternalAPI::terminate_all_APIs() ;
      MAC_Exec::check_for_remaining_objects() ;

      if( nb_ranks>1 ) // communicator() is destroyed...
      {
         std::cout << "Process #" << rank << " : terminated." << std::endl ;
      }
   }
   catch( MAC_Exceptions::Error )
   {
      if( VERBOSE )
      {
         MAC_Exec::out() << "Unable to achieve MAC termination !!!"
                         << std::endl ;
      }   
   }

   out().flush() ;
   err().flush() ;
   if( MAC_OUT!=0 )
   {
      delete MAC_OUT ; MAC_OUT = 0 ;
   }
   if( MAC_ERR!=0 )
   {
      delete MAC_ERR ; MAC_ERR = 0 ;
   }
   if( MAC_IN!=0 )
   {
       delete MAC_IN ; MAC_IN = 0 ;
   }
}




//------------------------------------------------------------------------
void
MAC_Exec:: check_for_remaining_objects( void )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec::check_for_remaining_objects" ) ;
   if( MAC_Object::GetNumberOf_MAC_objects() > 0 )
   {
      MAC_Exec::out() << std::endl << std::endl ;
      MAC_Exec::out() << "after cleaning up, number of remaining P-objects : "
                 << MAC_Object::GetNumberOf_MAC_objects()
                 << std::endl ;
      MAC_Object::TraceRemainingObjects( MAC_Exec::out() ) ;
      set_exit_code( 2 ) ;
   }
}




//----------------------------------------------------------------------
void
MAC_Exec::print( std::string const& data, std::ostream& os ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec:: print" ) ;
   
   os << endl ;
   os << "*** Operating system: " <<
      MAC_System::sysname() << endl  ;
   os << endl ;
   os << "*** Executable: " << name_of_exe() << endl << endl ;

   if( !dynamic_check_list().empty() )
   {
      os << "*** Built-in tests that are dynamically enabled:" << endl ; 
      os << "       " << dynamic_check_list() << endl << endl ;
   }

   os << "*** Data file: " ;
   if( data== "-stdin" )
   {
      os << "read from standard input" ;
   }
   else
   {
      os << data ;
   }
   os << endl << endl;
   
   os << "*** MAC library" << endl ;
   os << "       compiler          : " << MAC_System::compiler_name() << endl ;
   os << "       compilation date  : " << date_of_compilation() << endl ;
   os << "       compilation level : opt" << compilation_level() << endl ;
   os << endl ;
}




//-------------------------------------------------------------------------
void
MAC_Exec:: parse_arguments( stringVector& args )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Exec:: parse_arguments" ) ;

   EXE_FILE = args(0) ;
   args.remove_at( 0 ) ;
   size_t iArgs=0 ;
   while( iArgs<args.size() )
   {
      std::string const& argument = args(iArgs) ;
      bool done = true ;
      // Here take place loop for argument scanning
      if( argument == "-Cpost" )
      {
         MAC_Assertion::add_handled_check( MAC_Assertion::Postcondition ) ;
      }
      else if( argument == "-Cobjects" )
      {
         MAC_Assertion::add_handled_check( MAC_Assertion::Objects ) ;
      }
      else if( argument == "-catch" )
      {
         MAC_Assertion::add_handled_check( MAC_Assertion::Objects ) ;
         args.remove_at( iArgs ) ;
         std::istringstream is( args(iArgs) ) ;
         size_t r ;
         is >> r ;
         
         MAC_Object::catch_object_by_rank( r ) ;
      }
      else if( argument == "-build_pattern" )
      {
         args.remove_at( iArgs ) ;
         MAC_ModulePattern::build_pattern_base( args(iArgs).c_str() ) ;
         PATTERN_STATUS = MAC_ModuleExplorer::build ;
      }
      else if( argument == "-check_pattern" )
      {
         args.remove_at( iArgs ) ;
         MAC_ModulePattern::build_pattern_base( args(iArgs).c_str() ) ;
         PATTERN_STATUS = MAC_ModuleExplorer::verify ;
      }
      else if( argument == "-Call" )
      {
         MAC_Assertion::add_handled_check( MAC_Assertion::Postcondition ) ;
         MAC_Assertion::add_handled_check( MAC_Assertion::Invariant ) ;
         MAC_Assertion::add_handled_check( MAC_Assertion::Check ) ;         
      }
      else if( argument == "-v" )
      {
         VERBOSE = true ;
      }
      else if( argument == "-o" )
      {
         args.remove_at( iArgs ) ;
         std::string outfile = args(iArgs) ;
         if( communicator()->nb_ranks()>1 )
         {
            std::ostringstream os ;
            os << "." <<  (int)communicator()->rank() ;
            outfile += os.str() ;
         }
         MAC_OUT = new std::ofstream( outfile.c_str() ) ;
         if( !(*MAC_OUT) )
         {
            MAC_Error::object()->raise_plain(
               "Unable to open standard output file "+outfile ) ;
         }
      }
      else if( argument == "-H" )
      {
         ONLY_HELP = true ;
      }
      else if( argument == "-no_signal_handling" )
      {
         SIGNAL_HANDLING = false ;
      }
      else if( argument == "-no_external_API" )
      {
         // done
      }
      else if( argument == "-notify_parallel" )
      {  
         MAC_Context const* ct = MAC_Exec::execution_context() ;
         if( !ct->has_variable( MAC_Variable::object( "BS_with_MPI" ) ) )
         {
            MAC_Error::object()->raise_plain(
               "Parallel execution not allowed: no MPI tools defined" ) ;
         }
      }
      else
      {
         done = false ;
      }
      if( done )
      {
         args.remove_at( iArgs ) ;
      }
      else
      {
         iArgs++ ;
      }
   }
   if( args.size() == 0 )
   {
      ONLY_HELP = true ;
   }
}




// Create module from args NOT SUPPORTED for now
// //------------------------------------------------------------------------
// MAC_Module* 
// MAC_Exec:: create_module_from_args( MAC_Object* a_owner,
//                                     std::string const& application_name,
//                                     stringVector& args )
// //------------------------------------------------------------------------
// {
//    MAC_LABEL( "MAC_Exec:: create_module_from_args" ) ;
//    
//    std::stringstream str ;
//    str << "MODULE MAC_Application" << std::endl ;
//    str << "concrete_name=\"" << application_name << "\"" << std::endl ;
//    size_t i=0 ;
//    
//    while( i<args.size() )
//    {
//       std::string const& arg = args(i) ;
//       size_t idx = arg.find_first_of( "=" ) ;
//       if( idx<arg.length() )
//       {
//          str << arg << std::endl;
//          args.remove_at(i) ;
//       }
//       else
//       {
//          i++ ;
//       }
//    }
//    str << "END MODULE MAC_Application" << std::endl ;
//    
//    MAC_Module* result = MAC_Module::create( a_owner, "MAIN", str ) ;
//    
//    return( result ) ;
// }




//------------------------------------------------------------------------
void
MAC_Exec:: print_usage( std::ostream& os )
//------------------------------------------------------------------------
{
   os << "Usage: " << endl ;
   os << "   [executable] [-H] [-V] [-Cobjects] [-catch <obj_nb>] "
      << "[-Cpost|-Call] [ -build_pattern <file> ] ( data_file | -stdin )" 
      << endl << endl ;
   os << "Options:" << endl ;
   os << "   [-h] : Display help and exit" << endl ;
   os << "   [-v] : verbosity"  << endl ;
   os << "   [-Cobjects]    : orphean objects tracking" << endl ;
   os << "   [-catch <obj_nb>]: object creation catching" << endl ;
   os << "   [-Cpost|-Call] : assertion level modification" << endl ;
   os << "   [ -build_pattern <file> ] : extend given pattern" << endl ;
   os << "   [ -check_pattern <file> ] : check from given pattern" << endl ;
   os << "   [ -no_signal_handling ] : doesn't provide signal handling facility"
   	 << endl ;
   os << "   [ -no_external_API ] : doesn't provide extra libraries management"
   	 << endl ;
   os << endl << "Arguments:" << endl ;
   os << "   ( data_file | -stdin ) : data file or standard input"<< endl ;
}
