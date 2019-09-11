#ifndef MAC_EXEC_HH
#define MAC_EXEC_HH

#include <cstddef>
#include <iosfwd>

#include <MAC_ModuleExplorer.hh>

class stringVector ;

class MAC_Application ;
class MAC_Communicator ;
class MAC_Context ;
class MAC_ContextSimple ;
class MAC_Module ;
class MAC_Object ;
class MAC_Timer ;

class MAC_Exec
{
   public: //-----------------------------------------------------------------
      
   //-- MAC library self-description
      
      // Compilation level of MAC library.
      static size_t compilation_level( void ) ;

      // Date when was compiled MAC library currently used.
      static std::string const& date_of_compilation( void ) ;
      
      // Name of executable.
      static std::string const& name_of_exe( void ) ;

      // List of dynamic checks activated during execution.
      static std::string dynamic_check_list( void ) ;

   //-- MAC core execution

      // Initialize platform from argument list `args', 
      // and return 0 if successful.
      // Arguments not recognized by MAC are set in args.
      static int initialize( int argc, char * argv[],
                              stringVector& args ) ;

      // Terminate MAC session.
      static void terminate( void ) ;

      // Return to system exit code.
      // >=1 : exit on user-defined user error
      // 1 : exit on default user error
      // 0 : normal exit
      // -1 : exit on internal error
      // -2 : exit on internal error du to remaining objects after termination.
      static int exit_code( void ) ;
      
      // Modify exit code.
      static void set_exit_code( int status ) ;
      
   //-- `::MAC_Application' factory
      
      // Create application module from file whose name is given in `args'
      //  remaining argument list.
      static MAC_Module const* create_module_with_data( MAC_Object* a_owner,
                                                        stringVector& args ) ;

      // Create application from module `mod'
      static MAC_Application* create_application( MAC_Object* a_owner,
                                                  MAC_Module const* appli_mod,
                                                  stringVector& args ) ;
      
      // Launch application `appli'.
      static void run_application( MAC_Application* appli  ) ;

   //-- Parallel tools

      static MAC_Communicator const* communicator( void ) ;

   //-- Context of execution

      static MAC_Context const* execution_context( void ) ;
      static void add_variable_to_execution_context(
                          MAC_Variable const* a_variable, MAC_Data* a_value ) ;
      
   //-- Input - Output
      
      // Pelicans standard output stream.
      static std::ostream& out( void ) ;
      
      // Pelicans standard error stream.
      static std::ostream& err( void ) ;
      
      // Pelicans standard input stream.
      static std::istream& in( void ) ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      MAC_Exec( void ) ;
     ~MAC_Exec( void ) ;
      MAC_Exec( MAC_Exec const& other ) ;
      MAC_Exec& operator=( MAC_Exec const& other ) ;

      // Check for remaining objects.
      static void check_for_remaining_objects( void ) ;

      static void print( std::string const& data,
                         std::ostream& os ) ;
      
      static void parse_arguments( stringVector& args ) ;

//       Create module from args NOT SUPPORTED for now
//       static MAC_Module* create_module_from_args(
//                                      MAC_Object* a_owner,
//                                      std::string const& application_name,
//                                      stringVector& args ) ;
      
      static void print_usage( std::ostream& os ) ;
      
      static bool VERBOSE ;
      static bool ONLY_HELP ;
      static bool SIGNAL_HANDLING ;
      static bool EXTERNAL_API ;
      
      static MAC_Timer* APPLI_TIMER ;
      
      static std::string EXE_FILE ;

      static int EXIT_STATUS ;
      
      static MAC_ModuleExplorer::PatternStatus PATTERN_STATUS ;

      static MAC_ContextSimple* EXEC_CONTEXT ;

      static std::istream* MAC_IN ;
      static std::ostream* MAC_OUT ;
      static std::ostream* MAC_ERR ;
      
} ;

#endif
