#ifndef MAC_ERROR_HH
#define MAC_ERROR_HH

#include <MAC_Object.hh>
#include <MAC_Data.hh>

#include <iosfwd>
#include <string>

class MAC_Map ;
class MAC_ModuleExplorer ;
class stringVector ;

//---------------------------------------------------------------------------
//   Facilities for adapting the error handling mechanism
//---------------------------------------------------------------------------

class MAC_Error
{
   public: //----------------------------------------------------------------

      static MAC_Error* object( void ) ;

      void raise_plain( std::string message ) ;

      void raise_internal( std::string message ) ;

      void raise_read_syntax_error( std::string const& file,
                                    int line,
                                    std::string const& last,
                                    std::string const& nature ) ;

      void display_info( std::string message ) ;
      
   //-- Incomplete implementation

      void raise_not_tested( std::string file, int line, std::string test ) ;

      void raise_not_implemented( MAC_Object const* oo,
                                  std::string method ) ;

   //-- Simulation of virtual methods with respect to many objects

      void raise_bad_types( MAC_Object const* oo,
                            std::string method,
                            MAC_Object const* argument ) ;

      void raise_bad_types( MAC_Object const* oo,
                            std::string method,
                            MAC_Object const* argument_1, 
                            MAC_Object const* argument_2 ) ;

   //-- Registration

      void raise_missing_keyword( MAC_ModuleExplorer const* exp,
                                  std::string keyword ) ;

      void raise_missing_module( MAC_ModuleExplorer const* exp,
                                 std::string path_and_name ) ;
      
      void raise_module_error( MAC_ModuleExplorer const* exp,
                               std::string error ) ;

      void raise_bad_data_type( MAC_ModuleExplorer const* exp,
                                std::string keyword,
				MAC_Data::Type query_kind ) ;

      void raise_bad_file( MAC_ModuleExplorer const* exp,
                           std::string const& filename,
                           std::string const& access ) ;

      void raise_not_evaluable( MAC_ModuleExplorer const* exp,
                                std::string keyword,
                                stringVector const& undefined_variables ) ;

      void raise_bad_data_value( MAC_ModuleExplorer const* exp,
                                 std::string keyword,
				 std::string allowed_values ) ;
      
      void raise_data_error( MAC_ModuleExplorer const* exp,
                             std::string keyword,
                             std::string error ) ;

   //-- Contract

      void raise_precondition_violation( std::string test ) ;

      void raise_postcondition_violation( std::string file,
                                          int line,
                                          std::string test ) ;
      void raise_invariant_violation( std::string file,
                                      int line,
                                      std::string test ) ;
      void raise_assertion_violation( std::string file,
                                      int line,
                                      std::string test ) ;

      void raise_bad_object( std::string message,
                             MAC_Object const* a_object ) ;
      
   //-- File manipulation

      void raise_file_handling( std::string file,
                                std::string operation ) ;

   //-- Warnings

      void trace( std::string const& event ) ;
      
      void display_new_syntax( MAC_ModuleExplorer const* older_module,
                               MAC_ModuleExplorer const* result_module ) const ;
      
      static void notify( MAC_Module* list,
                          std::string const& message,
                          std::string const& where,
                          stringVector const& choices )  ;
      
      static MAC_Module* notify( MAC_Module* list,
                                  std::string const& message,
                                  std::string const& where )  ;
      
      void display_data_checking( MAC_ModuleExplorer* report ) const ;
      
   //-- Abnormal termination
      
      // Exits on error and return exit_status to shell.
      static void exit( size_t status = 1 ) ;
      
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

     ~MAC_Error( void ) ;
      MAC_Error( MAC_Error const& other ) ;
      MAC_Error const& operator=( MAC_Error const& other ) ;

      MAC_Error( void ) ;

      void print_invocated_methods( void ) ;
      void print_client_and_server( void ) ;
      
      void begin( void ) ;
      void end( void ) ;
      std::string const& hline( void ) const ;
      
      static void exit_with_internal_error( void ) ;
      
      //---------------------------------------------------------------------
      //   ATTRIBUTES
      //---------------------------------------------------------------------
      std::ostream& os ;

} ;

#endif
