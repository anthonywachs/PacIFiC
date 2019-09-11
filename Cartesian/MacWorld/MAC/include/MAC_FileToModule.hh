#ifndef MAC_FILE_TO_MODULE_HH
#define MAC_FILE_TO_MODULE_HH

#include <MAC_Object.hh>

class MAC_ModuleExplorer ;
class MAC_ObjectRegister ;

/*
FRAMEWORK INSTANTIATION
PUBLISHED
*/

class MAC_FileToModule : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_FileToModule const* object( std::string const& format ) ;
      
      static bool has( std::string const& format ) ;
      
      static void find_file_format( std::string const& a_file_name,
                                    std::string& a_format ) ;
      
      static std::string const& list_of_formats( void ) ;

   //-- Characteristics

      std::string const& format( void ) const ;
      
      std::string const& default_motif( void ) const ;

   //-- Module building
      
      virtual MAC_Module* create_from_file( 
                                 MAC_Object* a_owner,
                                 std::string const& module_name,
                                 std::string const& file_name ) const = 0 ;
      
   protected: //--------------------------------------------------------

   //-- Plug in

      virtual ~MAC_FileToModule( void ) ;

      MAC_FileToModule( std::string const& a_format,
                        std::string const& a_default_motif ) ;

   //-- Preconditions, Postconditions, Invariant

      bool create_from_file_PRE( MAC_Object* a_owner,
                                 std::string const& module_name,
                                 std::string const& file_name ) const ;
      
      bool create_from_file_POST( MAC_Module const* result,
                                  MAC_Object* a_owner,
                                  std::string const& module_name,
                                  std::string const& file_name ) const ;
      
   private: //----------------------------------------------------------

      MAC_FileToModule( void ) ;
      MAC_FileToModule( MAC_FileToModule const& other ) ;
      MAC_FileToModule& operator=( MAC_FileToModule const& other ) ;

      static MAC_ObjectRegister* plugins_map( void ) ;
      
      static std::string& formats( void ) ;
      
   //-- Attributes

      std::string MY_FORMAT ;
      std::string MY_MOTIF ;
} ;

#endif
