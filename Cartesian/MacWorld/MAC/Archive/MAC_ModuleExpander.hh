#ifndef MAC_MODULE_EXPANDER_HH
#define MAC_MODULE_EXPANDER_HH

#include <MAC_Application.hh>
#include <string>

class MAC_Module ;
class MAC_ModuleExplorer ;

/*
Builders of data structures expanded from two parts : the first one is
a skeleton data structure with missing pieces, the second one is a
complementary data structure providing these missing pieces.

The missing pieces of the skeleton are entries whose data expressions
implemented by the `MAC_ExtractionExp::' class.

PUBLISHED
*/

class MAC_ModuleExpander : public MAC_Application
{
   public: //-----------------------------------------------------------

      static MAC_ModuleExpander* create( MAC_Object* a_owner,
                                         MAC_ModuleExplorer const* exp ) ;
      
   //-- Program core execution

      virtual void run( void ) ;

   //-- Expander

      static MAC_Module const* create_expanded_module(
                          MAC_Object* a_owner,
                          MAC_Module const* input_mod,
                          std::string const& skeleton_file_name ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~MAC_ModuleExpander( void ) ; 
      MAC_ModuleExpander( MAC_ModuleExpander const& other ) ;
      MAC_ModuleExpander& operator=( MAC_ModuleExpander const& other ) ;

      MAC_ModuleExpander( MAC_Object* a_owner,
                          std::string const& skeleton_file,
                          std::string const& input_file,
                          std::string const& expanded_file,
                          std::string const& submod_name ) ;
      
   //-- Plug in
      
      MAC_ModuleExpander( void ) ;
      
      virtual MAC_ModuleExpander* create_replica( 
                            MAC_Object* a_owner,
                            MAC_ModuleExplorer const* exp ) const ;

      virtual MAC_ModuleExpander* create_replica_from_args(
                            MAC_Object* a_owner,
                            stringVector& args ) const ;
    //-- Command line

      virtual void print_usage( void ) const ;
      virtual void print_operands( void ) const ;

   //-- Expander

      static MAC_Module* create_skeleton_module(
                   MAC_Object* a_owner, std::string const& file_name ) ;

   //-- Class attribute
      
      static MAC_ModuleExpander const* PROTOTYPE ;
      
   //-- Attribute
      
      std::string const BASE ;
      std::string const OUTPUT ;
      std::string const INPUT ;
      std::string SUBMOD ;

} ;

#endif



