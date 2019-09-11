#ifndef MAC_MODULE_PATTERN_HH
#define MAC_MODULE_PATTERN_HH

#include <MAC_Object.hh>
#include <MAC_Data.hh>

#include <string>

class MAC_Module ;
class stringVector ;

/*
  Provides functionalities to analyse and/or verify data deck.
*/

class MAC_ModulePattern : public MAC_Object
{
   public: //----------------------------------------------------------------

  //-- Instance delivery and initialization

      // Create and return an instance attached to base pattern `a_name'.
      static MAC_ModulePattern* create( MAC_Object* a_owner,
                                        MAC_Module const* way,
                                        std::string const& a_name ) ;
      // Duplicate `self'.
      virtual MAC_ModulePattern* create_clone( MAC_Object* a_owner ) const ;

      // Search for description of class `class_name'.
      // Return NULL if no description is found.
      static MAC_ModulePattern* create_pattern_description(
                                        MAC_Object* a_owner,
                                        std::string const& class_name ) ;
      
  //-- Status report

      // Does current pattern is an indirection on sub-pattern ?
      static bool is_mutable( MAC_Module const* module ) ;
      
      
      // refering pattern description ( follow or the class name )
      static std::string const& indirection_to( MAC_Module const* module ) ;
      
      // does `module' refer to a module description ?
      static bool is_module_description( MAC_Module const* module ) ;
      
      // does `module' refer to an entry description ?
      static bool is_entry_description( MAC_Module const* module ) ;
      
      // does `module' refer to a variable description ?
      static bool is_variable_description( MAC_Module const* module ) ;
      
      // does `module' refer to a condition description ?
      static bool is_condition_description( MAC_Module const* module ) ;
      
      // valid choices for indirection
      static stringVector const& valid_indirections( MAC_Module const* module ) ;
      // needed module for current pattern
      stringVector const& mandatory_modules( MAC_Module const* way, bool first=true ) const ;
      
      // needed entries for current pattern
      stringVector const& mandatory_entries( MAC_Module const* way, bool first=true ) const ;
      
      // type of data
      std::string const& type_of_entry( std::string const& a_name ) const ;
      
      // name of `self'
      std::string const& name( void ) const ;
      
      // validity of `checked' with current pattern
      MAC_Module* validity( MAC_Object* a_owner,
                            MAC_Module const* checked,
                            bool recurse,
                            MAC_Module* result=0 ) const ;
      
      stringVector const&  provided_variables( void ) const ;
      
      static std::string const& variable_access( MAC_Module const* module,
                                                 MAC_Module const* way )  ;
      
   //-- Modifier
      
      enum Access { mandatory, optional, generic, unspecified } ;

      // Add module `path_and_name' to current pattern with accessibility
      //  `access' following `way'.
      void add_pattern( std::string const& path_and_name,
                        MAC_Module const* way,
                        Access access ) const ;
      
      // Add generic entry in current pattern following `way'.
      void add_generic_keyword( std::string const& a_keyword,
                                MAC_Module const* way,
                                MAC_Data::Type type ) ;
      
      // Add entry `a_keyword' in current pattern following `way'.
      void add_entry( std::string const& a_keyword,
                      MAC_Module const* way,
                      MAC_Data::Type type,
                      Access acc ) ;
      
      // Add provided entry `a_keyword' in current pattern following `way'.
      void add_provided_variable( std::string const& a_keyword,
                                  MAC_Module const* way,
                                  MAC_Data::Type type ) ;
      
      // Create subpattern of `self' following `way' at `path_and_name' position.
      MAC_ModulePattern* create_subpattern(
                                  MAC_Object* a_owner,
                                  MAC_Module const* way,
                                  std::string const& path_and_name ) const ;
      
      // pattern on generic module of current pattern following `way'
      MAC_ModulePattern* generic_pattern( MAC_Object* a_owner,
                                          MAC_Module const* way ) const ;
      
      static void expand_then_simplify( void ) ;
      
      void attach_help_data( std::string const& keyword,
                             std::string const& help ) ;
      
      void attach_default_data( std::string const& keyword,
                                std::string const& value ) ;
      
      void attach_verify_data( std::string const& keyword,
                               std::string const& expression ) ;
      
      void attach_list_of_valid_choices( std::string const& keyword,
                                         stringVector const& valid_choices ) ;
      
      void attach_list_of_dynamic_choices( std::string const& keyword,
                                           std::string const& regexp,
                                           std::string const& where ) ;
      
      void attach_file_extension( std::string const& filename,
                                  std::string const& mode ) ;
      
      // Is `a_name' an allowed module of current pattern ?
      MAC_Module* allowed_module( std::string const& a_name,
                                  MAC_Module const* way ) const ;
      
      MAC_Module* generic_module( MAC_Module const* way ) const ;
      
      MAC_Module* generic_entry( MAC_Module const* way ) const ;
      
      // Is `a_name' an allowed entry of current pattern ?
      MAC_Module const*  allowed_entry( std::string const& a_name,
                                        MAC_Module const* way  ) const ;
      
   //-- Module pattern base management

      // Attach module pattern base to file `filename' and start
      //  pattern recognition process.
      // If the file exists, it is read and inserted in base.
      static void build_pattern_base( std::string const& filename ) ;
      
      // Attach module pattern base to file `filename'.
      static void open_pattern_base( std::string const& filename ) ;
      
      // Close module pattern base.
      static void close_pattern_base( void ) ;
      
      // Save base to file.
      static void save_pattern( void ) ;
      
      // Has pattern to be built ?
      static bool build_pattern( void ) ;
      
      // Has pattern base ?
      static bool has_pattern( void ) ;
     
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------
      
      MAC_ModulePattern( void ) ;
     ~MAC_ModulePattern( void ) ;
      MAC_ModulePattern( MAC_ModulePattern const& other ) ;
      MAC_ModulePattern& operator=( MAC_ModulePattern const& other ) ;

      MAC_ModulePattern( MAC_Object* a_owner, MAC_Module* mm, MAC_Module const* a_father=0 ) ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Internal
      
      static void set_property( MAC_Module* mod,
                                std::string const& name,
                                std::string const& property ) ;
      
      static std::string const& property( MAC_Module const* mod,
                                          std::string const& name ) ;
      
      static bool is_module_matching( std::string const& module_name,
                                      std::string const& pattern_name ) ;
      
      void add_pattern_simple( std::string const& a_name,
                               MAC_Module const* way,
                               Access access ) const ;
      
      void add_entry_simple( std::string const& a_keyword,
                             MAC_Data::Type type,
                             Access access  ) ;
      
      static MAC_List* list_of_selected_conditional_modules(
                             MAC_Object * a_owner,
                             MAC_Module const* module,
                             MAC_Module const* way ) ;
      
      static MAC_List* list_of_conditional_modules(
                             MAC_Object * a_owner,
                             MAC_Module const* module ) ;
      
      static void check_name_validity( std::string const& a_name ) ;
      
      static MAC_Module*  module_of_pattern(
                             std::string const& class_name,
                             std::string const& instance_name,
                             bool build_if_not_exist ) ;
      
      static stringVector const& access_name( void ) ;
      
      MAC_Module* find_subpattern(
                             MAC_Module const* way,
                             std::string const& a_name,
                             MAC_Module *& a_father ) const ;
      
      MAC_Module* find_subpattern_simple(
                             MAC_Module const* way,
                             std::string const& a_name,
                             MAC_Module *& a_father ) const ;
      
      static void split( std::string const& a_name,
                         std::string& first_dir,
                         std::string& last ) ;
      
      static std::string class_part( std::string const& class_name ) ;
      
      void expand_then_simplify_one( void ) ;
      
      // Is value of `a_name' as an entry of `way' a valid
      //  with respect to `model' ?
      void check_allowed_value_for_entry(
                             MAC_Module* result,
                             MAC_Module const* model,
                             std::string const& a_name,
                             MAC_Module const* way ) const ;
      
      std::string inferred_generic_name( std::string const& a_name ) const ;
      
      static MAC_Module* child( MAC_Module const* module,
                                std::string const& name ) ;
      
      static void extend_entries( MAC_Module * module,
                                  std::string const entrie_name,
                                  std::string const& value ) ;
      
      static void complete_indirection_choices( MAC_Module * root ) ;

      static void check_evaluable( bool pattern_mod,
                                   MAC_Module const* mod,
                                   std::string const& keyword,
                                   MAC_Data::Type data_type ) ;
      
      static void check_in_spec( MAC_Module* result,
                                 MAC_Module const* model,
                                 std::string const& a_name,
                                 MAC_Module const* way ) ;
      
      static void check_vector_in_spec( MAC_Module* result,
                                        MAC_Module const* model,
                                        std::string const& a_name,
                                        MAC_Module const* way ) ;
      static void check_unique_spec( MAC_Module* result,
                                     MAC_Module const* model,
                                     std::string const& a_name,
                                     MAC_Module const* way ) ;
      static void merge( MAC_Module* current, MAC_Module const* to_add ) ;
      
   //-- Static attribute
      
      static MAC_Module* MP_File ;
      static MAC_Module* MP_Base ;
      static MAC_Module* MP_Pattern ;
      static std::string pattern_filename ;
      static bool HAS ;
      static bool BUILD ;

   //-- Attribute
      
      MAC_Module* const mod ;
      MAC_Module const* father ;
} ;

#endif
