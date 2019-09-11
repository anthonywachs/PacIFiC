#ifndef MAC_EXTRACTION_EXP_HH
#define MAC_EXTRACTION_EXP_HH

#include <MAC_TransferExp.hh>
#include <string>

class MAC_ModuleExplorer ;

/*
Expressions extracting data from an attached data structure that has
to be primarily specified by calling `::initialize'( `mod' ).
mod must not be modified until having called `::reset' method.

Before calling `::initialize', none of the expressions implemented here
can be used.

---
name     : has_data
argument : String
type     : Bool

The returned value is true if there exists, in the attached data structure,
an entry whose keyword match the argument.

Example:
  if( `has_data'( "/TOTO/my_entry" ) )
  MODULE titi
     ...
  END MODULE titi

---
name      : extracted_data
arguments : String, optional second argument
type      : everything

The returned value is the result of the evaluation of the data
whose keyword matches the first argument (in the attached data structure),
if such an entry exists (in the attached data structure). If not, the
second argument can be used to define a default value.

Example 1:

   toto = `extracted_data'( "/TOTO/my_entry" )
   
Example 2:

   toto = `extracted_data'( "/TOTO/my_entry", 3. )
   (toto is "/TOTO/my_entry" in `mod' if any, 3. elsewhere)
      
   titi = `extracted_data'( "/TOTO/my_entry", < "a" "b" > )
   (titi is "/TOTO/my_entry" in `mod' if any, < "a" "b" > elsewhere)

Example 3:

   For implementation reasons, the second argument becomes mandatory in
   "if" constructions, even when it is not used.

   if( has_data( "/TOTO/my_entry" ) )
   MODULE titi
      toto = `extracted_data'( "/TOTO/my_entry", 0. )
   END MODULE titi
   
---
name     : has_module
argument : String
type     : Bool

The returned value is true if there exists, in the attached data structure,
a module whose keyword match the argument.

Example:
  if( `has_module'( "/TOTO" ) )
  MODULE titi
     toto = `extracted_data'( "/TOTO/my_entry", 3. )
     ...
  END MODULE titi
  
---
name     : extracted_module
argument : String, String, optional third argument
type     : String

The returned value is the name of a temporary file containing the
MODULE called according to the first argument in the attached data structure.
The extracted module is renamed as the second argument of the function.

Example1:

    MODULE titi
       #include( `extracted_module'( "/TOTO", "TITI" ) )
    END MODULE titi

    The module of name "TOTO" is extracted from the database,
    and included as module "TITI".

Example2:

   For implementation reasons, a special syntax (with a dummy third argument)
   is required in "if" constructions.

   if( has_module( "/TOTO/module" ) )
   MODULE titi
      #include( `extracted_module'( "/TOTO/module", "module", "" ) )
   END MODULE titi
   
PUBLISHED
*/

class MAC_ExtractionExp : public MAC_TransferExp
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization
      
      static void initialize( MAC_Module const* mod ) ;
      static void reset( void ) ;

      static bool is_initialized( void ) ;
      static MAC_Module const* data_base( void ) ;
      
   //-- Context
      
      virtual void declare( MAC_List* lst ) const ;
      virtual bool context_has_required_variables( 
                                               MAC_Context const* ct ) const ;
      
   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual bool value_can_be_evaluated( MAC_Context const* ct ) const ;
      virtual stringVector const& undefined_variables(
                                           MAC_Context const* ct ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------
      
      MAC_ExtractionExp( void ) ;
     ~MAC_ExtractionExp( void ) ;
      MAC_ExtractionExp( MAC_ExtractionExp const& other ) ;
      MAC_ExtractionExp& operator=( MAC_ExtractionExp const& other ) ;
      
      enum ExtractionExp{ has_data, ext_data, has_mod, ext_mod } ;
      
      MAC_ExtractionExp( MAC_Object* a_owner,
                         ExtractionExp op,
                         std::string const& a_name,
                         MAC_Sequence const* argument_list ) ;

   //-- Plug in

      MAC_ExtractionExp( std::string const& a_name, ExtractionExp op ) ;

      virtual MAC_ExtractionExp* create_replica( 
                      MAC_Object* a_owner,
                      MAC_Sequence const* argument_list ) const ;

      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Transfer implementation
      
      virtual MAC_Data const* data( MAC_Context const* ct ) const ;
      
   //-- Private methods

      static std::string const& temporary_file( void ) ;

      static std::string const& data_name( std::string const& exp_name,
                                           MAC_Data const* d ) ;
      
      void extract_module( std::string const& file_name,
                           std::string const& d_name,
                           std::string const& m_name ) const ;

   //-- Class attributes
      
      static MAC_Module const* DB_MOD ;

      static MAC_ExtractionExp* PROTO_HAS_DATA ;
      static MAC_ExtractionExp* PROTO_DATA ;
      static MAC_ExtractionExp* PROTO_HAS_MOD ;
      static MAC_ExtractionExp* PROTO_EXTRACTED_MODULE ;
 
   //-- Attributes

      ExtractionExp const OP ;
      std::string TEMP_FILE_NAME ;
      MAC_Data const* SRC ;
} ;

#endif
