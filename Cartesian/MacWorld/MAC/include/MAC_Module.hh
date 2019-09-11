#ifndef MAC_MODULE_HH
#define MAC_MODULE_HH

#include <MAC_Object.hh>

#include <string>
#include <iosfwd>

class MAC_KeywordDataPair ;
class MAC_KeywordDataIterator ;
class MAC_List ;
class MAC_ModuleIterator ;
class MAC_ModuleExplorer ;
class MAC_Context ;
class MAC_ContextPair ;
class MAC_ContextSimple ;
class MAC_Data ;
class MAC_DataWithContext ;
class MAC_String ;
class MAC_Variable ;

/* 
Nodes of the MAC Hierarchical Data System.

The MAC Hierarchical Data System is described in `MAC_ModuleExplorer::'.
Note on context : all childs module share context from their father with
 their own context.
*/

class MAC_Module : public MAC_Object
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization
      
      // Creates and return an instance containing neither modules nor entries.
      static MAC_Module* create( MAC_Object* a_owner,
                                 std::string const& a_name ) ;
      
      // Create and return an instance being the root a the module hierarchy
      // stored in `file_name'.
      static MAC_Module* create( MAC_Object* a_owner,
                                 std::string const& a_name,
                                 std::string const& file_name,
                                 MAC_Context const* ct = 0 ) ;
      
      static MAC_Module* create( MAC_Object* a_owner,
                                 std::string const& a_name,
                                 std::istream& input_stream,
                                 MAC_Context const* ct = 0 ) ;

      // Create and return a clone of `self': modules and entries are clones,
      // all variables of `::father'() context needed is retrieved.
      virtual MAC_Module* create_clone( MAC_Object* a_owner ) const ;

      // Create and return a copy of `self': modules and entries are clones,
      // `::father'() context is lost.
      MAC_Module* create_copy( MAC_Object* a_owner ) const ;
      
      static MAC_Module* create_as_difference( MAC_Object* a_owner,
                                               std::string const& a_name,
                                               MAC_Module const* m1,
                                               MAC_Module const* m2,
                                               MAC_ModuleExplorer const* exp ) ;

   //-- Modifiers

      void modify_module_name( std::string const& a_name ) ;
      
      // Make `self' contain the keyword-data pair defined by `keyword'
      // and `data' as an entry. If `self' already contains an entry
      // with the same keyword, a fatal error is raised.
      void add_entry( std::string const& keyword, MAC_Data const* data ) ;

      // Make `self' contain `a_module'. 
      // If `self' already contains a module with the same name as `a_module', 
      // the two modules are merged (if it happens that two leaves have the 
      // same keyword, a fatal error is raised).
      void add_module( MAC_Module* a_module ) ;
      
      // Add all the modules and entries of `a_module' to `self'.
      // If it happens that two leaves have the same keyword, the value of
      // `a_module' is taken).
      void merge_module( MAC_Module* a_module ) ;

      // Remove the module identified by the `path_and_name' in the 
      // module hierarchy of root `self'.
      void remove_module( std::string const& path_and_name ) ;

      // Remove the leaf identified by the `path_and_name' in the 
      // module hierarchy of root `self'.
      void remove_entry( std::string const& path_and_name ) ;

   //-- Access
      
      // name of `self'
      std::string const& name( void ) const ;

      // absolute path-name of `self'
      std::string const& absolute_path_name( void ) const ;

      // Is `self' empty ?
      bool is_empty( void ) const ;

      // Is there a module in the 
      // module hierarchy of root `self' ? 
      bool has_module( void ) const ;

      // Is there a module identified by the path-name `path_and_name' in the 
      // module hierarchy of root `self' ? 
      bool has_module( std::string const& path_and_name ) const ;
      
      // module containing this if any (NULL otherwise)
      MAC_Module const* father( void ) const ;
      
      // first module whose path-name is `path_and_name' in the module
      // hierarchy of root `self' (if none, a fatal error is raised)
      MAC_Module* module( std::string const& path_and_name ) const ;
      
      // Is there an entry in the 
      // module hierarchy of root `self' ? 
      bool has_entry( void ) const ;

      // Is there an entry identified by the path-name `path_and_name' in the 
      // module hierarchy of root `self' ? 
      bool has_entry( std::string const& path_and_name ) const ;
      
      // data of the first entry identified by the path-name `path_and_name' 
      // in the module hierarchy of root `self' 
      // (if none, a fatal error is raised)
      MAC_Data const* data_of_entry( std::string const& path_and_name ) const ;

      // Make `new_data' be the data of the first  entry identified by the 
      // path-name `path_and_name' in the module hierarchy of root `self'
      // (if none, a fatal error is raised).
      void replace_data_of_entry( std::string const& path_and_name,
                                  MAC_Data const* new_data ) const ;

      // Create and return a list of entries such that the regular expression
      // `regexp' matches their path name relative to `self'.
      MAC_List* create_data_selection( MAC_Object* a_owner,
                                       std::string const& regexp,
                                       MAC_List* result,
                                       std::string const& where ) const ;
   //-- Context management
      
      // variable context associated to `self'
      MAC_Context const* context( void ) const ;
      
      // Add to `::context()' variable `variable' with default
      // `value'.
      void add_variable( MAC_Variable const* variable,
                         MAC_Data* value ) ;
      
      // Modify `::context()' variable `variable' with 
      // `value'.
      void modify_variable( MAC_Variable const* variable,
                            MAC_Data* value ) ;
      
   //-- Input - Output
      
      // Append the Module Hierarchy of root `self' to file whose name is
      // `file' with either "text" or "hybrid" `format'.
      // In hybrid format, whenever possible, each entry is replaced
      // by reference to an associated binary file named `file'.bin
      // The hybrid format is used :
      //   - to reduce file size ;
      //   - to ensure exact float recovering.
      // (binary files generated with hybrid format might not be readable on 
      // system different from these used to produce them).
      void write( std::string const& file,
                  std::string const& format ) const ;
      
      // IMPLEMENTATION : write the Module Hierarchy of root `self' with the
      // format of a MAC data file.
      virtual void print( std::ostream& os, size_t indent_width ) const ;

      virtual void display_info( std::ostream& os, size_t indent_width ) const ;

      // name of module beeing parsed: usefull in error treatment
      static std::string const& current_parsed_module_path_name( void ) ;

      // string representation of owned entry given by its path
      std::string data_as_string( std::string const& path_and_name,
                                  MAC_Context const* ct,
                                  bool & failed ) const ;
      
   //-- Iterators
      
      // Create and return an iterator on modules contained in `self'.
      MAC_ModuleIterator* create_module_iterator( MAC_Object* a_owner ) const ;

      // Create and return an iterator on entries contained in `self'.
      MAC_KeywordDataIterator* create_entry_iterator( MAC_Object* a_owner ) const ;

   //-- Path utilities

      // name of the module or entry defined by the path-name `path_and_name'
      // (similar to the unix command basename associated to file names)
      static std::string basename( std::string const& path_and_name ) ;
      
      // name of the module containing the module or entry defined by the
      // path-name `path_and_name'
      // (similar to the unix command dirname associated to file names)
      static std::string dirname( std::string const& path_and_name ) ;

      bool substitute_variables( std::string& replaced,
                                 MAC_Context const* ct ) const ;

      MAC_DataWithContext const* create_evaluation(
                                   MAC_Object * a_owner,
                                   std::string const& expression,
                                   MAC_Context const* ct ) const ;
   
   //-- Comparison

      // IMPLEMENTATION : `other' must be a `MAC_Module::' object
      // or a `MAC_String::' object
      virtual bool comparable( MAC_Object const* other ) const ;
      
      // IMPLEMENTATION : if `other' is a `MAC_Module::' object, it is equal to
      // to `self' if it has the same name ; if `other' is a `MAC_String::' 
      // object, it is equal to `self' if it represents its name
      virtual bool is_equal( MAC_Object const* other ) const ;
      
      // IMPLEMENTATION : compare `self' name to `other' name
      virtual int three_way_comparison( MAC_Object const* other ) const ;

      // IMPLEMENTATION : hash code of a `MAC_String::' object defined from
      // `self' name
      virtual size_t hash_code( void ) const ;

   protected: //-------------------------------------------------------
    	    	
   private: //-------------------------------------------------------

      MAC_Module( void ) ;
     ~MAC_Module( void ) ;
      MAC_Module( MAC_Module const& other ) ;
      MAC_Module const& operator=( MAC_Module const& other ) ;
  
      MAC_Module( MAC_Object* a_owner,
                  std::string const& a_name,
                  MAC_Context const* ct ) ;

      void complete_context( MAC_Module* root,
                             MAC_Module* dup,
                             MAC_Context const* ref_ctx ) const ;
      
      void complete_context( MAC_Module* root,
                             MAC_Data const* dat,
                             MAC_Context const* ref_ctx,
                             bool& has_modif ) const ;

      bool find( std::string const& nom,
                 MAC_Module*& theModule,
                 MAC_KeywordDataPair*& theAssignment ) const ;
      static void split( std::string const& nom,
                         std::string& token,
                         std::string& otherToken ) ;
      
      // Print entire tree and inserts tabulation for each level
      std::ostream& recursive_print( std::ostream& s,
                                     int n,
                                     bool hybrid,
                                     std::string const& bin_file ) const ;

   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;
      
   //-- Attributes

      MAC_String* const NAME ;

      // Children structure:
      MAC_Module const* FATHER ;
      MAC_List* const MODS ;     // List of MAC_Module*
      MAC_List* const ENTRIES ;  // List of MAC_KeywordDataPair*

      // Context:
      MAC_ContextSimple* const CTX ;
      MAC_ContextPair* const TMP_CTX ;

};

#endif 
