#ifndef MAC_MODULE_EXPLORER_HH
#define MAC_MODULE_EXPLORER_HH

#include <MAC_Object.hh>
#include <MAC_Data.hh>

#include <string>

class MAC_Application ;
class MAC_Container ;
class MAC_Context ;
class MAC_ContextPair ;
class MAC_DataWithContext ;
class MAC_List ;
class MAC_Module ;
class MAC_ModuleIterator ;
class MAC_ModulePattern ;
class MAC_KeywordDataIterator ;
class boolVector ;
class doubleVector ;
class doubleArray2D ;
class boolArray2D ;
class stringArray2D ;
class doubleArray3D ;
class intArray2D ;
class intArray3D ;
class intVector ;
class stringVector ;

/*
Navigators to interrogate a database of the MAC Hierarchical Data System.

Each `MAC_ModuleExplorer::' object is attached to a `MAC_Module::' object
and offers ad hoc navigational and query facilities for accessing the data
of the associated Module Hierarchy.
Moreover, MAC_ModuleExplorer objects can be linked with a pattern recognition
mecanism. Two modes are then available :
build mode : pattern describing current HDS is build ;
verify mode : current HDS is compared to match a pattern and `validity' feature
 become available.
Since each `MAC_Module::' object defines a Module Hierarchy, each 
`MAC_ModuleExplorer::' object is said to be attached to a Module Hierarchy.

DESCRIPTION OF THE MAC HIERARCHICAL DATA SYSTEM
----------------------------------------------------

Data are organized into a hierarchy. The nodes of the hierarchy are called 
modules (hence the name: Module Hierarchy) while the leaves are the 
data themselves.

Each module may contain :
   - other modules, and 
   - entries represented by a keyword-data pair (thus, within a module,
     data are uniquely identified by a keyword).

Each module defines a Module Hierarchy (a module contains other modules
themselves containing other modules and so on).
The upper module of a Module Hierarchy is denoted the root module.
 
Given the root module of a Module Hierarchy, modules and entries can be 
retrieved knowing a path-name. A path-name is the name of a module or 
of the keyword of an entry preceeded by a concatenation of 
one or several module names separated with "/" . 
The starting node of a path is the root module if the first character is "/" 
(absolute path-name), otherwise (relative path-name) that
starting node is unspecified  (thus path-name looks like a unix file name). 

Example : in the Module Hierarchy defined by :
   MODULE mod1
     MODULE mod2
       k = "v"
     END MODULE mod2
   END MODULE mod1 
the keyword-data pair  {k,"v"} can be retrieved in module mod1 with either of 
the following path-names : 
        "/mod1/mod2/k" (absolute path-name)
        "mod2/k"       (relative path-name)
        "k"            (relative path-name)

Note that modules and entries might not be uniquely defined by
relative path-names. They are uniquely defined by absolute path-names.

Each `MAC_ModuleExplorer::' object is attached to a `MAC_Module::' object.
*/

class MAC_ModuleExplorer : public MAC_Object
{
   public: //----------------------------------------------------------------

   //-- Module pattern

      enum PatternStatus { ignore, build, verify } ;
      
      // status of `self' relative to pattern recognition
      PatternStatus pattern_status( void ) const ;
      
      // pattern module used to describe pattern ( can be null if no pattern 
      // exist for associated module)
      MAC_ModulePattern* pattern( void ) const ;
      
      // compare `self' to underlying pattern and return diagnostic message
      //  if pattern is not respected
      // `recurse' is used to checked recursively all sub-modules
      MAC_ModuleExplorer* validity( MAC_Object* a_owner,
                                    bool recurse=true ) const ;
      
   //-- Instance delivery and initialization

      // Create and return an instance attached to `mm'.
      static MAC_ModuleExplorer* create( MAC_Object* a_owner,
                                         MAC_Module const* mm,
                                         PatternStatus status = ignore ) ;

      // Create and return an instance attached to the first module identified
      // by the path-name `path_and_name' in the attached Module Hierarchy.
      MAC_ModuleExplorer* create_subexplorer( 
                                 MAC_Object* a_owner,
                                 std::string const& path_and_name ) const ;

      virtual MAC_ModuleExplorer* create_clone( MAC_Object* a_owner ) const ;

   //-- Status report

      // name of the attached module
      std::string const& name( void ) const ;

      // absolute path-name of `self'
      std::string const& absolute_path_name( void ) const ;

      // owner of the attached module
      MAC_Object const* owner_of_attached_module( void ) const ;

      // Is the attached module empty ?
      bool is_empty( void ) const ;
      
      // Does the attached module contains a module of name `a_path_and_name' ?
      bool has_module( std::string const& a_path_and_name ) const ;

      // Does the attached module contains an entry whose keyword is 
      // `a_path_and_keyword' ?
      bool has_entry( std::string const& a_path_and_keyword ) const ;

  //-- Iterator on modules contained in the attached module
      
      // Move module iterator to first position.
      void start_module_iterator( void ) ;

      // Is module iterator position valid ?
      bool is_valid_module( void ) const ;

      // Move module iterator one position within the modules contained
      // in the attached module.
      void go_next_module( void ) ;

      // Create and return an instance attached to the module at
      // module iterator position.      
      MAC_ModuleExplorer* create_subexplorer( MAC_Object* a_owner ) const ;

  //-- Iterator on entries contained in the attached module
      
      // Move entry iterator to first position.
      void start_entry_iterator( void ) ;

      // Is entry iterator position valid ?
      bool is_valid_entry( void ) const ;

      // Move entry iterator one position within the entries contained in the
      // attached module.
      void go_next_entry( void ) ;

      // keyword of the entry at current entry iterator position
      std::string const& keyword( void ) const ;
      
      // data of the entry at current entry iterator position
      MAC_DataWithContext* data( MAC_Object* a_owner,
                                 MAC_Context const* ct = 0 ) const ;

  //-- Access to data contained in the attached module
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none)
      MAC_DataWithContext*
      abstract_data( MAC_Object* a_owner,
                     std::string const& a_path_and_keyword,
                     MAC_Context const* ct = 0 ) const  ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      bool bool_data( std::string const& a_path_and_keyword,
                      MAC_Context const* ct = 0  ) const ;

      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      int int_data( std::string const& a_path_and_keyword,
                    MAC_Context const* ct = 0  ) const ;

      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      double double_data( std::string const& a_path_and_keyword,
                          MAC_Context const* ct = 0  ) const ;

      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      std::string const& string_data( std::string const& a_path_and_keyword,
                                      MAC_Context const* ct = 0  ) const ;

      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      intVector const&  intVector_data( std::string const& a_path_and_keyword,
                                        MAC_Context const* ct = 0  ) const ;

      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      doubleVector const& doubleVector_data( std::string const& a_path_and_keyword,
                                             MAC_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      doubleArray2D const& doubleArray2D_data( std::string const& a_path_and_keyword,
                                               MAC_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      stringArray2D const& stringArray2D_data( std::string const& a_path_and_keyword,
                                               MAC_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      boolArray2D const& boolArray2D_data( std::string const& a_path_and_keyword,
                                               MAC_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      doubleArray3D const& doubleArray3D_data( std::string const& a_path_and_keyword,
                                               MAC_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      intArray2D const& intArray2D_data( std::string const& a_path_and_keyword,
                                         MAC_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      intArray3D const& intArray3D_data( std::string const& a_path_and_keyword,
                                         MAC_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      boolVector const& boolVector_data( std::string const& a_path_and_keyword,
                                         MAC_Context const* ct = 0  ) const ;
      
      // data called `a_path_and_keyword' in the attached module
      // (fatal error raised if none or if not of the proper type)
      stringVector const& stringVector_data( std::string const& a_path_and_keyword,
                                             MAC_Context const* ct = 0  ) const  ;

      MAC_Module* create_clone_of_attached_module( MAC_Object* a_owner ) const ;
      
   //-- Meta data information
      
      void set_help( std::string const& a_keyword,
                     std::string const& expression ) const ;
      void set_default( std::string const& a_keyword,
                        std::string const& expression ) const ;
      void test_data( std::string const& a_keyword,
                      std::string const& expression,
                      MAC_Context const* ct = 0 ) const ;
      void test_data_in( std::string const& a_keyword,
                         stringVector const& choices ) const ;
      void test_data_as( std::string const& a_keyword,
                         std::string const& regexp,
                         std::string const& where="true" ) const ;
      void test_file( std::string const& filename,
                      std::string const& mode ) const ;
      
   //-- Input - Output

      // IMPLEMENTATION : write the Module Hierarchy of root `self' with the
      // format of a MAC data file.
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
      void write( std::string const& file,
                  std::string const& format ) const ;

      
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Module pattern
      
      bool build_pattern( void ) const ;

      void declare_data( std::string const& a_path_and_keyword,
                         MAC_Data::Type type ) const ;
      
   //-- Access to context
      
      MAC_Context const* context( std::string const& a_path_and_name,
                                  MAC_Context const* ct ) const ;
      
      MAC_ModuleExplorer( void ) ;
     ~MAC_ModuleExplorer( void ) ;
      MAC_ModuleExplorer( MAC_ModuleExplorer const& other ) ;
      MAC_ModuleExplorer const& operator=(
          	          MAC_ModuleExplorer const& other ) ;

      MAC_ModuleExplorer( MAC_Object* a_owner,
                          MAC_Module const* mm,
                          PatternStatus status ,
                          MAC_ModulePattern* pat ) ;

      //---------------------------------------------------------------------
      //   ATTRIBUTES
      //---------------------------------------------------------------------
      MAC_Module const* mod ;
      MAC_ModuleIterator * submodule_iterator ;
      MAC_KeywordDataIterator * keyword_iterator ;
      mutable MAC_ContextPair* tmp_context ;
      
      MAC_ModulePattern* MP ;
      bool keyword_iterator_started ;
      PatternStatus my_status ;
} ;

#endif
