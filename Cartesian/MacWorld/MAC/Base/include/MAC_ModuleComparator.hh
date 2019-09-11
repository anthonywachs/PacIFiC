#ifndef MAC_MODULE_COMPARATOR_HH
#define MAC_MODULE_COMPARATOR_HH

#include <MAC_Object.hh>

#include <string>

class MAC_Module;
class MAC_ModuleExplorer ;
class MAC_Context;
class MAC_Data;

class MAC_ModuleComparator : public MAC_Object
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization
      
      // Creates and return an instance containing neither modules nor entries.
      static MAC_ModuleComparator* create( MAC_Object* a_owner, 
                                           MAC_ModuleExplorer const* exp ) ;

      int compare( MAC_Module const* m1, 
                   MAC_Module const* m2, 
                   MAC_Module* result ) ;

   protected: //-------------------------------------------------------

   private: //-------------------------------------------------------

      MAC_ModuleComparator( void ) ;
     ~MAC_ModuleComparator( void ) ;
      MAC_ModuleComparator( MAC_ModuleComparator const& other ) ;
      MAC_ModuleComparator const& operator=( 
                            MAC_ModuleComparator const& other ) ;

      MAC_ModuleComparator( MAC_Object* a_owner, 
                            MAC_ModuleExplorer const* exp ) ;

      bool is_valid_module( std::string const& name ) ;
      bool is_valid_data( std::string const& name, std::string const& abs_path_name ) ;
      bool is_verbose( void ) ;
      int internalCompare( MAC_Module const * m1, 
                           MAC_Module const * m2, 
                           MAC_Module* result, 
                           std::string const& path) ;
      int compare( std::string const& dataname, 
                   MAC_Module const* m1, 
                   MAC_Module const* m2, 
                   MAC_Module* result) ;
      int compare( double v1, double v2, double& status ) ;
      int compare( bool v1, bool v2 ) ;
      int compare( int v1, int v2 ) ;
      int compare( std::string const& v1, std::string const& v2 ) ;

   //-- Attributes

      std::string left;
      std::string right;
      bool verbose;
      class stringVector *valid_module;
      class stringVector *valid_data;
      class stringVector *ignore_data;
      double MY_DBL_EPS ;
      double MY_DBL_MIN ;
      
};

#endif
