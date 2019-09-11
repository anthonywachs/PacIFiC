#ifndef MAC_COMPARATOR_HH
#define MAC_COMPARATOR_HH

#include <MAC_Application.hh>

class MAC_ModuleExplorer ;

/*
PUBLISHED
*/

class MAC_Comparator : public MAC_Application
{
   public: //-----------------------------------------------------------------

   //-- File formats
      
      static void set_preferred_motifs_formats( 
                                           MAC_ModuleExplorer const* exp ) ;
         
      static void detect_file_format( MAC_ModuleExplorer const* exp,
                                      std::string const& file_name,
                                      std::string& format ) ;
         
   //-- Program core execution

      virtual void run( void ) ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~MAC_Comparator( void ) ;
      MAC_Comparator( MAC_Comparator const& other ) ;
      MAC_Comparator& operator=( MAC_Comparator const& other ) ;

      MAC_Comparator( MAC_Object* a_owner,
                      MAC_ModuleExplorer const* exp ) ;

      MAC_Comparator( MAC_Object* a_owner,
                      stringVector& args ) ;

   //-- Plug in

      MAC_Comparator( void ) ;

      virtual MAC_Comparator* create_replica(
                                        MAC_Object* a_owner,
                                        MAC_ModuleExplorer const* exp ) const ;

      virtual MAC_Comparator* create_replica_from_args(
                                        MAC_Object* a_owner,
                                        stringVector& args ) const ;
      
   //-- Command line

      virtual void print_usage( void ) const ;
      virtual void print_operands( void ) const ;
      virtual void print_exit_status( void ) const ;

   //-- Internals

      //????? changer ces noms
      void builtin_compare( std::ostream& os ) ;
      
      void do_diff( std::ostream& os ) ;
      
      void display_glob_info_diff( std::ostream& os,
                                   MAC_Module const* mod ) const ;

      void display_submodule_diff( std::ostream& os,
                                   MAC_Module const* mod ) const ;

   //-- Motifs
      
      static stringVector& pref_motifs( void ) ;
      static stringVector& pref_formats( void ) ;
      
   //-- Class attributes

      static MAC_Comparator const* PROTOTYPE ;
      
   //-- Attributes

      std::string FILE1 ;
      std::string FILE2 ;
      std::string FILE_OUT ;
      MAC_ModuleExplorer* EXP ;
      bool VERBOSE ;
      double MY_DBL_EPS ;
      double MY_DBL_MIN ;
      stringVector IGNORE ;
      std::string FORMAT ;
      
} ;

#endif
