#ifndef MAC_SimplifyPattern_HH
#define MAC_SimplifyPattern_HH

#include <MAC_Application.hh>
#include <string>
#include <stringVector.hh>

class MAC_Saver ;
class MAC_ModuleExplorer ;

/*
  Check validity of several data files from a pattern model
*/

class MAC_SimplifyPattern : public MAC_Application
{

   public: //-----------------------------------------------------------

      static MAC_SimplifyPattern* create( MAC_Object* a_owner,
                                           MAC_ModuleExplorer const* exp ) ;
      
   //-- Program core execution

      virtual void run( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~MAC_SimplifyPattern( void ) ; 
      MAC_SimplifyPattern( MAC_SimplifyPattern const& other ) ;
      MAC_SimplifyPattern & operator=(
                            MAC_SimplifyPattern const& other ) ;

      MAC_SimplifyPattern( MAC_Object* a_owner,
                            std::string const& pattern_file ) ;
      
   //-- Plug in
      
      MAC_SimplifyPattern( void ) ;
      
      virtual MAC_SimplifyPattern* create_replica( 
                            MAC_Object* a_owner,
                            MAC_ModuleExplorer const* exp ) const ;

      virtual MAC_SimplifyPattern* create_replica_from_args(
                            MAC_Object* a_owner,
                            stringVector& args ) const ;
      
   //-- Class attribute
      
      static MAC_SimplifyPattern const* PROTOTYPE ;
      
   //-- Attribute
      
      std::string const PATTERN ;
      
} ;

#endif



