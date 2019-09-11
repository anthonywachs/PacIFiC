#ifndef MAC_SYSTEM_EXP_HH
#define MAC_SYSTEM_EXP_HH

#include <MAC_Expression.hh>

/* System expressions.
   There are :
     - getcwd() : return current working directory
     - getenv( varname ) : return value of some environement variable
                           (it must exist)
     - join( path1, path2, .. ) : The return value is the concatenation
            of path1,  path2, and optionally path3, etc., with exactly one
            // slash ('/') inserted between components.

PUBLISHED
*/

class MAC_SystemExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
  //-- Value
      
      virtual std::string const& to_string( MAC_Context const* ct = 0 ) const ;
      
   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------

      MAC_SystemExp( void ) ;
     ~MAC_SystemExp( void ) ;
      MAC_SystemExp( MAC_SystemExp const& other ) ;
      MAC_SystemExp& operator=( MAC_SystemExp const& other ) ;

      MAC_SystemExp( MAC_Object* a_owner,
                     std::string const& a_name,
                     MAC_Sequence const* argument_list ) ;

   //-- Plug in

      MAC_SystemExp( std::string const& a_name ) ;

      virtual MAC_SystemExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes
            
      static MAC_SystemExp const* PROTOTYPE_PWD ;
      static MAC_SystemExp const* PROTOTYPE_GETENV ;
      static MAC_SystemExp const* PROTOTYPE_DIRNAME ;
      static MAC_SystemExp const* PROTOTYPE_BASENAME ;
      static MAC_SystemExp const* PROTOTYPE_SEPARATOR ;
      static MAC_SystemExp const* PROTOTYPE_JOIN ;
      static MAC_SystemExp const* PROTOTYPE_GETPID ;
      static MAC_SystemExp const* PROTOTYPE_UNAME ;
      static MAC_SystemExp const* PROTOTYPE_HOSTNAME ;

   //-- Attributes
      
      mutable std::string RESULT_STR ;
} ;

#endif
