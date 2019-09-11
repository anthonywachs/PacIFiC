#ifndef MAC_STRINGVECTOR_HH
#define MAC_STRINGVECTOR_HH

#include <MAC_Data.hh>
#include <stringVector.hh>

/* Implements list of string constant.
   The list can be built by adding simple value to it then finalizes it
   or from a stringVector.  */

class MAC_Container ;

class MAC_StringVector : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_StringVector* create( MAC_Object* a_owner,
                                       stringVector const& aStringVector ) ;
      
      static MAC_StringVector* create( MAC_Object* a_owner,
                                       MAC_Container const* aDataList ) ;
      
      virtual MAC_StringVector* create_clone( MAC_Object* a_owner ) const ;
      
   //-- Type

      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value

      virtual stringVector const& to_string_vector( 
                                            MAC_Context const* ct = 0 ) const ;  
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Modifier

      void set( stringVector const& other ) ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_StringVector( void ) ;
     ~MAC_StringVector( void ) ;
      MAC_StringVector( MAC_StringVector const & other ) ;
      MAC_StringVector& operator=( MAC_StringVector const& other ) ;

      MAC_StringVector( MAC_Object* a_owner,
                      stringVector const& aStringVector ) ;

      MAC_StringVector( MAC_Object* a_owner,
                        MAC_Container const* aDataList,
                        bool& error ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Attributes

      stringVector dv ;
} ;


#endif
