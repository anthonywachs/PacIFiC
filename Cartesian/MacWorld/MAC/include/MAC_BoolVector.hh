#ifndef MAC_BOOLVECTOR_HH
#define MAC_BOOLVECTOR_HH

#include <MAC_Data.hh>
#include <boolVector.hh>

/* Implements list of bool constant.
   The list can be built by adding simple value to it then finalizes it
   or from a boolVector.  */

class MAC_Container ;

class MAC_BoolVector : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_BoolVector* create( MAC_Object* a_owner,
                                       boolVector const& aBoolVector ) ;
      
      static MAC_BoolVector* create( MAC_Object* a_owner,
                                       MAC_Container const* aDataList ) ;
      
      virtual MAC_BoolVector* create_clone( MAC_Object* a_owner ) const ;
      
   //-- Type

      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value

      virtual boolVector const& to_bool_vector( 
                                            MAC_Context const* ct = 0 ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Modifier

      void set( boolVector const& other ) ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_BoolVector( void ) ;
     ~MAC_BoolVector( void ) ;
      MAC_BoolVector( MAC_BoolVector const & other ) ;
      MAC_BoolVector& operator=( MAC_BoolVector const& other ) ;

      MAC_BoolVector( MAC_Object* a_owner,
                      boolVector const& aBoolVector ) ;

      MAC_BoolVector( MAC_Object* a_owner,
                      MAC_Container const* aDataList,
                      bool& error ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Attributes

      boolVector dv ;
} ;


#endif
