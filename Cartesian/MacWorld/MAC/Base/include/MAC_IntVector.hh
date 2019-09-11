#ifndef MAC_INTVECTOR_HH
#define MAC_INTVECTOR_HH

#include <MAC_Data.hh>
#include <intVector.hh>

/* Implements list of int constant.
   The list can be built by adding simple value to it then finalizes it
   or from a intVector.  */

class MAC_Container ;

class MAC_IntVector : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_IntVector* create( MAC_Object* a_owner,
                                       intVector const& aIntVector ) ;
      
      static MAC_IntVector* create( MAC_Object* a_owner,
                                       MAC_Container const* aDataList ) ;
      
      virtual MAC_IntVector* create_clone( MAC_Object* a_owner ) const ;
      
   //-- Type

      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value

      virtual intVector const& to_int_vector( 
                                            MAC_Context const* ct = 0 ) const ;  
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Modifier

      void set( intVector const& other ) ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_IntVector( void ) ;
     ~MAC_IntVector( void ) ;
      MAC_IntVector( MAC_IntVector const & other ) ;
      MAC_IntVector& operator=( MAC_IntVector const& other ) ;

      MAC_IntVector( MAC_Object* a_owner,
                      intVector const& aIntVector ) ;

      MAC_IntVector( MAC_Object* a_owner,
                        MAC_Container const* aDataList,
                        bool& error ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Attributes

      intVector dv ;
} ;


#endif
