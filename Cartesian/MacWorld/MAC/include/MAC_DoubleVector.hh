#ifndef MAC_DOUBLEVECTOR_HH
#define MAC_DOUBLEVECTOR_HH

#include <MAC_Data.hh>
#include <doubleVector.hh>

/* Implements list of double constant.
   The list can be built by adding simple value to it then finalizes it
   or from a doubleVector.  */

class MAC_Container ;

class MAC_DoubleVector : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_DoubleVector* create( MAC_Object* a_owner,
                                       doubleVector const& aDoubleVector ) ;
      
      static MAC_DoubleVector* create( MAC_Object* a_owner,
                                       MAC_Container const* aDataList ) ;
      
      virtual MAC_DoubleVector* create_clone( MAC_Object* a_owner ) const ;
      
   //-- Type

      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value

      virtual doubleVector const& to_double_vector( 
                                            MAC_Context const* ct = 0 ) const ;
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Modifier

      void set( doubleVector const& other ) ;
      
   //-- Formal calculus
      
      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct ) const ;

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_DoubleVector( void ) ;
     ~MAC_DoubleVector( void ) ;
      MAC_DoubleVector( MAC_DoubleVector const & other ) ;
      MAC_DoubleVector& operator=( MAC_DoubleVector const& other ) ;

      MAC_DoubleVector( MAC_Object* a_owner,
                      doubleVector const& aDoubleVector ) ;

      MAC_DoubleVector( MAC_Object* a_owner,
                        MAC_Container const* aDataList,
                        bool& error ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Attributes

      doubleVector dv ;
} ;


#endif
