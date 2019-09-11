#ifndef MAC_VECTOR_EXP_HH
#define MAC_VECTOR_EXP_HH

#include <MAC_Expression.hh>

#include <boolArray2D.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <stringArray2D.hh>

/*

Operator to form array from simple vectors.

PUBLISHED
*/

class MAC_ArrayExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value

      virtual doubleArray2D const& to_double_array2D(
                                        MAC_Context const* ct = 0 ) const ;

      virtual intArray2D const& to_int_array2D( 
                                        MAC_Context const* ct = 0 ) const ;
      
      virtual boolArray2D const& to_bool_array2D( 
                                        MAC_Context const* ct = 0 ) const ;
      
      virtual stringArray2D const& to_string_array2D( 
                                        MAC_Context const* ct = 0 ) const ;
      
      virtual doubleArray3D const& to_double_array3D(
                                        MAC_Context const* ct = 0 ) const ;

      virtual intArray3D const& to_int_array3D( 
                                        MAC_Context const* ct = 0 ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

     ~MAC_ArrayExp( void ) ;
      MAC_ArrayExp( MAC_ArrayExp const& other ) ;
      MAC_ArrayExp& operator=( MAC_ArrayExp const& other ) ;
      
      MAC_ArrayExp( MAC_Object* a_owner,
                    MAC_Sequence const* argument_list ) ;
      
   //-- Plug in

      MAC_ArrayExp( void ) ;

      virtual MAC_ArrayExp* create_replica( 
                                  MAC_Object * a_owner,
                                  MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments(
                                 MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes
            
      static MAC_ArrayExp const* PROTOTYPE ;

   //-- Attributes

      mutable doubleArray2D RESULT_D2D ;
      mutable intArray2D    RESULT_I2D ;
      mutable boolArray2D   RESULT_B2D ;
      mutable stringArray2D RESULT_S2D ;
      mutable doubleArray3D RESULT_D3D ;
      mutable intArray3D    RESULT_I3D ;
} ;

#endif
