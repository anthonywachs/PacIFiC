#ifndef MAC_VECTOR_EXP_HH
#define MAC_VECTOR_EXP_HH

#include <MAC_Expression.hh>

#include <doubleVector.hh>
#include <intVector.hh>
#include <boolVector.hh>
#include <stringVector.hh>

class MAC_Iterator ;

/*
Operators to form vector from simple items:
   - vector of doubles:
       vector( 0., 3., 4., -1.E3 )              ->  < 0. 3. 4. -1.E3 >
   - vector of integers:
       vector( 0, 3, 4, -5 )                    ->  < 0 3 4 -5 >
   - vector of strings:
       vector( "MAC", "is", "beautiful" )  ->  < "MAC" "is" "beautiful" >
   - vector of booleans:
       vector( true, false )                    ->  < true, false >

   - nvector( 2, 1.0 )                          ->  < 1.0 1.0 >
   
   - nvector( 3, true )                         ->  < true true true >
   
Operator to retrieve the size of such vectors:
   size( vector( 0., 3., 4., -1.E3 ) )          ->  4

Operator to retrieve component from vector:
   First arg is a vector and second one an integer.
   component( vector( "MAC", "is", "beautiful") ), 1 ) -> "MAC"

Operator to test if the elements of a vector are in increasing order
   increasing( < 0.0 1.0 1.0 2.0 > )  -> true
   increasing( < 0 1 2 1 > )          -> false

Operator to test if all the elements of a vector are greater than or equal to
a given value
   greater( < -10 3 >, 1 )      -> false
   greater( < 1.0 3.0 >, 0.5 )  -> true

Operator to form vector from another one:
   apply( < 1. 4. 9. 16. >, $DS_x*$DS_x, "DS_x" ) -> < 1. 4. 9. 16. >
   apply( < 1  4  9  16  >, $IS_x*$IS_x, "IS_x" ) -> < 1  4  9  16  >
   apply( < true true false >, ! $BS_x, "BS_x" )  -> < false false true >
   apply( < "titi" "toto" >, $SS_x+"0", "SS_x" )  -> < "titi0" "toto0" >
   apply( < "titi" "toto" >, $SS_x+to_string($IS_ic), "SS_x", "SS_ic" )
                                                  -> < "titi0" "toto1" >

Operator to sum the elements of a vector :
   sum( < 0. 1.0 2.0 > )  ->  3.0
   
Operator to reverse order of vector :
   reverse( < 0. 1.0 2.0 > )  ->  < 2.0 1.0 0. >
   
PUBLISHED
*/

class MAC_VectorExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
 
   //-- Context
                  
      virtual void declare( MAC_List* lst ) const ;

      virtual bool context_has_required_variables( 
                                           MAC_Context const* ct ) const ;
      
   //-- Value
      
      virtual bool value_can_be_evaluated( MAC_Context const* ct ) const ;
      
      virtual stringVector const& undefined_variables(
                                           MAC_Context const* ct ) const ;
      
      virtual bool to_bool( MAC_Context const* ct = 0 ) const ;

      virtual double to_double( MAC_Context const* ct = 0 ) const ;

      virtual int to_int( MAC_Context const* ct = 0 ) const ;

      virtual std::string const& to_string( MAC_Context const* ct = 0 ) const ;
      
      virtual doubleVector const& to_double_vector(
                                       MAC_Context const* ct = 0 ) const ;
      
      virtual intVector const& to_int_vector( 
                                       MAC_Context const* ct = 0 ) const ;

      virtual stringVector const& to_string_vector( 
                                       MAC_Context const* ct = 0 ) const ;
      
      virtual boolVector const& to_bool_vector(
                                       MAC_Context const* ct = 0 ) const ;
      
   //-- Formal calculus
      
      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      MAC_VectorExp( void ) ;
     ~MAC_VectorExp( void ) ;
      MAC_VectorExp( MAC_VectorExp const& other ) ;
      MAC_VectorExp& operator=( MAC_VectorExp const& other ) ;

      enum VectorExp{ vect, size, component,
                      increasing, greater,
                      nvector, cond_vect,
                      apply, reverse, sum } ;

      MAC_VectorExp( MAC_Object* a_owner,
                     VectorExp exp_id,
                     std::string const& a_name,
                     MAC_Sequence const* argument_list ) ;

      MAC_Context const* create_apply_context(
                     MAC_Object* a_owner,
                     MAC_Context const* ctx,
                     MAC_Variable const* v,
                     MAC_Variable const* vic ) const ;
                     
   //-- Plug in

      MAC_VectorExp( VectorExp exp_id, std::string const& a_name ) ;
      
      virtual MAC_VectorExp* create_replica( 
                                    MAC_Object * a_owner,
                                    MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes      

      static MAC_VectorExp const* PROTOTYPE_CONDITIONAL_VECTOR ;
      static MAC_VectorExp const* PROTOTYPE_VECTOR ;
      static MAC_VectorExp const* PROTOTYPE_NVECTOR ;
      static MAC_VectorExp const* PROTOTYPE_SIZE ;
      static MAC_VectorExp const* PROTOTYPE_COMPO ;
      static MAC_VectorExp const* PROTOTYPE_INCREASING ;
      static MAC_VectorExp const* PROTOTYPE_GREATER ;
      static MAC_VectorExp const* PROTOTYPE_APPLY ;
      static MAC_VectorExp const* PROTOTYPE_REVERSE ;
      static MAC_VectorExp const* PROTOTYPE_SUM ;

   //-- Attributes

      VectorExp const OP ;
      
      mutable MAC_Iterator* MY_IT ; // To speed up evaluation
      
      mutable doubleVector RESULT_D ;
      mutable intVector    RESULT_I ;
      mutable boolVector   RESULT_B ;
      mutable stringVector RESULT_S ;
      
} ;

#endif
