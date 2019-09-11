#ifndef MAC_DATA_HH
#define MAC_DATA_HH

#include <MAC_Object.hh>

class boolVector ;
class boolArray2D ;
class doubleArray2D ;
class doubleArray3D ;
class doubleVector ;
class intArray2D ;
class intArray3D ;
class intVector ;
class stringVector ;
class stringArray2D ;
class MAC_Context ;
class MAC_List ;
class MAC_ContextSimple ;
class MAC_Variable ;

/*
Data of the MAC Hierarchical Data System.

A Data has two characteristics : 
   - a type, and
   - a value.

The value of a data might depend of a context (but not its type).

A context is  represented by a `MAC_Context::' object (possibly NULL in
which case the type and the value are context independant). 

To evaluate a data in a given context (resulting in its type and value) :
   1. the required variables should occur in the context (which can be
      achieved by calling `declare' ;
   2. the context should assign relevant values to these required variables.
*/

class MAC_Data : public MAC_Object
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      virtual MAC_Data* create_clone( MAC_Object* a_owner ) const = 0 ;
      
   //-- Context(10.)
                  
      // Ensure that `lst' contains all necessary instances of
      // `MAC_Variable::' for the evaluation of `self'.
      // IMPLEMENTATION : do nothing : `lst' is left unchanged.
      virtual void declare( MAC_List* lst ) const ;

      // Does `ct' contains all necessary variables for the evaluation 
      // of `self' (the variables of `ct' can not necessary be evaluated) ?
      // IMPLEMENTATION : true
      virtual bool context_has_required_variables( 
                                               MAC_Context const* ct ) const ;
      
   //-- Type(20.)

      enum Type { Undefined, Double, Int, Bool, String,
                  DoubleVector, IntVector, BoolVector, StringVector,
                  DoubleArray2D, IntArray2D, BoolArray2D, StringArray2D,
                  DoubleArray3D, IntArray3D } ;

      static std::string type_name( Type kind ) ;
      
      // type of self
      virtual Type data_type( void ) const = 0 ;
      
  //-- Value(30.)
      
      // Is it possible to evaluate the value according to `ct' ?
      // IMPLEMENTATION : true
      virtual bool value_can_be_evaluated( MAC_Context const* ct ) const ;

      // undefined variable names for the evaluation of the value according to `ct'
      // IMPLEMENTATION : empty vector
      virtual stringVector const& undefined_variables(
                                           MAC_Context const* ct ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual bool to_bool( MAC_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual double to_double( MAC_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual int to_int( MAC_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual std::string const& to_string( MAC_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual doubleVector const& to_double_vector(
                                       MAC_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual intVector const& to_int_vector( 
                                       MAC_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual stringVector const& to_string_vector(
                                        MAC_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual boolVector const& to_bool_vector(
                                        MAC_Context const* ct = 0 ) const ;
      
      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual doubleArray2D const& to_double_array2D(
                                        MAC_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual intArray2D const& to_int_array2D( 
                                        MAC_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual boolArray2D const& to_bool_array2D( 
                                        MAC_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual stringArray2D const& to_string_array2D( 
                                        MAC_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual doubleArray3D const& to_double_array3D(
                                        MAC_Context const* ct = 0 ) const ;

      // value evaluated according to `ct'
      // IMPLEMENTATION : a fatal error is raised
      virtual intArray3D const& to_int_array3D( 
                                        MAC_Context const* ct = 0 ) const ;

      // value represented as a string evaluated according to `ct'
      virtual std::string value_as_string( 
                                        MAC_Context const* ct = 0 ) const ;

   //-- Formal calculus

      // Return partial derivative of `self' with respect to `var'
      // expanding expressions from `ct'.
      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct  ) const ;

      // Create a simplified expression of `self'.
      MAC_Data* create_simplification( MAC_Object* a_owner ) const ;
      
      // Is `self' a constant value ?
      virtual bool is_constant( void ) const ;
      
      // Is `self' a raw data ?
      // IMPLEMENTATION : returns true
      virtual bool is_raw_data( void ) const ;
      
   protected: //-------------------------------------------------------
      
      virtual ~MAC_Data( void ) ;

      MAC_Data( MAC_Object* a_owner ) ;

   //-- Formal calculus

      // Called by `::create_simplification'.
      virtual MAC_Data* create_non_const_simplification(
                                               MAC_Object* a_owner ) const ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool declare_PRE( MAC_List const* lst ) const ;
      virtual bool declare_POST( MAC_List const* lst ) const ;
      virtual bool context_has_required_variables_PRE( 
                                             MAC_Context const* ct ) const ;
      virtual bool to_double_PRE( MAC_Context const* ct ) const ;
      virtual bool to_int_PRE( MAC_Context const* ct ) const ;
      virtual bool to_bool_PRE( MAC_Context const* ct ) const ;
      virtual bool to_string_PRE( MAC_Context const* ct ) const ;
      virtual bool to_double_vector_PRE( MAC_Context const* ct ) const ;
      virtual bool to_int_vector_PRE( MAC_Context const* ct ) const ;
      virtual bool to_string_vector_PRE( MAC_Context const* ct ) const ;
      virtual bool to_bool_vector_PRE( MAC_Context const* ct ) const ;
      virtual bool to_double_array2D_PRE( MAC_Context const* ct ) const ;
      virtual bool to_int_array2D_PRE( MAC_Context const* ct ) const ;
      virtual bool to_bool_array2D_PRE( MAC_Context const* ct ) const ;
      virtual bool to_string_array2D_PRE( MAC_Context const* ct ) const ;
      virtual bool to_double_array3D_PRE( MAC_Context const* ct ) const ;
      virtual bool to_int_array3D_PRE( MAC_Context const* ct ) const ;
      virtual bool create_derivative_PRE( MAC_Object* a_owner,
                                          MAC_Variable const* var,
                                          MAC_Context const* ct  ) const ;
      virtual bool create_derivative_POST( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Data const* result ) const ;
      virtual bool create_non_const_simplification_PRE(
                                           MAC_Object* a_owner ) const ;
      virtual bool create_non_const_simplification_POST(
                                           MAC_Object* a_owner,
                                           MAC_Data const* result ) const ;
     
      virtual bool invariant( void ) const ;
      virtual bool is_raw_data_POST( bool result ) const ;
      
   private: //-------------------------------------------------------

      MAC_Data( void ) ;
      MAC_Data( MAC_Data const& other ) ;
      MAC_Data& operator=( MAC_Data const& other ) ;

      void exitWithError( std::string const& mess ) const ;
};

#endif
