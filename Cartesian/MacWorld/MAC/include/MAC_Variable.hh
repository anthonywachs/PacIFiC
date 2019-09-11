#ifndef MAC_VARIABLE_HH
#define MAC_VARIABLE_HH

#include <MAC_Data.hh>

class MAC_Context ;
class MAC_ContextSimple ;

/*
Variables on which expressions depend.

The result of a variable evaluation in a context is its value
in that context.
Type of variable depends on its name (after $ symbol indicating that
identifier stands for a variable name):
 First character stands for element type ( D for Double, I for
  Integer, B for Boolean, S for String );
 Second one stands for dimension ( S for Scalar, V for Vector and A for
  2D-Array ).
  
Examples :
$SS_HW = "Hello, world!"   // A scalar string
$DVvar = < 0.0 1.0 2.0 >   // A vector of doubles
$BSistrue = true           // A scalar boolean
$IA_itab = array( < 0 0 > , < 1 1 > ) // An integer array
*/

class MAC_Variable : public MAC_Data
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      virtual MAC_Variable* create_clone( MAC_Object* a_owner ) const ;

      // number of instances
      static size_t nb_objects( void ) ;

      // an instance identified by `a_name'
      static MAC_Variable const* object( std::string const& a_name ) ;
      
      static MAC_Variable const* object( size_t id ) ;
      
   //-- Identification

      // number, uniquely determining `self'      
      size_t id_number( void ) const ;

      std::string const& name( void ) const ;

   //-- Context

      // IMPLEMENTATION : ensure that `lst' contains `self'.
      virtual void declare( MAC_List* lst ) const ;

      virtual bool context_has_required_variables( 
                                             MAC_Context const* ct ) const ;

   //-- Type

      // type of data of name `a_name' (raise an error in case of bad syntaxe)
      static MAC_Data::Type data_type( std::string const& a_name ) ;
      
      virtual MAC_Data::Type data_type( void ) const ;
      
  //-- Value
      
      // IMPLEMENTATION : true if `ct' associates a value to `self'
      virtual bool value_can_be_evaluated( MAC_Context const* ct ) const ;
      
      virtual stringVector const& undefined_variables(
                                           MAC_Context const* ct ) const ;

      virtual bool to_bool( MAC_Context const* ct ) const ;

      virtual double to_double( MAC_Context const* ct ) const ;

      virtual int to_int( MAC_Context const* ct ) const ;

      virtual std::string const& to_string( MAC_Context const* ct ) const ;

      virtual doubleVector const& to_double_vector(
                                            MAC_Context const* ct ) const ;

      virtual intVector const& to_int_vector( MAC_Context const* ct ) const ;

      virtual stringVector const& to_string_vector( 
                                            MAC_Context const* ct ) const ;

      virtual boolVector const& to_bool_vector( 
                                            MAC_Context const* ct ) const ;

      virtual doubleArray2D const& to_double_array2D( 
                                            MAC_Context const* ct ) const ;

      virtual intArray2D const& to_int_array2D( MAC_Context const* ct ) const ;
      
      virtual boolArray2D const& to_bool_array2D( 
                                        MAC_Context const* ct = 0 ) const ;
      
      virtual stringArray2D const& to_string_array2D( 
                                        MAC_Context const* ct = 0 ) const ;
      
      virtual doubleArray3D const& to_double_array3D( 
                                            MAC_Context const* ct ) const ;

      virtual intArray3D const& to_int_array3D( MAC_Context const* ct ) const ;
      
   //-- Formal calculus

      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct ) const ;
      virtual bool is_raw_data( void ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      MAC_Variable( void ) ;
     ~MAC_Variable( void ) ;
      MAC_Variable( MAC_Variable const& other ) ;
      MAC_Variable& operator=( MAC_Variable const& other ) ;

      MAC_Variable( MAC_Object* a_owner, std::string const& a_name ) ;

      MAC_Variable( MAC_Object* a_owner, MAC_Variable const* other ) ;

      // value of `self' in `ct'.
      MAC_Data const* data( MAC_Context const* ct ) const ;
      
      static MAC_List*  variable_list( void ) ;

   //-- Attributes

      std::string const NAME ;
      MAC_Data::Type const KIND ;
      size_t ID ;

      mutable bool EVALUATING ;
};

#endif
