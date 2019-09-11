#ifndef MAC_VARIABLE_EXP_HH
#define MAC_VARIABLE_EXP_HH

#include <MAC_Expression.hh>

/*
Operators on `MAC_Variable::' objects.

- check if a given variable in an HDS

  if( is_defined( "DS_value" ) )
  MODULE toto
     ...
  END MODULE toto

  Module toto is built only if a variable of name $DS_value has been
  defined.

- value of a variable with default value

  my_variable = value( "DS_value", 1.E-5 )

  my_variable is set with the value of the variable $DS_value this variable
  is defined, or with the default value 1.E-5 elsewhere.
  
*/

class MAC_VariableExp : public MAC_Expression
{
   public: //---------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual bool to_bool( MAC_Context const* ct ) const ;
      virtual double to_double( MAC_Context const* ct ) const ;
      virtual int to_int( MAC_Context const* ct ) const ;
      virtual std::string const& to_string( MAC_Context const* ct ) const ;
      virtual doubleVector const& to_double_vector( MAC_Context const* ct ) const ;
      virtual intVector const& to_int_vector( MAC_Context const* ct ) const ;
      virtual stringVector const& to_string_vector(
                                              MAC_Context const* ct ) const ;
      virtual boolVector const& to_bool_vector( MAC_Context const* ct ) const ;
      virtual doubleArray2D const& to_double_array2D( MAC_Context const* ct ) const ;
      virtual intArray2D const& to_int_array2D( MAC_Context const* ct ) const ;
      virtual boolArray2D const& to_bool_array2D( 
                                        MAC_Context const* ct = 0 ) const ;      
      virtual stringArray2D const& to_string_array2D( 
                                        MAC_Context const* ct = 0 ) const ;
 
   //-- Formal calculus

      virtual bool is_constant( void ) const ;

   protected: //-------------------------------------------------------

   private: //---------------------------------------------------------

      MAC_VariableExp( void ) ; 
     ~MAC_VariableExp( void ) ;
      MAC_VariableExp( MAC_VariableExp const& other ) ;
      MAC_VariableExp& operator=( MAC_VariableExp const& other ) ;

      enum VarExp{ var_def, var_value } ;

      MAC_VariableExp( MAC_Object* a_owner,
                       VarExp exp_id,
                       std::string const& a_name,
                       MAC_Sequence const* argument_list ) ;
      
   //-- Plug in

      MAC_VariableExp( VarExp exp_id, std::string const& a_name ) ;
      
      virtual MAC_VariableExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;
   
      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static MAC_VariableExp const* PROTOTYPE_var_def ;
      static MAC_VariableExp const* PROTOTYPE_var_value ;

   //-- Attributes

      VarExp const OP ;
} ;

#endif
