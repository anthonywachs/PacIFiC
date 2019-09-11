#include <MAC_VariableExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Context.hh>
#include <MAC_Error.hh>
#include <MAC_Sequence.hh>
#include <MAC_Variable.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <stringVector.hh>
#include <boolArray2D.hh>
#include <stringArray2D.hh>
#include <doubleArray2D.hh>
#include <intArray2D.hh>

MAC_VariableExp const*  
MAC_VariableExp::PROTOTYPE_var_def = 
                             new MAC_VariableExp( var_def, "is_defined" ) ;

MAC_VariableExp const*  
MAC_VariableExp::PROTOTYPE_var_value = 
                             new MAC_VariableExp( var_value, "value" ) ;

struct MAC_VariableExp_ERROR
{
   static void n0( std::string const& f_name, std::string const& op_name ) ;
   static void n1( std::string const& a_name ) ;
   static void n2( std::string const& a_name,
                   stringVector const& undef_var ) ;
} ;

//----------------------------------------------------------------------
MAC_VariableExp:: MAC_VariableExp( VarExp exp_id,
                                   std::string const& a_name  ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( exp_id )
{
   MAC_LABEL( "MAC_VariableExp:: MAC_VariableExp" ) ;
}

//----------------------------------------------------------------------
MAC_VariableExp*
MAC_VariableExp:: create_replica( MAC_Object* a_owner,
                                  MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_VariableExp* result = new MAC_VariableExp( a_owner,
                                                  OP,
                                                  name(),
                                                  argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_VariableExp:: MAC_VariableExp( MAC_Object* a_owner,
                                   VarExp exp_id,
                                   std::string const& a_name,
                                   MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP( exp_id )
{
   MAC_LABEL( "MAC_VariableExp:: MAC_VariableExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_VariableExp:: ~MAC_VariableExp( void ) 
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
std::string const& 
MAC_VariableExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   switch( OP )
   {
      case var_def :
         result = name() + "(SS)" ;
	 break ;
      case var_value :
         result = name() + "(SS,<default value>)" ;
	 break ;
      default :
         MAC_VariableExp_ERROR::n0( "usage", name() ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_VariableExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = false ;
   switch( OP )
   {
      case var_def :
         result = some_arguments->count()==1 ;
	 if( result )
	 {
	    Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
	    result = ( t0==String ) ;
	 }
	 break ;
      case var_value :
         result = some_arguments->count()==2 ;
	 if( result )
	 {
	    Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
	    result = ( t0==String ) ;
	 }
	 break ;
      default :
         MAC_VariableExp_ERROR::n0( "valid_arguments", name() ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_VariableExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: data_type" ) ;
   
   MAC_Data::Type result = MAC_Data::Undefined ;
   switch( OP )
   {
      case var_def :
         result = MAC_Data::Bool ;
         break ;
      case var_value :
         result = arg(1)->data_type() ;
	 break ;
      default :
         MAC_VariableExp_ERROR::n0( "data_type", name() ) ;
         break ;   
   }   
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_VariableExp:: is_constant( void ) const
//----------------------------------------------------------------------
{
   return( false ) ;
}

//----------------------------------------------------------------------
bool
MAC_VariableExp:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_bool" ) ;
   MAC_CHECK_PRE( to_bool_PRE( ct ) ) ;

   bool result = false ;
   switch( OP )
   {
      case var_def :
         if( ct != 0 )
         {
            MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
            result = ct->has_variable( var ) ;
         }
         break ;
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ; 
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }   
            result = var->to_bool( ct ) ;
         }
         else
         {
            result = arg(1)->to_bool( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_bool", name() ) ;
         break ;   
   }
   return( result ) ;
}
   
//----------------------------------------------------------------------
double
MAC_VariableExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE( ct ) ) ;

   double result = 0. ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }   
            result = var->to_double( ct ) ;
         }
         else
         {
            result = arg(1)->to_double( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_double", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
int
MAC_VariableExp:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE( ct ) ) ;

   int result = 0 ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_int( ct ) ;
         }
         else
         {
            result = arg(1)->to_int( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_int", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_VariableExp:: to_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_string" ) ;
   MAC_CHECK_PRE( to_string_PRE( ct ) ) ;

   static std::string result = "unexpected" ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_string( ct ) ;
         }
         else
         {
            result = arg(1)->to_string( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_string", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
boolVector const&
MAC_VariableExp:: to_bool_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_bool_vector" ) ;
   MAC_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;

   static boolVector result(0) ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_bool_vector( ct ) ;
         }
         else
         {
            result = arg(1)->to_bool_vector( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_bool_vector", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
doubleVector const&
MAC_VariableExp:: to_double_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_double_vector" ) ;
   MAC_CHECK_PRE( to_double_vector_PRE( ct ) ) ;

   static doubleVector result(0) ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_double_vector( ct ) ;
         }
         else
         {
            result = arg(1)->to_double_vector( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_double_vector", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
intVector const&
MAC_VariableExp:: to_int_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_int_vector" ) ;
   MAC_CHECK_PRE( to_int_vector_PRE( ct ) ) ;

   static intVector result(0) ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_int_vector( ct ) ;
         }
         else
         {
            result = arg(1)->to_int_vector( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_int_vector", name() ) ;
         break ;   
   }
   return( result ) ;
}
   
//----------------------------------------------------------------------
stringVector const&
MAC_VariableExp:: to_string_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_string_vector" ) ;
   MAC_CHECK_PRE( to_string_vector_PRE( ct ) ) ;

   static stringVector result(0) ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_string_vector( ct ) ;
         }
         else
         {
            result = arg(1)->to_string_vector( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_string_vector", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
MAC_VariableExp:: to_double_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_double_array2D" ) ;
   MAC_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;

   static doubleArray2D result(0,0) ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_double_array2D( ct ) ;
         }
         else
         {
            result = arg(1)->to_double_array2D( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_double_array2D", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
stringArray2D const&
MAC_VariableExp:: to_string_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_string_array2D" ) ;
   MAC_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;

   static stringArray2D result(0,0) ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_string_array2D( ct ) ;
         }
         else
         {
            result = arg(1)->to_string_array2D( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_string_array2D", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
boolArray2D const&
MAC_VariableExp:: to_bool_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_bool_array2D" ) ;
   MAC_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;

   static boolArray2D result(0,0) ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_bool_array2D( ct ) ;
         }
         else
         {
            result = arg(1)->to_bool_array2D( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_bool_array2D", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
intArray2D const&
MAC_VariableExp:: to_int_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VariableExp:: to_int_array2D" ) ;
   MAC_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;

   static intArray2D result(0,0) ;
   switch( OP )
   {
      case var_value :
      {
         MAC_Variable const* var =
                       MAC_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            MAC_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               MAC_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_int_array2D( ct ) ;
         }
         else
         {
            result = arg(1)->to_int_array2D( ct ) ;
         }
      }
      break ;
      default :
         MAC_VariableExp_ERROR::n0( "to_int_array2D", name() ) ;
         break ;   
   }
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
MAC_VariableExp_ERROR:: n0( std::string const& f_name,
                            std::string const& op_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** MAC_VariableExp::" + f_name +"\n" ;
   mesg += "    operation " + op_name + " not implemented." ;
   MAC_Error::object()->raise_internal( mesg ) ;
}

//internal--------------------------------------------------------------
void 
MAC_VariableExp_ERROR:: n1( std::string const& a_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** MAC_VariableExp:: value error\n" ;
   mesg += "    variable name \""+a_name+"\" and default value should have the same type." ;
   MAC_Error::object()->raise_plain( mesg ) ;
}

//internal--------------------------------------------------------------
void 
MAC_VariableExp_ERROR:: n2( std::string const& a_name,
                            stringVector const& undef_var )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** MAC_VariableExp:: value error\n" ;
   mesg += "    variable name \""+a_name+"\" cannot be evaluated.\n" ;
   if( undef_var.size() > 0 )
   {
      mesg += "    undefined variable(s):\n" ;
      for( size_t i=0 ; i<undef_var.size() ; ++i )
      {
         mesg += "       - \"" + undef_var(i) + "\"\n" ;
      }
   }
   MAC_Error::object()->raise_plain( mesg ) ;
}
