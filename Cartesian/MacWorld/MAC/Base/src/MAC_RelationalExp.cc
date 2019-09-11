#include <MAC_RelationalExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Sequence.hh>

#include <iostream>

MAC_RelationalExp const* 
MAC_RelationalExp::PROTOTYPE_LE = new MAC_RelationalExp( "<=", LE ) ;

MAC_RelationalExp const* 
MAC_RelationalExp::PROTOTYPE_GE = new MAC_RelationalExp( ">=", GE ) ;

MAC_RelationalExp const* 
MAC_RelationalExp::PROTOTYPE_LT = new MAC_RelationalExp( "<", LT ) ;

MAC_RelationalExp const*
MAC_RelationalExp::PROTOTYPE_GT = new MAC_RelationalExp( ">", GT ) ;

MAC_RelationalExp const* 
MAC_RelationalExp::PROTOTYPE_EQ = new MAC_RelationalExp( "=", EQ ) ;

MAC_RelationalExp const* 
MAC_RelationalExp::PROTOTYPE_NEQ = new MAC_RelationalExp( "!=", NEQ ) ;

//----------------------------------------------------------------------
MAC_RelationalExp*
MAC_RelationalExp:: create_replica( MAC_Object* a_owner,
                             MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_RelationalExp:: create_replica" ) ;
   return new MAC_RelationalExp( a_owner, name(), argument_list, OP ) ;
}

//----------------------------------------------------------------------
MAC_RelationalExp:: MAC_RelationalExp( std::string const& a_name,
                                       ComparisonOperator a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP(a_op)
   , ARG0(0)
   , ARG1(0)
{
   MAC_LABEL( "MAC_RelationalExp:: MAC_RelationalExp" ) ;
}

//----------------------------------------------------------------------
MAC_RelationalExp:: MAC_RelationalExp( MAC_Object* a_owner,
                                      std::string const& a_name,
                                      MAC_Sequence const* argument_list,
                                      ComparisonOperator a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP(a_op)
   , ARG0( arg(0) )
   , ARG1( arg(1) )
{
   MAC_LABEL( "MAC_RelationalExp:: MAC_RelationalExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_RelationalExp:: ~MAC_RelationalExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_RelationalExp:: ~MAC_RelationalExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_RelationalExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return( MAC_Data::Bool ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_RelationalExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   if( OP == EQ || OP == NEQ )
   {
      result = "IS|DS|SS|BS " + name() + " <same type>" ;
   }
   else
   {
      result = "IS|DS " + name() + " <same type>" ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_RelationalExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_RelationalExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()==2 ;
   if( result )
   {
      MAC_Data::Type k0 =  extract_arg( some_arguments, 0 )->data_type() ;
      MAC_Data::Type k1 =  extract_arg( some_arguments, 1 )->data_type() ;
      result = ( k0==k1 ) ;
      if( result )
      {
         if( OP == EQ || OP == NEQ )
         {
            result = k0==Double || k0==Int || k0==String || k0== Bool ;
         }
         else
         {
            result = k0==Double || k0==Int ;
         }
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_RelationalExp:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_RelationalExp:: to_bool" ) ;
   
   MAC_CHECK_PRE( to_bool_PRE( ct ) ) ;
   bool result = true ;
   MAC_Data::Type k0 = ARG0->data_type() ;
   
   if( OP==EQ && k0==String )
   {
      result = ( ARG0->to_string( ct ) == ARG1->to_string( ct ) ) ;
   }
   else if( OP==EQ && k0==Bool )
   {
      result = ( ARG0->to_bool( ct ) == ARG1->to_bool( ct ) ) ;
   }
   else if( OP==NEQ && k0==String )
   {
      result = ( ARG0->to_string( ct ) != ARG1->to_string( ct ) ) ;
   }
   else if( OP==NEQ && k0==Bool )
   {
      result = ( ARG0->to_bool( ct ) != ARG1->to_bool( ct ) ) ;
   }
   else
   {
      double v1 ;
      double v2 ;
      if( ARG0->data_type()==Double )
      {
         v1 = ARG0->to_double( ct ) ;
      }
      else
      {
         v1 = ARG0->to_int( ct ) ;
      }
      if( ARG1->data_type()==Double )
      {
         v2 = ARG1->to_double( ct ) ;
      }
      else
      {
         v2 = ARG1->to_int( ct ) ;
      }
      switch(OP)
      {
         case LE : result = ( v1 <= v2 ) ;
            break ;
         case GE : result = ( v1 >= v2 ) ;
            break ;
         case LT : result = ( v1 < v2 ) ;
            break ;
         case GT : result = ( v1 > v2 ) ;
            break ;
         case EQ : result = ( v1 == v2 ) ;
            break ;
         case NEQ : result = ( v1 != v2 ) ;
            break ;
         default :
            MAC_Error::object()->raise_plain( "Internal error" ) ;
      }
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_RelationalExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space ;
   if( external_brackets_are_set() ) os << "(" ;
   ARG0->print( os, 0 ) ;
   os << name() ;
   ARG1->print( os, 0 ) ;
   if( external_brackets_are_set() ) os << ")" ;
}
