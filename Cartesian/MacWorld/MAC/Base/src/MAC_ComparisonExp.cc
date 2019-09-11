#include <MAC_ComparisonExp.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Sequence.hh>

#include <iostream>

MAC_ComparisonExp const*
MAC_ComparisonExp::PROTOTYPE_min = new MAC_ComparisonExp( "min", min_op ) ;

MAC_ComparisonExp const*
MAC_ComparisonExp::PROTOTYPE_max = new MAC_ComparisonExp( "max", max_op ) ;

//----------------------------------------------------------------------
MAC_ComparisonExp:: MAC_ComparisonExp(
   std::string const& a_name,
   Function a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( a_op )
   , FIRST( 0 ) 
   , SECOND( 0 )
{
   MAC_LABEL( "MAC_ComparisonExp:: MAC_ComparisonExp" ) ;
}

//----------------------------------------------------------------------
MAC_ComparisonExp:: MAC_ComparisonExp( MAC_Object* a_owner,
				       std::string const& a_name,
				       MAC_Sequence const* argument_list,
				       Function a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP(a_op)
   , FIRST( arg(0) ) 
   , SECOND( arg(1) )
{
   MAC_LABEL( "MAC_ComparisonExp:: MAC_ComparisonExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ComparisonExp:: ~MAC_ComparisonExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ComparisonExp:: ~MAC_ComparisonExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ComparisonExp*
MAC_ComparisonExp:: create_replica( MAC_Object* a_owner,
				    MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ComparisonExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_ComparisonExp* result = new MAC_ComparisonExp( a_owner, 
						      name(), 
						      argument_list, 
						      OP ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_ComparisonExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   result = name() + "(DS,DS) " ;
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_ComparisonExp:: valid_arguments(
                              MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ComparisonExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()==2 ;
   if( result )
   {
      result &=
         ( extract_arg( some_arguments, 0 )->data_type() == MAC_Data::Double )&&
         ( extract_arg( some_arguments, 1 )->data_type() == MAC_Data::Double ) ;
   }
   return result ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_ComparisonExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return Double ;
}

//----------------------------------------------------------------------
double
MAC_ComparisonExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ComparisonExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE( ct ) ) ;
   
   double result = MAC::bad_double() ;
   double v1 = FIRST->to_double( ct ) ;
   double v2 = SECOND->to_double( ct ) ;
   if( OP == min_op ) result = MAC::min( v1, v2 ) ;
   else if ( OP == max_op ) result = MAC::max( v1, v2 ) ;
   return( result ) ;
}
