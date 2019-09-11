#include <MAC_BooleanExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Sequence.hh>

#include <iostream>

MAC_BooleanExp const*  
MAC_BooleanExp::PROTOTYPE_or  = new MAC_BooleanExp( OR, "||" ) ;

MAC_BooleanExp const*
MAC_BooleanExp::PROTOTYPE_and = new MAC_BooleanExp( AND, "&&" ) ;

MAC_BooleanExp const*
MAC_BooleanExp::PROTOTYPE_not = new MAC_BooleanExp( NOT, "!" ) ;

//----------------------------------------------------------------------
MAC_BooleanExp:: MAC_BooleanExp( BoolExp exp_id, 
                                 std::string const& a_name  ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( exp_id )
   , ARG0( 0 )
   , ARG1( 0 )
{
   MAC_LABEL( "MAC_BooleanExp:: MAC_BooleanExp" ) ;
}

//----------------------------------------------------------------------
MAC_BooleanExp*
MAC_BooleanExp:: create_replica( MAC_Object* a_owner,
                                 MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BooleanExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_BooleanExp* result = new MAC_BooleanExp( a_owner, 
                                                OP,
                                                name(),
                                                argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_BooleanExp:: MAC_BooleanExp( MAC_Object* a_owner,
                                 BoolExp exp_id,
                                 std::string const& a_name,
                                 MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list ),
     OP( exp_id ),
     ARG0( arg(0) ),
     ARG1( nb_arguments()>1 ? arg(1) : 0 )
{
   MAC_LABEL( "MAC_BooleanExp:: MAC_BooleanExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_BooleanExp:: ~MAC_BooleanExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BooleanExp:: ~MAC_BooleanExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_BooleanExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   switch( OP )
   {
      case OR :
         result = "BS || BS" ;
	 break ;
      case AND :
         result = "BS && BS" ;
	 break ;
      case NOT :
         result = "! BS" ;
	 break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_BooleanExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BooleanExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = false ;
   switch( OP )
   {
      case OR  :
      case AND :
         result = some_arguments->count()==2 ;
	 if( result )
	 {
	    Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
	    Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
	    result = result && ( t0==Bool && t1==Bool ) ;
	 }
	 break ;
      case NOT :
         result = ( some_arguments->count() == 1 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = result && ( t0 == Bool ) ;
         }
	 break ;
      default : result = false ;
   }
   return result ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_BooleanExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return MAC_Data::Bool ;
}

//----------------------------------------------------------------------
bool
MAC_BooleanExp:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BooleanExp:: to_bool" ) ;
   MAC_CHECK_PRE( to_bool_PRE( ct ) ) ;

   bool result ;
   switch( OP )
   {
      case AND :
         result = ( ARG0->to_bool(ct) && ARG1->to_bool(ct) ) ;
         break ;
      case OR :
         result = ( ARG0->to_bool(ct) || ARG1->to_bool(ct) ) ;
         break ;
      case NOT :
         result = !( ARG0->to_bool( ct ) ) ;
         break ;
      default :
	MAC_Error::object()->raise_internal( "Bad operator" ) ;
	result = false ;
   }
   return result ;
}

//----------------------------------------------------------------------
void
MAC_BooleanExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space ;
   if( external_brackets_are_set() ) os << "(" ;
   switch( OP )
   {
      case AND :
         ARG0->print( os, 0 ) ;
         os << " && " ;
         ARG1->print( os, 0 ) ;
         break ;
      case OR :
         ARG0->print( os, 0 ) ;
         os << " || " ;
         ARG1->print( os, 0 ) ;
         break ;
      case NOT :
         os << "! " ;
         ARG0->print( os, 0 ) ;
         break ;
      default :
	MAC_Error::object()->raise_internal( "Bad operator" ) ;
   }
   if( external_brackets_are_set() ) os << ")" ;
}
