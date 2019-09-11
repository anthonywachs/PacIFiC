#include <MAC_UnaryArithmeticExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_List.hh>
#include <MAC_Root.hh>
#include <MAC_Sequence.hh>
#include <MAC.hh>

#include <iostream>

//----------------------------------------------------------------------
MAC_UnaryArithmeticExp const* 
MAC_UnaryArithmeticExp:: PROTOTYPE_Minus = 
                          new MAC_UnaryArithmeticExp( "unary_minus", Minus ) ;
//----------------------------------------------------------------------

//----------------------------------------------------------------------
MAC_UnaryArithmeticExp:: MAC_UnaryArithmeticExp( std::string const& a_name,
                                                 UnaryOperator a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( a_op )
   , FIRST( 0 )
{
   MAC_LABEL( "MAC_UnaryArithmeticExp:: MAC_UnaryArithmeticExp" ) ;
}

//----------------------------------------------------------------------
MAC_UnaryArithmeticExp*
MAC_UnaryArithmeticExp:: create_replica( MAC_Object* a_owner,
                             MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_UnaryArithmeticExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_UnaryArithmeticExp* result = new MAC_UnaryArithmeticExp( a_owner, 
                                                          name(), 
                                                          argument_list, 
                                                          OP ) ;

   MAC_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_UnaryArithmeticExp:: MAC_UnaryArithmeticExp( MAC_Object* a_owner,
                                     std::string const& a_name,
                                     MAC_Sequence const* argument_list,
                                     UnaryOperator a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , FIRST( arg(0) )
{
   MAC_LABEL( "MAC_UnaryArithmeticExp:: MAC_UnaryArithmeticExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_UnaryArithmeticExp:: ~MAC_UnaryArithmeticExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_UnaryArithmeticExp:: ~MAC_UnaryArithmeticExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_UnaryArithmeticExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_UnaryArithmeticExp:: data_type" ) ;
   
   return FIRST->data_type() ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_UnaryArithmeticExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = name() + "<double|integer>" ;
   
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_UnaryArithmeticExp:: valid_arguments(
                              MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_UnaryArithmeticExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()==1 ;
   if( result )
   {
      MAC_Data::Type k0 =  extract_arg( some_arguments, 0 )->data_type() ;
      result = result && ( k0==Double || k0==Int ) ;
   }
   return result ;
}

//----------------------------------------------------------------------
double
MAC_UnaryArithmeticExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_UnaryArithmeticExp:: to_double" ) ;
   
   MAC_CHECK_PRE( to_double_PRE(ct) ) ;
   double result = MAC::max_double() ;
   double v = FIRST->to_double(ct) ;
   switch(OP)
   {
      case Minus : result = -v ;
         break ;
      default :
         MAC_Error::object()->raise_plain( "Internal error" ) ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
int
MAC_UnaryArithmeticExp:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_UnaryArithmeticExp:: to_int" ) ;
   
   MAC_CHECK_PRE( to_int_PRE(ct) ) ;
   int result = MAC::max_int() ;
   int v = FIRST->to_int(ct) ;
   switch(OP)
   {
      case Minus : result = -v ;
         break ;
      default :
         MAC_Error::object()->raise_plain( "Internal error" ) ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
MAC_Data*
MAC_UnaryArithmeticExp:: create_derivative( MAC_Object* a_owner,
                                            MAC_Variable const* var,
                                            MAC_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_UnaryArithmeticExp:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;

   MAC_List* list = MAC_List::create( 0 ) ;
   list->append( FIRST->create_derivative( list, var, ct ) ) ;
   
   MAC_Data* result = create_replica( a_owner, list ) ;
   list->set_owner( result ) ;
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_UnaryArithmeticExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space << "-" ;
   if( external_brackets_are_set() ) os << "(" ;
   FIRST->print( os, 0 ) ;
   if( external_brackets_are_set() ) os << ")" ;
}
