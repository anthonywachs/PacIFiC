#include <MAC_ArithmeticExp.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_Double.hh>
#include <MAC_Error.hh>
#include <MAC_List.hh>
#include <MAC_Sequence.hh>

#include <iostream>

MAC_ArithmeticExp const* 
MAC_ArithmeticExp::PROTOTYPE_M = new MAC_ArithmeticExp( "+", M ) ;

MAC_ArithmeticExp const* 
MAC_ArithmeticExp::PROTOTYPE_L = new MAC_ArithmeticExp( "-", L ) ;

MAC_ArithmeticExp const* 
MAC_ArithmeticExp::PROTOTYPE_T = new MAC_ArithmeticExp( "*", T ) ;

MAC_ArithmeticExp const* 
MAC_ArithmeticExp::PROTOTYPE_D = new MAC_ArithmeticExp( "/", D ) ;

MAC_ArithmeticExp const* 
MAC_ArithmeticExp::PROTOTYPE_MOD = new MAC_ArithmeticExp( "modulo", MOD ) ;

//----------------------------------------------------------------------
MAC_ArithmeticExp:: MAC_ArithmeticExp( std::string const& a_name,
                                       AlgebraicOperator a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( a_op )
   , ARG0( 0 )
   , ARG1( 0 )
{
   MAC_LABEL( "MAC_ArithmeticExp:: MAC_ArithmeticExp" ) ;
}

//----------------------------------------------------------------------
MAC_ArithmeticExp*
MAC_ArithmeticExp:: create_replica( MAC_Object* a_owner,
                                    MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArithmeticExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_ArithmeticExp* result = new MAC_ArithmeticExp( a_owner,
                                                      name(),
                                                      argument_list,
                                                      OP ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ArithmeticExp:: MAC_ArithmeticExp( MAC_Object* a_owner,
                                       std::string const& a_name,
                                       MAC_Sequence const* argument_list,
                                       AlgebraicOperator a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , ARG0( arg(0) )
   , ARG1( arg(1) )
{
   MAC_LABEL( "MAC_ArithmeticExp:: MAC_ArithmeticExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ArithmeticExp:: ~MAC_ArithmeticExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArithmeticExp:: ~MAC_ArithmeticExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_ArithmeticExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   if( OP==M )
   {
      result = "IS|DS|SS " + name() + " <same type> " ;
   }
   else if( OP==MOD )
   {
      result = "modulo(IS,IS)" ;
   }
   else 
   {
      result = "IS|DS " + name() + " <same type> " ;
   }
       
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_ArithmeticExp:: valid_arguments(
                              MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArithmeticExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = some_arguments->count()==2 ;
   if( result )
   {
      MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
      MAC_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
      result = ( k0==k1 ) ;
      if( OP==MOD )
      {
         result &= ( k0 == MAC_Data::Int ) ;
      }
      else if( OP==M )
      {
         result &= ( k0 == MAC_Data::Int ||
                     k0 == MAC_Data::Double ||
                     k0 == MAC_Data::String )  ;
      }
      else
      {
         result &= ( k0 == MAC_Data::Int ||
                     k0 == MAC_Data::Double )  ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_ArithmeticExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return arg(0)->data_type() ;
}

//----------------------------------------------------------------------
double
MAC_ArithmeticExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArithmeticExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE( ct ) ) ;
   
   double result = MAC::max_double() ;
   double v1 = ARG0->to_double( ct ) ;
   double v2 = ARG1->to_double( ct ) ;
   switch( OP )
   {
      case M : result = v1 + v2 ;
         break ;
      case L : result = v1 - v2 ;
         break ;
      case T : result = v1 * v2 ;
         break ;
      case D : result = v1 / v2 ;
         break ;
      default :
         MAC_Error::object()->raise_plain( "Internal error" ) ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
int
MAC_ArithmeticExp:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArithmeticExp:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE(ct) ) ;
   
   int result = MAC::max_int() ;
   int v1 = ARG0->to_int( ct ) ;
   int v2 = ARG1->to_int( ct ) ;

   switch( OP )
   {
      case M : result = v1 + v2 ;
         break ;
      case L : result = v1 - v2 ;
         break ;
      case T : result = v1 * v2 ;
         break ;
      case D : result = v1 / v2 ;
         break ;
      case MOD : result = v1 % v2 ;
         break ;
      default :
         MAC_Error::object()->raise_plain( "Internal error" ) ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
std::string const&
MAC_ArithmeticExp:: to_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArithmeticExp:: to_string" ) ;
   MAC_CHECK_PRE( to_string_PRE(ct) ) ;
   MAC_ASSERT( name()=="+" ) ;
   
   std::string const& v1 = ARG0->to_string( ct ) ;
   std::string const& v2 = ARG1->to_string( ct ) ;
   RESULT_STR = v1 + v2 ;
   return RESULT_STR ;
}

//----------------------------------------------------------------------
MAC_Data*
MAC_ArithmeticExp:: create_derivative( MAC_Object* a_owner,
                                       MAC_Variable const* var,
                                       MAC_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArithmeticExp:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   MAC_Data* result = 0 ;
   MAC_List* list = MAC_List::create( 0 ) ;
   MAC_Data* d1 = ARG0->create_derivative( list, var, ct ) ;
   MAC_Data* d2 = ARG1->create_derivative( list, var, ct ) ;
   
   if( OP==M || OP==L )
   {
      list->append(d1) ;
      list->append(d2) ;
      
      result = MAC_Expression::create( a_owner, name(), list ) ;
      list->set_owner( result ) ;
   }
   else if( OP==T )
   {
      list->append(ARG0->create_clone(list)) ;
      list->append(d2) ; 
      MAC_Expression* m1 = MAC_Expression::create( 0, "*", list ) ;
      list->set_owner( m1 ) ;
      
      list = MAC_List::create( 0 ) ;
      list->append(d1) ;
      list->append(ARG1->create_clone(list)) ; 
      MAC_Expression* m2 = MAC_Expression::create( 0, "*", list ) ;
      list->set_owner( m2 ) ;
      
      list = MAC_List::create( 0 )  ;
      list->append(m1) ; m1->set_owner( list ) ;
      list->append(m2) ; m2->set_owner( list ) ;
      result = MAC_Expression::create( a_owner, "+", list ) ;
      list->set_owner( result ) ;
   }
   else if( OP==D )
   {
      list->append(d1) ;
      list->append(ARG1->create_clone(list)) ; 
      MAC_Expression* m1 = MAC_Expression::create( 0, "*", list ) ;
      list->set_owner( m1 ) ;
      
      list = MAC_List::create( 0 ) ;
      list->append(ARG0->create_clone(list)) ;
      list->append(d2) ; 
      MAC_Expression* m2 = MAC_Expression::create( 0, "*", list ) ;
      list->set_owner( m2 ) ;
      
      list = MAC_List::create( 0 ) ;
      list->append(m1) ; m1->set_owner( list ) ;
      list->append(m2) ; m2->set_owner( list ) ;
      MAC_Expression* num = MAC_Expression::create( 0, "-", list ) ;
      list->set_owner( num ) ;
      
      list = MAC_List::create( 0 ) ;
      list->append(ARG1->create_clone(list)) ;
      MAC_Expression* den = MAC_Expression::create( 0, "sqr", list ) ;
      list->set_owner( den ) ;

      list = MAC_List::create( 0 ) ;
      list->append(num) ; num->set_owner( list ) ;
      list->append(den) ; den->set_owner( list ) ;
      result = MAC_Expression::create( a_owner, "/", list ) ;
      list->set_owner( result ) ;
   }
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_ArithmeticExp:: print( std::ostream& os, size_t indent_width ) const
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

//----------------------------------------------------------------------
MAC_Data*
MAC_ArithmeticExp:: create_operator_simplification( MAC_Object* a_owner )  
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArithmeticExp:: create_operator_simplification" ) ;

   MAC_Data* result = this ;
   
   bool first_const = ARG0->is_constant() ;
   bool second_const = ARG1->is_constant() ;
   MAC_CHECK( ! ( first_const && second_const ) ) ;
   if( first_const || second_const )
   {
      double v = ( first_const ? ARG0->to_double() : ARG1->to_double() ) ;
      MAC_Data const* non_const = ( first_const ? ARG1 : ARG0 ) ;
      bool null = v==0.0 ;
      bool unity = v==1.0 ;
      if( OP==M && null )
      {
         result = non_const->create_clone( a_owner ) ;
      }
      else if( OP==L && null && second_const )
      {
         result = non_const->create_clone( a_owner ) ;   
      }
      else if( OP==T && null )
      {
         result = MAC_Double::create( a_owner, 0.0 ) ;
      }
      else if( OP==T && unity )
      {
         result = non_const->create_clone( a_owner ) ;
      }
      else if( OP==D && null )
      {
         if( first_const )
         {
            result = MAC_Double::create( a_owner, 0.0 ) ;
         }
         else
         {
            MAC_Error::object()->raise_plain(
               "When simplifiing expression, null dividend found" ) ;
         }
      }
      else if( OP==D && unity &&  second_const )
      {
         
         result = non_const->create_clone( a_owner ) ;
      }
   }

   MAC_CHECK_POST( create_operator_simplification_POST( a_owner, result ) ) ;
   return result ;
}
