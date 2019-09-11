#include <MAC_DerivativeExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Double.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_Error.hh>
#include <MAC_List.hh>
#include <MAC_Root.hh>
#include <MAC_Sequence.hh>
#include <MAC_Variable.hh>
#include <MAC.hh>

#include <iostream>

MAC_DerivativeExp const* 
MAC_DerivativeExp:: PROTOTYPE_d = new MAC_DerivativeExp( "d", d ) ;

MAC_DerivativeExp const* 
MAC_DerivativeExp:: PROTOTYPE_dnum = new MAC_DerivativeExp( "dnum", dnum ) ;

//----------------------------------------------------------------------
MAC_DerivativeExp:: MAC_DerivativeExp( std::string const& a_name,
                                       OP_TYPE a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( a_op )
   , EXP( 0 )
   , VAR( 0 )
   , DERIVATIVE( 0 )
{
   MAC_LABEL( "MAC_DerivativeExp:: MAC_DerivativeExp" ) ;
}

//----------------------------------------------------------------------
MAC_DerivativeExp*
MAC_DerivativeExp:: create_replica( MAC_Object* a_owner,
                                    MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DerivativeExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_DerivativeExp* result = new MAC_DerivativeExp( a_owner,
                                                      name(),
                                                      OP,
                                                      argument_list ) ;

   MAC_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DerivativeExp:: MAC_DerivativeExp( MAC_Object* a_owner,
                                       std::string const& a_name,
                                       OP_TYPE a_op,
                                       MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , EXP( arg(0) )
   , VAR( MAC_Variable::object( arg(1)->to_string() ) )
   , DERIVATIVE( 0 )
{
   MAC_LABEL( "MAC_DerivativeExp:: MAC_DerivativeExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_DerivativeExp:: ~MAC_DerivativeExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DerivativeExp:: ~MAC_DerivativeExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_DerivativeExp*
MAC_DerivativeExp:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DerivativeExp:: create_clone" ) ;
   MAC_CHECK_PRE( !is_a_prototype() ) ;
   
   MAC_DerivativeExp* result =
      static_cast<MAC_DerivativeExp*>(
         MAC_Expression::create_clone( a_owner ) ) ;

   if( !is_a_prototype() && DERIVATIVE!=0 )
      result->DERIVATIVE = DERIVATIVE->create_clone( result ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_DerivativeExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   result = name() + "(<expression>,SS)" ;
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_DerivativeExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DerivativeExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = some_arguments->count()==2 ;
   if( result )
   {
      MAC_Data::Type k0 =  extract_arg( some_arguments, 0 )->data_type() ;
      result = result && ( k0==Double || k0==DoubleVector ) ;
      MAC_Data::Type k1 =  extract_arg( some_arguments, 1 )->data_type() ;
      result = result && ( k1==String ) ;
   }
   return result ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_DerivativeExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DerivativeExp:: data_type" ) ;
   
   return EXP->data_type() ;
}

//----------------------------------------------------------------------
double
MAC_DerivativeExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DerivativeExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE(ct) ) ;

   double result ;
   
   if( OP==d )
   {
      result = derivative(ct)->to_double(ct) ;
   }
   else
   {
      MAC_ASSERT( OP==dnum ) ;
      static double const eps = 1.0e-8 ;
      
      double val = ct->value( VAR )->to_double(ct) ;
      MAC_ContextSimple* ctx = MAC_ContextSimple::create( 0 ) ;
      ctx->extend( ct ) ;
      double dx = MAC::max( eps, val*eps ) ;
      
      ctx->set_value_of( VAR, MAC_Double::create( ctx, val+dx ) ) ;
      result = ( EXP->to_double( ctx ) - EXP->to_double( ct ) )/dx ;
      ctx->destroy() ; ctx=0 ;
   }
   return result ;
}

//----------------------------------------------------------------------
doubleVector const&
MAC_DerivativeExp:: to_double_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DerivativeExp:: to_double_vector" ) ;
   MAC_CHECK_PRE( to_double_vector_PRE( ct ) ) ;

   MAC_ASSERT( OP==d ) ;

   doubleVector const& result = derivative(ct)->to_double_vector(ct) ;
   return result ;
}

//----------------------------------------------------------------------
MAC_Data*
MAC_DerivativeExp:: create_derivative( MAC_Object* a_owner,
                                       MAC_Variable const* var,
                                       MAC_Context const* ct ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DerivativeExp:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;

   if( OP!=d )
   {
      MAC_Error::object()->raise_plain( "Unable to differentiate "+name()+" operator" ) ;
   }
   
   MAC_Data* result = derivative(ct)->create_derivative( a_owner, var, ct ) ;
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_DerivativeExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space << name() << "( " ;
   EXP->print( os, 0 ) ;
   os << ", \"" << VAR->name() << "\" )" ;
}

//----------------------------------------------------------------------
MAC_Data*
MAC_DerivativeExp:: create_non_const_simplification(
                                              MAC_Object* a_owner ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DerivativeExp:: create_non_const_simplification" ) ;
   MAC_CHECK( create_non_const_simplification_PRE( a_owner ) ) ;
   MAC_Data* result = 0 ;
   
   if( OP==d )
   {
      MAC_ASSERT( DERIVATIVE!=0 ) ;
      result = DERIVATIVE->create_simplification( a_owner ) ;
   }
   
   MAC_CHECK( create_non_const_simplification_POST( a_owner, result ) ) ;
   return result ;
}

//----------------------------------------------------------------------
MAC_Data const*
MAC_DerivativeExp:: derivative( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DerivativeExp:: derivative" ) ;
   MAC_CHECK( OP==d ) ;
   
   if( DERIVATIVE==0 )
   {
      MAC_Data* der = EXP->create_derivative( 0, VAR, ct ) ;
//       std::cout << "Derivative of " ;
//       EX->print( std::cout, 1 ) ;
//       std::cout << std::endl << " with respect to "  ;
//       VAR->print( std::cout, 1 ) ;
//       std::cout << std::endl << " is : " ;
//       der->print( std::cout, 0 ) ;
//       std::cout << std::endl << " Simplified in " ;
      DERIVATIVE = der->create_simplification(
         const_cast<MAC_DerivativeExp*>(this) ) ;
//       DERIVATIVE->print( std::cout, 0 ) ;
//      std::cout << std::endl ;
      
      der->destroy() ;
   }
   return DERIVATIVE ;
}
