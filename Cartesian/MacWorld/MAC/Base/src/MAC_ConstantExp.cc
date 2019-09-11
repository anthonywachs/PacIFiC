#include <MAC_ConstantExp.hh>

#include <MAC.hh>

#include <MAC_assertions.hh>
#include <MAC_Double.hh>
#include <MAC_Error.hh>
#include <MAC_Sequence.hh>

MAC_ConstantExp const* 
MAC_ConstantExp::PROTOTYPE_pi = new MAC_ConstantExp( "pi", MAC::pi() ) ;

MAC_ConstantExp const* 
MAC_ConstantExp::PROTOTYPE_e = new MAC_ConstantExp( "e", MAC::e() ) ;

MAC_ConstantExp const* 
MAC_ConstantExp::PROTOTYPE_eu = new MAC_ConstantExp( "euler", MAC::euler() ) ;

//----------------------------------------------------------------------
MAC_ConstantExp:: MAC_ConstantExp( std::string const& a_name,
                                   double a_val) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , VAL( a_val )
{
   MAC_LABEL( "MAC_ConstantExp:: MAC_ConstantExp" ) ;
}

//----------------------------------------------------------------------
MAC_ConstantExp*
MAC_ConstantExp:: create_replica( MAC_Object* a_owner,
                                  MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConstantExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_ConstantExp* result =
                  new MAC_ConstantExp( a_owner, name(), VAL, argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ConstantExp:: MAC_ConstantExp( MAC_Object* a_owner,
                                   std::string const& a_name,
                                   double a_val,
                                   MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , VAL( a_val )
{
   MAC_LABEL( "MAC_ConstantExp:: MAC_ConstantExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ConstantExp:: ~MAC_ConstantExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConstantExp:: ~MAC_ConstantExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_ConstantExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConstantExp:: data_type" ) ;
   MAC_Data::Type result = Double ;
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_ConstantExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConstantExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = some_arguments->count()==0 ;
   return result ;
}

//----------------------------------------------------------------------
std::string const&
MAC_ConstantExp:: usage( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConstantExp:: usage" ) ;
   static std::string result ;
   result = name()+"()" ;
   return result ;
}

//----------------------------------------------------------------------
double
MAC_ConstantExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConstantExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE(ct) ) ;
   return VAL ;
}

//----------------------------------------------------------------------
MAC_Data*
MAC_ConstantExp:: create_derivative( MAC_Object* a_owner,
                                     MAC_Variable const* var,
                                     MAC_Context const* ct ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConstantExp:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   MAC_Data* result = MAC_Double::create( a_owner, 0.0 ) ;
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

