#include <MAC_DataWithContextExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Context.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_Error.hh>
#include <MAC_List.hh>
#include <MAC_String.hh>
#include <MAC_Variable.hh>
#include <stringVector.hh>

#include <iostream>

MAC_DataWithContextExp const*
MAC_DataWithContextExp::PROTO =
                     new MAC_DataWithContextExp( "data_with_context" ) ;

//----------------------------------------------------------------------
MAC_DataWithContextExp:: MAC_DataWithContextExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : MAC_TransferExp( a_name )
   , DATA( 0 )
{
   MAC_LABEL( "MAC_DataWithContextExp:: MAC_DataWithContextExp" ) ;
}

//----------------------------------------------------------------------
MAC_DataWithContextExp:: MAC_DataWithContextExp(
                                    MAC_Object* a_owner,
                                    std::string const& a_name,
                                    MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_TransferExp( a_owner, a_name, argument_list )
   , DATA( 0 )
{
   MAC_LABEL( "MAC_DataWithContextExp:: MAC_DataWithContextExp" ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Data const* d = arg(0) ;
   MAC_ContextSimple* ct = MAC_ContextSimple::create( this ) ;
   for( size_t i=1 ; i<argument_list->count() ; )
   {
      std::string const& n = arg(i++)->to_string() ;
      MAC_Variable const* v = MAC_Variable::object( n ) ;
      ct->extend( v, arg(i++)->create_clone( ct ) ) ;
   }
   DATA = MAC_DataWithContext::create( this, d, ct ) ;
}

//----------------------------------------------------------------------
MAC_DataWithContextExp:: ~MAC_DataWithContextExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp:: ~MAC_DataWithContextExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_DataWithContextExp*
MAC_DataWithContextExp:: create(
                           MAC_Object* a_owner,
                           MAC_Data const* data, MAC_Context const* ct )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp:: create" ) ;
   MAC_CHECK_PRE( data != 0 ) ;

   MAC_List* argument_list = MAC_List::create( 0 ) ;
   argument_list->append( data->create_clone( argument_list ) ) ;
   if( ct != 0 )
   {
      for( size_t i=0 ; i<ct->nb_variables() ; ++i )
      {
         MAC_Variable const* v = ct->variable(i) ;
         argument_list->append( MAC_String::create( argument_list, v->name() ) ) ;
         argument_list->append( ct->value( v )->create_clone( argument_list ) ) ;
      }
   }
   
   MAC_DataWithContextExp* result =
                       PROTO->create_replica( a_owner, argument_list ) ;

   argument_list->set_owner( result ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DataWithContextExp*
MAC_DataWithContextExp:: create_replica(
                             MAC_Object* a_owner,
                             MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_DataWithContextExp* result =
      new MAC_DataWithContextExp( a_owner, name(), argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_DataWithContextExp:: declare( MAC_List* lst ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp:: declare" ) ;
   MAC_CHECK_PRE( declare_PRE( lst ) ) ;
   DATA->declare( lst ) ;
   MAC_CHECK_POST( declare_POST( lst ) ) ;
}

//----------------------------------------------------------------------
bool
MAC_DataWithContextExp:: context_has_required_variables(
                                           MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp::context_has_required_variables " ) ;
   MAC_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;
   return( DATA->context_has_required_variables( ct ) ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_DataWithContextExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp:: data_type" ) ;
   return( DATA->data_type() ) ;
}

//----------------------------------------------------------------------
bool
MAC_DataWithContextExp:: is_constant( void ) const
//----------------------------------------------------------------------
{
   return( false ) ;
}

//----------------------------------------------------------------------
bool
MAC_DataWithContextExp:: value_can_be_evaluated( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp:: value_can_be_evaluated" ) ;
   return( DATA->value_can_be_evaluated( ct ) ) ;
}

//----------------------------------------------------------------------
stringVector const& 
MAC_DataWithContextExp:: undefined_variables( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp:: undefined_variables" ) ;
   return( DATA->undefined_variables( ct ) ) ;
}

//----------------------------------------------------------------------
bool
MAC_DataWithContextExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   size_t const n = some_arguments->count() ;
   bool result = ( n%2 == 1 ) ;
   for( size_t i=1 ; result && i<n ; i += 2 )
   {
      result = ( extract_arg( some_arguments, i )->data_type() == String ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_DataWithContextExp:: usage( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp:: usage" ) ;
   static std::string result =
      name() + "(<expression>[,SS,<value>])" ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data const*
MAC_DataWithContextExp:: data( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContextExp:: data" ) ;
   MAC_CHECK( data_PRE( ct ) ) ;

   MAC_Data const* result = DATA ;

   MAC_CHECK_POST( data_POST( result, ct ) ) ;
   return( result ) ;
}
