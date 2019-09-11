#include <MAC_ContextSimple.hh>

#include <MAC_assertions.hh>
#include <MAC_Iterator.hh>
#include <MAC_List.hh>
#include <MAC_Vector.hh>
#include <MAC_Variable.hh>

//----------------------------------------------------------------------
MAC_ContextSimple*
MAC_ContextSimple:: create( MAC_Object* a_owner ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextSimple:: create" ) ;
   MAC_ContextSimple* result = new MAC_ContextSimple( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ContextSimple:: MAC_ContextSimple( MAC_Object* a_owner ) 
//----------------------------------------------------------------------
   : MAC_Context( a_owner )
   , VALUES( MAC_Vector::create( this, 0 ) )
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ContextSimple*
MAC_ContextSimple:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextSimple:: create_clone" ) ;

   MAC_ContextSimple* result = new MAC_ContextSimple( a_owner ) ;
   
   MAC_Iterator* it = VALUES->create_iterator( 0 ) ;
   for( size_t i=0 ; i<VALUES->index_limit() ; i++  )
   {
      MAC_Data const* dat = static_cast<MAC_Data const*>( VALUES->at(i) ) ;
      if( dat!=0 )
      {
         MAC_Variable const*var=MAC_Variable::object(i) ;
         result->extend( var, dat->create_clone( result ) ) ;
      }
   }
   it->destroy() ;

   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ContextSimple:: ~MAC_ContextSimple( void ) 
//----------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;

   notify_observers_of_my_destruction() ;
}


//----------------------------------------------------------------------
size_t
MAC_ContextSimple:: nb_variables( void ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextSimple:: nb_variables" ) ;
   return( VALUES->count() ) ;
}

//----------------------------------------------------------------------
MAC_Variable const* 
MAC_ContextSimple:: variable( size_t i ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextSimple:: variable" ) ;
   MAC_CHECK_PRE( variable_PRE(i) ) ;
   
   MAC_Variable const* result = 0 ;
   for( size_t j=0 ; j<VALUES->index_limit() ; j++ )
   {
      if( VALUES->at(j)!=0 )
      {
         if( i==0 )
         {
            result = MAC_Variable::object(j) ;
            break ;
         }
         i-- ;
      }
   }
   
   MAC_CHECK_POST( variable_POST( result ) ) ;
   return( result ) ;
      
}

//----------------------------------------------------------------------
bool
MAC_ContextSimple:: has_variable( MAC_Variable const* var ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextSimple:: has_variable" ) ;
   MAC_CHECK_PRE( has_variable_PRE( var ) ) ;

   size_t idx = var->id_number() ;
   bool result = ( VALUES->index_limit()>idx && VALUES->at( idx )!=0 ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data* 
MAC_ContextSimple:: value( MAC_Variable const* var ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextSimple:: value" ) ;
   MAC_CHECK_PRE( value_PRE( var ) ) ;
   
   MAC_Data* result = static_cast<MAC_Data* >(
      VALUES->at( var->id_number() ) ) ;

   MAC_CHECK_POST( value_POST( result, var ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_ContextSimple:: has_circular_definition( MAC_Variable const* var,
                                             MAC_Data const* a_value ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextSimple:: has_circular_definition" ) ;
   MAC_CHECK_PRE( var != 0 ) ;
   MAC_CHECK_PRE( a_value != 0 ) ;
   
   MAC_List* dummy_context = MAC_List::create( 0 ) ;
   a_value->declare( dummy_context ) ;
   bool result = dummy_context->has( var ) ;
   dummy_context->destroy() ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_ContextSimple:: extend( MAC_Variable const* var,
                            MAC_Data const* a_value )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextSimple:: extend" ) ;
   MAC_CHECK_PRE( var != 0 ) ;
   MAC_CHECK_PRE( a_value != 0 ) ;
   MAC_CHECK_PRE( !has_circular_definition( var, a_value ) ) ;
   MAC_CHECK_PRE( var->data_type()==a_value->data_type() ) ;
   MAC_CHECK_PRE( a_value->is_under_ownership_of( this ) ) ;
   
   size_t idx = var->id_number() ;
   if( VALUES->index_limit()<=idx )
   {
      VALUES->resize( idx+10 ) ;
   }
   VALUES->set_at( idx, const_cast<MAC_Data*>(a_value) ) ;
   update_observers() ;
   MAC_CHECK_POST( has_variable( var ) ) ;
   MAC_CHECK_POST( value(var) == a_value ) ;
}

//----------------------------------------------------------------------
void
MAC_ContextSimple:: extend( MAC_Context const* other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextSimple:: extend( MAC_Context const* )" ) ;
   MAC_CHECK_PRE( other != 0 ) ;
   
   for( size_t i=0 ; i<other->nb_variables() ; i++ )
   {
      MAC_Variable const* var = other->variable(i) ;
      extend( other->variable(i), other->value(var)->create_clone( this ) ) ;
   }
   update_observers() ;
}

//----------------------------------------------------------------------
void
MAC_ContextSimple:: set_value_of( MAC_Variable const* var,
                                  MAC_Data const* a_value )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextSimple:: set_value_of" ) ;
   MAC_CHECK_PRE( var !=0 ) ;
   MAC_CHECK_PRE( has_variable( var ) ) ;
   MAC_CHECK_PRE( a_value !=0 ) ;
   MAC_CHECK_PRE( var->data_type()==a_value->data_type() ) ;
   MAC_CHECK_PRE( !has_circular_definition(var,a_value ) ) ;
   MAC_CHECK_PRE( a_value->is_under_ownership_of( this ) ) ;
   
   VALUES->set_at( var->id_number(), const_cast<MAC_Data*>(a_value) ) ;
   update_observers() ;

   MAC_CHECK_POST( value( var )==a_value ) ;
   
}

//----------------------------------------------------------------------
void
MAC_ContextSimple:: update( void )
//----------------------------------------------------------------------
{
   update_observers() ;
}

//----------------------------------------------------------------------
void
MAC_ContextSimple:: update_for_destruction_of( MAC_Context const* subject )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
MAC_ContextSimple:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( VALUES!=0 ) ;
   return( true ) ;
}




