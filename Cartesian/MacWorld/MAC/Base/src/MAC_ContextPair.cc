#include <MAC_ContextPair.hh>

#include <MAC_assertions.hh>
#include <MAC_Vector.hh>
#include <MAC_Variable.hh>

//----------------------------------------------------------------------
MAC_ContextPair*
MAC_ContextPair:: create( MAC_Object* a_owner,
                          MAC_Context const* first,
                          MAC_Context const* second ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextPair:: create" ) ;
   MAC_CHECK_PRE( EQUIVALENT( first == 0, second == 0 ) ) ;

   MAC_ContextPair* result = new MAC_ContextPair( a_owner, first, second  ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( IMPLIES( first == 0 || second == 0,
                            result->nb_variables() == 0 ) ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; first != 0 && i<first->nb_variables() ; ++i ),
         result->has_variable( first->variable( i ) ) ) ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; second != 0 && i<second->nb_variables() ; ++i ),
         result->has_variable( second->variable( i ) ) ) ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<result->nb_variables() ; ++i ),
         first->has_variable( result->variable( i ) ) ||
         second->has_variable( result->variable( i ) ) ) ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<result->nb_variables() ; ++i ),
         IMPLIES(
            second->has_variable( result->variable( i ) ),
            result->value( result->variable( i ) ) == second->value( result->variable( i ) ) ) ) ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<result->nb_variables() ; ++i ),
         IMPLIES(
            !second->has_variable( result->variable( i ) ),
            first->has_variable( result->variable( i ) ) &&
            result->value( result->variable( i ) ) == first->value( result->variable( i ) ) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ContextPair:: MAC_ContextPair( MAC_Object* a_owner,
                                   MAC_Context const* first,
                                   MAC_Context const* second ) 
//----------------------------------------------------------------------
   : MAC_Context( a_owner )
   , CT1( 0 )
   , CT2( 0 )
   , VALUES( MAC_Vector::create( this, MAC_Variable::nb_objects() ) )
{
   MAC_CHECK( EQUIVALENT( first != 0, second != 0 ) ) ;
   if( first != 0 || second != 0 )
   {
      re_initialize( first, second ) ;
   }
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
MAC_ContextPair:: re_initialize( MAC_Context const* first,
                                 MAC_Context const* second )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextPair:: re_initialize" ) ;
   MAC_CHECK_PRE( first != 0 && second != 0 ) ;
   MAC_CHECK_PRE( first != second ) ;
   MAC_CHECK_PRE( first != this ) ;
   MAC_CHECK_PRE( this != second ) ;
   
   if( CT1!=first || CT2!=second )
   {
      if( (CT1!=0) && (CT1!=first) )
      {
         CT1->detach_observer( this ) ;
      }
      if( (CT2!=0) && (CT2!=second) )
      {
         CT2->detach_observer( this ) ;
      }
      CT1 = first ;
      CT2 = second ;
      if( CT1!=0 )
      {
         CT1->attach_observer( this ) ;
      }
      if( CT2!=0 )
      {
         CT2->attach_observer( this ) ;
      }
      update() ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<first->nb_variables() ; ++i ),
         has_variable( first->variable( i ) ) ) ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<second->nb_variables() ; ++i ),
         has_variable( second->variable( i ) ) ) ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<nb_variables() ; ++i ),
         first->has_variable( variable( i ) ) ||
         second->has_variable( variable( i ) ) ) ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<nb_variables() ; ++i ),
         IMPLIES(
            second->has_variable( variable( i ) ),
            value( variable( i ) ) == second->value( variable( i ) ) ) ) ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<nb_variables() ; ++i ),
         IMPLIES(
            !second->has_variable( variable( i ) ),
            first->has_variable( variable( i ) ) &&
            value( variable( i ) ) == first->value( variable( i ) ) ) ) ) ;
}

//----------------------------------------------------------------------
MAC_ContextPair*
MAC_ContextPair:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextPair:: create_clone" ) ;

   MAC_ContextPair* result = new MAC_ContextPair( a_owner, CT1, CT2  ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->nb_variables() == nb_variables() ) ;
   MAC_CHECK_POST(
      FORALL(
         ( size_t i=0 ; i<nb_variables() ; ++i ),
         result->has_variable( variable( i ) ) &&
         result->value( variable( i ) ) == value( variable( i ) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ContextPair:: ~MAC_ContextPair( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextPair:: ~MAC_ContextPair" ) ;
   
   MAC_CHECK_INV( invariant() ) ;
   if( CT1!=0 )
   {
      CT1->detach_observer( this ) ;
   }
   if( CT2!=0 )
   {
      CT2->detach_observer( this ) ;
   }
//   MAC_Context::update_dependencies(this) ;
   
   notify_observers_of_my_destruction() ;
}

//----------------------------------------------------------------------
size_t
MAC_ContextPair:: nb_variables( void ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextPair:: nb_variables" ) ;
   
   size_t result = VALUES->count() ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Variable const* 
MAC_ContextPair:: variable( size_t i ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextPair:: variable" ) ;
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
MAC_ContextPair:: has_variable( MAC_Variable const* var ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextPair:: has_variable" ) ;
   MAC_CHECK_PRE( has_variable_PRE( var ) ) ;
   
   size_t idx = var->id_number() ;
   bool result = ( VALUES->index_limit()>idx ) &&
                 ( VALUES->at( idx )!=0 ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data* 
MAC_ContextPair:: value( MAC_Variable const* var ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextPair:: value" ) ;
   MAC_CHECK_PRE( value_PRE( var ) ) ;
   
   MAC_Data* result = 
                   static_cast<MAC_Data*>( VALUES->at( var->id_number() ) ) ;

   MAC_CHECK_POST( value_POST( result, var ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_ContextPair:: update( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ContextPair:: update" ) ;
   
   VALUES->resize(MAC_Variable::nb_objects()) ;
   if( CT1!=0 )
   {
      for( size_t i=0 ; i<CT1->nb_variables() ; ++i )
      {
	 MAC_Variable const* var = CT1->variable( i ) ;
	 size_t idx = var->id_number() ;
	 VALUES->set_at( idx, CT1->value(var) ) ;
      }
   }
   if( CT2!=0 )
   {
      for( size_t i=0 ; i<CT2->nb_variables() ; ++i )
      {
	 MAC_Variable const* var = CT2->variable( i ) ;
	 size_t idx = var->id_number() ;
	 VALUES->set_at( idx, CT2->value( var ) ) ;
      }
   }
   update_observers() ;
}

//----------------------------------------------------------------------
void
MAC_ContextPair:: update_for_destruction_of( MAC_Context const* subject )
//----------------------------------------------------------------------
{
   if( subject == CT1 )
   {
      CT1 = 0 ;
   }
   else
   {
      MAC_ASSERT( subject == CT2 ) ;
      CT2 = 0 ;
   }
   VALUES->clear() ;
}

//----------------------------------------------------------------------
bool
MAC_ContextPair:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( VALUES != 0 ) ;
   return( true ) ;
}




