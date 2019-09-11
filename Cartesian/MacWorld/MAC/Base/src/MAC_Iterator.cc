#include <MAC_Iterator.hh>

#include <MAC_assertions.hh>
#include <MAC_Container.hh>


//----------------------------------------------------------------------
MAC_Iterator:: MAC_Iterator( MAC_Object* a_owner,
                             MAC_Container const* list )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , my_list( list )
{
   MAC_CHECK( list!=0 ) ;
   
   record_container_state_id() ;
   
   MAC_CHECK_POST( !container_has_been_modified() ) ;
}




//----------------------------------------------------------------------
void
MAC_Iterator::base_initialize( MAC_Container const* list )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Iterator::base_initialize" ) ;
   MAC_CHECK( list != 0 ) ;
   
   my_list = list ;
   record_container_state_id() ;
   
   MAC_CHECK_POST( !container_has_been_modified() ) ;
}




//----------------------------------------------------------------------
MAC_Iterator:: MAC_Iterator( MAC_Object* a_owner )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , my_list( 0 )
{
   MAC_CHECK_POST( container_has_been_modified() ) ;
}




//----------------------------------------------------------------------
MAC_Iterator:: ~MAC_Iterator( void )
//----------------------------------------------------------------------
{
}




//----------------------------------------------------------------------
bool
MAC_Iterator:: container_has_been_modified( void ) const
//----------------------------------------------------------------------
{
   return( my_list==0 || my_list->state_id()!=list_state ) ;
}




//----------------------------------------------------------------------
void
MAC_Iterator:: record_container_state_id( void )
//----------------------------------------------------------------------
{
   list_state = my_list->state_id() ;
}




//----------------------------------------------------------------------
bool
MAC_Iterator:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool 
MAC_Iterator:: item_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_valid() ) ;  
   return( true ) ;
}




//----------------------------------------------------------------------
bool 
MAC_Iterator:: is_valid_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( !container_has_been_modified() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool 
MAC_Iterator:: go_next_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_valid() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool 
MAC_Iterator:: start_POST( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( !container_has_been_modified() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Iterator:: item_POST( MAC_Object const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   return( true ) ;
}

