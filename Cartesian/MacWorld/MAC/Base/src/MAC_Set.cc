#include <MAC_Set.hh>

#include <MAC_assertions.hh>


//-----------------------------------------------------------------------------
MAC_Set:: MAC_Set( MAC_Object* a_owner )
//-----------------------------------------------------------------------------
   : MAC_Collection( a_owner )
{
}



//-----------------------------------------------------------------------------
MAC_Set:: ~MAC_Set( void ) 
//-----------------------------------------------------------------------------
{
}



//-----------------------------------------------------------------------------
bool
MAC_Set:: extend_POST( bool old_has, 
                       size_t old_count, 
                       MAC_Object const* object,
                       size_t old_state_id  ) const
//-----------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Collection::extend_POST( old_has, old_count, object, old_state_id ) ) ;

   MAC_ASSERT( !old_has || ( count() == old_count )   ) ;
   MAC_ASSERT(  old_has || ( count() == old_count+1 ) ) ;
   MAC_ASSERT( state_id() != old_state_id ) ;
   return( true ) ;
}



//-----------------------------------------------------------------------------
bool
MAC_Set:: remove_POST( size_t old_count,
                       MAC_Object const* object,
                       size_t old_state_id ) const
//-----------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Collection::remove_POST( old_count, object, old_state_id ) ) ;

   MAC_ASSERT( !has( object ) ) ;

   return( true ) ; 
}



