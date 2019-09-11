#include <MAC_Collection.hh>

#include <iostream>

#include <MAC_Iterator.hh>
#include <MAC_assertions.hh>

using std::cerr ;
using std::endl ;
using std::ostream ;


//----------------------------------------------------------------------
MAC_Collection:: MAC_Collection( MAC_Object* a_owner )
//----------------------------------------------------------------------
    : MAC_Container( a_owner )
{
   MAC_LABEL( "MAC_Collection:: MAC_Collection" ) ;
}



//----------------------------------------------------------------------
MAC_Collection:: ~MAC_Collection( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Collection:: ~MAC_Collection" ) ;
}



//----------------------------------------------------------------------
bool
MAC_Collection:: extend_PRE( MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( object != 0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Collection:: extend_POST( bool old_has,
                              size_t old_count,
                              MAC_Object const* object,
                              size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( has( object ) ) ;
   MAC_ASSERT( old_state_id!=state_id() ) ;
   
   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Collection:: remove_PRE( MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( has( object ) ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Collection:: remove_POST( size_t old_count,
                              MAC_Object const* object,
                              size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( count() == old_count-1 ) ;
   MAC_ASSERT( old_state_id!=state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Collection:: clear_POST( size_t old_state_id ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( count() == 0 ) ;
   MAC_ASSERT( old_state_id!=state_id() ) ;

   return( true ) ;
}



