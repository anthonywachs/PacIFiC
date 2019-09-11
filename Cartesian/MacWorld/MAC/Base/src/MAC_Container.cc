#include <MAC_Container.hh>

#include <iostream>

#include <MAC_Iterator.hh>
#include <MAC_assertions.hh>

using std::endl ;
using std::ostream ;

//-----------------------------------------------------------------------------
MAC_Container:: MAC_Container( MAC_Object* a_owner )
//-----------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , MY_STATE( 0 )
{
   MAC_LABEL( "MAC_Container:: MAC_Container" ) ;
}

//----------------------------------------------------------------------
MAC_Container:: ~MAC_Container( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Container:: ~MAC_Container" ) ;
   report_state_change() ;
}

//-------------------------------------------------------------------------
bool
MAC_Container:: matching_items( MAC_Object const* object1,
                                MAC_Object const* object2 ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Container:: matching_items" ) ;
   MAC_CHECK_PRE( matching_items_PRE( object1, object2 ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   return( object1->is_equal( object2 ) ) ;

   MAC_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------
size_t
MAC_Container:: state_id( void ) const
//-----------------------------------------------------------------------
{
   return MY_STATE ;
}

//-------------------------------------------------------------------------
bool
MAC_Container:: has( MAC_Object const* object ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Container:: has" ) ;
   MAC_CHECK_PRE( object != 0 ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   bool result = ( item( object ) != 0 ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( IMPLIES( result!=0 , count()!=0 ) ) ;
   MAC_CHECK_POST( OLD( state_id )==state_id() ) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
void
MAC_Container:: print( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Container:: print" ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   std::string const s( indent_width, ' ' ) ;
   os << s << "Collection of type " << type_name() ;
   if( count()==1 )
   {
      os << " contains " << count() << " item :" << endl ;
      MAC_Iterator* it = create_iterator( 0 ) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         MAC_Object const* obj = it->item() ;
         obj->display_info( os, indent_width+3 ) ;
         obj->print( os, indent_width+3 ) ;
         os << endl ;
      }
      it->destroy() ;
   }
   else if( count()!=0 )
   {
      os << " contains " << count() << " items :" << endl ;
      MAC_Iterator* it = create_iterator( 0 ) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         MAC_Object const* obj = it->item() ;
         os << s << " - " << endl ;
         obj->display_info( os, indent_width+3 ) ;
         obj->print( os, indent_width+3 ) ;
         os << endl ;
         os << s << "    -----------------------------" << endl ;
      }
      it->destroy() ;
   }
   else
   {
      os << " contains " << count() << " item" << endl ;
   }
   
   MAC_CHECK_POST( OLD( state_id )==state_id() ) ;

}

//----------------------------------------------------------------------
bool
MAC_Container:: create_clone_POST( MAC_Container* result,
                                   MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{ 
   MAC_ASSERT( MAC_Object::create_clone_POST( result, a_owner ) ) ;
   MAC_ASSERT( result->count() == count() ) ;
   MAC_Iterator* it1 = 0 ;
   MAC_Iterator* it2 = 0 ;
   MAC_ASSERT( 
      FORALL( ( it1 = create_iterator( 0 ),
                it2 = result->create_iterator( 0 ) ;
                it1->is_valid() && it2->is_valid() ;
                it1->go_next(), it2->go_next() ),
                it1->item()->has_same_address( it2->item() ) ) ) ;
   
   it1->destroy() ;
   it2->destroy() ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Container:: count_POST( size_t result ) const
//----------------------------------------------------------------------
{
   MAC_Iterator* it = create_iterator(0) ;
   size_t val = 0 ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      if( it->item()!=0 )
      {
         val++ ;
      }
   }
   it->destroy() ;
   MAC_ASSERT( result==val ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Container:: matching_items_PRE( MAC_Object const* object1,
                                    MAC_Object const* object2 ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( (object1 != 0) && (object2 != 0) ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Container:: item_PRE( MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( object != 0 ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Container:: item_POST( MAC_Object const* result,
                           MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result==0 || 
               ( count() != 0 &&  matching_items( result, object ) ) ) ;

  
   return( true ) ;
}

//-----------------------------------------------------------------------
bool
MAC_Container:: create_iterator_POST( MAC_Iterator* result,
                                      MAC_Object* a_owner ) const
//-----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   MAC_ASSERT( EQUIVALENT( count() != 0 , result->is_valid() ) ) ;

   return( true ) ;
}

//-----------------------------------------------------------------------
void
MAC_Container:: report_state_change( void ) 
//-----------------------------------------------------------------------
{
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   MY_STATE ++ ;
   MAC_CHECK_POST( OLD(state_id)!=state_id() ) ;
}

