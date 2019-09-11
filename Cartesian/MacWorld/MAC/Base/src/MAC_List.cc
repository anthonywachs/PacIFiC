#include <MAC_List.hh>

#include <iostream>

#include <MAC_assertions.hh>
#include <MAC_ListItem.hh>
#include <MAC_ListIterator.hh>

using std::ostream ;
using std::endl ;


//----------------------------------------------------------------------
MAC_List*
MAC_List:: create( MAC_Object* a_owner )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: create" ) ;
   MAC_List* result = new MAC_List( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->index_limit() == 0 ) ;

   return( result ) ;
}



//-------------------------------------------------------------------------
MAC_List:: MAC_List( MAC_Object* a_owner )
//-------------------------------------------------------------------------
   : MAC_Sequence( a_owner ),
     theList(0), 
     theLast(0), 
     nbElem( 0 )
{
   MAC_LABEL( "MAC_List:: MAC_List" ) ;
   MAC_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
MAC_List:: ~MAC_List( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: ~MAC_List" ) ;
   MAC_CHECK_INV( invariant() ) ;

   clear() ;
}



//-------------------------------------------------------------------------
MAC_List*
MAC_List:: create_clone( MAC_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: create_clone" ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_List* result = create( a_owner ) ;
   MAC_ListItem* ret = theList ;         
   while( ret!=0 )
   {
      result->append( ret->val() ) ;
      ret = ret->next() ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;  

   return( result ) ;
}



//-------------------------------------------------------------------------
void      
MAC_List:: copy( MAC_List const* other )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: copy" ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   MAC_CHECK_PRE( other !=0 ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( this != other )
   {
      clear() ;
      MAC_ListItem* ret = other->theList ;
         
      while( ret!=0 )
      {
         append( ret->val() ) ;
         ret = ret->next() ;
      }
   }
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( state_id() != OLD(state_id) ) ;
   MAC_CHECK_POST( index_limit() == other->index_limit() ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<index_limit() ; ++i ),
                           at(i) == other->at(i) ) ) ;
}



//-------------------------------------------------------------------------
size_t
MAC_List:: index_limit( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: index_limit" ) ;
   MAC_CHECK_INV( invariant() ) ;

   return( nbElem ) ;
}



//-------------------------------------------------------------------------
void
MAC_List:: append( MAC_Object* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: append" ) ;
   MAC_CHECK_PRE( append_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   MAC_ListItem* theNew = new MAC_ListItem( object ) ;
   
   if( nbElem==0 )
   {
      theList = theNew ;
   }
   else
   {
      theLast->doLink( theNew ) ;
   }
   theLast = theNew ;
   nbElem++ ;
   report_state_change() ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( append_POST( OLD( index_limit ),
                                OLD( count ),
                                object,
                                OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
void
MAC_List:: prepend( MAC_Object* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: prepend" ) ;
   MAC_CHECK_PRE( prepend_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   MAC_ListItem* theNew = new MAC_ListItem( object ) ;
   MAC_ListItem* old = theList ;

   theList = theNew ;

   if( old==0 )
   {
      theLast = theNew ;
   }
   else
   {
      theNew->doLink( old ) ;
   }
   
   nbElem++ ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( prepend_POST( OLD( index_limit ),
                                 OLD( count ),
                                 object,
                                 OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
void
MAC_List:: set_at( size_t i, MAC_Object* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: set_at" ) ;
   MAC_CHECK_PRE( set_at_PRE( i, object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   MAC_ListItem* it = theItem( i ) ;
   it->replaceValue(object) ;
   report_state_change() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_at_POST( OLD( index_limit ),
                                OLD( count ),
                                i,
                                object,
                                OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
void
MAC_List:: insert_at( size_t i, MAC_Object* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: insert_at" ) ;
   MAC_CHECK_PRE( insert_at_PRE( i, object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   MAC_ListItem* theNew = new MAC_ListItem( object ) ;
   MAC_ListItem* old ;
   
   if( i==0 )
   {
      old = theList ;
      theList = theNew ;
   }
   else
   {
      MAC_ListItem* ptr = theItem( i-1 ) ;
      old = ptr->next() ;
      ptr->doLink( theNew ) ;
   }
   theNew->doLink( old ) ;
   if( old==0 )
   {
      theLast=theNew ;
   }
   nbElem++ ;
   report_state_change() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( insert_at_POST( OLD( index_limit ),
                                   OLD( count ),
                                   i,
                                   object,
                                   OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
size_t
MAC_List:: count( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: count" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( count_POST( nbElem ) ) ;
   return( nbElem ) ;
}



//-------------------------------------------------------------------------
MAC_Object*
MAC_List:: item( MAC_Object const* object ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: item" ) ;
   MAC_CHECK_PRE( item_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_ListItem* ptr = theList ;
   MAC_Object* result = 0 ;
   while( ptr!=0 )
   {
      if( matching_items( ptr->val(), object ) )
      {
         result = ptr->val() ;
         break ;
      }
      ptr = ptr->next() ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( item_POST( result, object ) ) ;

   return( result ) ; 
}



//-------------------------------------------------------------------------
MAC_Object*
MAC_List:: at( size_t i ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: at" ) ;
   MAC_CHECK_PRE( at_PRE( i ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Object* result = theItem(i)->val() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( at_POST( result , i ) ) ;

   return( result ) ;
}



//-------------------------------------------------------------------------
size_t
MAC_List:: index_of( MAC_Object const* object ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: index_of" ) ;
   MAC_CHECK_PRE( index_of_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_ListItem* ptr = theList ;
   size_t ret = badIndex ;
   size_t cpt = 0 ;
   
   while( ptr!=0 )
   {
      if( matching_items( ptr->val(), object ) )
      {
         ret = cpt ;
         break ;
      }
      ptr = ptr->next() ;
      cpt++ ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( index_of_POST( ret, object ) ) ;

   return( ret ) ;
}



//-------------------------------------------------------------------------
MAC_ListIterator*
MAC_List:: create_iterator( MAC_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: create_iterator" ) ;
   MAC_ListIterator* result = MAC_ListIterator::create( a_owner, this ) ;

   MAC_CHECK_POST( create_iterator_POST( result, a_owner ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
void
MAC_List:: remove( MAC_Object const* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: remove" ) ;
   MAC_CHECK_PRE( remove_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   MAC_ListItem* ptr = theList ;
   MAC_ListItem* old = 0 ;
   MAC_Object * ret = 0 ;
   
   while( ptr!=0 )
   {
      if( matching_items( ptr->val(), object ) )
      {
         MAC_ListItem* nextItem = ptr->next() ;
         if( old==0 )
         {
            theList = nextItem ;
            if( theLast==ptr )
            {
               theLast=0 ;
            }
         }
         else
         {
            old->doLink( nextItem ) ;
            if( theLast==ptr )
            {
               theLast=old ;
            }
         }
         ret = ptr->val() ;
         delete ptr ;
         break ;
      }
      old = ptr ;
      ptr = ptr->next() ;
   }
   nbElem-- ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( remove_POST( OLD( count ), object, OLD( state_id ) ) ) ;

//   return( ret ) ;
}



//-------------------------------------------------------------------------
void 
MAC_List:: remove_at( size_t i )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: remove_at" ) ;
   MAC_CHECK_PRE( remove_at_PRE( i ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   MAC_ListItem* old ;
   
   if( i==0 )
   {
      old=theList ;
      theList = old->next() ;
      if( old == theLast )
      {
         theLast = 0 ;
      }
   }
   else
   {
      MAC_ListItem* ptr = theItem( i-1 ) ;
      old = ptr->next() ;
      ptr->doLink( old->next() ) ;
      if( old == theLast )
      {
         theLast = ptr ;
      }
   }
   
   delete old ;
   nbElem-- ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( remove_at_POST( OLD( index_limit ),
                                   OLD( count ),
                                   OLD( state_id ) ) ) ;
}



//---------------------------------------------------------------------
void
MAC_List:: remove_section( size_t iFirst, size_t length )
//---------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: remove_section" ) ;
   MAC_CHECK_PRE( remove_section_PRE( iFirst, length ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
      
   MAC_ListItem* old = theItem( iFirst+length-1 ) ;
   MAC_ListItem* end = old->next() ;
   MAC_ListItem* firstToDelete ;
   
   if( iFirst==0 )
   {
      firstToDelete = theList ;
      theList = end ;
      if( old == theLast )
      {
         theLast = 0 ;
      }
   }
   else
   {
      MAC_ListItem* prev = theItem( iFirst-1 ) ;
      firstToDelete = prev->next() ;
      prev->doLink( end ) ;
      if( old == theLast )
      {
         theLast = prev ;
      }

   }
   for( size_t cpt = 0 ; cpt<length ; cpt++ )
   {
      MAC_ListItem* ptr = firstToDelete->next() ;
      delete firstToDelete ;
      firstToDelete = ptr ;
   }
   
   MAC_CHECK( end == firstToDelete ) ;
   
   nbElem -= length ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( remove_section_POST( OLD( index_limit ), 
                                        OLD( count ), 
                                        length,
                                        OLD( state_id ) ) ) ;
}



//---------------------------------------------------------------------
void
MAC_List:: destroy_items_and_remove_section( size_t iFirst, size_t length )
//---------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: destroy_items_and_remove_section" ) ;
   MAC_CHECK_PRE( destroy_items_and_remove_section_PRE( iFirst, length ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   MAC_ListItem* deb = theItem( iFirst ) ;
   for( size_t cpt = 0 ; cpt<length ; cpt++ )
   {
      MAC_ListItem* nextItem = deb->next() ;
      deb->val()->destroy() ;
      deb=nextItem ;
   }
   remove_section( iFirst, length ) ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( remove_section_POST( OLD( index_limit ), 
                                        OLD( count ), 
                                        length,
                                        OLD( state_id ) ) ) ;
   MAC_CHECK_POST( index_limit() == OLD( index_limit ) - length ) ;
   MAC_CHECK_POST( count() == OLD( count ) - length ) ;
}



//----------------------------------------------------------------------
void
MAC_List:: destroy_items_and_clear( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: destroy_items_and_clear" ) ;
   MAC_CHECK_PRE( destroy_items_and_remove_section_PRE( 0, index_limit() ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   MAC_ListItem* ptr = theList ;
   while( ptr!=0 )
   {
      MAC_ListItem* nextItem = ptr->next() ;
      ptr->val()->destroy() ;
      ptr = nextItem ;
   }
   clear() ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( clear_POST( OLD( state_id ) ) ) ;

}



//-------------------------------------------------------------------------
void
MAC_List:: clear( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_List:: clear" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   MAC_ListItem* ptr = theList ;
   MAC_ListItem* old ;
   
   while( ptr!=0 )
   {
      old = ptr ;
      ptr = ptr->next() ;
      delete old ;
   }
   theList = 0 ;
   theLast = 0 ;
   nbElem = 0 ;     
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( clear_POST( OLD( state_id ) ) ) ;
}



//----------------------------------------------------------------------
// ostream&
// operator<<( ostream& out, MAC_List const& l )
//----------------------------------------------------------------------
// {
//    out << "List : " << endl ;
//    out << "  Nb elms = " << l.si() << endl ;
//    return( out ) ;
// }



//-------------------------------------------------------------------------
MAC_ListItem*
MAC_List:: theItem( size_t n ) const
//-------------------------------------------------------------------------
{
   size_t nb = n ;
   MAC_CHECK( n<nbElem ) ; //??????????????????????????
   MAC_ListItem* ret = theList ;
         
   while( nb>0 )
   {
      nb-- ;
      ret=ret->next() ;
   }
   return( ret ) ;
}



//-------------------------------------------------------------------------
bool
MAC_List:: invariant( void ) const
//-------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Sequence::invariant() ) ;
   MAC_ASSERT( count() == index_limit() ) ;

   return( true ) ;
}



//------------------------------------------------------------------------
bool
MAC_List:: set_at_POST( size_t old_index_limit,
                        size_t old_count,
                        size_t i,
                        MAC_Object const* object,
                        size_t old_state_id ) const
//------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Sequence::set_at_POST( old_index_limit, 
                                          old_count, 
                                          i, 
                                          object,
                                          old_state_id )  ) ;

   MAC_ASSERT( old_count == count() ) ;
   MAC_ASSERT( state_id() != old_state_id ) ;
   
   return( true ) ;
}



//-------------------------------------------------------------------------
bool
MAC_List:: at_POST( MAC_Object const* result, size_t i ) const
//-------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Sequence::at_POST( result, i ) ) ;
   MAC_ASSERT( result != 0 ) ;

   return( true ) ;
}


//----------------------------------------------------------------------
bool
MAC_List:: remove_at_POST( size_t old_index_limit,
                           size_t old_count,
                           size_t old_state_id ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Sequence::remove_at_POST( old_index_limit,
                                             old_count,
                                             old_state_id  ) ) ;

   MAC_ASSERT( index_limit() == old_index_limit - 1 ) ;
   MAC_ASSERT( count() == old_count - 1 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_List:: remove_section_POST( size_t old_index_limit,
                                size_t old_count,
                                size_t length,
                                size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Sequence::remove_section_POST( old_index_limit, 
                                                  old_count, 
                                                  length,
                                                  old_state_id ) ) ;

   MAC_ASSERT( index_limit() == old_index_limit - length ) ;
   MAC_ASSERT( count() == old_count - length ) ;

   return( true ) ;
}


