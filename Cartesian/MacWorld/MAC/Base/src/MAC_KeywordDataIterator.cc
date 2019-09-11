#include <MAC_KeywordDataIterator.hh>

#include <MAC_KeywordDataPair.hh>
#include <MAC_List.hh>
#include <MAC_ListItem.hh>
#include <MAC_assertions.hh>



//----------------------------------------------------------------------
MAC_KeywordDataIterator*
MAC_KeywordDataIterator:: create( MAC_Object* a_owner,
                                  MAC_List const* a_list )
//----------------------------------------------------------------------
{
   return ( new MAC_KeywordDataIterator( a_owner, a_list ) ) ;
}



//----------------------------------------------------------------------
MAC_KeywordDataIterator:: MAC_KeywordDataIterator( MAC_Object* a_owner,
                                                   MAC_List const* a_list )
//----------------------------------------------------------------------
   : MAC_ListIterator( a_owner, a_list )
{
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
MAC_KeywordDataIterator:: ~MAC_KeywordDataIterator( void )
//----------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
MAC_KeywordDataPair*
MAC_KeywordDataIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_PRE( item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ; 

   MAC_Object* obj = MAC_ListIterator::item() ;
   MAC_CHECK( dynamic_cast<MAC_KeywordDataPair *>( obj )!=0 ) ;
   
   MAC_KeywordDataPair* resu = static_cast<MAC_KeywordDataPair *>( obj ) ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( item_POST( resu ) ) ;

   return( resu ) ;
}



//---------------------------------------------------------------------
bool
MAC_KeywordDataIterator:: invariant( void ) const
//---------------------------------------------------------------------
{
   MAC_ASSERT( MAC_ListIterator::invariant() ) ;

   return( true ) ;
}
