#include <MAC_ModuleIterator.hh>

#include <MAC_assertions.hh>
#include <MAC_List.hh>
#include <MAC_ListItem.hh>
#include <MAC_Module.hh>



//----------------------------------------------------------------------
MAC_ModuleIterator*
MAC_ModuleIterator:: create( MAC_Object* a_owner, MAC_List const* a_list )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleIterator:: create" ) ;
   MAC_CHECK_PRE( a_list != 0 ) ;

   MAC_ModuleIterator* result = new MAC_ModuleIterator( a_owner, a_list ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_ModuleIterator:: MAC_ModuleIterator( MAC_Object* a_owner,
                                         MAC_List const* a_list )
//----------------------------------------------------------------------
   : MAC_ListIterator( a_owner, a_list )
{
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
MAC_ModuleIterator:: ~MAC_ModuleIterator( void )
//----------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
MAC_Module*
MAC_ModuleIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ; 

   MAC_Object* obj = MAC_ListIterator::item() ;
   
   MAC_Module* resu = dynamic_cast<MAC_Module*>( obj ) ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( item_POST( resu ) ) ;

   return( resu ) ;
}




//---------------------------------------------------------------------
bool
MAC_ModuleIterator:: invariant( void ) const
//---------------------------------------------------------------------
{
   MAC_ASSERT( MAC_ListIterator::invariant() ) ;

   return( true ) ;
}
