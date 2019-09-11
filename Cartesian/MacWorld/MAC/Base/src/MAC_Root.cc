#include <MAC_Root.hh>

#include <MAC_assertions.hh>
#include <MAC_Vector.hh>

#include <iostream>

//-------------------------------------------------------------------------
MAC_Object*
MAC_Root:: object( size_t i ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Root:: object" ) ;

   MAC_Vector* objs = objects() ;
   if( objs->index_limit()<=i )
   {
      objs->resize( i+1 ) ;
   }
   if( objs->at(i) == 0 )
   {
      objs->set_at( i, new MAC_Root( 0 ) ) ;
   }
   
   MAC_Object* result = objs->at(i) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
MAC_Root:: MAC_Root( MAC_Object* a_owner ) 
//-------------------------------------------------------------------------
   : MAC_Object( a_owner )
{
}

//-------------------------------------------------------------------------
MAC_Root:: ~MAC_Root( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
MAC_Root:: cleanup( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Root:: cleanup" ) ;
   
   MAC_Vector* objs = objects() ;
   objs->destroy_items_and_remove_section( 0, objs->index_limit() ) ;
   objects( true ) ;
}

//-------------------------------------------------------------------------
MAC_Vector*
MAC_Root:: objects( bool destroy )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Root:: objects" ) ;

   static MAC_Vector* result = 0 ;
   if( result == 0 )
   {
      result = MAC_Vector::create( 0, 0 ) ;
   }
   if( destroy )
   {
      result->destroy() ;
      result = 0 ;
   }
   
   MAC_CHECK_POST( IMPLIES( destroy, result == 0 ) ) ;
   MAC_CHECK_POST( IMPLIES( !destroy, result != 0 && result->owner() == 0 ) ) ;
   return( result ) ;
}


