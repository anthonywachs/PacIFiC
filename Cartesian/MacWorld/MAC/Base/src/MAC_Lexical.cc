#include <MAC_Lexical.hh>

#include <MAC_assertions.hh>
#include <MAC_Module.hh>
#include <MAC_Error.hh>
#include <MAC_List.hh>
#include <MAC_String.hh>

#include <iostream>

//----------------------------------------------------------------------
MAC_Lexical* MAC_Lexical::the_common_owner = 0 ;
//----------------------------------------------------------------------



//----------------------------------------------------------------------
MAC_Lexical*
MAC_Lexical:: common_owner( void )
//----------------------------------------------------------------------
{
   if( the_common_owner==0 )
   {
      the_common_owner = new MAC_Lexical() ;
   }
   return the_common_owner ;
}



//----------------------------------------------------------------------
void
MAC_Lexical:: remove_all_lexical( void )
//----------------------------------------------------------------------
{
   if( the_common_owner!=0 )
   {
      the_common_owner->destroy() ;
      the_common_owner=0 ;
   }
}



//----------------------------------------------------------------------
MAC_Lexical*
MAC_Lexical:: create( MAC_Module* aModule )
//----------------------------------------------------------------------
{
   return new MAC_Lexical( common_owner(), aModule ) ;
}



//----------------------------------------------------------------------
MAC_Lexical*
MAC_Lexical:: create( MAC_Data* aSimpleValue )
//----------------------------------------------------------------------
{
   return new MAC_Lexical( common_owner(), aSimpleValue ) ;
}



//----------------------------------------------------------------------
MAC_Lexical*
MAC_Lexical:: create( MAC_List* aSimpleValue )
//----------------------------------------------------------------------
{
   return new MAC_Lexical( common_owner(), aSimpleValue ) ;
}



//----------------------------------------------------------------------
MAC_Lexical:: MAC_Lexical( void )
//----------------------------------------------------------------------
      :   MAC_Object( 0 ),
          myModule( 0 ),
          myData( 0 ),
          myList( 0 )
{
}



//----------------------------------------------------------------------
MAC_Lexical:: MAC_Lexical( MAC_Object* anOwner,
                           MAC_Module* aModule )
//----------------------------------------------------------------------
      :   MAC_Object( anOwner ),
          myModule( aModule ),
          myData( 0 ),
          myList( 0 )
{
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
MAC_Lexical:: MAC_Lexical( MAC_Object* anOwner,
                           MAC_Data* aSimpleValue )
//----------------------------------------------------------------------
      :   MAC_Object( anOwner ),
          myModule( 0 ),
          myData( aSimpleValue ),
          myList( 0 )
{
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK( aSimpleValue->owner()==0 ) ;
   aSimpleValue->set_owner( this ) ;
}




//----------------------------------------------------------------------
MAC_Lexical:: MAC_Lexical( MAC_Object* anOwner,
                           MAC_List* aSimpleValue )
//----------------------------------------------------------------------
      :   MAC_Object( anOwner ),
          myModule( 0 ),
          myData( 0 ),
          myList( aSimpleValue )
{
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK( aSimpleValue->owner()==0 ) ;
   aSimpleValue->set_owner( this ) ;
}




//----------------------------------------------------------------------
MAC_Lexical:: ~MAC_Lexical( void ) 
//----------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------
bool
MAC_Lexical:: is_list( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Lexical:: is_list" ) ;
   MAC_CHECK_INV( invariant() ) ;
   return myList!=0 ;
}



//----------------------------------------------------------------------
bool
MAC_Lexical:: is_module( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Lexical:: is_module" ) ;
   MAC_CHECK_INV( invariant() ) ;
   return myModule!=0 ;
}



//----------------------------------------------------------------------
bool
MAC_Lexical:: is_data( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Lexical:: is_data" ) ;
   MAC_CHECK_INV( invariant() ) ;
   return myData!=0 ;
}



//-------------------------------------------------------------------------
MAC_List *
MAC_Lexical:: to_list( void ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Lexical:: to_module" ) ;
   MAC_CHECK_PRE( is_list() ) ;
   MAC_CHECK_INV( invariant() ) ;
   return myList ;
}



//-------------------------------------------------------------------------
MAC_Module *
MAC_Lexical:: to_module( void ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Lexical:: to_module" ) ;
   MAC_CHECK_PRE( is_module() ) ;
   MAC_CHECK_INV( invariant() ) ;
   return myModule ;
}



//-------------------------------------------------------------------------
MAC_Data *
MAC_Lexical:: to_data( void )  
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Lexical:: to_data" ) ;
   MAC_CHECK_PRE( is_data() ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_Data *resu = myData ;
   return resu ;
}



//----------------------------------------------------------------------
bool
MAC_Lexical:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   MAC_ASSERT( this==the_common_owner || myModule!=0 || myData!=0 || myList!=0 ) ;
   MAC_ASSERT( myModule==0 || myModule->owner()!=0 ) ;
   return( true ) ;
}





