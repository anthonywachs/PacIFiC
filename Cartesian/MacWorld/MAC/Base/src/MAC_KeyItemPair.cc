#include <MAC_KeyItemPair.hh>

#include <MAC_assertions.hh>


//-------------------------------------------------------------------------
MAC_KeyItemPair*
MAC_KeyItemPair:: create( MAC_Object* a_owner,
                          MAC_Object* a_key,
                          MAC_Object* a_item  ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeyItemPair:: create" ) ;
   return( new MAC_KeyItemPair( a_owner, a_key, a_item ) ) ;
}




//-------------------------------------------------------------------------
MAC_KeyItemPair:: MAC_KeyItemPair( MAC_Object* a_owner,
                                   MAC_Object* a_key,
                                   MAC_Object* a_item ) 
//-------------------------------------------------------------------------
   : MAC_Object( a_owner ), 
     the_key( a_key ),
     the_item( a_item )
{
   MAC_LABEL( "MAC_KeyItemPair:: MAC_KeyItemPair" ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//-------------------------------------------------------------------------
MAC_KeyItemPair:: ~MAC_KeyItemPair( void ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeyItemPair:: ~MAC_KeyItemPair" ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//-------------------------------------------------------------------------
size_t
MAC_KeyItemPair:: hash_code( void ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeyItemPair:: hash_code" ) ;
   return( key()->hash_code() ) ;  
}




//-------------------------------------------------------------------------
bool
MAC_KeyItemPair:: invariant( void ) const
//-------------------------------------------------------------------------
{
   MAC_ASSERT( key() != 0 ) ;

   return( true ) ;
} 




//-------------------------------------------------------------------------
MAC_Object*
MAC_KeyItemPair:: item( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeyItemPair:: item" ) ;
   return( the_item ) ;
} 




//-------------------------------------------------------------------------
void
MAC_KeyItemPair:: set_item( MAC_Object* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeyItemPair:: set_item" ) ;
   the_item = object ;
} 




//-------------------------------------------------------------------------
MAC_Object*
MAC_KeyItemPair:: key( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeyItemPair:: key" ) ;
   return( the_key ) ;
}




//-------------------------------------------------------------------------
bool
MAC_KeyItemPair:: is_equal( MAC_Object const* other ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeyItemPair:: is_equal" ) ;
   MAC_CHECK_PRE( is_equal_PRE( other ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   MAC_KeyItemPair const* other_as_pair =
                           dynamic_cast<MAC_KeyItemPair const*>( other ) ;
   bool resu ;
   if( other_as_pair != 0 )
   {
      resu = key()->is_equal( other_as_pair->key() ) ;
   }
   else
   {
      resu = key()->is_equal( other ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( is_equal_POST( resu, other ) ) ;

   return( resu ) ;
}




//-------------------------------------------------------------------------
bool
MAC_KeyItemPair:: comparable( MAC_Object const* other ) const 
//-------------------------------------------------------------------------
{
   return MAC_Object::comparable( other ) ||
      key()->comparable( other ) ;
}




//-------------------------------------------------------------------------
int
MAC_KeyItemPair:: three_way_comparison( MAC_Object const* other ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeyItemPair:: three_way_comparison" ) ;
   MAC_CHECK_PRE( three_way_comparison_PRE( other ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_KeyItemPair const* other_as_pair =
                             dynamic_cast<MAC_KeyItemPair const*>( other ) ;
   int resu ;
   if( other_as_pair != 0 )
   {
      resu = key()->three_way_comparison( other_as_pair->key() ) ;
   }
   else
   {
      resu = key()->three_way_comparison( other ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( three_way_comparison_POST( resu, other ) ) ;

   return( resu );   
}
