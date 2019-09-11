#include <MAC_KeywordDataPair.hh>

#include <doubleVector.hh>
#include <intVector.hh>

#include <MAC_Error.hh>
#include <MAC_KeyItemPair.hh>
#include <MAC_Root.hh>
#include <MAC_KeywordDataPair.hh>
#include <MAC_Data.hh>
#include <MAC_List.hh>
#include <MAC_Vector.hh>
#include <MAC_String.hh>

#include <iostream>


//----------------------------------------------------------------------
MAC_KeywordDataPair*
MAC_KeywordDataPair:: create( MAC_Object* a_owner,
                              MAC_String const* a_keyword, 
                              MAC_Data const* a_data )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeywordDataPair:: create" ) ;
   return( new MAC_KeywordDataPair( a_owner, a_keyword, a_data ) ) ;
} 




//----------------------------------------------------------------------
MAC_KeywordDataPair:: MAC_KeywordDataPair( MAC_Object* a_owner,
                                           MAC_String const* a_keyword, 
                                           MAC_Data const* a_data )
//----------------------------------------------------------------------
   : MAC_Object( a_owner ),
     pair( 0 )
{
   pair = MAC_KeyItemPair::create( this,
                                   const_cast<MAC_String*>(a_keyword), 
                                   const_cast<MAC_Data*>(a_data) ) ;

   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
std::string const&
MAC_KeywordDataPair:: keyword( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeywordDataPair:: keyword" ) ;
   return( static_cast<MAC_String*>( pair->key() )->to_string() ) ;
}




//----------------------------------------------------------------------
MAC_Data const* 
MAC_KeywordDataPair:: data( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeywordDataPair:: data" ) ;
   return( static_cast<MAC_Data const*>( pair->item() ) ) ;
}




//----------------------------------------------------------------------
void
MAC_KeywordDataPair:: replace_data( MAC_Data const* a_data )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeywordDataPair:: replace_data" ) ;
   MAC_CHECK_INV( invariant() ) ;

   pair->set_item( const_cast<MAC_Data*>(a_data) ) ;
}




//----------------------------------------------------------------------
MAC_KeywordDataPair:: ~MAC_KeywordDataPair( void ) 
//----------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}




//-------------------------------------------------------------------------
size_t
MAC_KeywordDataPair:: hash_code( void ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeywordDataPair:: hash_code" ) ;
   return( pair->hash_code() ) ;  
}




//----------------------------------------------------------------------
bool
MAC_KeywordDataPair:: comparable( MAC_Object const* other ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeywordDataPair:: comparable" ) ;
   return MAC_Object::comparable( other ) ||
      pair->comparable( other ) ;
}




//----------------------------------------------------------------------
bool
MAC_KeywordDataPair:: is_equal( MAC_Object const* other ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeywordDataPair:: is_equal" ) ;
   MAC_CHECK_PRE( is_equal_PRE( other ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   MAC_KeywordDataPair const* other_as_assignment =
                           dynamic_cast<MAC_KeywordDataPair const*>( other ) ;
   bool result ;
   if( other_as_assignment != 0 )
   {
      result = pair->is_equal( other_as_assignment->pair ) ;
   }
   else
   {
      result = pair->is_equal( other ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( is_equal_POST( result, other ) ) ;

   return( result ) ;
}




//-------------------------------------------------------------------------
int
MAC_KeywordDataPair:: three_way_comparison( MAC_Object const* other ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeywordDataPair:: three_way_comparison" ) ;
   MAC_Error::object()->raise_not_implemented( this, "three_way_comparison" ) ;
   return( 1 ) ;
}




//----------------------------------------------------------------------
void
MAC_KeywordDataPair:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_KeywordDataPair:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << keyword() ;
   os << " = " ;
   bool is_str = data()->data_type()==MAC_Data::String ;
   
   if( is_str )
   {
      os << "\"" ;
   }
   data()->print( os, 0 ) ;
   if( is_str )
   {
      os << "\"" ;
   }
}
