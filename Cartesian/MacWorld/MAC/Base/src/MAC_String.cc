#include <MAC_String.hh>

#include <iostream>
#include <numeric>

#include <MAC_assertions.hh>

//----------------------------------------------------------------------------
MAC_String*
MAC_String:: create( MAC_Object* a_owner, std::string const& a )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_String:: create" ) ;
   return( new MAC_String( a_owner, a ) ) ;
}



//----------------------------------------------------------------------------
MAC_String:: MAC_String( MAC_Object* a_owner, std::string const& a )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner ),
     str( a )
{
   MAC_LABEL( "MAC_String:: MAC_String" ) ;
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
MAC_String:: ~MAC_String( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_String:: ~MAC_String" ) ;
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
MAC_Data::Type
MAC_String:: data_type( void) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_String:: data_type" ) ;
   return( MAC_Data::String ) ;
}



//----------------------------------------------------------------------------
MAC_String*
MAC_String:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_String:: create_clone" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_String* resu = new MAC_String( a_owner, str ) ;
   MAC_CHECK_POST( resu!=0 ) ;
   return resu ;
}



//----------------------------------------------------------------------------
std::string const&
MAC_String:: to_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_String:: to_string" ) ;
   MAC_CHECK_PRE( to_string_PRE(ct) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( str ) ;
}



//----------------------------------------------------------------------------
bool
MAC_String:: is_equal( MAC_Object const* other ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_String:: is_equal" ) ;
   MAC_CHECK_PRE( is_equal_PRE( other ) ) ;
   
   bool result = three_way_comparison( other ) == 0 ;
   MAC_CHECK_POST( is_equal_POST( result, other ) ) ;
   
   return( result ) ;
}



//----------------------------------------------------------------------------
int
MAC_String:: three_way_comparison( MAC_Object const* other ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_String:: three_way_comparison" ) ;
   MAC_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   MAC_String const* tString = static_cast<MAC_String const*>( other ) ;

   int result = str.compare( tString->str ) ;
   MAC_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return result ;
}



//-------------------------------------------------------------------------
size_t
MAC_String:: hash_code( void ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_String:: hash_code" ) ;
   size_t resu = std::accumulate( str.begin(), str.end(), 0 ) ;
   return( resu ) ;
}



//----------------------------------------------------------------------
void
MAC_String:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_String:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   if(str.find('"')<str.length())
   {
      std::string dup(str) ;
      size_t idx ;
      while( (idx=dup.find('"'))<dup.length()) dup.replace( idx, 1, "'" ) ;
      os << space << "\"" << dup << "\"" ;
   }
   else
   {
      os << space << "\"" << str << "\"" ;
   }
   
}



//----------------------------------------------------------------------------
void
MAC_String:: set( std::string const& val )
//----------------------------------------------------------------------------
{
   str = val ;
}



