#include <MAC_IndexSet.hh>

#include <MAC_assertions.hh>
#include <MAC.hh>

#include <iostream>
#include <string>

//-------------------------------------------------------------------------
MAC_IndexSet*
MAC_IndexSet:: create( MAC_Object* a_owner ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IndexSet:: create" ) ;
   
   MAC_IndexSet* result = new MAC_IndexSet( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->elements().size() == 0 ) ;
   MAC_CHECK_POST( result->id() == MAC::bad_index() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
MAC_IndexSet:: MAC_IndexSet( MAC_Object* a_owner )
//-------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , ID( MAC::bad_index() )
   , SET( 0 )
{
}

//-------------------------------------------------------------------------
MAC_IndexSet*
MAC_IndexSet:: create( MAC_Object* a_owner,
                       size_t_vector const& vec,
                       size_t a_id ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IndexSet:: create" ) ;
   
   MAC_IndexSet* result = new MAC_IndexSet( a_owner, vec, a_id ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->elements().size() == vec.size() ) ;
   MAC_CHECK_POST( FORALL( (size_t i=0;i<vec.size()-1;i++),
                           result->elements()(i)<=result->elements()(i+1) )) ;
   MAC_CHECK_POST( FORALL( (size_t i=0;i<vec.size();i++),
                           vec.has(result->elements()(i) )) ) ;
   MAC_CHECK_POST( result->id() == a_id ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
MAC_IndexSet:: MAC_IndexSet( MAC_Object* a_owner,
                             size_t_vector const& vec,
                             size_t a_id )
//-------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , ID( a_id )
   , SET( vec )
{
   MAC_LABEL( "MAC_IndexSet:: MAC_IndexSet" ) ;
   SET.sort_increasingly() ;
}

//-------------------------------------------------------------------------
MAC_IndexSet:: ~MAC_IndexSet( void )
//-------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
MAC_IndexSet:: re_initialize( size_t_vector const& vec, 
                              size_t a_id )
//----------------------------------------------------------------------
{
   SET = vec ;
   SET.sort_increasingly() ;
   ID = a_id ;
}

//----------------------------------------------------------------------
size_t
MAC_IndexSet:: id( void ) const
//----------------------------------------------------------------------
{
   return( ID ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
MAC_IndexSet:: elements( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IndexSet:: elements" ) ;
   
   return( SET ) ;
}

//----------------------------------------------------------------------
bool
MAC_IndexSet:: is_equal( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IndexSet:: is_equal" ) ;
   MAC_CHECK_PRE( is_equal_PRE( other ) ) ;
   
   bool result = ( three_way_comparison( other ) == 0 ) ;
   
   MAC_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
int
MAC_IndexSet:: three_way_comparison( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IndexSet:: three_way_comparison" ) ;
   MAC_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   MAC_CHECK( dynamic_cast<MAC_IndexSet const*>( other ) != 0 ) ;
   MAC_IndexSet const* oo = static_cast<MAC_IndexSet const*>( other ) ;
                  
   size_t const ss = SET.size() ;
   size_t const os = oo->SET.size() ;

   int result = 0 ;
   if( ss<os )
   {
      result = -1 ;
   }
   else if( ss>os )
   {
      result = 1 ;
   }
   else if( ss>0 )
   {
      for( size_t i=0 ; i<ss ; i++ )
      {
         size_t si = SET(i) ;
         size_t osi = oo->SET(i) ;
         
         if( si < osi )
         {
            result = -1 ;
            break ;
         }
         else if( si > osi )
         {
            result = 1 ;
            break ;
         }
      }
      
   }
   
   MAC_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
MAC_IndexSet:: hash_code( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IndexSet:: hash_code" ) ;
   size_t result = 0 ;
   for( size_t i=0 ; i<SET.size() ; i++ )
   {
      result += SET(i) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_IndexSet:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IndexSet:: print" ) ;

   std::string const s( indent_width, ' ' ) ;
   os << s << "index: " << ID << ",  " << SET ;
}
