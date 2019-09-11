#include <MAC_BoolVector.hh>

#include <iostream>

#include <MAC_assertions.hh>
#include <MAC_Container.hh>
#include <MAC_Iterator.hh>
#include <boolVector.hh>

//----------------------------------------------------------------------------
MAC_BoolVector*
MAC_BoolVector:: create( MAC_Object* a_owner,
                         boolVector const& aBoolVector )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BoolVector:: create(boolVector)" ) ;

   MAC_BoolVector* result = new MAC_BoolVector( a_owner, aBoolVector ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->to_bool_vector() == aBoolVector ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_BoolVector:: MAC_BoolVector( MAC_Object* a_owner,
                                 boolVector const& aBoolVector )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner )
   , dv( aBoolVector )
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_BoolVector*
MAC_BoolVector:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BoolVector:: create_clone" ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_BoolVector* result = new MAC_BoolVector( a_owner, dv ) ;

   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_BoolVector:: MAC_BoolVector( MAC_Object* a_owner,
                                 MAC_Container const* aDataList,
                                 bool& error )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner )
   , dv( aDataList->count() )
{
   MAC_Iterator* list_iterator = aDataList->create_iterator( 0 ) ;
   size_t cpt=0 ;
   for( list_iterator->start() ;
        list_iterator->is_valid() && !error ;
        list_iterator->go_next() )
   {
      MAC_Data const* val =
         dynamic_cast<MAC_Data const*>( list_iterator->item() ) ;
      error |= ( val==0  ||
                 val->data_type()!=Bool ||
                 !val->value_can_be_evaluated( 0 ) ) ;
      if( !error ) dv(cpt++) = val->to_bool() ;
   }
   
   list_iterator->destroy() ; list_iterator=0 ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_BoolVector*
MAC_BoolVector:: create( MAC_Object* a_owner,
                         MAC_Container const* aDataList ) 
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BoolVector:: create(MAC_Container)" ) ;
   MAC_CHECK_PRE( aDataList!=0 ) ;
   MAC_CHECK_PRE( aDataList->count()!=0 ) ;

   bool error = false ;
   MAC_BoolVector* result = new MAC_BoolVector( 0, aDataList, error ) ;
   if( error )
   {
      result->destroy() ;
      result = 0 ;
   }
   else if( a_owner != 0 )
   {
      result->set_owner( a_owner ) ;
   }

   MAC_CHECK_POST( IMPLIES( result != 0, result->owner() == a_owner ) ) ;
   MAC_CHECK_POST( IMPLIES( result != 0, result->to_bool_vector().size() == aDataList->count() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_BoolVector:: ~MAC_BoolVector( void )
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
MAC_BoolVector:: set( boolVector const& other )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BoolVector:: set" ) ;
   dv=other ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<other.size() ; i++ ),
                           to_bool_vector()(i)==other(i) ) ) ;
}

//----------------------------------------------------------------------------
MAC_Data::Type
MAC_BoolVector:: data_type( void ) const
//----------------------------------------------------------------------------
{
   return( BoolVector ) ;
}

//----------------------------------------------------------------------------
boolVector const& 
MAC_BoolVector:: to_bool_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BoolVector:: to_bool_vector" ) ;
   MAC_CHECK_PRE( to_bool_vector_PRE(ct) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( dv ) ;
}

//----------------------------------------------------------------------------
void
MAC_BoolVector:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BoolVector:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   MAC_CHECK_INV( invariant() ) ;
   os << space << "< " ;
   for( size_t i=0 ; i<dv.size() ; i++ )
   {
      if( (i+1) % 10 == 0 )
      {
         os << std::endl << space ;
      }
      if( dv(i) )
      {
         os << "true " ;
      }
      else
      {
         os << "false " ;
      }
   }
   os << " > " << std::endl ;
}

//-------------------------------------------------------------------------
bool
MAC_BoolVector:: invariant( void ) const 
//-------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Data::invariant() ) ;
   return( true ) ;
}
