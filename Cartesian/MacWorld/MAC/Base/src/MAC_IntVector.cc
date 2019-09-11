#include <MAC_IntVector.hh>

#include <iostream>

#include <MAC_assertions.hh>
#include <MAC_Container.hh>
#include <MAC_Iterator.hh>
#include <intVector.hh>

//----------------------------------------------------------------------------
MAC_IntVector*
MAC_IntVector:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IntVector:: create_clone" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_IntVector* result = 0 ;
   result = new MAC_IntVector( a_owner, dv ) ;
   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_IntVector*
MAC_IntVector:: create( MAC_Object* a_owner,
                        intVector const& aIntVector )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IntVector:: create(intVector)" ) ;
   MAC_IntVector* result = new MAC_IntVector( a_owner, aIntVector ) ;
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->to_int_vector() == aIntVector ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_IntVector:: MAC_IntVector( MAC_Object* a_owner,
                                 intVector const& aIntVector )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner )
   , dv( aIntVector )
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_IntVector:: MAC_IntVector( MAC_Object* a_owner,
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
                 val->data_type()!=Int ||
                 !val->value_can_be_evaluated( 0 ) ) ;
      if( !error ) dv(cpt++) = val->to_int() ;
   }
   
   list_iterator->destroy() ; list_iterator=0 ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_IntVector*
MAC_IntVector:: create( MAC_Object* a_owner,
                        MAC_Container const* aDataList ) 
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IntVector:: create(MAC_Container)" ) ;
   MAC_CHECK_PRE( aDataList!=0 ) ;
   MAC_CHECK_PRE( aDataList->count()!=0 ) ;

   bool error = false ;
   MAC_IntVector* result = new MAC_IntVector( 0, aDataList, error ) ;
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
   MAC_CHECK_POST( IMPLIES( result != 0, result->to_int_vector().size() == aDataList->count() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_IntVector:: ~MAC_IntVector( void )
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
MAC_IntVector:: set( intVector const& other )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IntVector:: set" ) ;
   dv=other ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<other.size() ; i++ ),
                           to_int_vector()(i)==other(i) ) ) ;
}

//----------------------------------------------------------------------------
MAC_Data::Type
MAC_IntVector:: data_type( void ) const
//----------------------------------------------------------------------------
{
   return( IntVector ) ;
}

//----------------------------------------------------------------------------
intVector const& 
MAC_IntVector:: to_int_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IntVector:: to_int_vector" ) ;
   MAC_CHECK_PRE( to_int_vector_PRE(ct) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( dv ) ;
}

//----------------------------------------------------------------------------
void
MAC_IntVector:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IntVector:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   MAC_CHECK_INV( invariant() ) ;
   os << space << "< " ;
   for( size_t i=0 ; i<dv.size() ; i++ )
   {
       os << dv(i) << " " ;
   }
   os << " > " ;
}

//-------------------------------------------------------------------------
bool
MAC_IntVector:: invariant( void ) const 
//-------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Data::invariant() ) ;
   return true ;
}
