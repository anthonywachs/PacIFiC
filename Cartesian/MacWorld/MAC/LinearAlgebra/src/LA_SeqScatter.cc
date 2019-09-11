#include <LA_SeqScatter.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_DistributedPartition.hh>

#include <LA_SeqImplementation.hh>
#include <LA_SeqVector.hh>

//----------------------------------------------------------------------
LA_SeqScatter*
LA_SeqScatter:: create( MAC_Object* a_owner,
                        size_t a_nb_rows,
                        size_t_vector const& a_repatriated_items_table,
                        size_t_vector const& a_local_indices_table )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqScatter:: create" ) ;
   MAC_CHECK_PRE( a_repatriated_items_table.size() ==
                                        a_local_indices_table.size() ) ;
   MAC_CHECK_PRE(
      FORALL(
         ( size_t i=0 ; i<a_repatriated_items_table.size() ; ++i ),
         a_repatriated_items_table(i)<a_nb_rows ) ) ;

   LA_SeqScatter* result =
      new LA_SeqScatter( a_owner,
                         a_nb_rows,
                         a_repatriated_items_table, a_local_indices_table ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->size() == a_repatriated_items_table.size() ) ;
   MAC_CHECK_POST( result->repatriated_items() == a_repatriated_items_table ) ;
   MAC_CHECK_POST( result->local_indices() == a_local_indices_table ) ;
   MAC_CHECK_POST( result->distribution()->global_number() == a_nb_rows ) ;
   MAC_CHECK_POST( result->distribution()->local_number() == a_nb_rows ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqScatter:: LA_SeqScatter(
                     MAC_Object* a_owner,
                     size_t a_nb_rows,
                     size_t_vector const& a_repatriated_items_table,
                     size_t_vector const& a_local_indices_table )
//----------------------------------------------------------------------
   : LA_Scatter( a_owner )
   , NB_ROWS( a_nb_rows )
   , NEEDED( a_repatriated_items_table )
   , LOCAL( a_local_indices_table )
   , DIST( 0 )
{
}

//----------------------------------------------------------------------
LA_SeqScatter:: ~LA_SeqScatter( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_Implementation const*
LA_SeqScatter:: implementation( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqScatter:: implementation" ) ;
   
   static LA_Implementation const* result = LA_SeqImplementation::object() ;

   MAC_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_SeqScatter:: size( void ) const
//----------------------------------------------------------------------
{
   return( NEEDED.size() ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
LA_SeqScatter:: repatriated_items( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqScatter:: repatriated_items" ) ;
   
   size_t_vector const& result = NEEDED ;

   MAC_CHECK_POST( repatriated_items_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
LA_SeqScatter:: local_indices( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqScatter:: local_indices" ) ;
   
   size_t_vector const& result = LOCAL ;

   MAC_CHECK_POST( local_indices_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DistributedPartition const*
LA_SeqScatter:: distribution( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqScatter:: distribution" ) ;

   if( DIST == 0 )
   {
      DIST = MAC_DistributedPartition::create(
                                const_cast<LA_SeqScatter*>( this ) ) ;
      DIST->set_global_number( NB_ROWS ) ;
   }
   
   MAC_DistributedPartition const* result = DIST ;

   MAC_CHECK_POST( distribution_POST( result ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
void
LA_SeqScatter:: get( LA_Vector const* source,
                     LA_SeqVector* dest ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqScatter:: get" ) ;
   MAC_CHECK_PRE( get_PRE( source, dest) ) ;

   LA_SeqVector const* bsource = static_cast<LA_SeqVector const*>( source ) ;
   MAC_CHECK( dynamic_cast<LA_SeqVector const*>( source ) != 0 ) ;

   size_t n = size() ;
   
   for( size_t i=0 ; i<n ; ++i )
   {
      dest->set_item( LOCAL(i), bsource->item( NEEDED(i) ) ) ;
   }

   dest->synchronize() ;
   MAC_CHECK_POST( get_POST( source, dest) ) ;
}

//----------------------------------------------------------------------
void
LA_SeqScatter:: set( LA_SeqVector const* source,
                     LA_Vector* dest ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqScatter:: set" ) ;
   MAC_CHECK_PRE( set_PRE( source, dest) ) ;

   MAC_CHECK( dynamic_cast<LA_SeqVector const*>( dest ) != 0 ) ;

   size_t n = size() ;
   
   for( size_t i=0 ; i<n ; ++i )
   {
      dest->set_item( NEEDED(i), source->item( LOCAL(i) ) ) ;
   }
   
   dest->synchronize() ;
   MAC_CHECK_POST( set_POST( source, dest) ) ;

}

//----------------------------------------------------------------------
bool
LA_SeqScatter:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Scatter::implementation_POST( result ) ) ;
   MAC_ASSERT( result == LA_SeqImplementation::object() ) ;
   return( true ) ;
}
