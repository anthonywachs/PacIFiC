#include <LA_Scatter.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_DistributedPartition.hh>

#include <LA_SeqVector.hh>

#include <size_t_vector.hh>

//----------------------------------------------------------------------
LA_Scatter:: LA_Scatter( MAC_Object* a_owner )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
{
}

//----------------------------------------------------------------------
LA_Scatter:: ~LA_Scatter( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_Scatter:: implementation_POST( LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: repatriated_items_POST( size_t_vector const& result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result.size() == size() ) ;
   MAC_ASSERT(
      FORALL(
         ( size_t i=0 ; i<size() ; ++i ),
               result(i)<distribution()->global_number() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: local_indices_POST( size_t_vector const& result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result.size() == size() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: distribution_POST( MAC_DistributedPartition const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == this ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: get_PRE( LA_Vector const* source,
                      LA_SeqVector const* dest ) const 
//----------------------------------------------------------------------
{
   MAC_ASSERT( source != 0 ) ;
   MAC_ASSERT( source->is_synchronized() ) ;
   MAC_ASSERT( source->implementation() == implementation() ) ;
   MAC_ASSERT( source->row_distribution()->is_compatible( distribution() ) ) ;
   MAC_ASSERT( dest!=0 ) ;
   MAC_ASSERT(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              repatriated_items()(i) < source->nb_rows() ) ) ;
   MAC_ASSERT(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              local_indices()(i) < dest->nb_rows() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: get_POST( LA_Vector const* source,
                       LA_SeqVector const* dest ) const 
//----------------------------------------------------------------------
{
   MAC_ASSERT( dest->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: set_PRE( LA_SeqVector const* source,
                      LA_Vector const* dest ) const 
//----------------------------------------------------------------------
{
   MAC_ASSERT( source!=0 ) ;
   MAC_ASSERT( source->is_synchronized() ) ;
   MAC_ASSERT( dest!=0 ) ;
   MAC_ASSERT( dest->implementation() == implementation() ) ;
   MAC_ASSERT( dest->row_distribution()->is_compatible( distribution() ) ) ;
   MAC_ASSERT(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              local_indices()(i) < source->nb_rows() ) ) ;
   MAC_ASSERT(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              repatriated_items()(i) < dest->nb_rows() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Scatter:: set_POST( LA_SeqVector const* source,
                       LA_Vector const* dest ) const 
//----------------------------------------------------------------------
{
   MAC_ASSERT( dest->is_synchronized() ) ;
   return( true ) ;
}


