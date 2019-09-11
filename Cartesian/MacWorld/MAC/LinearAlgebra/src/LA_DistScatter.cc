#include <LA_DistScatter.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <doubleVector.hh>

#include <LA_CRSmatrix.hh>
#include <LA_DistImplementation.hh>
#include <LA_DistMatrix.hh>
#include <LA_DistVector.hh>
#include <LA_MatrixIterator.hh>
#include <LA_PelMatrix.hh>

//----------------------------------------------------------------------
LA_DistScatter*
LA_DistScatter:: create( MAC_Object* a_owner,
                         MAC_DistributedPartition const* partition )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: create" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( partition != 0 ) ;

   LA_DistScatter* result = new LA_DistScatter( a_owner, partition ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->distribution()->is_compatible( partition ) ) ;
   MAC_CHECK_POST( result->size() == 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DistScatter:: LA_DistScatter( MAC_Object* a_owner,
                                 MAC_DistributedPartition const* partition )
//----------------------------------------------------------------------
   : LA_Scatter( a_owner )
   , DIST( MAC_DistributedPartition::create( this ) )
   , COMM( MAC_Exec::communicator() )
   , RECEIVED_NB( 0 )
   , SENT_NB( 0 )
   , RECEIVED_IDX( 0 )
   , SENT_IDX( 0 )
   , ITEMS( 0 )
   , LOCAL_IDS( 0 )
   , REQUEST( 0 )
{
   MAC_LABEL( "LA_DistScatter:: LA_DistScatter" ) ;
   re_initialize( partition ) ;
}

//----------------------------------------------------------------------
LA_DistScatter:: ~LA_DistScatter( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_DistScatter:: re_initialize( MAC_DistributedPartition const* partition )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: re_initialize" ) ;
   MAC_CHECK_PRE( partition != 0 ) ;

   DIST->set( partition ) ;
   ITEMS.resize(0) ;
   LOCAL_IDS.resize(0) ;

   MAC_CHECK_POST( distribution()->is_compatible( partition ) ) ;
   MAC_CHECK_POST( size() == 0 ) ;
}

//----------------------------------------------------------------------
LA_Implementation const*
LA_DistScatter:: implementation( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: implementation" ) ;

   LA_Implementation const* result = LA_DistImplementation::object() ;

   MAC_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_DistScatter:: size( void ) const
//----------------------------------------------------------------------
{
   return( ITEMS.size() ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
LA_DistScatter:: repatriated_items( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: repatriated_items" ) ;

   size_t_vector const& result = ITEMS ;

   MAC_CHECK_POST( repatriated_items_POST( result ) ) ;
   MAC_CHECK_POST(
      FORALL( ( size_t i=1 ; i<result.size() ; ++i ),
              result(i)>result(i-1) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
LA_DistScatter:: local_indices( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: local_indices" ) ;

   size_t_vector const& result = LOCAL_IDS ;

   MAC_CHECK_POST( local_indices_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DistributedPartition const*
LA_DistScatter:: distribution( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: distribution" ) ;

   MAC_DistributedPartition const* result = DIST ;

   MAC_CHECK_POST( distribution_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistScatter:: set_sorted( size_t_vector const& a_repatriated_items_table,
                             size_t_vector const& a_local_indices_table )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: set_sorted" ) ;
   MAC_CHECK_PRE(
      FORALL(
         ( size_t i=0 ; i<a_repatriated_items_table.size() ; ++i ),
         a_repatriated_items_table(i)<distribution()->global_number() ) ) ;
   MAC_CHECK_PRE(
      FORALL(
         ( size_t i=1 ; i<a_repatriated_items_table.size() ; ++i ),
         a_repatriated_items_table(i)>a_repatriated_items_table(i-1) ) ) ;
   MAC_CHECK_PRE(
      a_local_indices_table.size()==a_repatriated_items_table.size() ) ;

   ITEMS = a_repatriated_items_table ;
   LOCAL_IDS = a_local_indices_table ;

   RECEIVED_NB.re_initialize( COMM->nb_ranks() ) ;
   SENT_NB.re_initialize( COMM->nb_ranks() ) ;

   if( RECEIVED_IDX.size()<a_repatriated_items_table.size() )
   {
      RECEIVED_IDX.resize( a_repatriated_items_table.size() ) ;
   }

   intVector const& start = DIST->start_of_partition() ;
   intVector const& partitioning = DIST->partitioning() ;

   size_t const rank = COMM->rank() ;
   int const my_start = start( rank )  ;
   int const my_up = my_start + partitioning( rank ) ;

   for( size_t i=0 ; i<ITEMS.size() ; i++ )
   {
      size_t idx = ITEMS( i ) ;
      RECEIVED_IDX( i ) = idx ;//???S : recopie ITEMS

      size_t curr = DIST->rank_of( idx ) ;
      RECEIVED_NB(curr )++ ;

      MAC_CHECK( idx >= (size_t) start( curr ) &&
                 idx <  (size_t) (start( curr )+partitioning(curr)) ) ;
   }

   COMM->all_to_all( RECEIVED_NB, SENT_NB ) ;
   SENT_NB( rank ) = 0 ;

   size_t send_size = SENT_NB.sum() ;

   //if( SENT_IDX.size() < send_size )
   SENT_IDX.resize( send_size ) ;

   size_t N_send_tot = 0 ;
   size_t N_rec_tot = 0 ;
   for( size_t r=0 ; r<COMM->rank() ; r++ )
   {
      size_t N = SENT_NB( r ) ;

      if( N>0 )
      {
         receive_index_from( N_rec_tot, r ) ;
         N_rec_tot += N ;
      }
   }

   for( size_t r=0 ; r<COMM->nb_ranks() ; r++ )
   {
      size_t N = RECEIVED_NB( r ) ;

      if( N>0 )
      {
         if( r!=COMM->rank() )
         {
            send_index_to( N_send_tot, r ) ;
         }
         N_send_tot += N ;
      }
   }
   for( size_t r=COMM->rank()+1 ; r<COMM->nb_ranks() ; r++ )
   {
      size_t N = SENT_NB( r ) ;
      if( N>0 )
      {
         receive_index_from( N_rec_tot, r ) ;
         N_rec_tot += N ;
      }
   }

   MAC_ASSERT( (size_t)my_start == DIST->first_local_index() ) ;
   MAC_ASSERT( (size_t)my_up == DIST->local_index_limit() ) ;
   MAC_ASSERT( FORALL( ( size_t i=0 ; i<send_size ; ++i ),
                        SENT_IDX(i)>=my_start && SENT_IDX(i)<my_up ) ) ;


   MAC_CHECK( FORALL( ( size_t i=0 ; i<send_size ; ++i ),
                      SENT_IDX(i)>=my_start && SENT_IDX(i)<my_up ) ) ;
   MAC_CHECK( N_rec_tot == send_size ) ;

   MAC_CHECK_POST( size() == a_repatriated_items_table.size() ) ;
   MAC_CHECK_POST( repatriated_items() == a_repatriated_items_table ) ;
   MAC_CHECK_POST( local_indices() == a_local_indices_table ) ;
}

//----------------------------------------------------------------------
void
LA_DistScatter:: set_unsorted( size_t_vector const& a_repatriated_items_table,
                               size_t_vector const& a_local_indices_table )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: set_unsorted" ) ;
   MAC_CHECK_PRE(
      FORALL(
         ( size_t i=0 ; i<a_repatriated_items_table.size() ; ++i ),
         a_repatriated_items_table(i)<distribution()->global_number() ) ) ;
   MAC_CHECK_PRE(
      a_repatriated_items_table.size() == a_local_indices_table.size() ) ;

   size_t N = a_repatriated_items_table.size() ;

   size_t_vector sorted_items( N ) ;
   size_t_vector sorted_indexes( N ) ;
   //global_to_local
   size_t_vector tmp( DIST->global_number(), MAC::bad_index() ) ;

   for( size_t i=0 ; i<N ; i++ )
   {
      tmp( a_repatriated_items_table(i) ) = a_local_indices_table(i) ;
   }
   size_t idx=0 ;

   size_t GN = DIST->global_number() ;

   for( size_t i=0 ; i<GN ; i++ )
   {
      if( tmp(i) != MAC::bad_index() )
      {
         sorted_items(idx) = i ;
         sorted_indexes(idx) = tmp(i) ;
         idx++ ;
      }
   }
   MAC_CHECK( idx==N ) ;

   /*
    * Cela trie les indices globaux par ordre croissant.
    * La fonction set_sorted utilise le fait que les indices sont rangé par
    * blocs correspondant au proc dans l'ordre croissant.
    */



   set_sorted( sorted_items, sorted_indexes ) ;

   MAC_CHECK_POST( size() == a_repatriated_items_table.size() ) ;
   MAC_CHECK_POST(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              repatriated_items().has( a_repatriated_items_table(i) ) ) ) ;
   MAC_CHECK_POST(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              local_indices()(i) == a_local_indices_table(
                 a_repatriated_items_table.index_of( repatriated_items()(i) ) ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistScatter:: get( LA_Vector const* source,
                      LA_SeqVector* dest ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: get" ) ;
   MAC_CHECK_PRE( get_PRE( source, dest) ) ;

   LA_DistVector const* dsource = static_cast<LA_DistVector const*>( source ) ;
   MAC_CHECK( dynamic_cast<LA_DistVector const*>( source ) != 0 ) ;

   size_t N_rec_tot = 0 ;
   size_t N_send_tot = 0 ;
   for( size_t r=0 ; r<COMM->rank() ; r++ )
   {
      receive_values_from( dest, N_rec_tot, r ) ;
      N_rec_tot += RECEIVED_NB( r ) ;
   }
   for( size_t r=0 ; r<COMM->nb_ranks() ; r++ )
   {
      if( r!=COMM->rank() )
         send_values_to( dsource, N_send_tot, r ) ;
      N_send_tot += SENT_NB( r ) ;
   }

   size_t n = RECEIVED_NB( COMM->rank() ) ;
   for( size_t j=0 ; j<n ; j++ )
   {
      size_t idx = RECEIVED_IDX(j+N_rec_tot) ;
      dest->set_item( LOCAL_IDS( j+N_rec_tot ), dsource->item(idx) ) ;
   }
   N_rec_tot += n ;

   for( size_t r=COMM->rank()+1 ; r<COMM->nb_ranks() ; r++ )
   {
      receive_values_from( dest, N_rec_tot, r ) ;
      N_rec_tot += RECEIVED_NB( r ) ;
   }
   dest->synchronize() ;

   MAC_CHECK_POST( get_POST( source, dest ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistScatter:: get_begin( LA_DistVector const* source ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: get_begin" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( source!=0 ) ;
   MAC_CHECK_PRE( source->is_synchronized() ) ;
   MAC_CHECK_PRE( distribution()->is_compatible(
                                        source->row_distribution() ) ) ;
   MAC_CHECK_PRE(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              repatriated_items()(i) < source->nb_rows() ) ) ;
   MAC_CHECK_PRE( !is_getting() ) ;

   static doubleVector BUFFER(0) ;

   size_t N_send_tot = 0 ;
   size_t req = SENT_NB.sum() ;
   if( BUFFER.size()<req ) BUFFER.resize( req ) ;
   double * ptr = const_cast<double*>( BUFFER.data() ) ;
   REQUEST = new void * [ COMM->nb_ranks() ] ;

   for( size_t i=0 ; i<COMM->nb_ranks() ; i++ )
   {
      if( i!= COMM->rank() )
      {
         REQUEST[i] = send_values_NB_to( source, N_send_tot, i, ptr ) ;
         ptr += SENT_NB(i) ;
         N_send_tot += SENT_NB(i) ;
      }
      else
      {
         REQUEST[i]  = 0 ;
      }
   }

   MAC_CHECK_POST( is_getting() ) ;
}

//----------------------------------------------------------------------
void
LA_DistScatter:: get_end( LA_SeqVector* dest ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: get_end" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( dest!=0 ) ;
   MAC_CHECK_PRE(
      FORALL( ( size_t i=0 ; i<size() ; ++i ),
              local_indices()(i) < dest->nb_rows() ) ) ;
   MAC_CHECK_PRE( is_getting() ) ;

   size_t N_rec_tot = 0 ;
   for( size_t i=0 ; i<COMM->nb_ranks() ; i++ )
   {
      if( i!=COMM->rank() )
      {
         receive_values_from( dest, N_rec_tot, i ) ;
      }
      N_rec_tot += RECEIVED_NB(i) ;
   }
   for( size_t i=0 ; i<COMM->nb_ranks() ; i++ )
   {
      if( REQUEST[i] != 0 ) COMM->wait( REQUEST[i] ) ;
   }
   dest->synchronize() ;

   delete [] REQUEST ;
   REQUEST = 0 ;

   MAC_CHECK_POST( !is_getting() ) ;
   MAC_CHECK_POST( dest->is_synchronized() ) ;
}

//----------------------------------------------------------------------
bool
LA_DistScatter:: is_getting( void ) const
//----------------------------------------------------------------------
{
   return( REQUEST != 0 ) ;
}

//----------------------------------------------------------------------
void
LA_DistScatter:: set( LA_SeqVector const* source,
                      LA_Vector* dest ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: set" ) ;
   MAC_CHECK_PRE( set_PRE( source, dest) ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const*>( dest ) != 0 ) ;
   LA_DistVector * ddest = static_cast<LA_DistVector *>( dest ) ;

   for( size_t i=0 ; i<ITEMS.size() ; i++ )
   {
      size_t idx = ITEMS(i) ;
      dest->set_item( idx, source->item( LOCAL_IDS( i) ) ) ;
   }
   ddest->synchronize() ;

   MAC_CHECK_POST( set_POST( source, dest) ) ;

}
//----------------------------------------------------------------------
void
LA_DistScatter:: get_rows( LA_DistMatrix const* source,
                           LA_SeqMatrix* dest ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: get_rows" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( source!=0 ) ;
   MAC_CHECK_PRE( source->is_synchronized() ) ;
   MAC_CHECK_PRE( distribution()->is_compatible(
                                       source->row_distribution() ) ) ;
   MAC_CHECK_PRE( dest!=0 ) ;
   MAC_CHECK_PRE( dest->state() != LA::NotSync_set ) ;
   MAC_CHECK_PRE( repatriated_items() == local_indices() ) ;

   LA_MatrixIterator* it = source->create_stored_item_iterator( 0 ) ;

   size_t N_rec_tot = 0 ;
   size_t N_send_tot = 0 ;
   for( size_t i=0 ; i<COMM->rank() ; i++ )
   {
      receive_rows_from( dest, N_rec_tot, i ) ;
   }

   size_t local_to_local = RECEIVED_NB(COMM->rank()) ;

   for( size_t i=0 ; i<local_to_local ; i++ )
   {
      size_t idx = ITEMS( N_rec_tot ) ;
      size_t jdx = LOCAL_IDS( N_rec_tot ) ;
      for( it->start_row_items( idx ) ; it->is_valid() ; it->go_next() )
      {
         dest->set_item( jdx, it->col(), it->item() ) ;
      }
      N_rec_tot++ ;
   }
   dest->synchronize() ;

   for( size_t i=0 ; i<COMM->nb_ranks() ; i++ )
   {
      send_rows_to( it, source, N_send_tot, i ) ;
   }
   for( size_t i=COMM->rank()+1 ; i<COMM->nb_ranks() ; i++ )
   {
      receive_rows_from( dest, N_rec_tot, i ) ;
   }
   dest->synchronize() ;

   it->destroy() ; it = 0 ;

   MAC_CHECK_POST( dest->is_synchronized() ) ;
}

//----------------------------------------------------------------------
void
LA_DistScatter:: send_values_to( LA_DistVector const* source,
                                 size_t start,
                                 size_t rank ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: send_values_to" ) ;
   MAC_CHECK( rank < COMM->nb_ranks() ) ;
   MAC_CHECK( rank != COMM->rank() ) ;

   static doubleVector tmp( 0 ) ;

   size_t N = SENT_NB( rank ) ;
   if( N > 0 )
   {
      if( tmp.size() < N ) tmp.resize( N ) ;
      for( size_t j=0 ; j<N ; j++ )
         tmp( j ) = source->item( SENT_IDX( start+j ) ) ;
      COMM->send( rank, tmp.data(), N ) ;
   }
}

//----------------------------------------------------------------------
void
LA_DistScatter:: send_index_to( size_t start, size_t rank ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: send_index_to" ) ;
   MAC_CHECK( rank < COMM->nb_ranks() ) ;

   size_t N = RECEIVED_NB( rank ) ;

   intVector waited_index( N ) ;
   for( size_t j=0 ; j<N ; j++ )
   {
      waited_index(j) = RECEIVED_IDX( j + start ) ;
   }
   COMM->send( rank, waited_index.data(), N ) ;
}

//----------------------------------------------------------------------
void
LA_DistScatter:: receive_index_from( size_t start, size_t rank )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: receive_index_from" ) ;
   MAC_CHECK( rank < COMM->nb_ranks() ) ;

   size_t N = SENT_NB( rank ) ;
   intVector sended_index( N ) ;
   COMM->receive( rank, const_cast<int*>(sended_index.data()), N ) ;
   for( size_t j=0 ; j<N ; j++ )
   {
      SENT_IDX( start+j ) = sended_index( j ) ;
   }
}


//----------------------------------------------------------------------
void*
LA_DistScatter:: send_values_NB_to( LA_DistVector const* source,
                                    size_t start,
                                    size_t rank,
                                    double * ptr ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: send_values_NB_to" ) ;
   MAC_CHECK( rank < COMM->nb_ranks() ) ;
   MAC_CHECK( rank != COMM->rank() ) ;

   void* result = 0 ;

   size_t N = SENT_NB( rank ) ;
   if( N>0 )
   {
      for( size_t j=0 ; j<N ; j++ )
         ptr[j] = source->item( SENT_IDX( start+j ) ) ;
      result = COMM->Isend( rank, ptr, N ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistScatter:: receive_values_from( LA_SeqVector* dest,
                                      size_t start,
                                      size_t rank ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: receive_values_from" ) ;
   MAC_CHECK( rank < COMM->nb_ranks() ) ;
   MAC_CHECK( rank != COMM->rank() ) ;

   static doubleVector tmp( 0 ) ;
   size_t N = RECEIVED_NB( rank ) ;

   if( N>0 )
   {
      if( tmp.size() < N ) tmp.resize( N ) ;
      double* ptr_tmp = const_cast<double*>( tmp.data() ) ;
      COMM->receive( rank, ptr_tmp, N ) ;
      double* ptr_dest = dest->data() ;
      for( size_t j=0 ; j<N ; ++j )
      {
         ptr_dest[ LOCAL_IDS( j+start ) ] = ptr_tmp[j] ;
      }
   }
}

//----------------------------------------------------------------------
void
LA_DistScatter:: send_rows_to( LA_MatrixIterator* it,
                               LA_DistMatrix const* source,
                               size_t& start,
                               size_t rank ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: send_rows_to" ) ;
   MAC_CHECK( source != 0 ) ;
   MAC_CHECK( it != 0 ) ;
   MAC_CHECK( it->matrix() == source ) ;
   MAC_CHECK( rank < COMM->nb_ranks() ) ;

   size_t N = SENT_NB( rank ) ;
   if( N>0 )
   {
      LA_PelMatrix* pmat =
         LA_PelMatrix::create( 0, source->nb_rows(), source->nb_cols() ) ;
      for( size_t j=0 ; j<N ; ++j )
      {
         for( it->start_row_items( SENT_IDX( start+j ) ) ;
              it->is_valid() ; it->go_next() )
         {
            pmat->set_item( it->row(), it->col(), it->item() ) ;
         }
      }
      pmat->synchronize() ;

      LA_CRSmatrix* cmat = LA_CRSmatrix::create( 0, pmat ) ;
      pmat->destroy() ;

      cmat->send( COMM, rank, false ) ;
      cmat->destroy() ;

   }
   start += N ;
}

//----------------------------------------------------------------------
void
LA_DistScatter:: receive_rows_from( LA_SeqMatrix* dest,
                                    size_t& start,
                                    size_t rank ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistScatter:: receive_rows_from" ) ;
   MAC_CHECK( dest != 0 ) ;
   MAC_CHECK( dest->state() != LA::NotSync_set ) ;
   MAC_CHECK( rank < COMM->nb_ranks() ) ;

   size_t N = RECEIVED_NB( rank ) ;
   if( N>0 )
   {
      LA_CRSmatrix const* mat = LA_CRSmatrix::receive( 0, COMM, rank ) ;
      dest->add_Mat( mat, 1.0 ) ;
      mat->destroy() ;
   }
   start += N ;
}

//----------------------------------------------------------------------
bool
LA_DistScatter:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Scatter::implementation_POST( result ) ) ;
   MAC_ASSERT( result == LA_DistImplementation::object() ) ;
   return( true ) ;
}
