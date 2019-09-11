#include <LA_DistVector.hh>

#include <LA_DistImplementation.hh>

#include <iomanip>
#include <fstream>
#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Int.hh>
#include <MAC_Vector.hh>

#ifdef OUTLINE
   #define inline
   #include <LA_DistVector.icc>
   #undef inline
#endif

//----------------------------------------------------------------------
LA_DistVector*
LA_DistVector:: create( MAC_Object* a_owner,
                        size_t a_nb_rows, size_t a_nb_local_rows,
                        LA::DistributionStrategy dist_strat )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: create" ) ;
   MAC_CHECK_PRE( dist_strat == LA::FromLocalSize ||
                  dist_strat == LA::FromGlobalSize ) ;

   LA_DistVector* result = new LA_DistVector( a_owner,
                                              a_nb_rows, a_nb_local_rows,
                                              dist_strat ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST(
       FORALL( ( size_t i=result->row_distribution()->first_local_index() ;
                 i<result->row_distribution()->local_index_limit() ; ++i ),
                 result->item(i) == 0.0 ) ) ;
   MAC_CHECK_POST( result->distribution_strategy() == dist_strat ) ;
   MAC_CHECK_POST( IMPLIES( result->distribution_strategy() == LA::FromGlobalSize,
                            result->nb_rows() == a_nb_rows  ) ) ;
   MAC_CHECK_POST( IMPLIES( result->distribution_strategy() == LA::FromLocalSize,
                            result->row_distribution()->local_index_limit()
                          - result->row_distribution()->first_local_index()
                          == a_nb_local_rows ) ) ;
   MAC_CHECK_POST( result->state() == LA::Sync ) ;
   MAC_CHECK_POST( result->is_resizable() ) ;
   MAC_CHECK_POST( result->is_synchronized() ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
LA_DistVector:: LA_DistVector( MAC_Object* a_owner,
                               size_t a_nb_rows,
                               size_t a_nb_local_rows,
                               LA::DistributionStrategy dist_strat )
//----------------------------------------------------------------------
   : LA_Vector( a_owner, 0 )
   , DIST( MAC_DistributedPartition::create( this ) )
   , GLOBAL_INDEXES( MAC_Exec::communicator()->nb_ranks(), intVector( 0 ) )
   , GLOBAL_VALUES( MAC_Exec::communicator()->nb_ranks(), doubleVector( 0 ) )
   , LOCAL_VECTOR( LA_SeqVector::create( this, 0 ) )
   , VERBOSE( false )
   , INITIALIZED( false )
   , FIRST( 0 )
   , LAST( 0 )
   , NB_RANKS( MAC_Exec::communicator()->nb_ranks() )
   , HAS_BEEN_NULLIFIED( true )
   , NB_EXTRA(0)
{
   MAC_LABEL( "LA_DistVector:: LA_DistVector" ) ;
   MAC_CHECK_PRE( dist_strat == LA::FromLocalSize ||
                  dist_strat == LA::FromGlobalSize ) ;

   set_distribution_strategy( dist_strat ) ;

   re_initialize( a_nb_rows, a_nb_local_rows ) ;
   INITIALIZED = true ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( distribution_strategy() == dist_strat ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: re_initialize( size_t a_nb_rows, size_t a_nb_local_rows )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: re_initialize" ) ;
   MAC_CHECK_PRE( re_initialize_PRE( a_nb_rows, a_nb_local_rows ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( distribution_strategy() == LA::FromLocalSize ||
       DIST->global_number() != a_nb_rows )
   {
      if( distribution_strategy() == LA::FromLocalSize )
      {
         DIST->set_local_number( a_nb_local_rows ) ;
      }
      else if( distribution_strategy() == LA::FromGlobalSize )
      {
         DIST->distribute_global_number( a_nb_rows ) ;
      }
      LOCAL_VECTOR->re_initialize( DIST->local_number() ) ;
      FIRST = DIST->first_local_index() ;
      LAST  = DIST->local_index_limit() ;
      set_rows_number( DIST->global_number() ) ;
   }
   set( 0.0 ) ;

   MAC_CHECK_POST( re_initialize_POST( a_nb_rows, a_nb_local_rows ) ) ;
}

//----------------------------------------------------------------------
LA_DistVector*
LA_DistVector:: create_vector( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: create_vector" ) ;
   MAC_CHECK_PRE( create_vector_PRE( a_owner ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_DistVector* result =
      new LA_DistVector( a_owner,
                         DIST->global_number(),
                         DIST->local_number(),
                         distribution_strategy() ) ;

   MAC_CHECK_POST( create_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DistScatter*
LA_DistVector:: create_scatter( MAC_Object* a_owner,
                                size_t_vector const& from,
                                size_t_vector const& to ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: create_scatter" ) ;
   MAC_CHECK_PRE( create_scatter_PRE( a_owner, from, to ) ) ;

   MAC_CHECK_INV( invariant() ) ;

   LA_DistScatter* result =
      LA_DistScatter::create( a_owner, row_distribution() ) ;
   result->set_unsorted( from, to ) ;

   MAC_CHECK_POST( create_scatter_POST( result, a_owner, from, to ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector*
LA_DistVector:: create_local_vector( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: create_local_vector" ) ;

   MAC_CHECK_PRE( create_local_vector_PRE( a_owner ) ) ;

   LA_SeqVector* result = LA_SeqVector::create( a_owner, nb_rows() ) ;

   for( size_t i=FIRST ; i<LAST; i++ )
   {
      result->set_item( i, LOCAL_VECTOR->item( i-FIRST ) ) ;
   }
   result->synchronize() ;

   MAC_CHECK_POST( create_local_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DistVector:: ~LA_DistVector( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
LA_Implementation const*
LA_DistVector:: implementation( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: implementation" ) ;

   LA_Implementation const* result = LA_DistImplementation::object() ;


   MAC_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set_synchronized( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: set_synchronized" ) ;
   MAC_CHECK_INV( invariant() ) ;
   NB_EXTRA=0 ;
   for( size_t i=0 ; i<NB_RANKS ; i++ )
   {
      GLOBAL_VALUES[i].re_initialize( 0 ) ;
      GLOBAL_INDEXES[i].re_initialize( 0 ) ;
   }
   LA_Vector::synchronize() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( is_synchronized() ) ;
}

//----------------------------------------------------------------------
bool
LA_DistVector:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: is_desynchronizable" ) ;

   bool result = true ;

   MAC_CHECK_POST( result == true ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set( double value )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: set LA_Vector" ) ;
   MAC_CHECK_PRE( set_PRE( value ) ) ;

   LOCAL_VECTOR->set( value ) ;
   set_synchronized() ;

   MAC_CHECK_POST( set_POST( value ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set( LA_Vector const* a )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: set LA_Vector" ) ;
   MAC_CHECK_PRE( set_PRE( a ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const*>( a ) !=0 ) ;
   LA_DistVector const* da = static_cast<LA_DistVector const*>( a ) ;
   LOCAL_VECTOR->set( da->LOCAL_VECTOR ) ;
   set_synchronized() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_POST( a ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set( LA_SeqVector const* a,
                     size_t_vector const& local_to_global )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: set" ) ;
   MAC_CHECK_PRE( a->nb_rows() == local_to_global.size() ) ;
   MAC_CHECK_PRE( state() != LA::NotSync_add ) ;
   MAC_CHECK_INV( invariant() ) ;

   set_unsynchronized_state( LA::NotSync_set ) ;
   for( size_t i=0 ; i<local_to_global.size() ; i++ )
   {
      set_item( local_to_global(i), a->item(i) ) ;
   }

   MAC_CHECK_POST( state() == LA::NotSync_set ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: add( LA_SeqVector const* a,
                     size_t_vector const& local_to_global )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: add" ) ;
   MAC_CHECK_PRE( a->nb_rows() == local_to_global.size() ) ;
   MAC_CHECK_PRE( state()!=LA::NotSync_set ) ;
   MAC_CHECK_INV( invariant() ) ;

   set_unsynchronized_state( LA::NotSync_add ) ;
   for( size_t i=0 ; i<local_to_global.size() ; i++ )
      if( a->item(i)!=0.0 )
         add_to_item( local_to_global(i), a->item(i) ) ;

   MAC_CHECK_POST( state() == LA::NotSync_add ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: scale( double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: scale" ) ;
   MAC_CHECK_PRE( scale_PRE( alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   LOCAL_VECTOR->scale( alpha ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: sum( LA_Vector const* a, double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: sum LA_Vector" ) ;
   MAC_CHECK_PRE( sum_PRE( a, alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const*>( a ) !=0 ) ;
   LA_DistVector const* da = static_cast<LA_DistVector const*>( a ) ;
   LOCAL_VECTOR->sum( da->LOCAL_VECTOR, alpha ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( sum_POST( a, alpha ) ) ;

}

//----------------------------------------------------------------------
double
LA_DistVector:: dot( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: dot LA_Vector" ) ;
   MAC_CHECK_PRE( dot_PRE( a ) ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const*>( a ) !=0 ) ;
   LA_DistVector const* da = static_cast<LA_DistVector const*>( a ) ;

   double result =
      ( DIST->local_number() > 0 ?
                      LOCAL_VECTOR->dot( da->LOCAL_VECTOR ) : 0.0 ) ;
   result = MAC_Exec::communicator()->sum( result ) ;

   MAC_CHECK_POST( dot_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_DistVector:: two_norm( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: two_norm" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( two_norm_PRE() ) ;

   double result =
      ( DIST->local_number() > 0 ? LOCAL_VECTOR->two_norm() : 0.0 ) ;
   result = MAC::sqrt( MAC_Exec::communicator()->sum( MAC::sqr( result ) ) ) ;

   MAC_CHECK_POST( two_norm_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_DistVector:: max_norm( void ) const
//----------------------------------------------------------------------
// max_norm = max( |v( i )| )
{
   MAC_LABEL( "LA_DistVector:: max_norm" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( max_norm_PRE() ) ;

   double result =
      ( DIST->local_number() > 0 ?
                    LOCAL_VECTOR->max_norm() : -MAC::max_double() ) ;
   result = MAC_Exec::communicator()->max( result ) ;

   MAC_CHECK_POST( max_norm_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: recover_global_vector( LA_SeqVector* vec ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: recover_global_vector" ) ;
   MAC_CHECK_PRE( vec!=0 ) ;
   MAC_CHECK_PRE( vec->nb_rows()==nb_rows() ) ;
   MAC_CHECK_PRE( is_synchronized() ) ;
   MAC_CHECK_INV( invariant() ) ;

   size_t dec = FIRST ;
   size_t first = FIRST ;
   size_t last = LAST ;

   for( size_t i=first ; i<last; i++ )
   {
      vec->set_item( i, LOCAL_VECTOR->item( i - dec ) ) ;
   }

   MAC_Exec::communicator()->all_gather_v( LOCAL_VECTOR->data(),
                                           DIST->local_number(),
                                           const_cast< double* >( vec->data() ),
                                           DIST->partitioning(),
                                           DIST->start_of_partition() ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: synchronize( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: synchronize" ) ;
   MAC_CHECK_PRE( synchronize_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK(
      MAC_Exec::communicator()->boolean_and( state()!=LA::NotSync_add ) ||
      MAC_Exec::communicator()->boolean_and( state()!=LA::NotSync_set ) ) ;

   LA::SyncState effective_mode =
      ( MAC_Exec::communicator()->boolean_and( state() != LA::NotSync_set ) ?
        LA::NotSync_add : LA::NotSync_set ) ;

   if( VERBOSE )
   {
      MAC::out() << "LA_DistVector:: synchronize" << std::endl ;
   }

   size_t const nb_ranks = MAC_Exec::communicator()->nb_ranks() ;
   size_t const rank = MAC_Exec::communicator()->rank() ;

   // C'est la que les choses serieuses vont commencer ...

   // Tout d'abords, on va commencer par s'envoyer le nombre d'elements
   //  a echanger entre process
   intVector exported_items( nb_ranks ) ;
   for( size_t i=0 ; i<nb_ranks ; ++i )
   {
      if( i != rank )
      {
         exported_items( i ) = GLOBAL_VALUES[i].size() ;
      }
   }

   intVector exchanged_items( nb_ranks ) ;
   MAC_Exec::communicator()->all_to_all( exported_items, exchanged_items ) ;

   // Maintenant, on va les envoyer et recevoir
   for( size_t i=0 ; i<rank ; ++i )
   {
      size_t N = (size_t) exchanged_items( i ) ;
      if( N > 0 )
      {
         receive_subvector( N, i, effective_mode ) ;
      }
   }

   for( size_t i=0 ; i<nb_ranks ; ++i )
   {
      int N = exported_items( i ) ;
      if( N > 0 )
      {
         send_subvector( N, i ) ;
      }

   }

   for( size_t i=rank+1 ; i<nb_ranks ; ++i )
   {
      size_t N = (size_t) exchanged_items( i ) ;
      if( N > 0 )
      {
         receive_subvector( N, i, effective_mode ) ;
      }
   }

   set_synchronized() ;

   HAS_BEEN_NULLIFIED = false ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( is_synchronized() ) ;
   MAC_CHECK_POST( local_vector()->is_synchronized() ) ;
}


//----------------------------------------------------------------------
void
LA_DistVector:: send_subvector( size_t N, size_t i ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: send_subvector" ) ;
   MAC_CHECK( N > 0 ) ;
   MAC_CHECK( i < MAC_Exec::communicator()->nb_ranks() ) ;
   MAC_CHECK( i != MAC_Exec::communicator()->rank() ) ;

   if( VERBOSE )
   {
      MAC::out() << "Process " << MAC_Exec::communicator()->rank()
                 << " send " << N << " values to proccess " << i ;
   }

   MAC_Exec::communicator()->send( i, GLOBAL_INDEXES[i].data(), N ) ;
   MAC_Exec::communicator()->send( i, GLOBAL_VALUES[i].data(), N ) ;

   if( VERBOSE )
   {
      MAC::out() << " ... OK " << std::endl ;
      MAC::out().flush() ;
   }
}

//----------------------------------------------------------------------
void
LA_DistVector:: receive_subvector( size_t N, size_t i,
                                   LA::SyncState effective_mode )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: receive_subvector" ) ;
   MAC_CHECK( N > 0 ) ;
   MAC_CHECK( i < MAC_Exec::communicator()->nb_ranks() ) ;
   MAC_CHECK( i != MAC_Exec::communicator()->rank() ) ;

   doubleVector values( N ) ;
   double* ptr_values = const_cast<double*>( values.data() ) ;

   intVector idx( N ) ;
   int* ptr_idx = const_cast<int*>( idx.data() ) ;

   double* ptr_vec = LOCAL_VECTOR->data() ;

   if( VERBOSE )
   {
      MAC::out() << "Process " << MAC_Exec::communicator()->rank()
                 << " receive " << N << " values from " << i ;
   }

   MAC_Exec::communicator()->receive( i, ptr_idx, N ) ;
   MAC_Exec::communicator()->receive( i, ptr_values, N ) ;

   if( effective_mode == LA::NotSync_add )
   {
      for( size_t j=0 ; j<N ; j++ )
      {
         MAC_CHECK( idx(j)>=(int)FIRST && idx(j)<(int)LAST ) ;
         ptr_vec[ptr_idx[j]-FIRST] += ptr_values[j] ;
      }

   }
   else
   {
      for( size_t j=0 ; j<N ; j++ )
      {
         MAC_CHECK( idx(j)>=(int)FIRST && idx(j)<(int)LAST ) ;
         MAC_CHECK(
            IMPLIES( HAS_BEEN_NULLIFIED ,
                     ( LOCAL_VECTOR->item(idx(j)-FIRST)==0.0
                       || LOCAL_VECTOR->item(idx(j)-FIRST)==values(j) ) ) ) ;
         ptr_vec[ptr_idx[j]-FIRST] = ptr_values[j] ;
      }
   }
   if( VERBOSE )
   {
      MAC::out() << " ... OK " << std::endl ;
      MAC::out().flush() ;
   }
}

//----------------------------------------------------------------------
LA_SeqVector const*
LA_DistVector:: local_vector( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: local_vector const" ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_SeqVector* result = LOCAL_VECTOR ;

   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->nb_rows()==row_distribution()->local_number() ) ;
   MAC_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqVector*
LA_DistVector:: local_vector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: local_vector" ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_SeqVector* result = LOCAL_VECTOR ;

   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->nb_rows()==row_distribution()->local_number() ) ;
   MAC_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: print_items( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: print_items" ) ;
   MAC_CHECK_PRE( print_items_PRE( os, indent_width ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( is_synchronized() ) ;

   std::string const space( indent_width, ' ' ) ;
   os << space << "nb_rows:" << nb_rows() << std::endl ;
   if( nb_rows()>0 )
   {
      std::ios_base::fmtflags original_flags = os.flags() ;
      os.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
      std::streamsize p = os.precision() ;
      os << std::setprecision( 7 ) ;
      for( size_t iRow = FIRST ; iRow<LAST ; ++iRow )
      {
         os << space << "Row n°" << iRow << "  " ;
         double const x = item(iRow) ;
         os << std::setw(15) << x << std::endl ;
      }
      os << std::setprecision(p) ;
      os.flags( original_flags ) ;
   }
}

//----------------------------------------------------------------------
void
LA_DistVector:: set_as_reciprocal( LA_Vector const* a,
                                   double smallest_inverted_item,
                                   double default_value )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: set_as_reciprocal LA_Vector" ) ;
   MAC_CHECK_PRE( set_as_reciprocal_PRE(
                         a, smallest_inverted_item, default_value ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const*>( a ) !=0 ) ;
   LA_DistVector const* da = static_cast<LA_DistVector const*>( a ) ;
   LOCAL_VECTOR->set_as_reciprocal(
           da->LOCAL_VECTOR, smallest_inverted_item, default_value ) ;
   set_synchronized() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_as_reciprocal_POST(a) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set_as_v_product( LA_Vector const* a, LA_Vector const* b )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: set_as_v_product LA_Vector" ) ;
   MAC_CHECK_PRE( set_as_v_product_PRE(a,b) ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const*>( a ) !=0 ) ;
   LA_DistVector const* da = static_cast<LA_DistVector const*>( a ) ;
   MAC_CHECK( dynamic_cast<LA_DistVector const*>( b ) !=0 ) ;
   LA_DistVector const* db = static_cast<LA_DistVector const*>( b ) ;

   LOCAL_VECTOR->set_as_v_product( da->LOCAL_VECTOR, db->LOCAL_VECTOR ) ;
   set_synchronized() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_as_v_product_POST( a, b ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: set_item( size_t i, double x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: set_item" ) ;
   MAC_CHECK_PRE( set_item_PRE( i ) ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;

   if( !only_local_modifs() )
   {
      set_unsynchronized_state( LA::NotSync_set ) ;
   }

   if( i<FIRST || i>=LAST )
   {
      size_t r = DIST->rank_of( i ) ;
      MAC_CHECK( r != MAC_Exec::communicator()->rank() ) ;
      intVector& idx = GLOBAL_INDEXES[r] ;
      doubleVector& val = GLOBAL_VALUES[r] ;
      bool found = false ;
      size_t n = idx.size() ;
      for( size_t j=0 ; j<n ; j++ )
      {
         if( idx(j)==(int)i )
         {
            val(j) = x ;
            found = true ;
            break ;
         }
      }
      if( !found )
      {
         idx.append( i ) ;
         val.append( x ) ;
         NB_EXTRA++ ;
      }
   }
   else
   {
      LOCAL_VECTOR->data()[i-FIRST] = x ;
   }
   MAC_CHECK_POST( set_item_POST( i, OLD( state ) ) ) ;

}

//----------------------------------------------------------------------
void
LA_DistVector:: add_to_item( size_t i, double x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: add_to_item" ) ;
   MAC_CHECK_PRE( add_to_item_PRE( i ) ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;

   if( !only_local_modifs() )
   {
      set_unsynchronized_state( LA::NotSync_add ) ;
   }

   if( x != 0.0 )
   {
      if( i<FIRST || i>=LAST )
      {
         size_t r = DIST->rank_of( i ) ;
         MAC_CHECK( r != MAC_Exec::communicator()->rank() ) ;
         intVector& idx = GLOBAL_INDEXES[r] ;
         doubleVector& val = GLOBAL_VALUES[r] ;
         bool found = false ;
         size_t n = idx.size() ;
         for( size_t j=0 ; j<n ; j++ )
         {
            if( idx(j)==(int)i )
            {
               val(j) += x ;
               found = true ;
               break ;
            }
         }
         if( !found )
         {
            idx.append( i ) ;
            val.append( x ) ;
            NB_EXTRA++ ;
         }
      }
      else
      {
         LOCAL_VECTOR->data()[i-FIRST] += x ;
      }
   }

   MAC_CHECK_POST( add_to_item_POST( i, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistVector:: write( std::string const& filename ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistVector:: write" ) ;
   MAC_CHECK_PRE( is_synchronized() ) ;

   size_t rank = MAC_Exec::communicator()->rank() ;
   size_t size = MAC_Exec::communicator()->nb_ranks() ;
   int dummy ;

   if( rank>0 ) MAC_Exec::communicator()->receive( rank-1, dummy ) ;

   std::ofstream out( filename.c_str(),
                      ( rank==0 ? std::ios::out
                                : std::ios::out|std::ios::app ) ) ;
   if( rank==0 ) out << nb_rows() ;

   for( size_t i=0 ; i<LOCAL_VECTOR->nb_rows() ; i++ )
   {
      out << std::endl << LOCAL_VECTOR->item(i) ;
   }

   out.close() ;

   if( rank!=size-1 )
   {
      MAC_Exec::communicator()->send( rank+1, dummy ) ;
      MAC_Exec::communicator()->receive( size-1, dummy ) ;
   }
   else
   {
      for(size_t i=0 ; i<size-1 ; i++ )
         MAC_Exec::communicator()->send( i, dummy ) ;
   }
}

//----------------------------------------------------------------------
bool
LA_DistVector:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Vector::invariant() ) ;
   if( INITIALIZED )
   {
      MAC_ASSERT( distribution_strategy() == LA::FromLocalSize ||
                  distribution_strategy() == LA::FromGlobalSize ) ;
      MAC_ASSERT( MAC_Exec::communicator()!=0 ) ;
      MAC_ASSERT( DIST!=0 ) ;
      MAC_ASSERT( LOCAL_VECTOR!=0 ) ;
      MAC_ASSERT( LOCAL_VECTOR->distribution_strategy() == 
                  LA::NoDistribution ) ;
      MAC_ASSERT( GLOBAL_VALUES.size() ==
                             MAC_Exec::communicator()->nb_ranks() ) ;
      MAC_ASSERT( GLOBAL_INDEXES.size() ==
                             MAC_Exec::communicator()->nb_ranks() ) ;
      MAC_ASSERT( DIST->local_number()==LOCAL_VECTOR->nb_rows() ) ;
      MAC_ASSERT( ( state() != LA::Sync ) ||
                  FORALL( (size_t i=0 ; i<GLOBAL_VALUES.size() ; i++),
                          GLOBAL_VALUES[i].size()==0 ) ) ;
      MAC_ASSERT( FORALL( (size_t i=0 ; i<GLOBAL_VALUES.size() ; i++),
                          GLOBAL_VALUES[i].size()==GLOBAL_INDEXES[i].size() ) ) ;
      MAC_ASSERT( FIRST == DIST->first_local_index() ) ;
      MAC_ASSERT( LAST == DIST->local_index_limit() ) ;
   }
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_DistVector:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Vector::implementation_POST( result ) ) ;
   MAC_ASSERT( result == LA_DistImplementation::object() ) ;
   return( true ) ;
}
