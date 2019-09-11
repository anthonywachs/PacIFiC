#include <LA_DistMatrix.hh>

#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Timer.hh>
#include <MAC_Vector.hh>
#include <MAC_assertions.hh>
#include <boolVector.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>

#include <LA_BlockSeqMatrix.hh>
#include <LA_CRSmatrix.hh>
#include <LA_PelMatrix.hh>
#include <LA_SeqMatrix.hh>
#include <LA_MatrixIterator.hh>
#include <LA_DistImplementation.hh>
#include <LA_DistVector.hh>
#include <LA_DistScatter.hh>
#include <LA_PairOfMatrixIterator.hh>
#include <LA_ShiftedIndexMatrixIterator.hh>

#include <iostream>

LA_DistMatrix const* LA_DistMatrix:: PROTOTYPE = new LA_DistMatrix() ;

//----------------------------------------------------------------------
LA_DistMatrix:: LA_DistMatrix( void )
//----------------------------------------------------------------------
   : LA_Matrix( "LA_DistMatrix" )
   , VERBOSE( false )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , ROW_DIST( 0 )
   , COL_DIST( 0 )
   , FIRST_LOCAL_ROW( MAC::bad_index() )
   , LAST_LOCAL_ROW( MAC::bad_index() )
   , LOCAL_DIAG_MATRIX( 0 )
   , LOCAL_NODIAG_MATRIX( 0 )
   , NON_LOCAL_MATRIX( 0 )
   , GLOBAL_VEC( 0 )
   , SCATTER_OK( false )
   , SCATTER( 0 )
   , SEQ_PROTO( 0 )
   , INITIAL_PROTO( 0 )
   , FINAL_PROTO( 0 )
   , INITIALIZED( false )
   , IN_PLACE( false )
   , SQUARE( false )
{
   MAC_LABEL( "LA_DistMatrix:: LA_DistMatrix( prototype )" ) ;

   MAC_CHECK_POST( is_a_prototype() ) ;
   MAC_CHECK_POST( name() == "LA_DistMatrix" ) ;
   MAC_CHECK_POST( !is_symmetric() ) ;
   MAC_CHECK_POST( is_resizable() ) ;
   MAC_CHECK_POST( state() == LA::NotSync_undef ) ;
   MAC_CHECK_POST( nb_rows() == 0 ) ;
   MAC_CHECK_POST( nb_cols() == 0 ) ;
   MAC_CHECK_POST( nb_stored_items() == 0 ) ;
   MAC_CHECK_POST( distribution_strategy() == LA::InvalidDistribution ) ;
}

//----------------------------------------------------------------------
LA_DistMatrix:: LA_DistMatrix( MAC_Object* a_owner,
                               size_t a_nb_rows, size_t a_nb_cols,
                               size_t a_nb_local_rows, size_t a_nb_local_cols,
                               LA::DistributionStrategy dist_strat,
                               LA_SeqMatrix* initial_prototype,
                               LA_SeqMatrix* final_prototype,
                               bool verbose )
//----------------------------------------------------------------------
   : LA_Matrix( a_owner, "LA_DistMatrix", dist_strat )
   , VERBOSE( verbose )
   , NB_ROWS( MAC::bad_index() )
   , NB_COLS( MAC::bad_index() )
   , ROW_DIST( MAC_DistributedPartition::create( this ) )
   , COL_DIST( MAC_DistributedPartition::create( this ) )
   , LOCAL_DIAG_MATRIX( 0 )
   , LOCAL_NODIAG_MATRIX( 0 )
   , NON_LOCAL_MATRIX( 0 )
   , GLOBAL_VEC( LA_SeqVector::create( this, 0 ) )
   , SCATTER_OK( false )
   , SCATTER( 0 )
   , SEQ_PROTO( initial_prototype )
   , INITIAL_PROTO( initial_prototype )
   , FINAL_PROTO( final_prototype )
   , HAS_SUBST_PROTO( false )
   , INITIALIZED( false )
   , IN_PLACE( false )
   , SQUARE( false )
{
   MAC_LABEL( "LA_DistMatrix:: LA_DistMatrix" ) ;
   MAC_CHECK_PRE( initial_prototype != 0 ) ;
   MAC_CHECK_PRE( initial_prototype->owner() == 0 ) ;
   MAC_CHECK_PRE( !initial_prototype->is_desynchronizable() ) ;
   MAC_CHECK_PRE( IMPLIES( final_prototype != 0,
                           final_prototype->owner() == 0 ) ) ;
   MAC_CHECK_PRE( IMPLIES( final_prototype != 0,
                           !final_prototype->is_desynchronizable() ) ) ;
   MAC_CHECK_PRE( dist_strat == LA::FromGlobalSize ||
                  dist_strat == LA::FromLocalSize ) ;

   initial_prototype->set_owner( this ) ;
   if( final_prototype != 0 ) final_prototype->set_owner( this ) ;

   re_initialize( a_nb_rows, a_nb_cols,
                  a_nb_local_rows, a_nb_local_cols ) ;
   INITIALIZED = true ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( distribution_strategy() == dist_strat ) ;
   MAC_CHECK_POST( !is_a_prototype() ) ;
   MAC_CHECK_POST( name() == "LA_DistMatrix" ) ;
   MAC_CHECK_POST( !is_symmetric() ) ;
   MAC_CHECK_POST( is_resizable() ) ;
   MAC_CHECK_POST( state() == LA::NotSync_undef ) ;
   MAC_CHECK_POST( initial_prototype->owner() == this ) ;
   MAC_CHECK_POST(
      IMPLIES( final_prototype != 0, final_prototype->owner() == this ) ) ;
}

//----------------------------------------------------------------------
LA_DistMatrix*
LA_DistMatrix:: create( MAC_Object* a_owner,
                        size_t a_nb_rows, size_t a_nb_cols,
                        size_t a_nb_local_rows, size_t a_nb_local_cols,
                        LA::DistributionStrategy dist_strat,
                        bool verbose )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: create" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( EQUIVALENT( a_nb_rows == 0, a_nb_cols == 0 ) ) ;
   MAC_CHECK_PRE( dist_strat == LA::FromGlobalSize ||
                  dist_strat == LA::FromLocalSize ) ;
   MAC_ASSERT( IMPLIES( dist_strat == LA::FromGlobalSize,
                        a_nb_rows != MAC::bad_index() ) ) ;
   MAC_ASSERT( IMPLIES( dist_strat == LA::FromGlobalSize,
                        a_nb_cols != MAC::bad_index() ) ) ;
   MAC_ASSERT(
     IMPLIES( dist_strat == LA::FromGlobalSize,
        MAC_Exec::communicator()->same_value_everywhere( (int) a_nb_rows ) ) ) ;
   MAC_ASSERT(
     IMPLIES( dist_strat == LA::FromGlobalSize,
        MAC_Exec::communicator()->same_value_everywhere( (int) a_nb_cols ) ) ) ;

   MAC_ASSERT( IMPLIES( dist_strat == LA::FromLocalSize,
                        a_nb_local_rows != MAC::bad_index() ) ) ;
   MAC_ASSERT( IMPLIES( dist_strat == LA::FromLocalSize,
                        a_nb_local_cols != MAC::bad_index() ) ) ;


   LA_SeqMatrix* initial_proto = LA_PelMatrix::create( 0, 0, 0 ) ;
   LA_DistMatrix* result = new LA_DistMatrix( a_owner,
                                              a_nb_rows, a_nb_cols,
                                              a_nb_local_rows, a_nb_local_cols,
                                              dist_strat,
                                              initial_proto, 0,
                                              verbose ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->name() == "LA_DistMatrix" ) ;
   MAC_CHECK_POST( result->is_resizable() ) ;
   MAC_CHECK_POST( !result->is_symmetric() ) ;
   MAC_CHECK_POST( result->state() == LA::NotSync_undef ) ;
   MAC_ASSERT( IMPLIES( result->is_desynchronizable(),
                        ! result->is_synchronized() ) ) ;
   MAC_CHECK_PRE( result->distribution_strategy() == dist_strat ) ;
   MAC_CHECK_POST(
             IMPLIES( result->distribution_strategy() == LA::FromGlobalSize,
                      result->nb_rows() == a_nb_rows ) ) ;
   MAC_CHECK_POST(
             IMPLIES( result->distribution_strategy() == LA::FromGlobalSize,
                      result->nb_cols() == a_nb_cols ) ) ;
   MAC_ASSERT(
          IMPLIES( result->distribution_strategy() == LA::FromLocalSize,
          result->row_distribution()->local_number() ==  a_nb_local_rows ) ) ;
   MAC_ASSERT(
          IMPLIES( result->distribution_strategy() == LA::FromLocalSize,
          result->col_distribution()->local_number() ==  a_nb_local_cols ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DistMatrix*
LA_DistMatrix:: create_matrix( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: create_matrix" ) ;
   MAC_CHECK_PRE( create_matrix_PRE( a_owner ) ) ;

   LA_DistMatrix* result =
      new LA_DistMatrix( a_owner,
                         ROW_DIST->global_number(), COL_DIST->global_number(),
                         ROW_DIST->local_number(), COL_DIST->local_number(),
                         distribution_strategy(),
                         INITIAL_PROTO->create_matrix(0),
                         ( FINAL_PROTO!=0 ? FINAL_PROTO->create_matrix(0) : 0 ),
                         VERBOSE ) ;

   MAC_CHECK_POST( create_matrix_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DistMatrix*
LA_DistMatrix:: create_replica( MAC_Object* a_owner,
                                MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   size_t a_nb_rows = 0 ;
   size_t a_nb_cols = 0 ;
   size_t a_nb_local_rows = 0 ;
   size_t a_nb_local_cols = 0 ;
   LA::DistributionStrategy dist_strat = LA::FromGlobalSize ;

   if( exp->has_module( "DistributionStrategy" ) )
   {
      MAC_ModuleExplorer* se =
                exp->create_subexplorer( 0, "DistributionStrategy" ) ;
      std::string const& type = se->string_data( "type" ) ;
      if( type == "from_global_sizes" )
      {
         dist_strat = LA::FromGlobalSize ;
         if( se->has_entry( "nb_rows" ) )
         {
            a_nb_rows = exp->int_data( "nb_rows" ) ;
            a_nb_cols = exp->int_data( "nb_cols" ) ;
            if( a_nb_rows <= 0 )
            {
               MAC_Error::object()->raise_bad_data_value( exp, "nb_rows",
                                    "A strictly positive value is expected" ) ;
            }
            if( a_nb_cols <= 0 )
            {
               MAC_Error::object()->raise_bad_data_value( exp, "nb_cols",
                                    "A strictly positive value is expected" ) ;
            }
         }
      }
      else if( type == "from_local_sizes" )
      {
         dist_strat = LA::FromLocalSize ;
         if( se->has_entry( "nb_local_rows" ) )
         {
            intVector nb_rows_ll = exp->intVector_data( "nb_local_rows" ) ;
            intVector nb_cols_ll = exp->intVector_data( "nb_local_cols" ) ;
            if( nb_rows_ll.size() != communicator()->nb_ranks() )
            {
               MAC_Error::object()->raise_bad_data_value( exp, "nb_local_rows",
                                    "A vector of size nb_ranks() is expected" ) ;
            }
            if( nb_cols_ll.size() != communicator()->nb_ranks() )
            {
               MAC_Error::object()->raise_bad_data_value( exp, "nb_local_cols",
                                    "A vector of size nb_ranks() is expected" ) ;
            }
            a_nb_local_rows = nb_rows_ll( communicator()->rank() ) ;
            a_nb_local_cols = nb_cols_ll( communicator()->rank() ) ;
            if( a_nb_local_rows <= 0 )
            {
               MAC_Error::object()->raise_bad_data_value( exp, "nb_local_rows",
                                    "A strictly positive value is expected" ) ;
            }
            if( a_nb_local_cols <= 0 )
            {
               MAC_Error::object()->raise_bad_data_value( exp, "nb_local_cols",
                                    "A strictly positive value is expected" ) ;
            }
         }
      }
      else
      {
         MAC_Error::object()->raise_bad_data_value( se, "type",
                                             "\"from_global_sizes\" or "
                                             "\"from_local_sizes\"" ) ;
      }
      se->destroy() ; se = 0 ;
   }

   LA_SeqMatrix* initial_proto = 0 ;
   if( exp->has_module( "initial_block_prototype" ) )
   {
      MAC_ModuleExplorer const* sexp =
         exp->create_subexplorer( 0, "initial_block_prototype" ) ;
      initial_proto = LA_SeqMatrix::make( 0, sexp ) ;
      sexp->destroy() ;
   }
   else
   {
      initial_proto = LA_PelMatrix::create( 0, 0, 0 ) ;
   }

   LA_SeqMatrix* final_proto = 0 ;
   if( exp->has_module( "final_block_prototype" ) )
   {
      MAC_ModuleExplorer const* sexp =
         exp->create_subexplorer( 0, "final_block_prototype" ) ;
      final_proto = LA_SeqMatrix::make( 0, sexp ) ;
      sexp->destroy() ;
   }

   bool verbose = false ;
   if( exp->has_entry( "verbose" ) )
   {
      verbose = exp->bool_data( "verbose" ) ;
      exp->set_default( "verbose", "false" ) ;
   }

   LA_DistMatrix* result =
      new LA_DistMatrix( a_owner,
                         a_nb_rows, a_nb_cols,
                         a_nb_local_rows, a_nb_local_cols,
                         dist_strat,
                         initial_proto, final_proto,
                         verbose ) ;

   MAC_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   MAC_CHECK_POST( ( result->distribution_strategy() == LA::FromGlobalSize ) ||
                   ( result->distribution_strategy() == LA::FromLocalSize ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_DistMatrix:: ~LA_DistMatrix( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: re_initialize( size_t a_nb_rows, size_t a_nb_cols,
                               size_t a_nb_local_rows, size_t a_nb_local_cols )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: re_initialize" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( re_initialize_PRE( a_nb_rows, a_nb_cols,
                                     a_nb_local_rows, a_nb_local_cols ) ) ;
   MAC_CHECK_PRE( EQUIVALENT( a_nb_rows == 0, a_nb_cols == 0 ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   SEQ_PROTO = INITIAL_PROTO ;
   HAS_SUBST_PROTO = false ;

   if( nb_rows() != MAC::bad_index()  )
   {
      destroy_possession( NON_LOCAL_MATRIX ) ; NON_LOCAL_MATRIX = 0 ;
      if( LOCAL_DIAG_MATRIX!=0 )
      {
         destroy_possession( LOCAL_DIAG_MATRIX ) ; LOCAL_DIAG_MATRIX=0 ;
      }
      if( LOCAL_NODIAG_MATRIX!=0 )
      {
         destroy_possession( LOCAL_NODIAG_MATRIX ) ; LOCAL_NODIAG_MATRIX=0 ;
      }
   }

   if( distribution_strategy() == LA::FromLocalSize )
   {
      ROW_DIST->set_local_number( a_nb_local_rows ) ;
      COL_DIST->set_local_number( a_nb_local_cols ) ;
   }
   else if( distribution_strategy() == LA::FromGlobalSize )
   {
      ROW_DIST->distribute_global_number( a_nb_rows ) ;
      COL_DIST->distribute_global_number( a_nb_cols ) ;
   }
   FIRST_LOCAL_ROW = ROW_DIST->first_local_index() ;
   LAST_LOCAL_ROW  = ROW_DIST->local_index_limit() ;

   size_t_vector global_col_partition( 1 ) ;
   global_col_partition( 0 ) = COL_DIST->global_number() ;
   NON_LOCAL_MATRIX = LA_BlockSeqMatrix::create( this,
                                                 ROW_DIST->partitioning(),
                                                 global_col_partition,
                                                 SEQ_PROTO ) ;

   size_t loc = ROW_DIST->local_number() ;
   NB_COLS = COL_DIST->global_number() ;
   NB_ROWS = ROW_DIST->global_number() ;
   SQUARE = ( NB_COLS==NB_ROWS ) ;

   if( SQUARE )
   {
      LOCAL_DIAG_MATRIX = SEQ_PROTO->create_matrix( this ) ;
      LOCAL_DIAG_MATRIX->re_initialize( loc, loc ) ;
   }

   LOCAL_NODIAG_MATRIX = SEQ_PROTO->create_matrix( this ) ;
   LOCAL_NODIAG_MATRIX->re_initialize( loc, NB_COLS  ) ;

   GLOBAL_VEC->re_initialize( NB_COLS ) ;

   SCATTER_OK = false ;
   if( SCATTER == 0 )
   {
      SCATTER = LA_DistScatter::create( this, COL_DIST ) ;
   }
   else
   {
      SCATTER->re_initialize( COL_DIST ) ;
   }

   set_unsynchronized_state( LA::NotSync_undef ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( re_initialize_POST( a_nb_rows, a_nb_cols,
                                       a_nb_local_rows, a_nb_local_cols ) ) ;
}

//----------------------------------------------------------------------
LA_DistVector*
LA_DistMatrix:: create_vector( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: create_vector" ) ;
   MAC_CHECK_PRE( create_vector_PRE( a_owner ) ) ;

   LA_DistVector* result =
      LA_DistVector::create( a_owner,
                             row_distribution()->global_number(),
                             row_distribution()->local_number(),
                             distribution_strategy() ) ;

   MAC_CHECK_POST( create_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_MatrixIterator*
LA_DistMatrix:: create_stored_item_iterator( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: create_stored_item_iterator" ) ;
   MAC_CHECK_PRE( create_stored_item_iterator_PRE( a_owner ) ) ;

   LA_MatrixIterator* result = 0 ;

   LA_MatrixIterator* first_internal =
      LOCAL_NODIAG_MATRIX->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* first_shift =
      LA_ShiftedIndexMatrixIterator::create(
         0, FIRST_LOCAL_ROW, 0, this, first_internal ) ;
   if( SQUARE )
   {
      LA_MatrixIterator* second_internal =
         LOCAL_DIAG_MATRIX->create_stored_item_iterator( 0 ) ;
      LA_MatrixIterator* second_shift =
         LA_ShiftedIndexMatrixIterator::create(
            0, FIRST_LOCAL_ROW, FIRST_LOCAL_ROW, this, second_internal ) ;
      result = LA_PairOfMatrixIterator::create(
                          a_owner, this, first_shift, second_shift ) ;
   }
   else
   {
      result = first_shift ;
      if( a_owner != 0 ) first_shift->set_owner( a_owner ) ;
   }

   MAC_CHECK_POST( create_stored_item_iterator_POST( a_owner, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix*
LA_DistMatrix:: create_local_matrix( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: create_local_matrix" ) ;
   MAC_CHECK_PRE( create_local_matrix_PRE( a_owner ) ) ;

   LA_PelMatrix* result = LA_PelMatrix::create( a_owner,
                                                nb_rows(), nb_cols() ) ;

   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
   for( it->start_all_items() ; it->is_valid() ; it->go_next() )
   {
      result->set_item( it->row() , it->col(), it->item() ) ;
   }
   it->destroy() ; it = 0 ;

   result->synchronize() ;

   MAC_CHECK_POST( create_local_matrix_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Implementation const*
LA_DistMatrix:: implementation( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: implementation" ) ;

   LA_Implementation const* result = LA_DistImplementation::object() ;

   MAC_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_DistMatrix:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: is_desynchronizable" ) ;

   bool result = true ;

   MAC_CHECK_POST( result == true ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Communicator const*
LA_DistMatrix:: communicator( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: communicator" ) ;
   static MAC_Communicator const* result = MAC_Exec::communicator() ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_DistMatrix:: nb_stored_items( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: nb_stored_items" ) ;
   MAC_CHECK_PRE( nb_stored_items_PRE() ) ;

   size_t result = LOCAL_NODIAG_MATRIX->nb_stored_items() ;
   if( SQUARE ) result += LOCAL_DIAG_MATRIX->nb_stored_items() ;

   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_DistMatrix:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: nb_rows" ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( NB_ROWS ) ;
}

//----------------------------------------------------------------------
size_t
LA_DistMatrix:: nb_cols( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: nb_cols" ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( NB_COLS ) ;
}

//----------------------------------------------------------------------
double
LA_DistMatrix:: item( size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: item" ) ;
   MAC_CHECK_PRE( i < LAST_LOCAL_ROW && i>=FIRST_LOCAL_ROW ) ;
   MAC_CHECK_PRE( j < nb_cols() ) ;
   MAC_CHECK_PRE( state() == LA::Sync ) ;
   MAC_CHECK_INV( invariant() ) ;

   return ( !SQUARE || j < FIRST_LOCAL_ROW || j >= LAST_LOCAL_ROW ?
            LOCAL_NODIAG_MATRIX->item( i-FIRST_LOCAL_ROW, j ) :
            LOCAL_DIAG_MATRIX->item( i-FIRST_LOCAL_ROW, j-FIRST_LOCAL_ROW ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: extract_diag( LA_Vector* diag ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: extract_diag" ) ;
   MAC_CHECK_PRE( extract_diag_PRE(diag) ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector* >( diag ) != 0 ) ;
   LA_DistVector* dvec = static_cast<LA_DistVector* >( diag ) ;

   diagonal_block_matrix()->extract_diag( dvec->local_vector() ) ;
   dvec->synchronize() ;

   MAC_CHECK_POST( extract_diag_POST( diag )  ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix*
LA_DistMatrix:: diagonal_block_matrix( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: diagonal_block_matrix" ) ;
   MAC_CHECK_PRE( nb_rows() > 0 ) ;
   MAC_CHECK_PRE( nb_rows() == nb_cols() ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_SeqMatrix* result = LOCAL_DIAG_MATRIX ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == this ) ;
   MAC_CHECK_POST( result->nb_rows() == result->nb_cols() ) ;
   MAC_CHECK_POST( result->nb_rows() == row_distribution()->local_number() ) ;
   MAC_CHECK_POST( ! result->is_desynchronizable() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_SeqMatrix const*
LA_DistMatrix:: block_prototype( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: block_prototype" ) ;

   LA_SeqMatrix const* result = SEQ_PROTO ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->nb_rows() == 0 ) ;
   MAC_CHECK_POST( result->nb_cols() == 0 ) ;
   MAC_CHECK_POST( result->state() == LA::Sync ) ;
   MAC_CHECK_POST( result->is_resizable() ) ;
   MAC_CHECK_POST( ! result->is_desynchronizable() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: set_block_prototype( LA_SeqMatrix const* a_proto )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: set_block_prototype" ) ;
   MAC_CHECK_PRE( a_proto != 0 ) ;
   MAC_CHECK_PRE( a_proto->owner() == this ) ;
   MAC_CHECK_PRE( a_proto->nb_rows() == 0 ) ;
   MAC_CHECK_PRE( a_proto->nb_cols() == 0 ) ;
   MAC_CHECK_PRE( ! a_proto->is_desynchronizable() ) ;
   MAC_CHECK_PRE( a_proto->is_resizable() ) ;

   if( VERBOSE )
   {
      MAC::out() << "Setting block prototype to "
                 << a_proto->name() << std::endl ;
   }

   SEQ_PROTO = a_proto ;

   if( LOCAL_DIAG_MATRIX!=0 )
   {
      LA_SeqMatrix* old = LOCAL_DIAG_MATRIX ;
      LOCAL_DIAG_MATRIX = SEQ_PROTO->create_copy( this, LOCAL_DIAG_MATRIX ) ;
      destroy_possession( old ) ;
   }

   LA_SeqMatrix* old = LOCAL_NODIAG_MATRIX ;
   LOCAL_NODIAG_MATRIX = SEQ_PROTO->create_copy( this, LOCAL_NODIAG_MATRIX ) ;
   destroy_possession( old ) ;

   NON_LOCAL_MATRIX->set_block_prototype(
                              SEQ_PROTO->create_matrix( NON_LOCAL_MATRIX ) ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( block_prototype() == a_proto ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: stop_local_modifs( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: stop_local_modifs" ) ;
   MAC_CHECK_PRE( stop_local_modifs_PRE() ) ;

   set_only_local_modifs_state( false ) ;
   SCATTER_OK = false ;

   MAC_CHECK_POST( stop_local_modifs_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: nullify( void  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: nullify" ) ;
   MAC_CHECK_PRE( nullify_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( LOCAL_DIAG_MATRIX!=0 ) LOCAL_DIAG_MATRIX->nullify() ;
   LOCAL_NODIAG_MATRIX->nullify() ;

   NON_LOCAL_MATRIX->set_stored_items( MAC::bad_double() ) ;
   set_unsynchronized_state( LA::NotSync_undef ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( nullify_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: scale( double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: scale" ) ;
   MAC_CHECK_PRE( scale_PRE( alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( alpha != 1.0 )
   {
      if( LOCAL_DIAG_MATRIX!=0 ) LOCAL_DIAG_MATRIX->scale( alpha ) ;
      LOCAL_NODIAG_MATRIX->scale( alpha)  ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: set_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: set_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;
   MAC_CHECK_PRE( set_item_PRE( i, j ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( !only_local_modifs() )
   {
      set_unsynchronized_state( LA::NotSync_set ) ;
   }

   if( i<LAST_LOCAL_ROW && i>=FIRST_LOCAL_ROW )
   {
      if( SQUARE && j<LAST_LOCAL_ROW && j>=FIRST_LOCAL_ROW )
      {
         LOCAL_DIAG_MATRIX->set_item( i-FIRST_LOCAL_ROW,j-FIRST_LOCAL_ROW,x ) ;
      }
      else
      {
         LOCAL_NODIAG_MATRIX->set_item( i-FIRST_LOCAL_ROW,j,x ) ;
      }

   }
   else
   {
      NON_LOCAL_MATRIX->set_item( i,j,x ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: add_to_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: add_to_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;
   MAC_CHECK_PRE( add_to_item_PRE( i, j ) ) ;
   MAC_CHECK_INV( invariant() ) ;

  if( !only_local_modifs() )
  {
     set_unsynchronized_state( LA::NotSync_add ) ;
  }

   if( i<LAST_LOCAL_ROW && i>=FIRST_LOCAL_ROW )
   {
      if( SQUARE && j<LAST_LOCAL_ROW && j>=FIRST_LOCAL_ROW )
      {
         LOCAL_DIAG_MATRIX->add_to_item(
                    i-FIRST_LOCAL_ROW, j-FIRST_LOCAL_ROW, x ) ;
      }
      else
      {
         LOCAL_NODIAG_MATRIX->add_to_item( i-FIRST_LOCAL_ROW, j, x ) ;
      }
   }
   else
   {
      if( NON_LOCAL_MATRIX->item( i,j) == MAC::bad_double() )
      {
         NON_LOCAL_MATRIX->nullify_item( i, j ) ;
      }
      NON_LOCAL_MATRIX->add_to_item( i, j, x ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_to_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: set( LA_Matrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: set" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( set_PRE( A, same_pattern ) ) ;

   MAC_CHECK( dynamic_cast<LA_DistMatrix const* >( A ) != 0 ) ;
   LA_DistMatrix const* dA = static_cast<LA_DistMatrix const* >( A ) ;

   if( LOCAL_DIAG_MATRIX!=0 )
   {
      LOCAL_DIAG_MATRIX->set( dA->LOCAL_DIAG_MATRIX, same_pattern ) ;
   }
   LOCAL_NODIAG_MATRIX->set( dA->LOCAL_NODIAG_MATRIX, same_pattern ) ;

   NON_LOCAL_MATRIX->set_stored_items( MAC::bad_double() ) ;
   if( !same_pattern )
   {
      SCATTER_OK = false ;
   }
   set_unsynchronized_state( LA::NotSync_undef ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_POST( A ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: multiply_vec_then_add( LA_Vector const* x,
                                       LA_Vector* y,
                                       double alpha, double beta ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: multiply_vec_then_add" ) ;
   MAC_CHECK_PRE( multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const* >( x ) != 0 ) ;
   LA_DistVector const* dx = static_cast<LA_DistVector const* >( x ) ;
   MAC_CHECK( dynamic_cast<LA_DistVector* >( y ) != 0 ) ;
   LA_DistVector* dy = static_cast<LA_DistVector* >( y ) ;

   multiply_vec_then_add_IMP( dx, dy, alpha, beta ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: multiply_vec_then_add_IMP(
                               LA_DistVector const* x, LA_DistVector* y,
                               double alpha, double beta ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: multiply_vec_then_add_IMP" ) ;
   MAC_CHECK_PRE( multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( !SCATTER_OK ) build_scatter() ;

   if( SQUARE )
   {
      SCATTER->get_begin( x ) ;

      LOCAL_DIAG_MATRIX->multiply_vec_then_add(
         x->local_vector(), y->local_vector(), alpha, beta ) ;

      SCATTER->get_end( const_cast<LA_SeqVector*>( GLOBAL_VEC ) ) ;

      LOCAL_NODIAG_MATRIX->multiply_vec_then_add(
                 GLOBAL_VEC, y->local_vector(), alpha, 1.0 ) ;
   }
   else
   {
      SCATTER->get( x, GLOBAL_VEC ) ;
      LOCAL_NODIAG_MATRIX->multiply_vec_then_add(
                     GLOBAL_VEC, y->local_vector(), alpha, beta ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: tr_multiply_vec_then_add( LA_Vector const* x, LA_Vector* y,
                                          double alpha, double beta ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: tr_multiply_vec_then_add" ) ;
   MAC_CHECK_PRE( tr_multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const* >( x ) != 0 ) ;
   LA_DistVector const* dx = static_cast<LA_DistVector const* >( x ) ;
   MAC_CHECK( dynamic_cast<LA_DistVector* >( y ) != 0 ) ;
   LA_DistVector* dy = static_cast<LA_DistVector* >( y ) ;

   tr_multiply_vec_then_add_IMP( dx, dy, alpha, beta ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( tr_multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: tr_multiply_vec_then_add_IMP(
                                LA_DistVector const* x, LA_DistVector* y,
                                double alpha, double beta ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: tr_multiply_vec_then_add_IMP" ) ;
   MAC_CHECK_PRE( tr_multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   y->scale( beta ) ;

   size_t N = nb_cols() ;

   LA_SeqVector* global_y = LA_SeqVector::create( 0, N ) ;
   LOCAL_NODIAG_MATRIX->tr_multiply_vec_then_add(
                             x->local_vector(), global_y, alpha, beta ) ;
   for( size_t i=0 ; i<N ; ++i )
   {
      y->add_to_item( i, global_y->item( i ) ) ;
   }
   global_y->destroy() ; global_y = 0 ;
   y->synchronize() ;

   if( LOCAL_DIAG_MATRIX!=0 )
   {
      LOCAL_DIAG_MATRIX->tr_multiply_vec_then_add(
         x->local_vector(), y->local_vector(), alpha, 1.0 ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( tr_multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: add_Mat( LA_Matrix const* A, double alpha,
                         bool same_pattern )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: add_Mat" ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;
   MAC_CHECK_PRE( add_Mat_PRE( A, alpha, same_pattern ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( dynamic_cast<LA_DistMatrix const* >( A ) != 0 ) ;
   LA_DistMatrix const* dA = static_cast<LA_DistMatrix const* >( A ) ;

   if( LOCAL_DIAG_MATRIX!=0 )
   {
      LOCAL_DIAG_MATRIX->add_Mat( dA->LOCAL_DIAG_MATRIX, alpha, same_pattern ) ;
   }
   LOCAL_NODIAG_MATRIX->add_Mat( dA->LOCAL_NODIAG_MATRIX, alpha, same_pattern ) ;

   if( !same_pattern )
   {
      SCATTER_OK = false ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                             double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: add_Mat_Mat" ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;
   MAC_CHECK_PRE( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( dynamic_cast<LA_DistMatrix const* >( A ) != 0 ) ;
   LA_DistMatrix const* dA = static_cast<LA_DistMatrix const* >( A ) ;
   MAC_CHECK( dynamic_cast<LA_DistMatrix const* >( B ) != 0 ) ;
   LA_DistMatrix const* dB = static_cast<LA_DistMatrix const* >( B ) ;

   LA_SeqMatrix* GLOBAL_MAT =
      LA_PelMatrix::create( 0, B->nb_rows(), B->nb_cols() ) ;
   if( !dA->SCATTER_OK ) dA->build_scatter() ;
   LA_DistScatter* a_scatter = dA->SCATTER ;
   size_t_vector const& receive_index = a_scatter->repatriated_items() ;
   if( dA->SQUARE )
   {
      a_scatter = LA_DistScatter::create( 0, B->row_distribution() ) ;
      size_t n = B->nb_rows() ;
      boolVector need( n ) ;
      need.set( false ) ;
      size_t N = receive_index.size() ;
      for( size_t i=0 ; i<N; i++ ) need( receive_index(i) ) = true ;
      size_t i1 = dA->row_distribution()->first_local_index() ;
      size_t i2 = dA->row_distribution()->local_index_limit() ;

      for( size_t i=i1 ; i<i2 ; i++ )
      {
         if( !need(i) ) N++ ;
         need(i) = true ;
      }

      size_t_vector rec( N ) ;
      size_t idx = 0 ;
      for( size_t i=0 ; i<n ; i++ )
         if( need(i) ) rec(idx++) = i ;
      MAC_CHECK( N == idx ) ;

      a_scatter->set_sorted( rec, rec) ;
   }

   a_scatter->get_rows( dB, GLOBAL_MAT ) ;

   LA_MatrixIterator* itA = dA->create_stored_item_iterator( 0 ) ;
   LA_MatrixIterator* itB = GLOBAL_MAT->create_stored_item_iterator( 0 ) ;
   start_local_modifs() ;
   for( itA->start_all_items() ; itA->is_valid() ; itA->go_next() )
   {
      for( itB->start_row_items(itA->col()) ; itB->is_valid() ; itB->go_next() )
      {
         add_to_item( itA->row(),
                      itB->col(),
                      alpha * itA->item() * itB->item() ) ;
      }
   }
   stop_local_modifs() ;//SCATTER_OF = false ;
   itA->destroy() ; itA = 0 ;
   itB->destroy() ; itB = 0 ;

   GLOBAL_MAT->destroy() ;
   if( a_scatter->owner()==0 ) a_scatter->destroy() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_Mat_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
size_t
LA_DistMatrix:: allocated_memory( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: allocated_memory" ) ;
   MAC_CHECK_INV( invariant() ) ;

   size_t result = NON_LOCAL_MATRIX->allocated_memory() ;
   if( LOCAL_DIAG_MATRIX!=0 )
   {
      result += LOCAL_DIAG_MATRIX->allocated_memory()
                            + LOCAL_NODIAG_MATRIX->allocated_memory() ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: synchronize( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: synchronize" ) ;
   MAC_CHECK_PRE( synchronize_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( communicator()->boolean_or( state() == LA::NotSync_add ) &&
       communicator()->boolean_or( state() == LA::NotSync_set ) )
   {
      MAC_Error::object()->raise_plain(
         "*** LA_DistMatrix synchronization error:\n"
         "    unable to synchronize matrices with\n"
         "    NotSync_add AND NotSync_set states." ) ;
   }

   if( !communicator()->boolean_and( state() == LA::NotSync_undef ) )
   {
      size_t rank = communicator()->rank() ;
      size_t nb_ranks = communicator()->nb_ranks() ;

      LA::SyncState effective_mode =
         ( communicator()->boolean_and( state()!=LA::NotSync_set ) ?
               LA::NotSync_add : LA::NotSync_set ) ;

      if( VERBOSE )
      {
         MAC::out() << "LA_DistMatrix:: synchronize" << std::endl ;
      }

      // Tout d'abord, on va commencer par s'envoyer le nombre d'elements
      //  a echanger entre process
      intVector exported_items( nb_ranks ) ;
      IN_PLACE = true ;

      for( size_t i=0 ; i<nb_ranks ; i++ )
      {
         if( i != rank && ROW_DIST->partitioning()(i) > 0 )
         {
            int N = 0 ;
            if( NON_LOCAL_MATRIX->has_submatrix( i, 0 ) )
            {
               LA_SeqMatrix const* sub_matrix =
                                       NON_LOCAL_MATRIX->submatrix( i, 0 ) ;
               MAC_CHECK( NON_LOCAL_MATRIX->state()==LA::Sync ) ;
               MAC_CHECK( sub_matrix->state()==LA::Sync ) ;
               N = sub_matrix->nb_stored_items() ;
               IN_PLACE &=
                  dynamic_cast<LA_CRSmatrix const*>( sub_matrix ) != 0 ;
            }
            exported_items( i ) = N ;
         }
      }

      intVector exchanged_items( nb_ranks ) ;
      communicator()->all_to_all( exported_items, exchanged_items ) ;

      if( VERBOSE )
      {
         MAC::out() << "Echanged item number to synchronize matrix : "
                    << std::endl ;
         MAC::out() << " Exported : " << exported_items << std::endl ;
         MAC::out() << " Imported : " << exchanged_items << std::endl ;
      }

      // Maintenant, on va les envoyer et recevoir
      if( IN_PLACE )
      {
         for( size_t i=0 ; i<nb_ranks ; i++ )
         {
            send_submatrix(exported_items( i ),i) ;
         }
         for( size_t i=0 ; i<nb_ranks ; i++ )
         {
            receive_submatrix( exchanged_items( i ), i, effective_mode ) ;
         }
      }
      else
      {
         for( size_t i=0 ; i<rank ; i++ )
         {
            receive_submatrix( exchanged_items( i ), i, effective_mode ) ;
         }

         for( size_t i=0 ; i<nb_ranks ; i++ )
         {
            send_submatrix(exported_items( i ), i ) ;
         }

         for( size_t i=rank+1 ; i<nb_ranks ; i++ )
         {
            receive_submatrix( exchanged_items( i ), i, effective_mode ) ;
         }
      }

      if( IN_PLACE )
      {
         for( size_t i=0 ; i<nb_ranks ; i++ )
         {
            wait_send_submatrix(exported_items( i ), i ) ;
         }
      }

      if( FINAL_PROTO!=0 && !HAS_SUBST_PROTO )
      {
         set_block_prototype( FINAL_PROTO ) ;
         HAS_SUBST_PROTO = true ;
      }
      NON_LOCAL_MATRIX->set_stored_items( MAC::bad_double() ) ;
      SCATTER_OK = false ;
   }
   LA_Matrix::synchronize() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( synchronize_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: send_submatrix( size_t N, size_t i ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: send_submatrix" ) ;

   if( N > 0 )
   {
      MAC_Timer* timer = 0 ;

      if( VERBOSE )
      {
         MAC::out() << "LA_DistMatrix:: send "<< N <<" items to " << i ;
         MAC::out().flush() ;
         timer = MAC_Timer::create( 0 ) ;
         timer->start() ;
      }

      MAC_CHECK( NON_LOCAL_MATRIX->has_submatrix(i,0) ) ;
      LA_SeqMatrix const* sub_matrix = NON_LOCAL_MATRIX->submatrix(i,0) ;
      MAC_ASSERT( sub_matrix->state() == LA::Sync ) ;

      LA_CRSmatrix const* crs_matrix = convert_if_needed( sub_matrix ) ;
      // MAC_ASSERT( IN_PLACE == (crs_matrix==sub_matrix) ) ;

      crs_matrix->send( communicator(), i, IN_PLACE ) ;

      if( crs_matrix!=sub_matrix ) crs_matrix->destroy() ;
      if( VERBOSE )
      {
         timer->stop() ;
         MAC::out() << " -> done in " ;
         timer->print( MAC::out(), 0 ) ;
         MAC::out() << std::endl ;
         MAC::out().flush() ;
         timer->destroy() ;
         timer = 0 ;
      }
   }
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: wait_send_submatrix( size_t N, size_t i ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: wait_send_submatrix" ) ;

   if( N > 0 )
   {
      if( VERBOSE )
      {
         MAC::out() << "LA_DistMatrix:: wait_send_submatrix to " << i << std::endl ;
      }
      MAC_CHECK( NON_LOCAL_MATRIX->has_submatrix(i,0) ) ;
      LA_SeqMatrix const* sub_matrix = NON_LOCAL_MATRIX->submatrix(i,0) ;
      LA_CRSmatrix const* crs_matrix = convert_if_needed( sub_matrix ) ;
      MAC_ASSERT( crs_matrix==sub_matrix ) ;

      crs_matrix->wait_send( communicator() ) ;
   }
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: receive_submatrix( size_t N,
                                   size_t i,
                                   LA::SyncState effective_mode )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: receive_submatrix" ) ;
   MAC_CHECK( effective_mode == LA::NotSync_add ||
              effective_mode == LA::NotSync_set ) ;
   if( N > 0 )
   {
      MAC_Timer* timer = 0 ;
      if( VERBOSE )
      {
         MAC::out() << "LA_DistMatrix:: receive "<< N <<" items from " << i  ;
         MAC::out().flush() ;
         timer = MAC_Timer::create( 0 ) ;
         timer->start() ;
      }
      LA_CRSmatrix* crs_matrix = LA_CRSmatrix::receive( 0, communicator(), i ) ;
      MAC_CHECK( crs_matrix->nb_stored_items()==N ) ;

      LA_MatrixIterator* it = crs_matrix->create_stored_item_iterator( crs_matrix ) ;

      if( effective_mode == LA::NotSync_add )
      {
         for( it->start_all_items() ; it->is_valid() ; it->go_next() )
         {
            double xx = it->item() ;
            if( xx!=MAC::bad_double() )
            {
               MAC_CHECK( it->row()<LAST_LOCAL_ROW-FIRST_LOCAL_ROW ) ;
               add_to_item( it->row()+FIRST_LOCAL_ROW, it->col(), xx ) ;
            }
         }
      }
      else
      {
         for( it->start_all_items() ; it->is_valid() ; it->go_next() )
         {
            double xx = it->item() ;
            if( xx!=MAC::bad_double() )
            {
               MAC_CHECK( it->row()<LAST_LOCAL_ROW-FIRST_LOCAL_ROW ) ;
               set_item( it->row()+FIRST_LOCAL_ROW, it->col(), xx ) ;
            }
         }
      }
      crs_matrix->destroy() ;
      if( VERBOSE )
      {
         timer->stop() ;
         MAC::out() << " -> done in " ;
         timer->print( MAC::out(), 0 ) ;
         MAC::out() << std::endl ;
         MAC::out().flush() ;
         timer->destroy() ;
         timer = 0 ;
      }
   }
}


//----------------------------------------------------------------------
void
LA_DistMatrix:: add_tMat( LA_Matrix const* A, double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: add_tMat" ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;
   MAC_CHECK_PRE( add_tMat_PRE( A, alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( dynamic_cast<LA_DistMatrix const* >( A ) != 0 ) ;
   LA_DistMatrix const* dA = static_cast<LA_DistMatrix const* >( A ) ;

   if(VERBOSE)
      MAC::out() << "LA_DistMatrix:: add_tMat" << std::endl ;

   MAC_DistributedPartition const* dist = ROW_DIST ;
   bool has_A_items_localy = dist->local_number() > 0 ;
   size_t_vector const& start = dist->start_of_partition() ;
   size_t_vector const& part = dist->partitioning() ;

   if( !dA->SCATTER_OK ) dA->build_scatter() ;
   LA_DistScatter* a_scatter = dA->SCATTER ;
   size_t_vector const& A_col_indexes = a_scatter->repatriated_items() ;

   size_t nb_ranks = communicator()->nb_ranks() ;
   size_t rank = communicator()->rank() ;
   size_t FIRST_LOCAL_ROW_A = dA->FIRST_LOCAL_ROW ;
   size_t_vector const& start_A = dA->ROW_DIST->start_of_partition() ;

   intVector exported_items( nb_ranks ) ;
   exported_items.set( 0 ) ;

   // First, lets determine from which processors matrices will be exchanged.
   // We use column scatter to detect local items to send to other processors.
   size_t curr = 0 ;
   size_t end = 0 ;
   for( size_t i=0 ; i<A_col_indexes.size() ; i++ )
   {
      size_t idx = A_col_indexes(i) ;
      if( end <= idx )
      {
         while( (end = start( curr ) + part( curr ) ) <= idx )
         {
            MAC_CHECK( curr < nb_ranks ) ;
            curr++ ;
         }
      }
      exported_items( curr ) = 1 ;
   }
   exported_items(rank) = 0 ;

   intVector exchanged_items( nb_ranks ) ;
   communicator()->all_to_all( exported_items, exchanged_items ) ;

   MAC_Vector * mat_vec = MAC_Vector::create( 0, nb_ranks ) ;
   for( size_t i=0 ; i<nb_ranks ; i++ )
      if( exported_items(i) )
         mat_vec->set_at( i, LA_PelMatrix::create( 0, dA->LOCAL_NODIAG_MATRIX->nb_rows(), part(i) ) ) ;

   if( has_A_items_localy )
   {
      // Now, we're going to fill matrices to send
      LA_MatrixIterator* it = dA->LOCAL_NODIAG_MATRIX->create_stored_item_iterator( 0 ) ;
      for( it->start_all_items() ; it->is_valid() ; it->go_next() )
      {
         size_t col = it->col() ;
         size_t dest = dist->rank_of(col) ;
         if( dest==rank )
         {
            MAC_CHECK( !SQUARE ) ;
            LOCAL_NODIAG_MATRIX->add_to_item( col-FIRST_LOCAL_ROW, it->row()+FIRST_LOCAL_ROW_A,it->item()*alpha ) ;
         }
         else
         {
            LA_PelMatrix * mat =
               static_cast<LA_PelMatrix*>(mat_vec->at(dest)) ;
              MAC_CHECK( mat!=0 ) ;
            mat->set_item(
               it->row(), it->col() - start(dest), it->item() ) ;
         }
      }
      it->destroy() ;

      // Then, we're start to send them
      for( size_t i=0 ; i<nb_ranks ; i++ )
      {
         if( exported_items( i ) > 0 )
         {
            LA_SeqMatrix * sub = static_cast<LA_SeqMatrix*>(mat_vec->at(i));
            sub->synchronize() ;
            LA_CRSmatrix const* crs_matrix = convert_if_needed( sub ) ;
            bool in_place = true ;
            crs_matrix->send( communicator(), i, in_place ) ;
            sub->destroy() ;
            mat_vec->set_at(i,const_cast<LA_CRSmatrix* >(crs_matrix)) ;
         }
      }

      // We deal with local diagonal part
      if( LOCAL_DIAG_MATRIX != 0 )
      {
         LOCAL_DIAG_MATRIX->add_tMat( dA->LOCAL_DIAG_MATRIX, alpha ) ;
      }
   }

   // We now receive sended matrices
   for( size_t i=0 ; i<nb_ranks ; i++ )
   {
      if( exchanged_items( i ) !=0 )
      {
         LA_CRSmatrix* crs_matrix = LA_CRSmatrix::receive( 0, communicator(), i ) ;
         LA_MatrixIterator* it = crs_matrix->create_stored_item_iterator( crs_matrix ) ;

         for( it->start_all_items() ; it->is_valid() ; it->go_next() )
            LOCAL_NODIAG_MATRIX->add_to_item( it->col(), it->row()+start_A(i), alpha*it->item() ) ;
         crs_matrix->destroy() ;
      }
   }

  // Then, we finalize matrix sending
   for( size_t i=0 ; i<nb_ranks ; i++ )
   {
      if( exported_items( i ) != 0 )
      {
         LA_CRSmatrix * crs = static_cast<LA_CRSmatrix*>(mat_vec->at(i));
         crs->wait_send( communicator() ) ;
         crs->destroy() ;
      }
   }
   mat_vec->destroy() ;

   SCATTER_OK = false ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_tMat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
MAC_DistributedPartition const*
LA_DistMatrix:: row_distribution( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: row_distribution" );
   MAC_CHECK_PRE( row_distribution_PRE() ) ;

   MAC_DistributedPartition const* result = ROW_DIST ;

   MAC_CHECK_POST( row_distribution_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DistributedPartition const*
LA_DistMatrix:: col_distribution( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: col_distribution" );
   MAC_CHECK_PRE( col_distribution_PRE() ) ;

   MAC_DistributedPartition const* result = COL_DIST ;

   MAC_CHECK_POST( col_distribution_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrix const*
LA_DistMatrix:: convert_if_needed( LA_SeqMatrix const* mat )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: convert_if_needed" ) ;
   MAC_CHECK_PRE( mat->state() == LA::Sync ) ;

   LA_CRSmatrix const* result = dynamic_cast<LA_CRSmatrix const*>(mat) ;

   if( result==0 ) result = LA_CRSmatrix::create( 0, mat ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->nb_rows() == mat->nb_rows() ) ;
   MAC_CHECK_POST( result->nb_cols() == mat->nb_cols() ) ;
   MAC_CHECK_POST(
      IMPLIES( dynamic_cast<LA_CRSmatrix const*>(mat) != 0,
               result == mat ) ) ;
   MAC_CHECK_POST( result->state() == LA::Sync ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: scale_as_diag_mat_mat( LA_Vector const* lvec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: scale_as_diag_mat_mat LA_Vector" ) ;
   MAC_CHECK_PRE( scale_as_diag_mat_mat_PRE( lvec ) ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const* >( lvec ) != 0 ) ;
   LA_DistVector const* dx = static_cast<LA_DistVector const* >( lvec ) ;

   if( LOCAL_DIAG_MATRIX!=0 )
      LOCAL_DIAG_MATRIX->scale_as_diag_mat_mat( dx->local_vector() ) ;
   LOCAL_NODIAG_MATRIX->scale_as_diag_mat_mat( dx->local_vector() ) ;

   MAC_CHECK_POST( scale_as_diag_mat_mat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: scale_as_mat_diag_mat( LA_Vector const* rvec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: scale_as_mat_diag_mat " ) ;
   MAC_CHECK_PRE( scale_as_mat_diag_mat_PRE( rvec ) ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const* >( rvec ) != 0 ) ;
   LA_DistVector const* dvec = static_cast<LA_DistVector const* >( rvec ) ;
   if( !SCATTER_OK ) build_scatter() ;
   SCATTER->get( rvec, GLOBAL_VEC ) ;
   LA_SeqVector const* loc = dvec->local_vector() ;

   start_local_modifs() ;
   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
   for( it->start_all_items() ; it->is_valid() ; it->go_next() )
   {
      size_t col = it->col() ;
      double d_j = ( SQUARE && col>=FIRST_LOCAL_ROW && col<LAST_LOCAL_ROW ?
                     loc->item( col-FIRST_LOCAL_ROW ) :
                     GLOBAL_VEC->item( col ) ) ;
      it->set_item( it->item() * d_j ) ;
   }
   stop_local_modifs() ;
   it->destroy() ;

   MAC_CHECK_POST( scale_as_mat_diag_mat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: add_to_diag( LA_Vector const* vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: add_to_diag LA_Vector" ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;
   MAC_CHECK_PRE( add_to_diag_PRE( vec ) ) ;

   MAC_CHECK( dynamic_cast<LA_DistVector const* >( vec ) != 0 ) ;
   LA_DistVector const* dx = static_cast<LA_DistVector const* >( vec ) ;
   LA_SeqVector const* lavec = dx->local_vector() ;

   diagonal_block_matrix()->add_to_diag( lavec ) ;

   MAC_CHECK_POST( add_to_diag_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: print" ) ;

   LA_Matrix::print( os, indent_width ) ;
   std::string s( indent_width+3, ' ' ) ;
   if( INITIAL_PROTO != 0 )
   {
      os << s << "initial prototype:" << std::endl ;
      INITIAL_PROTO->print( os, indent_width+6 ) ;
   }
   if( FINAL_PROTO != 0 )
   {
      os << s << "final prototype:" << std::endl ;
      FINAL_PROTO->print( os, indent_width+6 ) ;
   }
}

//----------------------------------------------------------------------
bool
LA_DistMatrix:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Matrix::invariant() ) ;

   if( INITIALIZED )
   {
      MAC_ASSERT( distribution_strategy() == LA::FromGlobalSize ||
                  distribution_strategy() == LA::FromLocalSize ) ;
       MAC_ASSERT( ROW_DIST!=0 ) ;
      MAC_ASSERT( COL_DIST!=0 ) ;
      MAC_ASSERT( ROW_DIST->global_number()==nb_rows() ) ;
      MAC_ASSERT( COL_DIST->global_number()==nb_cols() ) ;
      MAC_ASSERT(
         EQUIVALENT(
            SQUARE,
            nb_rows() == nb_cols() ) ) ;
      MAC_ASSERT(
         EQUIVALENT(
            SQUARE,
            LOCAL_DIAG_MATRIX!=0 &&
            LOCAL_DIAG_MATRIX->nb_rows()==ROW_DIST->local_number() &&
            LOCAL_DIAG_MATRIX->nb_cols()==ROW_DIST->local_number() ) ) ;
      MAC_ASSERT(
         LOCAL_NODIAG_MATRIX!=0 &&
         LOCAL_NODIAG_MATRIX->nb_rows()==ROW_DIST->local_number() &&
         LOCAL_NODIAG_MATRIX->nb_cols()==nb_cols() ) ;
      MAC_ASSERT( NON_LOCAL_MATRIX->nb_rows()==nb_rows() ) ;
      MAC_ASSERT( NON_LOCAL_MATRIX->nb_cols()==nb_cols() ) ;
      if( state() == LA::Sync )
      {
         LA_MatrixIterator* it =
                    NON_LOCAL_MATRIX->create_stored_item_iterator( 0 ) ;
         MAC_ASSERT(
            FORALL( ( it->start_all_items() ; it->is_valid() ; it->go_next() ),
                    ( it->item()==MAC::bad_double() ) ) ) ;
         it->destroy() ;
      }
   }
   return( true ) ;
}

//----------------------------------------------------------------------
double
LA_DistMatrix:: two_norm( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: two_norm" ) ;
   MAC_CHECK_PRE( is_synchronized() ) ;

   double result = 0.0 ;

   LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;

   for( it->start_all_items() ; it->is_valid() ; it->go_next() )
   {
      result += MAC::sqr( it->item() ) ;
   }
   it->destroy() ;

   result = MAC::sqrt( communicator()->sum( result ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_DistMatrix:: implementation_POST( LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Matrix::implementation_POST( result ) ) ;
   MAC_ASSERT( result == LA_DistImplementation::object() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
void
LA_DistMatrix:: build_scatter( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistMatrix:: build_scatter" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK( !SCATTER_OK ) ;
   MAC_CHECK( SCATTER != 0 ) ;

   size_t_vector receive_index( 0 ) ;
   if( COL_DIST->local_number() > 0 )
   {
      LA_MatrixIterator* it =
                 LOCAL_NODIAG_MATRIX->create_stored_item_iterator( 0 ) ;
      boolVector need( nb_cols(), false ) ;
      size_t N = 0 ;
      for( it->start_all_items() ; it->is_valid() ; it->go_next() )
      {
         size_t j = it->col() ;
         if( !need(j) ) N++ ;
         need( j ) = true ;
      }
      it->destroy() ; it=0 ;
      receive_index.resize( N ) ;
      size_t idx = 0 ;
      for( size_t i=0 ; i<need.size() ; i++ )
         if( need(i) ) receive_index(idx++) = i ;
      MAC_ASSERT( idx==N ) ;

   }
   SCATTER->set_sorted( receive_index, receive_index ) ;
   SCATTER_OK = true ;

   MAC_CHECK( SCATTER_OK ) ;
}
