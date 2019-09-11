#include <LA_MatrixIterator.hh>

#include <MAC_assertions.hh>
#include <MAC_DistributedPartition.hh>

#include <LA.hh>
#include <LA_Matrix.hh>

//----------------------------------------------------------------------
LA_MatrixIterator:: LA_MatrixIterator( MAC_Object* a_owner )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
{
}

//----------------------------------------------------------------------
LA_MatrixIterator:: ~LA_MatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_MatrixIterator:: unsynchronized_matrix_state_for_set( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_MatrixIterator:: unsynchronized_matrix_state_for_set" ) ;
   MAC_CHECK_PRE( matrix()->is_desynchronizable() ) ;
   MAC_CHECK_PRE( !matrix()->only_local_modifs() ) ;
   MAC_CHECK_PRE( !matrix()->only_local_modifs() ) ;
   MAC_CHECK_PRE( matrix()->state() != LA::NotSync_add ) ;

   LA_Matrix* m = const_cast<LA_Matrix*>( matrix() ) ;
   m->set_unsynchronized_state( LA::NotSync_set ) ;

   MAC_CHECK_POST( matrix()->state() == LA::NotSync_set ) ;
}

//----------------------------------------------------------------------
void
LA_MatrixIterator:: unsynchronized_matrix_state_for_add( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_MatrixIterator:: unsynchronized_matrix_state_for_add" ) ;
   MAC_CHECK_PRE( matrix()->is_desynchronizable() ) ;
   MAC_CHECK_PRE( !matrix()->only_local_modifs() ) ;
   MAC_CHECK_PRE( matrix()->state() != LA::NotSync_set ) ;

   LA_Matrix* m = const_cast<LA_Matrix*>( matrix() ) ;
   m->set_unsynchronized_state( LA::NotSync_add ) ;

   MAC_CHECK_POST( matrix()->state() == LA::NotSync_add ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: start_all_items_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( matrix()->state() ==  LA::Sync ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: start_row_items_PRE( size_t i_row ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( matrix()->state() == LA::Sync ) ;
   MAC_ASSERT( i_row < nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: go_next_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_valid() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: matrix_POST( LA_Matrix const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->nb_rows() == nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: row_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_valid() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: row_POST( size_t result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result < matrix()->nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: col_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_valid() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: col_POST( size_t result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result < matrix()->nb_cols() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
double
LA_MatrixIterator:: item_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_valid() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: set_item_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_valid() ) ;
   MAC_ASSERT(
       IMPLIES( matrix()->only_local_modifs(),
                row() >= matrix()->row_distribution()->first_local_index() &&
                row() < matrix()->row_distribution()->local_index_limit() ) ) ;
   MAC_ASSERT( matrix()->state() != LA::NotSync_add ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: set_item_POST( LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT(
   IMPLIES( matrix()->is_desynchronizable() && !matrix()->only_local_modifs(),
            matrix()->state() == LA::NotSync_set ) ) ;
   MAC_ASSERT(
   IMPLIES( !matrix()->is_desynchronizable() || matrix()->only_local_modifs(),
            matrix()->state() == old_state ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: add_to_item_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_valid() ) ;
   MAC_ASSERT(
       IMPLIES( matrix()->only_local_modifs(),
                row() >= matrix()->row_distribution()->first_local_index() &&
                row() < matrix()->row_distribution()->local_index_limit() ) ) ;
   MAC_ASSERT( matrix()->state() != LA::NotSync_set ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_MatrixIterator:: add_to_item_POST( LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT(
      IMPLIES( matrix()->is_desynchronizable() && !matrix()->only_local_modifs(),
               matrix()->state() == LA::NotSync_add ) ) ;
   MAC_ASSERT(
   IMPLIES( !matrix()->is_desynchronizable() || matrix()->only_local_modifs(),
            matrix()->state() == old_state ) ) ;
   return( true ) ;
}



