#include <LA_ShiftedIndexMatrixIterator.hh>

#include <MAC_assertions.hh>
#include <MAC.hh>
#include <LA_Matrix.hh>

//----------------------------------------------------------------------
LA_ShiftedIndexMatrixIterator*
LA_ShiftedIndexMatrixIterator:: create( MAC_Object* a_owner,
                                        size_t row_shift,
                                        size_t col_shift,
                                        LA_Matrix const* A,
                                        LA_MatrixIterator* internal )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_ShiftedIndexMatrixIterator:: create" ) ;
   MAC_CHECK_PRE( A != 0 ) ;
   MAC_CHECK_PRE( internal != 0 ) ;
   MAC_CHECK_PRE( internal->owner()==0 ) ;
   MAC_CHECK_PRE( internal->matrix()->nb_rows()+row_shift <= A->nb_rows() ) ;
   MAC_CHECK_PRE( internal->matrix()->nb_cols()+col_shift <= A->nb_cols() ) ;

   LA_ShiftedIndexMatrixIterator* result =
      new LA_ShiftedIndexMatrixIterator( a_owner,
                                         row_shift, col_shift, A,
                                         internal  ) ;
   internal->set_owner( result ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->matrix() == A ) ;
   MAC_CHECK_POST( result->is_valid() == false ) ;
   MAC_CHECK_PRE( internal->owner()==result ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_ShiftedIndexMatrixIterator:: LA_ShiftedIndexMatrixIterator(
                                            MAC_Object* a_owner,
                                            size_t row_shift,
                                            size_t col_shift,
                                            LA_Matrix const* A,
                                            LA_MatrixIterator* internal )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , ROW_SHIFT( row_shift )
   , COL_SHIFT( col_shift )
   , INTERNAL( internal )
   , NULLROW(false)
{
}

//----------------------------------------------------------------------
LA_ShiftedIndexMatrixIterator:: ~LA_ShiftedIndexMatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_ShiftedIndexMatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( !NULLROW && INTERNAL->is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_ShiftedIndexMatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_ShiftedIndexMatrixIterator:: start_all_items" ) ;
   MAC_CHECK_PRE( start_all_items_PRE() ) ;

   NULLROW = false ;
   INTERNAL->start_all_items() ;
}

//----------------------------------------------------------------------
void
LA_ShiftedIndexMatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_ShiftedIndexMatrixIterator:: start_row_items" ) ;
   MAC_CHECK_PRE( start_row_items_PRE( i_row ) ) ;

   NULLROW = i_row<ROW_SHIFT || i_row>=INTERNAL->nb_rows()+ROW_SHIFT ;
   if( !NULLROW )
   {
      INTERNAL->start_row_items( (size_t)( (int)i_row-(int)ROW_SHIFT ) );
   }
}

//----------------------------------------------------------------------
void
LA_ShiftedIndexMatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_ShiftedIndexMatrixIterator:: go_next" ) ;
   MAC_CHECK_PRE( go_next_PRE() ) ;

   INTERNAL->go_next() ;
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_ShiftedIndexMatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_ShiftedIndexMatrixIterator:: matrix" ) ;

   LA_Matrix const* result = MAT ;

   MAC_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_ShiftedIndexMatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
size_t
LA_ShiftedIndexMatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_ShiftedIndexMatrixIterator:: row" ) ;
   MAC_CHECK_PRE( row_PRE() ) ;

   size_t result = INTERNAL->row()+ROW_SHIFT ;

   MAC_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_ShiftedIndexMatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_ShiftedIndexMatrixIterator:: col" ) ;
   MAC_CHECK_PRE( col_PRE() ) ;

   size_t result = INTERNAL->col()+COL_SHIFT ;

   MAC_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_ShiftedIndexMatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_ShiftedIndexMatrixIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;

   return( INTERNAL->item() ) ;
}

//----------------------------------------------------------------------
void
LA_ShiftedIndexMatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_ShiftedIndexMatrixIterator:: set_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   MAC_CHECK_PRE( set_item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   INTERNAL->set_item( x ) ;
   if( MAT->is_desynchronizable() && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_set() ;
   }


   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_ShiftedIndexMatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_ShiftedIndexMatrixIterator:: add_to_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   MAC_CHECK_PRE( add_to_item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   INTERNAL->add_to_item( x ) ;
   if( MAT->is_desynchronizable() && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_add() ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}
