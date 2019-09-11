#include <LA_PairOfMatrixIterator.hh>

#include <MAC_assertions.hh>
#include <MAC.hh>
#include <LA_Matrix.hh>

//----------------------------------------------------------------------
LA_PairOfMatrixIterator*
LA_PairOfMatrixIterator:: create( MAC_Object* a_owner,
                                  LA_Matrix const* A,
                                  LA_MatrixIterator* first,
                                  LA_MatrixIterator* second )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: create" ) ;
   MAC_CHECK_PRE( A != 0 ) ;
   MAC_CHECK_PRE( first->owner() == 0 ) ;
   MAC_CHECK_PRE( first->nb_rows() == A->nb_rows() ) ;
   MAC_CHECK_PRE( second->owner()==0 ) ;
   MAC_CHECK_PRE( second->nb_rows() == A->nb_rows() ) ;


   LA_PairOfMatrixIterator* result =
         new LA_PairOfMatrixIterator( a_owner, A, first, second  ) ;
   first->set_owner( result ) ;
   second->set_owner( result ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->matrix() == A ) ;
   MAC_CHECK_POST( result->is_valid() == false ) ;
   MAC_CHECK_PRE( first->owner()==result ) ;
   MAC_CHECK_PRE( second->owner()==result ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_PairOfMatrixIterator:: LA_PairOfMatrixIterator( MAC_Object* a_owner,
                                                   LA_Matrix const* A,
                                                   LA_MatrixIterator* first,
                                                   LA_MatrixIterator* second )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , NB_ROWS( first->nb_rows() )
   , FIRST( first )
   , SECOND( second )
   , CURR( 0 )
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: LA_PairOfMatrixIterator" ) ;
}

//----------------------------------------------------------------------
LA_PairOfMatrixIterator:: ~LA_PairOfMatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_PairOfMatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( CURR!=0 && CURR->is_valid() ) ;
}

//----------------------------------------------------------------------
void
LA_PairOfMatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: start_all_items" ) ;
   MAC_CHECK_PRE( start_all_items_PRE() ) ;

   FIRST->start_all_items() ;
   SECOND->start_all_items() ;
   CURR=( FIRST->is_valid() ? FIRST : SECOND ) ;
}

//----------------------------------------------------------------------
void
LA_PairOfMatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: start_row_items" ) ;
   MAC_CHECK_PRE( start_row_items_PRE( i_row ) ) ;

   FIRST->start_row_items(i_row) ;
   SECOND->start_row_items(i_row) ;
   CURR=( FIRST->is_valid() ? FIRST : SECOND ) ;
}

//----------------------------------------------------------------------
void
LA_PairOfMatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: go_next" ) ;
   MAC_CHECK_PRE( go_next_PRE() ) ;

   CURR->go_next() ;
   if( !CURR->is_valid() && CURR!=SECOND )
   {
      CURR = SECOND ;
   }
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_PairOfMatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: matrix" ) ;

   LA_Matrix const* result = MAT ;

   MAC_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_PairOfMatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( NB_ROWS ) ;
}

//----------------------------------------------------------------------
size_t
LA_PairOfMatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: row" ) ;
   MAC_CHECK_PRE( row_PRE() ) ;

   size_t result = CURR->row() ;

   MAC_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_PairOfMatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: col" ) ;
   MAC_CHECK_PRE( col_PRE() ) ;

   size_t result = CURR->col() ;

   MAC_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_PairOfMatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;

   return( CURR->item() ) ;
}

//----------------------------------------------------------------------
void
LA_PairOfMatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: set_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   MAC_CHECK_PRE( set_item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   CURR->set_item(x) ;
   if( MAT->is_desynchronizable() && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_set() ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_PairOfMatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PairOfMatrixIterator:: add_to_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   MAC_CHECK_PRE( add_to_item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   CURR->add_to_item(x) ;
   if( MAT->is_desynchronizable() && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_add() ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}
