#include <LA_BlockSeqMatrixIterator.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <LA.hh>
#include <LA_BlockSeqMatrix.hh>

//----------------------------------------------------------------------
LA_BlockSeqMatrixIterator*
LA_BlockSeqMatrixIterator:: create(
                       MAC_Object* a_owner, LA_BlockSeqMatrix const* A )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: create" ) ;
   MAC_CHECK_PRE( A != 0 ) ;
   
   LA_BlockSeqMatrixIterator* result =
                   new LA_BlockSeqMatrixIterator( a_owner, A ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->matrix() == A ) ;
   MAC_CHECK_POST( !result->is_valid() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_BlockSeqMatrixIterator:: LA_BlockSeqMatrixIterator(
                       MAC_Object* a_owner, LA_BlockSeqMatrix const* A )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , MAT_IT( 0 )
   , IS_VALID( false )
   , IB( MAC::bad_index() )
   , JB( MAC::bad_index() )
   , STAY_IN_ROW( false )
   , IROW( MAC::bad_index() )
   , ROW_SHIFT( MAC::bad_index() )
   , COL_SHIFT( MAC::bad_index() )
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: LA_BlockSeqMatrixIterator" ) ;
}

//----------------------------------------------------------------------
LA_BlockSeqMatrixIterator:: ~LA_BlockSeqMatrixIterator( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: ~LA_BlockSeqMatrixIterator" ) ;
   
   if( MAT_IT!=0 )
   {
      MAT_IT->destroy() ;
      MAT_IT = 0 ;
   }   
}

//----------------------------------------------------------------------
bool
LA_BlockSeqMatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( IS_VALID ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: start_all_items" ) ;
   MAC_CHECK_PRE( start_all_items_PRE() ) ;
   
   IB = 0 ;
   JB = 0 ;
   if( MAT_IT != 0 )
   {
      MAT_IT->destroy() ;
      MAT_IT = 0 ;
   }
   STAY_IN_ROW = false ;
   IROW = MAC::bad_index() ;

   IS_VALID = next() ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: start_row_items" ) ;
   MAC_CHECK_PRE( start_row_items_PRE( i_row ) ) ;
   
   IB = MAT->ROW_ELEM_2_BLOCK( i_row ) ;
   JB = 0 ;
   if( MAT_IT != 0 )
   {
      MAT_IT->destroy() ;
      MAT_IT = 0 ;
   }
   STAY_IN_ROW = true ;
   IROW = i_row-MAT->ROW_BLOCK_2_FIRST_ELEM( IB ) ;

   IS_VALID = next() ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: go_next" ) ;
   MAC_CHECK_PRE( go_next_PRE() ) ;

   IS_VALID = next() ;
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_BlockSeqMatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: matrix" ) ;
   
   LA_Matrix const* result = MAT ;
   
   MAC_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: nb_rows" ) ;
   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: row" ) ;
   MAC_CHECK_PRE( row_PRE() ) ;

   MAC_CHECK( MAT_IT!=0 ) ;
   size_t result = MAT_IT->row()+ROW_SHIFT ;

   MAC_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_BlockSeqMatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: col" ) ;
   MAC_CHECK_PRE( col_PRE() ) ;
   
   MAC_CHECK( MAT_IT!=0 ) ;
   size_t result = MAT_IT->col() + COL_SHIFT ;
   
   MAC_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_BlockSeqMatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;

   MAC_CHECK( MAT_IT!=0 ) ;   
   MAC_CHECK( MAT_IT->item()==MAT->item( row(), col() ) ) ;
   
   return( MAT_IT->item() ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: set_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   MAC_CHECK_PRE( set_item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( MAT_IT!=0 ) ;
   MAT_IT->set_item( x ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_BlockSeqMatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: add_to_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   MAC_CHECK_PRE( add_to_item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( MAT_IT!=0 ) ;
   MAT_IT->add_to_item( x ) ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
bool
LA_BlockSeqMatrixIterator:: next( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_BlockSeqMatrixIterator:: next" ) ;
   
   bool result = false ;
   if( MAT_IT!=0 )
   {
      MAT_IT->go_next() ;
      result = MAT_IT->is_valid() ;
      if( !result )
      {
         MAT_IT->destroy() ;
         MAT_IT = 0 ;
         JB++ ;
      }
   }
   for( ; !result && IB<MAT->ROW_PARTITION.size() ; ++IB )
   {
      for( ; !result && JB<MAT->COL_PARTITION.size() ; ++JB )
      {
         if( MAT->has_submatrix( IB, JB ) )
         {
            LA_SeqMatrix* matS = MAT->submat( IB, JB ) ;
            if( !matS->is_synchronized() ) matS->synchronize() ;
            MAT_IT = matS->create_stored_item_iterator( 0 ) ;
            if( STAY_IN_ROW )
            {
               MAT_IT->start_row_items( IROW ) ;
            }
            else
            {
               MAT_IT->start_all_items() ;
            }
            result = MAT_IT->is_valid() ;
            if( !result )
            {
               MAT_IT->destroy() ;
               MAT_IT = 0 ;
            }
            else
            {
               ROW_SHIFT = MAT->ROW_BLOCK_2_FIRST_ELEM( IB ) ;
               COL_SHIFT = MAT->COL_BLOCK_2_FIRST_ELEM( JB ) ;
            }
         }
         if( result  )
         {
            break ;
         }
      }
      if( result || STAY_IN_ROW )
      {
         break ;
      }
      JB = 0 ;
   }
   MAC_CHECK( !result || MAT_IT!=0 ) ;
   return( result ) ;
}
   




