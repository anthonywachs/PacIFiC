#include <LA_PelMatrixIterator.hh>

#include <MAC_assertions.hh>
#include <MAC.hh>
#include <LA_PelMatrix.hh>

//----------------------------------------------------------------------
LA_PelMatrixIterator*
LA_PelMatrixIterator:: create( MAC_Object* a_owner,
                               LA_PelMatrix const* A )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: create" ) ;
   MAC_CHECK_PRE( A != 0 ) ;

   LA_PelMatrixIterator* result =  new LA_PelMatrixIterator( a_owner, A ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->matrix() == A ) ;
   MAC_CHECK_POST( result->is_valid() == false ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
LA_PelMatrixIterator:: LA_PelMatrixIterator( MAC_Object* a_owner,
                                             LA_PelMatrix const* A )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , CUR_ELEM( 0 )
   , I_ROW( 0 )
   , STAY_IN_ROW( false )
{
}

//----------------------------------------------------------------------
LA_PelMatrixIterator:: ~LA_PelMatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_PelMatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: is_valid" ) ;

   bool result = false ;
   if( MAT->ROW_TABLE != 0 )
   {
      result = ( CUR_ELEM != 0 ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: start_all_items" ) ;
   MAC_CHECK_PRE( start_all_items_PRE() ) ;

   I_ROW = 0 ;
   STAY_IN_ROW = false ;
   CUR_ELEM = 0 ;
   if( MAT->ROW_TABLE != 0 )
   {
      CUR_ELEM = MAT->ROW_TABLE[ I_ROW ] ;
      if( CUR_ELEM == 0 )
      {
         go_next_element() ;
      }
   }
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: start_row_items" ) ;
   MAC_CHECK_PRE( start_row_items_PRE( i_row ) ) ;

   I_ROW = i_row ;
   STAY_IN_ROW = true ;
   CUR_ELEM = 0 ;
   if( MAT->ROW_TABLE != 0 )
   {
      CUR_ELEM = MAT->ROW_TABLE[ I_ROW ] ;
   }
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: go_next" ) ;
   MAC_CHECK_PRE( go_next_PRE() ) ;

   go_next_element() ;
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_PelMatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: matrix" ) ;

   LA_Matrix const* result = MAT ;

   MAC_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_PelMatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
size_t
LA_PelMatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: row" ) ;
   MAC_CHECK_PRE( row_PRE() ) ;

   size_t result = I_ROW ;

   MAC_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_PelMatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: col" ) ;
   MAC_CHECK_PRE( col_PRE() ) ;

   size_t result =  CUR_ELEM->iCol ;

   MAC_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_PelMatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;

   return( CUR_ELEM->xVal ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: set_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   MAC_CHECK_PRE( set_item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   CUR_ELEM->xVal = x ;
   if( MAT->UNSYNCHRO && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_set() ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: add_to_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   MAC_CHECK_PRE( add_to_item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   CUR_ELEM->xVal += x ;
   if( MAT->UNSYNCHRO && !MAT->only_local_modifs() )
   {
      unsynchronized_matrix_state_for_add() ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrixIterator:: go_next_element( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrixIterator:: go_next_element" ) ;

   if( CUR_ELEM!=0 )
   {
      CUR_ELEM = CUR_ELEM->next ;
   }
   if( !STAY_IN_ROW )
   {
      while( I_ROW<MAT->nb_rows()-1 && CUR_ELEM==0 )
      {
         I_ROW++ ;
         CUR_ELEM = MAT->ROW_TABLE[ I_ROW ] ;
      }
   }
}
