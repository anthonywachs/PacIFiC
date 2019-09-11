#include <LA_CRSmatrixIterator.hh>

#include <MAC_assertions.hh>
#include <MAC.hh>
#include <LA_CRSmatrix.hh>

//----------------------------------------------------------------------
LA_CRSmatrixIterator*
LA_CRSmatrixIterator:: create( MAC_Object* a_owner,
                               LA_CRSmatrix const* A )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_CRSmatrixIterator:: create" ) ;
   MAC_CHECK_PRE( A != 0 ) ;

   LA_CRSmatrixIterator* result =  new LA_CRSmatrixIterator( a_owner, A ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->matrix() == A ) ;
   MAC_CHECK_POST( result->is_valid() == false ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_CRSmatrixIterator:: LA_CRSmatrixIterator( MAC_Object* a_owner,
                                             LA_CRSmatrix const* A )
//----------------------------------------------------------------------
   : LA_MatrixIterator( a_owner )
   , MAT( A )
   , VALID( false )
   , I_ROW( 0 )
   , STAY_IN_ROW( false )
   , I_CUR( MAC::bad_index() )
{
}

//----------------------------------------------------------------------
LA_CRSmatrixIterator:: ~LA_CRSmatrixIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
LA_CRSmatrixIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( VALID ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrixIterator:: start_all_items( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_CRSmatrixIterator:: start_all_items" ) ;
   MAC_CHECK_PRE( start_all_items_PRE() ) ;
   
   STAY_IN_ROW = false ;
   I_CUR = 0 ;
   I_ROW = 0 ;
   VALID = ( I_CUR< (int)MAT->NB_ELEMS ) ;
   if( VALID )
   {
      while( MAT->START(I_ROW+1)<=I_CUR ) I_ROW++ ;
      MAC_CHECK( I_CUR >=  MAT->START(I_ROW) ) ;
      MAC_CHECK( I_CUR <  MAT->START(I_ROW+1) ) ;
   }
}

//----------------------------------------------------------------------
void
LA_CRSmatrixIterator:: start_row_items( size_t i_row )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_CRSmatrixIterator:: start_row_items" ) ;
   MAC_CHECK_PRE( start_row_items_PRE( i_row ) ) ;

   I_ROW = i_row ;
   STAY_IN_ROW = true ;
   I_CUR = MAT->START(i_row) ;
   VALID = I_CUR < MAT->START(i_row+1) ;        
}

//----------------------------------------------------------------------
void
LA_CRSmatrixIterator:: go_next( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_CRSmatrixIterator:: go_next" ) ;
   MAC_CHECK_PRE( go_next_PRE() ) ;

   I_CUR++ ;
   if( STAY_IN_ROW )
   {
      VALID = ( I_CUR<MAT->START(I_ROW+1) ) ;
   }
   else
   {
      VALID = ( I_CUR<(int)MAT->NB_ELEMS ) ;
      if( VALID )
      {
         while( MAT->START(I_ROW+1)<=I_CUR ) I_ROW++ ;
         MAC_CHECK( I_CUR >=  MAT->START(I_ROW) ) ;
         MAC_CHECK( I_CUR <  MAT->START(I_ROW+1) ) ;
      }
   }
}

//----------------------------------------------------------------------
LA_Matrix const*
LA_CRSmatrixIterator:: matrix( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_CRSmatrixIterator:: matrix" ) ;
   
   LA_Matrix const* result = MAT ;
   
   MAC_CHECK_POST( matrix_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_CRSmatrixIterator:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( MAT->nb_rows() ) ;
}

//----------------------------------------------------------------------
size_t
LA_CRSmatrixIterator:: row( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_CRSmatrixIterator:: row" ) ;
   MAC_CHECK_PRE( row_PRE() ) ;

   size_t result = I_ROW  ;
   
   MAC_CHECK_POST( row_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_CRSmatrixIterator:: col( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_CRSmatrixIterator:: col" ) ;
   MAC_CHECK_PRE( col_PRE() ) ;

   size_t result = MAT->COL(I_CUR) ;

   MAC_CHECK_POST( col_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
LA_CRSmatrixIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_CRSmatrixIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;

   return( MAT->VALUES(I_CUR) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrixIterator:: set_item( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_CRSmatrixIterator:: set_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   MAC_CHECK_PRE( set_item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_CRSmatrix* dummy = const_cast<LA_CRSmatrix*>( MAT ) ;
   dummy->VALUES(I_CUR) = x ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_item_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_CRSmatrixIterator:: add_to_item( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_CRSmatrixIterator:: add_to_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, matrix()->state() ) ;
   MAC_CHECK_PRE( add_to_item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_CRSmatrix* dummy = const_cast<LA_CRSmatrix*>( MAT ) ;
   dummy->VALUES(I_CUR) += x ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_to_item_POST( OLD( state ) ) ) ;
}
