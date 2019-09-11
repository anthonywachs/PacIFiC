#include <LA_PelMatrix.hh>

#include <LA_PelMatrixIterator.hh>
#include <LA_SeqVector.hh>
#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_assertions.hh>

#include <iomanip>
#include <iostream>
#include <fstream>

LA_PelMatrix const* LA_PelMatrix::PROTOTYPE = new LA_PelMatrix() ;

//----------------------------------------------------------------------
LA_PelMatrix*
LA_PelMatrix:: create( MAC_Object* a_owner,
                       size_t a_nb_rows,
                       size_t a_nb_cols )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: create" ) ;

   LA_PelMatrix* result = new LA_PelMatrix( a_owner, a_nb_rows, a_nb_cols ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( !result->is_symmetric() ) ;
   MAC_CHECK_POST( result->is_resizable() ) ;
   MAC_CHECK_POST( !result->is_desynchronizable() ) ;
   MAC_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   MAC_CHECK_POST( result->state() == LA::Sync ) ;
   MAC_CHECK_POST( result->is_synchronized() ) ;
   MAC_CHECK_POST( result->nb_rows() == a_nb_rows ) ;
   MAC_CHECK_POST( result->nb_cols() == a_nb_cols ) ;
   MAC_CHECK_POST( result->nb_stored_items() == 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_PelMatrix:: LA_PelMatrix( MAC_Object* a_owner,
                             size_t a_nb_rows,
                             size_t a_nb_cols )
//----------------------------------------------------------------------
   : LA_SeqMatrix( a_owner, "LA_PelMatrix" )
   , UNSYNCHRO( false )
   , NB_ROWS( a_nb_rows )
   , NB_COLS( a_nb_cols )
   , ROW_TABLE( 0 )
   , FIRST_BUCKET( 0 )
   , CURRENT_BUCKET( 0 )
{
   MAC_LABEL( "LA_PelMatrix:: LA_PelMatrix" ) ;

   if( NB_ROWS > 0 )
   {
      ROW_TABLE = new RowElm* [ NB_ROWS ] ;
      for( size_t i=0 ; i<NB_ROWS ; ++i )
         ROW_TABLE[i] = 0 ;
   }

   init_allocating() ;

   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
LA_PelMatrix*
LA_PelMatrix:: create_matrix( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: create_matrix" ) ;
   MAC_CHECK_PRE( create_matrix_PRE( a_owner ) ) ;

   MAC_CHECK_INV( invariant() ) ;
   LA_PelMatrix* result = new LA_PelMatrix( a_owner, NB_ROWS, NB_COLS ) ;
   if( UNSYNCHRO ) result->make_desynchronizable() ;

   MAC_CHECK_POST( create_matrix_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_PelMatrix*
LA_PelMatrix:: create_replica( MAC_Object* a_owner,
                               MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   int a_nb_rows = 0 ;
   int a_nb_cols = 0 ;
   if( exp->has_entry( "nb_rows" ) )
   {
      a_nb_rows = exp->int_data( "nb_rows" ) ;
      exp->test_data( "nb_rows", "nb_rows>=0" ) ;
      exp->set_default( "nb_rows", "0" ) ;
   }
   if( exp->has_entry( "nb_cols" ) )
   {
      a_nb_cols = exp->int_data( "nb_cols" ) ;
      exp->test_data( "nb_cols", "nb_cols>=0" ) ;
      exp->set_default( "nb_cols", "0" ) ;
   }

   LA_PelMatrix* result =
           LA_PelMatrix::create( a_owner, a_nb_rows, a_nb_cols ) ;

   if(  exp->has_entry( "is_desynchronizable" ) )
   {
      bool const dist = exp->bool_data( "is_desynchronizable" ) ;
      exp->set_default( "is_desynchronizable", "false" ) ;
      if( dist )
      {
         if( MAC_Exec::communicator()->nb_ranks() != 1 )
         {
            MAC_Error::object()->raise_data_error(
               exp, "is_desynchronizable",
               "testing distributed behavior is only allowed on one processor" ) ;
         }
         result->make_desynchronizable() ;
      }
   }

   MAC_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_PelMatrix:: LA_PelMatrix( void )
//----------------------------------------------------------------------
   : LA_SeqMatrix( "LA_PelMatrix" )
   , UNSYNCHRO( false )
   , NB_ROWS( 0 )
   , NB_COLS( 0 )
   , ROW_TABLE( 0 )
   , FIRST_BUCKET( 0 )
   , CURRENT_BUCKET( 0 )
{
   MAC_LABEL( "LA_PelMatrix:: LA_PelMatrix" ) ;
   MAC_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
LA_PelMatrix:: ~LA_PelMatrix( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: ~LA_PelMatrix" ) ;
   MAC_CHECK_INV( invariant() ) ;

   destroy_all_elements() ;
   if( ROW_TABLE != 0 )
   {
      delete [] ROW_TABLE ;
      ROW_TABLE = 0 ;
   }
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: re_initialize_with_global_sizes( size_t a_nb_rows,
                                                size_t a_nb_cols )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: re_initialize_with_global_sizes" ) ;
   MAC_CHECK_PRE( re_initialize_with_global_sizes_PRE( a_nb_rows,
                                                       a_nb_cols ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   destroy_all_elements() ;
   init_allocating() ;

   NB_COLS = a_nb_cols ;

   if( a_nb_rows != NB_ROWS )
   {
      NB_ROWS = a_nb_rows ;

      if( ROW_TABLE != 0 )
      {
         delete [] ROW_TABLE ;
         ROW_TABLE = 0 ;
      }
      if( NB_ROWS > 0 )
      {
         ROW_TABLE = new RowElm* [ NB_ROWS ] ;
      }
   }
   for( size_t iRow=0 ; iRow<NB_ROWS ; iRow++ )
   {
      ROW_TABLE[iRow] = 0 ;
   }
   if( UNSYNCHRO )
   {
      set_unsynchronized_state( LA::NotSync_undef ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( re_initialize_with_global_sizes_POST( a_nb_rows,
                                                         a_nb_cols ) ) ;
}

//----------------------------------------------------------------------
size_t
LA_PelMatrix:: nb_stored_items( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: nb_stored_items" ) ;
   MAC_CHECK_PRE( nb_stored_items_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   size_t result = 0 ;
   for( size_t i=0 ; i<NB_ROWS ; i++ )
   {
      RowElm* elm = ROW_TABLE[i] ;
      while( elm!=0 )
      {
         ++result ;
         elm = elm->next ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_PelMatrix:: allocated_memory( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: allocated_memory" ) ;

   RowElmBucket * ptr = FIRST_BUCKET ;
   size_t result = 0 ;

   while( ptr!=0 )
   {
      ptr = ptr->next ;
      result += sizeof(RowElmBucket) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_PelMatrix:: nb_rows( void ) const
//----------------------------------------------------------------------
{
   return( NB_ROWS ) ;
}

//----------------------------------------------------------------------
size_t
LA_PelMatrix:: nb_cols( void ) const
//----------------------------------------------------------------------
{
   return( NB_COLS ) ;
}

//----------------------------------------------------------------------
double
LA_PelMatrix:: item( size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: item" ) ;
   MAC_CHECK_PRE( item_PRE( i, j ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   double result = 0. ;
   RowElm* elm = ROW_TABLE[i] ;
   while( elm!=0 && elm->iCol<=j )
   {
      if( elm->iCol==j )
      {
         result = elm->xVal ;
         break ;
      }
      elm = elm->next ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: make_desynchronizable( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: make_desynchronizable" ) ;
   MAC_CHECK_PRE( !is_desynchronizable() ) ;

   UNSYNCHRO = true ;
   set_unsynchronized_state( LA::NotSync_undef ) ;

   MAC_CHECK_POST( is_desynchronizable() ) ;
   MAC_CHECK_POST( state() == LA::NotSync_undef ) ;
   MAC_CHECK_POST( ! is_synchronized() ) ;
}

//----------------------------------------------------------------------
bool
LA_PelMatrix:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: is_desynchronizable" ) ;

   bool result = UNSYNCHRO ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_MatrixIterator*
LA_PelMatrix:: create_stored_item_iterator( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: create_stored_item_iterator" ) ;
   MAC_CHECK_PRE( create_stored_item_iterator_PRE( a_owner ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_PelMatrixIterator* result =
                           LA_PelMatrixIterator::create( a_owner, this ) ;

   MAC_CHECK_POST( create_stored_item_iterator_POST( a_owner, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: set_stored_items( double val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: set_stored_items" ) ;
   MAC_CHECK_INV( invariant() ) ;

   for( size_t iRow=0 ; iRow<NB_ROWS ; iRow++ )
   {
      RowElm* elm = ROW_TABLE[iRow] ;
      while( elm!=0 )
      {
         elm->xVal = val ;
         elm = elm->next ;
      }
   }
   if( is_desynchronizable() )
   {
      set_unsynchronized_state( LA::NotSync_undef ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_stored_items_POST( val ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: nullify_row( size_t i )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: nullify_row" ) ;
   MAC_CHECK_PRE( nullify_row_PRE( i ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   RowElm* elm = ROW_TABLE[i] ;
   while( elm!=0 )
   {
      elm->xVal = 0. ;
      elm = elm->next ;
   }
   if( UNSYNCHRO ) set_unsynchronized_state( LA::NotSync_undef ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( nullify_row_POST( i ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: set_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: set_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;
   MAC_CHECK_PRE( set_item_PRE( i, j ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( is_desynchronizable() && !only_local_modifs() )
   {
      set_unsynchronized_state( LA::NotSync_set ) ;
   }

   RowElm* elm = ROW_TABLE[i] ;
   RowElm* dummy = 0 ;
   bool found = false ;
   while( elm!=0 && elm->iCol<=j && !found )
   {
      if( elm->iCol==j )
      {
         found = true ;
      }
      else
      {
         dummy = elm ;
         elm = elm->next ;
      }
   }

   if( !found )
   {
      elm = create_element() ;
      MAC_CHECK( elm != 0 ) ;
      elm->iCol = j ;
      elm->xVal = x ;
      if( dummy == 0 )
      {
         elm->next = ROW_TABLE[i] ;
         ROW_TABLE[i] = elm ;
      }
      else
      {
         elm->next = dummy->next ;
         dummy->next  = elm ;
      }
   }
   else
   {
      elm->xVal = x ;
   }

   MAC_CHECK_POST( set_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: add_to_item( size_t i, size_t j, double x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: add_to_item" ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;
   MAC_CHECK_PRE( add_to_item_PRE( i, j ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( is_desynchronizable() && !only_local_modifs() )
   {
      set_unsynchronized_state( LA::NotSync_add ) ;
   }

   RowElm* elm = ROW_TABLE[i] ;
   RowElm* dummy = 0 ;
   bool found = false ;
   while( elm!=0 && elm->iCol<=j && !found )
   {
      if( elm->iCol==j )
      {
         found = true ;
      }
      else
      {
         dummy = elm ;
         elm = elm->next ;
      }
   }
   if( !found )
   {
         elm = create_element() ;
         MAC_ASSERT( elm != 0 ) ;
         elm->iCol = j ;
         elm->xVal = x ;
         if( dummy==0 )
         {
            elm->next = ROW_TABLE[i] ;
            ROW_TABLE[i] = elm ;
         }
         else
         {
            elm->next = dummy->next ;
            dummy->next  = elm ;
         }
   }
   else
   {
      elm->xVal += x ;
   }

   MAC_CHECK_POST( add_to_item_POST( i, j, OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: multiply_vec_then_add( LA_Vector const* x, LA_Vector* y,
				      double alpha, double beta ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: multiply_vec_then_add" ) ;
   MAC_CHECK_PRE( multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_SeqVector const* bx = static_cast<LA_SeqVector const* >( x ) ;
   MAC_CHECK( dynamic_cast<LA_SeqVector const* >( x ) != 0 ) ;
   LA_SeqVector* by = static_cast<LA_SeqVector* >( y ) ;
   MAC_CHECK( dynamic_cast<LA_SeqVector* >( y ) != 0 ) ;

   double* resu = by->data() ;
   double const* x_data = bx->data() ;
   for( size_t i=0 ; i<NB_ROWS ; i++ )
   {
      RowElm* elm = ROW_TABLE[i] ;
      *resu *= beta ;
      while( elm!=0 )
      {
         *resu += alpha * elm->xVal * x_data[ elm->iCol ] ;
         elm = elm->next ;
      }
      resu++ ;
   }
   y->synchronize() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: tr_multiply_vec_then_add( LA_Vector const* x, LA_Vector* y,
					 double alpha, double beta ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: tr_multiply_vec_then_add" ) ;
   MAC_CHECK_PRE( tr_multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   y->scale( beta ) ;

   LA_SeqVector const* bx = static_cast<LA_SeqVector const* >( x ) ;
   MAC_CHECK( dynamic_cast<LA_SeqVector const* >( x ) != 0 ) ;
   LA_SeqVector* by = static_cast<LA_SeqVector* >( y ) ;
   MAC_CHECK( dynamic_cast<LA_SeqVector* >( y ) != 0 ) ;

   double const* x_data = bx->data() ;
   double* y_data = by->data() ;
   for( size_t i=0 ; i<NB_ROWS ; i++ )
   {
      RowElm* elm = ROW_TABLE[i] ;
      while( elm!=0 )
      {
         MAC_CHECK( elm->iCol<by->nb_rows() ) ;
         y_data[elm->iCol] += alpha * elm->xVal*x_data[i] ;
         elm = elm->next ;
      }
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( tr_multiply_vec_then_add_POST( y ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: set( LA_Matrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: set" ) ;
   MAC_CHECK_PRE( set_PRE( A, same_pattern ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_PelMatrix const* pA = dynamic_cast<LA_PelMatrix const* >( A ) ;
   if( pA==0 )
   {
      LA_SeqMatrix::set( A, same_pattern ) ;
   }
   else
   {
      if( same_pattern )
      {
         for( size_t i=0 ; i<NB_ROWS ; i++ )
         {
            RowElm* elm = ROW_TABLE[i] ;
            for( RowElm* elmA = pA->ROW_TABLE[i] ;
                 elmA != 0 ;
                 elmA = elmA->next )
            {
               MAC_CHECK( elm!=0 ) ;
               MAC_CHECK( elm->iCol == elmA->iCol ) ;
               elm->xVal = elmA->xVal ;
               elm = elm->next ;
            }
         }
      }
      else
      {
         destroy_all_elements() ;
         init_allocating() ;
         for( size_t i=0 ; i<NB_ROWS ; i++ )
         {
            RowElm* elm = 0 ;
            bool first = true ;
            ROW_TABLE[i] = 0 ;
            for( RowElm* elmA = pA->ROW_TABLE[i] ;
                 elmA != 0 ;
                 elmA = elmA->next )
            {
               RowElm* elmNew = create_element() ;
               elmNew->iCol = elmA->iCol ;
               elmNew->xVal = elmA->xVal ;
               elmNew->next = 0 ;
               if( first )
               {
                  ROW_TABLE[i] = elmNew ;
                  elm = elmNew ;
                  first = false ;
               }
               else
               {
                  elm->next = elmNew ;
                  elm = elm->next ;
               }
            }
         }
      }
   }
   //synchronize() ;
   if( is_desynchronizable() ) set_unsynchronized_state( LA::NotSync_undef ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_POST( A ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: scale( double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: scale" ) ;
   MAC_CHECK_PRE( scale_PRE( alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( alpha==0.0 )
   {
      nullify() ;
   }
   else if( alpha!=1.0 )
   {
      for( size_t i=0 ; i<NB_ROWS ; i++ )
      {
         RowElm* elm = ROW_TABLE[i] ;
         while( elm!=0 )
         {
            elm->xVal *= alpha ;
            elm = elm->next ;
         }
      }
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( scale_POST( alpha ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                            double alpha )
//----------------------------------------------------------------------
// M <- alpha*A* B +beta*M
{
   MAC_LABEL( "LA_PelMatrix:: add_Mat_Mat" ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;
   MAC_CHECK_PRE( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_PelMatrix const* AMAC = dynamic_cast<LA_PelMatrix const*>( A ) ;
   LA_PelMatrix const* BMAC = dynamic_cast<LA_PelMatrix const*>( B ) ;

   if( ( AMAC != 0 ) && ( BMAC != 0 ) )
   {
      add_Mat_Mat_IMP( AMAC, BMAC, alpha ) ;
   }
   else
   {
      LA_SeqMatrix::add_Mat_Mat( A, B, alpha ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_Mat_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: add_Mat_tMat( LA_Matrix const* A, LA_Matrix const* B,
                             double alpha )
//----------------------------------------------------------------------
// M <- alpha*A*transpose( B )+beta*M
{
   MAC_LABEL( "LA_PelMatrix:: add_Mat_tMat" ) ;
   MAC_CHECK_PRE( add_Mat_tMat_PRE( A, B, alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_PelMatrix const* AMAC = dynamic_cast<LA_PelMatrix const*>( A ) ;
   LA_PelMatrix const* BMAC = dynamic_cast<LA_PelMatrix const*>( B ) ;

   if( ( AMAC != 0 ) && ( BMAC != 0 ) )
   {
      add_Mat_tMat_IMP( AMAC, BMAC, alpha ) ;
   }
   else
   {
      LA_SeqMatrix::add_Mat_tMat( A, B, alpha ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_Mat_tMat_POST() ) ;
}
//----------------------------------------------------------------------
void
LA_PelMatrix:: factorize_MILU0( bool modified, double piv_min )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: factorize_MILU0" ) ;
   MAC_CHECK_PRE( factorize_MILU0_PRE( modified, piv_min ) ) ;

   size_t const n = nb_rows() ;

   RowElm** diag = new RowElm* [ NB_ROWS ] ;
   for( size_t i=0 ; i<n ; ++i )
   {
      RowElm* elm = ROW_TABLE[i] ;
      RowElm* dummy = 0 ;
      bool found = false ;
      while( elm!=0 && elm->iCol<=i && !found )
      {
         if( elm->iCol==i )
         {
            found = true ;
         }
         else
         {
            dummy = elm ;
            elm = elm->next ;
         }
      }
      if( !found )
      {
         elm = create_element() ;
         MAC_CHECK( elm != 0 ) ;
         elm->iCol = i ;
         elm->xVal = 0 ;
         if( dummy == 0 )
         {
            elm->next = ROW_TABLE[i] ;
            ROW_TABLE[i] = elm ;
         }
         else
         {
            elm->next = dummy->next ;
            dummy->next  = elm ;
         }

      }
      diag[i] = elm ;
      MAC_CHECK( diag[i]->iCol == i ) ;
   }

   if( MAC::abs( diag[0]->xVal )<piv_min )
   {
      nullify_row( 0 ) ;
      diag[0]->xVal = 1. ;
      raise_MILU0_zero_pivot( 0 ) ;
   }

   for( size_t i=1 ; i<n ; ++i )
   {
      RowElm* elm_ik = ROW_TABLE[i] ;
      double drop = 0.0 ;
      while( elm_ik!=0 && elm_ik->iCol<i )
      {
         size_t const k = elm_ik->iCol ;
         double const akk = diag[k]->xVal ;
         double const aik = elm_ik->xVal/akk ;
         elm_ik->xVal = aik ;
         RowElm* elm_kj = diag[k]->next ;
         while( elm_kj != 0 )
         {
            size_t const j = elm_kj->iCol ;
            MAC_CHECK( j>k ) ;
            double const akj = elm_kj->xVal ;
            bool found = false ;
            RowElm* elm_ij = ROW_TABLE[i] ;
            while( elm_ij != 0 && elm_ij->iCol<=j )
            {
               if( elm_ij->iCol == j )
               {
                  found = true ;
                  elm_ij->xVal -= aik*akj ;
                  break ;
               }
               elm_ij = elm_ij->next ;
            }
            if( !found ) drop -= aik*akj ;
            elm_kj = elm_kj->next ;
         }
         elm_ik = elm_ik->next ;
      }
      double aii = diag[i]->xVal ;
      if( modified )
      {
         aii -= drop ;
         diag[i]->xVal -= drop ;
      }
      if( MAC::abs( aii )<piv_min )
      {
         nullify_row( i ) ;
         diag[i]->xVal = 1. ;
         raise_MILU0_zero_pivot( i ) ;
      }
   }
   synchronize() ;

   delete [] diag ; diag = 0 ;

   MAC_CHECK_POST( factorize_MILU0_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: solve_LU( LA_SeqVector const* rhs,
                         LA_SeqVector* sol ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: solve_LU" ) ;
   MAC_CHECK_PRE( solve_LU_PRE( rhs, sol ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   size_t const n = NB_ROWS ;

   RowElm* elm = 0 ;
   double const* ptr_rhs = rhs->data() ;
   double* ptr_sol = sol->data() ;

   for( size_t i=0 ; i<n ; i++ )
   {
      elm = ROW_TABLE[i] ;
      double r = ptr_rhs[i] ;
      while( elm!=0 )
      {
         size_t j = elm->iCol ;
         if( j>= i ) break ;
         r -= elm->xVal * ptr_sol[j] ;
         elm = elm->next ;
      }
      ptr_sol[i] = r ;
   }

   for( size_t i=n-1 ; i<n ; i-- )
   {
      elm = ROW_TABLE[i] ;
      double r = ptr_sol[i] ;
      while( elm!=0 && elm->iCol<i )
      {
         elm = elm->next ;
      }
      MAC_CHECK( elm->iCol==i ) ;
      double d = elm->xVal ;
      elm = elm->next ;
      while( elm!=0 )
      {
         r -=  elm->xVal * ptr_sol[elm->iCol] ;
         elm = elm->next ;
      }
      ptr_sol[i] = r/d ;
   }

   sol->synchronize() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( solve_LU_POST( rhs, sol ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: relax( double omega,
                      LA_SeqMatrix::relaxation_mode mode,
                      LA_SeqVector const* omega_inv_diag,
                      LA_SeqVector const* rhs,
                      LA_SeqVector* sol ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: relax" ) ;
   MAC_CHECK_PRE( relax_PRE( omega, mode, omega_inv_diag, rhs, sol ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   size_t const n = NB_ROWS ;

   RowElm* elm = 0 ;

   double const* ptr_rhs = rhs->data() ;
   double const* ptr_omega_inv_diag = omega_inv_diag->data() ;
   double* ptr_sol = sol->data() ;

   if( mode == forward || mode == symmetric )
   {
      for( size_t i=0 ; i<n ; ++i )
      {
         elm = ROW_TABLE[i] ;
         double xt =  ptr_rhs[i] ;
         while( elm!=0 )
         {
            xt -= elm->xVal * ptr_sol[elm->iCol] ;
            elm = elm->next ;
         }
         ptr_sol[i] += xt * ptr_omega_inv_diag[i] ;
      }
   }

   if( mode == backward || mode == symmetric )
   {
      for( size_t i=n-1 ; i<n ; --i )
      {
         elm = ROW_TABLE[i] ;
         double xt =  ptr_rhs[i] ;
         while( elm!=0 )
         {
            xt -= elm->xVal * ptr_sol[elm->iCol] ;
            elm = elm->next ;
         }
         ptr_sol[i] += xt * ptr_omega_inv_diag[i] ;
      }
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( relax_POST( omega, mode, omega_inv_diag, rhs, sol ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: add_Mat_Mat_IMP( LA_PelMatrix const* A,
                                LA_PelMatrix const* B,
                                double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_PelMatrix:: add_Mat_Mat_IMP" ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;
   MAC_CHECK_PRE( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   start_local_modifs() ;
   for( size_t i=0 ; i<A->NB_ROWS ; i++ )
   {
      RowElm* elmA = A->ROW_TABLE[i] ;
      while( elmA!=0 )
      {
         size_t k = elmA->iCol ;
         RowElm* elmB = B->ROW_TABLE[k] ;
         while( elmB != 0 )
         {
            size_t j = elmB->iCol ;
            double x = alpha * elmA->xVal * elmB->xVal ;
            add_to_item( i, j, x ) ; //??? trop lent
            elmB = elmB->next ;
         }
         elmA = elmA->next ;
      }
   }
   stop_local_modifs() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_Mat_Mat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: add_Mat_tMat_IMP( LA_PelMatrix const* A,
                                 LA_PelMatrix const* B,
                                 double alpha )
//----------------------------------------------------------------------
// M <- alpha*A*d*transpose( B )+beta*M
{
   MAC_LABEL( "LA_PelMatrix:: add_Mat_tMat_IMP" ) ;
   MAC_CHECK_PRE( add_Mat_tMat_PRE( A, B, alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   start_local_modifs() ;
   for( size_t i=0 ; i<NB_ROWS ; i++ )
   {
      for( size_t j=0 ; j<NB_COLS ; j++ )
      {
         RowElm* elmA = A->ROW_TABLE[i] ;
         while( elmA!=0 )
         {
            size_t k = elmA->iCol ;
            RowElm* elmB = B->ROW_TABLE[j] ;
            while( elmB!=0 )
            {
               if( elmB->iCol==k )
               {
                  double x = alpha*elmA->xVal*elmB->xVal ;
                  add_to_item( i, j, x ) ;
                  break ;
               }
               elmB = elmB->next ;
            }
            elmA = elmA->next ;
         }
      }
   }
   stop_local_modifs() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( add_Mat_tMat_POST() ) ;
}

//----------------------------------------------------------------------
LA_PelMatrix::RowElm*
LA_PelMatrix:: create_element( void )
//----------------------------------------------------------------------
{
   if( CURRENT_BUCKET->used == 128 )
   {
      CURRENT_BUCKET->next = new RowElmBucket ;
      CURRENT_BUCKET=CURRENT_BUCKET->next ;
      CURRENT_BUCKET->used=0 ;
      CURRENT_BUCKET->next=0 ;
   }
   return &(CURRENT_BUCKET->elem[ CURRENT_BUCKET->used++ ] ) ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: destroy_all_elements( void )
//----------------------------------------------------------------------
{
   RowElmBucket * ptr = FIRST_BUCKET ;
   RowElmBucket * next ;

   while( ptr!=0 )
   {
      next = ptr->next ;
      delete ptr ;
      ptr = next ;
   }
   FIRST_BUCKET = 0 ;
   CURRENT_BUCKET = 0 ;
}

//----------------------------------------------------------------------
void
LA_PelMatrix:: init_allocating( void )
//----------------------------------------------------------------------
{
   MAC_ASSERT( FIRST_BUCKET==0 ) ;
   MAC_ASSERT( CURRENT_BUCKET==0 ) ;

   FIRST_BUCKET = CURRENT_BUCKET = new RowElmBucket ;
   CURRENT_BUCKET->next = 0 ;
   CURRENT_BUCKET->used = 0 ;
}

//----------------------------------------------------------------------
bool
LA_PelMatrix:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( EQUIVALENT( NB_ROWS==0, ROW_TABLE==0 ) ) ;
   MAC_ASSERT( is_a_prototype() || (FIRST_BUCKET!=0 && CURRENT_BUCKET!=0) ) ;
   return true ;
}
