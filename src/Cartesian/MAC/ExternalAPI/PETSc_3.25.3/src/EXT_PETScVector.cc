#include <EXT_PETScVector.hh>

#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Int.hh>
#include <MAC_String.hh>
#include <EXT_PETScImplementation.hh>
#include <MAC_DoubleVector.hh>
#include <doubleVector.hh>

#include <iomanip>
#include <fstream>


//----------------------------------------------------------------------
EXT_PETScVector*
EXT_PETScVector:: create( MAC_Object* a_owner,
                          bool sequential,
                          size_t a_nb_rows )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: create" ) ;

   EXT_PETScVector* result = new EXT_PETScVector( a_owner, sequential,
                                                  a_nb_rows ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->state() == LA::Sync ) ;
   MAC_CHECK_POST( result->is_desynchronizable() ) ;
   MAC_CHECK_POST( result->is_resizable() ) ;
   MAC_CHECK_POST( result->is_synchronized() ) ;
   MAC_CHECK_POST( result->distribution_strategy() !=
                   LA::InvalidDistribution ) ;
   MAC_CHECK_POST( 
       EQUIVALENT( result->distribution_strategy() == LA::NoDistribution, 
                   sequential ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
EXT_PETScVector:: EXT_PETScVector( MAC_Object* a_owner,
                                   bool sequential,
                                   size_t a_nb_rows )
//----------------------------------------------------------------------
   : LA_Vector( a_owner, 0 )
   , VECTOR( 0 )
   , DIST( MAC_DistributedPartition::create( this ) )
   , SEQ( sequential )
{
   MAC_LABEL( "EXT_PETScVector:: EXT_PETScVector" ) ;

   if( sequential )
   {
      set_distribution_strategy( LA::NoDistribution ) ;
   }
   else
   {
      set_distribution_strategy( LA::FromGlobalSize ) ;
      //la strategie LA::FromLocalSize est disponible dans PETSc
      //mais pas implemente ici
   }
   re_initialize( a_nb_rows ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( 
       EQUIVALENT( distribution_strategy() == LA::NoDistribution, 
                   sequential ) ) ;
   MAC_CHECK_POST( nb_rows() == a_nb_rows ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: re_initialize( size_t a_nb_rows, size_t a_nb_local_rows )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: re_initialize" ) ;
   MAC_CHECK_PRE( re_initialize_PRE( a_nb_rows, a_nb_local_rows ) ) ;

   int nloc = PETSC_DECIDE ;
   size_t nglob = a_nb_rows ;

   if( VECTOR == 0 || a_nb_rows != nb_rows() )
   {
      if( VECTOR!=0 ) PETSc_do( VecDestroy( &VECTOR ) ) ;

      if( SEQ )
      {
         PETSc_do( VecCreateSeq( PETSC_COMM_SELF, PetscInt(a_nb_rows), 
	 	&VECTOR ) ) ;
      }
      else
      {
         PETSc_do( VecCreateMPI( PETSC_COMM_WORLD, nloc, PetscInt(nglob), 
	 	&VECTOR ) ) ;
      }
      int s ;
      PETSc_do( VecGetLocalSize( VECTOR, &s ) ) ;
      DIST->set_local_number( s ) ;
      int i0, i1 ;
      PETSc_do( VecGetOwnershipRange( VECTOR, &i0, &i1 ) ) ;
      FIRST = i0 ;
      LAST = i1 ;
      set_rows_number( a_nb_rows ) ;
   }
   set( 0.0 ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( re_initialize_POST( a_nb_rows, a_nb_local_rows ) ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=row_distribution()->first_local_index() ;
                             i<row_distribution()->local_index_limit() ; ++i ),
                             item( i ) == 0.0 ) ) ;
}




//----------------------------------------------------------------------
EXT_PETScVector*
EXT_PETScVector:: create_vector( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: create_vector" ) ;
   MAC_CHECK_PRE( create_vector_PRE( a_owner ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   EXT_PETScVector* result =
      new EXT_PETScVector( a_owner, SEQ, nb_rows() ) ;

   MAC_CHECK_POST( create_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
EXT_PETScScatter*
EXT_PETScVector:: create_scatter( MAC_Object* a_owner,
                                  size_t_vector const& from,
                                  size_t_vector const& to ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: create_scatter" ) ;
   MAC_CHECK_PRE( create_scatter_PRE( a_owner, from, to ) ) ;

   MAC_CHECK_INV( invariant() ) ;

   EXT_PETScScatter* result =
      EXT_PETScScatter::create( a_owner, this, from, to ) ;

   MAC_CHECK_POST( create_scatter_POST( result, a_owner, from, to ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
LA_SeqVector*
EXT_PETScVector:: create_local_vector( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: create_local_vector" ) ;

   MAC_CHECK_PRE( create_local_vector_PRE( a_owner ) ) ;

   LA_SeqVector* result = LA_SeqVector::create( a_owner, nb_rows() ) ;
   PetscScalar* values = 0 ;
   PETSc_do( VecGetArray( VECTOR, &values ) ) ;

   for( size_t i=FIRST ; i<LAST; i++ )
   {
      result->set_item(i, (double) values[i-FIRST]  ) ;
   }
   PETSc_do( VecRestoreArray( VECTOR, &values ) ) ;
   result->synchronize() ;

   MAC_CHECK_POST( create_local_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
double
EXT_PETScVector:: item( size_t i ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: item" ) ;
   MAC_CHECK_PRE( item_PRE( i ) ) ;

   PetscScalar values = PetscScalar(MAC::bad_index()) ;
   PetscInt ni = 1 ;
   PetscInt ix = PetscInt(i) ;
   PETSc_do( VecGetValues( VECTOR, ni, &ix, &values ) ) ;

   return( values ) ;
}




//----------------------------------------------------------------------------
LA_Implementation const*
EXT_PETScVector:: implementation( void ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: implementation" ) ;

   LA_Implementation const* result = EXT_PETScImplementation::object() ;

   MAC_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
EXT_PETScVector:: ~EXT_PETScVector( void )
//----------------------------------------------------------------------
{
   if( VECTOR!=0 ) PETSc_do( VecDestroy( &VECTOR ) ) ;
}




//----------------------------------------------------------------------
bool
EXT_PETScVector:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: is_desynchronizable" ) ;

   bool result = true ;

   MAC_CHECK_POST( result == true ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_DistributedPartition const*
EXT_PETScVector:: row_distribution( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: row_distribution" ) ;
   MAC_CHECK_PRE( row_distribution_PRE() ) ;

   MAC_DistributedPartition const* result = DIST ;

   MAC_CHECK_POST( row_distribution_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: start_local_modifs( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: start_local_modifs" ) ;
   MAC_CHECK_PRE( start_local_modifs_PRE() ) ;

   std::string mesg ;
   mesg += "*** EXT_PETScVector:: start_local_modifs" "\n" ;
   mesg += "    PETSc does not allow this functionality." ;
   MAC_Error::object()->raise_internal( mesg ) ;

   MAC_CHECK_POST( start_local_modifs_POST() ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: stop_local_modifs( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: stop_local_modifs" ) ;
   MAC_CHECK_PRE( stop_local_modifs_PRE() ) ;

   std::string mesg ;
   mesg += "*** EXT_PETScVector:: stop_local_modifs" "\n" ;
   mesg += "    PETSc does not allow this functionality." ;
   MAC_Error::object()->raise_internal( mesg ) ;

   MAC_CHECK_POST( stop_local_modifs_POST() ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: set( double value )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: set double" ) ;
   MAC_CHECK_PRE( set_PRE( value ) ) ;

   PETSc_do( VecSet( vector(), value ) ) ;
   synchronize() ;

   MAC_CHECK_POST( set_POST( value ) ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: set( LA_Vector const* a )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: set LA_Vector" ) ;
   MAC_CHECK_PRE( set_PRE( a ) ) ;
   MAC_CHECK( dynamic_cast<EXT_PETScVector const*>( a ) !=0 ) ;
   EXT_PETScVector const* pa = static_cast<EXT_PETScVector const*>( a ) ;
   MAC_CHECK_INV( invariant() ) ;

   PETSc_do( VecCopy( pa->vector(), vector() ) ) ;
   synchronize() ;

   MAC_CHECK_POST( set_POST( a ) ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: scale( double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: scale" ) ;
   MAC_CHECK_PRE( scale_PRE( alpha ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   PETSc_do( VecScale( VECTOR, alpha ) ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( scale_POST( alpha ) ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: sum( LA_Vector const* a, double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: sum LA_Vector" ) ;
   MAC_CHECK_PRE( sum_PRE( a, alpha ) ) ;
   MAC_CHECK( dynamic_cast<EXT_PETScVector const*>( a ) !=0 ) ;
   EXT_PETScVector const* pa = static_cast<EXT_PETScVector const*>( a ) ;

   PETSc_do( VecAXPY( VECTOR, alpha, pa->vector() ) ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( sum_POST( a, alpha ) ) ;

}




//----------------------------------------------------------------------
double
EXT_PETScVector:: dot( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: dot LA_Vector" ) ;
   MAC_CHECK_PRE( dot_PRE( a ) ) ;
   MAC_CHECK( dynamic_cast<EXT_PETScVector const*>( a ) !=0 ) ;
   EXT_PETScVector const* pa = static_cast<EXT_PETScVector const*>( a ) ;
   double result ;
   PETSc_do( VecDot( VECTOR, pa->vector(), &result ) ) ;

   MAC_CHECK_POST( dot_POST( result ) ) ;
   return result ;
}




//----------------------------------------------------------------------
double
EXT_PETScVector:: two_norm( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: two_norm" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( two_norm_PRE() ) ;

   double result ;
   PETSc_do( VecNorm( VECTOR, NORM_2, &result ) ) ;

   MAC_CHECK_POST( two_norm_POST( result ) ) ;

   return( result ) ;
}




//----------------------------------------------------------------------
double
EXT_PETScVector:: max_norm( void ) const
//----------------------------------------------------------------------
// max_norm = max( |v( i )| )
{
   MAC_LABEL( "EXT_PETScVector:: max_norm" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( max_norm_PRE() ) ;

   double result ;
   PETSc_do( VecNorm( VECTOR, NORM_INFINITY, &result ) ) ;

   MAC_CHECK_POST( max_norm_POST( result ) ) ;

   return( result ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: synchronize( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: synchronize" ) ;
   MAC_CHECK_INV( invariant() ) ;

   PETSc_do( VecAssemblyBegin( VECTOR ) ) ;
   PETSc_do( VecAssemblyEnd( VECTOR ) ) ;

   LA_Vector::synchronize() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( is_synchronized() ) ;
}




//----------------------------------------------------------------------
Vec&
EXT_PETScVector:: vector( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: vector" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK( VECTOR!=0 ) ;

   return const_cast<EXT_PETScVector*>(this)->VECTOR ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: print_items( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: print_items" ) ;
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
      double * values ;
      PETSc_do( VecGetArray( VECTOR, &values ) ) ;
      for( size_t iRow = FIRST ; iRow<LAST ; ++iRow )
      {
         os << space << "Row n°" << iRow << "  " ;
         double const x = values[iRow-FIRST] ;
         os << std::setw(15) << x << std::endl ;
      }
      os << std::setprecision( int(p) ) ;
      os.flags( original_flags ) ;
      PETSc_do( VecRestoreArray( VECTOR, &values ) ) ;
   }
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: set_as_reciprocal( LA_Vector const* a,
                                     double smallest_inverted_item,
                                     double default_value )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: set_as_reciprocal LA_Vector" ) ;
   MAC_CHECK_PRE( set_as_reciprocal_PRE(a,smallest_inverted_item,
   	default_value) ) ;

   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( dynamic_cast<EXT_PETScVector const*>( a ) !=0 ) ;
   EXT_PETScVector const* pa = static_cast<EXT_PETScVector const*>( a ) ;

   PETSc_do( VecCopy( pa->vector(), VECTOR ) ) ;

   if( smallest_inverted_item!=0.0 )
   {
      double * values ;
      PETSc_do( VecGetArray( pa->vector(), &values ) ) ;
      for( size_t iRow = FIRST ; iRow<LAST ; ++iRow )
      {
         if( MAC::abs( values[iRow-FIRST] )<smallest_inverted_item )
         {
            VecSetValue( VECTOR, PetscInt(iRow), 1.0/default_value, 
	    	INSERT_VALUES ) ;
         }
      }
      PETSc_do( VecRestoreArray( VECTOR, &values ) ) ;
   }

   PETSc_do( VecReciprocal( VECTOR ) ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_as_reciprocal_POST(a) ) ;

}




//----------------------------------------------------------------------
void
EXT_PETScVector:: set_as_v_product( LA_Vector const* a,
                                         LA_Vector const* b )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: set_as_v_product LA_Vector" ) ;
   MAC_CHECK_PRE( set_as_v_product_PRE(a,b) ) ;

   MAC_CHECK( dynamic_cast<EXT_PETScVector const*>( a ) !=0 ) ;
   MAC_CHECK( dynamic_cast<EXT_PETScVector const*>( b ) !=0 ) ;
   EXT_PETScVector const* pa = static_cast<EXT_PETScVector const*>( a ) ;
   EXT_PETScVector const* pb = static_cast<EXT_PETScVector const*>( b ) ;

   PETSc_do( VecPointwiseMult( VECTOR, pa->vector(), pb->vector() ) ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_as_v_product_POST( a, b ) ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: set_item( size_t i, double x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: set_item" ) ;
   MAC_CHECK_PRE( set_item_PRE( i ) ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;

   if( !only_local_modifs() && is_desynchronizable() )
   {
      set_unsynchronized_state( LA::NotSync_set ) ;
   }
   VecSetValue( VECTOR, PetscInt(i), x, INSERT_VALUES ) ;

   MAC_CHECK_POST( set_item_POST( i, OLD( state ) ) ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: add_to_item( size_t i, double x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: add_to_item" ) ;
   MAC_CHECK_PRE( add_to_item_PRE( i ) ) ;
   MAC_SAVEOLD( LA::SyncState, state, state() ) ;

   if( !only_local_modifs() && is_desynchronizable() )
   {
      set_unsynchronized_state( LA::NotSync_add ) ;
   }
   if( x != 0.0 )
   {
      VecSetValue( VECTOR, PetscInt(i), x, ADD_VALUES ) ;
   }
   MAC_CHECK_POST( add_to_item_POST( i, OLD( state ) ) ) ;
}




//----------------------------------------------------------------------
bool
EXT_PETScVector:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Vector::invariant() ) ;
   MAC_ASSERT( distribution_strategy() != LA::InvalidDistribution ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
EXT_PETScVector:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Vector::implementation_POST( result ) ) ;
   MAC_ASSERT( result == EXT_PETScImplementation::object() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: write( std::string const& filename ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: write" ) ;
   MAC_CHECK_PRE( is_synchronized() ) ;

   size_t rank = MAC_Exec::communicator()->rank() ;
   size_t size = MAC_Exec::communicator()->nb_ranks() ;
   int dummy ;

   if( rank>0 ) MAC_Exec::communicator()->receive( rank-1, dummy ) ;

   std::ofstream out( filename.c_str(), 
   	(rank==0 ? std::ios::out : std::ios::out|std::ios::app) ) ;
   if( rank==0 )
      out << nb_rows() ;

   double * values ;
   PETSc_do( VecGetArray( VECTOR, &values ) ) ;
   for( size_t iRow = FIRST ; iRow<LAST ; ++iRow )
   {
      out << std::endl << values[iRow-FIRST]  ;
   }
   PETSc_do( VecRestoreArray( VECTOR, &values ) ) ;
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
void
EXT_PETScVector:: save_state( MAC_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: save_state" ) ;
   MAC_CHECK_PRE( save_state_PRE( writer ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   writer->start_new_object( "LA_Vector" ) ;
   
   writer->add_entry( "name", MAC_String::create( 0, name() ) ) ;
   writer->add_entry( "first_global_index", MAC_Int::create( 0, int(FIRST) ) ) ;
   writer->add_entry( "last_global_index", MAC_Int::create( 0, int(LAST) ) ) ;

   // Saving values
   doubleVector localvalues( LAST - FIRST );
   
   double * values ;
   PETSc_do( VecGetArray( VECTOR, &values ) ) ;
   for( size_t iRow = FIRST ; iRow<LAST ; ++iRow )
     localvalues(iRow-FIRST) = values[iRow-FIRST];
   PETSc_do( VecRestoreArray( VECTOR, &values ) ) ;
   
   writer->add_entry( "values", MAC_DoubleVector::create( 0, localvalues ) ) ;
   
   writer->finalize_object() ;

   MAC_CHECK_POST( save_state_POST( writer ) ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScVector:: restore_state( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScVector:: restore_state" ) ;
   MAC_CHECK_PRE( restore_state_PRE( reader ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   reader->start_object_retrieval( "LA_Vector" ) ;

   // Retrieving stored datas :
   size_t const iFIRST = 
   	size_t(reader->data_of_entry( "first_global_index" )->to_int()) ;
   size_t const iLAST = 
   	size_t(reader->data_of_entry( "last_global_index" )->to_int()) ;

   // Does some checks
   MAC_ASSERT( FIRST==iFIRST ) ;
   MAC_ASSERT( LAST==iLAST ) ;

   // Retrieving values
   doubleVector localvalues 
   		= reader->data_of_entry( "values" )->to_double_vector() ;
   MAC_ASSERT( localvalues.size() == LAST - FIRST ) ;     

   // Set values in PETSC vector
   for( size_t iRow = FIRST ; iRow<LAST ; ++iRow )
     set_item( iRow, localvalues(iRow-FIRST) );       

   reader->end_object_retrieval() ;
   
   synchronize();
   
   MAC_CHECK_POST( restore_state_POST( reader ) ) ;
}
