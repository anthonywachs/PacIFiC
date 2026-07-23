#include <EXT_PETScScatter.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_DistributedPartition.hh>

#include <LA_SeqVector.hh>

#include <EXT_PETScImplementation.hh>
#include <EXT_PETScVector.hh>


//----------------------------------------------------------------------
EXT_PETScScatter*
EXT_PETScScatter:: create( MAC_Object* a_owner,
                           EXT_PETScVector const* global_vector,
                           size_t_vector const& a_repatriated_items_table,
                           size_t_vector const& a_local_indices_table ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScScatter:: create" ) ;
   MAC_CHECK_PRE(
      FORALL(
         ( size_t i=0 ; i<a_repatriated_items_table.size() ; ++i ),
         a_repatriated_items_table(i)<global_vector->nb_rows() ) ) ;
   MAC_CHECK_PRE(
      a_repatriated_items_table.size()==a_local_indices_table.size() ) ;
   
   EXT_PETScScatter* result = new EXT_PETScScatter( a_owner,
                                                    global_vector,
                                                    a_repatriated_items_table,
                                                    a_local_indices_table ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST(
      IMPLIES( global_vector->is_desynchronizable(),
               result->distribution() != 0 ) ) ;
   MAC_CHECK_POST(
      IMPLIES( result->distribution() != 0,
               result->distribution()->is_compatible(
                               global_vector->row_distribution() ) ) ) ;
   MAC_CHECK_POST( result->size() == a_repatriated_items_table.size() ) ;
   MAC_CHECK_POST( result->repatriated_items() ==
                                            a_repatriated_items_table ) ;
   MAC_CHECK_POST( result->local_indices() == a_local_indices_table ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
EXT_PETScScatter:: EXT_PETScScatter(
                           MAC_Object* a_owner,
                           EXT_PETScVector const* global_vector,
                           size_t_vector const& a_repatriated_items_table,
                           size_t_vector const& a_local_indices_table )
//----------------------------------------------------------------------
   : LA_Scatter( a_owner )
   , SIZE( a_repatriated_items_table.size() )
   , DIST( MAC_DistributedPartition::create( this ) )
   , NEEDED( a_repatriated_items_table )
   , LOCAL( a_local_indices_table )
{
   MAC_LABEL( "EXT_PETScScatter:: EXT_PETScScatter" ) ;

   DIST->set( global_vector->row_distribution() ) ;
   
   size_t s = 0 ;
   for( size_t i=0 ; i<SIZE ; ++i ) if( s<=LOCAL(i) ) s = LOCAL(i) ;
   
   int* idx_from = new int [ SIZE ] ;
   int* idx_to = new int [ SIZE ] ;
   for( size_t i=0 ; i<SIZE ; ++i )
   {
      idx_from[i] = NEEDED(i) ;
      idx_to[i] = LOCAL(i) ;
   }
   IS from, to ;
   PETSc_do( ISCreateGeneral( PETSC_COMM_SELF, SIZE, idx_from,
   	PETSC_COPY_VALUES, &from ) ) ;
   PETSc_do( ISCreateGeneral( PETSC_COMM_SELF, SIZE, idx_to, 
   	PETSC_COPY_VALUES, &to ) ) ;
   
   PETSc_do( VecCreateSeq( PETSC_COMM_SELF, s+1, &SEQ ) ) ;
   PETSc_do( VecScatterCreate(
                global_vector->vector(), from, SEQ, to, &SCATTER ) ) ;
   
   delete [] idx_from ;
   delete [] idx_to ;
   PETSc_do( ISDestroy( &from ) ) ;
   PETSc_do( ISDestroy( &to ) ) ;
}




//----------------------------------------------------------------------
EXT_PETScScatter:: ~EXT_PETScScatter( void )
//----------------------------------------------------------------------
{
   PETSc_do( VecScatterDestroy( &SCATTER ) ) ;
   PETSc_do( VecDestroy( &SEQ ) ) ;
}




//----------------------------------------------------------------------
LA_Implementation const*
EXT_PETScScatter:: implementation( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScScatter:: implementation" ) ;
   
   static LA_Implementation const* result = EXT_PETScImplementation::object() ;

   MAC_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
size_t
EXT_PETScScatter:: size( void ) const
//----------------------------------------------------------------------
{
   return( SIZE ) ;
}




//----------------------------------------------------------------------
size_t_vector const&
EXT_PETScScatter:: repatriated_items( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScScatter:: repatriated_items" ) ;
   
   size_t_vector const& result = NEEDED ;

   MAC_CHECK_POST( repatriated_items_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
size_t_vector const&
EXT_PETScScatter:: local_indices( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScScatter:: local_indices" ) ;
   
   size_t_vector const& result = LOCAL ;

   MAC_CHECK_POST( local_indices_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_DistributedPartition const*
EXT_PETScScatter:: distribution( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScScatter:: distribution" ) ;

   MAC_DistributedPartition const* result = DIST ;

   MAC_CHECK_POST( distribution_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScScatter:: get( LA_Vector const* source,
                        LA_SeqVector* dest ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScScatter:: get" ) ;
   MAC_CHECK_PRE( get_PRE( source, dest) ) ;

   EXT_PETScVector const* psource =
                        static_cast<EXT_PETScVector const*>( source ) ;
   MAC_CHECK( dynamic_cast<EXT_PETScVector const*>( source ) != 0 ) ;

   PETSc_do( VecScatterBegin( SCATTER, psource->vector(),
                              SEQ, INSERT_VALUES, SCATTER_FORWARD ) ) ;
   PETSc_do( VecScatterEnd( SCATTER, psource->vector(),
                            SEQ, INSERT_VALUES, SCATTER_FORWARD ) ) ;
   PetscScalar *ptr_x ;

   PETSc_do( VecGetArray( SEQ, &ptr_x ) ) ;
   for( size_t i=0 ; i<SIZE ; i++ )
   {
      dest->set_item( LOCAL(i), ptr_x[LOCAL(i)] ) ;
   }
   PETSc_do( VecRestoreArray(SEQ,&ptr_x) ) ;
   dest->synchronize() ;
   
   MAC_CHECK_POST( get_POST( source, dest) ) ;
}




//----------------------------------------------------------------------
void
EXT_PETScScatter:: set( LA_SeqVector const* source,
                        LA_Vector* dest ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScScatter:: set" ) ;
   MAC_CHECK_PRE( set_PRE( source, dest) ) ;

   MAC_CHECK( dynamic_cast<EXT_PETScVector const*>( dest ) != 0 ) ;
   EXT_PETScVector const* pdest = static_cast<EXT_PETScVector const*>( dest ) ;

   for( size_t i=0 ; i<SIZE ; i++ )
   {
      VecSetValue( SEQ, LOCAL(i), source->item(LOCAL(i)), INSERT_VALUES ) ;
   }
   PETSc_do( VecScatterBegin( SCATTER, SEQ, pdest->vector(),
                              INSERT_VALUES, SCATTER_REVERSE ) ) ;
   PETSc_do( VecScatterEnd( SCATTER, SEQ, pdest->vector(),
                            INSERT_VALUES, SCATTER_REVERSE ) ) ;
   static_cast<LA_Vector*>(dest)->synchronize() ;
   
   MAC_CHECK_POST( set_POST( source, dest) ) ;
}




//----------------------------------------------------------------------
bool
EXT_PETScScatter:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Scatter::implementation_POST( result ) ) ;
   MAC_ASSERT( result == EXT_PETScImplementation::object() ) ;
   return( true ) ;
}
