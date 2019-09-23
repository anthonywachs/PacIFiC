#include <LA_Vector.hh>

#include <LA_SeqVector.hh>
#include <LA_Scatter.hh>

#include <MAC_Communicator.hh>
#include <MAC_DistributedPartition.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <size_t_vector.hh>

#ifdef OUTLINE
#define inline
   #include <LA_Vector.icc>
#undef inline
#endif

#include <iostream>
#include <fstream>


//----------------------------------------------------------------------
LA_Vector:: LA_Vector( MAC_Object* a_owner, size_t a_nb_rows )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , NB_ROWS( a_nb_rows )
   , DIST_STATUS( LA::Sync )
   , RESIZABLE( true )
   , ONLY_LOCAL_MODIFS( false )
   , DIST_STRAT( LA::InvalidDistribution )
   , VECNAME( MAC::undefined_string )
{
   MAC_LABEL( "LA_Vector:: LA_Vector" ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK_POST( owner() == a_owner ) ;
   MAC_CHECK_POST( nb_rows() == a_nb_rows ) ;
   MAC_CHECK_POST( is_resizable() ) ;
   MAC_CHECK_POST( distribution_strategy() == LA::InvalidDistribution ) ;
}




//----------------------------------------------------------------------
LA_Vector:: ~LA_Vector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: ~LA_Vector" ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: is_resizable( void ) const
//----------------------------------------------------------------------
{
   return( RESIZABLE ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: make_non_resizable( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: make_non_resizable" ) ;

   RESIZABLE = false ;

   MAC_CHECK_POST( !is_resizable() ) ;
}




//----------------------------------------------------------------------
size_t
LA_Vector:: nb_local_rows( void ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: nb_local_rows" ) ;
   
   size_t result = row_distribution()->local_number() ;
   
   MAC_CHECK_POST( result == row_distribution()->local_number() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: only_local_modifs( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: only_local_modifs" ) ;

   bool result = ONLY_LOCAL_MODIFS ;

   return( result ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: set_only_local_modifs_state( bool only_local )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: set_only_local_modifs_state" ) ;

   ONLY_LOCAL_MODIFS = only_local ;

   MAC_CHECK_POST( only_local_modifs() == only_local ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: start_local_modifs( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: start_local_modifs" ) ;
   MAC_CHECK_PRE( start_local_modifs_PRE() ) ;

   set_only_local_modifs_state( true ) ;

   MAC_CHECK_POST( start_local_modifs_POST() ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: stop_local_modifs( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: stop_local_modifs" ) ;
   MAC_CHECK_PRE( stop_local_modifs_PRE() ) ;

   set_only_local_modifs_state( false ) ;

   MAC_CHECK_POST( stop_local_modifs_POST() ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: set_distribution_strategy(
                              LA::DistributionStrategy dist_strat ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: set_distribution_strategy" ) ;
   MAC_CHECK_PRE( distribution_strategy() == LA::InvalidDistribution ) ;
   MAC_CHECK_PRE( dist_strat != LA::InvalidDistribution ) ;

   DIST_STRAT = dist_strat ;

   MAC_CHECK_POST( distribution_strategy() == dist_strat ) ;
}




//----------------------------------------------------------------------
LA::DistributionStrategy
LA_Vector:: distribution_strategy( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: distribution_strategy" ) ;

   return( DIST_STRAT ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: is_synchronized( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: is_synchronized" ) ;
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;

   bool result = ( is_desynchronizable() ?
                      MAC_Exec::communicator()->boolean_and(
                                                        state() == LA::Sync ) :
                      true ) ;

   MAC_CHECK_POST(
      IMPLIES( !is_desynchronizable(), result == true ) ) ;
   MAC_CHECK_POST(
      IMPLIES( is_desynchronizable(),
               result == MAC_Exec::communicator()->boolean_and(
                                               state() == LA::Sync ) ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: nullify( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: nullify" ) ;
   MAC_CHECK_PRE( nullify_PRE() ) ;

   set( 0.0 ) ;

   MAC_CHECK_POST( nullify_POST() ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string s( indent_width, ' ' ) ;
   os << s << "vector: \"" << type_name() << "\"" << std::endl ;
}




//----------------------------------------------------------------------
void
LA_Vector:: read( std::string const& filename )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: read" ) ;
   MAC_CHECK_PRE( read_PRE( filename ) ) ;

   if( !is_resizable() )
   {
      MAC_Error::object()->raise_internal(
         "Reading of non resizable vector is not implemented" ) ;
   }

   size_t rank = MAC_Exec::communicator()->rank() ;

   std::ifstream in( filename.c_str() ) ;
   if( !in )
   {
      MAC_Error::object()->raise_file_handling( filename, "open" ) ;
   }
   size_t nbrows ;

   in >> nbrows ;
   re_initialize( nbrows ) ;
   size_t i = 0 ;

   if( rank==0 || !is_desynchronizable() )
   {
      while( !in.eof() )
      {
         if( i>=nbrows )
         {
            MAC_Error::object()->raise_plain( "Invalid vector (rows>nb_rows) in "
                                              +filename ) ;
         }
         double v ;
         in >> v ;
         set_item( i, v ) ;
         i++ ;
      }
      if( i!=nbrows )
      {
         MAC_Error::object()->raise_plain( "Invalid vector (rows<nb_rows) in "
                                           +filename ) ;
      }
   }
   in.close() ;
   synchronize() ;

   MAC_CHECK_POST( read_POST() ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: set_rows_number( size_t a_nb_rows )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: set_rows_number" ) ;
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;

   NB_ROWS = a_nb_rows ;

   MAC_CHECK_POST( nb_rows() == a_nb_rows ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: re_initialize_PRE( size_t a_nb_rows,
                               size_t a_nb_local_rows ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( is_resizable() ) ;

   MAC_ASSERT( IMPLIES( distribution_strategy() != LA::FromLocalSize,
                        a_nb_rows != MAC::bad_index() ) ) ;
   MAC_ASSERT(
     IMPLIES( distribution_strategy() != LA::FromLocalSize && is_desynchronizable(),
       MAC_Exec::communicator()->same_value_everywhere( (int) a_nb_rows ) ) ) ;

   MAC_ASSERT( IMPLIES( distribution_strategy() == LA::FromLocalSize,
                           a_nb_local_rows != MAC::bad_index() ) ) ;

   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: re_initialize_POST( size_t a_nb_rows,
                                size_t a_nb_local_rows ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( IMPLIES( distribution_strategy() != LA::FromLocalSize,
                        nb_rows() == a_nb_rows ) ) ;
   MAC_ASSERT(
       IMPLIES( distribution_strategy() == LA::FromLocalSize,
                row_distribution()->local_number() ==  a_nb_local_rows ) ) ;

   MAC_ASSERT( state() == LA::Sync ) ;
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: create_vector_PRE( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: create_vector_POST(
                          LA_Vector* result, MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   MAC_ASSERT( result->implementation() == implementation() ) ;
   MAC_ASSERT( result->is_resizable() == is_resizable() ) ;
   MAC_ASSERT( result->nb_rows() == nb_rows() ) ;
   MAC_ASSERT( result->is_desynchronizable() == is_desynchronizable() ) ;
   MAC_ASSERT( result->state() == LA::Sync ) ;
   MAC_ASSERT( result->row_distribution()->is_compatible( row_distribution() ) ) ;
   MAC_ASSERT( result->is_synchronized() ) ;
   MAC_ASSERT( result->distribution_strategy() == distribution_strategy() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: implementation_POST( LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: start_local_modifs_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( !only_local_modifs() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: start_local_modifs_POST( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( only_local_modifs() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: stop_local_modifs_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( only_local_modifs() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: stop_local_modifs_POST( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( !only_local_modifs() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: synchronize_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: synchronize_POST( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: row_distribution_PRE( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: row_distribution_POST(
                          MAC_DistributedPartition const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == this ) ;
   MAC_ASSERT( result->partitioning().size()==
                               MAC_Exec::communicator()->nb_ranks() ) ;
   MAC_ASSERT( result->global_number()==nb_rows() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: create_scatter_PRE( MAC_Object* a_owner,
                                size_t_vector const& from,
                                size_t_vector const& to ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( from.size() == to.size() ) ;
   MAC_ASSERT( FORALL( ( size_t i=0 ; i<from.size() ; ++i ),
                       from(i) < nb_rows() ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: create_scatter_POST( LA_Scatter const* result,
                                 MAC_Object const* a_owner,
                                 size_t_vector const& from,
                                 size_t_vector const& to ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   MAC_ASSERT( result->implementation() == implementation() ) ;
   MAC_ASSERT( result->size() == from.size() ) ;
   MAC_ASSERT(
      FORALL( ( size_t i=0 ; i<result->size() ; ++i ),
              result->repatriated_items().has( from(i) ) ) ) ;
   MAC_ASSERT(
      FORALL( ( size_t i=0 ; i<result->size() ; ++i ),
              result->local_indices()(i) == to(
                 from.index_of( result->repatriated_items()(i) ) ) ) ) ;
   MAC_ASSERT( result->distribution()->is_compatible(
                                                row_distribution() ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: create_local_vector_PRE( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: create_local_vector_POST(
                       LA_SeqVector* result, MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   MAC_ASSERT( result->nb_rows() == nb_rows() ) ;
   MAC_ASSERT( result->is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: item_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( i < row_distribution()->local_index_limit()  &&
               i >= row_distribution()->first_local_index() ) ;
   MAC_ASSERT( state() == LA::Sync ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: set_item_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( i < nb_rows() ) ;
   MAC_ASSERT( IMPLIES( only_local_modifs(),
                        i >= row_distribution()->first_local_index() &&
                        i < row_distribution()->local_index_limit() ) ) ;
   MAC_ASSERT( state() != LA::NotSync_add ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: set_item_POST( size_t i, LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT(
   IMPLIES( is_desynchronizable() && !only_local_modifs(),
            state() == LA::NotSync_set ) ) ;
   MAC_ASSERT(
   IMPLIES( !is_desynchronizable() || only_local_modifs(),
            state() == old_state ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: add_to_item_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( i < nb_rows() ) ;
   MAC_ASSERT( IMPLIES( only_local_modifs(),
                        i >= row_distribution()->first_local_index() &&
                        i < row_distribution()->local_index_limit() ) ) ;
   MAC_ASSERT( state() != LA::NotSync_set ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: add_to_item_POST( size_t i, LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT(
      IMPLIES( is_desynchronizable() && !only_local_modifs(),
               state() == LA::NotSync_add ) ) ;
   MAC_ASSERT(
   IMPLIES( !is_desynchronizable() || only_local_modifs(),
            state() == old_state ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: dot_PRE( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( is_synchronized() ) ;
   MAC_ASSERT( a != 0 ) ;
   MAC_ASSERT( a->nb_rows() == nb_rows() ) ;
   MAC_ASSERT( a->implementation() == implementation() ) ;
   MAC_ASSERT( a->is_synchronized() ) ;
   MAC_ASSERT( a->row_distribution()->is_compatible( row_distribution() ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: dot_POST( double result ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: two_norm_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( nb_rows() > 0 ) ;
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: two_norm_POST( double result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result >= 0.0 ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: max_norm_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( nb_rows() > 0 ) ;
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: max_norm_POST( double result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result >= 0.0 ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: nullify_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: nullify_POST( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: scale_PRE( double value ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( is_synchronized() ) ;
   MAC_ASSERT(
      IMPLIES( is_desynchronizable(),
               MAC_Exec::communicator()->same_value_everywhere( value ) ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: scale_POST( double value ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: set_PRE( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( a != 0 ) ;
   MAC_ASSERT( a->nb_rows() == nb_rows() ) ;
   MAC_ASSERT( a->implementation() == implementation() ) ;
   MAC_ASSERT( a->is_synchronized() ) ;
   MAC_ASSERT( a->row_distribution()->is_compatible( row_distribution() ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: set_POST( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: set_PRE( double value ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT(
      IMPLIES( is_desynchronizable(),
               MAC_Exec::communicator()->same_value_everywhere( value ) ) )  ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: set_POST( double value ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: set_as_v_product_PRE(
                          LA_Vector const* a, LA_Vector const* b ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( a != 0 ) ;
   MAC_ASSERT( a->nb_rows() == nb_rows() ) ;
   MAC_ASSERT( a->implementation() == implementation() ) ;
   MAC_ASSERT( a->is_synchronized() ) ;
   MAC_ASSERT( a->row_distribution()->is_compatible( row_distribution() ) ) ;
   MAC_ASSERT( b != 0 ) ;
   MAC_ASSERT( b->nb_rows() == nb_rows() ) ;
   MAC_ASSERT( b->implementation() == implementation() ) ;
   MAC_ASSERT( b->is_synchronized() ) ;
   MAC_ASSERT( b->row_distribution()->is_compatible( row_distribution() ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: set_as_v_product_POST(
                          LA_Vector const* a, LA_Vector const* b ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_synchronized() ) ;
   MAC_ASSERT( FORMAL( FORALL( ( size_t i=0 ; i<nb_rows() ; ++i ),
                               item(i) == a->item(i)*b->item(i) ) ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: set_as_reciprocal_PRE( LA_Vector const* a,
                                   double smallest_inverted_item,
                                   double default_value ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( a != 0 ) ;
   MAC_ASSERT( a->nb_rows() == nb_rows() ) ;
   MAC_ASSERT( a->implementation() == implementation() ) ;
   MAC_ASSERT( a->is_synchronized() ) ;
   MAC_ASSERT( a->row_distribution()->is_compatible( row_distribution() ) ) ;
   MAC_ASSERT( smallest_inverted_item >= 0. ) ;
   MAC_ASSERT(
      IMPLIES( is_desynchronizable(),
               MAC_Exec::communicator()->same_value_everywhere( smallest_inverted_item ) ) ) ;
   MAC_ASSERT(
      IMPLIES( is_desynchronizable(),
               MAC_Exec::communicator()->same_value_everywhere( default_value ) ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: set_as_reciprocal_POST( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: sum_PRE( LA_Vector const* a, double alpha ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( is_synchronized() ) ;
   MAC_ASSERT( a != 0 ) ;
   MAC_ASSERT( a->nb_rows() == nb_rows() ) ;
   MAC_ASSERT( a->implementation() == implementation() ) ;
   MAC_ASSERT( a->is_synchronized() ) ;
   MAC_ASSERT( a->row_distribution()->is_compatible( row_distribution() ) ) ;
   MAC_ASSERT(
      IMPLIES( is_desynchronizable(),
               MAC_Exec::communicator()->same_value_everywhere( alpha ) ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: sum_POST(  LA_Vector const* a, double alpha ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: print_items_PRE( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
//   MAC_ASSERT( os ) ; // Not accepted from gcc-9.x.x
   MAC_ASSERT( os.good() ) ;
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: write_PRE( std::string const& filename ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( is_synchronized() ) ;
   MAC_ASSERT( !filename.empty() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: write_POST( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: read_PRE( std::string const& filename ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   MAC_ASSERT( !filename.empty() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_Vector:: read_POST( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: save_state( MAC_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: save_state" ) ;
   MAC_CHECK_PRE( save_state_PRE( writer ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Error::object()->raise_not_implemented( this, "save_state" ) ;

   MAC_CHECK_POST( save_state_POST( writer ) ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: restore_state( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: restore_state" ) ;
   MAC_CHECK_PRE( restore_state_PRE( reader ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   MAC_Error::object()->raise_not_implemented( this, "restore_state" ) ;
   
   MAC_CHECK_POST( restore_state_POST( reader ) ) ;
}




//----------------------------------------------------------------------
void
LA_Vector:: read_state_nonrestored( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: restore_state" ) ;

   reader->start_object_retrieval( "LA_Vector" ) ;
   reader->end_object_retrieval() ;
}




//----------------------------------------------------------------------
std::string const&
LA_Vector:: name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: name" ) ;
   MAC_CHECK_INV( invariant() ) ;

   return( VECNAME );
}




//----------------------------------------------------------------------
void
LA_Vector:: set_name( std::string const& name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Vector:: set_name" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   VECNAME = name ;
}
