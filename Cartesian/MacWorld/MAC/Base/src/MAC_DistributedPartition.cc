#include <MAC_DistributedPartition.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>

#include <doubleVector.hh>

#ifdef OUTLINE
   #define inline
   #include <MAC_DistributedPartition.icc>
   #undef inline
#endif


//----------------------------------------------------------------------
MAC_DistributedPartition*
MAC_DistributedPartition:: create( MAC_Object* a_owner )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DistributedPartition:: create" ) ;

   MAC_DistributedPartition* result =
                               new MAC_DistributedPartition( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->communicator() == MAC_Exec::communicator() ) ;
   MAC_CHECK_POST( result->global_number() == 0 ) ;
   MAC_CHECK_POST( result->local_number() == 0 ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_DistributedPartition:: MAC_DistributedPartition( MAC_Object* a_owner )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , COMM( MAC_Exec::communicator() )
   , SIZE( MAC_Exec::communicator()->nb_ranks() )
   , RANK( MAC_Exec::communicator()->rank() )
   , FIRST( 0 )
   , LAST( 0 )
   , GLOBAL_NB( 0 )
   , LOCAL_NB( 0 )
   , PARTITION( 0 )
   , START( 0 )
{
   MAC_LABEL( "MAC_DistributedPartition:: MAC_DistributedPartition" ) ;

   PARTITION.re_initialize( SIZE ) ;
   START.re_initialize( SIZE ) ;
}




//----------------------------------------------------------------------
MAC_DistributedPartition:: ~MAC_DistributedPartition( void )
//----------------------------------------------------------------------
{
}




//----------------------------------------------------------------------
void
MAC_DistributedPartition:: set( MAC_DistributedPartition const* other ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DistributedPartition:: set" ) ;
   MAC_CHECK_PRE( other!=0 ) ;

   PARTITION = other->PARTITION ;
   START = other->START ;
   FIRST = other->FIRST ;
   LAST = other->LAST ;
   LOCAL_NB = other->LOCAL_NB ;
   GLOBAL_NB = other->GLOBAL_NB ;

   MAC_CHECK_POST( first_local_index() == other->first_local_index() ) ;
   MAC_CHECK_POST( local_index_limit() == other->local_index_limit() ) ;
   MAC_CHECK_POST( global_number() == other->global_number() ) ;
   MAC_CHECK_POST( local_number() == other->local_number() ) ;
   MAC_CHECK_POST( partitioning() == other->partitioning() ) ;
   MAC_CHECK_POST( start_of_partition() ==  other->start_of_partition() ) ;
}




//----------------------------------------------------------------------
bool
MAC_DistributedPartition:: is_compatible( 
                              MAC_DistributedPartition const* other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DistributedPartition:: is_compatible" ) ;
   MAC_CHECK_PRE( other!=0 ) ;
   
   bool result = ( PARTITION == other->PARTITION ) ;

   MAC_CHECK_POST(
             EQUIVALENT( result, other->partitioning()==partitioning() ) ) ; 
   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC_DistributedPartition:: set_local_number( size_t a_local_number )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DistributedPartition:: set_local_number" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   COMM->all_gather( a_local_number, PARTITION ) ;
   GLOBAL_NB = PARTITION.sum() ;
   START( 0 ) = 0 ;
   for( size_t i=1 ; i<SIZE ; ++i )
   {
      START( i ) = START( i-1 ) + PARTITION( i-1 ) ;
   }   
   FIRST = START( RANK ) ;
   LAST = FIRST + PARTITION( RANK ) ;
   LOCAL_NB = PARTITION( RANK ) ;

   MAC_CHECK_POST( local_number() == a_local_number ) ;
   MAC_CHECK_POST( global_number() ==
                   (  size_t) communicator()->sum( (double) a_local_number ) ) ;
   MAC_CHECK_POST( partitioning()( communicator()->rank() ) == 
   	(int) a_local_number ) ;
   MAC_CHECK_POST( partitioning().sum() == (int) global_number() ) ;
}




//----------------------------------------------------------------------
void
MAC_DistributedPartition:: distribute_global_number( size_t a_global_number )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DistributedPartition:: distribute_global_number" ) ;
   MAC_CHECK_COLLECTIVE( false ) ;

   GLOBAL_NB = a_global_number ;
   size_t const dim = GLOBAL_NB / SIZE ;
   size_t const r   = GLOBAL_NB % SIZE ;
   for( size_t i=0 ; i<SIZE ; ++i )
   {
      PARTITION(i) = ( i<r ? dim+1 : dim ) ;
   }
   START(0) = 0 ;
   for( size_t i=1 ; i<SIZE ; ++i )
   {
      START(i) = START(i-1)+PARTITION(i-1) ;
   }   
   FIRST = START( RANK ) ;
   LAST = FIRST + PARTITION( RANK ) ;
   LOCAL_NB = PARTITION( RANK ) ;

   MAC_CHECK_POST( global_number() == a_global_number ) ;
   MAC_CHECK_POST( partitioning().sum() == (int) global_number() ) ;
   MAC_CHECK_POST(
      local_number() == a_global_number/communicator()->nb_ranks() ||
      local_number() == a_global_number/communicator()->nb_ranks()+1 ) ;
}




//----------------------------------------------------------------------
void
MAC_DistributedPartition:: set_global_number( size_t a_global_number )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DistributedPartition:: set_global_number" ) ;
   MAC_CHECK_COLLECTIVE( false ) ;

   GLOBAL_NB = a_global_number ;
   PARTITION.set( 0 ) ;
   PARTITION( RANK ) = GLOBAL_NB ;
   START(0) = 0 ;
   for( size_t i=1 ; i<SIZE ; ++i )
   {
      START(i) = START(i-1)+PARTITION(i-1) ;
   }   
   FIRST = 0 ;
   LAST =  GLOBAL_NB ;
   LOCAL_NB = GLOBAL_NB ;

   MAC_CHECK_POST( global_number() == a_global_number ) ;
   MAC_CHECK_POST( partitioning().sum() == (int) global_number() ) ;
   MAC_CHECK_POST( partitioning()( communicator()->rank() ) == 
   	a_global_number ) ;
   MAC_CHECK_POST( local_number() == a_global_number ) ;
   MAC_CHECK_POST( first_local_index() == 0 ) ;
   MAC_CHECK_POST( local_index_limit() == a_global_number ) ;
}
