#ifndef MAC_DISTRIBUTED_PARTITION_HH
#define MAC_DISTRIBUTED_PARTITION_HH

#include <MAC_Object.hh>
#include <intVector.hh>

class MAC_Communicator ;

/* 
Distribution of a set of indices on processes.

PUBLISHED
*/

class MAC_DistributedPartition : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_DistributedPartition* create( MAC_Object* a_owner ) ;
      
   //-- Distribution built

      void set( MAC_DistributedPartition const* other ) ;

      void set_local_number( size_t a_local_number ) ;
      
      void distribute_global_number( size_t a_global_number ) ;

      void set_global_number( size_t a_global_number ) ;

   //-- Distribution comparison
     
      // Has `other' the same distributed structure than `self' ?
      bool is_compatible( MAC_DistributedPartition const* other ) const ;
      
   //-- Access

      // communicator
      MAC_Communicator const* communicator( void ) const ;
      
      // global number of items owned by all processes
      size_t global_number( void ) const ;
      
      // local number of items owned by current process  
      size_t local_number( void ) const ;
      
      /* first index of local owned items:
           `::first_local_index()' = `::start_of_partition'( `::communicator'()->rank() ) */
      size_t first_local_index( void ) const ;
      
      /* last index + 1 of local owned items:
           `::local_index_limit()' = `::first_local_index()'+`::local_number'() */
      size_t local_index_limit( void ) const ;

      /* table of partition:
           `::partitioning'( `::communicator'()->rank() ) = `::local_number'()
           `::partitioning'().sum() = `::global_number'() */
      intVector const& partitioning( void ) const ;

      // table of the first index of owned by each process
      intVector const& start_of_partition( void ) const ;

      // rank of the process which owned the index `i'
      size_t rank_of( size_t i ) const ;
      
   protected: //--------------------------------------------------------
     
   private: //----------------------------------------------------------

      MAC_DistributedPartition( MAC_Object* a_owner ) ;
     ~MAC_DistributedPartition( void ) ;
      
      MAC_DistributedPartition( void ) ;
      MAC_DistributedPartition( MAC_DistributedPartition const& other ) ;
      MAC_DistributedPartition& operator=( 
                                MAC_DistributedPartition const& other ) ;

   //-- Attributes
      
      MAC_Communicator const* const COMM ;
      
      size_t const SIZE ;
      size_t const RANK ;      

      size_t FIRST ;
      size_t LAST ;
      size_t GLOBAL_NB ;
      size_t LOCAL_NB ;
      intVector PARTITION ;
      intVector START ;
} ;

#ifndef OUTLINE
   #include <MAC_DistributedPartition.icc>
#endif

#endif



