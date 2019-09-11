#ifndef LA_DISTRIBUTED_SCATTER_HH
#define LA_DISTRIBUTED_SCATTER_HH

#include <LA_Scatter.hh>

#include <intVector.hh>
#include <size_t_vector.hh>

class MAC_Communicator ;
class MAC_DistributedPartition ;

class LA_DistMatrix ;
class LA_DistVector ;
class LA_Matrix ;
class LA_MatrixIterator ;
class LA_SeqMatrix ;
class LA_SeqVector ;

/*
  `LA_Scatter::' server for `LA_DistImplementation::' objects.

  PUBLISHED
*/

class LA_DistScatter : public LA_Scatter
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static LA_DistScatter* create(
                        MAC_Object* a_owner,
                        MAC_DistributedPartition const* partition ) ;

      // Reinitialize the internal state, as if `::create' was just completed.
      void re_initialize( MAC_DistributedPartition const* partition ) ;

   //-- Characteristics

      virtual LA_Implementation const* implementation( void ) const ;

      virtual size_t size( void ) const ;

      virtual size_t_vector const& repatriated_items( void ) const ;

      virtual size_t_vector const& local_indices( void ) const ;

      virtual MAC_DistributedPartition const* distribution( void ) const ;

   //-- Repatriate items

      // Set the items which have to be repatriated
      // (`a_repatriated_items_table' is supposed to be sorted).
      void set_sorted( size_t_vector const& a_repatriated_items_table,
                       size_t_vector const& a_local_indices_table ) ;

      // Set the items which have to be repatriated.
      void set_unsorted( size_t_vector const& a_repatriated_items_table,
                         size_t_vector const& a_local_indices_table ) ;

      virtual void get( LA_Vector const* source,
                        LA_SeqVector* dest ) const ;

      // Operator `::get'( source, dest ) is split in
      //      - `::get_begin'( source )
      //      - `::get_end'( dest )
      // hence, operations which not need repatriated items could
      // be done as items are repatriated for time optimisation.
      void get_begin( LA_DistVector const* source ) const ;
      void get_end( LA_SeqVector* dest ) const ;
      bool is_getting( void ) const ;

      virtual void set( LA_SeqVector const* source,
                        LA_Vector* dest ) const ;

      // Repatriate distributed items in distributed `source' to local `dest'.
      void get_rows( LA_DistMatrix const* source,
                     LA_SeqMatrix* dest ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      LA_DistScatter( void ) ;
     ~LA_DistScatter( void ) ;
      LA_DistScatter( LA_DistScatter const& other ) ;
      LA_DistScatter& operator=( LA_DistScatter const& other ) ;

      LA_DistScatter( MAC_Object* a_owner,
                      MAC_DistributedPartition const* partition ) ;

   //-- Distributed processing

      void send_index_to( size_t start, size_t rank ) const ;

      void receive_index_from( size_t start, size_t rank ) ;

      void send_values_to( LA_DistVector const* source,
                           size_t start,
                           size_t rank ) const ;

      void* send_values_NB_to( LA_DistVector const* source,
                               size_t start,
                               size_t rank,
                               double * ptr ) const ;

      void receive_values_from( LA_SeqVector* dest,
                                size_t start,
                                size_t rank ) const ;

      void send_rows_to( LA_MatrixIterator* it,
                         LA_DistMatrix const* source,
                         size_t& start,
                         size_t rank ) const ;

      void receive_rows_from(
                         LA_SeqMatrix* dest,
                         size_t& start,
                         size_t rank ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool implementation_POST(
                            LA_Implementation const* result ) const ;

   //-- Attributes

      MAC_DistributedPartition* DIST ;
      MAC_Communicator const* COMM ;

      intVector RECEIVED_NB ;  // Nb elements / process to receive
      intVector SENT_NB ;      // Nb elements / process to send

      intVector RECEIVED_IDX ; // Indices of received elements
      intVector SENT_IDX ;     // Indices of sent elements

      size_t_vector ITEMS ;
      size_t_vector LOCAL_IDS ;
      mutable void ** REQUEST ;
} ;

#endif
