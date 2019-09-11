#ifndef MAC_SEQUENCE_HH
#define MAC_SEQUENCE_HH

#include <MAC_Collection.hh>

/*
Collections whose items are in one-to-one correspondance with
integers of a contiguous interval with zero as lower bound.
items can be the NULL object.
*/

class MAC_Sequence : public MAC_Collection
{

   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      virtual MAC_Sequence* create_clone( MAC_Object* a_owner ) const = 0 ;

   //-- Comparison

      virtual bool comparable( MAC_Object const* other ) const ;

      virtual size_t hash_code( void ) const ;

      // IMPLEMENTATION : true if three_way_comparison() returns 0
      virtual bool is_equal( MAC_Object const* other ) const ;

      // IMPLEMENTATION : lexicographic order
      virtual int three_way_comparison( MAC_Object const* other ) const ;

   //-- Measurement

      // exclusive upper limit for indices
      virtual size_t index_limit( void ) const = 0 ;

   //-- Access
      
      /** Class constant returned by index when no element of self match
          is_equal relation. */
      static size_t const badIndex ;

      // item at index `i'
      virtual MAC_Object* at( size_t i ) const = 0 ;

      // index of the first item matching `object' if any ; a number out
      // of the valid index interval if none
      virtual size_t index_of( MAC_Object const* object ) const = 0 ;

   //-- Sorting

      virtual void sort( void ) ;

   //-- Element change

      virtual void extend( MAC_Object* object ) ;

      // Add `object' to end.
      virtual void append( MAC_Object* object ) = 0 ;

      // Add `object' to beginning.
      virtual void prepend( MAC_Object* object ) = 0  ;

      // Assign `object' to the `i'-th entry.
      virtual void set_at( size_t i, MAC_Object* object ) = 0  ;

      // Increase by one the index of all items starting from the `i'-th 
      // position and add `object' to the `i'-th position.
      virtual void insert_at( size_t i, MAC_Object* object ) = 0  ;

   //-- Removal

      virtual void remove( MAC_Object const* object ) ;

      virtual void destroy_item_and_remove( MAC_Object const* object ) ;
      
      // Remove item at index `i'.
      virtual void remove_at( size_t i ) = 0 ;

      // Terminate and remove item at index `i'.
      void destroy_item_and_remove_at( size_t i ) ;

      // Remove `length' contiguous items, starting from the one of index
      // `iFirst'.
      virtual void remove_section( size_t iFirst, size_t length ) = 0 ;

      // Terminate and remove `length' contiguous items, starting from the 
      // one of index `iFirst'.
      virtual void destroy_items_and_remove_section( size_t iFirst,
                                                     size_t length ) = 0 ;

      virtual void destroy_items_and_clear( void ) ;
      
      virtual void clear( void ) = 0 ;
      
   //-- Statics

   protected: //------------------------------------------------------------

      virtual ~MAC_Sequence( void ) ;

      MAC_Sequence( MAC_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

      bool create_clone_POST( MAC_Sequence* result,
                              MAC_Object* a_owner ) const ;

      virtual bool append_PRE( MAC_Object const* object ) const ;
      virtual bool append_POST( size_t old_index_limit,
                                size_t old_count,
                                MAC_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool prepend_PRE( MAC_Object const* object ) const ;
      virtual bool prepend_POST( size_t old_index_limit,
                                 size_t old_count,
                                 MAC_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool set_at_PRE( size_t i, MAC_Object const* object ) const ;
      virtual bool set_at_POST( size_t old_index_limit,
                                size_t old_count,
                                size_t i,
                                MAC_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool insert_at_PRE( size_t i, MAC_Object const* object ) const ;
      virtual bool insert_at_POST( size_t old_index_limit,
                                   size_t old_count,
                                   size_t i,
                                   MAC_Object const* object,
                                   size_t old_state_id ) const ;

      virtual bool remove_at_PRE( size_t i ) const ;
      virtual bool remove_at_POST( size_t old_index_limit,
                                   size_t old_count,
                                   size_t old_state_id ) const ;

      virtual bool remove_section_PRE( size_t iFirst, size_t length ) const ;
      virtual bool remove_section_POST( size_t old_index_limit,
                                        size_t old_count,
                                        size_t length,
                                        size_t old_state_id ) const ;

      virtual bool at_PRE( size_t i ) const ;
      virtual bool at_POST( MAC_Object const* resu, size_t i ) const ;

      virtual bool index_of_PRE( MAC_Object const* object ) const ;
      virtual bool index_of_POST( size_t result, MAC_Object const* object ) const ;
      
      virtual bool destroy_items_and_remove_section_PRE( size_t iFirst,
                                                         size_t length ) const ;

      virtual bool clear_POST( size_t old_state_id ) const ;

   private: //--------------------------------------------------------------

      MAC_Sequence( void ) ;
      MAC_Sequence( MAC_Sequence const& other ) ;
      MAC_Sequence const& operator=( MAC_Sequence const& other ) ;

} ;

#endif
