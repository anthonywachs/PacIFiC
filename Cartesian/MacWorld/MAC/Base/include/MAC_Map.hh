#ifndef MAC_MAP_HH
#define MAC_MAP_HH

#include <MAC_Collection.hh>

#include <MAC_MapIterator.hh>

//----------------------------------------------------------------------------
// Data structures, used to store non NULL items identified by keys that may
// be hashed into an integer index
//----------------------------------------------------------------------------
// Implemented using a hash table of key/item pairs
//----------------------------------------------------------------------------
// To find an item associated to a given key, that key is first hashed into
// an integer index to identify the bucket in which the key/item is stored.
// Then, the key is searched within that bucket with a comparing criteria 
// based on the is_equal method.
//
// The number of buckets is set at creation. A number of (key/value) 
// pairs greatly exceeding the number of buckets will lead to a
// significant loss of efficiency when retrieving an item associated to a
// given key since the searching algorithm within a bucket is linear.
//----------------------------------------------------------------------------


class MAC_HashTableSet ;

class MAC_Map : public MAC_Container
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return a new instance.
      static MAC_Map* create( MAC_Object* a_owner, size_t size = 20 ) ; 

      virtual MAC_Map* create_clone( MAC_Object* a_owner ) const ;

   //-- Measurement

      // number of buckets of the underlying hash table
      size_t nb_buckets( void ) const ;

      virtual size_t count( void ) const ;

   //-- Access

      // Is there an item with a key comparing equal to `key' ?
      bool has_key( MAC_Object const* key ) const ;

      virtual MAC_Object* item( MAC_Object const* object ) const ;

      // item of key comparing equal to `key' if any, 0 otherwise 
      MAC_Object* item_at( MAC_Object const* key ) const ;

      virtual MAC_MapIterator* create_iterator( MAC_Object* a_owner ) const ;

   //-- Element change

      // Ensure that the item of key comparing equal to `key' is `a_item'.
      void set_item_at( MAC_Object* key, MAC_Object* a_item ) ;
   
   //-- Removal         

      // Remove the pair key/item of key comparing equal to `key'.
      void remove_at( MAC_Object const* key ) ;

      // Remove all pairs key/item.
      void clear( void ) ;

      
   protected: //-----------------------------------------------------------

   private: //-------------------------------------------------------------

      friend class MAC_MapIterator ;

      MAC_Map( void ) ;
     ~MAC_Map( void ) ;
      MAC_Map( MAC_Map const& other ) ;
      MAC_Map const& operator=( MAC_Map const& other ) ;

      MAC_Map( MAC_Object* a_owner, size_t aSize ) ;
      MAC_Map( MAC_Object* a_owner, MAC_Map const* other ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

   //-- Attributes

      MAC_HashTableSet* hList ;
} ;


#endif
