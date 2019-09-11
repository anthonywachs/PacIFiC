#ifndef MAC_HASH_TABLE_SET_ITERATOR_HH
#define MAC_HASH_TABLE_SET_ITERATOR_HH

#include <MAC_Iterator.hh>

class MAC_HashTableSet ;

//---------------------------------------------------------------------------
// Iterators on items of MAC_HashTableSet collections
//---------------------------------------------------------------------------
// The sequence of access is the internal ordering of items within the
// table : the first bucket is traversed, then the second and so on.
//---------------------------------------------------------------------------

class MAC_HashTableSetIterator : public MAC_Iterator
{
   public: //---------------------------------------------------------------

   //-- Initialization

      // Create and return an iterator over items of `a_table'.
      static MAC_HashTableSetIterator* create( MAC_Object* a_owner,
                                              MAC_HashTableSet const* a_table ) ;

   //-- Status report

      virtual bool is_valid( void ) const ;

   //-- Access

      virtual MAC_Object* item( void ) const ;

   //-- Cursor movement

      virtual void start( void ) ;

      virtual void go_next( void ) ;
       
   protected: //------------------------------------------------------------

      virtual ~MAC_HashTableSetIterator( void ) ;

      MAC_HashTableSetIterator( MAC_Object* a_owner, 
                             MAC_HashTableSet const* hTable ) ;
      
   private: //-------------------------------------------------------------- 

      MAC_HashTableSetIterator( void ) ;
      MAC_HashTableSetIterator( MAC_HashTableSetIterator const& other ) ;
      MAC_HashTableSetIterator const& operator=( 
                                MAC_HashTableSetIterator const& other ) ;
      
      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      MAC_HashTableSet const* table ;
      MAC_Iterator* where ;
      bool init ;
      size_t bucket ;
} ;

#endif
