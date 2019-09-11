#ifndef MAC_HASH_TABLE_SET_HH
#define MAC_HASH_TABLE_SET_HH

#include <MAC_Set.hh>

#include <MAC_HashTableSetIterator.hh>
#include <MAC_List.hh>

//---------------------------------------------------------------------------
//  Sets implemented as hash tables
//---------------------------------------------------------------------------

class MAC_Vector ;

class MAC_HashTableSet : public MAC_Set
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return a new instance.
      static MAC_HashTableSet* create( MAC_Object* a_owner,
                                       size_t size = 20 ) ; 

      virtual MAC_HashTableSet* create_clone( MAC_Object* a_owner ) const ;

   //-- Measurement   

      // number of buckets
      size_t nb_buckets( void ) const ;

      virtual size_t count( void ) const ;

   //-- Access

      virtual MAC_Object* item( MAC_Object const* object ) const ;

      virtual MAC_HashTableSetIterator* create_iterator( 
                                             MAC_Object* a_owner ) const ;

   //-- Element change

      virtual void extend( MAC_Object* object ) ;

   //-- Removal
      
      virtual void remove( MAC_Object const* object ) ;

      virtual void destroy_item_and_remove( MAC_Object const* object ) ;

      virtual void clear( void ) ;

      virtual void destroy_items_and_clear( void ) ;

   protected: //------------------------------------------------------------

      virtual ~MAC_HashTableSet( void ) ;

      MAC_HashTableSet( MAC_Object* a_owner, size_t size ) ;

   //-- Preconditions, Postconditions, Invariant     
 
      virtual bool invariant( void ) const ;

   private: //--------------------------------------------------------------

      friend class MAC_HashTableSetIterator ;

      MAC_HashTableSet( void ) ;
      MAC_HashTableSet( MAC_HashTableSet const& other ) ;
      MAC_HashTableSet const& operator=( MAC_HashTableSet const& other ) ;
      
      MAC_List * getList( size_t entry ) ;
      MAC_List const* getList( size_t entry ) const ;
      
      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      MAC_Vector* hTable ;
      size_t the_size ;
      
} ;


#endif
