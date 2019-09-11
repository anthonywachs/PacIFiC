#ifndef MAC_MAP_ITERATOR_HH
#define MAC_MAP_ITERATOR_HH

#include <MAC_Iterator.hh>

class MAC_Map ;
class MAC_Iterator ;

//--------------------------------------------------------------------------
// Iterators on items of MAC_Map containers
//--------------------------------------------------------------------------
// The sequence of access is the internal ordering of items within the map :
// the first bucket is traversed first, then the second and so on.
//--------------------------------------------------------------------------

class MAC_MapIterator : public MAC_Iterator
{
   public: //---------------------------------------------------------------

   //-- Initialization

      // Create an return an iterator over items of `map'
      static MAC_MapIterator* create( MAC_Object* a_owner,
                                      MAC_Map const* map ) ;

   //-- Status report

      virtual bool is_valid( void ) const ;

   //-- Access

      virtual MAC_Object* item( void ) const ;

      // key of the item at current position
      MAC_Object* key( void ) const ;

   //-- Cursor movement

      virtual void start( void ) ;

      virtual void go_next( void ) ;
       
   protected: //------------------------------------------------------------

      virtual ~MAC_MapIterator( void ) ;

      MAC_MapIterator( MAC_Object* a_owner, 
                       MAC_Map const* hTable ) ;
      
   private: //-------------------------------------------------------------- 

      MAC_MapIterator( void  ) ;
      MAC_MapIterator( MAC_MapIterator const& other ) ;
      MAC_MapIterator const& operator=( MAC_MapIterator const& other ) ;

      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      MAC_Map const* table ;
      MAC_Iterator* where ;
} ;

#endif
