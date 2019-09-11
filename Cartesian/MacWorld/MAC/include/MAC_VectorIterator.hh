#ifndef MAC_VECTOR_ITERATOR_HH
#define MAC_VECTOR_ITERATOR_HH

#include <MAC_Iterator.hh>

class MAC_Vector ;

/*
Iterators on items of MAC_Vector collections.

The sequence of access is the order of increasing indices of items.
*/

class MAC_VectorIterator : public MAC_Iterator
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an iterator over non NULL items of `vector'.
      static MAC_VectorIterator* create( MAC_Object* a_owner,
                                         MAC_Vector const* vector ) ;

   //-- Status report

      virtual bool is_valid( void ) const ;

   //-- Access

      size_t index_of_item( void ) const ;

      virtual MAC_Object* item( void ) const ;

   //-- Cursor movement

      void go_i_th( size_t i ) ;

      virtual void start( void ) ;

      virtual void go_next( void ) ;
             
   protected: //------------------------------------------------------------

      virtual ~MAC_VectorIterator( void ) ;

      MAC_VectorIterator( MAC_Object* a_owner,
                          MAC_Vector const* vector ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   private: //--------------------------------------------------------------

      MAC_VectorIterator( void ) ;
      MAC_VectorIterator( MAC_VectorIterator const& other ) ;
      MAC_VectorIterator const& operator=( MAC_VectorIterator const& other ) ;

   //-- Attributes

      MAC_Vector const* vec ;
      bool valid ;
      size_t current_i ;
      MAC_Object* current_item ;
} ;

#endif
