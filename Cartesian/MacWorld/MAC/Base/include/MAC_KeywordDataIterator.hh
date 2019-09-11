#ifndef MAC_KEYWORD_DATA_ITERATOR_HH
#define MAC_KEYWORD_DATA_ITERATOR_HH

#include <MAC_ListIterator.hh>
#include <MAC_KeywordDataPair.hh>

class MAC_List ;

/*
   Iterators on MAC_KeywordDataPair objects.

   Return the `MAC_KeywordDataPair::' objects stored in a `MAC_List::'.
*/

class MAC_KeywordDataIterator : public MAC_ListIterator
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance for traversing `a_list' whose
      // items are `MAC_KeywordDataPair::' objects.
      static MAC_KeywordDataIterator* create( MAC_Object* a_owner,
                                              MAC_List const* a_list ) ;
      
   //-- Access

      virtual MAC_KeywordDataPair* item( void ) const ;
      
   protected: //------------------------------------------------------------
      
   private: //--------------------------------------------------------------

      MAC_KeywordDataIterator( void ) ;
     ~MAC_KeywordDataIterator( void ) ;
      MAC_KeywordDataIterator( MAC_KeywordDataIterator const& other ) ;
      MAC_KeywordDataIterator const& operator=(
                               MAC_KeywordDataIterator const& other ) ;
      MAC_KeywordDataIterator( MAC_Object* a_owner,
                               MAC_List const* a_list ) ;
      
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;      
} ;

#endif
