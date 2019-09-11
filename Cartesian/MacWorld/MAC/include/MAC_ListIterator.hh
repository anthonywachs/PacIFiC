#ifndef MAC_LIST_ITERATOR_HH
#define MAC_LIST_ITERATOR_HH

#include <MAC_Iterator.hh>

class MAC_List ;
class MAC_ListItem ;

/*
Iterators on items of MAC_List collections

The sequence of access is the order in which items have been inserted
in the list
*/

class MAC_ListIterator : public MAC_Iterator
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an iterator over items of `a_list'.
      static MAC_ListIterator* create( MAC_Object* a_owner,
                                       MAC_List const* a_list ) ;
      
      static MAC_ListIterator* create( MAC_Object* a_owner ) ;
      
      void re_initialize( MAC_List const* a_list ) ;

   //-- Status report

      virtual bool is_valid( void ) const ;

   //-- Access

      virtual MAC_Object* item( void ) const ;

   //-- Cursor movement

      virtual void start( void ) ;
  
      virtual void go_next( void ) ;       
      
   protected: //------------------------------------------------------------

      virtual ~MAC_ListIterator( void ) ;

      MAC_ListIterator( MAC_Object* a_owner, MAC_List const* aList ) ;
      
      MAC_ListIterator( MAC_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;
      
   private: //--------------------------------------------------------------

      MAC_ListIterator( void ) ;
      MAC_ListIterator( MAC_ListIterator const& other ) ;
      MAC_ListIterator& operator=( MAC_ListIterator const& other ) ;

   //-- Attributes
      
      MAC_List const* list ;
      MAC_ListItem* where ;
      bool init ;
} ;

#endif
