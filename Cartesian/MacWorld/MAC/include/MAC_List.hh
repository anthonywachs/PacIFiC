#ifndef MAC_LIST_HH
#define MAC_LIST_HH

#include <MAC_Sequence.hh>

#include <MAC_ListIterator.hh>

class MAC_ListItem ;

/*
Sequential, one-way linked lists
*/

class MAC_List : public MAC_Sequence
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance.
      static MAC_List* create( MAC_Object* a_owner ) ;

      virtual MAC_List* create_clone( MAC_Object* a_owner ) const ;

      // Reinitialize by copying all the items of `other'.
      void copy( MAC_List const* other ) ;

  //-- Measurement

      virtual size_t index_limit( void ) const ;

      virtual size_t count( void ) const ;

   //-- Access

      virtual MAC_Object* item( MAC_Object const* object ) const ;

      virtual MAC_Object* at( size_t i ) const ;
         
      virtual size_t index_of( MAC_Object const* object ) const ;

      virtual MAC_ListIterator* create_iterator( MAC_Object* a_owner ) const ;

   //-- Element change

      virtual void append( MAC_Object* object ) ;
      
      virtual void prepend( MAC_Object* object ) ;

      virtual void set_at( size_t i, MAC_Object* object ) ;

      virtual void insert_at( size_t i, MAC_Object* object ) ;

   //-- Removal

      virtual void remove( MAC_Object const* object ) ;
         
      virtual void remove_at( size_t i ) ;
      
      virtual void remove_section( size_t iFirst, size_t length ) ;
      
      virtual void destroy_items_and_clear( void ) ;
      
      virtual void destroy_items_and_remove_section( size_t iFirst, 
                                                 size_t length ) ;

      virtual void clear( void ) ;

   protected: //------------------------------------------------------------

      virtual ~MAC_List( void ) ;

      MAC_List( MAC_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

      virtual bool set_at_POST( size_t old_index_limit,
                                size_t old_count,
                                size_t i,
                                MAC_Object const* object,
                                size_t old_state_id ) const ;
      virtual bool at_POST( MAC_Object const* result, size_t i ) const ;

      virtual bool remove_at_POST( size_t old_index_limit,
                                   size_t old_count,
                                   size_t old_state_id ) const ;

      virtual bool remove_section_POST( size_t old_index_limit,
                                        size_t old_count,
                                        size_t length,
                                        size_t old_state_id ) const ;
      
   private: //--------------------------------------------------------------

      friend class MAC_ListIterator ;

      MAC_List( void ) ;
      MAC_List( MAC_List const& other ) ;
      MAC_List const& operator=( MAC_List const& other ) ;

      MAC_ListItem* theItem( size_t n ) const ;
      
      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      MAC_ListItem* theList ;
      MAC_ListItem* theLast ;
      size_t nbElem ;
      
} ;


#endif
