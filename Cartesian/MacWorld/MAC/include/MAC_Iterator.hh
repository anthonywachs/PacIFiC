#ifndef MAC_ITERATOR_HH
#define MAC_ITERATOR_HH

#include <MAC_Object.hh>

class MAC_Container ;

/*
Iterators to traverse data structures and visit every non NULL items 
exactly once

Items are accessed sequentially, one-way.
The sequence of access depends on the underlying data structure.
*/

class MAC_Iterator : public MAC_Object
{
   public: //---------------------------------------------------------------

   //-- Status report

      // Has the traversed container been modified since last call 
      // to `::start' ?
      virtual bool container_has_been_modified( void ) const ;

      // Is current position valid ?
      virtual bool is_valid( void ) const = 0 ;

   //-- Access

      // item at the current position
      virtual MAC_Object* item( void ) const = 0 ;

   //-- Cursor movement

      // Move to first position.
      virtual void start( void ) = 0 ;

      // Move to next position.
      virtual void go_next( void ) = 0 ;
      
   protected: //------------------------------------------------------------
      
      virtual ~MAC_Iterator( void ) ;

      MAC_Iterator( MAC_Object* a_owner, MAC_Container const* list ) ;
      
      MAC_Iterator( MAC_Object* a_owner ) ;
      
      void base_initialize( MAC_Container const* list ) ;

      void record_container_state_id( void ) ;
      
   //-- Preconditions, Postconditions, Invariant   
   
      virtual bool invariant( void ) const ;

      virtual bool start_POST( void ) const ;
      
      virtual bool item_PRE( void ) const ;
      virtual bool is_valid_PRE( void ) const ;
      virtual bool go_next_PRE( void ) const ;
      virtual bool item_POST( MAC_Object const* result ) const ;
      
   private: //--------------------------------------------------------------

      MAC_Iterator( void ) ;
      MAC_Iterator( MAC_Iterator const& other ) ;
      MAC_Iterator const& operator=( MAC_Iterator const& other ) ;

   //-- Attributes

      MAC_Container const* my_list ;
      size_t list_state ;

} ;

#endif
