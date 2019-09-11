#ifndef MAC_CONTAINER_HH
#define MAC_CONTAINER_HH

#include <MAC_Object.hh>

class MAC_Iterator ;

/*
Data structures
  - used to hold zero or more items, where item denotes a possibly NULL
    object ;
  - with a finite item count ;
  - for which there exist a traversal policy that will visit every item
    exactly once.
The items are referred to, accessed and manipulated exclusively through
pointers. The NULL object is thus represented by the 0 pointer.
*/

class MAC_Container : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      virtual MAC_Container* create_clone( MAC_Object* a_owner ) const = 0 ;

   //-- Status

      // Do `object1' and `object2' compare equal for search operations ?
      // IMPLEMENTATION : use of is_equal()
      virtual bool matching_items( MAC_Object const* object1,
                                   MAC_Object const* object2 ) const ;

      // identification of the object state (two different states have
      // different identifiers)
      size_t state_id( void ) const ;

   //-- Measurement

      // number of non NULL items
      virtual size_t count( void ) const = 0 ;

   //-- Access

      // Does self include an item matching `object' ?
      bool has( MAC_Object const* object ) const ;

      // first item matching object if any ; O otherwise
      virtual MAC_Object* item( MAC_Object const* object ) const = 0 ;

      // Create and return an iterator on non NULL items.
      virtual MAC_Iterator* create_iterator( MAC_Object* a_owner ) const = 0 ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;      
      
   protected: //--------------------------------------------------------

      virtual ~MAC_Container( void ) ;

      MAC_Container( MAC_Object* a_owner ) ;
      
      bool create_clone_POST( MAC_Container* result,
                              MAC_Object* a_owner ) const ;

      virtual bool count_POST( size_t result ) const ;

      virtual bool matching_items_PRE( MAC_Object const* object1,
                                       MAC_Object const* object2 ) const ;

      virtual bool item_PRE( MAC_Object const* object ) const ;
      virtual bool item_POST( MAC_Object const* result, 
                              MAC_Object const* object ) const ;

      virtual bool create_iterator_POST( MAC_Iterator* result,
                                         MAC_Object* a_owner ) const ;

      void report_state_change( void ) ;
      
   private: //----------------------------------------------------------

      MAC_Container( void ) ;
      MAC_Container( MAC_Container const& other ) ;
      MAC_Container const& operator=( MAC_Container const& other ) ;

   //-- Attributes
      
      size_t MY_STATE ;
} ;

#endif
