#ifndef MAC_COLLECTION_HH
#define MAC_COLLECTION_HH

#include <MAC_Container.hh>

/*
Containers for which basic insertion and removal of items is  meaningfull
*/

class MAC_Collection : public MAC_Container
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      virtual MAC_Collection* create_clone( MAC_Object* a_owner ) const = 0 ;
           
   //-- Element change 

      // Ensure that self includes an item matching `object'.
      virtual void extend( MAC_Object* object ) = 0 ;

   //-- Removal

      // Remove the first item matching `object'.
      virtual void remove( MAC_Object const* object ) = 0 ;

      // Remove and terminate the first item matching `object'.  */
      virtual void destroy_item_and_remove( MAC_Object const* object ) = 0 ;
      
      // Remove all items.
      virtual void clear( void ) = 0 ;

      // Remove and terminate all items.
      virtual void destroy_items_and_clear( void ) = 0 ;

   protected: //--------------------------------------------------------

      virtual ~MAC_Collection( void ) ;

      MAC_Collection( MAC_Object* a_owner ) ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool extend_PRE( MAC_Object const* object ) const ;
      virtual bool extend_POST( bool old_has,
                                size_t old_count,
                                MAC_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool remove_PRE( MAC_Object const* object ) const ;
      virtual bool remove_POST( size_t old_count,
                                MAC_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool clear_POST( size_t old_state_id ) const ;
      
   private: //----------------------------------------------------------

      MAC_Collection( void ) ;
      MAC_Collection( MAC_Collection const& other ) ;
      MAC_Collection const& operator=( MAC_Collection const& other ) ;
      
} ;

#endif
