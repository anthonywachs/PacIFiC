#ifndef MAC_SET_HH
#define MAC_SET_HH

#include <MAC_Collection.hh>

//---------------------------------------------------------------------------
//   Collections with no NULL items and no matching different items.
//---------------------------------------------------------------------------

class MAC_Set : public MAC_Collection
{
   public : //----------------------------------------------------------

   //-- Instance delivery and initialization

      virtual MAC_Set* create_clone( MAC_Object* a_owner ) const = 0 ;

   //-- Element change 

      virtual void extend( MAC_Object* object ) = 0 ;
      
   //-- Removal

      virtual void remove( MAC_Object const* object ) = 0 ;

   protected : //-------------------------------------------------------

      virtual ~MAC_Set( void ) ;

      MAC_Set( MAC_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool extend_POST( bool old_has, 
                                size_t old_count, 
                                MAC_Object const* object,
                                size_t old_state_id ) const ;

      virtual bool remove_POST( size_t old_count,
                                MAC_Object const* object,
                                size_t old_state_id ) const ;

   private: //----------------------------------------------------------

      MAC_Set( void ) ;
      MAC_Set( MAC_Set const& other ) ;
      MAC_Set const& operator=( MAC_Set const& other ) ;
} ;

#endif
