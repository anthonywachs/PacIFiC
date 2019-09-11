#ifndef MAC_LIST_IDENTITY_HH
#define MAC_LIST_IDENTITY_HH

#include <MAC_List.hh>

//--------------------------------------------------------------------------
//  Lists whose items are matching provided that they have the same adress
//  (i.e. the same identity).
//--------------------------------------------------------------------------


class MAC_ListIdentity : public MAC_List
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance.
      static MAC_ListIdentity* create( MAC_Object* a_owner ) ;

      virtual MAC_ListIdentity* create_clone( MAC_Object* a_owner ) const ;

   //-- Status

      // IMPLEMENTATION : true if `object1' and 'object2' 
      //                  have the same adress
      virtual bool matching_items( MAC_Object const* object1,
                                   MAC_Object const* object2 ) const ;

   protected: //------------------------------------------------------------

      virtual ~MAC_ListIdentity( void ) ;

      MAC_ListIdentity( MAC_Object* a_owner ) ;
            
   private: //--------------------------------------------------------------

      MAC_ListIdentity( void ) ;
      MAC_ListIdentity( MAC_ListIdentity const& other ) ;
      MAC_ListIdentity const& operator=( MAC_ListIdentity const& other ) ;
      
} ;

#endif
