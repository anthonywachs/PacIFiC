#ifndef MAC_KEY_ITEM_PAIR_HH
#define MAC_KEY_ITEM_PAIR_HH

#include <MAC_Object.hh>

class MAC_KeyItemPair : public MAC_Object
{
      
   public: //---------------------------------------------------------------


   //-- Initialization

      static MAC_KeyItemPair* create( MAC_Object* a_owner,
                                      MAC_Object* a_key,
                                      MAC_Object* a_item ) ; 

   //-- Access

      virtual size_t hash_code( void ) const ;

      MAC_Object* key( void ) const ;

      MAC_Object* item( void ) const ;

   //-- Element change

      void set_item( MAC_Object* object ) ;

   //-- Comparison

      virtual bool comparable( MAC_Object const* other ) const ;

      virtual bool is_equal( MAC_Object const* other ) const ;

      virtual int three_way_comparison( MAC_Object const* other ) const ;

      
   protected: //----------------------------------------------------

      virtual bool invariant() const ;

   private: //------------------------------------------------------

      MAC_KeyItemPair( void ) ;
     ~MAC_KeyItemPair( void ) ;
      MAC_KeyItemPair( MAC_KeyItemPair const& other ) ;
      MAC_KeyItemPair const& operator=( MAC_KeyItemPair const& other ) ;
                      
      MAC_KeyItemPair( MAC_Object* a_owner,
                       MAC_Object* a_key,
                       MAC_Object* a_item ) ;
      
      //------------------------------------------------------------
      //   ATTRIBUTES
      //------------------------------------------------------------
      MAC_Object* the_key ;
      MAC_Object* the_item ;

} ;


#endif
