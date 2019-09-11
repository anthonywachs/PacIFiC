#ifndef MAC_LIST_ITEM_HH
#define MAC_LIST_ITEM_HH

#include <MAC_Object.hh>

/**  Item for MAC_List implementation. */
class MAC_ListItem
{
   public: //---------------------------------------------------------
      MAC_ListItem( MAC_Object* a_val ) : value(a_val), theNext(0) 
      {
      }
            
      ~MAC_ListItem( void ) {}

      void doLink( MAC_ListItem* another ) 
      {
         theNext = another ;
      }
            
      MAC_ListItem * next() 
      {
         return theNext ;
      }

      MAC_Object* val() 
      {
         return value ;
      }
            
      MAC_Object* replaceValue( MAC_Object* newValue ) 
      {
         MAC_Object* old = value ;
         value = newValue ;
         return old ;
      }
            
   protected: //--------------------------------------------------------
   private: //----------------------------------------------------------
      //----------------------------------------------------------------
      // ATTRIBUTES for class MAC_ListItem
      //----------------------------------------------------------------
      MAC_Object* value ;
      MAC_ListItem * theNext ;
} ;


#endif
