#ifndef MAC_OBJECT_REGISTER_HH
#define MAC_OBJECT_REGISTER_HH

#include <MAC_Object.hh>

#include <string>

/*
Registers of `MAC_Object::' instances, eg usable in pluggable factories.

PUBLISHED   
*/

class MAC_Iterator ;
class MAC_Map ;

class MAC_ObjectRegister : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization
      
      static MAC_ObjectRegister* create( MAC_Object* a_owner,
                                         std::string const& a_register_name ) ;
      
   //-- Characteristics

      // name of `self'
      std::string const& register_name( void ) const ;
      
   //-- Access

      // Is there an object registered under the name `a_name' ?
      bool has( std::string const& a_name ) const ;

      // object registered under the name `a_name'
      // (fatal error raised if none)
      MAC_Object* item( std::string const& a_name ) const ;

      // Create and return an iterator on registered objects.
      MAC_Iterator* create_iterator( MAC_Object* a_owner ) const ;
      
   //-- Registration
      
      // Register `an_item' under the name `a_name'.
      // (fatal error raised if such a registration name already exists).
      void register_item( std::string const& a_name, MAC_Object* an_item ) ;

      // Suppress the object that was registered under the name `a_name'
      // (fatal error raised if none).
      void unregister_item( std::string const& a_name ) ;

      // Suppress the registered object `an_item'
      // (fatal error raised if none).
      void unregister_item( MAC_Object* an_item ) ;

   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      MAC_ObjectRegister( void ) ;
     ~MAC_ObjectRegister( void ) ;
      MAC_ObjectRegister( MAC_ObjectRegister const& other ) ;
      MAC_ObjectRegister& operator=( MAC_ObjectRegister const& other ) ;

      MAC_ObjectRegister( MAC_Object* a_owner,
                          std::string const& a_register_name ) ;
      
   //-- Attributes
      
      MAC_Map* const REGISTER ;
      std::string const NAME ;
} ;


#endif

