#ifndef MAC_LEXICAL_H
#define MAC_LEXICAL_H

#include <MAC_Object.hh>

#include <string>

/* Provides a common interface to Pelicans data structures for reading. */

class MAC_KeywordDataPair ;
class MAC_Module ;
class MAC_List ;
class MAC_Data ;
class MAC_String ;

class MAC_Lexical : public MAC_Object
{
   public: //-------------------------------------------------------
      
      // Initialization
      static MAC_Lexical* create( MAC_Module* aModule ) ;
      static MAC_Lexical* create( MAC_Data * aSimpleValue ) ;
      static MAC_Lexical* create( MAC_List * aList ) ;
      
   //-- Access

      // Is self a module ?
      bool is_module( void ) const;
      
      // Is self a simple value ?
      bool is_data( void ) const ;
      
      // Is self a list ?
      bool is_list( void ) const ;
      
   //-- Retrieving
      
      // Retrieves self as a module.
      MAC_Module * to_module( void ) ;
      
       // Retrieves self as a simple value.
      MAC_Data* to_data( void ) ;
      
       // Retrieves self as a list.
      MAC_List * to_list( void ) ;
      
  //-- Removing   
      static void remove_all_lexical( void ) ;
      
   protected: //-------------------------------------------------------
      
      // Destructor
      virtual ~MAC_Lexical( void ) ;
      
   private: //-------------------------------------------------------
      
      MAC_Lexical( MAC_Lexical const& other  ) ;
      MAC_Lexical const& operator=( MAC_Lexical const& other  ) ;     

      MAC_Lexical( void ) ;
      MAC_Lexical( MAC_Object* a_owner, MAC_Module* aModule ) ;
      MAC_Lexical( MAC_Object* a_owner, MAC_KeywordDataPair * anAssignment ) ;
      MAC_Lexical( MAC_Object* a_owner, MAC_Data * aSimpleValue ) ;
      MAC_Lexical( MAC_Object* a_owner, MAC_List * aSimpleValue ) ;
      
      bool invariant( void ) const ;
      static MAC_Lexical* common_owner( void ) ;
      static MAC_Lexical* the_common_owner ;
      
      //-------------------------------------------------------------
      //   ATTRIBUTES
      //-------------------------------------------------------------
      MAC_Module* myModule ;
      MAC_Data* myData ;
      MAC_List* myList ;
};

typedef MAC_Lexical* MAC_LexicalPtr ;

#endif 
