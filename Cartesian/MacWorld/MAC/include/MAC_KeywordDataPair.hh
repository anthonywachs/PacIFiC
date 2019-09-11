#ifndef MAC_KEYWORD_DATA_PAIR_HH
#define MAC_KEYWORD_DATA_PAIR_HH

#include <MAC_Object.hh>

#include <string>

class MAC_KeyItemPair ;
class MAC_Data ;
class MAC_String ;
class MAC_Vector ;

class doubleVector ;
class intVector ;

/*
Keyword-Value pairs.
*/

class MAC_KeywordDataPair : public MAC_Object
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_KeywordDataPair* create( MAC_Object* a_owner,
                                          MAC_String const* a_keyword, 
                                          MAC_Data const* a_data ) ;

   //-- Access
      
      // keyword of the pair
      std::string const& keyword( void ) const ;

      // value of the pair
      MAC_Data const* data( void ) const ;

      // Make `a_data' be the value of the pair.
      void replace_data( MAC_Data const* a_data ) ;

   //-- Comparison

      virtual size_t hash_code( void ) const ;

      virtual bool comparable( MAC_Object const* other ) const ;

      virtual bool is_equal( MAC_Object const* other ) const ;

      virtual int three_way_comparison( MAC_Object const* other ) const ;

      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------

      MAC_KeywordDataPair( void ) ;
     ~MAC_KeywordDataPair( void ) ;
      MAC_KeywordDataPair( MAC_KeywordDataPair const& other ) ;
      MAC_KeywordDataPair const& operator=( MAC_KeywordDataPair const& other ) ;
      
     // Constructor : node with a single child
      MAC_KeywordDataPair( MAC_Object* a_owner,
                           MAC_String const* a_keyword, 
                           MAC_Data const* a_data ) ;
    	
      
   //-- Attributes

      MAC_KeyItemPair* pair ;
};

#endif 
