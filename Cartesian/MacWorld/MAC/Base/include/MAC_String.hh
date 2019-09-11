#ifndef MAC_String_HH
#define MAC_String_HH

#include <MAC_Data.hh>

#include <string>

class MAC_String : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization
      
      static MAC_String* create( MAC_Object* a_owner, std::string const& a ) ;
      
      virtual MAC_String* create_clone( MAC_Object* a_owner ) const ;

   //-- Type
      
       virtual MAC_Data::Type data_type( void ) const ;

   //-- Value
       
      virtual std::string const& to_string( MAC_Context const* ct = 0 ) const ;

   //-- Comparison

      virtual bool is_equal( MAC_Object const* other ) const ;
      
      virtual int three_way_comparison( MAC_Object const* other ) const ;

      virtual size_t hash_code( void ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   //-- Modifier
      
      void set( std::string const& val ) ;
      
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_String( void ) ;
     ~MAC_String( void ) ;
      MAC_String( MAC_String const& other ) ;
      MAC_String const& operator=( MAC_String const& other ) ;

      MAC_String( MAC_Object* a_owner, std::string const& a ) ;
            
   //--- Attributes

      std::string str ;
} ;

#endif
