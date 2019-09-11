#ifndef MAC_Bool_HH
#define MAC_Bool_HH

#include <MAC_Data.hh>

class MAC_Bool : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization
      
      static MAC_Bool* create( MAC_Object* a_owner,
                               bool val ) ;

      virtual MAC_Bool* create_clone( MAC_Object* a_owner ) const ;

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value
      
      virtual bool to_bool( MAC_Context const* ct = 0 ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Modifier
      
      void set( bool val ) ;      
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_Bool( void ) ;
     ~MAC_Bool( void ) ;
      MAC_Bool( MAC_Bool const & other ) ;
      MAC_Bool& operator=( MAC_Bool const& other ) ;

      MAC_Bool( MAC_Object* a_owner, bool val ) ;
            
      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      bool myValue ;
      
} ;


#endif
