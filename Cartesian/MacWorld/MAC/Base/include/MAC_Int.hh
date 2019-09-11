#ifndef MAC_Int_HH
#define MAC_Int_HH

#include <MAC_Data.hh>

class MAC_Int : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization
      
      static MAC_Int* create( MAC_Object* a_owner, int val ) ;

      virtual MAC_Int* create_clone( MAC_Object* a_owner ) const ;
      
   //-- Type 
      
      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value
      
      virtual int to_int( MAC_Context const* ct = 0 ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   //-- Modifier
      
      void set( int val ) ;      
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_Int( void ) ;
     ~MAC_Int( void ) ;
      MAC_Int( MAC_Int const & other ) ;
      MAC_Int& operator=( MAC_Int const& other ) ;

      MAC_Int( MAC_Object* a_owner, int val ) ;
            
      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      int myValue ;
      
} ;


#endif
