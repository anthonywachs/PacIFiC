#ifndef MAC_IntArray3D_HH
#define MAC_IntArray3D_HH

#include <MAC_Data.hh>

#include <intArray3D.hh>

class MAC_IntArray3D : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_IntArray3D* create( MAC_Object* a_owner,
                                    intArray3D const& val ) ;

      virtual MAC_IntArray3D* create_clone( MAC_Object* a_owner ) const ;

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value
      
      virtual intArray3D const& to_int_array3D( MAC_Context const* ct = 0 ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_IntArray3D( void ) ;
     ~MAC_IntArray3D( void ) ;
      MAC_IntArray3D( MAC_IntArray3D const & other ) ;
      MAC_IntArray3D& operator=( MAC_IntArray3D const& other ) ;

      MAC_IntArray3D( MAC_Object* a_owner, intArray3D const& val ) ;
            
   //-- Attributes

      intArray3D myValue ;
      
} ;


#endif
