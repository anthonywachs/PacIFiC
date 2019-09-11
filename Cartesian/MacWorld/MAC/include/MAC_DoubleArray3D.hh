#ifndef MAC_DOUBLEARRAY3D_HH
#define MAC_DOUBLEARRAY3D_HH

#include <MAC_Data.hh>

#include <doubleArray3D.hh>

class MAC_DoubleArray3D : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_DoubleArray3D* create( MAC_Object* a_owner,
                                      doubleArray3D const& val ) ;

      virtual MAC_DoubleArray3D* create_clone( MAC_Object* a_owner ) const ;

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value
      
      virtual doubleArray3D const& to_double_array3D( MAC_Context const* ct = 0 ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_DoubleArray3D( void ) ;
     ~MAC_DoubleArray3D( void ) ;
      MAC_DoubleArray3D( MAC_DoubleArray3D const & other ) ;
      MAC_DoubleArray3D& operator=( MAC_DoubleArray3D const& other ) ;

      MAC_DoubleArray3D( MAC_Object* a_owner, doubleArray3D const& val ) ;
            
   //-- Attributes

      doubleArray3D myValue ;
      
} ;


#endif
