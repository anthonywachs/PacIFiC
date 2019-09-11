#ifndef MAC_DOUBLEARRAY2D_HH
#define MAC_DOUBLEARRAY2D_HH

#include <MAC_Data.hh>

#include <doubleArray2D.hh>

class MAC_DoubleArray2D : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_DoubleArray2D* create( MAC_Object* a_owner,
                                      doubleArray2D const& val ) ;

      // Build array from list of vectors.
      // Returns 0 if entry vectors are inconsistent.
      static MAC_DoubleArray2D* create( MAC_Object* a_owner,
                                        MAC_List const* list ) ;
      
      virtual MAC_DoubleArray2D* create_clone( MAC_Object* a_owner ) const ;

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value
      
      virtual doubleArray2D const& to_double_array2D( MAC_Context const* ct = 0 ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_DoubleArray2D( void ) ;
     ~MAC_DoubleArray2D( void ) ;
      MAC_DoubleArray2D( MAC_DoubleArray2D const & other ) ;
      MAC_DoubleArray2D& operator=( MAC_DoubleArray2D const& other ) ;

      MAC_DoubleArray2D( MAC_Object* a_owner, doubleArray2D const& val ) ;
            
   //-- Attributes

      doubleArray2D myValue ;
      
} ;


#endif
