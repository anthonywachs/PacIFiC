#ifndef MAC_BOOL_ARRAY_2D_HH
#define MAC_BOOL_ARRAY_2D_HH

#include <MAC_Data.hh>

#include <boolArray2D.hh>

class MAC_BoolArray2D : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_BoolArray2D* create( MAC_Object* a_owner,
                                      boolArray2D const& val ) ;

      // Build array from list of vectors.
      // Returns 0 if entry vectors are inconsistent.
      static MAC_BoolArray2D* create( MAC_Object* a_owner,
                                      MAC_List const* list ) ;
      
      virtual MAC_BoolArray2D* create_clone( MAC_Object* a_owner ) const ;

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value
      
      virtual boolArray2D const& to_bool_array2D(
                                         MAC_Context const* ct = 0 ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_BoolArray2D( void ) ;
     ~MAC_BoolArray2D( void ) ;
      MAC_BoolArray2D( MAC_BoolArray2D const & other ) ;
      MAC_BoolArray2D& operator=( MAC_BoolArray2D const& other ) ;

      MAC_BoolArray2D( MAC_Object* a_owner, boolArray2D const& val ) ;
            
   //-- Attributes

      boolArray2D MY_VALUE ;
} ;


#endif
