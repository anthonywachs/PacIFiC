#ifndef MAC_STRING_ARRAY_2D_HH
#define MAC_STRING_ARRAY_2D_HH

#include <MAC_Data.hh>

#include <stringArray2D.hh>

class MAC_StringArray2D : public MAC_Data
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_StringArray2D* create( MAC_Object* a_owner,
                                        stringArray2D const& val ) ;

      // Build array from a list of vectors.
      // Returns 0 if entry vectors are inconsistent.
      static MAC_StringArray2D* create( MAC_Object* a_owner,
                                        MAC_List const* list ) ;
      
      virtual MAC_StringArray2D* create_clone( MAC_Object* a_owner ) const ;

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value
      
      virtual stringArray2D const& to_string_array2D( 
                                      MAC_Context const* ct = 0 ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_StringArray2D( void ) ;
     ~MAC_StringArray2D( void ) ;
      MAC_StringArray2D( MAC_StringArray2D const & other ) ;
      MAC_StringArray2D& operator=( MAC_StringArray2D const& other ) ;

      MAC_StringArray2D( MAC_Object* a_owner, stringArray2D const& val ) ;
            
   //-- Attributes

      stringArray2D MY_VALUE ;
} ;


#endif
