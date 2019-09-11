#ifndef MAC_IntArray2D_HH
#define MAC_IntArray2D_HH

#include <MAC_Data.hh>

#include <intArray2D.hh>

class size_t_array2D ;

class MAC_IntArray2D : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_IntArray2D* create( MAC_Object* a_owner,
                                     intArray2D const& val ) ;

      static MAC_IntArray2D* create( MAC_Object* a_owner,
                                     size_t_array2D const& val ) ;

      // Build array from list of vectors.
      // Returns 0 if entry vectors are inconsistent.
      static MAC_IntArray2D* create( MAC_Object* a_owner,
                                        MAC_List const* list ) ;
      
      virtual MAC_IntArray2D* create_clone( MAC_Object* a_owner ) const ;

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value
      
      virtual intArray2D const& to_int_array2D( MAC_Context const* ct = 0 ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_IntArray2D( void ) ;
     ~MAC_IntArray2D( void ) ;
      MAC_IntArray2D( MAC_IntArray2D const & other ) ;
      MAC_IntArray2D& operator=( MAC_IntArray2D const& other ) ;

      MAC_IntArray2D( MAC_Object* a_owner, intArray2D const& val ) ;

      MAC_IntArray2D( MAC_Object* a_owner, size_t_array2D const& val ) ;
            
   //-- Attributes

      intArray2D myValue ;
      
} ;


#endif
