#ifndef MAC_DOUBLE_HH
#define MAC_DOUBLE_HH

#include <MAC_Data.hh>

class MAC_DoubleComparator ;

class MAC_Double : public MAC_Data
{
      
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_Double* create( MAC_Object* a_owner, double val ) ;

      virtual MAC_Double* create_clone( MAC_Object* a_owner ) const ;
      
   //-- Type

      virtual MAC_Data::Type data_type( void ) const ;

   //-- Value
      
      virtual double to_double( MAC_Context const* ct = 0 ) const ;
     
   //-- Comparison
      
      static MAC_DoubleComparator const* double_comparator( void ) ;
      
      virtual size_t hash_code( void ) const ;
      
      virtual bool is_equal( MAC_Object const* other ) const ;

      virtual int three_way_comparison( MAC_Object const* other ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Modifier
      
      void set( double val ) ;
     
   //-- Formal calculus

      virtual MAC_Double* create_derivative( MAC_Object* a_owner,
                                             MAC_Variable const* var,
                                             MAC_Context const* ct ) const ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      MAC_Double( void ) ;
     ~MAC_Double( void ) ;
      MAC_Double( MAC_Double const & other ) ;
      MAC_Double& operator=( MAC_Double const& other ) ;

      MAC_Double( MAC_Object* a_owner, double val ) ;
      
   //-- Attributes
      
      double VALUE ;
      
} ;


#endif
