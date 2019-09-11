#ifndef MAC_DATA_WITH_CONTEXT_HH
#define MAC_DATA_WITH_CONTEXT_HH

#include <MAC_Data.hh>

class MAC_Context ;
class MAC_ContextPair ;

class MAC_DataWithContext : public MAC_Data
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization
      
      static MAC_DataWithContext* create( MAC_Object* a_owner,
                                          MAC_Data const* data,
                                          MAC_Context const* ct ) ;

      virtual MAC_DataWithContext* create_clone( MAC_Object* a_owner ) const ;
      
   //-- Context
                  
      virtual void declare( MAC_List* lst ) const ;

   //-- Type

      virtual Type data_type( void ) const ;
      
   //-- Value
      
      MAC_Context const* context( MAC_Context const* ct = 0 ) const ;
      
      virtual bool value_can_be_evaluated( MAC_Context const* ct = 0 ) const ;
      
      virtual stringVector const& undefined_variables(
                                           MAC_Context const* ct = 0 ) const ;

      virtual bool to_bool( MAC_Context const* ct = 0 ) const ;

      virtual double to_double( MAC_Context const* ct = 0 ) const ;

      virtual int to_int( MAC_Context const* ct = 0 ) const ;

      virtual std::string const& to_string(
                                        MAC_Context const* ct = 0 ) const ;

      virtual doubleVector const& to_double_vector(
                                        MAC_Context const* ct = 0 ) const ;

      virtual intVector const& to_int_vector(
                                        MAC_Context const* ct = 0 ) const ;

      virtual stringVector const& to_string_vector(
                                        MAC_Context const* ct = 0 ) const ;

      virtual boolVector const& to_bool_vector(
                                        MAC_Context const* ct = 0 ) const ;

      virtual doubleArray2D const& to_double_array2D(
                                        MAC_Context const* ct = 0 ) const ;

      virtual intArray2D const& to_int_array2D(
                                        MAC_Context const* ct = 0 ) const ;

      virtual boolArray2D const& to_bool_array2D( 
                                        MAC_Context const* ct = 0 ) const ;
      
      virtual stringArray2D const& to_string_array2D( 
                                        MAC_Context const* ct = 0 ) const ;
      
      virtual doubleArray3D const& to_double_array3D(
                                        MAC_Context const* ct = 0 ) const ;

      virtual intArray3D const& to_int_array3D(
                                        MAC_Context const* ct = 0 ) const ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Formal calculus

      virtual bool is_raw_data( void ) const ;

   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------
      
      MAC_DataWithContext( void ) ;
     ~MAC_DataWithContext( void ) ;
      MAC_DataWithContext( MAC_DataWithContext const& other ) ;
      MAC_DataWithContext& operator=( MAC_DataWithContext const& other ) ;

      MAC_DataWithContext( MAC_Object* a_owner,
                           MAC_Data const* data,
                           MAC_Context const* ct ) ;

   //-- Attributes

      MAC_Data const* DATA ;
      MAC_Context* CTX ;
      MAC_ContextPair* TMP_CTX ;
      
      
};

#endif
