#ifndef MAC_TRANSFER_EXP_HH
#define MAC_TRANSFER_EXP_HH

#include <MAC_Expression.hh>

/* Abstract class provide to transfer evaluation to context-dependant
   data.
   
   FRAMEWORK INSTANTIATION :
      Implement `::data' method.

PUBLISHED
*/

class MAC_TransferExp : public MAC_Expression
{
   public: //----------------------------------------------------------
      
   //-- Value
      
      virtual bool to_bool( MAC_Context const* ct ) const ;
      virtual double to_double( MAC_Context const* ct ) const ;
      virtual int to_int( MAC_Context const* ct ) const ;
      virtual std::string const& to_string( MAC_Context const* ct ) const ;
      virtual doubleVector const& to_double_vector( MAC_Context const* ct ) const ;
      virtual intVector const& to_int_vector( MAC_Context const* ct ) const ;
      virtual stringVector const& to_string_vector(
                                              MAC_Context const* ct ) const ;
      virtual boolVector const& to_bool_vector( MAC_Context const* ct ) const ;
      virtual doubleArray2D const& to_double_array2D( MAC_Context const* ct ) const ;
      virtual intArray2D const& to_int_array2D( MAC_Context const* ct ) const ;
      virtual boolArray2D const& to_bool_array2D( 
                                        MAC_Context const* ct = 0 ) const ;      
      virtual stringArray2D const& to_string_array2D( 
                                        MAC_Context const* ct = 0 ) const ;
      virtual doubleArray3D const& to_double_array3D( MAC_Context const* ct ) const ;
      virtual intArray3D const& to_int_array3D( MAC_Context const* ct ) const ;
      
   protected: //-------------------------------------------------------
      
      virtual ~MAC_TransferExp( void ) ;

      MAC_TransferExp( std::string const& a_name ) ;

      MAC_TransferExp( MAC_Object* a_owner,
                      std::string const& a_name,
                      MAC_Sequence const* argument_list ) ;
     
   //-- Transfer implementation
      
      // `MAC_Data::' object depending on context `ct'
      virtual MAC_Data const* data( MAC_Context const* ct ) const = 0 ;
     
   //-- Preconditions, Postconditions, Invariant

      virtual bool data_PRE( MAC_Context const* ct ) const ;
      
      virtual bool data_POST( MAC_Data const* result,
                              MAC_Context const* ct ) const ;
      
   private: //-------------------------------------------------------

      MAC_TransferExp( void ) ;
      MAC_TransferExp( MAC_TransferExp const& other ) ;
      MAC_TransferExp& operator=( MAC_TransferExp const& other ) ;
      
};

#endif
