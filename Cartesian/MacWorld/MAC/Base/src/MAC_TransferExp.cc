#include <MAC_TransferExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Context.hh>
#include <MAC_Sequence.hh>


//----------------------------------------------------------------------
MAC_TransferExp:: MAC_TransferExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
{
}




//----------------------------------------------------------------------
MAC_TransferExp:: MAC_TransferExp( MAC_Object* a_owner,
                                   std::string const& a_name,
                                   MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
{
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
MAC_TransferExp:: ~MAC_TransferExp( void ) 
//----------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
bool
MAC_TransferExp:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_bool" ) ;
   MAC_CHECK_PRE( to_bool_PRE( ct ) ) ;
   return( data( ct )->to_bool( ct ) ) ;
}




//----------------------------------------------------------------------
double
MAC_TransferExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE( ct ) ) ;
   return( data( ct )->to_double( ct ) ) ;
}




//----------------------------------------------------------------------
int
MAC_TransferExp:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE( ct ) ) ;
   return( data( ct )->to_int( ct ) ) ;
}




//----------------------------------------------------------------------
std::string const&
MAC_TransferExp:: to_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_string" ) ;
   MAC_CHECK_PRE( to_string_PRE( ct ) ) ;
   return( data( ct )->to_string( ct ) ) ;
}




//----------------------------------------------------------------------
doubleVector const& 
MAC_TransferExp:: to_double_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_double_vector" ) ;
   MAC_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   return( data( ct )->to_double_vector( ct ) ) ;
}




//----------------------------------------------------------------------
intVector const&
MAC_TransferExp:: to_int_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_int_vector" ) ;
   MAC_CHECK_PRE( to_int_vector_PRE( ct ) ) ;
   return( data( ct )->to_int_vector( ct ) ) ;
}




//----------------------------------------------------------------------
stringVector const& 
MAC_TransferExp:: to_string_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_string_vector" ) ;
   MAC_CHECK_PRE( to_string_vector_PRE( ct ) ) ;
   return( data( ct )->to_string_vector( ct ) ) ;
}




//----------------------------------------------------------------------
boolVector const&
MAC_TransferExp:: to_bool_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_bool_vector" ) ;
   MAC_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;
   return( data( ct )->to_bool_vector( ct ) ) ;
}




//----------------------------------------------------------------------
doubleArray2D const&
MAC_TransferExp:: to_double_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_double_array2D" ) ;
   MAC_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   return( data( ct )->to_double_array2D( ct ) ) ;
}




//----------------------------------------------------------------------
boolArray2D const&
MAC_TransferExp:: to_bool_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_bool_array2D" ) ;
   MAC_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   return( data( ct )->to_bool_array2D( ct ) ) ;
}




//----------------------------------------------------------------------
stringArray2D const&
MAC_TransferExp:: to_string_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_string_array2D" ) ;
   MAC_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;
   return( data( ct )->to_string_array2D( ct ) ) ;
}




//----------------------------------------------------------------------
intArray2D const&
MAC_TransferExp:: to_int_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_int_array2D" ) ;
   MAC_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;
   return( data( ct )->to_int_array2D( ct ) ) ;
}




//----------------------------------------------------------------------
doubleArray3D const&
MAC_TransferExp:: to_double_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_double_array3D" ) ;
   MAC_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   return( data( ct )->to_double_array3D( ct ) ) ;
}




//----------------------------------------------------------------------
intArray3D const&
MAC_TransferExp:: to_int_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TransferExp:: to_int_array3D" ) ;
   MAC_CHECK_PRE( to_int_array3D_PRE( ct ) ) ;
   return( data( ct )->to_int_array3D( ct ) ) ;
}




//----------------------------------------------------------------------
bool
MAC_TransferExp:: data_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( !is_a_prototype() ) ;
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_TransferExp:: data_POST( MAC_Data const* result,
                             MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->value_can_be_evaluated( ct ) ) ;
   return( true ) ;
}
