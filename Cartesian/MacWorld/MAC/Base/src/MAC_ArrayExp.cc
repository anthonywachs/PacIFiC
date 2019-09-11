#include <MAC_ArrayExp.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Sequence.hh>

#include <iostream>

MAC_ArrayExp const* MAC_ArrayExp::PROTOTYPE = new MAC_ArrayExp() ;

//----------------------------------------------------------------------
MAC_ArrayExp:: MAC_ArrayExp( void ) 
//----------------------------------------------------------------------
   : MAC_Expression( "array" )
   , RESULT_D2D( 0, 0 )
   , RESULT_I2D( 0, 0 )
   , RESULT_B2D( 0, 0 )
   , RESULT_S2D( 0, 0 )
   , RESULT_D3D( 0, 0, 0 )
   , RESULT_I3D( 0, 0, 0 )
{
   MAC_LABEL( "MAC_ArrayExp:: MAC_ArrayExp" ) ;
}

//----------------------------------------------------------------------
MAC_ArrayExp:: MAC_ArrayExp( MAC_Object* a_owner,
                             MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, "array", argument_list )
   , RESULT_D2D( 0, 0 )
   , RESULT_I2D( 0, 0 )
   , RESULT_B2D( 0, 0 )
   , RESULT_S2D( 0, 0 )
   , RESULT_D3D( 0, 0, 0 )
   , RESULT_I3D( 0, 0, 0 )
{
   MAC_LABEL( "MAC_ArrayExp:: MAC_ArrayExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ArrayExp:: ~MAC_ArrayExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArrayExp:: ~MAC_ArrayExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ArrayExp*
MAC_ArrayExp:: create_replica( MAC_Object* a_owner,
                               MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArrayExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_ArrayExp* result = new MAC_ArrayExp( a_owner, argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_ArrayExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "array(<list of IV|DV|SV|BV>)" ;
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_ArrayExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArrayExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()>0 ;
   bool prem = true ;
   Type k = Undefined ;
   for( size_t idx=0 ; result && idx<some_arguments->count() ; ++idx )
   {
      if( prem )
      {
         k = extract_arg( some_arguments, idx )->data_type() ;
         if( k!=IntVector &&  k!=DoubleVector &&
             k!=IntArray2D &&  k!=DoubleArray2D &&
             k!=BoolVector &&  k!=StringVector )
         {
            result = false ;
         }
         prem = false ;
      }
      else
      {
         
         result &= ( k==extract_arg( some_arguments, idx )->data_type() ) ;
      }
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_ArrayExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArrayExp:: data_type" ) ;
   MAC_Data::Type my_type = Undefined ;
   if( arg(0)->data_type()==IntVector )
   {
      my_type = IntArray2D ;
   }
   else if( arg(0)->data_type()==DoubleVector ) 
   {
      my_type = DoubleArray2D ;
   }
   else if( arg(0)->data_type()==IntArray2D ) 
   {
      my_type = IntArray3D ;
   }
   else if( arg(0)->data_type()==DoubleArray2D ) 
   {
      my_type = DoubleArray3D ;
   }
   else if( arg(0)->data_type()==BoolVector ) 
   {
      my_type = BoolArray2D ;
   }
   else if( arg(0)->data_type()==StringVector ) 
   {
      my_type = StringArray2D ;
   }
   else
   {
      MAC_Error::object()->raise_internal( "Bad type " ) ;
   }
   return my_type ;
}

//----------------------------------------------------------------------
doubleArray2D const&
MAC_ArrayExp:: to_double_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArrayExp:: to_double_array2D" ) ;
   MAC_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   
   size_t nb_col = 0 ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      nb_col = MAC::max( nb_col, arg(idx)->to_double_vector(ct).size() ) ;
   }
   RESULT_D2D.re_initialize( nb_arguments(), nb_col ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      doubleVector const& vec = arg(idx)->to_double_vector( ct ) ;
      for( size_t j=0 ; j<vec.size() ; ++j )
      {
         RESULT_D2D( idx, j ) = vec( j ) ;
      }
   }
   return( RESULT_D2D ) ;
}

//----------------------------------------------------------------------
intArray2D const&
MAC_ArrayExp:: to_int_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArrayExp:: to_int_array2D" ) ;
   MAC_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;
   
   size_t nb_col = 0 ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      nb_col = MAC::max( nb_col, arg(idx)->to_int_vector(ct).size() ) ;
   }
   RESULT_I2D.re_initialize( nb_arguments(), nb_col ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      intVector const& vec = arg(idx)->to_int_vector( ct ) ;
      for( size_t j=0 ; j<vec.size() ; ++j )
      {
         RESULT_I2D( idx, j ) = vec( j ) ;
      }
   }
   return( RESULT_I2D ) ;
}


//----------------------------------------------------------------------
boolArray2D const&
MAC_ArrayExp:: to_bool_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArrayExp:: to_bool_array2D" ) ;
   MAC_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   
   size_t nb_col = 0 ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      nb_col = MAC::max( nb_col, arg(idx)->to_bool_vector(ct).size() ) ;
   }
   RESULT_B2D.re_initialize( nb_arguments(), nb_col ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      boolVector const& vec = arg(idx)->to_bool_vector( ct ) ;
      for( size_t j=0 ; j<vec.size() ; ++j )
      {
         RESULT_B2D( idx, j ) = vec( j ) ;
      }
   }
   return( RESULT_B2D ) ;
}


//----------------------------------------------------------------------
stringArray2D const&
MAC_ArrayExp:: to_string_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArrayExp:: to_string_array2D" ) ;
   MAC_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;
   
   size_t nb_col = 0 ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      nb_col = MAC::max( nb_col, arg(idx)->to_string_vector(ct).size() ) ;
   }
   RESULT_S2D.re_initialize( nb_arguments(), nb_col ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      stringVector const& vec = arg(idx)->to_string_vector( ct ) ;
      for( size_t j=0 ; j<vec.size() ; ++j )
      {
         RESULT_S2D( idx, j ) = vec( j ) ;
      }
   }
   return( RESULT_S2D ) ;
}

//----------------------------------------------------------------------
doubleArray3D const&
MAC_ArrayExp:: to_double_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArrayExp:: to_double_array3D" ) ;
   MAC_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   
   size_t dim1 = 0 ;
   size_t dim2 = 0 ;

   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      doubleArray2D const& a = arg(idx)->to_double_array2D(ct) ;
      dim1 = MAC::max( dim1, a.index_bound(0) ) ;
      dim2 = MAC::max( dim2, a.index_bound(1) ) ;
   }
   RESULT_D3D.re_initialize( nb_arguments(), dim1, dim2 ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      doubleArray2D const& a = arg(idx)->to_double_array2D(ct) ;
      for( size_t j=0 ; j<a.index_bound(0) ; ++j )
      {
         for( size_t k=0 ; k<a.index_bound(1) ; ++k )
         {
            RESULT_D3D( idx, j, k ) = a( j, k ) ;
         }
      }
   }
   return( RESULT_D3D ) ;
}

//----------------------------------------------------------------------
intArray3D const&
MAC_ArrayExp:: to_int_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ArrayExp:: to_int_array3D" ) ;
   MAC_CHECK_PRE( to_int_array3D_PRE( ct ) ) ;
   
   size_t dim1 = 0 ;
   size_t dim2 = 0 ;
   
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      intArray2D const& a = arg(idx)->to_int_array2D(ct) ;
      dim1 = MAC::max( dim1, a.index_bound(0) ) ;
      dim2 = MAC::max( dim2, a.index_bound(1) ) ;
   }
   RESULT_I3D.re_initialize( nb_arguments(), dim1, dim2 ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      intArray2D const& a = arg(idx)->to_int_array2D(ct) ;
      for( size_t j=0 ; j<a.index_bound(0) ; j++ )
      {
         for( size_t k=0 ; k<a.index_bound(1) ; k++ )
         {
            RESULT_I3D( idx, j, k ) = a( j, k ) ;
         }
      }
   }
   return( RESULT_I3D ) ;
}
