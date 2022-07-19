#include <MAC_Communicator.hh>

#include <MAC.hh>
#include <MAC_DoubleComparator.hh>
#include <MAC_Error.hh>
#include <MAC_NumberedDoubleVectors.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_Root.hh>
#include <MAC_assertions.hh>

#include <boolArray2D.hh>
#include <boolVector.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <intArray2D.hh>
#include <intVector.hh>
#include <longLongIntVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <stringVector.hh>

#include <string>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;

struct MAC_Communicator_ERROR
{
   static void n0( void ) ;
   static void n1( void ) ;
   static void n2( void ) ;
   static void n3( void ) ;
} ;

bool MAC_Communicator:: TRACE = false ;

//----------------------------------------------------------------------
MAC_Communicator:: MAC_Communicator( std::string const& a_name )
//----------------------------------------------------------------------
   : MAC_Object( plugins_map() )
   , MY_NAME( a_name )
{
   MAC_LABEL( "MAC_Communicator:: MAC_Communicator" ) ;

   plugins_map()->register_item( a_name, this ) ;

   MAC_CHECK_POST( is_under_ownership_of( MAC_Root::object() ) ) ;
}

//----------------------------------------------------------------------
MAC_Communicator:: ~MAC_Communicator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
MAC_Communicator const*
MAC_Communicator:: object( std::string const& a_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: object" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;

   MAC_Communicator const* result =
      static_cast<MAC_Communicator const*>(
                               plugins_map()->item( a_name ) ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->is_under_ownership_of( MAC_Root::object() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_Communicator:: name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: name" ) ;
   return MY_NAME ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, size_t value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send size_t" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int val = (int) value ;
   send( dest, &val, 1 ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, size_t& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive size_t" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int val = MAC::bad_int() ;
   receive( src, &val, 1 ) ;
   if( val != (int) MAC::bad_index() && val < 0 )
      MAC_Communicator_ERROR:: n2() ;
   value = (size_t) val ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, size_t_vector const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send size_t_vector" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;
   MAC_CHECK_PRE( FORALL( ( size_t i=0 ; i<value.size() ; ++i ),
         value(i)==MAC::bad_index() || (int) value(i)<MAC::max_int() ) ) ;

   if( trace() ) print_method() ;

   intVector dummy(0) ;
   convert( value, dummy ) ;
   send( dest, dummy ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, size_t_vector& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive size_t_vector" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   intVector dummy(0) ;
   receive( src, dummy ) ;
   convert( dummy, value ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, size_t_array2D const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send size_t_array2D" ) ;
   MAC_CHECK_PRE( dest<nb_ranks() ) ;
   MAC_CHECK_PRE( dest!=rank() ) ;
   MAC_CHECK_PRE(
      FORALL( ( size_t i=0 ; i<value.index_bound(0) ; ++i ),
         FORALL(
            ( size_t j=0 ; j<value.index_bound(1) ; ++j ),
            value(i,j)==MAC::bad_index() || (int) value(i,j)<MAC::max_int() ) ) ) ;

   if( trace() ) print_method() ;

   intArray2D dummy(0,0) ;
   convert( value, dummy ) ;
   send( dest, dummy ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, size_t_array2D& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive size_t_array2D" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   intArray2D dummy(0,0) ;
   receive( src, dummy ) ;
   convert( dummy, value ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, int value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send int" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   send( dest, &value, 1 ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, int& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive int" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   receive( src, &value, 1 ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, intVector const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send intVector" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int dim = (int) value.size() ;
   send( dest, &dim, 1 ) ;
   if( dim>0 )
   {
      send( dest, value.data(), dim ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, longLongIntVector const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send longLongIntVector" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int dim = (int) value.size() ;
   send( dest, &dim, 1 ) ;
   if( dim>0 )
   {
      send( dest, value.data(), dim ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, intVector& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive intVector" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int dim = MAC::bad_int() ;
   receive( src, &dim, 1 ) ;
   if( dim<0 ) MAC_Communicator_ERROR::n0() ;
   value.re_initialize( (size_t) dim ) ;
   if( dim>0 )
   {
      receive( src, const_cast<int*>( value.data() ), dim ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, longLongIntVector& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive longLongIntVector" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int dim = MAC::bad_int() ;
   receive( src, &dim, 1 ) ;
   if( dim<0 ) MAC_Communicator_ERROR::n0() ;
   value.re_initialize( (size_t) dim ) ;
   if( dim>0 )
   {
      receive( src, const_cast<long long int*>( value.data() ), dim ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, intArray2D const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send intArray2D" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int D[2] ;
   D[0] = (int) value.index_bound( 0 ) ;
   D[1] = (int) value.index_bound( 1 ) ;
   send( dest, D, 2 ) ;
   if( D[0]*D[1]>0 )
   {
      send( dest, value.data(), D[0]*D[1] ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, intArray2D& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive intArray2D" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int D[2] = { MAC::bad_int(), MAC::bad_int() } ;
   receive( src, D, 2 ) ;
   if( D[0]<0 || D[1]<0 ) MAC_Communicator_ERROR::n0() ;
   value.re_initialize( (size_t) D[0], (size_t) D[1] ) ;
   if( D[0]*D[1]>0 )
   {
      receive( src, const_cast<int*>( value.data() ) , D[0]*D[1] ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, double value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send double" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   send( dest, &value, 1 ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, double& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive double" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   receive( src, &value, 1 ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, doubleVector const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send doubleVector" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int dim = (int) value.size() ;
   send( dest, &dim, 1 ) ;
   if( dim>0 )
   {
      send( dest, value.data(), dim ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, doubleVector& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive doubleVector" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int dim = MAC::bad_int() ;
   receive( src, &dim, 1 ) ;
   if( dim<0 ) MAC_Communicator_ERROR::n0() ;
   value.re_initialize( (size_t) dim ) ;
   if( dim>0 )
   {
      receive( src, const_cast<double*>( value.data() ), dim ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, doubleArray2D const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send doubleArray2D" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int D[2] ;
   D[0] = (int) value.index_bound( 0 ) ;
   D[1] = (int) value.index_bound( 1 ) ;
   send( dest, D, 2 ) ;
   if( D[0]*D[1]>0 )
   {
      send( dest, value.data(), D[0]*D[1] ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, doubleArray2D& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive doubleArray2D" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int D[2] = { MAC::bad_int(), MAC::bad_int() } ;
   receive( src, D, 2 ) ;
   if( D[0]<0 || D[1]<0 ) MAC_Communicator_ERROR::n0() ;
   value.re_initialize( (size_t) D[0], (size_t) D[1] ) ;
   if( D[0]*D[1]>0 )
   {
      receive( src, const_cast<double*>( value.data() ), D[0]*D[1] ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, bool value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send bool" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int i = ( value ? 1 : 0 ) ;
   send( dest, &i, 1 ) ;

}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, bool& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive bool" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int i ;
   receive( src, &i, 1 ) ;
   if( i != 0 && i != 1 ) MAC_Communicator_ERROR:: n3() ;
   value = ( i==1 ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, boolVector const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send boolVector" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   intVector dummy(0) ;
   convert( value, dummy ) ;
   send( dest, dummy ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, boolVector& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive boolVector" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   intVector dummy(0) ;
   receive( src, dummy ) ;
   convert( dummy, value ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, boolArray2D const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send boolArray2D" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   intArray2D dummy(0,0) ;
   convert( value, dummy ) ;
   send( dest, dummy ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, boolArray2D& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive boolArray2D" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   intArray2D dummy(0,0) ;
   receive( src, dummy ) ;
   convert( dummy, value ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, std::string const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send string" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int dim = (int) (value.size()+1) ;
   send( dest, &dim, 1 ) ;
   send( dest, value.c_str(), dim ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, std::string& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive string" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int dim = MAC::bad_int() ;
   receive( src, &dim, 1 ) ;
   if( dim<1 ) MAC_Communicator_ERROR::n0() ;
   char* val = new char[dim] ;
   receive( src, val, dim ) ;
   if( val[dim-1] != '\0' ) MAC_Communicator_ERROR::n1() ;
   value = std::string( val ) ;
   delete [] val ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: send( size_t dest, stringVector const& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: send stringVector" ) ;
   MAC_CHECK_PRE( dest < nb_ranks() ) ;
   MAC_CHECK_PRE( dest != rank() ) ;

   if( trace() ) print_method() ;

   int dim = (int) value.size() ;
   send( dest, &dim, 1 ) ;
   if( dim>0 )
   {
      for( size_t i=0 ; i<value.size() ; ++i )
      {
         send( dest, value(i) ) ;
      }
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: receive( size_t src, stringVector& value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: receive stringVector" ) ;
   MAC_CHECK_PRE( src < nb_ranks() ) ;
   MAC_CHECK_PRE( src != rank() ) ;

   if( trace() ) print_method() ;

   int dim = MAC::bad_int() ;
   receive( src, &dim, 1 ) ;
   if( dim<0 ) MAC_Communicator_ERROR::n0() ;
   value.re_initialize( dim ) ;
   if( dim>0 )
   {
      for( size_t i=0 ; i<value.size() ; ++i )
      {
         receive( src, value(i) ) ;
      }
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: wait( void* request ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( size_t& value, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(size_t)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      int dummy = ( rank() != root ? MAC::bad_int() : (int) value ) ;
      broadcast( &dummy, (size_t)1, root ) ;
      if( dummy != (int) MAC::bad_index() && dummy < 0 )
         MAC_Communicator_ERROR:: n2() ;
      value = (size_t) dummy ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( size_t_vector& values, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(size_t_vector)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      intVector dummy(0) ;
      if( rank() == root ) convert( values, dummy ) ;
      broadcast( dummy, root ) ;
      if( rank() != root ) convert( dummy, values ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( size_t_array2D& values, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(size_t_array2D)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      intArray2D dummy(0,0) ;
      if( rank() == root ) convert( values, dummy ) ;
      broadcast( dummy, root ) ;
      if( rank() != root ) convert( dummy, values ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( int& value, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(int)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   broadcast( &value, (size_t)1, root ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( intVector& values, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(intVector)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   int dim = ( rank() != root ? MAC::bad_int() : (int) values.size() ) ;
   broadcast( &dim , (size_t) 1, root ) ;
   if( rank() != root ) values.re_initialize( dim ) ;
   broadcast( const_cast<int*>( values.data() ), (size_t) dim, root ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( intArray2D& values, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(intArray2D)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   int dim0 = ( rank() != root ? MAC::bad_int() : (int) values.index_bound(0) ) ;
   int dim1 = ( rank() != root ? MAC::bad_int() : (int) values.index_bound(1) ) ;
   broadcast( &dim0, (size_t) 1, root ) ;
   broadcast( &dim1, (size_t) 1, root ) ;
   if( rank() != root ) values.re_initialize( dim0, dim1 ) ;
   broadcast( const_cast<int*>( values.data() ), (size_t) dim0*dim1, root ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( double& value, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(double)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   broadcast( &value, (size_t)1, root ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( doubleVector& values, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(doubleVector)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   int dim = ( rank() != root ? MAC::bad_int() : (int) values.size() ) ;
   broadcast( &dim , (size_t) 1, root ) ;
   if( rank() != root ) values.re_initialize( dim ) ;
   broadcast( const_cast<double*>( values.data() ), (size_t) dim, root ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( doubleArray2D& values, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(doubleArray2D)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   int dim0 = ( rank() != root ? MAC::bad_int() : (int) values.index_bound(0) ) ;
   int dim1 = ( rank() != root ? MAC::bad_int() : (int) values.index_bound(1) ) ;
   broadcast( &dim0, (size_t) 1, root ) ;
   broadcast( &dim1, (size_t) 1, root ) ;
   if( rank() != root ) values.re_initialize( dim0, dim1 ) ;
   broadcast( const_cast<double*>( values.data() ), (size_t) dim0*dim1, root ) ;
}


//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( bool& value, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(bool)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      int dummy = ( rank() != root ? MAC::bad_int() : ( value ? 1 : 0 ) ) ;
      broadcast( &dummy, (size_t)1, root ) ;
      if( dummy != 0 && dummy != 1 )
         MAC_Communicator_ERROR:: n2() ;
      value = ( dummy == 1 ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( boolVector& values, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(boolVector)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      intVector dummy(0) ;
      if( rank() == root ) convert( values, dummy ) ;
      broadcast( dummy, root ) ;
      if( rank() != root ) convert( dummy, values ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( boolArray2D& values, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(boolArray2D)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   if( nb_ranks() > 1 )
   {
      intArray2D dummy(0,0) ;
      if( rank() == root ) convert( values, dummy ) ;
      broadcast( dummy, root ) ;
      if( rank() != root ) convert( dummy, values ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( std::string& value, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(string)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   int dim = MAC::bad_int() ;
   if( rank() == root ) dim = int(value.size()+1) ;
   broadcast( dim, root ) ;

   char* tab = 0 ;
   if( rank() == root )
   {
      tab = const_cast<char*>( value.c_str() ) ;
   }
   else
   {
      tab = new char[dim] ;
   }
   broadcast( tab, dim, root ) ;
   if( rank() != root )
   {
      value = std::string( tab ) ;
      delete [] tab ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: broadcast( stringVector& value, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: broadcast(stringVector)" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;

   int dim = MAC::bad_int() ;
   if( rank() == root ) dim = (int) value.size() ;
   broadcast( dim, root ) ;
   if( rank() != root ) value.re_initialize( (size_t) dim ) ;
   for( size_t i=0 ; i<value.size() ; ++i )
   {
      broadcast( value(i), root ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: gather( double value, doubleVector& result,
                           size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: gather(double)" ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;
   MAC_CHECK_PRE( root!=rank() || result.size()==nb_ranks() ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   gather( &value, 1, const_cast<double*>( result.data() ), root ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: gather( int value, intVector& result,
                           size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: gather(double)" ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;
   MAC_CHECK_PRE( root!=rank() || result.size()==nb_ranks() ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   gather( &value, 1, const_cast<int*>( result.data() ), root ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: gather( doubleVector const& values,
                           doubleVector& result,
                           size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: gather(intVector)" ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;
   MAC_CHECK_PRE( root!=rank() || result.size()==nb_ranks()*values.size() ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   gather( values.data(), values.size(),
           const_cast<double*>( result.data() ), root ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: gather( intVector const& values,
                           intVector& result,
                           size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: gather(intVector)" ) ;
   MAC_CHECK_PRE( root < nb_ranks() ) ;
   MAC_CHECK_PRE( root!=rank() || result.size()==nb_ranks()*values.size() ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   gather( values.data(), values.size(),
           const_cast<int*>( result.data() ), root ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: boolean_and( bool value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: boolean_and" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   bool result = value ;
   if( rank()==0 )
   {
      bool recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result &= recep ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: boolean_or( bool value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: boolean_or" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   bool result = value ;
   if( rank()==0 )
   {
      bool recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result |= recep ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
double
MAC_Communicator:: sum( double value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: sum" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   double result = value ;
   if( rank()==0 )
   {
      double recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result += recep ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: sum_vector( doubleVector& vec ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: sum_vector" ) ;
   MAC_CHECK_PRE( same_value_everywhere( (int) vec.size() ) ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   MAC_SAVEOLD( size_t, vec_size, vec.size() ) ;

   if( trace() ) print_method() ;

   int const nb = (int) vec.size() ;
   if( nb_ranks()>1 && nb>0 )
   {
      sum_vector( const_cast<double*>( vec.data() ), nb ) ;
   }
   
   MAC_CHECK_POST( vec.size() == OLD( vec_size ) ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: sum_array( doubleArray2D& array ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: sum_array" ) ;
   MAC_CHECK_PRE( same_value_everywhere( (int) array.index_bound(0) ) ) ;
   MAC_CHECK_PRE( same_value_everywhere( (int) array.index_bound(1) ) ) ;
   MAC_SAVEOLD( size_t, array_index_bound0, array.index_bound(0) ) ;
   MAC_SAVEOLD( size_t, array_index_bound1, array.index_bound(1) ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   int const nb = (int) array.index_bound(0)*array.index_bound(1) ;
   if( nb_ranks()>1 && nb>0 )
   {
      sum_vector( const_cast<double*>( array.data() ), nb ) ;
   }

   MAC_CHECK_POST( array.index_bound(0) == OLD( array_index_bound0 ) ) ;
   MAC_CHECK_POST( array.index_bound(1) == OLD( array_index_bound1 ) ) ;
}


//----------------------------------------------------------------------
void
MAC_Communicator:: merge( MAC_DoubleComparator const* dbl_comp,
                          doubleArray2D& coord,
                          size_t_vector& idx ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: merge" ) ;
   MAC_CHECK_PRE( dbl_comp != 0 ) ;
   MAC_SAVEOLD( doubleArray2D, coord, coord ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   size_t dim = coord.index_bound( 0 ) ;

   size_t_array2D idx_all( 0, nb_ranks() ) ;
   size_t_vector nbvert( nb_ranks() ) ;
   doubleArray2D coords( 0, 0 ) ;

   if( rank()>0 )
   {
      send( 0, coord ) ;
      receive( 0, idx ) ;
   }
   else
   {
      MAC_NumberedDoubleVectors* numvec =
                  MAC_NumberedDoubleVectors::create( 0, dbl_comp, dim ) ;
      doubleVector vec( dim ) ;
      for( size_t i=0 ; i<nb_ranks() ; i++ )
      {
         if( i>0 )
         {
            receive( i, coords ) ;
         }
         else
         {
            coords = coord ;
         }
         MAC_CHECK( coords.index_bound( 0 ) == dim ) ;
         size_t n = nbvert( i ) = coords.index_bound( 1 ) ;
         if( idx_all.index_bound( 0 )<nbvert( i ) )
            idx_all.raise_first_index_bound( nbvert( i ) ) ;
         for( size_t j=0 ; j<n; j++ )
         {
            for( size_t d=0 ; d<dim; d++ )
               vec( d ) = coords( d, j ) ;
            numvec->extend( vec ) ;
            idx_all( j, i ) = numvec->index( vec ) ;
         }
      }

      // Re-ordering nodes to be independent from splitting
      size_t_vector const& perm = numvec->order() ;

      for( size_t i=nb_ranks()-1 ; i<nb_ranks() ; i-- )
      {
         idx.re_initialize( nbvert( i ) ) ;
         for( size_t j=0 ; j<idx.size() ; j++ )
         {
            idx( j ) = perm( idx_all( j, i ) ) ;
         }
         if( i>0 )
         {
            send( i, idx ) ;
         }
      }
      coord = numvec->ordered_items() ;
      numvec->destroy() ;
   }

   MAC_CHECK_POST( IMPLIES( rank()!=0, coord==OLD(coord) ) ) ;
   MAC_CHECK_POST( rank() != 0 ||
      FORALL( ( size_t i=0 ; i<coord.index_bound(1)-1 ; ++i ),
         dbl_comp->three_way_comparison( coord(0,i), coord(0,i+1) ) <= 0 ) ) ;
   MAC_CHECK_POST( coord.index_bound( 0 ) == OLD(coord).index_bound( 0 ) ) ;
   MAC_CHECK_POST( idx.size() == OLD(coord).index_bound( 1 ) ) ;
}

//----------------------------------------------------------------------
double
MAC_Communicator:: min( double value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: min" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   double result = value ;
   if( rank() == 0 )
   {
      double recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result = MAC::min( result, recep ) ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   MAC_CHECK_POST( min_POST( result, value ) ) ;
   return( result );
}

//----------------------------------------------------------------------
size_t
MAC_Communicator:: min( size_t value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: min" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   size_t result = value ;
   if( rank() == 0 )
   {
      size_t recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result = MAC::min( result, recep ) ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   MAC_CHECK_POST( min_POST( result, value ) ) ;
   return( result );
}

//----------------------------------------------------------------------
double
MAC_Communicator:: max( double value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: max" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   double result = value ;
   if( rank() == 0 )
   {
      double recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result = MAC::max( result, recep ) ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   MAC_CHECK_POST( max_POST( result, value ) ) ;
   return( result );
}

//----------------------------------------------------------------------
size_t
MAC_Communicator:: max( size_t value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: max" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   if( trace() ) print_method() ;

   size_t result = value ;
   if( rank() == 0 )
   {
      size_t recep ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep ) ;
         result = MAC::max( result, recep ) ;
      }
   }
   else
   {
      send( 0, value ) ;
   }
   broadcast( result, 0 ) ;

   MAC_CHECK_POST( max_POST( result, value ) ) ;
   return( result );
}

//----------------------------------------------------------------------
void
MAC_Communicator:: all_gather( int value, intVector& result ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: all_gather" ) ;
   MAC_CHECK_PRE( result.size()==nb_ranks() ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   all_gather( &value, 1, const_cast<int*>(result.data()) ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: all_gather( intVector const& values, intVector& result ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: all_gather" ) ;
   MAC_CHECK_PRE( values.size()!=0 ) ;
   MAC_CHECK_PRE( result.size()==nb_ranks()*values.size() ) ;
   MAC_CHECK_COLLECTIVE( true ) ;

   all_gather( values.data(), values.size(), const_cast<int*>(result.data()) ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: all_to_all( intVector const& values,
                               intVector& result ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: all_to_all" ) ;
   MAC_CHECK_PRE( values.size() % nb_ranks() == 0 ) ;
   MAC_CHECK_PRE( result.size()==values.size() ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   int N = values.size() / nb_ranks() ;

   all_to_all( values.data(), N, const_cast<int*>(result.data()) ) ;
}

//----------------------------------------------------------------------
MAC_ObjectRegister*
MAC_Communicator:: plugins_map( void )
//----------------------------------------------------------------------
{
   static MAC_ObjectRegister* result =
           MAC_ObjectRegister::create( MAC_Root::object(),
                                       "MAC_Communicator descendant" ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: trace( void )
//----------------------------------------------------------------------
{
   return( TRACE ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: print_method( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: print_method" ) ;
   if( MAC_Marker::nb_labels()>1 )
   {
      MAC::out() << "Methods : " ;
      for( size_t i=0 ; i<MAC_Marker::nb_labels()-1 ; ++i )
      {
         MAC::out() << MAC_Marker::label(i) << ";" ;
      }
      MAC::out() << std::endl ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: do_trace( void )
//----------------------------------------------------------------------
{
   TRACE = true ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: convert(
                       size_t_vector const& src, intVector& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.size() ) ;
   for( size_t i=0 ; i<src.size() ; ++i )
   {
      dest(i) = (int) src(i) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: convert(
                       intVector const& src, size_t_vector& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.size() ) ;
   for( size_t i=0 ; i<src.size() ; ++i )
   {
      if( src(i)!= (int) MAC::bad_index() && src(i) < 0 )
         MAC_Communicator_ERROR:: n2() ;
      dest(i) = (size_t) src(i) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: convert(
                     size_t_array2D const& src, intArray2D& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.index_bound(0), src.index_bound(1) ) ;
   for( size_t i=0 ; i<src.index_bound(0) ; ++i )
   {
      for( size_t j=0 ; j<src.index_bound(1) ; ++j )
      {
         dest(i,j) = (int) src(i,j) ;
      }
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: convert(
                     intArray2D const& src, size_t_array2D& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.index_bound(0), src.index_bound(1) ) ;
   for( size_t i=0 ; i<src.index_bound(0) ; ++i )
   {
      for( size_t j=0 ; j<src.index_bound(1) ; ++j )
      {
         if( src(i,j) != (int) MAC::bad_index() && src(i,j) < 0 )
            MAC_Communicator_ERROR:: n2() ;
         dest(i,j) = (size_t) src(i,j) ;
      }
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: convert(
                          boolVector const& src, intVector& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.size() ) ;
   for( size_t i=0 ; i<src.size() ; ++i )
   {
      dest(i) = ( src(i) ? 1 : 0 ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: convert(
                          intVector const& src, boolVector& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.size() ) ;
   for( size_t i=0 ; i<src.size() ; ++i )
   {
      if( src(i) != 0 && src(i) != 1 )
         MAC_Communicator_ERROR:: n3() ;
      dest(i) = ( src(i) == 1 ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: convert(
                        boolArray2D const& src, intArray2D& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.index_bound(0), src.index_bound(1) ) ;
   for( size_t i=0 ; i<src.index_bound(0) ; ++i )
   {
      for( size_t j=0 ; j<src.index_bound(1) ; ++j )
      {
         dest(i,j) = ( src(i,j) ? 1 : 0 ) ;
      }
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: convert(
                        intArray2D const& src, boolArray2D& dest ) const
//----------------------------------------------------------------------
{
   dest.re_initialize( src.index_bound(0), src.index_bound(1) ) ;
   for( size_t i=0 ; i<src.index_bound(0) ; ++i )
   {
      for( size_t j=0 ; j<src.index_bound(1) ; ++j )
      {
         if( src(i,j) != 0 && src(i,j) != 1 )
            MAC_Communicator_ERROR:: n3() ;
         dest(i,j) = ( src(i,j) == 1 ) ;
      }
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: sum_vector( double* values, int nb ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: sum_vector(double*)" ) ;
   MAC_CHECK_PRE( sum_vector_PRE( values, nb ) ) ;
   
   if( rank()==0 )
   {
      int s = nb ;
      double* recep = new double[nb] ;
      for( size_t i=1 ; i<nb_ranks() ; ++i )
      {
         receive( i, recep, s ) ;
         MAC_ASSERT( s == nb ) ;
         for( int j=0 ; j<nb ; ++j )
            values[j] += recep[j] ;
      }
      delete [] recep ;
   }
   else
   {
      send( 0, values, nb ) ;
   }
   broadcast( values, nb, 0 ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: nb_ranks_POST( size_t result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result > 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: rank_POST( size_t result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result < nb_ranks() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: send_PRE( size_t dest, int const* value, int nb ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( dest < nb_ranks() ) ;
   MAC_ASSERT( dest != rank() ) ;
   MAC_ASSERT( nb > 0 ) ;
   MAC_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: send_PRE( size_t dest, long long int const* value, int nb ) 
	const
//----------------------------------------------------------------------
{
   MAC_ASSERT( dest < nb_ranks() ) ;
   MAC_ASSERT( dest != rank() ) ;
   MAC_ASSERT( nb > 0 ) ;
   MAC_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: send_PRE( size_t dest, double const* value, int nb ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( dest < nb_ranks() ) ;
   MAC_ASSERT( dest != rank() ) ;
   MAC_ASSERT( nb > 0 ) ;
   MAC_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: send_PRE( size_t dest, char const* value, int nb ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( dest < nb_ranks() ) ;
   MAC_ASSERT( dest != rank() ) ;
   MAC_ASSERT( nb > 0 ) ;
   MAC_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: receive_PRE( size_t src, int const* value, int nb ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( src < nb_ranks() ) ;
   MAC_ASSERT( src != rank() ) ;
   MAC_ASSERT( nb > 0 ) ;
   MAC_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: receive_PRE( size_t src, long long int const* value, int nb )
	const
//----------------------------------------------------------------------
{
   MAC_ASSERT( src < nb_ranks() ) ;
   MAC_ASSERT( src != rank() ) ;
   MAC_ASSERT( nb > 0 ) ;
   MAC_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: receive_PRE( size_t src, double const* value, int nb ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( src < nb_ranks() ) ;
   MAC_ASSERT( src != rank() ) ;
   MAC_ASSERT( nb > 0 ) ;
   MAC_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: receive_PRE( size_t src, char const* value, int nb ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( src < nb_ranks() ) ;
   MAC_ASSERT( src != rank() ) ;
   MAC_ASSERT( nb > 0 ) ;
   MAC_ASSERT( value != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: min_POST( double result, double value ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result <= value ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: min_POST( size_t result, size_t value ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result <= value ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: sum_vector_PRE( double const* values, int nb ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( same_value_everywhere( nb ) ) ;
   MAC_ASSERT( nb > 0 ) ;
   MAC_ASSERT( nb_ranks() > 1 ) ;
   MAC_ASSERT( values != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: max_POST( double result, double value ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result >= value ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Communicator:: max_POST( size_t result, size_t value ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result >= value ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void
MAC_Communicator_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << std::endl
        << "*** MAC_Communicator:" << std::endl
        << "    negative table bound received" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
MAC_Communicator_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << std::endl
        << "*** MAC_Communicator:" << std::endl
        << "    pointer to a no null-terminated array of characters received." ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
MAC_Communicator_ERROR:: n2( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << std::endl
        << "*** MAC_Communicator:" << std::endl
        << "    negative size_t received." ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void
MAC_Communicator_ERROR:: n3( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << std::endl
        << "*** MAC_Communicator:" << std::endl
        << "    bad boolean value received." ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: all_gather( doubleVector const& values,
                           doubleVector& result) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: all_vgather(doubleVector)" ) ;
   MAC_CHECK_PRE(result.size()==nb_ranks()*values.size() ) ;
   
   all_gather( values.data(), values.size(),
           const_cast<double*>( result.data() )) ;
}

//----------------------------------------------------------------------
size_t
MAC_Communicator:: sum( size_t value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: sum" ) ;

   if( trace() ) print_method() ;
   
   size_t result = value ;
   if( rank()==0 )
   {
      size_t recep ;
      for( size_t i=1 ; i<nb_ranks() ; i++ )
      {
         receive( i, recep ) ;
         result += recep ;
      }
      for( size_t i=1 ; i<nb_ranks() ; i++ )
      {
         send( i, result ) ;
      }
   }
   else
   {
      send( 0, value ) ;
      receive( 0, result ) ;
   }
   return result ;
}

//----------------------------------------------------------------------
void
MAC_Communicator:: reduce_vector( doubleVector& vec, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: reduce_vector" ) ;

   if( trace() ) print_method() ;
   
   if( rank()==root )
   {
      doubleVector recep(vec) ;
      for( size_t i=1 ; i<nb_ranks() ; i++ )
      {
         receive( i, recep ) ;
         MAC_ASSERT( recep.size() == vec.size() ) ;
         for( size_t j=0 ; j<vec.size() ; j++ )
            vec(j) += recep(j) ;
      }
   }
   else
   {
      send( root, vec ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_Communicator:: reduce_vector_max( doubleVector& vec, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: reduce_vector" ) ;

   if( trace() ) print_method() ;
   
   if( rank()==root )
   {
      doubleVector recep(vec) ;
      for( size_t i=1 ; i<nb_ranks() ; i++ )
      {
         receive( i, recep ) ;
         MAC_ASSERT( recep.size() == vec.size() ) ;
         for( size_t j=0 ; j<vec.size() ; j++ )
            vec(j) = MAC::max( recep(j), vec(j) ) ;
      }
   }
   else
   {
      send( root, vec ) ;
   }
}

//----------------------------------------------------------------------
unsigned long long int
MAC_Communicator:: sum( unsigned long long int value ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Communicator:: sum" ) ;

   if( trace() ) print_method() ;
   
   unsigned long long int result = value ;
   return( result );
}
