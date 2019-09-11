#include <doubleVector.hh>

#ifdef OUTLINE
#define inline
#include <doubleVector.icc>
#undef inline
#endif

#include <MAC.hh>
#include <MAC_Error.hh>
#include <MAC_System.hh>

#include <iostream>

struct doubleVector_ERROR
{
   static void n0( void ) ;
} ;



//----------------------------------------------------------------------
doubleVector:: doubleVector( size_t dim, double val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   MAC_LABEL( "doubleVector:: doubleVector( size_t )" ) ;
   
   allocate( LENGTH ) ;
   set( val ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == dim ) ;
   MAC_CHECK_POST( capacity() == dim ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}




//----------------------------------------------------------------------
doubleVector:: doubleVector( int dim, double val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( static_cast<size_t>( dim ) )
   , CAPACITY( 0 )
{
   MAC_LABEL( "doubleVector:: doubleVector( int )" ) ;
   MAC_CHECK_PRE( dim >= 0 ) ;

   allocate( LENGTH ) ;
   set( val ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == static_cast<size_t>(dim) ) ;
   MAC_CHECK_POST( capacity() == size() ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}




//----------------------------------------------------------------------
doubleVector:: ~doubleVector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: ~doubleVector" ) ;
   MAC_CHECK_INV( invariant() ) ;
   deallocate() ;
}




//----------------------------------------------------------------------
doubleVector:: doubleVector( doubleVector const& other )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( other.LENGTH )
   , CAPACITY( 0 )
{
   MAC_LABEL( "doubleVector:: doubleVector( doubleVector const& )" ) ;

   allocate( other.capacity() ) ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      VECTOR[i] = other.VECTOR[i] ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == other.size() ) ;
   MAC_CHECK_POST( capacity() == other.capacity() ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == other(i) ) ) ;
}




//----------------------------------------------------------------------
doubleVector const&
doubleVector:: operator=( doubleVector const& other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: operator=" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   if( this != &other )
   {
      LENGTH = other.LENGTH ;
      if( LENGTH > CAPACITY )
      {
         deallocate() ;
         allocate( LENGTH ) ;
      }
      for( size_t i=0 ; i<LENGTH ; ++i )
      {
         VECTOR[i] = other.VECTOR[i] ;
      }
   }
   doubleVector const& result = *this ;
  
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( result.size() == other.size() ) ;
   MAC_CHECK_POST( IMPLIES( other.size() <= OLD(capacity),
                            result.capacity() == OLD(capacity) ) ) ;
   MAC_CHECK_POST( IMPLIES( other.size() > OLD(capacity),
                            result.capacity() == other.size() ) ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   result(i) == other(i) ) ) ;
   return( result ) ;
}




//---------------------------------------------------------------------
void
doubleVector:: re_initialize( size_t dim, double val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: re_initialize" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   LENGTH = dim ;
   if( CAPACITY != dim )
   {
      deallocate() ;
      allocate( LENGTH ) ;
   }
   set( val ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == dim ) ;
   MAC_CHECK_POST( capacity() == dim ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
            operator()(i) == val ) ) ;
}




//---------------------------------------------------------------------
void
doubleVector:: re_initialize( size_t ninter, double const& lower, 
      	double const& upper ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: re_initialize" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   LENGTH = ninter + 1 ;
   if( CAPACITY != ( ninter + 1 ) )
   {
      deallocate() ;
      allocate( LENGTH ) ;
   }

   VECTOR[0] = lower;
   double dx = ( upper - lower ) / ninter ;
   for( size_t i=1 ; i<LENGTH-1 ; ++i )
   {
      VECTOR[i] = lower + i*dx;
   }
   VECTOR[LENGTH-1] = upper;          

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == ( ninter + 1 ) ) ;
   MAC_CHECK_POST( capacity() == ( ninter + 1 ) ) ;
}




//----------------------------------------------------------------------
bool
doubleVector:: operator==( doubleVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: operator==" ) ;

   bool result = ( this == &other ) ;
   if( !result )
   {
      result = ( LENGTH==other.LENGTH ) ;
      for( size_t i=0 ; result && i<size() ; ++i )
      {
         result = ( operator()(i)==other(i) ) ;
      }
   }

   MAC_CHECK_POST( IMPLIES( result, size()==other.size() ) ) ;
   MAC_CHECK_POST( !result || FORALL( ( size_t i=0 ; i<size() ; ++i ),
                                       operator()(i) == other(i) ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
doubleVector:: operator!=( doubleVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: operator!=" ) ;

   bool result = !operator==( other ) ;

   MAC_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
size_t
doubleVector:: capacity( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: capacity" ) ;
   
   size_t result = CAPACITY ;
   
   MAC_CHECK_POST( result >= size() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
size_t
doubleVector:: index_of( double val, double a_dbl_eps, double a_dbl_min ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: index_of" ) ;
   
   size_t result = MAC::bad_index() ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      if( MAC::double_equality( VECTOR[i], val, a_dbl_eps, a_dbl_min ) )
      {
         result = i ;
         break ;
      }
   }
   
   MAC_CHECK_POST( 
      result == MAC::bad_index() ||
      MAC::double_equality( operator()(result) , val, 
                            a_dbl_eps, a_dbl_min ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
doubleVector:: has( double val, double a_dbl_eps, double a_dbl_min ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: has" ) ;
   
   size_t ifound = index_of( val, a_dbl_eps, a_dbl_min ) ;
   bool result = ( ifound != MAC::bad_index() ) ;

   MAC_CHECK_POST(
      IMPLIES( 
         result, 
         MAC::double_equality( operator()(
                                  index_of( val, a_dbl_eps, a_dbl_min ) ), 
                               val, a_dbl_eps, a_dbl_min ) ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
doubleVector:: resize( size_t dim, double val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: resize" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   if( dim > CAPACITY ) 
   {
      size_t const newCapacity = MAC_System::new_block_size( CAPACITY, dim ) ;
      double const* oldVector = VECTOR ;
      VECTOR = 0 ;
      allocate( newCapacity ) ;
      if( oldVector!=0 )
      {
         for( size_t i=0 ; i<LENGTH ; ++i )
         {
            VECTOR[i] = oldVector[i] ;
         }
         delete[] oldVector ;
      }
   }
   if( dim > LENGTH )
   {
      for( size_t i=LENGTH ; i<dim ; ++i )
      {
         VECTOR[i] = val ;
      }
   }
   LENGTH = dim ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( IMPLIES( dim <= OLD(capacity) ,
                            capacity() == OLD(capacity) ) ) ;
   MAC_CHECK_POST( 
      IMPLIES( dim > OLD(capacity) ,
               capacity() == MAC_System::new_block_size( OLD(capacity),dim ))) ;
   MAC_CHECK_POST( size() == dim ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=OLD(size) ; i<size() ; ++i ),
                           operator()(i) == val ) ) ;
}




//----------------------------------------------------------------------
void
doubleVector:: set( double val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: set" ) ;
   MAC_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      VECTOR[i] = val ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}




//----------------------------------------------------------------------
void
doubleVector:: append( double val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: append" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   size_t n = LENGTH ;
   resize( n+1 ) ;
   operator()( n ) = val ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( 
      IMPLIES( OLD(size)+1 <= OLD(capacity),
               capacity() == OLD( capacity ) ) ) ;
   MAC_CHECK_POST( 
      IMPLIES( OLD(size)+1 > OLD(capacity),
               capacity() == MAC_System::new_block_size( OLD(capacity),
                                                         OLD(size)+1 ) ) ) ;
   MAC_CHECK_POST( size() == OLD(size)+1 ) ;
   MAC_CHECK_POST( operator()(size()-1) == val ) ;
}




//----------------------------------------------------------------------
void
doubleVector:: remove_at( size_t idx )
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleVector:: remove_at" ) ;
   MAC_CHECK_PRE( idx<size() ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   for( size_t j=idx ; j<size()-1 ; j++ )
      VECTOR[j]=VECTOR[j+1] ;
   LENGTH-- ;

   MAC_CHECK_POST( size()==OLD(size)-1 ) ;
   MAC_CHECK_POST( capacity() == OLD(capacity) ) ;
}




//----------------------------------------------------------------------
std::ostream& operator<<( std::ostream& out, 
                          doubleVector const& vec ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, doubleVector const& )" ) ;

   size_t j=0 ;
   out << "< " ;
   for( size_t i=0 ; i<vec.size() ; ++i )
   {
      if( ++j>5 )
      {
         j = 1 ;
         out << std::endl << "  " ;
      }
      MAC::print_double( out, vec(i) ) ;
      out << " " ;
   }
   out << ">" ;
   return out ;
}




//----------------------------------------------------------------------
std::istream&
operator>>( std::istream& in, doubleVector& vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator>>( std::ostream&, doubleVector const& )" ) ;

   doubleVector aux( 0 ) ;
   
   std::string c ;
   double val = 0. ;
   
   if( !( in >> c ) || c != "<" )
   {
      doubleVector_ERROR::n0() ;
   }
   while( in >> val )
   {
      aux.append( val ) ;
   }
   in.clear() ;
   if( !( in >> c ) || c != ">" )
   {
      doubleVector_ERROR::n0() ;
   }

   vec = aux ;
   
   return( in ) ;
}




//----------------------------------------------------------------------
void
doubleVector:: allocate( size_t a_size )
//----------------------------------------------------------------------
{
   MAC_CHECK( VECTOR == 0 ) ;

   if( a_size > 0 )
   {
      VECTOR = new double [a_size] ;
      MAC_CHECK( VECTOR != 0 ) ;
   }
   CAPACITY = a_size ;
}




//----------------------------------------------------------------------
void
doubleVector:: deallocate( void )
//----------------------------------------------------------------------
{
   if( VECTOR != 0 )
   {
      delete [] VECTOR ;
      VECTOR = 0 ;
   }
   CAPACITY = 0 ;
}




//-----------------------------------------------------------------------
bool
doubleVector:: invariant( void ) const
//-----------------------------------------------------------------------
{
   MAC_ASSERT( LENGTH==0 || VECTOR != 0 ) ;
   MAC_ASSERT( LENGTH<=CAPACITY ) ;
   return( true ) ;
}




//internal--------------------------------------------------------------
void 
doubleVector_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   MAC_Error::object()->raise_plain(
      "*** Syntax error reading vector of \"double\":\n"
      "    < val1 val2 ... > is expected." ) ;
}

