#include <intVector.hh>

#ifdef OUTLINE
#define inline
#include <intVector.icc>
#undef inline
#endif

#include <MAC.hh>
#include <MAC_Error.hh>
#include <MAC_System.hh>

#include <iostream>

struct intVector_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
intVector:: intVector( size_t dim, int val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   MAC_LABEL( "intVector:: intVector( size_t )" ) ;

   allocate( LENGTH ) ;
   set( val ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == dim ) ;
   MAC_CHECK_POST( capacity() == dim ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
intVector:: intVector( int dim, int val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( static_cast<size_t>(dim) )
   , CAPACITY( 0 )
{
   MAC_LABEL( "intVector:: intVector( int )" ) ;
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
intVector:: intVector( intVector const& other )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( other.LENGTH )
   , CAPACITY( 0 )
{
   MAC_LABEL( "intVector:: intVector( intVector const& )" ) ;
   
   allocate( other.capacity() ) ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      VECTOR[i] = other.VECTOR[i] ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( *this == other ) ;
}

//----------------------------------------------------------------------
intVector:: ~intVector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: ~intVector" ) ;
   MAC_CHECK_INV( invariant() ) ;
   deallocate() ;
}

//----------------------------------------------------------------------
intVector const&
intVector:: operator=( intVector const& other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: operator=" ) ;
   MAC_CHECK_INV( invariant() ) ;

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
   intVector const& result = *this ;
  
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( result == other ) ;
   return( result ) ;
}

//---------------------------------------------------------------------
void
intVector:: re_initialize( size_t dim, int val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: re_initialize" ) ;
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

//----------------------------------------------------------------------
bool
intVector:: operator==( intVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: operator==" ) ;

   bool result = ( this == &other ) ;
   if( !result )
   {
      result = (LENGTH==other.LENGTH) ;
      for( size_t i=0 ; result && i<size() ; ++i )
      {
         result = operator()(i)==other(i) ;
      }
   }

   MAC_CHECK_POST( IMPLIES( result, size()==other.size() ) ) ;
   MAC_CHECK_POST( !result || FORALL( ( size_t i=0 ; i<size() ; ++i ),
                                       operator()(i) == other(i) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
intVector:: operator!=( intVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: operator!=" ) ;

   bool result = !operator==( other ) ;

   MAC_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
intVector:: capacity( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: capacity" ) ;
   
   size_t result = CAPACITY ;
   
   MAC_CHECK_POST( result >= size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
intVector:: index_of( int val ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: index_of" ) ;

   size_t result = MAC::bad_index() ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      if ( VECTOR[i] == val )
      {
         result = i ;
         break ;
      }
   }

   MAC_CHECK_POST( result == MAC::bad_index() || operator()(result) == val ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
intVector:: has( int val ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: has" ) ;

   size_t ifound = index_of( val ) ;
   bool result = ( ifound != MAC::bad_index() ) ;

   MAC_CHECK_POST( IMPLIES( result, operator()(index_of(val))==val ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
int
intVector:: sum( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: sum" ) ;
   int resu = 0 ;
   for( size_t i = 0 ; i<LENGTH ; ++i )
   {
      resu += operator()( i ) ;
   }
   return resu ;
}

//----------------------------------------------------------------------
void
intVector:: resize( size_t dim, int val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: resize" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   if( dim > CAPACITY ) 
   {
      size_t const newCapacity = MAC_System::new_block_size( CAPACITY, dim ) ;
      int const* oldVector = VECTOR ;
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
intVector:: set( int val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: set" ) ;
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
intVector:: append( int val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: append" ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   int n = LENGTH ;
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
intVector:: extend( int val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: extend" ) ;
   MAC_SAVEOLD( bool, has, has(val) ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   size_t i = index_of( val ) ;
   if( i == MAC::bad_index() )
   {
      append( val ) ;
   }

   MAC_CHECK_POST( has( val ) ) ;
   MAC_CHECK_POST( !OLD(has) || size()==OLD(size) ) ;
   MAC_CHECK_POST( IMPLIES( OLD(has), size()==OLD(size) ) ) ;
   MAC_CHECK_POST( IMPLIES( OLD(has), capacity() == OLD(capacity) ) ) ;
   MAC_CHECK_POST( 
         IMPLIES( !OLD(has), 
                  size()==OLD(size)+1 && operator()(size()-1)==val ) ) ;
   MAC_CHECK_POST( 
      IMPLIES( !OLD(has) && ( OLD(size)+1 <= OLD(capacity) ),
               capacity() == OLD( capacity ) ) ) ;
   MAC_CHECK_POST( 
      IMPLIES( !OLD(has) && ( OLD(size)+1 > OLD(capacity) ),
               capacity() == MAC_System::new_block_size( OLD(capacity),
                                                         OLD(size)+1 ) ) ) ;
}

//----------------------------------------------------------------------
void
intVector:: remove_at( size_t idx )
//----------------------------------------------------------------------
{
   MAC_LABEL( "intVector:: remove_at" ) ;
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
std::ostream&
operator<<( std::ostream& out, intVector const& vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, intVector const& )" ) ;

   size_t j=0 ;
   out << "< " ; 
   for( size_t i = 0 ; i<vec.size() ; ++i )
   {
      if( ++j>5 )
      {
         j = 1 ;
         out << std::endl << "  " ;
      }
      out << vec( i ) << " " ;
   }
   out << ">" ;
   return( out ) ;
}

//----------------------------------------------------------------------
std::istream&
operator>>( std::istream& in, intVector& vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator>>( std::ostream&, intVector const& )" ) ;

   intVector aux( 0 ) ;
   
   std::string c ;
   int val = 0 ;
   
   if( !( in >> c ) || c != "<" )
   {
      intVector_ERROR::n0() ;
   }
   while( in >> val )
   {
      aux.append( val ) ;
   }
   in.clear() ;
   if( !( in >> c ) || c != ">" )
   {
      intVector_ERROR::n0() ;
   }

   vec = aux ;
   
   return( in ) ;
}

//----------------------------------------------------------------------
void
intVector:: allocate( size_t a_size )
//----------------------------------------------------------------------
{
   MAC_CHECK( VECTOR == 0 ) ;

   if( a_size > 0 )
   {
      VECTOR = new int [a_size] ;
      MAC_CHECK( VECTOR != 0 ) ;
   }
   CAPACITY = a_size ;
}

//----------------------------------------------------------------------
void
intVector:: deallocate( void )
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
intVector:: invariant( void ) const
//-----------------------------------------------------------------------
{
   MAC_ASSERT( LENGTH==0 || VECTOR != 0 ) ;
   MAC_ASSERT( LENGTH<=CAPACITY ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
intVector_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   MAC_Error::object()->raise_plain(
      "*** Syntax error reading vector of \"int\":\n"
      "    < val1 val2 ... > is expected." ) ;
}
