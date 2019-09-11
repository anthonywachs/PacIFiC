#include <boolVector.hh>

#ifdef OUTLINE
#define inline
#include <boolVector.icc>
#undef inline
#endif

#include <MAC.hh>
#include <MAC_Error.hh>
#include <MAC_System.hh>

#include <iostream>

struct boolVector_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
boolVector:: boolVector( size_t dim, bool val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   MAC_LABEL( "boolVector:: boolVector( size_t )" ) ;
   
   allocate( LENGTH ) ;
   set( val ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == dim ) ;
   MAC_CHECK_POST( capacity() == dim ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
boolVector:: boolVector( int dim, bool val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   MAC_LABEL( "boolVector:: boolVector( int )" ) ;
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
boolVector:: boolVector( boolVector const& other )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( other.LENGTH )
   , CAPACITY( 0 )
{
   MAC_LABEL( "boolVector:: boolVector( boolVector const& )" ) ;
   
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
boolVector:: ~boolVector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: ~boolVector" ) ;
   MAC_CHECK_INV( invariant() ) ;   
   deallocate() ;
}

//----------------------------------------------------------------------
boolVector const&
boolVector:: operator=( boolVector const& other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: operator=" ) ;
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
   boolVector const& result = *this ;
  
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
boolVector:: re_initialize( size_t dim, bool val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: re_initialize" ) ;
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
boolVector:: operator==( boolVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: operator==" ) ;

   bool result = ( this == &other ) ;
   if( !result )
   {
      result = ( LENGTH == other.LENGTH ) ;
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
boolVector:: operator!=( boolVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: operator!=" ) ;

   bool result = !operator==( other ) ;

   MAC_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
boolVector:: capacity( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: capacity" ) ;
   
   size_t result = CAPACITY ;
   
   MAC_CHECK_POST( result >= size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
boolVector:: has( bool val ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: has" ) ;
   bool result = false ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      if ( VECTOR[i] == val )
      {
         result = true ;
         break ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
void
boolVector:: resize( size_t dim, bool val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: resize" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   if( dim > CAPACITY ) 
   {
      size_t const newCapacity = MAC_System::new_block_size( CAPACITY, dim ) ;
      bool const* oldVector = VECTOR ;
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
boolVector:: set( bool val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: set" ) ;
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
boolVector:: append( bool val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: append" ) ;
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
boolVector:: remove_at( size_t idx )
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolVector:: remove_at" ) ;
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
operator<<( std::ostream& out, boolVector const& vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, boolVector const& )" ) ;
   
   size_t j=0 ;
   out << "< " ;
   for( size_t i = 0 ; i<vec.size() ; ++i )
   {
      if( ++j>5 )
      {
         j = 1 ;
         out << std::endl << "  " ;
      }
      out << ( vec( i ) ? "true" : "false" ) << " " ;
   }
   out << ">" ;
   return( out ) ;
}

//----------------------------------------------------------------------
std::istream&
operator>>( std::istream& in, boolVector& vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator>>( std::ostream&, boolVector const& )" ) ;

   boolVector aux( 0 ) ;
   
   std::string c ;
   
   if( !( in >> c ) || c != "<" )
   {
      boolVector_ERROR::n0() ;
   }
   while( in >> c && c != ">" )
   {
      if( c == "true" )
      {
         aux.append( true ) ;
      }
      else if( c == "false" )
      {
         aux.append( false ) ;
      }
      else
      {
         boolVector_ERROR::n0() ;
      }
   }
   if( c != ">" )
   {
      boolVector_ERROR::n0() ;
   }
   in.clear() ;

   vec = aux ;
   
   return( in ) ;
}

//----------------------------------------------------------------------
void
boolVector:: allocate( size_t a_size )
//----------------------------------------------------------------------
{
   MAC_CHECK( VECTOR == 0 ) ;

   if( a_size > 0 )
   {
      VECTOR = new bool [a_size] ;
      MAC_CHECK( VECTOR != 0 ) ;
   }
   CAPACITY = a_size ;
}


//----------------------------------------------------------------------
void
boolVector:: deallocate( void )
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
boolVector:: invariant( void ) const
//-----------------------------------------------------------------------
{
   MAC_ASSERT( LENGTH==0 || VECTOR != 0 ) ;
   MAC_ASSERT( LENGTH<=CAPACITY ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
boolVector_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   MAC_Error::object()->raise_plain(
      "*** Syntax error reading vector of \"bool\":\n"
      "    < val1 val2 ... > is expected." ) ;
}
