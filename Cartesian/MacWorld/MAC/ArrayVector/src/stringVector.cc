#include <stringVector.hh>

#ifdef OUTLINE
#define inline
#include <stringVector.icc>
#undef inline
#endif

#include <MAC.hh>
#include <MAC_Error.hh>
#include <MAC_System.hh>
#include <iostream>

struct stringVector_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
stringVector:: stringVector( size_t dim )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   MAC_LABEL( "stringVector:: stringVector( size_t )" ) ;

   allocate( LENGTH ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == dim ) ;
   MAC_CHECK_POST( capacity() == dim ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; i++ ),
			   operator()(i) == "" ) ) ;
}

//----------------------------------------------------------------------
stringVector:: stringVector( int dim )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   MAC_LABEL( "stringVector:: stringVector( int )" ) ;
   MAC_CHECK_PRE( dim >= 0 ) ;

   allocate( LENGTH ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == static_cast<size_t>(dim) ) ;
   MAC_CHECK_POST( capacity() == size() ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; i++ ),
			   operator()(i) == "" ) ) ;
}

//----------------------------------------------------------------------
stringVector:: stringVector( std::string const& elements,
                             char const separator )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( 0 )
   , CAPACITY( 0 )
{
   MAC_LABEL( "stringVector:: stringVector( string, char )" ) ;

   build_from_string( elements, separator ) ;
}

//----------------------------------------------------------------------
stringVector:: stringVector( char const* elements,
                             char const separator )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( 0 )
   , CAPACITY( 0 )
{
   MAC_LABEL( "stringVector:: stringVector( const char*, char )" ) ;

   build_from_string( elements, separator ) ;
}

//----------------------------------------------------------------------
stringVector:: stringVector( stringVector const& other )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( other.LENGTH )
   , CAPACITY( 0 )
{
   MAC_LABEL( "stringVector:: stringVector( stringVector const& )" ) ;

   allocate( other.capacity() ) ;
   if( LENGTH != 0 )
   {
      for( size_t i=0 ; i<LENGTH ; ++i )
      {
         VECTOR[i] = other(i) ;
      }
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == other.size() ) ;
   MAC_CHECK_POST( capacity() == other.capacity() ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
                           operator()(i) == other(i) ) ) ;
}

//----------------------------------------------------------------------
stringVector:: ~stringVector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: ~stringVector" ) ;
   MAC_CHECK_INV( invariant() ) ;
   deallocate() ;
}

//----------------------------------------------------------------------
stringVector const&
stringVector:: operator=( stringVector const& other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: operator=" ) ;
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
      if( LENGTH != 0 )
      {
         for( size_t i=0 ; i<LENGTH ; ++i )
         {
            VECTOR[i] = other(i) ;
         }
      }
   }
   stringVector const& result = *this ;

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
stringVector:: re_initialize( size_t dim )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: re_initialize" ) ;
   MAC_CHECK_INV( invariant() ) ;

   LENGTH = dim ;
   if( CAPACITY != dim )
   {
      deallocate() ;
      allocate( LENGTH ) ;
   }
   else if( LENGTH != 0 )
   {
      for( size_t i=0 ; i<LENGTH ; ++i )
      {
         VECTOR[i] = "" ;
      }
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == dim ) ;
   MAC_CHECK_POST( capacity() == dim ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == "" ) ) ;
}

//----------------------------------------------------------------------
bool
stringVector:: operator==( stringVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: operator==" ) ;

   bool result = ( this == &other ) ;
   if( !result )
   {
      result = (LENGTH == other.LENGTH) ;
      for( size_t i=0 ; result && i<size() ; i++ )
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
stringVector:: operator!=( stringVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: operator!=" ) ;

   bool result = !operator==( other ) ;

   MAC_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
int
stringVector:: three_way_comparison( stringVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: three_way_comparison" ) ;

   size_t min_size = MAC::min( size(), other.size() ) ;
   int result = 0 ;

   size_t i=0 ;
   for( ; i < min_size ; ++i )
   {
      int comp = operator()(i).compare( other.operator()(i) ) ;
      if( comp != 0 )
      {
         result = comp ;
         break ;
      }
   }

   if( result == 0 )
   {
      if( i == size() && i< other.size() )
      {
         result = -1 ;
      }
      else if( i < size() && i == other.size() )
      {
         result = 1 ;
      }
   }

  return( result ) ;
}

//----------------------------------------------------------------------
bool
stringVector:: operator>( stringVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: operator>" ) ;

  return( three_way_comparison( other ) > 0 ) ;
}

//----------------------------------------------------------------------
bool
stringVector:: operator<( stringVector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: operator<" ) ;

   return( three_way_comparison( other ) < 0 ) ;
}

//----------------------------------------------------------------------
void stringVector:: build_from_string( std::string const& elements,
                                       char const separator )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: build_from_string" ) ;

   if( !elements.empty() )
   {
      size_t start = 0 ;
      size_t idx = 0 ;
      size_t end = elements.length()  ;
      
      while( ( idx = elements.find( separator, start ) ) < end )
      {
         std::string a_choice = elements.substr( start, idx-start ) ;
         MAC::remove_enclosing_characters(a_choice,'\"') ;
         append( a_choice ) ;
         start = idx+1 ;
      }
      std::string a_choice = elements.substr( start, elements.length()-start ) ;
      MAC::remove_enclosing_characters(a_choice,'\"') ;
      append( a_choice ) ;
   }
}

//----------------------------------------------------------------------
size_t
stringVector:: capacity( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: capacity" ) ;

   size_t result = CAPACITY ;

   MAC_CHECK_POST( result >= size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
stringVector:: index_of( std::string const& val ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: index_of" ) ;

   size_t result = MAC::bad_index() ;
   for( size_t i=0 ; i<LENGTH ; i++ )
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
stringVector:: has( std::string const& val ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: has" ) ;

   size_t ifound = index_of( val ) ;
   bool result = ( ifound != MAC::bad_index() ) ;

   MAC_CHECK_POST( IMPLIES( result, operator()(index_of(val))==val ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
stringVector:: resize( size_t dim, std::string const& val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: resize" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   if( dim > CAPACITY )
   {
      size_t const newCapacity = MAC_System::new_block_size( CAPACITY, dim ) ;
      std::string const* oldVector = VECTOR ;
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
      for( size_t i=LENGTH ; i<dim ; i++ )
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
stringVector:: set( std::string const& val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: set" ) ;
   MAC_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<LENGTH ; i++ )
   {
      VECTOR[i] = val ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
void
stringVector:: append( std::string const& val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: append" ) ;
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
stringVector:: extend( std::string const& val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: extend" ) ;
   MAC_SAVEOLD( bool, has, has(val) ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   size_t i = index_of( val ) ;
   if( i >= LENGTH )
   {
      append( val ) ;
   }

   MAC_CHECK_POST( has( val ) ) ;
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
stringVector:: remove_at( size_t idx )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: remove_at" ) ;
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
void
stringVector:: sort( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: sort" ) ;

   // Tri "à bulles"
   for( size_t i=0 ; i<size() ; ++i )
   {
      for( size_t j=i+1 ; j<size() ; ++j )
      {
         if( operator()(i).compare( operator()(j) )>0 )
         {
            std::string dummy = operator()(i) ;
            operator()(i) = operator()(j) ;
            operator()(j) = dummy ;
         }
      }
   }
}

//----------------------------------------------------------------------
stringVector stringVector:: sorted_copy( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringVector:: sorted_copy" ) ;

   stringVector result( *this ) ;
   result.sort() ;

   return( result ) ;
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, stringVector const& vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, stringVector const& )" ) ;

   size_t j=0 ;
   out << "< " ;
   for( size_t i = 0 ; i<vec.size() ; i++ )
   {
      if( ++j>5 )
      {
         j = 1 ;
         out << std::endl << "  " ;
      }
      out << "\"" << vec( i ) << "\" " ;
   }
   out << ">" ;
   return( out ) ;
}

//----------------------------------------------------------------------
std::istream&
operator>>( std::istream& in, stringVector& vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator>>( std::ostream&, stringVector const& )" ) ;

   stringVector aux( 0 ) ;

   std::string c ;

   if( !( in >> c ) || c != "<" )
   {
      stringVector_ERROR::n0() ;
   }
   while( in >> c && c != ">" )
   {
      if( c[0] != '\"' || c[c.size()-1] != '\"' )
      {
         stringVector_ERROR::n0() ;
      }
      std::string c2( c, 1, c.size()-2 ) ;
      aux.append( c2 ) ;
   }
   if( c != ">" )
   {
      stringVector_ERROR::n0() ;
   }
   in.clear() ;

   vec = aux ;

   return( in ) ;
}

//----------------------------------------------------------------------
void
stringVector:: allocate( size_t a_size )
//----------------------------------------------------------------------
{
   MAC_CHECK( VECTOR == 0 ) ;

   if( a_size > 0 )
   {
      VECTOR = new std::string [a_size] ;
      MAC_CHECK( VECTOR != 0 ) ;
      for( size_t i=0 ; i<a_size ; ++i )
      {
         VECTOR[i] = "" ;
      }
   }
   CAPACITY = a_size ;
}

//----------------------------------------------------------------------
void
stringVector:: deallocate( void )
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
stringVector:: invariant( void ) const
//-----------------------------------------------------------------------
{
   MAC_ASSERT( LENGTH==0 || VECTOR != 0 ) ;
   MAC_ASSERT( LENGTH<=CAPACITY ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void
stringVector_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   MAC_Error::object()->raise_plain(
      "*** Syntax error reading vector of \"string\":\n"
      "    < val1 val2 ... > is expected." ) ;
}

