#include <size_t_vector.hh>

#ifdef OUTLINE
#define inline
#include <size_t_vector.icc>
#undef inline
#endif

#include <MAC.hh>
#include <MAC_Error.hh>
#include <MAC_System.hh>
#include <intVector.hh>

#include <iostream>

struct size_t_vector_ERROR
{
   static void n0( void ) ;
} ;

//----------------------------------------------------------------------
size_t_vector:: size_t_vector( size_t dim, size_t val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( dim )
   , CAPACITY( 0 )
{
   MAC_LABEL( "size_t_vector:: size_t_vector( size_t )" ) ;

   allocate( LENGTH ) ;
   set( val ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == dim ) ;
   MAC_CHECK_POST( capacity() == dim ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == val ) ) ;
}

//----------------------------------------------------------------------
size_t_vector:: size_t_vector( int dim, size_t val )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( static_cast<size_t>( dim ) )
   , CAPACITY( 0 )
{
   MAC_LABEL( "size_t_vector:: size_t_vector( int )" ) ;
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
size_t_vector:: size_t_vector( intVector const& ivec )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( ivec.size() )
   , CAPACITY( 0 )
{
   MAC_LABEL( "size_t_vector:: size_t_vector( intVector )" ) ;
   MAC_CHECK_PRE( FORALL( (size_t i=0 ;i<ivec.size();++i ),
                          ivec(i)>=0 ) ) ;

   allocate( ivec.capacity() ) ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      VECTOR[i] = (size_t)ivec(i) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( size() == ivec.size() ) ;
   MAC_CHECK_POST( capacity() == ivec.capacity() ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size() ; ++i ),
			   operator()(i) == (size_t)ivec(i) ) ) ;
}

//----------------------------------------------------------------------
size_t_vector:: size_t_vector( size_t_vector const& other )
//----------------------------------------------------------------------
   : VECTOR( 0 )
   , LENGTH( other.LENGTH )
   , CAPACITY( 0 )
{
   MAC_LABEL( "size_t_vector:: size_t_vector( size_t_vector const& )" ) ;

   allocate( other.capacity() ) ;
   for( size_t i=0 ; i<LENGTH ; ++i )
   {
      VECTOR[i] = other.VECTOR[i] ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( *this == other ) ;
}

//----------------------------------------------------------------------
size_t_vector:: ~size_t_vector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: ~size_t_vector" ) ;
   MAC_CHECK_INV( invariant() ) ;
   deallocate() ;
}

//----------------------------------------------------------------------
size_t_vector const&
size_t_vector:: operator=( size_t_vector const& other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: operator=" ) ;
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
   size_t_vector const& result = *this ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( result == other ) ;
   return( result ) ;
}

//---------------------------------------------------------------------
void
size_t_vector:: re_initialize( size_t dim, size_t val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: re_initialize" ) ;
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
size_t_vector:: operator==( size_t_vector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: operator==" ) ;

   bool result = ( this == &other ) ;
   if( !result )
   {
      result = (LENGTH == other.LENGTH) ;
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
size_t_vector:: operator!=( size_t_vector const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: operator!=" ) ;

   bool result = !operator==( other ) ;

   MAC_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
size_t_vector:: capacity( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: capacity" ) ;
   
   size_t result = CAPACITY ;
   
   MAC_CHECK_POST( result >= size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
size_t_vector:: index_of( size_t val ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: index_of" ) ;
   
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
size_t_vector:: has( size_t val ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: has" ) ;
   
   size_t ifound = index_of( val ) ;
   bool result = ( ifound != MAC::bad_index() ) ;

   MAC_CHECK_POST( IMPLIES( result, operator()(index_of(val))==val ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
size_t_vector:: sum( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: sum" ) ;
   size_t resu = 0 ;
   for( size_t i = 0 ; i<LENGTH ; ++i )
   {
      resu += operator()( i ) ;
   }
   return resu ;
}

//----------------------------------------------------------------------
void
size_t_vector:: resize( size_t dim, size_t val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: resize" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;
    
   if( dim > CAPACITY ) 
   {
      size_t const newCapacity = MAC_System::new_block_size( CAPACITY, dim ) ;
      size_t const* oldVector = VECTOR ;
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
size_t_vector:: set( size_t val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: set" ) ;
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
size_t_vector:: append( size_t val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: append" ) ;
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
size_t_vector:: extend( size_t val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: extend" ) ;
   MAC_SAVEOLD( bool, has, has(val) ) ;
   MAC_SAVEOLD( size_t, size, size() ) ;
   MAC_SAVEOLD( size_t, capacity, capacity() ) ;

   size_t i = index_of( val ) ;
   if( i == MAC::bad_index() )
   {
      append( val ) ;
   }

   MAC_CHECK_POST( has( val ) ) ;
   MAC_CHECK_POST( IMPLIES( OLD(has), size() == OLD(size) ) ) ;
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
size_t_vector:: remove_at( size_t idx )
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: remove_at" ) ;
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
size_t_vector:: sort_increasingly( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "size_t_vector:: sort_increasingly" ) ;

   // BUBBLE SORT
   bool modif = true ;
   while(modif)
   {
      modif = false ;
      for( size_t i=0 ; i<LENGTH-1 ; ++i )
      {
         if(VECTOR[i]>VECTOR[i+1])
         {
            size_t tmp = VECTOR[i] ;
            VECTOR[i] = VECTOR[i+1] ;
            VECTOR[i+1] = tmp ;
            modif = true ;
         }
      }
   }

   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<size()-1 ; ++i ),
            operator()(i) <= operator()(i+1) ) ) ;
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, size_t_vector const& vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, size_t_vector const& )" ) ;
   
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
operator>>( std::istream& in, size_t_vector& vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator>>( std::istream&, size_t_vector& )" ) ;

   size_t_vector aux( 0 ) ;
   
   std::string c ;
   int val = 0 ;
   
   if( !( in >> c ) || c != "<" )
   {
      size_t_vector_ERROR::n0() ;
   }
   while( in >> val )
   {
      if( val<0 ) size_t_vector_ERROR::n0() ;
      aux.append( (size_t) val ) ;
   }
   in.clear() ;
   if( !( in >> c ) || c != ">" )
   {
      size_t_vector_ERROR::n0() ;
   }
   
   vec = aux ;
   
   return( in ) ;
}

//----------------------------------------------------------------------
void
size_t_vector:: allocate( size_t a_size )
//----------------------------------------------------------------------
{
   MAC_CHECK( VECTOR == 0 ) ;

   if( a_size > 0 )
   {
      VECTOR = new size_t [a_size] ;
      MAC_CHECK( VECTOR != 0 ) ;
   }
   CAPACITY = a_size ;
}

//----------------------------------------------------------------------
void
size_t_vector:: deallocate( void )
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
size_t_vector:: invariant( void ) const
//-----------------------------------------------------------------------
{
   MAC_ASSERT( LENGTH==0 || VECTOR != 0 ) ;
   MAC_ASSERT( LENGTH<=CAPACITY ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------
void 
size_t_vector_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   MAC_Error::object()->raise_plain(
      "*** Syntax error reading vector of \"size_t\":\n"
      "    < val1 val2 ... > is expected." ) ;
}


