#include <boolArray2D.hh>

#include <iostream>

#ifdef OUTLINE
#define inline
#include <boolArray2D.icc>
#undef inline
#endif

//----------------------------------------------------------------------
boolArray2D:: boolArray2D( size_t dim0, size_t dim1, bool val )
//----------------------------------------------------------------------
   : vector( dim0*dim1, val )
   , d0( dim0 )
   , d1( dim1 )
{
   MAC_LABEL( "boolArray2D:: boolArray2D( size_t, size_t )" ) ;
   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
boolArray2D:: ~boolArray2D( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
boolArray2D:: boolArray2D( boolArray2D const& other )
//----------------------------------------------------------------------
   : vector( other.vector )
   , d0( other.d0 )
   , d1( other.d1 )
{
   MAC_LABEL( "boolArray2D:: boolArray2D( boolArray2D const& )" ) ;
   MAC_CHECK_POST( *this == other ) ;
}

//----------------------------------------------------------------------
boolArray2D const&
boolArray2D:: operator=( boolArray2D const& other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray2D:: operator=" ) ;

   if( this != &other )
   {
      vector = other.vector ;
      d0 = other.d0 ;
      d1 = other.d1 ;
   }
   boolArray2D const& result = *this ;

   MAC_CHECK_POST( result == other ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
boolArray2D:: re_initialize( size_t dim0, size_t dim1, bool val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray2D:: re_initialize" ) ;
   
   d0 = dim0 ;
   d1 = dim1 ;

   vector.re_initialize( d0*d1, val ) ;

   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
bool
boolArray2D:: operator==( boolArray2D const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray2D:: operator==" ) ;

   bool result = ( d0==other.d0 &&
                   d1==other.d1 &&
                   vector==other.vector ) ;

   MAC_CHECK_POST(
      IMPLIES( result,
               index_bound(0)==other.index_bound(0) &&
               index_bound(1)==other.index_bound(1) ) ) ;
   MAC_CHECK_POST(
      !result ||
      FORALL( ( size_t i=0 ; i<index_bound(0) ; ++i ),
              FORALL( ( size_t j=0 ; j<index_bound(1) ; ++j ),
                      operator()(i,j) == other(i,j) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
boolArray2D:: operator!=( boolArray2D const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray2D:: operator!=" ) ;
  
   bool result = !operator==( other ) ;
   
   MAC_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;   
}

//----------------------------------------------------------------------
void
boolArray2D:: raise_first_index_bound( size_t dim0, bool val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray2D:: raise_first_index_bound" ) ;
   MAC_CHECK_PRE( dim0 > index_bound(0) ) ;
   MAC_CHECK_PRE( index_bound(1) > 0 ) ;
   MAC_SAVEOLD( size_t, index_bound, index_bound(0) ) ;

   d0 = dim0 ;
   vector.resize( d0*d1, val ) ;

   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( 
      FORALL( ( size_t i0=OLD(index_bound) ; i0<index_bound(0) ; ++i0 ),
         FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
            operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
void
boolArray2D:: set( bool val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray2D:: set" ) ;

   vector.set( val ) ;

   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
void
boolArray2D:: set_section( size_t an_index, size_t index_value, 
                          boolVector const& x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray2D:: set_section" ) ;
   MAC_CHECK_PRE( an_index < 2 ) ;
   MAC_CHECK_PRE( index_value < index_bound(an_index) ) ;
   MAC_CHECK_PRE( IMPLIES( an_index==0, x.size() == index_bound(1) ) ) ;
   MAC_CHECK_PRE( IMPLIES( an_index==1, x.size() == index_bound(0) ) ) ;
   
   if( an_index == 0 )
   {
      for( size_t k=0 ; k<d1 ; k++)
      {
         operator()(index_value,k) = x(k) ;
      }
   }
   else if( an_index == 1 )
   {
      for( size_t k=0 ; k<d0 ; k++ )
      {
         operator()(k,index_value) = x(k) ;
      }
   }

   MAC_CHECK_POST( an_index==1 ||
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                      operator()(index_value,i1) == x(i1) ) ) ;
   MAC_CHECK_POST( an_index==0 ||
                   FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                      operator()(i0,index_value) == x(i0) ) ) ;
}

//----------------------------------------------------------------------
void
boolArray2D:: extract_section( size_t an_index, size_t index_value,
                               boolVector& x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray2D:: extract_section" ) ;
   MAC_CHECK_PRE( an_index < 2 ) ;
   MAC_CHECK_PRE( index_value < index_bound(an_index) ) ;
   MAC_CHECK_PRE( IMPLIES( an_index==0, x.size() == index_bound(1) ) ) ;
   MAC_CHECK_PRE( IMPLIES( an_index==1, x.size() == index_bound(0) ) ) ;
   
   if( an_index == 0 )
   {
      for( size_t k=0 ; k<d1 ; k++)
      {
         x(k) = operator()(index_value,k) ;
      }
   }
   else if( an_index == 1 )
   {
      for( size_t k=0 ; k<d0 ; k++ )
      {
         x(k) = operator()(k,index_value) ;
      }
   }

   MAC_CHECK_POST( an_index==1 || 
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                           x(i1) == operator()(index_value,i1) ) ) ;
   MAC_CHECK_POST( an_index==0 || 
                   FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                           x(i0) == operator()(i0,index_value) ) ) ;
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, boolArray2D const& a )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, boolArray2D const& )" ) ;
   for( size_t i = 0 ; i<a.d0 ; i++ )
   {
      for( size_t j = 0 ; j<a.d1 ; j++ )
      {
         out << a( i, j ) << " " ;
      }
      out << std::endl ;
   }
   return( out ) ;
}
