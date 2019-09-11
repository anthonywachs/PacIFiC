#include <boolArray3D.hh>

#include <iostream>

#ifdef OUTLINE
#define inline
#include <boolArray3D.icc>
#undef inline
#endif

//----------------------------------------------------------------------
boolArray3D:: boolArray3D( size_t dim0, size_t dim1, size_t dim2,
                           bool val )
//----------------------------------------------------------------------
   : vector( dim0*dim1*dim2, val )
   , d0( dim0 )
   , d1( dim1 )
   , d2( dim2 )
{
   MAC_LABEL( "boolArray3D:: boolArray3D" ) ;

   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( index_bound(2) == dim2 ) ;
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                      operator()(i0,i1,i2) == val )))) ;
}

//----------------------------------------------------------------------
boolArray3D:: boolArray3D( boolArray3D const& other )
//----------------------------------------------------------------------
   : vector( other.d0*other.d1*other.d2 )
   , d0( other.d0 )
   , d1( other.d1 )
   , d2( other.d2 )
{
   MAC_LABEL( "boolArray3D:: boolArray3D" ) ;
   for( size_t i0=0 ; i0<index_bound(0) ; ++i0 )
      for( size_t i1=0 ; i1<index_bound(1) ; ++i1 )
         for( size_t i2=0 ; i2<index_bound(2) ; ++i2 )
            operator()(i0,i1,i2) = other(i0,i1,i2) ;
   

   MAC_CHECK_POST( *this == other ) ;
}

//----------------------------------------------------------------------
boolArray3D const&
boolArray3D:: operator=( boolArray3D const& other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray3D:: operator=" ) ;

   if( this != &other )
   {
      vector = other.vector ;
      d0 = other.d0 ;
      d1 = other.d1 ;
      d2 = other.d2 ;
   }
   boolArray3D const& result = *this ;

   MAC_CHECK_POST( result == other ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
boolArray3D:: ~boolArray3D( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
boolArray3D:: re_initialize( size_t dim0, size_t dim1, size_t dim2,
                             bool val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray3D:: re_initialize" ) ;

   d0 = dim0 ;
   d1 = dim1 ;
   d2 = dim2 ;
   vector.re_initialize( d0*d1*d2, val ) ;
   
   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( index_bound(2) == dim2 ) ;
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                      operator()(i0,i1,i2) == val )))) ;
}

//----------------------------------------------------------------------
bool
boolArray3D:: operator==( boolArray3D const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray3D:: operator==" ) ;

   bool result = ( d0==other.d0 &&
                   d1==other.d1 &&
                   d2==other.d2 &&
                   vector==other.vector ) ;

   MAC_CHECK_POST(
      IMPLIES( result,
               index_bound(0)==other.index_bound(0) &&
               index_bound(1)==other.index_bound(1) &&
               index_bound(2)==other.index_bound(2) ) ) ;
   MAC_CHECK_POST(
      !result ||
      FORALL( ( size_t i=0 ; i<index_bound(0) ; ++i ),
              FORALL( ( size_t j=0 ; j<index_bound(1) ; ++j ),
                      FORALL( ( size_t k=0 ; k<index_bound(2) ; ++k ),
                              operator()(i,j,k) == other(i,j,k) ) ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
boolArray3D:: operator!=( boolArray3D const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray3D:: operator!=" ) ;
  
   bool result = !operator==( other ) ;
   
   MAC_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;   
}

//----------------------------------------------------------------------
void
boolArray3D:: raise_first_index_bound( size_t dim0, bool val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray3D:: raise_first_index_bound" ) ;
   MAC_CHECK_PRE( dim0 > index_bound(0) ) ;
   MAC_CHECK_PRE( index_bound(1) > 0 &&
                  index_bound(2) > 0 ) ;
   MAC_SAVEOLD( size_t, index_bound, index_bound(0) ) ;

   d0 = dim0 ;

   vector.resize( d0*d1*d2, val ) ;

   MAC_CHECK_POST( index_bound(0)==dim0 ) ;
   MAC_CHECK_POST( 
         FORALL( ( size_t i0=OLD(index_bound) ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                      operator()(i0,i1,i2) == val )))) ;
}

//----------------------------------------------------------------------
void
boolArray3D:: set( bool val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "boolArray3D:: set" ) ;

   vector.set( val ) ;

   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                      operator()(i0,i1,i2) == val )))) ;
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, boolArray3D const& a )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, boolArray3D const& )" ) ;
   out << a.d0 << " " << a.d1 << " " << a.d2 << " : " << std::endl ;
   
   for( size_t i = 0 ; i<a.d0 ; i++ )
   {
      for( size_t j = 0 ; j<a.d1 ; j++ )
      {
         for( size_t k = 0 ; k<a.d2 ; k++ )
         {
            out << i << " " << j << " " << k
                << " " << a( i, j, k ) << std::endl ;
         }
      }
   }
   return( out ) ;
}
