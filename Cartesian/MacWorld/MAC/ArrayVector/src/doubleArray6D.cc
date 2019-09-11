#include <doubleArray6D.hh>

#include <iostream>

#ifdef OUTLINE
#define inline
#include <doubleArray6D.icc>
#undef inline
#endif

//----------------------------------------------------------------------
doubleArray6D:: doubleArray6D( size_t dim0, size_t dim1, size_t dim2,
                               size_t dim3, size_t dim4, size_t dim5,
                               double val )
//----------------------------------------------------------------------
   : vector( dim0*dim1*dim2*dim3*dim4*dim5, val )
   , d0( dim0 )
   , d1( dim1 )
   , d2( dim2 )
   , d3( dim3 )
   , d4( dim4 )
   , d5( dim5 )
{
   MAC_LABEL( "doubleArray6D:: doubleArray6D" ) ;
   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( index_bound(2) == dim2 ) ;
   MAC_CHECK_POST( index_bound(3) == dim3 ) ;
   MAC_CHECK_POST( index_bound(4) == dim4 ) ;
   MAC_CHECK_POST( index_bound(5) == dim5 ) ;
#ifndef _WIN32
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                   FORALL( ( size_t i4=0 ; i4<index_bound(4) ; ++i4 ),
                   FORALL( ( size_t i5=0 ; i5<index_bound(5) ; ++i5 ),
                      operator()(i0,i1,i2,i3,i4,i5) == val ))))))) ;
#endif
}

//----------------------------------------------------------------------
doubleArray6D:: ~doubleArray6D( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
doubleArray6D:: re_initialize( size_t dim0, size_t dim1, size_t dim2,
                               size_t dim3, size_t dim4, size_t dim5,
                               double val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray6D:: re_initialize" ) ;
   
   d0 = dim0 ;
   d1 = dim1 ;
   d2 = dim2 ;
   d3 = dim3 ;
   d4 = dim4 ;
   d5 = dim5 ;
   vector.re_initialize( d0*d1*d2*d3*d4*d5, val ) ;

   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( index_bound(2) == dim2 ) ;
   MAC_CHECK_POST( index_bound(3) == dim3 ) ;
   MAC_CHECK_POST( index_bound(4) == dim4 ) ;
   MAC_CHECK_POST( index_bound(5) == dim5 ) ;
#ifndef _WIN32
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                   FORALL( ( size_t i4=0 ; i4<index_bound(4) ; ++i4 ),
                   FORALL( ( size_t i5=0 ; i5<index_bound(5) ; ++i5 ),
                      operator()(i0,i1,i2,i3,i4,i5) == val ))))))) ;
#endif
}

//----------------------------------------------------------------------
void
doubleArray6D:: raise_first_index_bound( size_t dim0, double val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray6D:: raise_first_index_bound" ) ;
   MAC_CHECK_PRE( dim0 > index_bound(0) ) ;
   MAC_CHECK_PRE( index_bound(1) > 0 &&
                  index_bound(2) > 0 &&
                  index_bound(3) > 0 &&
                  index_bound(4) > 0 &&
                  index_bound(5) > 0 ) ;
   MAC_SAVEOLD( size_t, index_bound, index_bound(0) ) ;

   d0 = dim0 ;
   vector.resize( d0*d1*d2*d3*d4*d5, val ) ;

   MAC_CHECK_POST( index_bound(0)==dim0 ) ;
#ifndef _WIN32
   MAC_CHECK_POST( 
         FORALL( ( size_t i0=OLD(index_bound) ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                   FORALL( ( size_t i4=0 ; i4<index_bound(4) ; ++i4 ),
                   FORALL( ( size_t i5=0 ; i5<index_bound(5) ; ++i5 ),
                      operator()(i0,i1,i2,i3,i4,i5) == val ))))))) ;
#endif
}

//----------------------------------------------------------------------
size_t
doubleArray6D:: index_bound( size_t an_index ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray6D:: index_bound" ) ;
   MAC_CHECK_PRE( an_index < 6 ) ;
   size_t result ;
   
   switch( an_index )
   {
      case 0 : 
         result = ( d0 ) ;
         break ;
      case 1 : 
         result = ( d1 ) ;
         break ;
      case 2 : 
         result = ( d2 ) ;
         break ;
      case 3 : 
         result = ( d3 ) ;
         break ;
      case 4 : 
         result = ( d4 ) ;
         break ;
      default : 
         result = ( d5 ) ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
void
doubleArray6D:: set( double val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray6D:: set" ) ;

   vector.set( val ) ;

#ifndef _WIN32
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                   FORALL( ( size_t i4=0 ; i4<index_bound(4) ; ++i4 ),
                   FORALL( ( size_t i5=0 ; i5<index_bound(5) ; ++i5 ),
                      operator()(i0,i1,i2,i3,i4,i5) == val ))))))) ;
#endif
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, doubleArray6D const& a )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, doubleArray6D const& )" ) ;
   out << a.d0 << " " << a.d1 << " " << a.d2
       << " " << a.d3 << " " << a.d4 << " " << a.d5 << " : " << std::endl ;
   
   for( size_t i = 0 ; i<a.d0 ; i++ )
   {
      for( size_t j = 0 ; j<a.d1 ; j++ )
      {
         for( size_t k = 0 ; k<a.d2 ; k++ )
         {
            for( size_t l = 0 ; l<a.d3 ; l++ )
            {
               for( size_t m = 0 ; m<a.d4 ; m++ )
               {
                  for( size_t n = 0 ; n<a.d5 ; n++ )
                  {
                     out << i << " " << j << " " << k
                         << " " << l << " " << m
                         << " " << n << " " << a( i, j, k, l, m, n )
                         << std::endl ;
                  }
               }
            }
         }
      }
   }
   return( out ) ;
}
