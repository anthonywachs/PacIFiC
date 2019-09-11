#include <doubleArray4D.hh>

#include <iostream>

#ifdef OUTLINE
#define inline
#include <doubleArray4D.icc>
#undef inline
#endif

//----------------------------------------------------------------------
doubleArray4D:: doubleArray4D( size_t dim0, size_t dim1, size_t dim2,
                               size_t dim3, double val )
//----------------------------------------------------------------------
   : vector( dim0*dim1*dim2*dim3, val )
   , d0( dim0 )
   , d1( dim1 )
   , d2( dim2 )
   , d3( dim3 )
{
   MAC_LABEL( "doubleArray4D:: doubleArray4D" ) ;
   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( index_bound(2) == dim2 ) ;
   MAC_CHECK_POST( index_bound(3) == dim3 ) ;
#ifndef _WIN32
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                      operator()(i0,i1,i2,i3) == val ))))) ;
#endif
}

//----------------------------------------------------------------------
doubleArray4D:: ~doubleArray4D( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
doubleArray4D:: re_initialize( size_t dim0, size_t dim1, size_t dim2,
                               size_t dim3, double val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray4D:: re_initialize" ) ;
   
   d0 = dim0 ;
   d1 = dim1 ;
   d2 = dim2 ;
   d3 = dim3 ;
   vector.re_initialize( d0*d1*d2*d3, val ) ;

   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( index_bound(2) == dim2 ) ;
   MAC_CHECK_POST( index_bound(3) == dim3 ) ;
#ifndef _WIN32
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                      operator()(i0,i1,i2,i3) == val ))))) ;
#endif
}

//----------------------------------------------------------------------
void
doubleArray4D:: raise_first_index_bound( size_t dim0, double val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray4D:: raise_first_index_bound" ) ;
   MAC_CHECK_PRE( dim0 > index_bound(0) ) ;
   MAC_CHECK_PRE( index_bound(1) > 0 &&
                  index_bound(2) > 0 &&
                  index_bound(3) > 0 ) ;
   MAC_SAVEOLD( size_t, index_bound, index_bound(0) ) ;

   d0 = dim0 ;
   vector.resize( d0*d1*d2*d3, val ) ;

   MAC_CHECK_POST( index_bound(0)==dim0 ) ;
#ifndef _WIN32
   MAC_CHECK_POST( 
         FORALL( ( size_t i0=OLD(index_bound) ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                      operator()(i0,i1,i2,i3) == val ))))) ;
#endif
}

//----------------------------------------------------------------------
size_t
doubleArray4D:: index_bound( size_t an_index ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray4D:: index_bound" ) ;
   MAC_CHECK_PRE( an_index < 4 ) ;
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
      default : 
         result = ( d3 ) ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
void
doubleArray4D:: set( double val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray4D:: set" ) ;

   vector.set( val ) ;

#ifndef _WIN32
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   FORALL( ( size_t i2=0 ; i2<index_bound(2) ; ++i2 ),
                   FORALL( ( size_t i3=0 ; i3<index_bound(3) ; ++i3 ),
                      operator()(i0,i1,i2,i3) == val ))))) ;
#endif
}

//----------------------------------------------------------------------
std::ostream&
operator<<( std::ostream& out, doubleArray4D const& a )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, doubleArray4D const& )" ) ;
   out << a.d0 << " " << a.d1 << " " << a.d2
       << " " << a.d3 << " : " << std::endl ;
   
   for( size_t i = 0 ; i<a.d0 ; i++ )
   {
      for( size_t j = 0 ; j<a.d1 ; j++ )
      {
         for( size_t k = 0 ; k<a.d2 ; k++ )
         {
            for( size_t l = 0 ; l<a.d3 ; l++ )
            {
               out << i << " " << j << " " << k
                   << " " << l << " " << a( i, j, k, l ) << std::endl ;
            }
         }
      }
   }
   return( out ) ;
}
