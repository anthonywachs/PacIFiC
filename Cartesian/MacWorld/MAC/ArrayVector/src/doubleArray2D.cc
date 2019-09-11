#include <doubleArray2D.hh>

#ifdef OUTLINE
#define inline
#include <doubleArray2D.icc>
#undef inline
#endif

#include <ios>
#include <iostream>
#include <iomanip>

using std::endl ; 
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

//----------------------------------------------------------------------
doubleArray2D:: doubleArray2D( size_t dim0, size_t dim1, double val )
//----------------------------------------------------------------------
   : vector( dim0*dim1, val )
   , d0( dim0 )
   , d1( dim1 )
{
   MAC_LABEL( "doubleArray2D:: doubleArray2D" ) ;
   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
doubleArray2D:: ~doubleArray2D( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
doubleArray2D:: doubleArray2D( doubleArray2D const& other )
//----------------------------------------------------------------------
   : vector( other.vector )
   , d0( other.d0 )
   , d1( other.d1 )
{
   MAC_LABEL( "doubleArray2D:: doubleArray2D( doubleArray2D const& )" ) ;
   MAC_CHECK_POST( *this == other ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
doubleArray2D:: operator=( doubleArray2D const& other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray2D:: operator=" ) ;
   
   if( this != &other )
   {
      vector = other.vector ;
      d0 = other.d0 ;
      d1 = other.d1 ;
   }

   doubleArray2D const& result = *this ;

   MAC_CHECK_POST( result == other ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
doubleArray2D:: re_initialize( size_t dim0, size_t dim1, double val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray2D:: re_initialize" ) ;
   
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
doubleArray2D:: operator==( doubleArray2D const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray2D:: operator==" ) ;

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
doubleArray2D:: operator!=( doubleArray2D const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray2D:: operator!=" ) ;

   bool result = !operator==( other ) ;

   MAC_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
doubleArray2D:: raise_first_index_bound( size_t dim0, double val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray2D:: raise_first_index_bound" ) ;
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
doubleArray2D:: set( double val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray2D:: set" ) ;
   vector.set( val ) ;

   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
void
doubleArray2D:: set_section( size_t an_index, size_t index_value, 
                            doubleVector const& x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray2D:: set_section" ) ;
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
doubleArray2D:: extract_section( size_t an_index, size_t index_value,
                                 doubleVector& x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "doubleArray2D:: extract_section" ) ;
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
operator<<( std::ostream& out, doubleArray2D const& a )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, doubleArray2D const& )" ) ;

   ios_base::fmtflags original_flags = out.flags() ;
   out.setf( ios_base::uppercase | ios_base::scientific ) ;
   std::streamsize p = out.precision() ;
   out << setprecision( 6 )  ;

   out << setw( 5 ) << "i" << setw( 5 ) << "j" 
       << setw( 15 ) << "item(i,j)" << endl ;
   for( size_t i = 0 ; i<a.d0 ; i++ )
   {
      for( size_t j = 0 ; j<a.d1 ; j++ )
      {
         out << setw( 5 ) << i << setw( 5 ) << j ;
         out << setw( 15 ) << a( i, j ) << endl ;
      }
   }

   out << std::setprecision(p) ;
   out.flags( original_flags ) ;
   return( out ) ;
}
