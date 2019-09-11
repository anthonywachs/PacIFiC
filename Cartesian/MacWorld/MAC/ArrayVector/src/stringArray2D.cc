#include <stringArray2D.hh>

#ifdef OUTLINE
#define stringArray2D_HH
#include <stringArray2D.icc>
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
stringArray2D:: stringArray2D( size_t dim0, size_t dim1, string val )
//----------------------------------------------------------------------
   : vector( dim0*dim1 )
   , d0( dim0 )
   , d1( dim1 )
{
   MAC_LABEL( "stringArray2D:: stringArray2D" ) ;
   vector.set( val ) ;
   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
stringArray2D:: ~stringArray2D( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
stringArray2D:: stringArray2D( stringArray2D const& other )
//----------------------------------------------------------------------
   : vector( other.vector )
   , d0( other.d0 )
   , d1( other.d1 )
{
   MAC_LABEL( "stringArray2D:: stringArray2D( stringArray2D const& )" ) ;
   MAC_CHECK_POST( *this == other ) ;
}

//----------------------------------------------------------------------
stringArray2D const&
stringArray2D:: operator=( stringArray2D const& other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringArray2D:: operator=" ) ;
   
   if( this != &other )
   {
      vector = other.vector ;
      d0 = other.d0 ;
      d1 = other.d1 ;
   }

   stringArray2D const& result = *this ;

   MAC_CHECK_POST( result == other ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
stringArray2D:: re_initialize( size_t dim0, size_t dim1, string val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringArray2D:: re_initialize" ) ;
   
   d0 = dim0 ;
   d1 = dim1 ;
   vector.re_initialize( d0*d1 ) ;
   vector.set( val ) ;

   MAC_CHECK_POST( index_bound(0) == dim0 ) ;
   MAC_CHECK_POST( index_bound(1) == dim1 ) ;
   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
bool
stringArray2D:: operator==( stringArray2D const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringArray2D:: operator==" ) ;

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
stringArray2D:: operator!=( stringArray2D const& other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringArray2D:: operator!=" ) ;

   bool result = !operator==( other ) ;

   MAC_CHECK_POST( EQUIVALENT( result, !operator==(other) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
stringArray2D:: raise_first_index_bound( size_t dim0, string val ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringArray2D:: raise_first_index_bound" ) ;
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
stringArray2D:: set( string val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringArray2D:: set" ) ;
   vector.set( val ) ;

   MAC_CHECK_POST( FORALL( ( size_t i0=0 ; i0<index_bound(0) ; ++i0 ),
                   FORALL( ( size_t i1=0 ; i1<index_bound(1) ; ++i1 ),
                   operator()(i0,i1) == val ) ) ) ;
}

//----------------------------------------------------------------------
void
stringArray2D:: set_section( size_t an_index, size_t index_value, 
                            stringVector const& x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringArray2D:: set_section" ) ;
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
stringArray2D:: extract_section( size_t an_index, size_t index_value,
                                 stringVector& x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "stringArray2D:: extract_section" ) ;
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
operator<<( std::ostream& out, stringArray2D const& a )
//----------------------------------------------------------------------
{
   MAC_LABEL( "operator<<( std::ostream&, stringArray2D const& )" ) ;

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
