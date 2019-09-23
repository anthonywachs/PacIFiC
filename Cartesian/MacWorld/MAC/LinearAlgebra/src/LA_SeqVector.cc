#include <LA_SeqVector.hh>

#include <LA_SeqImplementation.hh>
#include <LA_SeqScatter.hh>

#include <MAC.hh>
#include <MAC_DistributedPartition.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Error.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Int.hh>
#include <MAC_Vector.hh>

#include <fstream>
#include <iomanip>

#ifdef OUTLINE
#define inline
#include <LA_SeqVector.icc>
#undef inline
#endif


//----------------------------------------------------------------------
LA_SeqVector*
LA_SeqVector:: create( MAC_Object* a_owner, size_t a_nb_rows )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: create" ) ;

   LA_SeqVector* result = new LA_SeqVector( a_owner, a_nb_rows ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->nb_rows() == a_nb_rows ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_rows() ; ++i ),
                           result->item(i) == 0.0 ) ) ;
   MAC_CHECK_POST( result->state() == LA::Sync ) ;
   MAC_CHECK_POST( result->is_resizable() ) ;
   MAC_CHECK_POST( ! result->is_desynchronizable() ) ;
   MAC_CHECK_POST( result->is_synchronized() ) ;
   MAC_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
LA_SeqVector*
LA_SeqVector:: create( MAC_Object* a_owner, doubleVector const& dvec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: create" ) ;
   MAC_CHECK_PRE( dvec.size()>0 ) ;

   LA_SeqVector* result = new LA_SeqVector( a_owner, dvec.size() ) ;
   for( size_t i=0 ; i<dvec.size() ; ++i )
   {
      result->set_item( i, dvec(i) ) ;
   }

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->nb_rows() == dvec.size() ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<result->nb_rows() ; ++i ),
                           result->item(i) == dvec(i) ) ) ;
   MAC_CHECK_POST( result->state() == LA::Sync ) ;
   MAC_CHECK_POST( result->is_resizable() ) ;
   MAC_CHECK_POST( ! result->is_desynchronizable() ) ;
   MAC_CHECK_POST( result->is_synchronized() ) ;
   MAC_CHECK_POST( result->distribution_strategy() == LA::NoDistribution ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
LA_SeqVector:: LA_SeqVector( MAC_Object* a_owner, size_t a_nb_rows )
//----------------------------------------------------------------------
   : LA_Vector( a_owner, a_nb_rows )
   , UNSYNCHRO( false )
   , ROW_DIST( 0 )
   , OWNS_DATA( true )
   , DATA( 0 )
{
   MAC_LABEL( "LA_SeqVector:: LA_SeqVector" ) ;

   if( a_nb_rows>0)
   {
      DATA = new double [ a_nb_rows ] ;
      for( size_t i=0 ; i<a_nb_rows ; ++i ) DATA[i] = 0.0 ;
   }

   set_distribution_strategy( LA::NoDistribution ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( distribution_strategy() == LA::NoDistribution ) ;
}




//----------------------------------------------------------------------
LA_SeqVector*
LA_SeqVector:: create_vector( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: create_vector" ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_SeqVector* result = new LA_SeqVector( a_owner, nb_rows() ) ;
   if( UNSYNCHRO ) result->make_desynchronizable() ;

   MAC_CHECK_POST( create_vector_POST( result, a_owner ) ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_rows() ; ++i ),
                           result->item(i) == 0. ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
LA_SeqVector:: ~LA_SeqVector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: ~LA_SeqVector" ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( OWNS_DATA && DATA != 0 )
   {
      delete [] DATA ;
   }
   DATA = 0 ;
}




//----------------------------------------------------------------------
LA_Implementation const*
LA_SeqVector:: implementation( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: implementation" ) ;

   LA_Implementation const* result = LA_SeqImplementation::object() ;

   MAC_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
LA_SeqVector:: is_desynchronizable( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: is_desynchronizable" ) ;

   bool result = UNSYNCHRO ;

   return( result ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: make_desynchronizable( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: make_desynchronizable" ) ;
   MAC_CHECK_PRE( ! is_desynchronizable() ) ;

   UNSYNCHRO = true ;

   MAC_CHECK_POST( is_desynchronizable() ) ;
   MAC_CHECK_POST( state() == LA::Sync ) ;
   MAC_CHECK_POST( is_synchronized() ) ;
}




//----------------------------------------------------------------------
MAC_DistributedPartition const*
LA_SeqVector:: row_distribution( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: row_distribution" ) ;
   MAC_CHECK_PRE( row_distribution_PRE() ) ;

   if( ROW_DIST == 0 )
   {
      ROW_DIST = MAC_DistributedPartition::create(
                                 const_cast<LA_SeqVector*>( this ) ) ;
   }
   if( ROW_DIST->global_number() != nb_rows() )
   {
      ROW_DIST->set_global_number( nb_rows() ) ;
   }
   MAC_DistributedPartition const* result = ROW_DIST ;

   MAC_CHECK_POST( row_distribution_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: re_initialize( size_t a_nb_rows, size_t a_nb_local_rows )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: re_initialize" ) ;
   MAC_CHECK_PRE( re_initialize_PRE( a_nb_rows, a_nb_local_rows ) ) ;

   re_initialize_with_global_sizes( a_nb_rows ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( re_initialize_POST( a_nb_rows, a_nb_local_rows ) ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_rows() ; ++i ),
                           item(i) == 0.0 ) ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: re_initialize_with_global_sizes( size_t a_nb_rows )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: re_initialize_with_global_sizes(size_t)" ) ;
   MAC_CHECK_PRE( re_initialize_with_global_sizes_PRE( a_nb_rows ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( OWNS_DATA ) ;

   if( a_nb_rows != nb_rows() )
   {
      if( DATA != 0 )
      {
         delete [] DATA ;
         DATA = 0 ;
      }
      if( a_nb_rows>0)
      {
         DATA = new double [ a_nb_rows ] ;
         for( size_t i=0 ; i<a_nb_rows ; ++i ) DATA[i] = 0.0 ;
      }
      set_rows_number( a_nb_rows ) ;
      if( UNSYNCHRO )
      {
         synchronize() ;
      }
   }
   else
   {
      nullify() ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( re_initialize_with_global_sizes_POST( a_nb_rows ) ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<nb_rows() ; ++i ),
                           item(i) == 0.0 ) ) ;
}




//----------------------------------------------------------------------
LA_SeqScatter*
LA_SeqVector:: create_scatter( MAC_Object* a_owner,
                               size_t_vector const& from,
                               size_t_vector const& to ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: create_scatter" ) ;
   MAC_CHECK_PRE( create_scatter_PRE( a_owner, from, to ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   LA_SeqScatter* result =
                LA_SeqScatter::create( a_owner, nb_rows(), from, to ) ;

   MAC_CHECK_POST( create_scatter_POST( result, a_owner, from, to ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
LA_SeqVector*
LA_SeqVector:: create_local_vector( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: create_local_vector" ) ;
   MAC_CHECK_PRE( create_local_vector_PRE( a_owner ) ) ;

   LA_SeqVector* result = create_vector( a_owner ) ;
   result->set( this ) ;

   MAC_CHECK_POST( create_local_vector_POST( result, a_owner ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: set( LA_Vector const* a )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: set" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( set_PRE( a ) ) ;

   MAC_CHECK( dynamic_cast<LA_SeqVector const*>(a) != 0 ) ;
   LA_SeqVector const* ba = static_cast<LA_SeqVector const*>(a) ;

   if( this != a )
   {
      for( size_t i=0 ; i<nb_rows() ; ++i ) DATA[i] = ba->DATA[i] ;
   }
   synchronize() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_POST( a ) ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: set( double value )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: set double" ) ;
   MAC_CHECK_PRE( set_PRE( value ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<nb_rows() ; ++i ) DATA[i] = value ;
   synchronize() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_POST( value ) ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: nullify( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: nullify" ) ;
   MAC_CHECK_PRE( nullify_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<nb_rows() ; ++i ) DATA[i] = 0.0 ;
   synchronize() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( nullify_POST() ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: scale( double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: scale" ) ;
   MAC_CHECK_PRE( scale_PRE( alpha ) ) ;

   if( alpha == 0. )
   {
      for( size_t i=0 ; i<nb_rows() ; ++i )
      {
         DATA[i] = 0. ;
      }
   }
   else if( alpha != 1. )
   {
      for( size_t i=0 ; i<nb_rows() ; ++i )
      {
         DATA[i] *= alpha ;
      }
   }

   MAC_CHECK_POST( scale_POST( alpha ) ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: sum( LA_Vector const* a, double alpha )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: sum" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( sum_PRE( a, alpha ) ) ;

   MAC_CHECK( dynamic_cast<LA_SeqVector const*>( a ) != 0 ) ;
   LA_SeqVector const* ba = static_cast<LA_SeqVector const*>( a ) ;

   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      DATA[i] += alpha*ba->DATA[i] ;
   }

   MAC_CHECK_POST( sum_POST( a, alpha ) ) ;
}




//----------------------------------------------------------------------
double
LA_SeqVector:: dot( LA_Vector const* a ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: dot" ) ;
   MAC_CHECK_PRE( dot_PRE(a) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK( dynamic_cast<LA_SeqVector const*>(a) != 0 ) ;
   LA_SeqVector const* ba = static_cast<LA_SeqVector const*>(a) ;

   double result = 0. ;
   for( size_t i=0 ; i<a->nb_rows() ; ++i )
   {
      result += DATA[i] * ba->DATA[i] ;
   }

   MAC_CHECK_POST( dot_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
double
LA_SeqVector:: two_norm( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: two_norm" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( two_norm_PRE() ) ;

   double result = 0. ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      result += DATA[i] * DATA[i] ;
   }
   result = MAC::sqrt( result ) ;

   MAC_CHECK_POST( two_norm_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
double
LA_SeqVector:: max_norm( void ) const
//----------------------------------------------------------------------
// max_norm = max( |v( i )| )
{
   MAC_LABEL( "LA_SeqVector:: max_norm" ) ;
   MAC_CHECK_PRE( max_norm_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   double result = 0. ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      if( DATA[i]>result )
      {
         result = DATA[i] ;
      }
      else if( -DATA[i]>result )
      {
         result = -DATA[i] ;
      }
   }
   MAC_CHECK_POST( max_norm_POST( result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
double
LA_SeqVector:: mean( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: mean" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( nb_rows() > 0 ) ;

   double result = 0. ;
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      result += DATA[i] ;
   }
   return( result/double( nb_rows() ) ) ;
}




//----------------------------------------------------------------------
double
LA_SeqVector:: standard_deviation( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: standard_deviation" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( nb_rows() > 1 ) ;
   double result = 0.0 ;

   if( nb_rows() > 1 )
   {
      double variance = 0.0 ;
      double moyenne =  mean() ;
      double elem ;
      for( size_t i=0 ; i<nb_rows() ; ++i )
      {
         elem = DATA[i] - moyenne ;
         variance += elem*elem ;
      }
      variance /= double( nb_rows() )  ;
      result = MAC::sqrt( variance ) ;
   }

   return( result ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: print_items( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: print_items" ) ;
   MAC_CHECK_PRE( print_items_PRE( os, indent_width ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   std::string const space( indent_width, ' ' ) ;
   os << space << "nb_rows:" << nb_rows() << std::endl ;
   if( nb_rows()>0 )
   {
      std::ios_base::fmtflags original_flags = os.flags() ;
      os.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
      std::streamsize p = os.precision() ;
      os << std::setprecision( 7 ) ;
      for( size_t iRow = 0 ; iRow<nb_rows() ; ++iRow )
      {
         os << space << "Row n°" << iRow << "  " ;
         double const x = DATA[iRow] ;
         os << std::setw(15) << x << std::endl ;
      }
      os << std::setprecision(p) ;
      os.flags( original_flags ) ;
   }
}




//----------------------------------------------------------------------
void
LA_SeqVector:: write( std::string const& filename ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: write" ) ;
   MAC_CHECK_INV( invariant() ) ;
   std::ofstream file( filename.c_str() ) ;
//   MAC_ASSERT( file ) ; // Not accepted from gcc-9.x.x
   MAC_ASSERT( file.is_open() ) ;
   file.precision( 15 ) ;

   file << nb_rows() ;
   for ( size_t iRow=0 ; iRow<nb_rows() ; iRow++ )
   {
      file << std::endl << DATA[iRow] ;
   }
}




//----------------------------------------------------------------------
void
LA_SeqVector:: save( MAC_Module* mod ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: save" ) ;
   MAC_CHECK_PRE( is_synchronized() ) ;
   MAC_CHECK_PRE( mod != 0 ) ;

   mod->add_entry( "nb_rows",
                   MAC_Int::create( mod, (int)nb_rows() ) ) ;
   doubleVector dvec( nb_rows() );
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      dvec(i)=DATA[i] ;
   }

   mod->add_entry( "data", MAC_DoubleVector::create( mod, dvec ) ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: restore( MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: restore" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   re_initialize_with_global_sizes( exp->int_data("nb_rows") ) ;
   doubleVector const& v = exp->doubleVector_data("data");
   for( size_t i=0 ; i<nb_rows() ; ++i )
   {
      DATA[i]=v(i) ;
   }
   synchronize() ;

   MAC_CHECK_POST( is_synchronized() ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: set_as_v_product( LA_Vector const* a, LA_Vector const* b )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: set_as_v_product" ) ;
   MAC_CHECK_PRE( set_as_v_product_PRE( a, b ) ) ;

   MAC_CHECK( dynamic_cast<LA_SeqVector const*>(a) != 0 ) ;
   LA_SeqVector const* ba = static_cast<LA_SeqVector const*>(a) ;
   MAC_CHECK( dynamic_cast<LA_SeqVector const*>(b) != 0 ) ;
   LA_SeqVector const* bb = static_cast<LA_SeqVector const*>(b) ;

   for ( size_t iRow=0 ; iRow<nb_rows() ; iRow++ )
   {
      DATA[iRow] = ba->DATA[iRow] * bb->DATA[iRow] ;
   }
   synchronize() ;

   MAC_CHECK_POST( set_as_v_product_POST( a, b ) ) ;
}




//----------------------------------------------------------------------
void
LA_SeqVector:: set_as_reciprocal( LA_Vector const* a,
                                  double smallest_inverted_item,
                                  double default_value )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqVector:: set_as_reciprocal" ) ;
   MAC_CHECK_PRE( set_as_reciprocal_PRE( a, smallest_inverted_item, 
   	default_value ) ) ;

   MAC_CHECK( dynamic_cast<LA_SeqVector const*>(a) != 0 ) ;
   LA_SeqVector const* ba = static_cast<LA_SeqVector const*>(a) ;

   for ( size_t iRow=0 ; iRow<nb_rows() ; iRow++ )
   {
      DATA[iRow] = ( MAC::abs( ba->DATA[iRow] ) < smallest_inverted_item ?
                     default_value :
                     1.0 / ba->DATA[iRow] ) ;
   }
   synchronize() ;

   MAC_CHECK_POST( set_as_reciprocal_POST( a ) ) ;
}




//----------------------------------------------------------------------
bool
LA_SeqVector:: re_initialize_with_global_sizes_PRE( size_t a_nb_rows ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_resizable() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_SeqVector:: re_initialize_with_global_sizes_POST( size_t a_nb_rows ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( nb_rows() == a_nb_rows ) ;
   MAC_ASSERT( state() == LA::Sync ) ;
   MAC_ASSERT( is_synchronized() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_SeqVector:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Vector::invariant() ) ;
   MAC_ASSERT( distribution_strategy() == LA::NoDistribution ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
LA_SeqVector:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( LA_Vector::implementation_POST( result ) ) ;
   MAC_ASSERT( result == LA_SeqImplementation::object() ) ;
   return( true ) ;
}
