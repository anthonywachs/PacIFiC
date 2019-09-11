#include <MAC_DoubleArray2D.hh>
#include <MAC_List.hh>
#include <doubleVector.hh>

#include <iostream>
#include <numeric>

#include <MAC_assertions.hh>

//----------------------------------------------------------------------------
MAC_DoubleArray2D*
MAC_DoubleArray2D:: create( MAC_Object* a_owner,
                          doubleArray2D const& val )
//----------------------------------------------------------------------------
{
   return( new MAC_DoubleArray2D( a_owner, val ) ) ;
}

//----------------------------------------------------------------------------
MAC_DoubleArray2D*
MAC_DoubleArray2D:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   return( new MAC_DoubleArray2D( a_owner, myValue ) ) ;
}

//----------------------------------------------------------------------------
MAC_DoubleArray2D*
MAC_DoubleArray2D:: create( MAC_Object* a_owner,
                            MAC_List const* list )
//----------------------------------------------------------------------------
{
   MAC_CHECK_PRE( list!=0 ) ;
   MAC_CHECK_PRE( list->count()!=0 ) ;
   MAC_CHECK_PRE( FORALL( (size_t i=0 ; i<list->count() ; i++ ),
                          dynamic_cast<MAC_Data const*>( list->at(i) )!=0 ) ) ;

   size_t nb_rows = list->count() ;
   size_t nb_cols = MAC::bad_index() ;
   bool error = nb_rows==0 ;
   doubleArray2D array(0,0);
   if( !error )
   {   
      for( size_t i=0 ; !error && i<nb_rows ; i++ )
      {
         MAC_Data const* val =
            static_cast<MAC_Data const*>( list->at(i) ) ;
         error = val==0 || !val->value_can_be_evaluated(0) || val->data_type()!=DoubleVector ;
         if( !error )
         {
            doubleVector const& vec = val->to_double_vector() ;
            if( i==0 )
            {
               nb_cols = vec.size() ;
               error = nb_cols==0 ;
               if( !error )
               {
                  array.re_initialize(nb_rows,nb_cols);
               }
            }
            error = error || vec.size()!=nb_cols ;
            for( size_t j=0 ; !error && j<nb_cols ; j++ )
            {
               array(i,j) = vec(j) ;
            }
         }
      }
   }

   return ( error ? 0 : new MAC_DoubleArray2D( a_owner, array ) ) ;
   
}

//----------------------------------------------------------------------------
MAC_DoubleArray2D:: MAC_DoubleArray2D( MAC_Object* a_owner,
                                   doubleArray2D const& val )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner ),
     myValue( val )
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_DoubleArray2D:: ~MAC_DoubleArray2D( void )
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_Data::Type
MAC_DoubleArray2D:: data_type( void ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   return( MAC_Data::DoubleArray2D ) ;
}

//----------------------------------------------------------------------------
doubleArray2D const&
MAC_DoubleArray2D:: to_double_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( myValue ) ;
}

//----------------------------------------------------------------------------
void
MAC_DoubleArray2D:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   std::string space( indent_width, ' ' ) ;
   std::ios::fmtflags oldoptions = os.flags( std::ios::scientific ) ;
   MAC_CHECK_INV( invariant() ) ;
   os << space << "[ " << std::endl << space ;
   for( size_t i=0 ; i<myValue.index_bound(0) ; i++ )
   {
      os << "< " ;
      for( size_t j=0 ; j<myValue.index_bound(1) ; j++ )
      {
         if( (j+1) % 10 == 0 )
         {
	   os << std::endl << space ;
         }
         MAC::print_double( os, myValue(i,j) ) ;
         os << " " ;
      }
      os << ">" ;
      if( i!=myValue.index_bound(0)-1 )
      {
         os << "," ;
      }
      os << std::endl << space ;
   }
   os.flags( oldoptions ) ;
   os << " ]" ;
}
