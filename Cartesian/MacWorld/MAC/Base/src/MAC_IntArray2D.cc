#include <MAC_IntArray2D.hh>

#include <iostream>
#include <numeric>

#include <MAC_assertions.hh>
#include <intVector.hh>
#include <MAC.hh>
#include <MAC_List.hh>

//----------------------------------------------------------------------------
MAC_IntArray2D*
MAC_IntArray2D:: create( MAC_Object* a_owner,
                         intArray2D const& val )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IntArray2D:: create(MAC_Object*,intArray2D const&)" ) ;

   MAC_IntArray2D* result = new MAC_IntArray2D( a_owner, val ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_IntArray2D*
MAC_IntArray2D:: create( MAC_Object* a_owner,
                         size_t_array2D const& val )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_IntArray2D:: create(MAC_Object*,intArray2D const&)" ) ;

   MAC_IntArray2D* result = new MAC_IntArray2D( a_owner, val ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_IntArray2D*
MAC_IntArray2D:: create( MAC_Object* a_owner,
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
   intArray2D array(0,0);
   if( !error )
   {   
      for( size_t i=0 ; !error && i<nb_rows ; i++ )
      {
         MAC_Data const* val =
            static_cast<MAC_Data const*>( list->at(i) ) ;
         error = val==0 || !val->value_can_be_evaluated(0) || val->data_type()!=IntVector ;
         if( !error )
         {
            intVector const& vec = val->to_int_vector() ;
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

   return ( error ? 0 : new MAC_IntArray2D( a_owner, array ) ) ;
   
}

//----------------------------------------------------------------------------
MAC_IntArray2D*
MAC_IntArray2D:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   return( new MAC_IntArray2D( a_owner, myValue ) ) ;
}

//----------------------------------------------------------------------------
MAC_IntArray2D:: MAC_IntArray2D( MAC_Object* a_owner,
                                 intArray2D const& val )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner )
   , myValue( val )
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_IntArray2D:: MAC_IntArray2D( MAC_Object* a_owner,
                                 size_t_array2D const& val )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner )
   , myValue( 0, 0 )
{
   myValue.set( val ) ;

   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_IntArray2D:: ~MAC_IntArray2D( void )
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
MAC_Data::Type
MAC_IntArray2D:: data_type( void ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   return( MAC_Data::IntArray2D ) ;
}



//----------------------------------------------------------------------------
intArray2D const&
MAC_IntArray2D:: to_int_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_PRE( to_int_array2D_PRE(ct) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( myValue ) ;
}



//----------------------------------------------------------------------------
void
MAC_IntArray2D:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   std::string space( indent_width, ' ' ) ;
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
         os << myValue(i,j) << " " ;
      }
      os << "  >" ;
      if( i!=myValue.index_bound(0)-1 )
      {
         os << "," ;
      }
      os << std::endl << space ;
   }
   os << " ]" ;
}
