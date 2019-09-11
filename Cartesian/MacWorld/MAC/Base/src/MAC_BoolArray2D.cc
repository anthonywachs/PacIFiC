#include <MAC_BoolArray2D.hh>

#include <iostream>
#include <numeric>

#include <boolVector.hh>
#include <MAC_assertions.hh>
#include <MAC.hh>
#include <MAC_List.hh>

//----------------------------------------------------------------------------
MAC_BoolArray2D*
MAC_BoolArray2D:: create( MAC_Object* a_owner,
                          boolArray2D const& val )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BoolArray2D:: create" ) ;

   MAC_BoolArray2D* result = new MAC_BoolArray2D( a_owner, val ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->to_bool_array2D() == val ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_BoolArray2D:: MAC_BoolArray2D( MAC_Object* a_owner,
                                   boolArray2D const& val )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner )
   , MY_VALUE( val )
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_BoolArray2D*
MAC_BoolArray2D:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BoolArray2D:: create_clone" ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_BoolArray2D* result = new MAC_BoolArray2D( a_owner, MY_VALUE ) ;

   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_BoolArray2D:: ~MAC_BoolArray2D( void )
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_BoolArray2D*
MAC_BoolArray2D:: create( MAC_Object* a_owner,
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
   boolArray2D array(0,0);
   if( !error )
   {   
      for( size_t i=0 ; !error && i<nb_rows ; i++ )
      {
         MAC_Data const* val =
            static_cast<MAC_Data const*>( list->at(i) ) ;
         error = val==0 || !val->value_can_be_evaluated(0) || val->data_type()!=BoolVector ;
         if( !error )
         {
            boolVector const& vec = val->to_bool_vector() ;
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

   return ( error ? 0 : new MAC_BoolArray2D( a_owner, array ) ) ;
   
}

//----------------------------------------------------------------------------
MAC_Data::Type
MAC_BoolArray2D:: data_type( void ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   return( MAC_Data::BoolArray2D ) ;
}

//----------------------------------------------------------------------------
boolArray2D const&
MAC_BoolArray2D:: to_bool_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( MY_VALUE ) ;
}

//----------------------------------------------------------------------------
void
MAC_BoolArray2D:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   std::string space( indent_width, ' ' ) ;
   std::ios::fmtflags oldoptions = os.flags( std::ios::scientific ) ;
   MAC_CHECK_INV( invariant() ) ;
   os << space << "[ " << std::endl << space ;
   for( size_t i=0 ; i<MY_VALUE.index_bound(0) ; i++ )
   {
      os << "< " ;
      for( size_t j=0 ; j<MY_VALUE.index_bound(1) ; j++ )
      {
         if( (j+1) % 10 == 0 )
         {
            os << std::endl << space ;
         }
         if( MY_VALUE(i,j) )
         {
            os << "true " ;
         }
         else
         {
            os << "false " ;
         }
         os << " " ;
      }
      os << ">" ;
      if( i!=MY_VALUE.index_bound(0)-1 )
      {
         os << "," ;
      }
      os << std::endl << space ;
   }
   os.flags( oldoptions ) ;
   os << " ]" ;
}
