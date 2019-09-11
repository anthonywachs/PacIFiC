#include <MAC_IntArray3D.hh>

#include <iostream>
#include <numeric>

#include <MAC_assertions.hh>

//----------------------------------------------------------------------------
MAC_IntArray3D*
MAC_IntArray3D:: create( MAC_Object* a_owner,
                          intArray3D const& val )
//----------------------------------------------------------------------------
{
   return( new MAC_IntArray3D( a_owner, val ) ) ;
}



//----------------------------------------------------------------------------
MAC_IntArray3D*
MAC_IntArray3D:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   return( new MAC_IntArray3D( a_owner, myValue ) ) ;
}



//----------------------------------------------------------------------------
MAC_IntArray3D:: MAC_IntArray3D( MAC_Object* a_owner,
                             intArray3D const& val )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner ),
     myValue( val )
{
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
MAC_IntArray3D:: ~MAC_IntArray3D( void )
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
MAC_Data::Type
MAC_IntArray3D:: data_type( void ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   return( MAC_Data::IntArray3D ) ;
}



//----------------------------------------------------------------------------
intArray3D const&
MAC_IntArray3D:: to_int_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_PRE( to_int_array3D_PRE(ct) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( myValue ) ;
}



//----------------------------------------------------------------------------
void
MAC_IntArray3D:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   std::string space( indent_width, ' ' ) ;
   MAC_CHECK_INV( invariant() ) ;
   os << space << "array( " << std::endl << space ;
   for( size_t i=0 ; i<myValue.index_bound(0) ; i++ )
   {
      os << space << "array( " << std::endl << space ;
      for( size_t j=0 ; j<myValue.index_bound(1) ; j++ )
      {
         os << "  < " ;
         for( size_t k=0 ; k<myValue.index_bound(2) ; k++ )
         {
            if( (k+1) % 10 == 0 )
            {
               os << std::endl << space ;
            }
            os << myValue(i,j,k) << " " ;
         }
         os << "  >" ;
         if( j!=myValue.index_bound(1)-1 )
         {
            os << "," ;
         }
      }
      os << "  ) " ;
      if( i!=myValue.index_bound(0)-1 )
      {
         os << "," ;
      }
      os << std::endl << space ;
   }
   os << " )" ;
}
