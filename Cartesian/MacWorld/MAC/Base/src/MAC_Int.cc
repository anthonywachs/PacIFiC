#include <MAC_Int.hh>

#include <numeric>
#include <iostream>

#include <MAC_assertions.hh>

//----------------------------------------------------------------------------
MAC_Int*
MAC_Int:: create( MAC_Object* a_owner, int val )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Int:: create" ) ;
   return( new MAC_Int( a_owner, val ) ) ;
}



//----------------------------------------------------------------------------
MAC_Int*
MAC_Int:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Int:: create_clone" ) ;
   return( new MAC_Int( a_owner, myValue ) ) ;
}



//----------------------------------------------------------------------------
MAC_Int:: MAC_Int( MAC_Object* a_owner, int val )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner ),
     myValue( val )
{
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
MAC_Int:: ~MAC_Int( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Int:: ~MAC_Int" ) ;
   MAC_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
MAC_Data::Type
MAC_Int:: data_type( void ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Int:: data_type" ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( MAC_Data::Int ) ;
}



//----------------------------------------------------------------------------
int
MAC_Int:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Int:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE(ct) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( myValue ) ;
}



//----------------------------------------------------------------------------
void
MAC_Int:: set( int val )
//----------------------------------------------------------------------------
{
   myValue = val ;
}



//----------------------------------------------------------------------------
void
MAC_Int:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Int:: print" ) ;
   MAC_CHECK_INV( invariant() ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << myValue ;
}
