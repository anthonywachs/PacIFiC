#include <MAC_Double.hh>

#include <iostream>

#include <MAC_assertions.hh>
#include <MAC.hh>
#include <MAC_DoubleComparatorExact.hh>

//----------------------------------------------------------------------------
MAC_Double*
MAC_Double:: create( MAC_Object* a_owner, double val )
//----------------------------------------------------------------------------
{
   return( new MAC_Double( a_owner, val ) ) ;
}

//----------------------------------------------------------------------------
MAC_Double*
MAC_Double:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   return( new MAC_Double( a_owner, VALUE ) ) ;
}

//----------------------------------------------------------------------------
MAC_Double:: MAC_Double( MAC_Object* a_owner, double val )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner )
   , VALUE( val )
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_Double:: ~MAC_Double( void )
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_Double*
MAC_Double:: create_derivative( MAC_Object* a_owner,
                                MAC_Variable const* var,
                                MAC_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Double:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   MAC_Double* result = create( a_owner, 0.0 ) ;
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DoubleComparator const*
MAC_Double:: double_comparator( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Double:: double_comparator" ) ;
   static MAC_DoubleComparator const* result =
                                   MAC_DoubleComparatorExact::object() ;
   MAC_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_Double:: is_equal( MAC_Object const* other ) const 
//----------------------------------------------------------------------
{
   MAC_CHECK_PRE( is_equal_PRE( other ) ) ;

   bool result = ( three_way_comparison( other ) == 0 ) ;

   MAC_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
int
MAC_Double:: three_way_comparison( MAC_Object const* other ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   static MAC_DoubleComparator const* DBL_COMP = double_comparator() ;
   MAC_Double const* otherDouble = static_cast<MAC_Double const*>( other ) ;
   int result = DBL_COMP->three_way_comparison( VALUE, otherDouble->VALUE )  ;

   MAC_CHECK_POST( three_way_comparison_POST( result, other ) ) ;      
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
MAC_Double:: hash_code( void ) const 
//----------------------------------------------------------------------
{
   return( (size_t) VALUE ) ;
}

//----------------------------------------------------------------------------
MAC_Data::Type
MAC_Double:: data_type( void ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   return( MAC_Data::Double ) ;
}

//----------------------------------------------------------------------------
double
MAC_Double:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_PRE( to_double_PRE( ct ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( VALUE ) ;
}

//----------------------------------------------------------------------------
void
MAC_Double:: set( double val )
//----------------------------------------------------------------------------
{
   VALUE = val ;
}

//----------------------------------------------------------------------------
void
MAC_Double:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   std::string space( indent_width, ' ' ) ;
   os << space ;
   MAC::print_double( os, VALUE ) ;
}
