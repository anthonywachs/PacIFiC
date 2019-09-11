#include <MAC_DoubleVector.hh>

#include <iostream>

#include <MAC_assertions.hh>
#include <MAC_Container.hh>
#include <MAC_Iterator.hh>
#include <doubleVector.hh>

//----------------------------------------------------------------------------
MAC_DoubleVector*
MAC_DoubleVector:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleVector:: create_clone" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_DoubleVector* result = 0 ;
   result = new MAC_DoubleVector( a_owner, dv ) ;
   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------------
MAC_DoubleVector*
MAC_DoubleVector:: create( MAC_Object* a_owner,
                           doubleVector const& aDoubleVector )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleVector:: create(doubleVector)" ) ;
   MAC_DoubleVector* result = new MAC_DoubleVector( a_owner, aDoubleVector ) ;
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->to_double_vector().size() == aDoubleVector.size() ) ;
   MAC_CHECK_POST(
      FORALL( ( size_t i=0 ; i<aDoubleVector.size() ; ++i ),
              result->to_double_vector()(i) == aDoubleVector(i) ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------------
MAC_DoubleVector:: MAC_DoubleVector( MAC_Object* a_owner,
                                     doubleVector const& aDoubleVector )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner )
   , dv( aDoubleVector )
{
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------------
MAC_DoubleVector:: MAC_DoubleVector( MAC_Object* a_owner,
                                     MAC_Container const* aDataList,
                                     bool& error )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner )
   , dv( aDataList->count() )
{
   MAC_Iterator* list_iterator = aDataList->create_iterator( 0 ) ;
   size_t cpt=0 ;
   for( list_iterator->start() ;
        list_iterator->is_valid() && !error ;
        list_iterator->go_next() )
   {
      MAC_Data const* val =
         dynamic_cast<MAC_Data const*>( list_iterator->item() ) ;
      error |= ( val==0  ||
                 val->data_type()!=Double ||
                 !val->value_can_be_evaluated( 0 ) ) ;
      if( !error ) dv(cpt++) = val->to_double() ;
   }
   
   list_iterator->destroy() ; list_iterator=0 ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------------
MAC_DoubleVector*
MAC_DoubleVector:: create( MAC_Object* a_owner,
                           MAC_Container const* aDataList ) 
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleVector:: create(MAC_Container)" ) ;
   MAC_CHECK_PRE( aDataList!=0 ) ;
   MAC_CHECK_PRE( aDataList->count()!=0 ) ;

   bool error = false ;
   MAC_DoubleVector* result = new MAC_DoubleVector( 0, aDataList, error ) ;
   if( error )
   {
      result->destroy() ;
      result = 0 ;
   }
   else if( a_owner != 0 )
   {
      result->set_owner( a_owner ) ;
   }

   MAC_CHECK_POST( IMPLIES( result != 0, result->owner() == a_owner ) ) ;
   MAC_CHECK_POST( IMPLIES( result != 0, result->to_double_vector().size() 
   	== aDataList->count() ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------------
MAC_DoubleVector:: ~MAC_DoubleVector( void )
//----------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------------
void
MAC_DoubleVector:: set( doubleVector const& other )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleVector:: set" ) ;
   dv=other ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<other.size() ; i++ ),
                           to_double_vector()(i)==other(i) ) ) ;
}




//----------------------------------------------------------------------------
MAC_Data::Type
MAC_DoubleVector:: data_type( void ) const
//----------------------------------------------------------------------------
{
   return( DoubleVector ) ;
}




//----------------------------------------------------------------------------
doubleVector const& 
MAC_DoubleVector:: to_double_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleVector:: to_double_vector" ) ;
   MAC_CHECK_PRE( to_double_vector_PRE(ct) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( dv ) ;
}




//----------------------------------------------------------------------------
void
MAC_DoubleVector:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleVector:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   std::ios::fmtflags oldoptions = os.flags( std::ios::scientific ) ;
   MAC_CHECK_INV( invariant() ) ;
   os << space << "< " ;
   for( size_t i=0 ; i<dv.size() ; i++ )
   {
      if( (i+1) % 10 == 0 ) os << std::endl << space ;
      MAC::print_double( os, dv(i) ) ;
      os << " " ;
   }
   os << "> " ;
   os.flags( oldoptions ) ;
}




//----------------------------------------------------------------------
MAC_Data*
MAC_DoubleVector:: create_derivative( MAC_Object* a_owner,
                                      MAC_Variable const* var,
                                      MAC_Context const* ct ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleVector:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;

   doubleVector adv(dv.size()) ;
   adv.set( 0.0 ) ;
   
   MAC_DoubleVector* result = create( a_owner, adv ) ;
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
bool
MAC_DoubleVector:: invariant( void ) const 
//-------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Data::invariant() ) ;
   return( true ) ;
}
