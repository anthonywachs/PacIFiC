#include <MAC_ConditionalExp.hh>

#include <MAC_assertions.hh>
#include <MAC_List.hh>
#include <MAC_Vector.hh>
#include <MAC_Root.hh>
#include <MAC_Sequence.hh>
#include <MAC_String.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <iostream>

MAC_ConditionalExp const*
MAC_ConditionalExp::PROTOTYPE = new MAC_ConditionalExp() ;

//----------------------------------------------------------------------
MAC_ConditionalExp:: MAC_ConditionalExp( void ) 
//----------------------------------------------------------------------
      : MAC_TransferExp( "(?:)" )
      , DEFAULT(0)
{
   MAC_LABEL( "MAC_ConditionalExp:: MAC_ConditionalExp" ) ;
}

//----------------------------------------------------------------------
MAC_ConditionalExp:: MAC_ConditionalExp( MAC_Object* a_owner,
                                         MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
      : MAC_TransferExp( a_owner, "(?:)", argument_list )
      , DEFAULT( arg( nb_arguments()-1 ) )
{
   MAC_LABEL( "MAC_ConditionalExp:: MAC_ConditionalExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ConditionalExp:: ~MAC_ConditionalExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConditionalExp:: ~MAC_ConditionalExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ConditionalExp*
MAC_ConditionalExp:: create_replica( MAC_Object* a_owner,
                             MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConditionalExp:: create_replica" ) ;
   MAC_CHECK_PRE( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_ConditionalExp* result = new MAC_ConditionalExp( a_owner, 
                                                      argument_list ) ;

   MAC_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data const*
MAC_ConditionalExp:: data( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConditionalExp:: data" ) ;
   MAC_CHECK( data_PRE( ct ) ) ;
   
   MAC_Data const* result = 0 ;
   size_t nb_tests = (nb_arguments()-1)/2 ;
   for( size_t i=0 ; i<nb_tests ; i++ )
   {
      size_t idx = 2*i ;
      if( arg(idx)->to_bool(ct) ) result=arg(idx+1) ;
   }
   if( result==0 ) result=DEFAULT ;
   
   MAC_CHECK( data_POST( result, ct ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data*
MAC_ConditionalExp:: create_derivative( MAC_Object* a_owner,
                                        MAC_Variable const* var,
                                        MAC_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConditionalExp:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   MAC_Data* result = 0 ;
   MAC_List* list = MAC_List::create( 0 ) ;
   size_t nb_tests = (nb_arguments()-1)/2 ;
   for( size_t i=0 ; i<nb_tests ; i++ )
   {
      size_t idx = 2*i ;
      list->append( arg(idx)->create_clone( list ) ) ;
      list->append( arg(idx+1)->create_derivative( list, var, ct ) ) ;
   }
   list->append( DEFAULT->create_derivative( list, var, ct ) ) ;
      
   result = MAC_Expression::create( a_owner, name(), list ) ;
   list->set_owner( result ) ;
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_ConditionalExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result =
      "( BS ? <val1> : ... BS ? <valN> : <defaultVal> )" ;
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_ConditionalExp:: valid_arguments(
                              MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConditionalExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()>=3 &&
      ( some_arguments->count()%2 ) == 1;
   if( result )
   {
      MAC_Data::Type k_default =
         extract_arg( some_arguments, some_arguments->count()-1 )->data_type() ;
      size_t nb_tests = (some_arguments->count()-1)/2 ;
      for( size_t i=0 ; i<nb_tests ; i++ )
      {
         size_t idx = 2*i ;
         MAC_Data::Type k0 =  extract_arg( some_arguments, idx )->data_type() ;
         MAC_Data::Type k1 =  extract_arg( some_arguments, idx+1 )->data_type() ;
         result = result && k0==Bool && k1==k_default ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_ConditionalExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   Type result = DEFAULT->data_type() ;
   return result ;
}

//----------------------------------------------------------------------
void
MAC_ConditionalExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space << "( " ;
   size_t nb_tests = (nb_arguments()-1)/2 ;
   for( size_t i=0 ; i<nb_tests ; i++ )
   {
      size_t idx = 2*i ;
      arg(idx)->print(os,0) ;
      os << " ? " ;
      arg(idx+1)->print(os,0) ;
      os << " : " ;
   }
   DEFAULT->print( os, 0 ) ;
   os << " ) " ;
}
