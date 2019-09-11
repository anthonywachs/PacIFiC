#include <MAC_SortExp.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_BalancedBinaryTree.hh>
#include <MAC_Double.hh>
#include <MAC_Error.hh>
#include <MAC_Vector.hh>
#include <MAC_Root.hh>
#include <MAC_Sequence.hh>
#include <MAC_String.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <iostream>

MAC_SortExp const* MAC_SortExp::PROTOTYPE = new MAC_SortExp() ;

//----------------------------------------------------------------------
MAC_SortExp:: MAC_SortExp( void ) 
//----------------------------------------------------------------------
   : MAC_Expression( "sort" )
   , VECTOR( 0 )
   , GROWING( false )
   , RESULT( 0 )
{
   MAC_LABEL( "MAC_SortExp:: MAC_SortExp" ) ;
}

//----------------------------------------------------------------------
MAC_SortExp*
MAC_SortExp:: create_replica( MAC_Object* a_owner,
                                 MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SortExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_SortExp* result = new MAC_SortExp( a_owner, argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_SortExp:: MAC_SortExp( MAC_Object* a_owner,
                                 MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, "sort", argument_list )
   , VECTOR( arg(0) )
   , GROWING( arg(1)->to_string() == "<" )
   , RESULT( 0 )
{
   MAC_LABEL( "MAC_SortExp:: MAC_SortExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_SortExp:: ~MAC_SortExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SortExp:: ~MAC_SortExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_SortExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "sort(DV,\">\"|\"<\")" ;
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_SortExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SortExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()==2 &&
      extract_arg( some_arguments, 0 )->data_type()==DoubleVector &&
      extract_arg( some_arguments, 1 )->data_type()==String ;
   if( result )
   {
      std::string const& opt = extract_arg( some_arguments, 1 )->to_string() ;
      result = opt=="<" || opt==">" ;
   }
   return result ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_SortExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SortExp:: data_type" ) ;
   return MAC_Data::DoubleVector ;
}

//----------------------------------------------------------------------
doubleVector const&
MAC_SortExp:: to_double_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SortExp:: to_double_vector" ) ;
   MAC_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   doubleVector const& initial = VECTOR->to_double_vector( ct ) ;
   MAC_BalancedBinaryTree* tree =  MAC_BalancedBinaryTree::create( 0 ) ;
   for( size_t j=0 ; j<initial.size() ; j++ )
   {
      MAC_Double* d = MAC_Double::create( tree, initial(j) ) ;
      tree->extend( d) ;
   }
   doubleVector& result = RESULT ;
   size_t n = tree->count() ;
   result.re_initialize( n ) ;
   MAC_Iterator* it = tree->create_iterator( tree ) ;
   size_t i = 0 ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      size_t idx = ( GROWING ? i : n-1-i ) ;
      result( idx ) = static_cast<MAC_Data*>( it->item() )->to_double() ;
      i++ ;
   }
   tree->destroy() ;
   return result ;
}
