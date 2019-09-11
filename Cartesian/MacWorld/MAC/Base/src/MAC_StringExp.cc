#include <MAC_StringExp.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_Vector.hh>
#include <MAC_Sequence.hh>
#include <MAC_String.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <iostream>
#include <sstream>

MAC_StringExp const*
MAC_StringExp:: PROTOTYPE_EMPTY = new MAC_StringExp( "empty" ) ;

MAC_StringExp const*
MAC_StringExp:: PROTOTYPE_TO_STRING = new MAC_StringExp( "to_string" ) ;


//----------------------------------------------------------------------
MAC_StringExp:: MAC_StringExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
{
   MAC_LABEL( "MAC_StringExp:: MAC_StringExp" ) ;
}

//----------------------------------------------------------------------
MAC_StringExp:: MAC_StringExp( MAC_Object* a_owner,
                               std::string const& a_name,
                               MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
{
   MAC_LABEL( "MAC_StringExp:: MAC_StringExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_StringExp:: ~MAC_StringExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_StringExp:: ~MAC_StringExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_StringExp*
MAC_StringExp:: create_replica( MAC_Object* a_owner,
                                MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_StringExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_StringExp* result = 
                  new MAC_StringExp( a_owner, name(), argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_StringExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_StringExp:: data_type" ) ;
   MAC_Data::Type result = ( name()=="empty" ? Bool : String ) ;
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_StringExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_StringExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   if( name()=="to_string" )
   {
      result = ( some_arguments->count() == 1 ) ;
      if( result )
      {
         MAC_Data::Type t = extract_arg( some_arguments, 0 )->data_type() ;
         result = result && ( t==Int || t==Double ) ;
      }
   }
   else if( name()=="empty" )
   {
      result = ( some_arguments->count() == 1 ) ;
      if( result )
      {
         result = result &&
         ( extract_arg( some_arguments, 0 )->data_type() == String ) ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
std::string const&
MAC_StringExp:: usage( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_StringExp:: usage" ) ;

   static std::string result ;
   if( name()=="to_string" )
   {
      result = "to_string(DS|IS)" ;
   }
   else if( name()=="empty" )
   {
      result = "empty(SS)" ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_StringExp:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_StringExp:: to_bool" ) ;
   MAC_CHECK_PRE( to_bool_PRE( ct ) ) ;
   bool result = false ;
   
   if( name()=="empty" )
   {
      std::string const& str = arg(0)->to_string( ct ) ;
      result = str.empty() ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
std::string const&
MAC_StringExp:: to_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_StringExp:: to_string" ) ;
   MAC_CHECK_PRE( to_string_PRE(ct) ) ;
   if( name()=="to_string" )
   {
      std::ostringstream os ;
      MAC_Data::Type t = arg(0)->data_type() ;
      if( t==Double )
      {
         os << std::scientific << arg(0)->to_double( ct ) ;
      }
      else
      {
         os << arg(0)->to_int( ct ) ;
      }
      STR_RESULT = os.str() ;
   }
   
   return STR_RESULT ;
}
