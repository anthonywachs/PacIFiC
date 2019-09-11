#include <MAC_ConvertTypeExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Root.hh>
#include <MAC_Sequence.hh>
#include <MAC.hh>

#include <iostream>
#include <cmath>

MAC_ConvertTypeExp const* 
MAC_ConvertTypeExp::PROTOTYPE_DOUBLE = new MAC_ConvertTypeExp( "double" ) ;

MAC_ConvertTypeExp const* 
MAC_ConvertTypeExp::PROTOTYPE_INT = new MAC_ConvertTypeExp( "int" ) ;


//----------------------------------------------------------------------
MAC_ConvertTypeExp:: MAC_ConvertTypeExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , P0( 0 )
{
   MAC_LABEL( "MAC_ConvertTypeExp:: MAC_ConvertTypeExp" ) ;
}

//----------------------------------------------------------------------
MAC_ConvertTypeExp*
MAC_ConvertTypeExp:: create_replica( MAC_Object* a_owner,
                             MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConvertTypeExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_ConvertTypeExp* result = new MAC_ConvertTypeExp( a_owner,
                                                        name(),
                                                        argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ConvertTypeExp:: MAC_ConvertTypeExp( MAC_Object* a_owner,
                                         std::string const& a_name,
                                         MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , P0( arg(0) )
{
   MAC_LABEL( "MAC_ConvertTypeExp:: MAC_ConvertTypeExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ConvertTypeExp:: ~MAC_ConvertTypeExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConvertTypeExp:: ~MAC_ConvertTypeExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_ConvertTypeExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConvertTypeExp:: data_type" ) ;
   
   return ( name()=="double" ? MAC_Data::Double : MAC_Data::Int ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_ConvertTypeExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   result = ( (name()=="double") ? "double(IS)" : "int(DS)" ) ;
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_ConvertTypeExp:: valid_arguments(
                              MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConvertTypeExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = some_arguments->count()==1 ;
   if( result )
   {
      MAC_Data::Type k0 =  extract_arg( some_arguments, 0 )->data_type() ;
      result = result && ( ( name() == "double" && k0 == Int     ) ||
                           ( name() == "int"    && k0 == Double  ) ) ;
   }
   return result ;
}

//----------------------------------------------------------------------
double
MAC_ConvertTypeExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConvertTypeExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE(ct) ) ;

   double result = (double) P0->to_int( ct ) ;

   return result ;
}

//----------------------------------------------------------------------
int
MAC_ConvertTypeExp:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ConvertTypeExp:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE(ct) ) ;

   int result = (int) P0->to_double(ct) ;

   return result ;
}
