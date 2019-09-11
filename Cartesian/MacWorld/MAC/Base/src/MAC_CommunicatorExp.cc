#include <MAC_CommunicatorExp.hh>

#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Sequence.hh>
#include <MAC_assertions.hh>


MAC_CommunicatorExp const*
MAC_CommunicatorExp::PROTOTYPE_RANK = new MAC_CommunicatorExp( "rank" ) ;

MAC_CommunicatorExp const*
MAC_CommunicatorExp::PROTOTYPE_NB_RANKS = new MAC_CommunicatorExp( "nb_ranks" ) ;

//----------------------------------------------------------------------
MAC_CommunicatorExp:: MAC_CommunicatorExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , COM( 0 )
{
   MAC_LABEL( "MAC_CommunicatorExp:: MAC_CommunicatorExp" ) ;
}

//----------------------------------------------------------------------
MAC_CommunicatorExp:: MAC_CommunicatorExp(
                               MAC_Object* a_owner,
                               std::string const& a_name,
                               MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , COM( MAC_Exec::communicator() )
{
   MAC_LABEL( "MAC_CommunicatorExp:: MAC_CommunicatorExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_CommunicatorExp:: ~MAC_CommunicatorExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CommunicatorExp:: ~MAC_CommunicatorExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_CommunicatorExp*
MAC_CommunicatorExp:: create_replica(
                                MAC_Object* a_owner,
                                MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CommunicatorExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_CommunicatorExp* result =
             new MAC_CommunicatorExp( a_owner, name(), argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_CommunicatorExp:: valid_arguments(
                              MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CommunicatorExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = some_arguments->count()==0 ;

   return result ;
}

//----------------------------------------------------------------------
std::string const&
MAC_CommunicatorExp:: usage( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CommunicatorExp:: usage" ) ;

   static std::string result ;
   result = name()+"()" ;
   return result ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_CommunicatorExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CommunicatorExp:: data_type" ) ;

   return Int ;
}

//----------------------------------------------------------------------
int
MAC_CommunicatorExp:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CommunicatorExp:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE( ct ) ) ;
   
   int result ;
   if( name()=="rank" )
   {
      result = COM->rank() ;
   }
   else
   {
      result = COM->nb_ranks() ;
   }
   return result ;
}
