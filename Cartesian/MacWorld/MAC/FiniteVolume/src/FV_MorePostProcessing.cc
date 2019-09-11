#include <FV_MorePostProcessing.hh>
#include <FV_StepByStepProgression.hh>
#include <FV.hh>
#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_Root.hh>
#include <MAC_assertions.hh>
#include <iostream>

using std::cout ;
using std::endl ;
using std::string ;

FV_MorePostProcessing const* 
FV_MorePostProcessing:: PROTOTYPE = new FV_MorePostProcessing() ;

//-------------------------------------------------------------------------
FV_MorePostProcessing:: FV_MorePostProcessing( void )
//-------------------------------------------------------------------------
   : MAC_Application( "FV_MorePostProcessing" )
   , READER( 0 )
   , APPLI( 0 )
{
}



   
//---------------------------------------------------------------------------
FV_MorePostProcessing*
FV_MorePostProcessing:: create_replica( 
        MAC_Object* a_owner,
        MAC_ModuleExplorer const* exp,
	double const& initial_time ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_MorePostProcessing:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   FV_MorePostProcessing* result = new FV_MorePostProcessing( a_owner,
	exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
FV_MorePostProcessing:: FV_MorePostProcessing( 
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : MAC_Application( a_owner, exp )
   , READER( 0 )
   , APPLI( 0 )
   , ROOT_EXP( exp )
{
   MAC_ModuleExplorer* ee = 0 ;

   ee = exp->create_subexplorer( 0, "MAC_ObjectReader" ) ;
   READER = MAC_ObjectReader::create( this, ee ) ;
   ee->destroy() ;

   size_t cycle_number = 0 ;
   if( exp->has_entry( "cycle_number" ) )
   {
      cycle_number = exp->int_data( "cycle_number" ) ;
   }
 
   MAC_Module* mod = create_modified_data_deck_module() ;
   ee = MAC_ModuleExplorer::create( 0, mod ) ;
   
   APPLI = FV_StepByStepProgression::create( this, ee ) ;
   ee->destroy() ;

   READER->seek_cycle( cycle_number ) ;
   APPLI->restore_registered_objects( READER ) ;
   READER->close_cycle() ;

   if ( MAC_Exec::communicator()->rank() == 0 )
   {
      FV::out() << endl << "++++++ FV_MorePostProcessing " 
              << " ++++++" << endl << endl ;   
   }

   MAC_CHECK_POST( APPLI!=0 ) ;
}




//---------------------------------------------------------------------------
FV_MorePostProcessing:: ~FV_MorePostProcessing( void )
//---------------------------------------------------------------------------
{
}




//-------------------------------------------------------------------------
void
FV_MorePostProcessing:: run( std::string const& inputRestartFileName ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_MorePostProcessing:: run" ) ;

   APPLI->do_more_post_processing( ROOT_EXP ) ;
}




//---------------------------------------------------------------------------
MAC_Module*
FV_MorePostProcessing:: create_modified_data_deck_module( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_MorePostProcessing:: create_modified_data_deck_module" ) ;

   MAC_Module* m = 0 ;

   MAC_Module* header = READER->header_module() ;
   if( header->has_module( "MAC_Application" ) )
   {
      m = header->module( "MAC_Application" ) ;
   }
   else
   {
      MAC_Error::object()->raise_plain( "invalid restart file" ) ; 
   }
   MAC_ASSERT( m != 0 ) ;
   MAC_Module* result = m->create_clone( this ) ;

   result->remove_entry( "concrete_name" ) ;
   result->remove_module( "MAC_ObjectWriter" ) ;
   result->remove_module( "FV_TimeIterator" ) ;
   result->remove_entry( "graphics_output_times" ) ;
   result->remove_entry( "state_saving_times" ) ;
   if ( result->has_entry( "INITIAL_CYCLE_NUMBER" ) )
     result->remove_entry( "INITIAL_CYCLE_NUMBER" ) ; 
   result->remove_module( "FV_DomainAndFields/FV_ResultSaver" ) ;     
   
   change_owner( MAC_Root::object(), result ) ;
   
   MAC_CHECK( result != 0 ) ;
   MAC_CHECK( result->owner() == MAC_Root::object() ) ;
   return( result ) ;
}
