#include <MAC_ObjectTest.hh>

#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_Root.hh>
#include <MAC_Timer.hh>
#include <MAC_assertions.hh>
#include <MAC.hh>

#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>

using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

bool MAC_ObjectTest::FAILURE = false ;

//-------------------------------------------------------------------------
MAC_ObjectTest*
MAC_ObjectTest:: object( std::string const& a_name )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectTest:: object" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;

   MAC_ObjectTest* result =
      static_cast<MAC_ObjectTest*>( plugins_map()->item( a_name ) ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->is_under_ownership_of( MAC_Root::object() ) ) ;
   MAC_CHECK_POST( result->registration_name() == a_name ) ; 
   MAC_CHECK_POST( !result->has_data_deck_explorer() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
MAC_ObjectTest:: ~MAC_ObjectTest( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
MAC_ObjectTest:: MAC_ObjectTest( std::string const& tested_class_name,
                                 std::string const& my_name )
//-------------------------------------------------------------------------
   : MAC_Object( plugins_map() )
   , NAME( my_name )
   , TESTED_CLASS( tested_class_name )
   , EXP( 0 )
   , NB_ELEMENTARY_TESTS( 0 )
   , NB_ELEMENTARY_TESTS_OK( 0 )
{
   MAC_LABEL( "MAC_ObjectTest:: MAC_ObjectTest" ) ;
   
   plugins_map()->register_item( my_name, this ) ;

   MAC_CHECK_POST( owner() == MAC_Root::object() ) ;
   MAC_CHECK_POST( registration_name() == my_name ) ;
   MAC_CHECK_POST( !has_data_deck_explorer() ) ;
}

//-------------------------------------------------------------------------
std::string const&
MAC_ObjectTest:: registration_name( void ) const
//-------------------------------------------------------------------------
{
   return( NAME ) ;
}

//-------------------------------------------------------------------------
void
MAC_ObjectTest:: set_data_deck_explorer( MAC_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectTest:: set_data_deck_explorer" ) ;
   MAC_CHECK_PRE( !has_data_deck_explorer() ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   MAC_CHECK_PRE( exp->string_data("concrete_name") == registration_name() ) ;

   EXP = exp->create_clone( this ) ;
   
   MAC_CHECK_POST( has_data_deck_explorer() ) ;
}

//-------------------------------------------------------------------------
bool
MAC_ObjectTest:: has_data_deck_explorer( void ) const
//-------------------------------------------------------------------------
{
   return( EXP!=0 ) ;
}

//-------------------------------------------------------------------------
MAC_ModuleExplorer const*
MAC_ObjectTest:: data_deck_explorer( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectTest:: data_deck_explorer" ) ;
   MAC_CHECK_PRE( has_data_deck_explorer() ) ;

   MAC_ModuleExplorer const* result = EXP ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
MAC_ObjectTest:: run_all_tests( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectTest:: run_all_tests" ) ;

   MAC_Timer* timer = MAC_Timer::create( 0 ) ;

   out() << endl ;
   out() << "----------------------------------------------------"    << endl ;
   out() << "|  Unit tests performed on class " << TESTED_CLASS << " :"
         << endl ;
   out() << "===================================================="    << endl ;
      
   reset_all_tests() ;
      
   timer->reset() ;
   timer->start() ;

   process_all_tests() ;

   timer->stop() ;
   out() << "===================================================="    << endl ;

   if( NB_ELEMENTARY_TESTS != NB_ELEMENTARY_TESTS_OK )
   {
      out() << "|  Class test FAILED !!!!!!!!!!!!!" << endl ;
   }
   out() << "|  End of " << NB_ELEMENTARY_TESTS << " test" ;
   if( NB_ELEMENTARY_TESTS> 1 )
   {
      out() << "s" ;
   }
   out() << " of class " <<  TESTED_CLASS 
         << " in " << timer->time() << " s " << endl ;
   out() << "----------------------------------------------------"    << endl ;

   timer->destroy() ;
}

//-------------------------------------------------------------------------
bool
MAC_ObjectTest:: tests_of_all_instances_are_successful( void )
//-------------------------------------------------------------------------
{
   return( !FAILURE ) ;
}

//-------------------------------------------------------------------------
void
MAC_ObjectTest:: reset_all_tests( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
MAC_ObjectTest:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectTest:: process_all_tests" ) ;

   if( EXP == 0 )
   {
      std::ostringstream m ;
      m << registration_name() << endl << endl ;
      m << "   The default implementation of \"process_all_tests\""  << endl ;
      m << "   provided by \"MAC_ObjectTest\" assumes that there is" << endl ;
      m << "   a specific data deck for \"" << registration_name() << "\"."
	<< endl << endl ;
      m << "   If not, the member function \"process_all_tests\""    << endl ;
      m << "   should be overridden." ;
      MAC_Error::object()->raise_plain( m.str() ) ;
   }
   
   EXP->start_module_iterator() ;
   for( ; EXP->is_valid_module() ; EXP->go_next_module() )
   {
      MAC_ModuleExplorer* ee = EXP->create_subexplorer( 0 ) ;
      process_one_test( ee ) ;
      ee->destroy() ;
   }
}

//-------------------------------------------------------------------------
void
MAC_ObjectTest:: process_one_test( MAC_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectTest:: process_one_test" ) ;

   std::ostringstream m ;
   m << registration_name() << endl << endl ;
   m << "   member function \"process_one_test\" is not implemented"  << endl ;
   MAC_Error::object()->raise_plain( m.str() ) ;
   
}

//-------------------------------------------------------------------------
void
MAC_ObjectTest:: notify_one_test_result( std::string const& displayed_name, 
                                         bool success )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectTest:: notify_one_test_result" ) ;
   MAC_CHECK_COLLECTIVE( true ) ;
   
   success = MAC_Exec::communicator()->boolean_and(success) ;
   
   out() << "| ... " << std::setw(40) << displayed_name << " :" ;
   if( success )
   {
      out() << "  OK" ;
      NB_ELEMENTARY_TESTS_OK++ ;
   }
   else
   {
      out() << " FAIL" ;
      FAILURE = true ;
   }
   out()  << endl ;
   NB_ELEMENTARY_TESTS++ ;
}

//----------------------------------------------------------------------------
void
MAC_ObjectTest:: print_time_result( std::string const& name, 
                                    double tt ) const
//----------------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = MAC::out().flags() ;
   MAC::out().setf( ios_base::uppercase | ios_base::scientific ) ;

   MAC::out() << "| ... " << setw( 50 ) << name << " CPU "
              << setprecision( 6 ) << setw( 15 ) << tt << endl ;

   MAC::out().flags( original_flags ) ;
}

//----------------------------------------------------------------------------
void
MAC_ObjectTest:: print_memory_result( std::string const& name, 
                                      size_t memory_size ) const
//----------------------------------------------------------------------------
{
   MAC::out() << "| ... " << setw( 50 ) << name << " MEM "
              << setw( 15 ) << memory_size << endl ;
}

//-------------------------------------------------------------------------
std::ostream &
MAC_ObjectTest:: out( void )
//-------------------------------------------------------------------------
{
   return( MAC::out() ) ;
}

//----------------------------------------------------------------------
MAC_ObjectRegister*
MAC_ObjectTest:: plugins_map( void )
//----------------------------------------------------------------------
{
   static MAC_ObjectRegister* result =
                MAC_ObjectRegister::create( MAC_Root::object(),
                                            "MAC_ObjectTest descendant" ) ;
   return( result ) ;
}

