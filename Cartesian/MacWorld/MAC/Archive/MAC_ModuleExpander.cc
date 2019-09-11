#include <MAC_ModuleExpander.hh>

#include <MAC.hh>
#include <MAC_Context.hh>
#include <MAC_Exec.hh>
#include <MAC_ExtractionExp.hh>
#include <MAC_Error.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_Root.hh>
#include <MAC_String.hh>
#include <MAC_Variable.hh>
#include <MAC_assertions.hh>

#include <fstream>
#include <iostream>


MAC_ModuleExpander const*
MAC_ModuleExpander::PROTOTYPE = new MAC_ModuleExpander() ;

//----------------------------------------------------------------------
MAC_ModuleExpander:: MAC_ModuleExpander( void )
//----------------------------------------------------------------------
   : MAC_Application( "MAC_ModuleExpander" )
   , BASE( "" )
   , OUTPUT( "" )
   , INPUT( "" )
{
}

//----------------------------------------------------------------------
MAC_ModuleExpander* 
MAC_ModuleExpander:: create_replica(
              MAC_Object* a_owner, MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExpander:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   std::string submod = ( exp->has_entry( "submodule_to_expand" ) ?
                          exp->string_data( "submodule_to_expand" ) :
                          "" ) ;
   
   MAC_ModuleExpander* result =
      new MAC_ModuleExpander( a_owner,
                              exp->string_data( "skeleton_file" ),
                              exp->string_data( "input_file" ),
                              exp->string_data( "expanded_file" ),
                              submod ) ;

   MAC_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ModuleExpander* 
MAC_ModuleExpander:: create_replica_from_args(
                         MAC_Object* a_owner, stringVector& args ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExpander:: create_replica_from_args" ) ;
   
   if( args.size() != 3 ) notify_error_in_arguments() ;
   
   MAC_ModuleExpander* result =
      new MAC_ModuleExpander( a_owner, args(0), args(1), args(2), "" ) ;
   args.remove_at(0) ;
   args.remove_at(0) ;
   args.remove_at(0) ;

   MAC_CHECK_POST( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_ModuleExpander:: run( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExpander:: run" ) ;

   // Read input file:
   MAC_Module* mod = MAC_Module::create( 0, "Root", INPUT,
                                         MAC_Exec::execution_context() ) ;
   MAC_ModuleIterator* it = mod->create_module_iterator( mod ) ;
   it->start() ;
   MAC_Module* input_mod = it->item()->create_clone( this ) ;
   if( !SUBMOD.empty() )
   {
      if( !input_mod->has_module( SUBMOD ) )
      {
         MAC_Error::object()->raise_plain(
            "Submodule " + SUBMOD + " is not a submodule name in " + INPUT ) ;
      }
      input_mod = input_mod->module( SUBMOD ) ;
   }
   
   MAC_Module const* expanded_module =
                           create_expanded_module( mod, input_mod, BASE ) ;

   std::ofstream out( OUTPUT.c_str() ) ;
   out.close() ;
   expanded_module->write( OUTPUT, "text" ) ;

   mod->destroy() ; mod = 0 ; it = 0 ; expanded_module = 0 ;
}

//---------------------------------------------------------------------------
void
MAC_ModuleExpander:: print_usage( void ) const
//---------------------------------------------------------------------------
{
   MAC::out() << usage_title( "macsdd" )  ;
   MAC::out() << "<base_dir> <input_file> <expanded_file>"
              << std::endl << std::endl ;
   MAC::out() << "     Create the full data deck associated to a simplified"
              << std::endl
              << "     one and a data base."
              << std::endl ;
}

//---------------------------------------------------------------------------
void
MAC_ModuleExpander:: print_operands( void ) const
//---------------------------------------------------------------------------
{
   MAC::out() << operands_title() ;
   MAC::out() << "     <skeleton_file>" << std::endl
              << "          skeleton data deck"
              << std::endl << std::endl ;
   MAC::out() << "     <input_file>" << std::endl
              << "          simplified data deck"
              << std::endl << std::endl ;
   MAC::out() << "     <expanded_file>" << std::endl
              << "          full data deck created"
              << std::endl << std::endl ;
}

//----------------------------------------------------------------------
MAC_ModuleExpander:: MAC_ModuleExpander( MAC_Object* a_owner,
                                         std::string const& skeleton_file,
                                         std::string const& input_file,
                                         std::string const& expanded_file,
                                         std::string const& submod_name  )
//----------------------------------------------------------------------
   : MAC_Application( a_owner, 0 )
   , BASE( skeleton_file )
   , OUTPUT( expanded_file )
   , INPUT( input_file )
   , SUBMOD( submod_name )
{
   MAC_LABEL( "MAC_ModuleExpander:: MAC_ModuleExpander" ) ;
}

//----------------------------------------------------------------------
MAC_ModuleExpander:: ~MAC_ModuleExpander( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
MAC_Module const*
MAC_ModuleExpander:: create_expanded_module(
                                 MAC_Object* a_owner,
                                 MAC_Module const* input_mod,
                                 std::string const& skeleton_file_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExpander:: create_expanded_module" ) ;
   MAC_CHECK_PRE( input_mod != 0 ) ;
   MAC_CHECK_PRE( ! skeleton_file_name.empty() ) ;

   // Initialize data-base:
   MAC_ExtractionExp::initialize( input_mod ) ;

   // Read skeleton file:
   MAC_Module* result = create_skeleton_module( a_owner, skeleton_file_name ) ;

   MAC_ExtractionExp::reset() ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->is_under_ownership_of( a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Module*
MAC_ModuleExpander:: create_skeleton_module( MAC_Object* a_owner,
                                             std::string const& file_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExpander:: create_skeleton_module" ) ;
   MAC_CHECK( !file_name.empty() ) ;
   MAC_CHECK( MAC_ExtractionExp::is_initialized() ) ;
   
   MAC_Module* mod = MAC_Module::create( a_owner, "Root", file_name,
                                         MAC_Exec::execution_context() ) ;
   MAC_ModuleIterator* it = mod->create_module_iterator( 0 ) ;
   it->start() ;
   MAC_ASSERT( it->is_valid() ) ;
   MAC_Module* result = it->item() ;
   it->go_next() ;
   MAC_ASSERT( !it->is_valid() ) ;
   it->destroy() ; it = 0 ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->is_under_ownership_of( a_owner ) ) ;
   return( result ) ;
}

