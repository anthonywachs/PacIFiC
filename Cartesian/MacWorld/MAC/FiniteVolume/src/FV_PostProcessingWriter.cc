#include <FV_PostProcessingWriter.hh>
#include <FV_DomainAndFields.hh>
#include <FV_Mesh.hh>
#include <FV.hh>
#include <FV_DiscreteField.hh>
#include <MAC_Communicator.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_Root.hh>
#include <MAC_Error.hh>
#include <list>
using std::list ;

/* Create and initialize an instance of FV_PostProcessingWriter
-----------------------------------------------------------------*/
FV_PostProcessingWriter* FV_PostProcessingWriter:: make( 
         MAC_Object* a_owner,
	 std::string const& a_name,
	 MAC_ModuleExplorer const* exp,
	 MAC_Communicator const* com,
	 list< FV_DiscreteField const* > a_fields,
	 FV_Mesh const* a_primary_mesh,
	 bool a_binary )
{
   MAC_LABEL( "FV_PostProcessingWriter:: make" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   
   FV_PostProcessingWriter const* proto =
      static_cast<FV_PostProcessingWriter const*>( 
      	plugins_map()->item( a_name ) ) ;
   MAC_ASSERT( proto->is_a_prototype() ) ;

   FV_PostProcessingWriter* result = proto->create_replica( a_owner, 
         exp, com, a_fields, a_primary_mesh, a_binary ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
   
}



/* In the constructor called by `::create_replica': initialization
   of the base class subobject
------------------------------------------------------------------*/
FV_PostProcessingWriter:: FV_PostProcessingWriter( MAC_Object* a_owner )
   : MAC_Object( a_owner )
   , IS_PROTO( false )
{
   MAC_LABEL( "FV_PostProcessingWriter:: FV_PostProcessingWriter" ) ;
   MAC_CHECK_POST( owner() == a_owner ) ;
   MAC_CHECK_POST( !is_a_prototype() ) ;
   
}




/* Registration of an instance
------------------------------*/
FV_PostProcessingWriter:: FV_PostProcessingWriter( 
	std::string const& a_name )
   : MAC_Object( plugins_map() )
   , IS_PROTO( true )
{
   MAC_LABEL( "FV_PostProcessingWriter:: FV_PostProcessingWriter" ) ;
   
   plugins_map()->register_item( a_name, this ) ;

   MAC_CHECK_POST( owner() == plugins_map() ) ;
   MAC_CHECK_POST( is_a_prototype() ) ;
   
}




/* Destructor
-------------*/
FV_PostProcessingWriter:: ~FV_PostProcessingWriter( void )
{}




/* Test invariance
------------------*/
bool FV_PostProcessingWriter:: invariant( void ) const
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   return( true ) ;
   
}




/* Postcondition of create
--------------------------*/
bool FV_PostProcessingWriter:: create_replica_POST( 
	 FV_PostProcessingWriter const* result,
	 MAC_Object* a_owner ) const
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   MAC_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
   
}




/* Return if this is a prototype
--------------------------------*/
bool FV_PostProcessingWriter:: is_a_prototype( void ) const
{
   return( IS_PROTO ) ;
   
}




/* Return a pointer to the object register
------------------------------------------*/
MAC_ObjectRegister* FV_PostProcessingWriter:: plugins_map( void )
{
   static MAC_ObjectRegister* result =
      MAC_ObjectRegister::create( MAC_Root::object(),
                                  "FV_PostProcessingWriter descendant" ) ;
   return( result ) ;
   
}







