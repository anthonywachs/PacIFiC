#include <MAC_DoubleComparator.hh>

#include <MAC_assertions.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_Root.hh>

//----------------------------------------------------------------------
MAC_DoubleComparator const*
MAC_DoubleComparator:: make( MAC_Object* a_owner,
                             MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleComparator:: make" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   std::string name = exp->string_data( "concrete_name" ) ;
   MAC_DoubleComparator const* proto =
      static_cast<MAC_DoubleComparator const*>(
                                    plugins_map()->item( name ) ) ;
   MAC_ASSERT( proto->is_a_prototype() ) ;
      
   MAC_DoubleComparator const* result = proto->create_replica( a_owner, exp ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
MAC_DoubleComparator:: MAC_DoubleComparator( MAC_Object* a_owner )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , IS_PROTO( false )
{
   MAC_CHECK_POST( owner() == a_owner ) ;
   MAC_CHECK_POST( !is_a_prototype() ) ;
}

//----------------------------------------------------------------------
MAC_DoubleComparator:: MAC_DoubleComparator( std::string const& name )
//----------------------------------------------------------------------
   : MAC_Object( plugins_map() )
   , IS_PROTO( true )
{
   MAC_LABEL( "MAC_DoubleComparator:: MAC_DoubleComparator" ) ;
   
   plugins_map()->register_item( name, this ) ;
   
   MAC_CHECK_POST( is_under_ownership_of( plugins_map() ) ) ;
   MAC_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
MAC_DoubleComparator:: ~MAC_DoubleComparator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
MAC_DoubleComparator:: create_replica_PRE(
                                   MAC_Object* a_owner,
                                   MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_a_prototype() ) ;
   MAC_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_DoubleComparator:: create_replica_POST(
                                   MAC_DoubleComparator const* result,
                                   MAC_Object* a_owner,
                                   MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   MAC_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}
//----------------------------------------------------------------------
bool
MAC_DoubleComparator:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_DoubleComparator:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------
MAC_ObjectRegister*
MAC_DoubleComparator:: plugins_map( void )
//----------------------------------------------------------------------
{
   static MAC_ObjectRegister* result =
            MAC_ObjectRegister::create( MAC_Root::object(),
                                        "MAC_DoubleComparator descendant" ) ;
   return( result ) ;
}
