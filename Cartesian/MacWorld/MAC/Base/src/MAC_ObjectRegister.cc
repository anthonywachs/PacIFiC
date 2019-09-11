#include <MAC_ObjectRegister.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Map.hh>
#include <MAC_MapIterator.hh>
#include <MAC_String.hh>

#include <stringVector.hh>

#include <sstream>

struct MAC_ObjectRegister_ERROR
{
   static void n0( std::string const& register_name,
                   std::string const& item_name ) ;
   static void n1( std::string const& register_name,
                   std::string const& item_name,
                   MAC_Map const* items ) ;
   static void n3( std::string const& register_name,
                   std::string const& item_name ) ;
   static void n4( std::string const& register_name ) ;
} ;

//----------------------------------------------------------------------
MAC_ObjectRegister*
MAC_ObjectRegister:: create( MAC_Object* a_owner,
                             std::string const& a_register_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectRegister:: create" ) ;
   MAC_CHECK_PRE( !a_register_name.empty() ) ;

   MAC_ObjectRegister* result =
                    new MAC_ObjectRegister( a_owner, a_register_name ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->register_name() == a_register_name ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ObjectRegister:: MAC_ObjectRegister( MAC_Object* a_owner,
                                         std::string const& a_register_name )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , REGISTER( MAC_Map::create( this ) )
   , NAME( a_register_name )
{
}

//----------------------------------------------------------------------
MAC_ObjectRegister:: ~MAC_ObjectRegister( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
std::string const&
MAC_ObjectRegister:: register_name( void ) const
//----------------------------------------------------------------------
{
   return( NAME ) ;
}

//----------------------------------------------------------------------
bool
MAC_ObjectRegister:: has( std::string const& a_name ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectRegister:: has" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;

   MAC_String const* nn = MAC_String::create( 0, a_name ) ;
   bool result = REGISTER->has_key( nn ) ;
   nn->destroy() ; nn = 0 ;

   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Object*
MAC_ObjectRegister:: item( std::string const& a_name ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectRegister:: item" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;

   MAC_String const* nn = MAC_String::create( 0, a_name ) ;
   MAC_Object* result = REGISTER->item_at( nn ) ;
   nn->destroy() ; nn = 0 ;
   if( result == 0 )
   {
      MAC_ObjectRegister_ERROR::n1( register_name(), a_name, REGISTER ) ;
   }

   MAC_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Iterator*
MAC_ObjectRegister:: create_iterator( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectRegister:: create_iterator" ) ;

   MAC_Iterator* result = REGISTER->create_iterator( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_ObjectRegister:: register_item( std::string const& a_name,
                                    MAC_Object* an_item )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectRegister:: register_item" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   MAC_CHECK_PRE( an_item != 0 ) ;

   MAC_String* nn = MAC_String::create( an_item, a_name ) ;
   if( REGISTER->has_key( nn ) )
   {
      MAC_ObjectRegister_ERROR::n0( register_name(), a_name ) ;
   }
   REGISTER->set_item_at( nn, an_item ) ;

   MAC_CHECK_POST( has( a_name ) ) ;
   MAC_CHECK_POST( item( a_name ) == an_item ) ;
}

//----------------------------------------------------------------------
void
MAC_ObjectRegister:: unregister_item( std::string const& a_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectRegister:: unregister_item" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;

   if( !has( a_name ) )
   {
      MAC_ObjectRegister_ERROR::n3( register_name(), a_name ) ;
   }
   
   MAC_String const* nn = MAC_String::create( 0, a_name ) ;
   REGISTER->remove_at( nn ) ;
   nn->destroy() ; nn = 0 ;

   MAC_CHECK_POST( !has( a_name ) ) ;
}

//----------------------------------------------------------------------
void
MAC_ObjectRegister:: unregister_item( MAC_Object* an_item )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectRegister:: unregister_item" ) ;
   MAC_CHECK_PRE( an_item != 0 ) ;

   MAC_Object const* nn = 0 ;
   MAC_MapIterator* it = REGISTER->create_iterator( 0 ) ;
   for( it->start() ; nn == 0 && it->is_valid() ; it->go_next() )
   {
      if( it->item() == an_item ) nn = it->key() ;
   }
   it->destroy() ; it = 0 ;
   
   if( nn == 0 )
   {
      MAC_ObjectRegister_ERROR::n4( register_name() ) ;
   }
   REGISTER->remove_at( nn ) ;
}

//internal--------------------------------------------------------------
void 
MAC_ObjectRegister_ERROR:: n0( std::string const& register_name,
                               std::string const& item_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "Attempt to register : "   << register_name << std::endl ;
   msg << "           of  name : \"" << item_name << "\"" << std::endl ;
   msg << "which has already been registered" << std::endl ;
   MAC_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
MAC_ObjectRegister_ERROR:: n1( std::string const& register_name,
                               std::string const& item_name ,
                               MAC_Map const* items )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "Request for non registered : "   << register_name << std::endl ;
   msg << "                  of  name : \"" << item_name << "\""
       << std::endl << std::endl ;
   if( items->count() == 0 )
   {
      msg << "No " + register_name + " registered" << std::endl ;
   }
   else
   {  
      msg << items->count() << " registered " + register_name + "(s) of name : " 
          << std::endl ;
      stringVector ref_table( items->count() ) ;
      MAC_MapIterator* it = items->create_iterator( 0 ) ;
      for( size_t i=0 ; it->is_valid() ; it->go_next(), ++i )
      {
         MAC_String const* ss = dynamic_cast<MAC_String const*>( it->key() ) ;
         ref_table(i) = ss->to_string() ;
      }
      it->destroy() ; it = 0 ;
      ref_table.sort() ;
      for( size_t i=0 ; i<ref_table.size() ; ++i )
      {
         msg << "   - \"" << ref_table(i) << "\"" << std::endl ;
      }
   }
   MAC_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
MAC_ObjectRegister_ERROR:: n3( std::string const& register_name,
                               std::string const& item_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "Attempt to unregister : "   << register_name << std::endl ;
   msg << "             of  name : \"" << item_name << "\"" << std::endl ;
   msg << "which is not registered" << std::endl ;
   MAC_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
MAC_ObjectRegister_ERROR:: n4( std::string const& register_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "Attempt to unregister : "   << register_name << std::endl ;
   msg << "a no registered object" << std::endl ;
   MAC_Error::object()->raise_plain( msg.str() ) ;
}
