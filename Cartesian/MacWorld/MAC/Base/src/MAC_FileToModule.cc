#include <MAC_FileToModule.hh>

#include <MAC_Error.hh>
#include <MAC_Iterator.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_Root.hh>
#include <MAC_assertions.hh>

#include <sstream>

using std::endl ;

struct MAC_FileToModule_ERROR
{
   static void n0( std::string const& a_format,
                   std::string const& a_default_motif,
                   MAC_FileToModule const* other ) ;
} ;

//----------------------------------------------------------------------------
MAC_FileToModule const*
MAC_FileToModule:: object( std::string const& format )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_FileToModule:: object" ) ;

   MAC_FileToModule const* result =
      static_cast<MAC_FileToModule const*>(
                                    plugins_map()->item( format ) ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->is_under_ownership_of( MAC_Root::object() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_FileToModule:: MAC_FileToModule( std::string const& a_format,
                                     std::string const& a_default_motif )
//----------------------------------------------------------------------------
   : MAC_Object( plugins_map() )
   , MY_FORMAT( a_format )
   , MY_MOTIF( a_default_motif )
{
   MAC_LABEL( "MAC_FileToModule:: MAC_FileToModule" ) ;
   
   MAC_Iterator* it = plugins_map()->create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      MAC_FileToModule const* oo = 
                       static_cast< MAC_FileToModule const* >( it->item() ) ;
      if( oo->default_motif() == a_default_motif )
         MAC_FileToModule_ERROR::n0( a_format, a_default_motif, oo ) ;
   }
   it->destroy() ; it = 0 ;

   plugins_map()->register_item( a_format, this ) ;
   if( formats().empty() )
   {
      formats() = a_format ;
   }
   else
   {
      formats() += "," + a_format ;
   }
}

//----------------------------------------------------------------------------
MAC_FileToModule:: ~MAC_FileToModule( void  )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
bool
MAC_FileToModule:: has( std::string const& format )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_FileToModule:: has" ) ;
   
   bool result = plugins_map()->has( format ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
MAC_FileToModule:: find_file_format( std::string const& a_file_name,
                                     std::string& a_format )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_FileToModule:: find_file_format" ) ;
   
   MAC_Iterator* it = plugins_map()->create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      MAC_FileToModule const* oo = 
                       static_cast< MAC_FileToModule const* >( it->item() ) ;
      if( a_file_name.find( oo->default_motif() ) < a_file_name.length() )
      {
         a_format = oo->format() ;
         break ;
      }
   }
   it->destroy() ; it = 0 ;
}

//----------------------------------------------------------------------------
std::string const&
MAC_FileToModule:: list_of_formats( void )
//----------------------------------------------------------------------------
{
   return( formats() ) ; 
}

//----------------------------------------------------------------------------
std::string const&
MAC_FileToModule:: format( void ) const
//----------------------------------------------------------------------------
{
   return( MY_FORMAT ) ;
}

//----------------------------------------------------------------------------
std::string const&
MAC_FileToModule:: default_motif( void ) const
//----------------------------------------------------------------------------
{
   return( MY_MOTIF ) ;
}

//----------------------------------------------------------------------------
MAC_ObjectRegister*
MAC_FileToModule:: plugins_map( void )
//----------------------------------------------------------------------------
{
   static MAC_ObjectRegister* result =
      MAC_ObjectRegister::create( MAC_Root::object(),
                                  "MAC_FileToModule descendant" ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
std::string&
MAC_FileToModule:: formats( void )
//----------------------------------------------------------------------------
{
   static std::string result ;
   return( result ) ;
}

//----------------------------------------------------------------------------
bool
MAC_FileToModule:: create_from_file_PRE( MAC_Object* a_owner,
                                         std::string const& module_name,
                                         std::string const& file_name ) const
//----------------------------------------------------------------------------
{
   MAC_ASSERT( !module_name.empty() ) ;
   MAC_ASSERT( !file_name.empty() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------------
bool
MAC_FileToModule:: create_from_file_POST( MAC_Module const* result,
                                          MAC_Object* a_owner,
                                          std::string const& module_name,
                                          std::string const& file_name ) const
//----------------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   MAC_ASSERT( result->name() == module_name ) ;
   return( true ) ;
}

//internal--------------------------------------------------------------------
void 
MAC_FileToModule_ERROR:: n0( std::string const& a_format,
                             std::string const& a_default_motif,
                             MAC_FileToModule const* other )
//internal--------------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** MAC_FileToModule error:" << endl << endl ;
   mesg << "    Attempt to register the format: \"" << a_format 
        << "\"" << endl ;
   mesg << "    associated to the filename motif: \"" 
        << a_default_motif << "\"" << endl << endl ;
   mesg << "    The filename motif \""
        << a_default_motif << "\" is already handled" << endl ;
   mesg << "    by the format: \"" << other->format() << "\"" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}

