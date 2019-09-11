#include <MAC_ExternalAPI.hh>

#include <MAC_Error.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_assertions.hh>

#include <stringVector.hh>

#include <iostream>

//----------------------------------------------------------------------
MAC_ExternalAPI:: MAC_ExternalAPI( std::string const& a_name,
                                   size_t a_priority_level )
//----------------------------------------------------------------------
   : MAC_Object( 0 )
   , MY_NAME( a_name )
   , MY_PRIORITY( a_priority_level )
{
   MAC_LABEL( "MAC_ExternalAPI:: MAC_ExternalAPI" ) ;

   if( a_priority_level>9 )
   {
      MAC_Error::object()->raise_plain(
         "MAC_ExternalAPI error :\n"
         "   external API of name : \n"+a_name+"\"\n"
         "   is defined with a priority level greater than 9" ) ;
   }
   
   plugins_map()->register_item( a_name, this ) ;
   plugins_names().append( a_name ) ;
}




//----------------------------------------------------------------------
MAC_ExternalAPI:: ~MAC_ExternalAPI( void )
//----------------------------------------------------------------------
{
}




//----------------------------------------------------------------------
void
MAC_ExternalAPI:: initialize_all_APIs( int& argc, char**& argv )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExternalAPI:: initialize_all_APIs" ) ;

   size_t const nb_APIs = plugins_names().size() ;
   for( size_t level=9 ; level<10 ; --level )
   {
      for( size_t i=0 ; i<nb_APIs ; ++i )
      {
         std::string const& name = plugins_names()(i) ;
         MAC_ExternalAPI* api =
            static_cast<MAC_ExternalAPI*>( plugins_map()->item( name ) ) ;
         if( api->MY_PRIORITY==level )
         {
            api->initialize( argc, argv ) ;
         }
      }
   }
}




//----------------------------------------------------------------------
void
MAC_ExternalAPI:: terminate_all_APIs( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExternalAPI:: terminate_all_APIs" ) ;
   
   size_t const nb_APIs = plugins_names().size() ;
   for( size_t level=0 ; level<10 ; ++level )
   {
      for( size_t i=nb_APIs-1 ; i<nb_APIs ; --i )
      {
         std::string const& name = plugins_names()(i) ;
         if( plugins_map()->has( name ) )
         {
            MAC_ExternalAPI* api =
               static_cast<MAC_ExternalAPI*>( plugins_map()->item( name ) ) ;
            if( api->MY_PRIORITY==level )
            {
               plugins_map()->unregister_item( name ) ;
               api->destroy() ;
            }
         }
      }
   }
   plugins_map()->destroy() ;
}




//----------------------------------------------------------------------
MAC_ObjectRegister*
MAC_ExternalAPI:: plugins_map( void )
//----------------------------------------------------------------------
{
   static MAC_ObjectRegister* result =
         MAC_ObjectRegister::create( 0, "MAC_ExternalAPI descendant" ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
stringVector&
MAC_ExternalAPI:: plugins_names( void )
//----------------------------------------------------------------------
{
   static stringVector result(0) ;
   return( result ) ;
}
