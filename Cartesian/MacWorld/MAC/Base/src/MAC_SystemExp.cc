#include <MAC_SystemExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Module.hh>
#include <MAC_Sequence.hh>
#include <MAC_System.hh>

#include <iostream>
#include <sstream>


MAC_SystemExp const*
MAC_SystemExp::PROTOTYPE_PWD = new MAC_SystemExp( "getcwd" ) ;

MAC_SystemExp const* 
MAC_SystemExp::PROTOTYPE_GETENV = new MAC_SystemExp( "getenv" ) ;

MAC_SystemExp const* 
MAC_SystemExp::PROTOTYPE_JOIN = new MAC_SystemExp( "join" ) ;

MAC_SystemExp const* 
MAC_SystemExp::PROTOTYPE_SEPARATOR = new MAC_SystemExp( "path_name_separator" ) ;

MAC_SystemExp const* 
MAC_SystemExp::PROTOTYPE_DIRNAME = new MAC_SystemExp( "dirname" ) ;

MAC_SystemExp const* 
MAC_SystemExp::PROTOTYPE_BASENAME = new MAC_SystemExp( "basename" ) ;

MAC_SystemExp const* 
MAC_SystemExp::PROTOTYPE_GETPID = new MAC_SystemExp( "getpid" ) ;

MAC_SystemExp const* 
MAC_SystemExp::PROTOTYPE_UNAME = new MAC_SystemExp( "uname" ) ;

MAC_SystemExp const* 
MAC_SystemExp::PROTOTYPE_HOSTNAME = new MAC_SystemExp( "host_name" ) ;

//----------------------------------------------------------------------
MAC_SystemExp:: MAC_SystemExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
{
   MAC_LABEL( "MAC_SystemExp:: MAC_SystemExp" ) ;
   MAC_CHECK( a_name=="getcwd" || a_name=="getenv"
              || a_name=="join" || a_name=="getpid"
              || a_name=="uname" || a_name=="host_name"
              || a_name=="path_name_separator"
              || a_name=="dirname" || a_name=="basename" ) ;
}

//----------------------------------------------------------------------
MAC_SystemExp:: MAC_SystemExp( MAC_Object* a_owner,
                               std::string const& a_name,
                               MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
{
   MAC_LABEL( "MAC_SystemExp:: MAC_SystemExp" ) ;
   MAC_CHECK( a_name=="getcwd" || a_name=="getenv"
              || a_name=="join" || a_name=="getpid"
              || a_name=="uname" || a_name=="host_name"
              || a_name=="path_name_separator"
              || a_name=="dirname" || a_name=="basename"  ) ;   
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_SystemExp:: ~MAC_SystemExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SystemExp:: ~MAC_SystemExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_SystemExp*
MAC_SystemExp:: create_replica( MAC_Object* a_owner,
                                MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SystemExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_SystemExp* result = new MAC_SystemExp( a_owner, name(), argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_SystemExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SystemExp:: data_type" ) ;
   MAC_Data::Type result = String ;
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_SystemExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SystemExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = false ;
   if( name()=="getcwd" || name()=="getpid" || name()=="uname" 
       || name()=="host_name" || name()=="path_name_separator" )
   {
      result = some_arguments->count()==0 ;
   }
   else if( name()=="getenv" || name()=="dirname" || name()=="basename" )
   {
      result = some_arguments->count()==1 &&
         extract_arg( some_arguments, 0 )->data_type() == String ;
   }
   else 
   {
      MAC_ASSERT( name()=="join" ) ;
      result = some_arguments->count()>1 ;
      for( size_t i=0 ; i<some_arguments->index_limit() ; i++ )
      {
         result = result &&
            extract_arg( some_arguments, i )->data_type() == String ;
      }
   }

   return result ;
}

//----------------------------------------------------------------------
std::string const&
MAC_SystemExp:: usage( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SystemExp:: usage" ) ;
   static std::string result ;
   if( name()=="getcwd" || name()=="getpid" || name()=="path_name_separator" ||
       name()=="uname" || name()=="host_name" )
   {
      result = name()+"()" ;
   }
   else if( name()=="getenv" || name()=="dirname" || name()=="basename" )
   {
      result = name() + "(SS)" ;
   }
   else if( name()=="join" )
   {
      result = "join(<list of SS>)" ;
   }
   return result ;
}

//----------------------------------------------------------------------
std::string
const&
MAC_SystemExp:: to_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SystemExp:: to_string" ) ;
   MAC_CHECK_PRE( to_string_PRE(ct) ) ;
   
   static char sep = MAC_System::path_name_separator() ;
   
   RESULT_STR = "" ;
  
   if( name()=="getcwd" )
   {
      RESULT_STR = MAC_System::working_directory() ;
   }
   else if( name()=="getpid" )
   {
      std::ostringstream os ;
      os << MAC_System::process_id() ;
      RESULT_STR = os.str() ;
   }
   else if( name()=="getenv" )
   {
      RESULT_STR = MAC_System::getenv( arg(0)->to_string( ct ) ) ;
   }
   else if( name()=="dirname" )
   {
      RESULT_STR = MAC_System::dirname( arg(0)->to_string( ct ) ) ;
   }
   else if( name()=="basename" )
   {
      RESULT_STR = MAC_System::basename( arg(0)->to_string( ct ) ) ;
   }
   else if( name()=="path_name_separator" )
   {
      RESULT_STR = MAC_System::path_name_separator() ;
   }
   else if( name()=="join" )
   {
      RESULT_STR = "" ;
      for( size_t i=0 ; i<nb_arguments() ; i++ )
      {
         std::string const& item = arg(i)->to_string( ct ) ;
         if( i!=0 && item.length()>0 && item[0]!=sep ) RESULT_STR += sep ;
         RESULT_STR += item ;
      }
   }
   else if( name()=="uname" )
   {
      RESULT_STR = MAC_System::sysname() ;
   }
   else if( name()=="host_name" )
   {
      RESULT_STR = MAC_System::host_name() ;
   }
   
   return RESULT_STR ;   
}
