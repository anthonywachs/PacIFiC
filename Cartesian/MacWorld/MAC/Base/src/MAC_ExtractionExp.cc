#include <MAC_ExtractionExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Bool.hh>
#include <MAC_Communicator.hh>
#include <MAC_Context.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_List.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Sequence.hh>
#include <MAC_String.hh>
#include <MAC_System.hh>
#include <MAC_Variable.hh>
#include <stringVector.hh>

#include <fstream>
#include <iostream>
#include <sstream>

MAC_Module const* MAC_ExtractionExp::DB_MOD = 0 ;
MAC_ExtractionExp* MAC_ExtractionExp::PROTO_HAS_DATA =
   new MAC_ExtractionExp( "has_data", MAC_ExtractionExp::has_data ) ;
MAC_ExtractionExp* MAC_ExtractionExp::PROTO_DATA =
   new MAC_ExtractionExp( "extracted_data", MAC_ExtractionExp::ext_data ) ;
MAC_ExtractionExp* MAC_ExtractionExp::PROTO_HAS_MOD =
   new MAC_ExtractionExp( "has_module", MAC_ExtractionExp::has_mod ) ;
MAC_ExtractionExp* MAC_ExtractionExp::PROTO_EXTRACTED_MODULE =
   new MAC_ExtractionExp( "extracted_module", MAC_ExtractionExp::ext_mod ) ;

//----------------------------------------------------------------------
MAC_ExtractionExp:: MAC_ExtractionExp( std::string const& a_name,
                                       ExtractionExp op ) 
//----------------------------------------------------------------------
   : MAC_TransferExp( a_name )
   , OP( op )
   , TEMP_FILE_NAME( "" )
   , SRC( 0 )
{
   MAC_LABEL( "MAC_ExtractionExp:: MAC_ExtractionExp" ) ;
}

//----------------------------------------------------------------------
MAC_ExtractionExp:: MAC_ExtractionExp( MAC_Object* a_owner,
                                       ExtractionExp op,
                                       std::string const& a_name,
                                       MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_TransferExp( a_owner, a_name, argument_list )
   , OP( op )
   , TEMP_FILE_NAME( "" )
   , SRC( 0 )
{
   MAC_LABEL( "MAC_ExtractionExp:: MAC_ExtractionExp" ) ;
   MAC_ASSERT( is_initialized() ) ;
   
   // Data name:
   std::string const& n = data_name( name(), arg(0) ) ;

   if( OP == ext_data )
   {
      MAC_ASSERT( name() == "extracted_data" ) ;
      if( DB_MOD->has_entry( n ) )
      {
         MAC_Data const* a_data = DB_MOD->data_of_entry( n ) ;
         MAC_List* l = MAC_List::create( 0 ) ;
         a_data->declare( l ) ;
         if( l->count() == 0 )
         {
            SRC = a_data ;
         }
         else
         {
            std::string const& dirname = MAC_Module::dirname( n )  ;
            MAC_Module const* m =
               ( dirname.empty() ? DB_MOD : DB_MOD->module( dirname ) ) ;
            MAC_Data* d =
                 MAC_DataWithContext::create( 0, a_data, m->context() ) ;
            SRC = d->create_simplification( this ) ;
            d->destroy() ; d = 0 ;
         }
         l->destroy() ; l = 0 ;
      }
      else if( argument_list->count() == 2 )
      {
         SRC = arg(1) ;
      }
      else
      {
         raise_error( "    missing entry: "+n ) ;
      }
   }
   else if( OP == has_data ) 
   {
      MAC_ASSERT( name() == "has_data" ) ;
      SRC = MAC_Bool::create( this, DB_MOD->has_entry( n ) ) ;
   }
   else if( OP == has_mod ) 
   {
      MAC_ASSERT( name() == "has_module" ) ;
      SRC = MAC_Bool::create( this, DB_MOD->has_module( n ) ) ;
   }
   else if( OP == ext_mod ) 
   {
      MAC_ASSERT( name() == "extracted_module" ) ;
      if( DB_MOD->has_module( n ) )
      {
         TEMP_FILE_NAME = temporary_file() ;
         extract_module( TEMP_FILE_NAME, n, data_name( name(), arg(1) ) ) ;
         SRC = MAC_String::create( this, TEMP_FILE_NAME ) ;
      }
      else if( argument_list->count() == 3 )
      {
         SRC = arg(2) ;
      }
      else
      {
         raise_error( "    missing module: "+n ) ;
      }
   }
   
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_ExtractionExp:: ~MAC_ExtractionExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: ~MAC_ExtractionExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( OP == ext_mod )
   {
      MAC_System::erase( TEMP_FILE_NAME ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_ExtractionExp:: initialize( MAC_Module const* mod )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: initialize" ) ;
   MAC_CHECK_PRE( mod != 0 ) ;
   
   DB_MOD = mod ;
   
   MAC_CHECK_POST( is_initialized() ) ;
   MAC_CHECK_POST( data_base() == mod ) ;
}

//----------------------------------------------------------------------
void
MAC_ExtractionExp:: reset( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: reset" ) ;
   
   DB_MOD = 0 ;
   
   MAC_CHECK_POST( !is_initialized() ) ;
}

//----------------------------------------------------------------------
bool
MAC_ExtractionExp:: is_initialized( void )
//----------------------------------------------------------------------
{
   return( DB_MOD != 0 ) ;
}

//----------------------------------------------------------------------
MAC_Module const*
MAC_ExtractionExp:: data_base( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: data_base" ) ;
   MAC_CHECK_PRE( is_initialized() ) ;
   MAC_Module const* result = DB_MOD ;
   MAC_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ExtractionExp*
MAC_ExtractionExp:: create_replica( MAC_Object* a_owner, 
                                    MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   if( !is_initialized() )
   {
      raise_error(
         "*** MAC_ExtractionExp: error\n"
         "    Can't evaluate expression before initialize self" ) ;
   }
   
   MAC_ExtractionExp* result =
           new MAC_ExtractionExp( a_owner, OP, name(), argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_ExtractionExp:: declare( MAC_List* lst ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: declare" ) ;
   MAC_CHECK_PRE( declare_PRE( lst ) ) ;
   
   SRC->declare( lst ) ;
   
   MAC_CHECK_POST( declare_POST( lst ) ) ;
}

//----------------------------------------------------------------------
bool
MAC_ExtractionExp:: context_has_required_variables(
                                           MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: context_has_required_variables" ) ;
   MAC_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;
   
   return( SRC->context_has_required_variables( ct ) ) ;
}

//----------------------------------------------------------------------
stringVector const&
MAC_ExtractionExp:: undefined_variables( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: undefined_variables" ) ;
   
   return( SRC->undefined_variables( ct ) ) ;
}

//----------------------------------------------------------------------
bool
MAC_ExtractionExp:: value_can_be_evaluated( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: value_can_be_evaluated" ) ;
   
   return( SRC->value_can_be_evaluated( ct ) ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_ExtractionExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: data_type" ) ;
   return( SRC->data_type() ) ;
}

//----------------------------------------------------------------------
bool
MAC_ExtractionExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = true ;
   if( OP == has_data || OP == has_mod )
   {
      result = some_arguments->count()==1 &&
               extract_arg( some_arguments, 0 )->data_type() == String ;
   }
   else if( OP == ext_mod )
   {
      result = ( some_arguments->count()==2 || some_arguments->count()==3 ) &&
               extract_arg( some_arguments, 0 )->data_type() == String &&
               extract_arg( some_arguments, 1 )->data_type() == String ;
      if( result && some_arguments->count()==3  )
      {
         result = ( extract_arg( some_arguments, 2 )->data_type() == String ) ;
      }
   }
   else
   {
      result = ( some_arguments->count()==1 || some_arguments->count()==2 ) &&
               extract_arg( some_arguments, 0 )->data_type() == String ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_ExtractionExp:: usage( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: usage" ) ;
   static std::string result ;
   if( OP == has_data )
   {
      result = name() + "(SS)" ;
   }
   else if( OP == has_mod )
   {
      result = name() + "(SS)" ;
   }
   else if( OP == ext_mod )
   {
      result = name() + "(SS,SS,[,SS])" ;
   }
   else
   {
      result = name() + "(SS[,SS])" ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_ExtractionExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: print" ) ;

   if( !is_a_prototype() )
   {
      SRC->print( os, indent_width ) ;
   }
   else
   {
      MAC_Expression::print( os, indent_width ) ;
   }
}

//----------------------------------------------------------------------
MAC_Data const*
MAC_ExtractionExp:: data( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: data" ) ;
   MAC_CHECK( data_PRE( ct ) ) ;

   MAC_Data const* result = SRC ;

   MAC_CHECK_POST( data_POST( result, ct ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_ExtractionExp:: temporary_file( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: temporary_file" ) ;

   MAC_Communicator const* com = MAC_Exec::communicator() ;
   
   static int EXTRACTED_IDX = 0 ;
   std::ostringstream ss ;
   ss << MAC_System::working_directory()
      << MAC_System::path_name_separator()
      << "temporary_" << EXTRACTED_IDX++ ;
   if( com->nb_ranks() > 1 )
   {
      ss << "#" << com->rank() ;
   }
   ss << ".mac" ;
   static std::string result = "" ;
   result = ss.str() ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_ExtractionExp:: data_name( std::string const& exp_name,
                               MAC_Data const* d )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: data_name" ) ;
   MAC_CHECK( !exp_name.empty() ) ;
   MAC_CHECK( d != 0 ) ;
   
   if( ! d->value_can_be_evaluated( 0 ) )
   {
      std::ostringstream msg ;
      msg << "*** " << exp_name << " expression error:" << std::endl
          << "    the entry name cannot be defined from variables " << std::endl ;
      stringVector const& undef_vars = d->undefined_variables( 0 ) ;
      if( undef_vars.size() > 0 )
      {
         msg << "    unexpected variable(s): " << std::endl ;
         for( size_t i=0 ; i<undef_vars.size() ; ++i )
         {
            msg << "        - \"" << undef_vars(i) << "\"" << std::endl ;
         }
      }
      MAC_Error::object()->raise_plain( msg.str() ) ;
   }   
   return( d->to_string() ) ;
}

//----------------------------------------------------------------------
void
MAC_ExtractionExp:: extract_module( std::string const& file_name,
                                    std::string const& d_name,
                                    std::string const& m_name ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ExtractionExp:: extract_module" ) ;
   MAC_CHECK( OP == ext_mod ) ;
   MAC_CHECK( ! file_name.empty() ) ;
   MAC_CHECK( ! d_name.empty() ) ;
   MAC_CHECK( ! m_name.empty() ) ;
   MAC_CHECK( DB_MOD != 0 ) ;
   
   std::ofstream out( file_name.c_str(),
                      std::ios::out | std::ios::trunc ) ;
   if( !out )
   {
     raise_error( "   unable to create temporary file \""+file_name+"\"" ) ;
   }
   MAC_Module* mod = DB_MOD->module( d_name ) ;
   MAC_Module* dup = mod->create_clone(0) ;
   dup->modify_module_name( m_name ) ;
   dup->print( out, 0 ) ;
   out.close() ;
   dup->destroy() ; dup = 0 ;
}
