#include <MAC_Error.hh>

#include <MAC_Exec.hh>
#include <MAC_Exceptions.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_Int.hh>
#include <MAC_Map.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Root.hh>
#include <MAC_String.hh>
#include <MAC_StringVector.hh>
#include <MAC_assertions.hh>
#include <stringVector.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;


//---------------------------------------------------------------------------
MAC_Error*
MAC_Error:: object( void )
//---------------------------------------------------------------------------
{
   static MAC_Error UNIQUE_INSTANCE ;
   return( &UNIQUE_INSTANCE ) ;
}




//---------------------------------------------------------------------------
MAC_Error:: MAC_Error( void )
//---------------------------------------------------------------------------
   : os( MAC_Exec::out() )
{
}




//---------------------------------------------------------------------------
MAC_Error:: ~MAC_Error( void )
//---------------------------------------------------------------------------
{
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_internal( std::string message )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Internal Error" << endl << endl ;
   os << message << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_plain( std::string message )
//---------------------------------------------------------------------------
{
   begin() ;
   os << message << endl << endl ;
   end() ;
   // print_invocated_methods() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: display_info( std::string message )
//---------------------------------------------------------------------------
{
   os << endl << hline() << endl ;
   os << message << endl ;
   os << hline() << endl ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_not_tested( std::string file, int line, std::string test )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Situation never tested : " << test << endl << endl ;
   os << "in file : " << file << endl ;
   os << "at line : " << line << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_not_implemented( MAC_Object const* oo,
                                   std::string method )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Call to not implemented method : " << method << endl ;
   os << "         for an object of type : " << oo->type_name()
      << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_bad_types( MAC_Object const* oo,
                             std::string method,
                             MAC_Object const* argument )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Uncovered combination of run-time types" << endl ;
   os << "in call to method     : " << method << endl ;
   os << "for an object of type : " << oo->type_name() << endl ;
   os << "with argument attached to object of type : " << endl ;
   os << "   " << argument->type_name() << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_bad_types( MAC_Object const* oo,
                             std::string method,
                             MAC_Object const* argument_1,
                             MAC_Object const* argument_2 )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Uncovered combination of run-time types" << endl ;
   os << "in call to method     : " << method << endl ;
   os << "for an object of type : " << oo->type_name() << endl ;
   os << "with arguments attached to objects of type : " << endl ;
   os << "   " << argument_1->type_name() << endl ;
   os << "   " << argument_2->type_name() << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_read_syntax_error( std::string const& file,
                                    int line,
                                    std::string const& last,
                                    std::string const& nature )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Error while reading data file." << endl << endl ;
   os << nature << endl << endl ;
   os << "Last line number " << line << " read in file " << file << " was :" << endl ;
   os << ">> " << last << " <<" << endl ;
   end() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_missing_keyword( MAC_ModuleExplorer const* exp,
                                   std::string keyword )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   os << "the following keyword is requested: " << endl
      << "   \"" << keyword << "\"" << endl ;
   end() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_missing_module( MAC_ModuleExplorer const* exp,
                                  std::string path_and_name )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   os << "the following MODULE is requested: " << endl
      << "   \"" << path_and_name << "\"" << endl ;
   end() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_module_error( MAC_ModuleExplorer const* exp,
                                std::string error )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   os << error << endl ;
   end() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_bad_data_type( MAC_ModuleExplorer const* exp,
                                 std::string keyword,
                                 MAC_Data::Type query_kind )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module : " << exp->absolute_path_name() << endl ;
   os << "the data of keyword : " << keyword << endl ;
   os << "should be of type : "
      << MAC_Data::type_name(query_kind) << endl ;
   end() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_bad_file( MAC_ModuleExplorer const* exp,
                            std::string const& filename,
                            std::string const& access )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module : " << exp->absolute_path_name() << endl ;
   os << " file : " << exp->string_data( filename ) << endl ;

   os << " specified by the data of keyword : " << filename << endl ;

   os << " can't be accessed in mode : " << access << endl ;
   end() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_not_evaluable( MAC_ModuleExplorer const* exp,
                                 std::string keyword,
                                 stringVector const& undefined_variables )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module : " << exp->absolute_path_name() << endl ;
   os << "the data of keyword : " << keyword << endl ;
   os << "cannot be evaluated " << endl ;
   if( undefined_variables.size() > 0 )
   {
      os << "  undefined variable(s): " << endl ;
      for( size_t i=0 ; i<undefined_variables.size() ; ++i )
      {
         os << "      - \"" << undefined_variables(i) << "\"" << endl ;
      }
   }
   end() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_bad_data_value( MAC_ModuleExplorer const* exp,
                                  std::string keyword,
                                  std::string allowed_values )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   os << "the data of keyword: " << keyword << endl ;
   os << "should have one of the following values: " << endl ;
   os << allowed_values << endl ;
   end() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_data_error( MAC_ModuleExplorer const* exp,
                              std::string keyword,
                              std::string error )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   os << "error in the data of keyword: " << keyword << endl ;
   os << error << endl ;
   end() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_bad_object( std::string message,
                              MAC_Object const* a_object )
//---------------------------------------------------------------------------
{
   begin() ;
   os << message << endl ;
   a_object->print( os, 1 ) ;
   end() ;
   exit() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_precondition_violation( std::string test )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "PRECONDITION assertion violation :" << endl ;
   os << "   " << test << endl << endl ;
   print_invocated_methods() ;
   os << endl ;
   print_client_and_server() ;
   os << endl ;
   os << "HINTS : " << endl ;
   os << "   - Meeting a PRECONDITION is the CLIENT\'s responsibility."
      << endl ;
   os << "   - A PRECONDITION violation is the manifestation of " << endl
      << "     a bug in the CLIENT." << endl ;
   end() ;
   exit_with_internal_error() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_postcondition_violation( std::string file,
                                           int line,
                                           std::string test )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "POSTCONDITION assertion violation :" << endl ;
   os << "   " << test << endl ;
   os << "in file : " << file << endl ;
   os << "at line : " << line << endl << endl ;
   print_invocated_methods() ;
   os << endl ;
   print_client_and_server() ;
   os << endl ;
   os << "HINTS : " << endl ;
   os << "   - Meeting a POSTCONDITION is the SUPPLIER\'s responsibility."
      << endl ;
   os << "   - A POSTCONDITION violation is the manifestation of " << endl
      << "     a bug in the SUPPLIER." << endl ;
   end() ;
   exit_with_internal_error() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_invariant_violation( std::string file,
                                       int line,
                                       std::string test )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "INVARIANT assertion violation :" << endl ;
   os << "   " << test << endl ;
   os << "in file : " << file << endl ;
   os << "at line : " << line << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_assertion_violation( std::string file,
                                       int line,
                                       std::string test )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "ASSERTION violation :" << endl ;
   os << "   " << test << endl  ;
   os << "in file : " << file << endl ;
   os << "at line : " << line << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: raise_file_handling( std::string file,
                                 std::string operation )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Failure of operation : " << operation << endl ;
   os << "for file : " << file << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: begin( void )
//---------------------------------------------------------------------------
{
   os << endl << endl << hline() << endl ;
   os << "              FATAL ERROR" << endl << hline() << endl ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: end( void )
//---------------------------------------------------------------------------
{
   os << hline() << endl ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: print_invocated_methods( void )
//---------------------------------------------------------------------------
{
   if( MAC_Marker::nb_labels() > 0 )
   {
      os << "Stack of invocated methods :" << endl ;
      for( size_t i=MAC_Marker::nb_labels()-1 ; i<MAC_Marker::nb_labels() ; i-- )
      {
         os << "   " << MAC_Marker::label(i) << endl ;
      }
      os << "WARNING : only the methods labelled with \'MAC_LABEL\' "
         << "are listed." << endl ;
   }
   else
   {
      os << "Stack of invocated methods unavailable" << endl ;
      os << "(none of them have been labelled with \'MAC_LABEL\')."
         << endl ;
   }
}




//---------------------------------------------------------------------------
void
MAC_Error:: print_client_and_server( void )
//---------------------------------------------------------------------------
{
   if( MAC_Marker::nb_labels() > 1 )
   {
      os << "CLIENT   : " << MAC_Marker::label( MAC_Marker::nb_labels()-2 ) << endl ;
   }
   if( MAC_Marker::nb_labels() > 0 )
   {
      os << "SUPPLIER : " << MAC_Marker::label( MAC_Marker::nb_labels()-1 ) << endl ;
      os << "(inferred from the above stack of invocated methods)."
         << endl ;
   }
}




//---------------------------------------------------------------------------
std::string const&
MAC_Error:: hline( void ) const
//---------------------------------------------------------------------------
{
   static std::string const h = string( 50, '-' ) ;
   return( h ) ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: trace( std::string const& event )
//---------------------------------------------------------------------------
{
   os << endl << hline() << endl ;
   os << " EVENT TRACKING " ;
   os << event << endl << hline() << endl ;
   print_invocated_methods() ;
   os << hline() << endl << endl ;
}




//-------------------------------------------------------------------------
void
MAC_Error:: display_new_syntax( MAC_ModuleExplorer const* older_module,
                                MAC_ModuleExplorer const* result_module ) const
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Error:: display_new_syntax" ) ;

   os << " ********* MAC information " << endl ;
   os << " Syntax for module : " << endl ;
   os << "<==============================" << endl ;
   older_module->print( os, 0 ) ;
   os << "<==============================" << endl ;
   os << " won't be supported anymore. " << endl ;
   os << " Please let replace with following :" << endl ;
   os << ">==============================" << endl ;
   result_module->print( os, 0 ) ;
   os << ">==============================" << endl ;
   os << endl << endl ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: notify( MAC_Module* list,
                    std::string const& message,
                    std::string const& where,
                    stringVector const& choices )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Error:: notify list" ) ;
   MAC_CHECK_PRE( list!=0 ) ;
   MAC_CHECK_PRE( !message.empty() ) ;
   MAC_CHECK_PRE( !where.empty() ) ;

   MAC_Module* result = notify( list, message, where ) ;
   if(choices.size()>0)
      result->add_entry( "valid_choices",
                         MAC_StringVector::create( result, choices ) ) ;

}




//---------------------------------------------------------------------------
MAC_Module*
MAC_Error:: notify( MAC_Module* list,
                    std::string const& message,
                    std::string const& where )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Error:: notify simple" ) ;
   MAC_CHECK_PRE( list!=0 ) ;
   MAC_CHECK_PRE( !message.empty() ) ;
   MAC_CHECK_PRE( !where.empty() ) ;

   size_t nb = 0 ;
   if( list->has_entry( "nb_items" ) )
      nb = list->data_of_entry( "nb_items" )->to_int() ;
   else
      list->add_entry( "nb_items", MAC_Int::create( list, nb ) ) ;

   std::ostringstream is ;
   is << "Error" << nb ;
   MAC_Module* result = MAC_Module::create( list, is.str() ) ;
   result->add_entry( "message", MAC_String::create( result, message ) ) ;
   result->add_entry( "where", MAC_String::create( result, where ) ) ;

   nb++ ;

   const_cast<MAC_Int*>(
      static_cast<MAC_Int const*>(list->data_of_entry( "nb_items" )))
      ->set( (int)nb ) ;
   list->add_module( result ) ;

   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->has_entry( "message" ) ) ;
   MAC_CHECK_POST( result->has_entry( "where" ) ) ;

   return result ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: display_data_checking( MAC_ModuleExplorer* report ) const
//---------------------------------------------------------------------------
{
   if( !report->is_empty() )
   {
      os << "Data deck checking failed : " << endl ;
      os << "*************************" << endl ;
      os << endl ;

      for( report->start_module_iterator() ;
           report->is_valid_module() ;
           report->go_next_module() )
      {
         MAC_ModuleExplorer const* error =
            report->create_subexplorer( 0 ) ;
         os << "   " << error->string_data( "message" ) << endl ;
         os << "   in module " << error->string_data( "where" ) << endl ;
         if( error->has_entry( "valid_choices" ) )
         {
            stringVector val = error->stringVector_data( "valid_choices" ) ;
            val.sort() ;
            os << "   Valid ones are :" << endl ;
            for( size_t j=0 ; j<val.size() ; ++j )
            {
               os << "      - \"" << val(j) << "\"" << endl ;
            }
         }
         os << endl ;
         error->destroy() ; error = 0 ;
      }
   }
}




//---------------------------------------------------------------------------
void
MAC_Error:: exit_with_internal_error( void )
//---------------------------------------------------------------------------
{
   MAC_Exec::set_exit_code( -1 ) ;
   throw MAC_Exceptions::InternalError() ;
}




//---------------------------------------------------------------------------
void
MAC_Error:: exit( size_t status )
//---------------------------------------------------------------------------
{
   MAC_Exec::set_exit_code( status ) ;
   throw MAC_Exceptions::UserError() ;
}



