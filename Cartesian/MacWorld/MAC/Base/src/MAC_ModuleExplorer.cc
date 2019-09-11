#include <MAC_ModuleExplorer.hh>

#include <MAC_assertions.hh>
#include <MAC_Context.hh>
#include <MAC_ContextPair.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_Error.hh>
#include <MAC_KeywordDataPair.hh>
#include <MAC_KeywordDataIterator.hh>
#include <MAC_List.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_ModulePattern.hh>
#include <MAC_String.hh>
#include <MAC_System.hh>
#include <MAC_Variable.hh>
#include <boolVector.hh>
#include <boolArray2D.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <stringVector.hh>
#include <stringArray2D.hh>

#include <fstream>
#include <sstream>
using std::string ;

//--------------------------------------------------------------------------
MAC_ModuleExplorer*
MAC_ModuleExplorer:: create( MAC_Object* a_owner,
                             MAC_Module const* mm,
                             PatternStatus status )
//--------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: create" ) ;
   MAC_CHECK_PRE( mm != 0 ) ;
   MAC_CHECK_PRE( IMPLIES( status==verify || status==build,
                           MAC_ModulePattern::has_pattern() ) ) ;
   MAC_CHECK_PRE( IMPLIES( status==build,
                           MAC_ModulePattern::build_pattern() ) ) ;

   MAC_ModulePattern* pat = 0 ;
   if( status!=ignore )
   {
      pat = MAC_ModulePattern::create( 0, mm, mm->name() ) ;
   }
   
   MAC_ModuleExplorer* result =
                       new MAC_ModuleExplorer( a_owner, mm, status, pat ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->pattern_status() == status ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
MAC_ModuleExplorer:: MAC_ModuleExplorer( MAC_Object* a_owner,
                                         MAC_Module const* mm,
                                         PatternStatus status,
                                         MAC_ModulePattern* pat )
//---------------------------------------------------------------------------
   : MAC_Object( a_owner ),
     mod( mm ),
     submodule_iterator( mm->create_module_iterator( this) ),
     keyword_iterator( mm->create_entry_iterator( this) ),
     tmp_context( MAC_ContextPair::create( this,0,0) ),
     MP( pat ),
     keyword_iterator_started(false),
     my_status( status )
{
   if( pat!=0 ) pat->set_owner( this ) ;
   
   MAC_CHECK_INV( invariant() ) ;
}




//--------------------------------------------------------------------------
MAC_ModuleExplorer*
MAC_ModuleExplorer:: create_subexplorer( MAC_Object* a_owner,
                                         std::string const& path_and_name ) const
//--------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: create_subexplorer( MAC_Object*, std::string const& )" ) ;
   MAC_CHECK_PRE( !path_and_name.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   if( !mod->has_module( path_and_name ) ) 
   {
      MAC_Error::object()->raise_missing_module( this, path_and_name ) ;
   }
   MAC_Module const* mm = mod->module( path_and_name ) ;
   MAC_ModulePattern* pat = 0 ;
   if( pattern()!=0 )
   {
      pat = pattern()->create_subpattern( 0, mod, path_and_name ) ;
   }
   
   MAC_ModuleExplorer* result =
                      new MAC_ModuleExplorer( a_owner, mm, my_status, pat ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->pattern_status() == pattern_status() ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK( result->mod == mm ) ;
   MAC_CHECK_POST( result->name() == MAC_Module::basename( path_and_name ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
MAC_ModuleExplorer*
MAC_ModuleExplorer:: create_clone( MAC_Object* a_owner ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: create_clone" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_ModulePattern* pat = 0 ;
   if( pattern()!=0 )
   {
      pat = pattern()->create_clone( 0 ) ;
   }
   MAC_ModuleExplorer* result =
                     new MAC_ModuleExplorer( a_owner, mod, my_status, pat ) ;
   
   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   MAC_CHECK_POST( result->pattern_status() == pattern_status() ) ;
   MAC_CHECK_INV( invariant() ) ;
   return result ;
}




//---------------------------------------------------------------------------
MAC_ModuleExplorer:: ~MAC_ModuleExplorer( void )
//---------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   mod = 0 ;
   MP = 0 ;
}




//--------------------------------------------------------------------------
MAC_ModuleExplorer::PatternStatus
MAC_ModuleExplorer:: pattern_status( void ) const
//--------------------------------------------------------------------------
{
   if( my_status != ignore && !MAC_ModulePattern::has_pattern() )
   {
      MAC_Error::object()->raise_internal(
         "*** MAC_ModuleExplorer error:\n"
         "    pattern file has already been closed" ) ;
   }
   return( my_status ) ;
}




//---------------------------------------------------------------------------
MAC_ModulePattern*
MAC_ModuleExplorer:: pattern( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: pattern" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_ModulePattern* result = MP ;
   
   MAC_CHECK_POST( IMPLIES( build_pattern(), result!=0 ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( result ) ;
}




//--------------------------------------------------------------------------
bool
MAC_ModuleExplorer:: build_pattern( void ) const
//--------------------------------------------------------------------------
{
   bool result = my_status == build &&
                 !MAC_Assertion::is_checking() ;
   if( my_status == build && !MAC_ModulePattern::build_pattern() )
   {
      MAC_Error::object()->raise_internal(
         "*** MAC_ModuleExplorer error:\n"
         "    pattern file has already been closed" ) ;
   }
   
   return( result ) ;
}




//--------------------------------------------------------------------------
MAC_ModuleExplorer*
MAC_ModuleExplorer:: validity( MAC_Object* a_owner,
                               bool recurse ) const
//--------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: validity" ) ;
   MAC_CHECK_PRE( pattern_status()==verify ) ;
   MAC_CHECK_PRE( MAC_ModulePattern:: has_pattern() ) ;

   MAC_ModuleExplorer* result  = 0 ;
   MAC_ModulePattern const* pat =
      MAC_ModulePattern::create_pattern_description( 0, mod->name() ) ;
   if( pat == 0 )
   {
      MAC_Error::object()->raise_plain(
         "*** pattern verification failed\n"
         "    unable to find description of:\n"
         "       module: "+mod->name()+"\n"
         "    the control file may not exist or may be corrupted" ) ;
   }
   else
   {
      MAC_Module* m = pat->validity( 0, mod, recurse ) ;
      result  = MAC_ModuleExplorer::create( a_owner, m ) ;
      m->set_owner( result ) ;
      pat->destroy() ; pat = 0 ;
   }

   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->owner()==a_owner ) ;   
   return( result ) ;
}




//---------------------------------------------------------------------------
bool 
MAC_ModuleExplorer:: has_module( std::string const& a_path_and_name ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: has_module" ) ;
   MAC_CHECK_PRE( !a_path_and_name.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;

   bool result = mod->has_module( a_path_and_name ) ;

   if( build_pattern() && result )
   {
      std::string dir_name = MAC_Module::dirname( a_path_and_name ) ;
      std::string base_name = MAC_Module::basename( a_path_and_name ) ;
      MAC_ASSERT( !base_name.empty() ) ;
      if( dir_name.empty() )
      {
         pattern()->add_pattern( a_path_and_name,
                                 mod->module( a_path_and_name ),
                                 MAC_ModulePattern::optional ) ;
      }
      else
      {
         MAC_ModuleExplorer const* exp = create_subexplorer( 0, dir_name ) ;
	 MAC_ASSERT( exp->has_module( base_name ) ) ;
	 exp->destroy() ; exp = 0 ;
      }
   }
   
   MAC_CHECK_INV( invariant() ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
bool
MAC_ModuleExplorer:: is_empty( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: is_empty" ) ;
   return( mod->is_empty() ) ;
}




//---------------------------------------------------------------------------
MAC_Context const*
MAC_ModuleExplorer:: context( std::string const& a_path_and_name,
                              MAC_Context const* ct ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: context" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_Module const* m = mod ;
   if( !a_path_and_name.empty() )
   {
      m = mod->module( a_path_and_name ) ;
   }
   MAC_Context const* result = m->context() ;
   if( ct!=0 )
   {
      if( build_pattern() )
      {
         stringVector vars( ct->nb_variables() ) ;
         for( size_t i=0 ; i<ct->nb_variables() ; ++i )
         {
            MAC_Variable const* var = ct->variable(i) ;
            vars(i) = var->name() ;
         }
         vars.sort() ;
         for( size_t i=0 ; i<vars.size() ; ++i )
         {
            std::string const& str = vars(i) ;
            MAC_Variable const* var = MAC_Variable::object( str ) ;
            pattern()->add_provided_variable( str, m, var->data_type() ) ;
         }
      }
      tmp_context->re_initialize( m->context(), ct ) ;
      result = tmp_context ;
   }
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( result!=0 ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
bool 
MAC_ModuleExplorer:: has_entry( std::string const& a_path_and_keyword ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: has_entry" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;

   bool result = mod->has_entry( a_path_and_keyword ) ;
   if( result && build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::Undefined ) ;      
   }
   
   MAC_CHECK_INV( invariant() ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: start_module_iterator( void ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: start_module_iterator" ) ;
   MAC_CHECK_INV( invariant() ) ;
   submodule_iterator->start() ;
   MAC_CHECK_INV( invariant() ) ;
}




//---------------------------------------------------------------------------
bool
MAC_ModuleExplorer:: is_valid_module( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: is_valid_module" ) ;
   MAC_CHECK_INV( invariant() ) ;
   return submodule_iterator->is_valid() ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: go_next_module( void ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: go_next_module" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( is_valid_module() ) ;
   submodule_iterator->go_next() ;
}




//---------------------------------------------------------------------------
MAC_ModuleExplorer*
MAC_ModuleExplorer:: create_subexplorer( MAC_Object* a_owner ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: create_subexplorer( MAC_Object* )" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( is_valid_module() ) ;

   MAC_ModulePattern* pat = 0 ;
   if( pattern()!=0 )
   {
      pat = pattern()->generic_pattern( 0, submodule_iterator->item() ) ;
   }
   
   MAC_ModuleExplorer* result =
      new MAC_ModuleExplorer( a_owner, submodule_iterator->item(), my_status, pat ) ;
   
   MAC_CHECK_POST( has_module( result->name() ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   return result ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: start_entry_iterator( void ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: start_entry_iterator" ) ;
   MAC_CHECK_INV( invariant() ) ;
   keyword_iterator->start() ;
   keyword_iterator_started = true ;
}




//---------------------------------------------------------------------------
bool
MAC_ModuleExplorer:: is_valid_entry( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: is_valid_entry" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   return keyword_iterator->is_valid() ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: go_next_entry( void ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: go_next_entry" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_PRE( is_valid_entry() ) ;
   keyword_iterator->go_next() ;
}




//---------------------------------------------------------------------------
std::string const&
MAC_ModuleExplorer:: keyword( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: keyword" ) ;
   MAC_CHECK_PRE( is_valid_entry() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   string const& result = keyword_iterator->item()->keyword() ;
   
   MAC_CHECK_POST( has_entry( result ) ) ;
   return result ;
}




//---------------------------------------------------------------------------
MAC_DataWithContext*
MAC_ModuleExplorer:: data( MAC_Object* a_owner,
                           MAC_Context const* ct ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: data" ) ;
   MAC_CHECK_PRE( is_valid_entry() ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Data const* a_data = keyword_iterator->item()->data() ;
   MAC_DataWithContext* result =
      MAC_DataWithContext::create( a_owner,
                                   a_data,
                                   context( "", ct ) ) ;
   if( build_pattern() )
   {
      pattern()->add_generic_keyword( keyword_iterator->item()->keyword(),  mod, a_data->data_type() ) ;
   }
   
   MAC_CHECK_POST( result !=0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
MAC_DataWithContext*
MAC_ModuleExplorer:: abstract_data( MAC_Object* a_owner,
                                    std::string const& a_path_and_keyword,
                                    MAC_Context const* ct ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: abstract_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   MAC_Data const* a_data = mod->data_of_entry( a_path_and_keyword ) ;
   
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, a_data->data_type() ) ;
   }
   MAC_DataWithContext* result =
      MAC_DataWithContext::create(
         a_owner,
         a_data,
         context( MAC_Module::dirname( a_path_and_keyword ), ct ) ) ;
   
   MAC_CHECK_POST( result !=0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
bool
MAC_ModuleExplorer:: bool_data( std::string const& a_path_and_keyword,
                                MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: bool_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   MAC_Data const* val = mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::Bool ) ;
   }
   if( val->data_type()!=MAC_Data::Bool )
   {
      MAC_Error::object()->raise_bad_data_type(
                              this, a_path_and_keyword, MAC_Data::Bool  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
         this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   return( val->to_bool( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
int
MAC_ModuleExplorer:: int_data( std::string const& a_path_and_keyword,
                               MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: int_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   MAC_Data const* val = mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::Int ) ;
   }
   if( val->data_type()!=MAC_Data::Int )
   {
      MAC_Error::object()->raise_bad_data_type(
                                this, a_path_and_keyword, MAC_Data::Int ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
         this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_int( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
double
MAC_ModuleExplorer:: double_data( std::string const& a_path_and_keyword,
                                  MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: double_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   MAC_Data const* val = mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::Double ) ;
   }
   if( val->data_type()!=MAC_Data::Double )
   {
      MAC_Error::object()->raise_bad_data_type(
                               this, a_path_and_keyword, MAC_Data::Double ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
         this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
      
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_double( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
std::string const&
MAC_ModuleExplorer:: string_data( std::string const& a_path_and_keyword,
                                  MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: string_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::String ) ;
   }
   if( val->data_type()!=MAC_Data::String )
   {
      MAC_Error::object()->raise_bad_data_type(
                                this, a_path_and_keyword, MAC_Data::String  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
         this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_string( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
intVector const&
MAC_ModuleExplorer:: intVector_data( std::string const& a_path_and_keyword,
                                     MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: intVector_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::IntVector ) ;
   }
   if( val->data_type()!=MAC_Data::IntVector )
   {
      MAC_Error::object()->raise_bad_data_type(
                            this, a_path_and_keyword, MAC_Data::IntVector  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
              this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_int_vector( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
doubleVector const&
MAC_ModuleExplorer:: doubleVector_data( std::string const& a_path_and_keyword,
                                        MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: doubleVector_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::DoubleVector ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=MAC_Data::DoubleVector )
   {
      MAC_Error::object()->raise_bad_data_type(
                         this, a_path_and_keyword, MAC_Data::DoubleVector  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ), ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
              this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_double_vector( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
doubleArray2D const&
MAC_ModuleExplorer:: doubleArray2D_data( std::string const& a_path_and_keyword,
                                         MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: doubleArray2D_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::DoubleArray2D ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=MAC_Data::DoubleArray2D )
   {
      MAC_Error::object()->raise_bad_data_type(
                        this, a_path_and_keyword, MAC_Data::DoubleArray2D  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
            this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_double_array2D( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
boolArray2D const&
MAC_ModuleExplorer:: boolArray2D_data( std::string const& a_path_and_keyword,
                                         MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: boolArray2D_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::BoolArray2D ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=MAC_Data::BoolArray2D )
   {
      MAC_Error::object()->raise_bad_data_type(
                        this, a_path_and_keyword, MAC_Data::BoolArray2D  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
            this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_bool_array2D( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
stringArray2D const&
MAC_ModuleExplorer:: stringArray2D_data( std::string const& a_path_and_keyword,
                                         MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: stringArray2D_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::StringArray2D ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=MAC_Data::StringArray2D )
   {
      MAC_Error::object()->raise_bad_data_type(
                        this, a_path_and_keyword, MAC_Data::StringArray2D  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
            this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_string_array2D( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
intArray2D const&
MAC_ModuleExplorer:: intArray2D_data( std::string const& a_path_and_keyword,
                                      MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: intArray2D_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::IntArray2D ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=MAC_Data::IntArray2D )
   {
      MAC_Error::object()->raise_bad_data_type(
                           this, a_path_and_keyword, MAC_Data::IntArray2D  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
            this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_int_array2D( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
doubleArray3D const&
MAC_ModuleExplorer:: doubleArray3D_data( std::string const& a_path_and_keyword,
                                         MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: doubleArray3D_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::DoubleArray3D ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=MAC_Data::DoubleArray3D )
   {
      MAC_Error::object()->raise_bad_data_type(
                      this, a_path_and_keyword, MAC_Data::DoubleArray3D  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
            this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_double_array3D( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
intArray3D const&
MAC_ModuleExplorer:: intArray3D_data( std::string const& a_path_and_keyword,
                                      MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: intArray3D_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::IntArray3D ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( val->data_type()!=MAC_Data::IntArray3D )
   {
      MAC_Error::object()->raise_bad_data_type(
                          this, a_path_and_keyword, MAC_Data::IntArray3D  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
         this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_int_array3D( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
boolVector const&
MAC_ModuleExplorer:: boolVector_data( std::string const& a_path_and_keyword,
                                      MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: boolVector_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::BoolVector ) ;
   }
   if( val->data_type()!=MAC_Data::BoolVector )
   {
      MAC_Error::object()->raise_bad_data_type(
                           this, a_path_and_keyword, MAC_Data::BoolVector  ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
              this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_bool_vector( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
stringVector const&
MAC_ModuleExplorer:: stringVector_data( std::string const& a_path_and_keyword,
                                        MAC_Context const* ct  ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: stringVector_data" ) ;
   MAC_CHECK_PRE( !a_path_and_keyword.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_path_and_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_path_and_keyword ) ;
   }
   MAC_Data const* val =  mod->data_of_entry( a_path_and_keyword ) ;
   if( build_pattern() )
   {
      declare_data( a_path_and_keyword, MAC_Data::StringVector ) ;
   }
   if( val->data_type()!=MAC_Data::StringVector )
   {
      MAC_Error::object()->raise_bad_data_type(
                          this, a_path_and_keyword, MAC_Data::StringVector ) ;
   }
   MAC_Context const* a_ctx =
                  context( MAC_Module::dirname( a_path_and_keyword ),ct ) ;
   if( !val->value_can_be_evaluated( a_ctx ) )
   {
      MAC_Error::object()->raise_not_evaluable(
               this, a_path_and_keyword, val->undefined_variables( a_ctx ) ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   return( val->to_string_vector( a_ctx ) ) ;
}




//---------------------------------------------------------------------------
MAC_Module*
MAC_ModuleExplorer:: create_clone_of_attached_module(
   MAC_Object* a_owner ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: create_clone_of_attached_module" ) ;

   MAC_Module* result = mod->create_clone( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   
   return result ;
}




//---------------------------------------------------------------------------
std::string const&
MAC_ModuleExplorer:: name( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: name" ) ;
   MAC_CHECK_INV( invariant() ) ;
   std::string const& result = mod->name() ;
   MAC_CHECK_POST( !result.empty() ) ;
   MAC_CHECK_POST( MAC_Module::dirname( result ).empty() ) ;
   MAC_CHECK_POST( result==MAC_Module::basename( absolute_path_name() ) ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
std::string const&
MAC_ModuleExplorer:: absolute_path_name( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: absolute_path_name" ) ;
   MAC_CHECK_INV( invariant() ) ;
   std::string const& result = mod->absolute_path_name() ;
   MAC_CHECK_POST( !result.empty() ) ;
   MAC_CHECK_POST( MAC_Module::basename( result )==name() ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
MAC_Object const*
MAC_ModuleExplorer:: owner_of_attached_module( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: owner_of_attached_module" ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( mod->owner() ) ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: print( std::ostream& os, size_t indent_width ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: print" ) ;
   MAC_CHECK_INV( invariant() ) ;
   mod->print( os, indent_width ) ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: write( std::string const& file,
                            std::string const& format ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: write" ) ;
   MAC_CHECK_INV( invariant() ) ;
   mod->write( file, format ) ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: test_data( std::string const& a_keyword,
                                std::string const& expression,
                                MAC_Context const* ct ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: test_data" ) ;
   MAC_CHECK_PRE( !a_keyword.empty() ) ;
   MAC_CHECK_PRE( a_keyword.find("/") > a_keyword.length() ) ;
   MAC_CHECK_PRE( !expression.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   if( !mod->has_entry( a_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_keyword ) ;
   }

   MAC_Data const* eval = mod->create_evaluation( 0, expression, ct ) ;
   if( eval==0 || eval->data_type()!=MAC_Data::Bool )
   {
      MAC_Error::object()->raise_data_error(
         this, a_keyword,
         "   condition ("+expression+")\n"
         "   can't be evaluated or is not a boolean expression" ) ;
   }

   if( !eval->to_bool() )
   {
      std::ostringstream msg ;
      msg << "   condition ("+expression+")\n" ;
      msg << "   is not verified:\n" ;
      msg << "      (" ;
      eval->print( msg, 0 ) ;
      msg << ")" ;
      MAC_Error::object()->raise_data_error(
                                     this , a_keyword, msg.str() ) ;
   }
   eval->destroy() ;

   if( build_pattern() )
   {
      pattern()->attach_verify_data( a_keyword, expression ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: test_data_in( std::string const& a_keyword,
                                   stringVector const& choices ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: test_data_in" ) ;
   MAC_CHECK_PRE( !a_keyword.empty() ) ;
   MAC_CHECK_PRE( a_keyword.find("/") > a_keyword.length() ) ;
   MAC_CHECK_PRE( choices.size() > 0 ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( !mod->has_entry( a_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_keyword ) ;
   }

   MAC_Data const* eval = mod->data_of_entry( a_keyword ) ;

   if( eval->data_type() != MAC_Data::String &&
       eval->data_type() != MAC_Data::Bool &&
       eval->data_type() != MAC_Data::Int  &&
       eval->data_type() != MAC_Data::StringVector &&
       eval->data_type() != MAC_Data::BoolVector &&
       eval->data_type() != MAC_Data::IntVector )
   {
      // Double are forbidden because the equality between two doubles
      // is not clear...
      MAC_Error::object()->raise_data_error(
         this, a_keyword,
         "*** MAC_ModuleExplorer:: test_data_in error\n"
         "    test only available for:\n"
         "        - string values\n"
         "        - boolean values\n"
         "        - integer values\n"
         "        - string vector values\n"
         "        - boolean vector values\n"
         "        - integer vector values" ) ;
   }

   std::string str = eval->value_as_string( mod->context() ) ;
   MAC::remove_enclosing_characters(str,'\"') ;
   
   if( !choices.has( str ) )
   {
      std::string mess ;
      for( size_t i=0 ; i<choices.size() ; i++ )
      {
         mess += "   - \""+choices(i)+"\"\n" ;
      }
      MAC_Error::object()->raise_bad_data_value(
         this, a_keyword, mess ) ;
   }

   if( build_pattern() && a_keyword != "type" && a_keyword != "concrete_name" )
   {
      pattern()->attach_list_of_valid_choices( a_keyword, choices ) ;
   }
   
   
   MAC_CHECK_INV( invariant() ) ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: test_data_as( std::string const& a_keyword,
                                   std::string const& regexp,
                                   std::string const& where ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: test_data_as" ) ;
   MAC_CHECK_PRE( !regexp.empty() ) ;
   MAC_CHECK_PRE( !a_keyword.empty() ) ;
   MAC_CHECK_PRE( a_keyword.find("/") > a_keyword.length() ) ;
   
   MAC_CHECK_INV( invariant() ) ;
   
   MAC_List* list = mod->create_data_selection( 0, regexp, 0, where ) ;
         
   stringVector choices(0) ;
   MAC_ListIterator * it = list->create_iterator( list ) ;
   for( it->start() ; it->is_valid() ; it->go_next() ) 
   {
      MAC_Data const* adata = static_cast<MAC_Data*>( it->item() ) ;
      if( adata->data_type()!=MAC_Data::String )
      {
         MAC_Error::object()->raise_bad_data_type(
            this, regexp, MAC_Data::String ) ;
      }
      choices.append( adata->to_string() ) ;
   }
   
   list->destroy() ;
   if( !mod->has_entry( a_keyword ) ) 
   {
      MAC_Error::object()->raise_missing_keyword( this, a_keyword ) ;
   }

   MAC_Data const* eval = mod->data_of_entry( a_keyword ) ;
   if( eval->data_type()!=MAC_Data::String )
   {
      MAC_Error::object()->raise_bad_data_type(
         this, a_keyword, MAC_Data::String ) ;
   }

   std::string const& str = eval->to_string( mod->context() ) ;
   
   if( !choices.has( str ) )
   {
      std::string mess ;
      for( size_t i=0 ; i<choices.size() ; i++ )
      {
         mess += "   - \""+choices(i)+"\"\n" ;
      }
      MAC_Error::object()->raise_bad_data_value(
         this, a_keyword, mess ) ;
   }

   if( build_pattern() )
   {
      pattern()->attach_list_of_dynamic_choices( a_keyword, regexp, where ) ;
   }
   
   
   MAC_CHECK_INV( invariant() ) ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: set_default( std::string const& a_keyword,
                                  std::string const& expression ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: set_default" ) ;
   MAC_CHECK_PRE( !a_keyword.empty() ) ;
   MAC_CHECK_PRE( a_keyword.find("/") > a_keyword.length() ) ;
   MAC_CHECK_PRE( !expression.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( build_pattern() )
   {
      if( !mod->has_entry( a_keyword ) ) 
      {
         MAC_Error::object()->raise_missing_keyword( this, a_keyword ) ;
      }
      pattern()->attach_default_data( a_keyword, expression ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: set_help( std::string const& a_keyword,
                               std::string const& expression ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: set_help" ) ;
   MAC_CHECK_PRE( !a_keyword.empty() ) ;
   MAC_CHECK_PRE( a_keyword.find("/") > a_keyword.length() ) ;
   MAC_CHECK_PRE( !expression.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( build_pattern() )
   {
      if( !mod->has_entry( a_keyword ) ) 
      {
         MAC_Error::object()->raise_missing_keyword( this, a_keyword ) ;
      }
      pattern()->attach_help_data( a_keyword, expression ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: test_file( std::string const& filename,
                                std::string const& mode ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: test_file" ) ;
   MAC_CHECK_PRE( !filename.empty() ) ;
   MAC_CHECK_PRE( filename.find("/") > filename.length() ) ;
   MAC_CHECK_PRE( mode=="read" || mode=="write" ) ;
   MAC_CHECK_INV( invariant() ) ;

   std::string const& fn = string_data( filename ) ;
   
   if( ( mode=="read" && !( MAC_System::can_read( fn ) ) ) ||
       ( mode=="write" && !( MAC_System::can_write( fn ) ) ) )
   {
      MAC_Error::object()->raise_bad_file( this, filename, mode ) ;
   }  
   
   if( build_pattern() )
   {
      if( !mod->has_entry( filename ) ) 
      {
         MAC_Error::object()->raise_missing_keyword( this, filename ) ;
      }
      pattern()->attach_file_extension( filename, mode ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
}




//---------------------------------------------------------------------------
void
MAC_ModuleExplorer:: declare_data( std::string const& a_path_and_keyword,
                                   MAC_Data::Type type ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleExplorer:: declare_data" ) ;
   MAC_CHECK( !a_path_and_keyword.empty() ) ;
   MAC_CHECK( my_status == build ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_ModulePattern::Access acc =  MAC_ModulePattern::mandatory ;
   if( type==MAC_Data::Undefined )
   {
      acc=MAC_ModulePattern::optional ;
   }
   
   if( keyword_iterator_started &&
       is_valid_entry() &&
       keyword()==MAC_Module::basename( a_path_and_keyword ) )
   {
      acc=MAC_ModulePattern::generic ;
   }
   
   pattern()->add_entry( a_path_and_keyword, mod, type, acc ) ;
}




//---------------------------------------------------------------------------
bool
MAC_ModuleExplorer:: invariant( void ) const
//---------------------------------------------------------------------------
{
   MAC_ASSERT( mod!=0 ) ;
   MAC_ASSERT( IMPLIES( my_status==build, MP!=0  ) ) ;
   MAC_ASSERT( IMPLIES( MP!=0, my_status!=ignore ) ) ;
   return( true ) ;
}
