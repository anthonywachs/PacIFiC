#include <MAC_Module.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <stringVector.hh>

#include <MAC_BinStored.hh>
#include <MAC_Context.hh>
#include <MAC_ContextPair.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_Error.hh>
#include <MAC_Int.hh>
#include <MAC_KeywordDataPair.hh>
#include <MAC_Root.hh>
#include <MAC_KeywordDataPair.hh>
#include <MAC_KeywordDataIterator.hh>
#include <MAC_ModuleComparator.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_List.hh>
#include <MAC_String.hh>
#include <MAC_System.hh>
#include <MAC_Variable.hh>
#include <MAC_Vector.hh>

#include <fstream>
#include <iomanip>
#include <sstream>

// Interface with Yacc
bool MAC_readFile( MAC_Module * top,
                   std::istream* input_stream,
                   std::string const& name,
                   bool debug ) ;
std::string MAC_current_parsed_module_path_name( void ) ;



//----------------------------------------------------------------------
MAC_Module*
MAC_Module:: create( MAC_Object* a_owner,
                     std::string const& a_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: create(a_owner,a_name)" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   MAC_CHECK_PRE( a_name.find("/") >= a_name.size() ) ;
   MAC_CHECK_PRE( a_name.find(" ") >= a_name.size() ) ;

   MAC_Module* result = new MAC_Module( a_owner, a_name, 0 ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->name() == a_name ) ;
   MAC_CHECK_POST( result->is_empty() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Module*
MAC_Module:: create( MAC_Object* a_owner, 
                     std::string const& a_name,
                     std::string const& file_name,
                     MAC_Context const* ct )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: create(a_owner,a_name,file_name)" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   MAC_CHECK_PRE( a_name.find("/") >= a_name.size() ) ;
   MAC_CHECK_PRE( a_name.find(" ") >= a_name.size() ) ;
   MAC_CHECK_PRE( !file_name.empty() ) ;
   
   MAC_Module* result = new MAC_Module( a_owner, a_name, ct ) ;
   bool ok = MAC_readFile( result, 0, file_name, false ) ;
   
   if( !ok )
   {
      MAC_Error::object()->raise_plain(
         "*** MAC_Module error:\n"
         "    unable to complete data deck reading \""+file_name+"\"" ) ;
   }
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->name() == a_name ) ;
   MAC_CHECK_POST( !result->is_empty() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Module*
MAC_Module:: create( MAC_Object* a_owner, 
                     std::string const& a_name,
                     std::istream& input_stream,
                     MAC_Context const* ct )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: create(a_owner,a_name,input_stream)" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   MAC_CHECK_PRE( a_name.find("/") >= a_name.size() ) ;
   MAC_CHECK_PRE( a_name.find(" ") >= a_name.size() ) ;

   MAC_Module* result = new MAC_Module( a_owner, a_name, ct ) ;
   bool ok = MAC_readFile( result, &input_stream, "", false ) ;

   if( !ok )
   {
      MAC_Error::object()->raise_plain(
         "*** MAC_Module error:\n"
         "    unable to complete data deck reading" ) ;
   }
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->name() == a_name ) ;
   MAC_CHECK_POST( !result->is_empty() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Module*
MAC_Module:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: create_clone" ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Module* result = create_copy( a_owner ) ;
   if( FATHER != 0 )
   {
      complete_context( result, result, FATHER->context() ) ;
   }  

   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   MAC_CHECK_POST( result->name() == name() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Module*
MAC_Module:: create_as_difference( MAC_Object* a_owner,
				   std::string const& a_name,
                                   MAC_Module const* m1,
                                   MAC_Module const* m2,
				   MAC_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: create_as_difference" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   MAC_CHECK_PRE( a_name.find("/") >= a_name.size() ) ;
   MAC_CHECK_PRE( a_name.find(" ") >= a_name.size() ) ;
   MAC_CHECK_PRE( m1 != 0 ) ;
   MAC_CHECK_PRE( m2 != 0 ) ;
   
   MAC_Module* result = MAC_Module::create( a_owner, a_name ) ;
   MAC_ModuleComparator* cmp = MAC_ModuleComparator::create( a_owner, exp ) ;
   int nb_err = cmp->compare( m1, m2, result ) ;
   result->add_entry( "nb_differences", MAC_Int::create( result, nb_err ) ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->name() == a_name ) ;
   MAC_CHECK_POST( result->has_entry( "nb_differences" ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Module:: MAC_Module( MAC_Object* a_owner,
                         std::string const& a_name,
                         MAC_Context const* ct )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , NAME( MAC_String::create( this, a_name ) )
   , FATHER( 0 )
   , MODS( MAC_List::create( this ) )
   , ENTRIES( MAC_List::create( this ) )
   , CTX( MAC_ContextSimple::create( this ) )
   , TMP_CTX( MAC_ContextPair::create( this, 0, 0) )
{
   MAC_LABEL( "MAC_Module:: MAC_Module" ) ;
   
   if( ct!=0 ) CTX->extend( ct ) ;
   
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
MAC_Module:: ~MAC_Module( void ) 
//----------------------------------------------------------------------
{
}




//----------------------------------------------------------------------
void
MAC_Module:: modify_module_name( std::string const& a_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: modify_module_name" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   MAC_CHECK_PRE( a_name.find("/") >= a_name.size() ) ;
   MAC_CHECK_PRE( a_name.find(" ") >= a_name.size() ) ;
   MAC_CHECK_PRE( IMPLIES( a_name != name() && father() != 0,
                           !father()->has_module( a_name ) ) ) ;
   MAC_CHECK_PRE( IMPLIES( a_name != name() && father() != 0,
                           !father()->has_entry( a_name ) ) ) ;
   
   NAME->set( a_name ) ;

   MAC_CHECK_POST( name() == a_name ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: add_entry( std::string const& keyword,
                        MAC_Data const* data )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: add_entry" ) ;
   MAC_CHECK_PRE( !keyword.empty() ) ;
   MAC_CHECK_PRE( !has_module( keyword ) ) ;
   MAC_CHECK_PRE( keyword.find("/") >= keyword.size() ) ;
   MAC_CHECK_PRE( keyword.find(" ") >= keyword.size() ) ;
   MAC_CHECK_PRE( data != 0 ) ;
   MAC_CHECK_PRE( data->is_under_ownership_of(this) ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( has_entry( keyword ) )
   {
      std::string mess = keyword ;
      mess += " is already the keyword of a data in Module " ;
      mess += name() ;
      MAC_Error::object()->raise_plain( mess ) ;
   }

   MAC_String* key_string = MAC_String::create( 0, keyword ) ;
   MAC_KeywordDataPair* a_entry = MAC_KeywordDataPair::create( this,
                                                               key_string,
                                                               data ) ;
   key_string->set_owner( a_entry ) ;

   if( keyword == "type" || keyword == "concrete_name" )
   {
      ENTRIES->prepend( a_entry ) ;
   }
   else
   {
      ENTRIES->append( a_entry ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( has_entry( keyword ) ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: add_module( MAC_Module* a_module )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: add_module" ) ;
   MAC_CHECK_PRE( a_module!=0 ) ;
   MAC_CHECK_PRE( !has_module( a_module->name() ) ) ;
   MAC_CHECK_PRE( !has_entry( a_module->name() ) ) ;
   MAC_CHECK_PRE( a_module->is_under_ownership_of(this) ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   MODS->append( a_module ) ;
   a_module->FATHER = this ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( has_module( a_module->name() ) ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: merge_module( MAC_Module* a_module )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: merge_module" ) ;
   MAC_CHECK_PRE( a_module!=0 ) ;
   MAC_CHECK_PRE( a_module->name()==name() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   CTX->extend( a_module->CTX ) ;
   
   MAC_ModuleIterator* iteratorOtherMod =
         a_module->create_module_iterator( 0 ) ;
   for( iteratorOtherMod->start() ;
        iteratorOtherMod->is_valid() ;
        iteratorOtherMod->go_next() )
   {
      if( has_module( iteratorOtherMod->item()->name() ) )
      {
         MAC_Module* m = module( iteratorOtherMod->item()->name() ) ;
         m->merge_module( iteratorOtherMod->item() ) ;
      }
      else
      {
         add_module( iteratorOtherMod->item()->create_clone( this ) ) ;
      }
   }
   iteratorOtherMod->destroy() ;
   MAC_KeywordDataIterator* iteratorOtherAss =
         a_module->create_entry_iterator( 0 ) ;
   for( iteratorOtherAss->start() ;
        iteratorOtherAss->is_valid() ;
        iteratorOtherAss->go_next() )
   {
      if( !has_entry( iteratorOtherAss->item()->keyword() ) )
      {
         add_entry( iteratorOtherAss->item()->keyword(),
                    iteratorOtherAss->item()->data()->create_clone( this ) ) ;
      }
      else
      {
         replace_data_of_entry(
                    iteratorOtherAss->item()->keyword(),
                    iteratorOtherAss->item()->data()->create_clone( this ) ) ;
      }
   }
   iteratorOtherAss->destroy() ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: remove_module( std::string const& path_and_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: remove_module" ) ;
   MAC_CHECK_PRE( !path_and_name.empty() ) ;
   MAC_CHECK_PRE( has_module( path_and_name ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Module* m = module( path_and_name ) ;
   std::string dir = dirname( path_and_name ) ;
   if( dir=="" )
   {
      MAC_CHECK( MODS->has( m ) ) ;
      MODS->remove( m ) ;
   }
   else
   {
      MAC_Module* mod = 0 ;
      MAC_KeywordDataPair* assign = 0 ;
      bool ok = find( dirname( path_and_name ), mod, assign ) ;
      MAC_CHECK( ok && mod!=0 && assign==0 ) ;
      MAC_CHECK( mod->MODS->has( m ) ) ;
      mod->MODS->remove( m ) ;
   }
   if( m->FATHER!=0 )
   {
      m->CTX->extend( m->FATHER->context() ) ;
   }
   m->FATHER=0 ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( !has_module( path_and_name ) ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: remove_entry( std::string const& path_and_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: remove_entry" ) ;
   MAC_CHECK_PRE( !path_and_name.empty() ) ;
   MAC_CHECK_PRE( has_entry( path_and_name ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Module* a_mod = 0 ;
   MAC_KeywordDataPair* key_data = 0 ;
   bool ok = find( path_and_name, a_mod, key_data ) ;
   MAC_CHECK( ok && a_mod==0 && key_data!=0 ) ;
   std::string dir = dirname( path_and_name ) ;
   if( dir=="" )
   {
      MAC_CHECK( ENTRIES->has( key_data ) ) ;
      ENTRIES->remove( key_data ) ;
   }
   else
   {
      MAC_Module* mod = 0 ;
      MAC_KeywordDataPair* assign = 0 ;
      ok = find( dirname( path_and_name ), mod, assign ) ;
      MAC_CHECK( ok && mod!=0 && assign==0 ) ;
      MAC_CHECK( mod->ENTRIES->has( key_data ) ) ;
      mod->ENTRIES->remove( key_data ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( !has_entry( path_and_name ) ) ;
}




//----------------------------------------------------------------------
std::string const&
MAC_Module:: name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: name" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   std::string const& result = NAME->to_string() ;

   MAC_CHECK_POST( !result.empty() ) ;
   MAC_CHECK_POST( result.find( "/" ) >= result.size() ) ;
   MAC_CHECK_POST( result.find( " " ) >= result.size() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
std::string const&
MAC_Module:: absolute_path_name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: absolute_path_name" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   static std::string result ;
   result = "" ;
   MAC_Module const* chain = this ;
   while( chain!=0 )
   {
      if( !result.empty() ) result = "/" + result ;
      
      result = chain->name() + result ;
      chain = chain->FATHER ;
   }
   result = "/" + result ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
MAC_Module:: is_empty( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: is_empty" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   bool result = !has_module() && !has_entry() ;
   
   MAC_CHECK_POST( EQUIVALENT( result, !has_module() && !has_entry() ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
MAC_Module:: has_module( void ) const
//----------------------------------------------------------------------
{
   return( MODS->count()>0 ) ;
}




//----------------------------------------------------------------------
bool
MAC_Module:: has_module( std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: has_module" ) ;
   MAC_CHECK_PRE( !path_and_name.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   MAC_Module* a_module = 0 ;
   MAC_KeywordDataPair* assign = 0 ;
   bool result = find( path_and_name, a_module, assign ) && a_module!=0 ;
   
   MAC_CHECK_INV( invariant() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Module const*
MAC_Module:: father( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: father" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_Module const* result = FATHER ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( IMPLIES( result!=0, result->module( name() )==this ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Module*
MAC_Module:: module( std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: module" ) ;
   MAC_CHECK_PRE( !path_and_name.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_Module* result = 0 ;
   MAC_KeywordDataPair* assign = 0 ;
   bool ok = find( path_and_name, result, assign ) && result!=0 ;
   if( !ok )
   {
      std::string mess = "Can't find module " ;
      mess += path_and_name ;
      mess += " in module " ;
      mess += name() ;
      MAC_Error::object()->raise_plain( mess ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   MAC_CHECK_POST( result->name() == basename(path_and_name) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
MAC_Module:: has_entry( void ) const
//----------------------------------------------------------------------
{
   return( ENTRIES->count()>0 ) ;
}




//----------------------------------------------------------------------
bool
MAC_Module:: has_entry( std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: has_entry" ) ;
   MAC_CHECK_PRE( !path_and_name.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   MAC_Module* a_module = 0 ;
   MAC_KeywordDataPair* assign = 0 ;
   bool result = find( path_and_name, a_module, assign ) && assign!=0 ;
   
   MAC_CHECK_INV( invariant() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Data const* 
MAC_Module:: data_of_entry( std::string const& path_and_name ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: data_of_entry" ) ;
   MAC_CHECK_PRE( !path_and_name.empty() ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_Module* a_module = 0 ;
   MAC_KeywordDataPair* key_data = 0 ;
   bool ok = find( path_and_name, a_module, key_data ) && key_data!=0 ;
   if( !ok )
   {
      std::string mess = "Can't find entry " ;
      mess += path_and_name ;
      mess += " in module " ;
      mess += name() ;
      MAC_Error::object()->raise_plain( mess ) ;
   }
   MAC_CHECK( key_data->keyword() == basename(path_and_name) ) ;
   // Comment: in case of reading from a binary file, the pointer result points
   // to an object of actual type MAC_BinStored 
   MAC_Data const* result = key_data->data() ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( result!=0 ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: replace_data_of_entry( std::string const& path_and_name,
                                    MAC_Data const* new_data ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: replace_data_of_entry" ) ;
   MAC_CHECK_PRE( !path_and_name.empty() ) ;
   MAC_CHECK_PRE( new_data != 0 ) ;
   MAC_CHECK_PRE( new_data->is_under_ownership_of( this ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Module* a_mod = 0 ;
   MAC_KeywordDataPair* key_data = 0 ;
   bool ok = find( path_and_name, a_mod, key_data ) && key_data!=0 ;
   if( !ok )
   {
      std::string mess = "Can't find entry " ;
      mess += path_and_name ;
      mess += " in module " ;
      mess += name() ;
      MAC_Error::object()->raise_plain( mess ) ;
   }
   MAC_CHECK( key_data->keyword() == basename(path_and_name) ) ;
   key_data->replace_data( new_data ) ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( data_of_entry( path_and_name ) == new_data ) ;
}




//----------------------------------------------------------------------
MAC_List*
MAC_Module:: create_data_selection( MAC_Object* a_owner,
                                    std::string const& regexp,
                                    MAC_List* result,
                                    std::string const& where ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: create_data_selection" ) ;
   MAC_CHECK_PRE( !regexp.empty() ) ;
   MAC_CHECK_PRE( result==0 || result->owner()==a_owner ) ;
   
   std::string where_cp = where ;
   bool failed = substitute_variables( where_cp, 0 ) ;
   if( failed )
   {
      MAC_Error::object()->raise_plain(
         "Bad subsitution in chain "+where_cp ) ;
   }
   if(result==0)
   {
      result = MAC_List::create(a_owner) ;
      if( regexp[0]=='/' ) 
      {
         MAC_Module const* root = this ;
         while( root->FATHER!=0 ) root = root->FATHER ;
         return root->create_data_selection( a_owner, regexp, result, 
	 	where_cp ) ;
      }
   }
   
   std::string token, otherToken ;
   split(regexp, token, otherToken) ;
   
   if( token==".." && FATHER!=0 )
   {
      FATHER->create_data_selection( a_owner, otherToken, result, where_cp ) ;
   }
   else if( has_module(token) )
   {
      module(token)->create_data_selection( a_owner, otherToken, result, 
      	where_cp ) ;
   }
   else if( token=="*" && !otherToken.empty() )
   {
      MAC_ModuleIterator* iteratorOtherMod = create_module_iterator( 0 ) ;
      for( iteratorOtherMod->start() ;
           iteratorOtherMod->is_valid() ;
           iteratorOtherMod->go_next() )
      {
         iteratorOtherMod->item()->create_data_selection( a_owner,
                                                          otherToken,
                                                          result, where_cp ) ;
      }
      iteratorOtherMod->destroy() ;
   }
   else if( token=="*" && otherToken.empty() ) 
   {
      MAC_KeywordDataIterator* iteratorOtherAss =
         create_entry_iterator( 0 ) ;
      for( iteratorOtherAss->start() ;
           iteratorOtherAss->is_valid() ;
           iteratorOtherAss->go_next() )
      {
         MAC_Data* clone_data =
            static_cast<MAC_Data*>
            ( iteratorOtherAss->item()->data()->create_clone( result ) ) ;
         result->append( clone_data ) ;
      }
      iteratorOtherAss->destroy() ;
   }
   else if( token=="$key" && otherToken.empty() ) 
   {
      MAC_KeywordDataIterator* iteratorOtherAss =
         create_entry_iterator( 0 ) ;
      for( iteratorOtherAss->start() ;
           iteratorOtherAss->is_valid() ;
           iteratorOtherAss->go_next() )
      {
         result->append( MAC_String::create( result,
			iteratorOtherAss->item()->keyword() ) ) ;
      }
      iteratorOtherAss->destroy() ;
   }
   else if( has_entry(token) && otherToken.empty() )
   {
      size_t idx ;
      while( ( idx = where_cp.find( "##(" ) ) < where_cp.length() )
      {
         where_cp.replace( idx, 3, "#(" ) ;
      }
      
      MAC_Data const* eval = create_evaluation( 0, where_cp, 0 ) ;
      if( eval->to_bool() )
      {
         result->append( MAC_DataWithContext::create( result,
                                                      data_of_entry(token),
                                                      context() ) ) ;
      }
      eval->destroy() ;
   }

   MAC_CHECK_POST( result!=0 && result->owner()==a_owner ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<result->index_limit() ; i++ ),
                           dynamic_cast<MAC_Data *>(result->at(i))!=0 ) ) ;
   return( result );
   
}




//----------------------------------------------------------------------
MAC_Context const*
MAC_Module:: context( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: context" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_Context const* result = CTX ;
   if( FATHER!=0 )
   {
      TMP_CTX->re_initialize( FATHER->context(), CTX ) ;
      result = TMP_CTX ;
   }
   MAC_CHECK_POST( result!=0 ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: add_variable( MAC_Variable const* variable,
                           MAC_Data* value ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: add_variable" ) ;
   MAC_CHECK_PRE( !context()->has_variable( variable ) ) ;
   MAC_CHECK_PRE( variable->data_type()==value->data_type() ) ;
   MAC_CHECK_PRE( value->owner()==0 ) ;
   MAC_CHECK_INV( invariant() ) ;
   value->set_owner( CTX ) ;
   CTX->extend( variable, value ) ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( context()->has_variable( variable ) ) ;
   MAC_CHECK_POST( context()->value( variable ) == value ) ;
   MAC_CHECK_POST( value->is_under_ownership_of( this ) ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: modify_variable( MAC_Variable const* variable,
                              MAC_Data* value ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: modify_variable" ) ;
   MAC_CHECK_PRE( context()->has_variable( variable ) ) ;
   MAC_CHECK_PRE( variable->data_type()==value->data_type() ) ;
   MAC_CHECK_PRE( value->owner()==0 ) ;
   MAC_CHECK_INV( invariant() ) ;
   value->set_owner( CTX ) ;
   if( CTX->has_variable( variable ) )
   {
      CTX->set_value_of( variable, value ) ;
   }
   else
   {
      CTX->extend( variable, value ) ;
   }
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( context()->has_variable( variable ) ) ;
   MAC_CHECK_POST( context()->value( variable ) == value ) ;
   MAC_CHECK_POST( value->is_under_ownership_of( this ) ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: write( std::string const& file,
                    std::string const& format ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: write" ) ;
   MAC_CHECK_PRE( !file.empty() ) ;
   MAC_CHECK_PRE( format=="text" || format=="hybrid" ) ;   
   MAC_CHECK_INV( invariant() ) ;
   
   std::ofstream out( file.c_str(), std::ios::out | std::ios::app ) ;
   if( !out )
   {
      std::string mess = "MAC_Module \"" ;
      mess += name() ;
      mess += "\" writing failure : \n   Unable to open file \"" ;
      mess += file ;
      mess += "\"" ;
      MAC_Error::object()->raise_plain( mess ) ;
   }
   bool hybrid = format=="hybrid" ;
   std::string const binary_file = file + ".bin" ;
   if( hybrid &&
       !MAC_BinStored::is_valid_binary_file( binary_file ) )
   {
      MAC_BinStored::init_binary_file( binary_file ) ;
   }
   recursive_print( out, 0, hybrid, binary_file ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: print( std::ostream& os, size_t indent_width ) const 
//----------------------------------------------------------------------
{	
   MAC_LABEL( "MAC_Module:: print" ) ;
   MAC_CHECK_INV( invariant() ) ;
   recursive_print( os, (int)indent_width, false, "" ) ;
}




//-------------------------------------------------------------------------
void
MAC_Module:: display_info( std::ostream& os, size_t indent_width ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: display_info" ) ;
   MAC_Object::display_info( os, indent_width ) ;
   std::string const s( indent_width, ' ' ) ;
   os << s << "module name : " << NAME->to_string() << std::endl ;
}




//----------------------------------------------------------------------
std::string const&
MAC_Module:: current_parsed_module_path_name( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: current_parsed_module_path_name" ) ;
   static std::string result = MAC_current_parsed_module_path_name() ;
   return( result ) ;
}




//---------------------------------------------------------------------------
std::string
MAC_Module:: data_as_string( std::string const& path_and_name,
                             MAC_Context const* ct,
                             bool& failed ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: data_as_string" ) ;
   MAC_CHECK_PRE( !path_and_name.empty() ) ;
   
   if( ct==0 ) ct=context() ;
   
   std::string result ;
   
   size_t idx = path_and_name.find("/") ;
   
   if( idx < path_and_name.size() ) 
   {
      MAC_Module const* root = 0 ;
      std::string path = path_and_name.substr(idx+1) ; ;
      
      if( idx==0 ) // Absolute path
      {
         root = this ;
         
         while( root->FATHER!=0 ) root=root->FATHER ;
      }
      else if( idx==2 && path_and_name[0]=='.' && path_and_name[1]=='.' ) // ../path
      {
         root = FATHER ;
      }
      else 
      {
         std::string sub_mod = path_and_name.substr(0,idx) ;
         
         if( has_module( sub_mod ) ) 
         {
            root = module( sub_mod ) ;
         }
      }
      if( root!=0 ) 
      {         
         result = root->data_as_string( path, ct, failed ) ;
      }
      
   }
   else
   {
      failed = !has_entry( path_and_name ) ;
      if( !failed ) 
      {         
         MAC_Data const* data = data_of_entry( path_and_name ) ;
         failed = !data->value_can_be_evaluated( ct ) ;
         
         std::ostringstream str ;
         if( !failed )
         {
            result = data->value_as_string( ct ) ;
         }
      }
   }

   MAC_CHECK_POST( IMPLIES( !failed, !result.empty() ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_ModuleIterator*
MAC_Module:: create_module_iterator( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: create_module_iterator" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   MAC_ModuleIterator* result =
                           MAC_ModuleIterator::create( a_owner, MODS ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_KeywordDataIterator*
MAC_Module:: create_entry_iterator( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: create_entry_iterator" ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_KeywordDataIterator* result =
                  MAC_KeywordDataIterator::create( a_owner, ENTRIES ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
std::string
MAC_Module:: basename( std::string const& path_and_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: basename" ) ;
   const char separator = '/' ;
   size_t idx = path_and_name.find_last_of( separator ) ;
   std::string result ;
   if( idx<path_and_name.length() )
   {
      result = path_and_name.substr( idx+1, path_and_name.length()-idx-1 ) ;
   }
   else
   {
      result = path_and_name ;
   }
   MAC_CHECK_POST( result.find( separator ) >= result.size() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
std::string
MAC_Module:: dirname( std::string const& path_and_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: dirname" ) ;
   const char separator = '/' ;
   size_t idx = path_and_name.find_last_of( separator ) ;
   std::string result ;
   if( idx>0 && idx<path_and_name.length() )
   {
      result = path_and_name.substr( 0, idx ) ;
   }
   else
   {
      result = "" ;
   }
   return( result ) ;
}




//---------------------------------------------------------------------------
bool
MAC_Module:: substitute_variables( std::string& replaced,
                                   MAC_Context const* ct ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: substitute_variables" ) ;
   MAC_CHECK_PRE( !replaced.empty() ) ;
   bool failed = false ;

   MAC_Context const* ctx = ( ct==0 ? context() :
                              MAC_ContextPair::create(0,context(),ct) ) ;
   // Adressage indirect des sous-modules
   size_t start = 0 ;
   size_t varidx = replaced.find("#(", start ) ;
   bool new_notation =  varidx < replaced.length() ;
   
   if( new_notation )
   {
      while( !failed &&
             varidx< replaced.length() )
      {
         if( varidx==0 || replaced[varidx-1] != '#' )
         {
            size_t endvaridx = replaced.find(")",varidx+2) ;
            if( endvaridx < replaced.length())
            {
               size_t i1 = varidx+2 ;
               size_t i2 = endvaridx-varidx-2 ;
               std::string default_value ;
            
               std::string varname =  replaced.substr( i1, i2 ) ;
               size_t i3 = varname.find( "," ) ;
               if( i3 < varname.length() ) 
               {
                  default_value = varname.substr( i3+1, varname.length()-i3 ) ;
                  varname = varname.substr( 0, i3 ) ;
               }

               std::string val = data_as_string(varname,ctx,failed) ;
               
               if( val.empty() )
               {
                  if( !default_value.empty() )
                  {
                     failed = false ;                  
                     replaced.replace( varidx, endvaridx-varidx+1, 
		     	default_value ) ;
                  }
               
               }
               else
               {            
                  std::string valuestr = " " + val + " " ;
            
                  replaced.replace( varidx, endvaridx-varidx+1, valuestr ) ;
               }
            }
            else
            {
               MAC_Error::object()->raise_plain(
                  "Bad syntax in variable substitution in chain "+replaced ) ;
            }
            start = 0 ;
         }
         else
         {
            start = varidx+1 ; 
         }
         varidx = replaced.find("#(", start ) ;
      
      }
   }
   else
   {
      
      MAC_KeywordDataIterator* it = create_entry_iterator( 0 ) ;
   
      for( it->start() ; !failed && it->is_valid() ; it->go_next() )
      {
         std::string prefix ;
         MAC_KeywordDataPair const* pair = it->item() ;
         std::string a_keyword = pair->keyword() ;

         std::string valuestr ;
      
         size_t idx=0 ;
      
         while( !failed && (idx=replaced.find(a_keyword,idx))<replaced.length())
         {
            bool to_be_replaced = true ;
            if( idx > 0 )
            {
               char c = replaced[idx-1] ;
               to_be_replaced = !( ( ('a'<=c) && (c<='z') ) ||
                                   ( ('A'<=c) && (c<='Z') ) ||
                                   ( ('0'<=c) && (c<='9') ) ||
                                   ( c=='_' ) ||
                                   ( c=='"' ) ||
                                   ( c=='#' )
                  ) ; 
               
            }
            if( idx > 1 )
            {
               to_be_replaced = to_be_replaced &&
                  !( replaced[idx-2]=='#' && replaced[idx-1]=='(' ) ;
               
            }
            if( to_be_replaced && idx+a_keyword.length()<replaced.length() )
            {
               char c = replaced[idx+a_keyword.length()] ;
               to_be_replaced = !( ( ('a'<=c) && (c<='z') ) ||
                                   ( ('A'<=c) && (c<='Z') ) ||
                                   ( ('0'<=c) && (c<='9') ) ||
                                   ( c=='_' ) ||
                                   ( c=='"' ) ||
                                   ( c=='#' )
                  ) ; 
            }
         
            if( to_be_replaced )
            {
               if( valuestr.empty() )
               {
                  valuestr = " " ;
                  valuestr += data_as_string(a_keyword,ctx,failed) ;
                  valuestr += " " ;
               }
            
               if( !failed )
               {
                  replaced.replace( idx, a_keyword.length(), valuestr ) ;
               }
            }
            else idx++ ;
      
         }
         
      }
      it->destroy() ;
   }
   
   if( ct!=0 ) ctx->destroy() ;
   
   return( failed ) ;
}




//---------------------------------------------------------------------------
MAC_DataWithContext const*
MAC_Module:: create_evaluation( MAC_Object * a_owner,
                                std::string const& expression,
                                MAC_Context const* ct ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: create_evaluation" ) ;
   MAC_CHECK_PRE( !expression.empty() ) ;
   
   std::ostringstream full_str ;
   full_str << "MODULE Root" << std::endl ;   
   std::string replaced = expression ;
   bool failed = substitute_variables( replaced, ct ) ;
   

   MAC_DataWithContext * result = 0 ;
   if( !failed )
   {
      full_str << "_expression_result = ( " << replaced << " )" << std::endl ;
      full_str << "END MODULE Root" << std::endl ;
   
      std::istringstream is( full_str.str() ) ;
   
      MAC_Module * built = MAC_Module::create( 0, "Test", is ) ;
      MAC_Module const* root = built->module( "Root" ) ;
      result = MAC_DataWithContext::create(
                     0,
                     root->data_of_entry( "_expression_result" ),
                     root->context() ) ;
      built->set_owner( result ) ;
      
      if( !result->value_can_be_evaluated( ct ) )
      {
         result->destroy() ; result = 0 ;
      }
      else if( a_owner != 0 )
      {
         result->set_owner( a_owner ) ;
      }
   }
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( IMPLIES( result != 0, result->owner() == a_owner ) ) ;
   MAC_CHECK_POST( IMPLIES( result != 0, 
   	result->value_can_be_evaluated( ct ) ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------------
bool
MAC_Module:: comparable( MAC_Object const* other ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: comparable" ) ;
   return( MAC_Object::comparable( other ) ||
                    dynamic_cast<MAC_String const*>( other ) != 0 ) ;
}




//----------------------------------------------------------------------------
bool
MAC_Module:: is_equal( MAC_Object const* other ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: is_equal" ) ;
   // less restrictive than MAC_Object::is_equal_PRE
   MAC_CHECK_PRE( is_equal_PRE( other ) ) ;

   MAC_Module const* lex = dynamic_cast<MAC_Module const* >( other ) ;
   bool result ;
   
   if( lex!=0 )
   {
      result = NAME->is_equal( lex->NAME ) ;
   }
   else
   {
      MAC_String const* ss = dynamic_cast<MAC_String const* >( other ) ;
      MAC_ASSERT( ss != 0 ) ;
      result = NAME->is_equal( ss ) ;
   }
   MAC_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------------
int
MAC_Module:: three_way_comparison( MAC_Object const* other ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: three_way_comparison" ) ;
   MAC_CHECK_PRE( three_way_comparison_PRE( other ) ) ;
   MAC_Module const* lex = dynamic_cast<MAC_Module const* >( other ) ;
   int result ;
   
   if( lex!=0 )
   {
      result = NAME->three_way_comparison( lex->NAME ) ;
   }
   else
   {
      MAC_String const* ss = dynamic_cast<MAC_String const* >( other ) ;
      MAC_ASSERT( ss != 0 ) ;
      result = NAME->three_way_comparison( ss ) ;
   }
   MAC_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
size_t
MAC_Module:: hash_code( void ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: hash_code" ) ;
   return( NAME->hash_code() ) ;
}




//----------------------------------------------------------------------
MAC_Module*
MAC_Module:: create_copy( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: create_copy" ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Module* result = new MAC_Module( a_owner, name(), CTX ) ;

   MAC_ModuleIterator* iteratorOtherMod = create_module_iterator( 0 ) ;
   for( iteratorOtherMod->start() ;
        iteratorOtherMod->is_valid() ;
        iteratorOtherMod->go_next() )
   {
      MAC_Module* clone_mod =
                        iteratorOtherMod->item()->create_copy( result ) ;
      result->add_module( clone_mod ) ;
   }
   iteratorOtherMod->destroy() ;
   MAC_KeywordDataIterator* iteratorOtherAss = create_entry_iterator( 0 ) ;
   for( iteratorOtherAss->start() ;
        iteratorOtherAss->is_valid() ;
        iteratorOtherAss->go_next() )
   {
      MAC_Data* clone_data =
         static_cast<MAC_Data*>
            ( iteratorOtherAss->item()->data()->create_clone( result ) ) ;
      result->add_entry( iteratorOtherAss->item()->keyword(),
                         clone_data ) ;
   }
   iteratorOtherAss->destroy() ;   

   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   MAC_CHECK_POST( result->name() == name() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: complete_context( MAC_Module* root,
                               MAC_Module* dup,
                               MAC_Context const* ref_ctx ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: complete_context" ) ;
   MAC_CHECK( root != 0 ) ;
   MAC_CHECK( dup != 0 ) ;
   MAC_CHECK( ref_ctx != 0 ) ;

   bool has_modif = false ;
   MAC_KeywordDataIterator* assIt = dup->create_entry_iterator( 0 ) ;
   for( assIt->start() ; assIt->is_valid() ; assIt->go_next() )
   {
      MAC_KeywordDataPair const* ass = assIt->item() ;
      MAC_Data const* dat = ass->data() ;
      complete_context( root, dat, ref_ctx, has_modif ) ;
      
   }
   assIt->destroy() ; assIt = 0 ;

   MAC_Context const* root_ctx = root->context() ;
   for( size_t i=0 ; i<root_ctx->nb_variables() ; i++ )
   {
      MAC_Variable const* var = root_ctx->variable(i) ;
      MAC_Data const* dat = root_ctx->value(var) ;
      complete_context( root, dat, ref_ctx, has_modif ) ;
   }
   
   MAC_Context const* dup_ctx = dup->context() ;
   for( size_t i=0 ; i<dup_ctx->nb_variables() ; i++ )
   {
      MAC_Variable const* var = dup_ctx->variable(i) ;
      MAC_Data const* dat = dup_ctx->value(var) ;
      complete_context( root, dat, ref_ctx, has_modif ) ;
   }
   
   if( has_modif )
   {
      complete_context( root, dup, ref_ctx ) ;
   }
   else
   {
      MAC_ModuleIterator* modIt = dup->create_module_iterator( 0 ) ;
      for( modIt->start() ; modIt->is_valid() ; modIt->go_next() )
      {
         MAC_Module* a_module = modIt->item() ;
         complete_context( root, a_module, ref_ctx ) ;
      }
      modIt->destroy() ; modIt = 0 ;
   }
}




//----------------------------------------------------------------------
void
MAC_Module:: complete_context( MAC_Module* root,
                               MAC_Data const* dat,
                               MAC_Context const* ref_ctx,
                               bool & has_modif ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: complete_context" ) ;
   MAC_CHECK( root != 0 ) ;
   MAC_CHECK( dat != 0 ) ;
   MAC_CHECK( ref_ctx != 0 ) ;

   MAC_Context const* ctx = root->context() ;
   if( !dat->context_has_required_variables(ctx) )
   { 
      MAC_List* needed = MAC_List::create( 0 ) ;
      dat->declare( needed ) ;
      MAC_ListIterator* it = needed->create_iterator( needed ) ;
      for( it->start() ;
           it->is_valid() ;
           it->go_next() )
      {
         MAC_Variable const* var =
            static_cast<MAC_Variable const*>( it->item() ) ;
         if( ref_ctx->has_variable( var ) && ! ctx->has_variable( var ) )
         {
            has_modif = true ;
            root->add_variable( var, 
	    	ref_ctx->value( var )->create_clone( 0 ) ) ;
         }   
      }
      needed->destroy() ; needed = 0 ;
   }
}




//----------------------------------------------------------------------
bool
MAC_Module:: find( std::string const& nom,
                   MAC_Module*& theModule,
                   MAC_KeywordDataPair*& theAssignment ) const
//----------------------------------------------------------------------
{
   MAC_CHECK( theModule==0 && theAssignment==0 ) ;
   MAC_CHECK_INV( invariant() ) ;
   bool result = false ;
   
   // The name is null : we have found the name
   if( nom.empty() )
   {
      theModule = const_cast<MAC_Module*>( this ) ;
      result = true ;
   }
   else
   {
      std::string token, otherToken ;
      split(nom, token, otherToken) ;
      
      // The name begins with "/" : we must find the child
      // root = ( nom[0]=='/' ) ;
      MAC_String * str = MAC_String::create( 0, token ) ;
      MAC_Object* obj = MODS->item( str ) ;
      if( obj!=0 )
      {
         // Find a module
         MAC_Module* a_module = static_cast<MAC_Module* >(obj) ;
         if( otherToken.empty() )
         {
            theModule = a_module ;
            result = true ;
         }
         else
         {
            result = a_module->find( otherToken, theModule, theAssignment ) ;
         }
      }
      else
      {
         // Searching for an assignment
         obj = ENTRIES->item( str ) ;
         theAssignment = static_cast<MAC_KeywordDataPair* >(obj) ;
         MAC_CHECK( IMPLIES( obj!=0, 
	 	dynamic_cast<MAC_KeywordDataPair* >(obj) != 0 ) ) ;
         
         result = theAssignment!=0 && otherToken.empty() ;
      }
      str->destroy() ;

   }
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( !result || theModule!=0 || theAssignment!=0 ) ;
   MAC_CHECK_POST( !( theModule!=0 && theAssignment!=0 ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC_Module:: split( std::string const& nom,
                    std::string& token,
                    std::string& otherToken ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: split" ) ;
   
   // si nom=="/toto/titi/tutu" ou nom=="toto/titi/tutu" alors 
   //     token="toto" et otherToken = "/titi/tutu"
   // si nom=="/toto" ou nom = "toto" alors
   //     token="toto" et otherToken=""
   int n = (int)nom.length() ;
      
   int start = (int)nom.find_first_not_of( "/" ) ;
   if( ( start >=0 ) && ( start < n ) )
   {
      int stop = (int)nom.find_first_of( "/", start ) ;
      if( ( stop >= 0 ) && ( stop < n ) ) 
      {
         token      = nom.substr( start, stop-start ) ;
         otherToken = nom.substr( stop, n ) ;
      }
      else
      {
         token = nom.substr( start, n-start ) ;
         otherToken = "" ;
      }
   }
   else
   {
      token = nom ;
      otherToken = "" ;
   }
}




//----------------------------------------------------------------------
std::ostream& 
MAC_Module:: recursive_print( std::ostream& s,
                              int n,
                              bool hybrid,
                              std::string const& bin_file) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Module:: recursive_print" ) ;
   MAC_CHECK( IMPLIES( hybrid, !bin_file.empty() ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   std::string bl( n, ' ' ) ;
   
   s << bl << "MODULE " << NAME->to_string()  << std::endl ;

   // print contexts:
   {
      MAC_Context const* father_ctx = 0 ;
      if( FATHER!=0 )
      {
         father_ctx = FATHER->context() ;
      }
      CTX->print( s, n+2, father_ctx ) ;
   }
   
   // print entries:
   {
      MAC_KeywordDataIterator* assIt = create_entry_iterator( 0 ) ;
      for( assIt->start() ; assIt->is_valid() ; assIt->go_next() )
      {
         MAC_KeywordDataPair const* ass = assIt->item() ;
      
         s << bl << "  " << ass->keyword() << " = " ;
         if( !hybrid ||
             !MAC_BinStored::is_type_supported( ass->data()->data_type() ) ||
             !ass->data()->is_constant() )
         {
            ass->data()->print( s, 0 ) ;
         }
         else
         {
            MAC_Data const* to_print =
               MAC_BinStored::create_reference( 0, ass->data(),
                                                bin_file, true ) ;
            to_print->print( s, 0 ) ;
            to_print->destroy() ;
         }
         s << std::endl ;
      }
      assIt->destroy() ;
   }

   // print modules:
   {
      MAC_ModuleIterator* modIt = create_module_iterator( 0 ) ;
      for( modIt->start() ; modIt->is_valid() ; modIt->go_next() )
      {
         MAC_Module const* a_module = modIt->item() ;
         a_module->recursive_print( s, n+2, hybrid, bin_file ) ;
      }
      modIt->destroy() ;
   }

   
   s << bl << "END MODULE " << name() << std::endl ;
   MAC_CHECK_INV( invariant() ) ;		
   return(s);
}




//----------------------------------------------------------------------
bool
MAC_Module:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   MAC_ASSERT( NAME!=0 ) ;
   MAC_ASSERT( MODS!=0 ) ;
   MAC_ASSERT( ENTRIES!=0 ) ;
   // To assume that FATHER hasn't been destroy before self
   MAC_ASSERT( FATHER==0 || FATHER->MODS!=0 ) ;
   return( true ) ;
}
