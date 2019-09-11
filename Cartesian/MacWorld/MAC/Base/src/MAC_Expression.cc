#include <MAC_Expression.hh>

#include <MAC_assertions.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_Error.hh>
#include <MAC_Iterator.hh>
#include <MAC_List.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_Root.hh>
#include <MAC_Sequence.hh>

#include <stringVector.hh>

#include <iostream>
#include <sstream>
#include <set>


//----------------------------------------------------------------------
MAC_Expression*
MAC_Expression:: create( MAC_Object* a_owner,
                         std::string const& a_name,
                         MAC_Sequence const* argument_list,
                         std::string const& a_comment ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: create" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   MAC_CHECK_PRE( argument_list!=0 ) ;
   MAC_CHECK_PRE( 
      FORALL( ( size_t i=0 ; i<argument_list->index_limit() ; ++i ),
         dynamic_cast<MAC_Data*>( argument_list->at(i) ) != 0 ) ) ;
   MAC_CHECK_PRE( valid_arguments_of( a_name, argument_list ) ) ;

   MAC_Expression const* proto =
      static_cast<MAC_Expression const*>(
                                    plugins_map()->item( a_name ) ) ;
   MAC_ASSERT( proto->is_a_prototype() ) ;
   
   MAC_Expression* result = proto->create_replica( a_owner, argument_list ) ;
   result->COMMENT = a_comment ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->name() == a_name ) ;
   MAC_CHECK_POST( result->external_brackets_are_set() ) ;
   MAC_CHECK_POST( result->comment() == a_comment ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
stringVector const&
MAC_Expression:: registered_expressions( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: registered_expressions" ) ;

   return( plugins_name() ) ;
}




//----------------------------------------------------------------------
MAC_Expression:: MAC_Expression( MAC_Object* a_owner,
                                 std::string const& a_name,
                                 MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Data( a_owner )
   , NAME( a_name )
   , HAS_BRACKETS( true )
   , COMMENT( "" )
   , ARGUMENTS( argument_list )
   , IT( argument_list->create_iterator( this ) )
{
   MAC_LABEL( "MAC_Expression:: MAC_Expression" ) ;
   MAC_CHECK_PRE( argument_list != 0 ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( owner() == a_owner ) ;
   MAC_CHECK_POST( name() == a_name ) ;
   MAC_CHECK_POST( !is_a_prototype() ) ;
}




//----------------------------------------------------------------------
MAC_Expression:: MAC_Expression( std::string const& a_name )
//----------------------------------------------------------------------
   : MAC_Data( plugins_map() )
   , NAME( a_name )
   , HAS_BRACKETS( true )
   , COMMENT( "" )
   , ARGUMENTS( 0 )
   , IT( 0 )
{
   MAC_LABEL( "MAC_Expression:: MAC_Expression" ) ;
   
   plugins_map()->register_item( a_name, this ) ;
   plugins_name().append( a_name ) ;
   
   MAC_CHECK_POST( is_a_prototype() ) ;
}




//----------------------------------------------------------------------
MAC_Expression*
MAC_Expression:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: create_clone" ) ;
   MAC_Expression* result = 0 ;
   if( is_a_prototype() )
   {
      result = const_cast<MAC_Expression*>(this) ;
   }
   else
   {
      MAC_Sequence const* args = ARGUMENTS ;
      MAC_Sequence* l = args->create_clone( 0 ) ;
      l->clear() ;
      for( size_t i=0 ; i<args->index_limit() ; ++i )
      {
         MAC_Object* obj = args->at(i) ;
         l->append( obj->create_clone( l ) ) ;
      }
      result = MAC_Expression::create( a_owner, name(), l, comment() ) ;
      l->set_owner( result ) ;
      if( HAS_BRACKETS )
         result->set_external_brackets() ;
      else
         result->unset_external_brackets() ;
   }
   MAC_CHECK_POST( IMPLIES( is_a_prototype(), result==this ) ) ;
   MAC_CHECK_POST( IMPLIES( !is_a_prototype(),
                            create_clone_POST( result, a_owner ) ) ) ;
   MAC_CHECK_POST( result->name()==name() ) ;
   MAC_CHECK_POST( result->comment()==comment() ) ;
   MAC_CHECK_POST( result->external_brackets_are_set()==external_brackets_are_set() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Expression:: ~MAC_Expression( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: ~MAC_Expression" ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
std::string const&
MAC_Expression:: name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: name" ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( NAME ) ;
}




//----------------------------------------------------------------------
std::string const&
MAC_Expression:: usage_of( std::string const& a_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: usage_of" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   
   MAC_Expression const* proto =
      static_cast<MAC_Expression const*>(
                                    plugins_map()->item( a_name ) ) ;
   MAC_ASSERT( proto->is_a_prototype() ) ;

   return( proto->usage() ) ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: valid_arguments_of( std::string const& a_name,
                                     MAC_Sequence const* argument_list )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: valid_arguments_of" ) ;
   MAC_CHECK_PRE( !a_name.empty() ) ;
   MAC_CHECK_PRE( argument_list!=0 ) ;
   MAC_CHECK_PRE( 
      FORALL( ( size_t i=0 ; i<argument_list->index_limit() ; ++i ),
         dynamic_cast<MAC_Data*>( argument_list->at(i) ) != 0 ) ) ;
   
   MAC_Expression const* proto =
      static_cast<MAC_Expression const*>(
                                    plugins_map()->item( a_name ) ) ;
   MAC_ASSERT( proto->is_a_prototype() ) ;

   return( proto->valid_arguments( argument_list ) ) ;
}




//----------------------------------------------------------------------
void
MAC_Expression:: declare( MAC_List * lst ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: declare" ) ;
   MAC_CHECK_PRE( declare_PRE( lst ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      MAC_Data const* data = static_cast<MAC_Data const*>( IT->item() ) ;
      data->declare( lst ) ;
      if( !IT->is_valid() ) raise_error( "circular definition" ) ; // Cycle
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( declare_POST( lst ) ) ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: context_has_required_variables( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: context_has_required_variables" ) ;
   MAC_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   bool result = true ;
   for( IT->start() ; result && IT->is_valid() ; IT->go_next() )
   {
      MAC_Data const* data = static_cast<MAC_Data const*>( IT->item() ) ;
      result = data->context_has_required_variables(ct) ;
      if( !IT->is_valid() ) raise_error( "circular definition" ) ; // Cycle
   }

   MAC_CHECK_POST( 
      FORALL( ( size_t i=0 ; i<nb_arguments() ; ++i ),
              !result || arg(i)->context_has_required_variables(ct) ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
size_t
MAC_Expression:: nb_arguments( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: nb_arguments" ) ;
   MAC_CHECK( !is_a_prototype() ) ;
   return( ARGUMENTS->count() ) ;
}




//----------------------------------------------------------------------
MAC_Data const*
MAC_Expression:: arg( size_t idx ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: arguments" ) ;
   MAC_CHECK( !is_a_prototype() ) ;
   MAC_CHECK( idx<nb_arguments() ) ;

   MAC_Data const* result = extract_arg( ARGUMENTS, idx) ;

   MAC_CHECK_POST( result != 0 ) ;   
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Data const*
MAC_Expression:: extract_arg( MAC_Sequence const* some_arguments,
                              size_t idx )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: extract_arg" ) ;
   MAC_CHECK( idx<some_arguments->count() ) ;
   MAC_CHECK_PRE( dynamic_cast<MAC_Data*>( some_arguments->at(idx) ) != 0 ) ;

   MAC_Data const* result = static_cast<MAC_Data*>( some_arguments->at(idx) ) ;

   MAC_CHECK_POST( result != 0 ) ;   
   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC_Expression:: raise_error( std::string const& message ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: raise_error" ) ;
   MAC_CHECK_PRE( !message.empty() ) ;

   std::ostringstream msg ;
   msg << "*** Expression " << name() << " error:" << std::endl ;
   msg << "    usage: " << usage() << std::endl ;
   msg << "    " << message << std::endl ;
   if( !COMMENT.empty() ) msg << "    " << COMMENT << std::endl ;
   MAC_Error::object()->raise_plain( msg.str() ) ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: value_can_be_evaluated( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: value_can_be_evaluated " ) ;
   MAC_CHECK_PRE( !is_a_prototype() ) ;
   bool result = true ;
   for( IT->start() ; IT->is_valid() && result ; IT->go_next() )
   {
      MAC_Data const* data = static_cast<MAC_Data const*>( IT->item() ) ;
      result = data->value_can_be_evaluated( ct ) ;
      if( !IT->is_valid() )
      {
         result = false ; // Cycle
         break ;
      }
   }
   return result ;
}




//----------------------------------------------------------------------
stringVector const&
MAC_Expression:: undefined_variables( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: undefined_variables" ) ;
   MAC_CHECK_PRE( !is_a_prototype() ) ;
   
   static stringVector result(0) ;
   stringVector undef_var(0) ;
   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      MAC_Data const* data = static_cast<MAC_Data const*>( IT->item() ) ;
      stringVector const& s = data->undefined_variables( ct ) ;
      for( size_t i=0 ; i<s.size() ; ++i )
      {
         undef_var.extend( s(i) ) ;
      }
      if( !IT->is_valid() )
      {
         break ; // Cycle
      }
   }
   result = undef_var ;
   result.sort() ;
   return result ;
}




//----------------------------------------------------------------------
MAC_Data*
MAC_Expression:: create_derivative( MAC_Object* a_owner,
                                    MAC_Variable const* var,
                                    MAC_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   MAC_Error::object()->raise_plain(
      "No create_derivative method implemented for expression "+name() ) ;
   MAC_Expression* result = 0 ;
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: external_brackets_are_set( void ) const
//----------------------------------------------------------------------
{
   return( HAS_BRACKETS ) ;
}




//----------------------------------------------------------------------
void
MAC_Expression:: set_external_brackets( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: set_external_brackets" ) ;
   HAS_BRACKETS = true ;
   MAC_CHECK_POST( external_brackets_are_set() ) ;
}




//----------------------------------------------------------------------
void
MAC_Expression:: unset_external_brackets( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: unset_external_brackets" ) ;
   HAS_BRACKETS = false ;
   MAC_CHECK_POST( !external_brackets_are_set() ) ;
}




//----------------------------------------------------------------------
void
MAC_Expression:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << name() << "(" ;
   if( ARGUMENTS!=0 )
   {
      bool prem = true ;
      for( IT->start() ; IT->is_valid() ; IT->go_next() )
      {
         if( !prem )
         {
            os << ", " ;
         }
         prem = false ;
         MAC_Data const* data = static_cast<MAC_Data const*>( IT->item() ) ;
         data->print( os, 0 ) ;
         if( !IT->is_valid() ) raise_error( "circular definition" ) ; // Cycle
      }
   }
   else
   {
      os << "prototype" ;
   }
   os << ")" ;
}




//----------------------------------------------------------------------
void
MAC_Expression:: print_prototypes( std::ostream& os, size_t indent_width )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: print_prototypes" ) ;

   std::string space( indent_width, ' ' ) ;
   MAC_Iterator* it = plugins_map()->create_iterator( 0 ) ;
   std::set< std::string > names ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      MAC_Expression const* expr = 
                     static_cast< MAC_Expression* >( it->item() ) ;
      names.insert( expr->name() ) ;
   }
   std::set< std::string >::const_iterator itn = names.begin() ;
   for( ; itn != names.end() ; ++itn ) 
   {
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         MAC_Expression const* expr = 
                     static_cast< MAC_Expression* >( it->item() ) ;
         if( expr->name() == (*itn) )
         {
            os << "     (\"" << expr->usage() << "\" \"\" \"\")" << std::endl ;
            break ;
         }
      }
   }
   it->destroy() ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( ARGUMENTS==0 ) ;  
}




//----------------------------------------------------------------------
std::string const&
MAC_Expression:: comment( void ) const
//----------------------------------------------------------------------
{
   return( COMMENT ) ;  
}




//----------------------------------------------------------------------
MAC_Data*
MAC_Expression:: create_non_const_simplification( MAC_Object* a_owner ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: create_non_const_simplification" ) ;
   MAC_CHECK_PRE( create_non_const_simplification_PRE( a_owner ) ) ;

   MAC_List* new_list = MAC_List::create( 0 ) ;
   for( IT->start(); IT->is_valid() ; IT->go_next() )
   {
      MAC_Data const* data = static_cast<MAC_Data const*>( IT->item() ) ;
      new_list->append( data->create_simplification( new_list ) ) ;
      if( !IT->is_valid() ) raise_error( "circular definition" ) ; // Cycle
   }
   MAC_Expression* tmp = create_replica( 0, new_list ) ;
   new_list->set_owner( tmp ) ;
   
   MAC_Data* result =0 ;
   
   if( tmp->is_constant() )
   {
      result = tmp->create_simplification( a_owner ) ;
   }
   else
   {
      result = tmp->create_operator_simplification( a_owner ) ;
   }
   if( result!=tmp )
   {
      tmp->destroy() ;
   }
   else
   {
      tmp->set_owner( a_owner ) ;
   }
   
   MAC_CHECK_POST( create_non_const_simplification_POST( a_owner, result ) ) ;
   return result ;
}




//----------------------------------------------------------------------
MAC_Data*
MAC_Expression:: create_operator_simplification( MAC_Object* a_owner )  
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: create_operator_simplification" ) ;
   MAC_Data* result = this ;
   MAC_CHECK_POST( create_operator_simplification_POST( a_owner, result ) ) ;
   return result ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: is_raw_data( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Expression:: is_raw_data" ) ;

   MAC_CHECK_POST( is_raw_data_POST( false ) ) ;
   
   return false ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: create_replica_PRE( MAC_Object const* a_owner,
                                     MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( argument_list != 0 ) ;
   MAC_ASSERT( argument_list->count() == argument_list->index_limit() ) ;
   MAC_ASSERT(
      FORALL( ( size_t i=0 ; i<argument_list->count() ; ++i ),
              dynamic_cast<MAC_Data const*>( argument_list->at(i) ) != 0 ) ) ;
   MAC_ASSERT( valid_arguments( argument_list ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: create_replica_POST( MAC_Expression const* result,
                                      MAC_Object const* a_owner,
                                      MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   MAC_ASSERT( result->name() == name() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: valid_arguments_PRE( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( some_arguments != 0 ) ;
   MAC_ASSERT( some_arguments->count() == some_arguments->index_limit() ) ;
   MAC_ASSERT(
      FORALL( ( size_t i=0 ; i<some_arguments->count() ; ++i ),
              dynamic_cast<MAC_Data const*>( some_arguments->at(i) ) != 0 ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: create_operator_simplification_POST(
               MAC_Object const* a_owner, MAC_Data const* result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result==this || result->owner()==a_owner ) ;
   return true ;
}




//----------------------------------------------------------------------
bool
MAC_Expression:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Data::invariant() ) ;
   return true ;  
}




//----------------------------------------------------------------------
stringVector&
MAC_Expression:: plugins_name( void )
//----------------------------------------------------------------------
{
   static stringVector result(0) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_ObjectRegister*
MAC_Expression:: plugins_map( void )
//----------------------------------------------------------------------
{
   static MAC_ObjectRegister* result =
             MAC_ObjectRegister::create( MAC_Root::object(),
                                         "MAC_Expression descendant" ) ;
   return( result ) ;
}
