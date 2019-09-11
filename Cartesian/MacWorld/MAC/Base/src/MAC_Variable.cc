#include <MAC_Variable.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Double.hh>
#include <MAC_Context.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_Root.hh>
#include <MAC_List.hh>
#include <stringVector.hh>

#include <iostream>

//----------------------------------------------------------------------
MAC_Variable:: MAC_Variable( MAC_Object* a_owner,
                             std::string const& a_name ) 
//----------------------------------------------------------------------
   : MAC_Data( a_owner )
   , NAME( a_name )
   , KIND( MAC_Variable::data_type( a_name ) )
   , EVALUATING( false )
{
   variable_list()->append( this ) ;
   ID = variable_list()->count()-1 ;
  
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_Variable*
MAC_Variable:: create_clone( MAC_Object * a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: create_clone" ) ;
   MAC_Variable* result = new MAC_Variable( a_owner, this ) ;
   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return result ;
}

//----------------------------------------------------------------------
MAC_Variable:: MAC_Variable( MAC_Object* a_owner,
                             MAC_Variable const* other ) 
//----------------------------------------------------------------------
      : MAC_Data( a_owner )
      , NAME( other->NAME )
      , KIND( other->KIND )
      , ID( other->ID )
      , EVALUATING( false )
{
   MAC_LABEL( "MAC_Variable:: MAC_Variable" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_Variable:: ~MAC_Variable( void ) 
//----------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
size_t
MAC_Variable:: nb_objects( void )
//----------------------------------------------------------------------
{
   return variable_list()->index_limit() ;
}

//----------------------------------------------------------------------
MAC_Variable const *
MAC_Variable:: object( std::string const& a_name ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: object" ) ;

   // Check name :
   MAC_Variable:: data_type( a_name ) ;
   
   MAC_Variable const* result = 0 ;
   MAC_Iterator* it = variable_list()->create_iterator(0) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      MAC_Variable const* var = static_cast<MAC_Variable const*>(it->item()) ;
      if( var->name()==a_name )
      {
         result=var ;
      }
   }
   it->destroy() ; it = 0 ;
   if( result==0 )
   {
      result = new MAC_Variable( variable_list(), a_name ) ;
   }
   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->name() == a_name ) ;
   return result ;
}

//----------------------------------------------------------------------
MAC_Variable const*
MAC_Variable:: object( size_t id )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: object" ) ;
   MAC_CHECK_PRE( id < nb_objects() ) ;
   
   MAC_Variable const* result=
             static_cast<MAC_Variable const*>( variable_list()->at( id ) ) ;

   MAC_CHECK_POST( result->id_number() == id ) ;
   return result ;
}

//----------------------------------------------------------------------
size_t
MAC_Variable:: id_number( void ) const
//----------------------------------------------------------------------
{
   return ID ;
}

//----------------------------------------------------------------------
std::string const&
MAC_Variable:: name( void ) const
//----------------------------------------------------------------------
{
   return NAME ;
}

//----------------------------------------------------------------------
void
MAC_Variable:: declare( MAC_List* lst ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: declare" ) ;
   MAC_CHECK_PRE( declare_PRE( lst ) ) ;
   
   lst->extend( const_cast<MAC_Variable*>(this) ) ;
   
   MAC_CHECK_POST( declare_POST( lst ) ) ;
}

//----------------------------------------------------------------------
bool
MAC_Variable:: context_has_required_variables( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: context_has_required_variables" ) ;
   MAC_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;

   bool result = ct->has_variable( this ) ;

   MAC_CHECK_POST( EQUIVALENT( result, ct->has_variable(this) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_Variable:: data_type( std::string const& a_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: data_type(a_name)" ) ;
   
   MAC_Data::Type result = MAC_Data::Undefined ;

   if( a_name.length() >= 2 )
   {
      std::string typ = a_name.substr( 0, 2 ) ;
      if( typ=="IS" )
      {
         result = Int ;
      }
      else if( typ=="IV" )
      {
         result = IntVector ;
      }
      else if( typ=="IA" )
      {
         result = IntArray2D ;
      }
      else if( typ=="BA" )
      {
         result = BoolArray2D ;
      }
      else if( typ=="SA" )
      {
         result = StringArray2D ;
      }
      else if( typ=="DS" )
      {
         result = Double ;
      }
      else if( typ=="DV" )
      {
         result = DoubleVector ;
      }
      else if( typ=="DA" )
      {
         result = DoubleArray2D ;
      }
      else if( typ=="BS" )
      {
         result = Bool ;
      }
      else if( typ=="BV" )
      {
         result = BoolVector ;
      }
      else if( typ=="SS" )
      {
         result = String ;
      }
      else if( typ=="SV" )
      {
         result = StringVector ;
      }
   }
   if( result == MAC_Data::Undefined )
   {
      std::string msg = "\""+a_name+"\" is not a valid variable name\n" ;
      msg += "A valid name is \"XY_name\"\n" ;
      msg += "   where \"X\" is the scalar type of the variable :\n" ;
      msg += "       - \"I\" : integer\n" ;
      msg += "       - \"D\" : double\n" ;
      msg += "       - \"B\" : boolean\n" ;
      msg += "       - \"S\" : string\n" ;
      msg += "   and \"Y\" defined its dimension :\n" ;
      msg += "       - \"S\" : simple (only one element)\n" ;
      msg += "       - \"V\" : vector\n" ;
      msg += "       - \"A\" : array2D\n" ;
      msg += "Examples : \"DV_coordinates\", \"SS_name\", \"IA_connectivity\"\n" ;
      MAC_Error::object()->raise_plain( msg ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_Variable:: data_type( void ) const
//----------------------------------------------------------------------
{
   return KIND ;
}

//----------------------------------------------------------------------
bool
MAC_Variable:: value_can_be_evaluated( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: value_can_be_evaluated" ) ;

   bool result = ! EVALUATING ;
   if( result )
   {
      EVALUATING = true ;
      result =  ( ct!=0 &&
                  ct->has_variable(this) &&
                  ct->value(this)!=0 &&
                  ct->value(this)->value_can_be_evaluated(ct) ) ;
      EVALUATING = false ;
   }

   MAC_CHECK_POST(
      IMPLIES(
         result,
         ct!=0 &&
            ct->has_variable(this) &&
            ct->value(this) != 0 &&
            ct->value(this)->value_can_be_evaluated(ct) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
stringVector const&
MAC_Variable::  undefined_variables( MAC_Context const* ct ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: undefined_variables" ) ;

   static stringVector result(0) ;
   result.re_initialize(0) ;
   if( EVALUATING ||
       ct == 0 || !ct->has_variable(this) || ct->value(this) == 0 )
   {
      result.append( name() ) ;
   }
   else
   {
      EVALUATING = true ;
      result = ct->value(this)->undefined_variables( ct ) ;
      EVALUATING = false ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_Variable:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_bool" ) ;
   MAC_CHECK_PRE( to_bool_PRE( ct ) ) ;
   return data( ct )->to_bool( ct ) ;
}

//----------------------------------------------------------------------
double
MAC_Variable:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE( ct ) ) ;
   return data( ct )->to_double( ct ) ;
}

//----------------------------------------------------------------------
int
MAC_Variable:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE( ct ) ) ;   
   return data( ct )->to_int( ct ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_Variable:: to_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_string" ) ;
   MAC_CHECK_PRE( to_string_PRE( ct ) ) ;
   return data( ct )->to_string( ct ) ;
}

//----------------------------------------------------------------------
doubleVector const& 
MAC_Variable:: to_double_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_double_vector" ) ;
   MAC_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   return data( ct )->to_double_vector( ct ) ;
}

//----------------------------------------------------------------------
intVector const&
MAC_Variable:: to_int_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_int_vector" ) ;
   MAC_CHECK_PRE( to_int_vector_PRE( ct ) ) ;
   return data( ct )->to_int_vector( ct ) ;
}

//----------------------------------------------------------------------
stringVector const&
MAC_Variable:: to_string_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_string_vector" ) ;
   MAC_CHECK_PRE( to_string_vector_PRE( ct ) ) ;
   return data( ct )->to_string_vector(ct ) ;
}

//----------------------------------------------------------------------
boolVector const&
MAC_Variable:: to_bool_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_bool_vector" ) ;
   MAC_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;
   return data( ct )->to_bool_vector( ct ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
MAC_Variable:: to_double_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_double_array2D" ) ;
   MAC_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   return data( ct )->to_double_array2D( ct ) ;
}

//----------------------------------------------------------------------
boolArray2D const&
MAC_Variable:: to_bool_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_bool_array2D" ) ;
   MAC_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   return data( ct )->to_bool_array2D( ct ) ;
}

//----------------------------------------------------------------------
stringArray2D const&
MAC_Variable:: to_string_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_string_array2D" ) ;
   MAC_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;
   return data( ct )->to_string_array2D( ct ) ;
}

//----------------------------------------------------------------------
intArray2D const&
MAC_Variable:: to_int_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_int_array2D" ) ;
   MAC_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;
   return data( ct )->to_int_array2D( ct ) ;
}

//----------------------------------------------------------------------
doubleArray3D const&
MAC_Variable:: to_double_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_double_array3D" ) ;
   MAC_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   return data( ct )->to_double_array3D( ct ) ;
}

//----------------------------------------------------------------------
intArray3D const&
MAC_Variable:: to_int_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: to_int_array3D" ) ;
   MAC_CHECK_PRE( to_int_array3D_PRE( ct ) ) ;
   return data( ct )->to_int_array3D( ct ) ;
}

//----------------------------------------------------------------------
MAC_Data const*
MAC_Variable:: data( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_CHECK( ct!=0 && ct->has_variable(this) ) ;

   MAC_Data const* result = ct->value( this ) ;

   MAC_CHECK( result!=0 ) ;   
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_Variable:: is_raw_data( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: is_raw_data" ) ;

   MAC_CHECK_POST( is_raw_data_POST( false ) ) ;
   
   return false ;
}

//----------------------------------------------------------------------
void
MAC_Variable:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << "$" << NAME ;
}

//----------------------------------------------------------------------
MAC_Data*
MAC_Variable:: create_derivative( MAC_Object* a_owner,
                                  MAC_Variable const* var,
                                  MAC_Context const* ct ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Variable:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   MAC_Data* result =
      ( var->name()==name() ?
        MAC_Double::create( a_owner, 1.0 ) :
        ct->value(this)->create_derivative( a_owner, var, ct ) ) ;
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_List*
MAC_Variable:: variable_list( void ) 
//----------------------------------------------------------------------
{
   static MAC_List* result = MAC_List::create( MAC_Root::object() ) ;
   return result ;
}
