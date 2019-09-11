#include <MAC_DataWithContext.hh>

#include <MAC_assertions.hh>
#include <MAC_Context.hh>
#include <MAC_ContextPair.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_DataWithContextExp.hh>
#include <MAC_List.hh>
#include <MAC_Sequence.hh>
#include <MAC_Variable.hh>

#include <iostream>

//----------------------------------------------------------------------------
MAC_DataWithContext*
MAC_DataWithContext:: create( MAC_Object* a_owner,
                              MAC_Data const* data,
                              MAC_Context const* ct )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: create" ) ;
   MAC_CHECK_PRE( data != 0 ) ;
   MAC_CHECK_PRE( ct != 0 ) ;

   MAC_DataWithContext* result = new MAC_DataWithContext( a_owner, data, ct ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->data_type() == data->data_type() ) ;
   MAC_CHECK_POST( result->is_raw_data() == data->is_raw_data() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_DataWithContext:: MAC_DataWithContext( MAC_Object* a_owner,
                                           MAC_Data const* data,
                                           MAC_Context const* ct )
//----------------------------------------------------------------------------
   : MAC_Data( a_owner )
   , DATA( data )
   , CTX( ct->create_clone( this ) )
   , TMP_CTX( MAC_ContextPair::create( this,0,0 ) )
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
MAC_DataWithContext*
MAC_DataWithContext:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: create_clone" ) ;

   MAC_DataWithContext* result = new MAC_DataWithContext( a_owner, DATA, CTX ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->data_type() == data_type() ) ;
   MAC_CHECK_POST( result->is_constant() == is_constant() ) ;
   MAC_CHECK_POST( result->is_raw_data() == is_raw_data() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DataWithContext:: ~MAC_DataWithContext( void ) 
//----------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
MAC_DataWithContext:: declare( MAC_List* lst ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: declare" ) ;
   MAC_CHECK_PRE( declare_PRE( lst ) ) ;

   MAC_List* ctx = MAC_List::create( 0 ) ;
   DATA->declare( ctx ) ;
   for( size_t i=0 ; i<ctx->count() ; i++ )
   {
      MAC_Variable* var = static_cast<MAC_Variable*>( ctx->at(i) ) ;
      if( !CTX->has_variable( var ) )
      {
         lst->extend( var ) ;
      }
      else
      {
         MAC_DataWithContext* data =
            MAC_DataWithContext::create( 0,
                                         CTX->value( var ),
                                         CTX ) ;
         data->declare( lst ) ;
         data->destroy() ; data = 0 ;
      }
   }
   ctx->destroy() ; ctx=0 ;
   
   MAC_CHECK_POST( declare_POST( lst ) ) ;
}

//----------------------------------------------------------------------------
MAC_Data::Type
MAC_DataWithContext:: data_type( void ) const
//----------------------------------------------------------------------------
{
   return( DATA->data_type() ) ;
}

//----------------------------------------------------------------------
MAC_Context const*
MAC_DataWithContext:: context( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: context" ) ;
   MAC_Context const* result = CTX ;
   if( ct!=0 && ct!=CTX )
   {
      TMP_CTX->re_initialize( CTX, ct ) ;
      result = TMP_CTX ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_DataWithContext:: value_can_be_evaluated( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: value_can_be_evaluated" ) ;
   return( DATA->value_can_be_evaluated( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
stringVector const&
MAC_DataWithContext:: undefined_variables( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: undefined_variables" ) ;
   return( DATA->undefined_variables( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
bool
MAC_DataWithContext:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_bool" ) ;
   MAC_CHECK_PRE( to_bool_PRE( ct ) ) ;
   return( DATA->to_bool( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
double
MAC_DataWithContext:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE( ct ) ) ;
   return( DATA->to_double( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
int
MAC_DataWithContext:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE( ct ) ) ;
   return( DATA->to_int( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_DataWithContext:: to_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_string" ) ;
   MAC_CHECK_PRE( to_string_PRE( ct ) ) ;
   return( DATA->to_string( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
doubleVector const& 
MAC_DataWithContext:: to_double_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_double_vector" ) ;
   MAC_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   return( DATA->to_double_vector( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
intVector const&
MAC_DataWithContext:: to_int_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_int_vector" ) ;
   MAC_CHECK_PRE( to_int_vector_PRE( ct ) ) ;
   return( DATA->to_int_vector( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
stringVector const& 
MAC_DataWithContext:: to_string_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_string_vector" ) ;
   MAC_CHECK_PRE( to_string_vector_PRE( ct ) ) ;
   return( DATA->to_string_vector( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
boolVector const&
MAC_DataWithContext:: to_bool_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_bool_vector" ) ;
   MAC_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;
   return( DATA->to_bool_vector( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
MAC_DataWithContext:: to_double_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_double_array2D" ) ;
   MAC_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   return( DATA->to_double_array2D( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
boolArray2D const&
MAC_DataWithContext:: to_bool_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_bool_array2D" ) ;
   MAC_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   return( DATA->to_bool_array2D( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
stringArray2D const&
MAC_DataWithContext:: to_string_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_string_array2D" ) ;
   MAC_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;
   return( DATA->to_string_array2D( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
intArray2D const&
MAC_DataWithContext:: to_int_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_int_array2D" ) ;
   MAC_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;
   return( DATA->to_int_array2D( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
doubleArray3D const&
MAC_DataWithContext:: to_double_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_double_array3D" ) ;
   MAC_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   return( DATA->to_double_array3D( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
intArray3D const&
MAC_DataWithContext:: to_int_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: to_int_array3D" ) ;
   MAC_CHECK_PRE( to_int_array3D_PRE( ct ) ) ;
   return( DATA->to_int_array3D(context( ct ) ) ) ;
}

//----------------------------------------------------------------------
bool
MAC_DataWithContext:: is_raw_data( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: is_raw_data" ) ;
   bool result = DATA->is_raw_data() ;
   MAC_CHECK_POST( is_raw_data_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_DataWithContext:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DataWithContext:: print" ) ;

   if( CTX->nb_variables() != 0 )
   {
      MAC_DataWithContextExp const* exp =
              MAC_DataWithContextExp::create( 0, DATA, CTX ) ;
      exp->print( os, indent_width ) ;
      exp->destroy() ; exp = 0 ;
   }
   else
   {
      DATA->print( os, indent_width ) ;
   }
}
