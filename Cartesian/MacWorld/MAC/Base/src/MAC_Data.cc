#include <MAC_Data.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_Bool.hh>
#include <MAC_BoolArray2D.hh>
#include <MAC_Error.hh>
#include <MAC_Int.hh>
#include <MAC_IntVector.hh>
#include <MAC_IntArray2D.hh>
#include <MAC_IntArray3D.hh>
#include <MAC_List.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_Double.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_DoubleArray2D.hh>
#include <MAC_DoubleArray3D.hh>
#include <MAC_String.hh>
#include <MAC_StringArray2D.hh>
#include <MAC_StringVector.hh>
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

#include <sstream>
#include <iostream>



//----------------------------------------------------------------------
MAC_Data:: MAC_Data( MAC_Object* a_owner ) 
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
{
   MAC_LABEL( "MAC_Data:: MAC_Data" ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
MAC_Data:: ~MAC_Data( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: ~MAC_Data" ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
void
MAC_Data:: declare( MAC_List* lst ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: declare" ) ;
   MAC_CHECK_PRE( declare_PRE( lst ) ) ;
   MAC_CHECK_POST( declare_POST( lst ) ) ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: is_constant( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: is_constant" ) ;
   MAC_List* lst = MAC_List::create( 0 ) ;
   declare( lst ) ;
   bool result = lst->count()==0 ;
   lst->destroy() ;
   
   return result ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: is_raw_data( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: is_raw_data" ) ;

   MAC_CHECK_POST( is_raw_data_POST( true ) ) ;
   
   return true ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: is_raw_data_POST( bool result ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( IMPLIES( result, is_constant() ) ) ;
   return true ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: context_has_required_variables( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: context_has_required_variables" ) ;
   MAC_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
std::string
MAC_Data:: type_name( Type kind ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: type_name" ) ;
   std::string result ;
   switch( kind )
   {
      case Double :
         result = "Double" ;
         break ;
      case Int :
         result = "Int" ;
         break ;
      case Bool :
         result = "Bool" ;
         break ;
      case String :
         result = "String" ;
         break ;
      case DoubleVector :
         result = "DoubleVector" ;
         break ;
      case IntVector :
         result = "IntVector" ;
         break ;
      case BoolVector :
         result = "BoolVector" ;
         break ;
      case DoubleArray2D :
         result = "DoubleArray2D" ;
         break ;
      case BoolArray2D :
         result = "BoolArray2D" ;
         break ;
      case StringArray2D :
         result = "StringArray2D" ;
         break ;
      case DoubleArray3D :
         result = "DoubleArray3D" ;
         break ;
      case IntArray2D :
         result = "IntArray2D" ;
         break ;
      case IntArray3D :
         result = "IntArray3D" ;
         break ;
      case StringVector :
         result = "StringVector" ;
         break ;
      default :
         MAC_Error::object()->raise_plain( "Bad kind" ) ;
         break ;
   }
   return result ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: value_can_be_evaluated( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   return true ;
}




//----------------------------------------------------------------------
stringVector const&
MAC_Data:: undefined_variables( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: undefined_variables" ) ;
   static stringVector result(0) ;
   return result ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_bool" ) ;
   MAC_CHECK_PRE( to_bool_PRE( ct ) ) ;
   exitWithError( "to_bool" ) ;
   return 0 ;
}




//----------------------------------------------------------------------
double
MAC_Data:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE( ct ) ) ;
   exitWithError( "to_double" ) ;
   return 0 ;
}




//----------------------------------------------------------------------
int
MAC_Data:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE( ct ) ) ;
   exitWithError( "to_int" ) ;
   return 0 ;
}




//----------------------------------------------------------------------
std::string
const&
MAC_Data:: to_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_string" ) ;
   MAC_CHECK_PRE( to_string_PRE( ct ) ) ;
   exitWithError( "to_string" ) ;
   static std::string dummy ;
   return dummy ;
}




//----------------------------------------------------------------------
doubleVector
const& 
MAC_Data:: to_double_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_double_vector" ) ;
   MAC_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   exitWithError( "to_double_vector" ) ;
   static doubleVector ret(0) ;
   return ret ;
}




//----------------------------------------------------------------------
intVector
const&
MAC_Data:: to_int_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_int_vector" ) ;
   MAC_CHECK_PRE( to_int_vector_PRE( ct ) ) ;
   exitWithError( "to_int_vector" ) ;
   static intVector ret(0) ;
   return ret ;
}




//----------------------------------------------------------------------
stringVector const& 
MAC_Data:: to_string_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_string_vector" ) ;
   MAC_CHECK_PRE( to_string_vector_PRE( ct ) ) ;
   exitWithError( "to_string_vector" ) ;
   static stringVector ret( 0 ) ;
   return ret ;     
}




//----------------------------------------------------------------------
boolVector const&
MAC_Data:: to_bool_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_bool_vector" ) ;
   MAC_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;
   exitWithError( "to_bool_vector" ) ;
   static boolVector ret( 0 ) ;
   return ret ;  
}




//----------------------------------------------------------------------
doubleArray2D const&
MAC_Data:: to_double_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_double_array2D" ) ;
   MAC_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   exitWithError( "to_double_array2D" ) ;
   static doubleArray2D ret(0,0) ;  
   return ret ;
}




//----------------------------------------------------------------------
boolArray2D const&
MAC_Data:: to_bool_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_bool_array2D" ) ;
   MAC_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   exitWithError( "to_bool_array2D" ) ;
   static boolArray2D ret(0,0) ;  
   return ret ;
}




//----------------------------------------------------------------------
stringArray2D const&
MAC_Data:: to_string_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_string_array2D" ) ;
   MAC_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;
   exitWithError( "to_string_array2D" ) ;
   static stringArray2D ret(0,0) ;  
   return ret ;
}




//----------------------------------------------------------------------
intArray2D const&
MAC_Data:: to_int_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_int_array2D" ) ;
   MAC_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;
   exitWithError( "to_int_array2D" ) ;
   static intArray2D ret(0,0) ;  
   return( ret ) ;
}




//----------------------------------------------------------------------
doubleArray3D const&
MAC_Data:: to_double_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_double_array3D" ) ;
   MAC_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   exitWithError( "to_double_array3D" ) ;
   static doubleArray3D ret(0,0,0) ;  
   return ret ;
}




//----------------------------------------------------------------------
intArray3D const&
MAC_Data:: to_int_array3D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: to_int_array3D" ) ;
   MAC_CHECK_PRE( to_int_array3D_PRE( ct ) ) ;
   exitWithError( "to_int_array3D" ) ;
   static intArray3D ret(0,0,0) ;  
   return( ret ) ;
}




//----------------------------------------------------------------------
std::string
MAC_Data:: value_as_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: value_as_string" ) ;
   MAC_CHECK_PRE( value_can_be_evaluated( ct ) ) ;
   
   std::ostringstream strout ;
   strout.precision(16) ;
   
   if( data_type()==String ) 
   {
      strout << "\"" << to_string( ct ) << "\"" ;
   }
   else if( data_type()==Bool ) 
   {
      strout <<  ( to_bool(ct) ? "true" : "false" ) ;
   }
   else if( data_type()==Int ) 
   {
      strout << to_int(ct) ;
   }
   else if( data_type()==Double )
   {
      MAC::print_double(strout, to_double(ct)) ;
   }
   else if( data_type()==StringVector )
   {
      stringVector const& val = to_string_vector(ct) ;
      strout << "< " ;
      for( size_t i=0 ; i<val.size() ; ++i )
      {
         strout << "\"" << val(i) << "\" "  ;
      }
      strout << ">" ;
   }
   else if( data_type()==BoolVector )
   {
      boolVector const& val = to_bool_vector(ct) ;
      strout << "< " ;
      for( size_t i=0 ; i<val.size() ; ++i )
      {
         strout << ( val(i) ? "true" : "false" ) <<  " "  ;
      }
      strout << ">" ;
   }
   else if( data_type()==IntVector )
   {
      intVector const& val = to_int_vector(ct) ;
      strout << "< " ;
      for( size_t i=0 ; i<val.size() ; ++i )
      {
         strout << val(i) << " "  ;
      }
      strout << ">" ;
   }
   else if( data_type()==DoubleVector )
   {
      doubleVector const& val = to_double_vector(ct) ;
      strout << "< " ;
      for( size_t i=0 ; i<val.size() ; ++i )
      {
         MAC::print_double(strout, val(i) ) ;
         strout << " "  ;
      }
      strout << ">" ;
   }
   else if( data_type()==BoolArray2D )
   {
      boolArray2D const& val = to_bool_array2D(ct) ;
      strout << "[ " ;
      for( size_t i=0 ; i<val.index_bound(0) ; ++i )
      {
         strout << "< " ;
         for( size_t j=0 ; j<val.index_bound(1) ; ++j )
         {
            strout << ( val(i,j) ? "true" : "false" ) ;
            strout << " " ;
         }
         strout << ">" ;
         if( i!=val.index_bound(0)-1 ) strout << "," ;
      }
      strout << " ]" ;
   }
   else if( data_type()==IntArray2D )
   {
      intArray2D const& val = to_int_array2D(ct) ;
      strout << "[ " ;
      for( size_t i=0 ; i<val.index_bound(0) ; ++i )
      {
         strout << "< " ;
         for( size_t j=0 ; j<val.index_bound(1) ; ++j )
         {
            strout << val(i,j) ;
            strout << " " ;
         }
         strout << ">" ;
         if( i!=val.index_bound(0)-1 ) strout << "," ;
      }
      strout << " ]" ;
   }
   else if( data_type()==StringArray2D )
   {
      stringArray2D const& val = to_string_array2D(ct) ;
      strout << "[ " ;
      for( size_t i=0 ; i<val.index_bound(0) ; ++i )
      {
         strout << "< " ;
         for( size_t j=0 ; j<val.index_bound(1) ; ++j )
         {
            strout << "\"" << val(i,j) << "\"" ;
            strout << " " ;
         }
         strout << ">" ;
         if( i!=val.index_bound(0)-1 ) strout << "," ;
      }
      strout << " ]" ;
   }
   else if( data_type()==DoubleArray2D )
   {
      doubleArray2D const& val = to_double_array2D(ct) ;
      strout << "[ " ;
      for( size_t i=0 ; i<val.index_bound(0) ; ++i )
      {
         strout << "< " ;
         for( size_t j=0 ; j<val.index_bound(1) ; ++j )
         {
            MAC::print_double(strout, val(i,j) ) ;
            strout << " " ;
         }
         strout << ">" ;
         if( i!=val.index_bound(0)-1 ) strout << "," ;
      }
      strout << " ]" ;
   }
   else
   {
      MAC_Error::object()->raise_internal(
         "*** MAC_Data:: value_as_string error\n"
         "    unexpected data type \""+type_name( data_type() )+"\"" ) ;
   }      
   return( strout.str() ) ;
}




//----------------------------------------------------------------------
MAC_Data*
MAC_Data:: create_derivative( MAC_Object* a_owner,
                              MAC_Variable const* var,
                              MAC_Context const* ct ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   MAC_Error::object()->raise_plain(
      "No create_derivative method implemented for expression of class "
      + MAC_Object::type_name() ) ;
   MAC_Data* result = 0 ;
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: create_derivative_PRE( MAC_Object* a_owner,
                                  MAC_Variable const* var,
                                  MAC_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   MAC_ASSERT( var!=0 ) ;
   MAC_ASSERT( data_type()==Double || data_type()==DoubleVector ) ;
   MAC_ASSERT( var->data_type()==Double ) ;
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: create_derivative_POST( MAC_Object* a_owner,
                                   MAC_Variable const* var,
                                   MAC_Data const* result ) const 
//----------------------------------------------------------------------
{
   MAC_ASSERT( result!=0 ) ;
   MAC_ASSERT( result->owner()==a_owner ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
MAC_Data*
MAC_Data:: create_simplification( MAC_Object* a_owner ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: create_simplification" ) ;
   MAC_CHECK_PRE( data_type()!=Undefined ) ;

   MAC_Data* result = 0 ;
   
   if( is_constant() )
   {
      switch( data_type() ) 
      {
         case Int :
            result = MAC_Int::create( a_owner, to_int() ) ;
            break ;
         case IntVector :
            result = MAC_IntVector::create( a_owner, to_int_vector() ) ;
            break ;
         case IntArray2D :
            result = MAC_IntArray2D::create( a_owner, to_int_array2D() ) ;
            break ;
         case BoolArray2D :
            result = MAC_BoolArray2D::create( a_owner, to_bool_array2D() ) ;
            break ;
         case StringArray2D :
            result = MAC_StringArray2D::create( a_owner, to_string_array2D() ) ;
            break ;
         case IntArray3D :
            result = MAC_IntArray3D::create( a_owner, to_int_array3D() ) ;
            break ;
         case Double :
            result = MAC_Double::create( a_owner, to_double() ) ;
            break ;
         case DoubleVector :
            result = MAC_DoubleVector::create( a_owner, to_double_vector() ) ;
            break ;
         case DoubleArray2D :
            result = MAC_DoubleArray2D::create( a_owner,to_double_array2D() ) ;
            break ;
         case DoubleArray3D :
            result = MAC_DoubleArray3D::create( a_owner,to_double_array3D() ) ;
            break ;
         case String :
            result = MAC_String::create( a_owner, to_string() ) ;
            break ;
         case StringVector :
            result = MAC_StringVector::create( a_owner, to_string_vector() ) ;
            break ;
         case Bool :
            result = MAC_Bool::create( a_owner, to_bool() ) ;
            break ;
        default :
            MAC_Error::object()->raise_internal( "Bad type" ) ;
      }
   }
   else
   {
      result = create_non_const_simplification( a_owner ) ;
   }
   
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_Data*
MAC_Data:: create_non_const_simplification( MAC_Object* a_owner ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Data:: create_non_const_simplification" ) ;
   MAC_CHECK_PRE( create_non_const_simplification_PRE( a_owner ) ) ;

   MAC_Data* result = create_clone( a_owner ) ;

   MAC_CHECK_POST( create_non_const_simplification_POST( a_owner, result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: declare_PRE( MAC_List const* lst ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( lst!=0 ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: declare_POST( MAC_List const* lst ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( FORALL( ( size_t i=0 ; i<lst->count(); i++),
                       dynamic_cast<MAC_Variable*>(lst->at(i))!=0 ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: context_has_required_variables_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( ct != 0 ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_double_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==Double ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_int_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==Int ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_bool_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==Bool ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_string_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==String ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_double_vector_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==DoubleVector ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_int_vector_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==IntVector ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_bool_vector_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==BoolVector ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_string_vector_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==StringVector ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_double_array2D_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==DoubleArray2D ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_bool_array2D_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==BoolArray2D ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_string_array2D_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==StringArray2D ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_int_array2D_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==IntArray2D ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_double_array3D_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==DoubleArray3D ) ;
   return true ;  
}




//----------------------------------------------------------------------
bool
MAC_Data:: to_int_array3D_PRE( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( value_can_be_evaluated( ct ) ) ;
   MAC_ASSERT( data_type()==IntArray3D ) ;
   return true ;  
}




//----------------------------------------------------------------------
void
MAC_Data:: exitWithError( std::string const& mess ) const
//----------------------------------------------------------------------
{
   MAC::out() << "Current lexical object : " << std::endl ;
   display_info( MAC::out(), 5 ) ;
   MAC::out() << "can't be converted through " << mess << " method ! " ;
   MAC_Error::object()->raise_plain(
      "Error in MAC_Data derived class in method "+mess ) ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: create_non_const_simplification_PRE( MAC_Object* a_owner ) const 
//----------------------------------------------------------------------
{
   MAC_ASSERT( data_type()!=Undefined ) ;
   MAC_ASSERT( !is_constant() ) ;
   
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: create_non_const_simplification_POST( MAC_Object* a_owner,
                                                 MAC_Data const* result ) const 
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
MAC_Data:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   return true ;  
}

