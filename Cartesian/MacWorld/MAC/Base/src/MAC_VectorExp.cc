#include <MAC_VectorExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Bool.hh>
#include <MAC_ContextPair.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_Double.hh>
#include <MAC_Error.hh>
#include <MAC_Int.hh>
#include <MAC_Iterator.hh>
#include <MAC_List.hh>
#include <MAC_Sequence.hh>
#include <MAC_String.hh>
#include <MAC_Variable.hh>

#include <iostream>

MAC_VectorExp const*
MAC_VectorExp::PROTOTYPE_NVECTOR = new MAC_VectorExp(
                                      MAC_VectorExp::nvector, "nvector") ;

MAC_VectorExp const*
MAC_VectorExp::PROTOTYPE_CONDITIONAL_VECTOR = new MAC_VectorExp(
                         MAC_VectorExp::cond_vect, "conditional_vector") ;
MAC_VectorExp const*
MAC_VectorExp::PROTOTYPE_VECTOR = new MAC_VectorExp(
                                          MAC_VectorExp::vect, "vector") ;
MAC_VectorExp const*
MAC_VectorExp::PROTOTYPE_SIZE = new MAC_VectorExp(
                                            MAC_VectorExp::size, "size") ;
MAC_VectorExp const*
MAC_VectorExp::PROTOTYPE_COMPO = new MAC_VectorExp(
                                  MAC_VectorExp::component, "component") ;

MAC_VectorExp const*
MAC_VectorExp::PROTOTYPE_INCREASING = new MAC_VectorExp(
                                MAC_VectorExp::increasing, "increasing") ;

MAC_VectorExp const*
MAC_VectorExp::PROTOTYPE_GREATER = new MAC_VectorExp(
                                      MAC_VectorExp::greater, "greater") ;

MAC_VectorExp const*
MAC_VectorExp::PROTOTYPE_APPLY = new MAC_VectorExp(
                                          MAC_VectorExp::apply, "apply") ;

MAC_VectorExp const*
MAC_VectorExp::PROTOTYPE_REVERSE = new MAC_VectorExp(
                                          MAC_VectorExp::reverse, "reverse") ;

MAC_VectorExp const*
MAC_VectorExp::PROTOTYPE_SUM = new MAC_VectorExp(
                                          MAC_VectorExp::sum, "sum") ;

//----------------------------------------------------------------------
MAC_VectorExp:: MAC_VectorExp( VectorExp exp_id, std::string const& a_name ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( exp_id )
   , MY_IT( 0 )
   , RESULT_D( 0 )
   , RESULT_I( 0 )
   , RESULT_B( 0 )
   , RESULT_S( 0 )
{
   MAC_LABEL( "MAC_VectorExp:: MAC_VectorExp" ) ;
}

//----------------------------------------------------------------------
MAC_VectorExp:: MAC_VectorExp( MAC_Object* a_owner,
                               VectorExp exp_id,
                               std::string const& a_name,
                               MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP( exp_id )
   , MY_IT( argument_list->create_iterator( this ) )
   , RESULT_D( 0 )
   , RESULT_I( 0 )
   , RESULT_B( 0 )
   , RESULT_S( 0 )
{
   MAC_LABEL( "MAC_VectorExp:: MAC_VectorExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_VectorExp:: ~MAC_VectorExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: ~MAC_VectorExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
   if( is_a_prototype() )
   {
      if( OP == vect )
      {
         PROTOTYPE_VECTOR = 0 ;
      }
      if( OP == size )
      {
         PROTOTYPE_SIZE = 0 ;
      }
      if( OP == component )
      {
         PROTOTYPE_COMPO = 0 ;
      }
      if( OP == increasing )
      {
         PROTOTYPE_INCREASING = 0 ;
      }
      if( OP == reverse )
      {
         PROTOTYPE_REVERSE = 0 ;
      }
      if( OP == sum )
      {
         PROTOTYPE_SUM = 0 ;
      }
       if( OP == greater )
      {
         PROTOTYPE_GREATER = 0 ;
      }
      if( OP == apply )
      {
         PROTOTYPE_APPLY = 0 ;
      }
   }
}

//----------------------------------------------------------------------
MAC_VectorExp*
MAC_VectorExp:: create_replica( MAC_Object* a_owner,
                                MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_VectorExp* result =
                   new MAC_VectorExp( a_owner, OP, name(), argument_list ) ;

   MAC_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_VectorExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "undefined" ;
   switch( OP )
   {
      case nvector :
         result = "nvector(IS,<scalar type>)" ;
         break ;
      case vect :
         result = "vector(<list of values with the same type>)" ;
         break ;
      case cond_vect :
         result = "conditional_vector(<list of pairs: BS,value with all values of the same type>)" ;
         break ;
      case size :
         result = "size(DV|IV|BV|SV)" ;
         break ;
      case component :
         result = "component(DV|IV|BV|SV,IS)" ;
         break ;
      case increasing :
         result = "increasing(DV|IV)" ;
         break ;
      case reverse :
         result = "reverse(DV|IV|SV|BV)" ;
         break ;
      case sum :
         result = "sum(DV|IV)" ;
         break ;
      case greater :
         result = "greater(DV|IV,DS|IS)" ;
         break ;
      case apply :
         result = "apply(source DV|IV|BV|SV,applied function DS|IS|BS|SS,variable name SS,[component name SS])" ;
         break ;
      default :
         MAC_Error::object()->raise_internal( "not implemented" ) ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_VectorExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = true ;
   switch( OP )
   {
      case vect :
         result = some_arguments->count()>0 ;
         if( result )
         {
            Type k = Undefined ;
            MAC_Iterator* iter = some_arguments->create_iterator( 0 ) ;
            for( iter->start() ; result && iter->is_valid() ; iter->go_next() )
            {
               MAC_Data const* data =
                  static_cast<MAC_Data const*>( iter->item() ) ;
               if( k == Undefined )
               {
                  k = data->data_type() ;
               }
               else
               {
                  result &= ( k==data->data_type() ) ;
               }
            }
            iter->destroy() ; iter = 0 ;
         }
         break ;
      case cond_vect :
         result = some_arguments->count()>0
            && some_arguments->count() % 2 == 0;
         if( result )
         {
            Type k = Undefined ;
            MAC_Iterator* iter = some_arguments->create_iterator( 0 ) ;
            for( iter->start() ; result && iter->is_valid() ; iter->go_next() )
            {
               MAC_Data const* a_bool =
                  static_cast<MAC_Data const*>( iter->item() ) ;
               iter->go_next() ;
               MAC_Data const* data =
                  static_cast<MAC_Data const*>( iter->item() ) ;
               if(  k == Undefined )
               {
                  k = data->data_type() ;
               }
               else
               {
                  result &= ( a_bool->data_type()==MAC_Data::Bool &&
                              ( k==data->data_type() ) ) ;
               }
            }
            iter->destroy() ; iter = 0 ;
         }
         break ;
      case size :
         result = some_arguments->count() == 1 ;
         if( result )
         {
            MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = ( k0==DoubleVector || k0==IntVector ||
                       k0==BoolVector || k0==StringVector ) ;
         }
         break ;
      case component :
         result = some_arguments->count() == 2 ;
         if( result )
         {
            MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
            MAC_Data::Type k1 =  extract_arg( some_arguments, 1 )->data_type() ;
            result = ( k0==DoubleVector || k0==IntVector ||
                       k0==BoolVector || k0==StringVector ) && ( k1==Int ) ;
         }
         break ;
      case increasing :
         result = some_arguments->count() == 1 ;
         if( result )
         {
            MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = ( k0==DoubleVector || k0==IntVector ) ;
         }
         break ;
      case reverse :
         result = some_arguments->count() == 1 ;
         if( result )
         {
            MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = ( k0==DoubleVector || k0==IntVector ||
                       k0==BoolVector || k0==StringVector ) ;
         }
         break ;
      case sum :
         result = some_arguments->count() == 1 ;
         if( result )
         {
            MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = ( k0==DoubleVector || k0==IntVector ) ;
         }
         break ;
      case greater :
         result = some_arguments->count() == 2 ;
         if( result )
         {
            MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
            MAC_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
            result = ( k0==DoubleVector && k1==Double ) ||
                     ( k0==IntVector    && k1==Int ) ;
         }
         break ;
      case nvector :
         result = some_arguments->count() == 2 ;
         if( result )
         {
            MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
            MAC_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
            result = k0==Int &&
               ( k1==Int || k1==Double || k1==String || k1==Bool  ) ;
         }
         break ;
      case apply :
         result = ( some_arguments->count() == 3 ||
                    some_arguments->count() == 4 ) ;
         if( result )
         {
            MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
            MAC_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
            MAC_Data::Type k2 = extract_arg( some_arguments, 2 )->data_type() ;
            result = ( k0 == DoubleVector && k1==Double ) ||
                     ( k0 == IntVector && k1==Int ) ||
                     ( k0 == BoolVector && k1==Bool ) ||
                     ( k0 == StringVector && k1==String ) ;
            result &= ( k2 == String ) ;
            if( result )
            {
               MAC_Data::Type k =
                  MAC_Variable::data_type(
                     extract_arg( some_arguments, 2 )->to_string() ) ;
               result = ( k==k1 ) ;
            }
            if( result && some_arguments->count() == 4 )
            {
               MAC_Data::Type k3 = extract_arg( some_arguments, 3 )->data_type() ;
               result = ( k3 == String ) ;
               if( result )
               {
                  MAC_Data::Type k =
                     MAC_Variable::data_type(
                        extract_arg( some_arguments, 3 )->to_string() ) ;
                  result = ( k==Int ) ;
               }
            }     
         }
         break ;        
      default :
         MAC_Error::object()->raise_internal( "not implemented" ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_VectorExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: data_type" ) ;

   Type result = Undefined ;
   
   switch( OP )
   {
      case cond_vect :
      case vect :
      case nvector :
         {
            size_t idx = ( OP==vect ? 0 : 1 ) ;
            
            Type tmp = arg(idx)->data_type() ;
            if( tmp==Double )
            {
               result = DoubleVector ;
            }
            else if( tmp==Int )
            {
               result = IntVector ;
            }
            else if( tmp==Bool )
            {
               result = BoolVector ;
            }
            else if( tmp==String )
            {
               result = StringVector ;
            }
         }
         break ;
      case size :
         result = Int ;
         break ;
      case component :
         {
            Type tmp = arg(0)->data_type() ;
            if( tmp==DoubleVector )
            {
               result = Double ;
            }
            else if( tmp==IntVector )
            {
               result = Int ;
            }
            else if( tmp==BoolVector )
            {
               result = Bool ;
            }
            else if( tmp==StringVector )
            {
               result = String ;
            }
         }
         break ;
      case increasing :
      case greater :
         result = Bool ;
         break ;
      case sum :
         {
            Type tmp = arg(0)->data_type() ;
            if( tmp==DoubleVector )
            {
               result = Double ;
            }
            else if( tmp==IntVector )
            {
               result = Int ;
            }
         }
         break ;
      case reverse :
         result = arg(0)->data_type() ;
         break ;
      case apply :
         result = arg(0)->data_type() ;
         break ;
      default :
         MAC_Error::object()->raise_internal( "not implemented" ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_VectorExp:: declare( MAC_List* lst ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: declare" ) ;
   MAC_CHECK_PRE( declare_PRE( lst ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Expression::declare( lst ) ;
   if( OP == apply )
   {
      MAC_Variable const* v = MAC_Variable::object( arg(2)->to_string() ) ;
      if( lst->has( v ) ) lst->remove( v ) ;
      if( nb_arguments()==4 )
      {
         MAC_Variable const* vic =
                              MAC_Variable::object( arg(3)->to_string() ) ;
         if( lst->has( vic ) ) lst->remove( vic ) ;
      }
   }
}

//----------------------------------------------------------------------
bool
MAC_VectorExp:: context_has_required_variables( 
                                           MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: context_has_required_variables" ) ;
   MAC_CHECK_PRE( context_has_required_variables_PRE( ct ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   bool result = false ;
   if( OP == apply )
   {
      MAC_Variable const* v = MAC_Variable::object( arg(2)->to_string() ) ;
      MAC_Variable const* vic =
         ( nb_arguments()==4 ?
                        MAC_Variable::object( arg(3)->to_string() ) : 0 ) ;
      MAC_Context const* ctx = create_apply_context( 0, ct, v, vic ) ;
      result = MAC_Expression::context_has_required_variables( ctx ) ;
      ctx->destroy() ;
   }
   else
   {
      result = MAC_Expression::context_has_required_variables( ct ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_VectorExp:: value_can_be_evaluated( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: value_can_be_evaluated" ) ;

   bool result = false ;
   if( OP == apply )
   {
      MAC_Variable const* v = MAC_Variable::object( arg(2)->to_string() ) ;
      MAC_Variable const* vic =
         ( nb_arguments()==4 ?
                        MAC_Variable::object( arg(3)->to_string() ) : 0 ) ;
      MAC_Context const* ctx = create_apply_context( 0, ct, v, vic ) ;
      result = MAC_Expression::value_can_be_evaluated( ctx ) ;
      ctx->destroy() ;
   }
   else
   {
      result = MAC_Expression::value_can_be_evaluated( ct ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
stringVector const&
MAC_VectorExp:: undefined_variables( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: undefined_variables" ) ;

   static stringVector result(0) ;
   if( OP == apply )
   {
      MAC_Variable const* v = MAC_Variable::object( arg(2)->to_string() ) ;
      MAC_Variable const* vic =
         ( nb_arguments()==4 ?
                        MAC_Variable::object( arg(3)->to_string() ) : 0 ) ;
      MAC_Context const* ctx = create_apply_context( 0, ct, v, vic ) ;
      result = MAC_Expression::undefined_variables( ctx ) ;
      ctx->destroy() ;
   }
   else
   {
      result = MAC_Expression::undefined_variables( ct ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_VectorExp:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: to_bool" ) ;
   MAC_CHECK_PRE( to_bool_PRE(ct) ) ;
   bool result = false ;
   if( OP==component )
   {
      result = arg(0)->to_bool_vector(ct)( arg(1)->to_int(ct) ) ;
   }
   else if( OP == increasing )
   {
      MAC_Data const* VECTOR = arg(0) ;
      result = true ;
      MAC_Data::Type dt = VECTOR->data_type() ;
      
      if( dt==IntVector )
      {
         intVector const& vec =  VECTOR->to_int_vector(ct) ;
         for( size_t i=0 ; result && i<vec.size()-1 ; i++ )
            result = result && vec(i) <= vec(i+1) ;
      }
      else if( dt==DoubleVector )
      {
         doubleVector const& vec =  VECTOR->to_double_vector(ct) ;
         for( size_t i=0 ; result && i<vec.size()-1 ; i++ )
            result = result && vec(i) <= vec(i+1) ;
      }
      else
         MAC_Error::object()->raise_internal( "not implemented" ) ;
   }
   else if( OP == greater )
   {
      MAC_Data const* VECTOR = arg(0) ;
      MAC_Data const* VAL = arg(1) ;
      result = true ;
      MAC_Data::Type dt = VECTOR->data_type() ;
      
      if( dt==IntVector )
      {
         intVector const& vec =  VECTOR->to_int_vector(ct) ;
         int val = VAL->to_int(ct) ;
         for( size_t i=0 ; result && i<vec.size() ; i++ )
            result = result && vec(i) >= val ;
      }
      else if( dt==DoubleVector )
      {
         doubleVector const& vec =  VECTOR->to_double_vector(ct) ;
         double val  = VAL->to_double(ct) ;
         for( size_t i=0 ; result && i<vec.size() ; i++ )
            result = result && vec(i) >= val ;
      }
      else 
         MAC_Error::object()->raise_internal( "not implemented" ) ;
   }
   else
      MAC_Error::object()->raise_internal( "not implemented" ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
double
MAC_VectorExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE(ct) ) ;
   double result = MAC::bad_double() ;
   if( OP == component )
   {
      result = arg(0)->to_double_vector(ct)( arg(1)->to_int(ct) ) ;
   }
   else
   {
      MAC_CHECK( OP==sum ) ;
      doubleVector const& v = arg(0)->to_double_vector(ct) ;
      result = 0.0 ;
      for( size_t i=0 ; i<v.size() ; i++ )
         result += v(i) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
int
MAC_VectorExp:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE(ct) ) ;
   
   int result = -1 ;
   
   switch( OP )
   {
      case size :
         {
            MAC_Data const* VECTOR = arg(0) ;
            switch( VECTOR->data_type() )
            {
               case IntVector :
                  result =  VECTOR->to_int_vector(ct).size() ;
                  break ;
               case DoubleVector :
                  result =  VECTOR->to_double_vector(ct).size() ;
                  break ;
               case BoolVector :
                  result =  VECTOR->to_bool_vector(ct).size() ;
                  break ;
               case StringVector :
                  result =  VECTOR->to_string_vector(ct).size() ;
                  break ;
               default :
                  MAC_Error::object()->raise_internal( "not implemented" ) ;
            }
         }
         break ;
      case component :
         result = arg(0)->to_int_vector(ct)( arg(1)->to_int(ct) ) ;
         break ;
      case sum :
         {
            intVector const& iv = arg(0)->to_int_vector(ct) ;
            result = iv.sum() ;
         }
         break ;   
      default :
         MAC_Error::object()->raise_internal( "not implemented" ) ;
         break ;
   }   
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_VectorExp:: to_string( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: to_string" ) ;
   MAC_CHECK_PRE( to_string_PRE(ct) ) ;
   std::string const& result =
                    arg(0)->to_string_vector(ct)( arg(1)->to_int(ct) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
doubleVector const&
MAC_VectorExp:: to_double_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: to_double_vector" ) ;
   MAC_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   if( OP==vect )
   {
      RESULT_D.re_initialize( nb_arguments() ) ;
      size_t idx = 0 ;
      for( MY_IT->start(); MY_IT->is_valid(); MY_IT->go_next() )
      {
         RESULT_D(idx++) =
            static_cast<MAC_Data const*>( MY_IT->item() )->to_double(ct) ;
      }
   }
   else if( OP==nvector )
   {
      RESULT_D.re_initialize( arg(0)->to_int(ct) ) ;
      RESULT_D.set( arg(1)->to_double(ct) ) ;
   }
   else if( OP==cond_vect )
   {
      RESULT_D.re_initialize( 0 ) ;
      for( MY_IT->start(); MY_IT->is_valid(); MY_IT->go_next() )
      {
         bool test = static_cast<MAC_Data const*>( MY_IT->item() )->to_bool(ct) ;
         MY_IT->go_next() ;
         if( test )
         {   
            RESULT_D.append(
               static_cast<MAC_Data const*>( MY_IT->item() )->to_double(ct) ) ;
         }
               
      }
   }
   else if( OP==apply )
   {
      RESULT_D = arg(0)->to_double_vector( ct ) ;
      MAC_Variable const* v = MAC_Variable::object( arg(2)->to_string() ) ;
      MAC_Variable const* vic =
         ( nb_arguments()==4 ?
                        MAC_Variable::object( arg(3)->to_string() ) : 0 ) ;
      MAC_Context const* ctx = create_apply_context( 0, ct, v, vic ) ;
      MAC_Double* d = static_cast<MAC_Double*>( ctx->value( v ) ) ;
      MAC_Int* ic = 
         ( vic!=0 ? static_cast<MAC_Int*>( ctx->value( vic ) ) : 0 ) ;
      for( size_t i=0 ; i<RESULT_D.size() ; ++i )
      {
         d->set( RESULT_D(i) ) ;
         if( ic != 0 ) ic->set( i ) ;
         RESULT_D(i) = arg(1)->to_double( ctx ) ;
      }
      ctx->destroy() ;
   }
   else if( OP==reverse )
   {
      doubleVector const& v = arg(0)->to_double_vector( ct ) ;
      size_t n = v.size() ;
      
      RESULT_D.resize(n) ;
      for( size_t i=0 ; i<n ; i++ )
         RESULT_D(i) = v(n-i-1) ;
   }
   
   return( RESULT_D ) ;
}

//----------------------------------------------------------------------
intVector const&
MAC_VectorExp:: to_int_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: to_int_vector" ) ;
   MAC_CHECK_PRE( to_int_vector_PRE( ct ) ) ;
   
   if( OP==vect )
   {
      RESULT_I.re_initialize( nb_arguments() ) ;
      size_t idx = 0 ;
      for( MY_IT->start(); MY_IT->is_valid(); MY_IT->go_next() )
      {
         RESULT_I(idx++) =
            static_cast<MAC_Data const*>( MY_IT->item() )->to_int(ct) ;
      }
   }
   else if( OP==nvector )
   {
      RESULT_I.re_initialize( arg(0)->to_int(ct) ) ;
      RESULT_I.set( arg(1)->to_int(ct) ) ;
   }
   else if( OP==cond_vect )
   {
      RESULT_I.re_initialize( 0 ) ;
      for( MY_IT->start(); MY_IT->is_valid(); MY_IT->go_next() )
      {
         bool test = static_cast<MAC_Data const*>( MY_IT->item() )->to_bool(ct) ;
         MY_IT->go_next() ;
         if( test )
         {   
            RESULT_I.append(
               static_cast<MAC_Data const*>( MY_IT->item() )->to_int(ct) ) ;
         }
               
      }
      
   }
   else if( OP==apply )
   {
      RESULT_I = arg(0)->to_int_vector( ct ) ;
      MAC_Variable const* v = MAC_Variable::object( arg(2)->to_string() ) ;
      MAC_Variable const* vic =
         ( nb_arguments()==4 ?
                        MAC_Variable::object( arg(3)->to_string() ) : 0 ) ;
      MAC_Context const* ctx = create_apply_context( 0, ct, v, vic ) ;
      MAC_Int* d = static_cast<MAC_Int*>( ctx->value( v ) ) ;
      MAC_Int* ic = 
         ( vic!=0 ? static_cast<MAC_Int*>( ctx->value( vic ) ) : 0 ) ;
      for( size_t i=0 ; i<RESULT_I.size() ; ++i )
      {
         d->set( RESULT_I(i) ) ;
         if( ic != 0 ) ic->set( i ) ;
         RESULT_I(i) = arg(1)->to_int( ctx ) ;
      }
      ctx->destroy() ;
   }
   else if( OP==reverse )
   {
      intVector const& v = arg(0)->to_int_vector( ct ) ;
      size_t n = v.size() ;
      
      RESULT_I.resize(n) ;
      for( size_t i=0 ; i<n ; i++ )
         RESULT_I(i) = v(n-i-1) ;
   }
   
   return( RESULT_I ) ;
}

//----------------------------------------------------------------------
boolVector const&
MAC_VectorExp:: to_bool_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: to_bool_vector" ) ;
   MAC_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;

   if( OP==vect )
   {
      RESULT_B.re_initialize( nb_arguments() ) ;
      size_t idx = 0 ;
      for( MY_IT->start(); MY_IT->is_valid(); MY_IT->go_next() )
      {
         RESULT_B(idx++) =
            static_cast<MAC_Data const*>( MY_IT->item() )->to_bool(ct) ;
      }
   }
   else if( OP==nvector )
   {
      RESULT_B.re_initialize( arg(0)->to_int(ct) ) ;
      RESULT_B.set( arg(1)->to_bool(ct) ) ;
   }
   else if( OP==cond_vect )
   {
      RESULT_B.re_initialize( 0 ) ;
      for( MY_IT->start(); MY_IT->is_valid(); MY_IT->go_next() )
      {
         bool test = static_cast<MAC_Data const*>( MY_IT->item() )->to_bool(ct) ;
         MY_IT->go_next() ;
         if( test )
         {   
            RESULT_B.append(
               static_cast<MAC_Data const*>( MY_IT->item() )->to_bool(ct) ) ;
         }
               
      }
   }
   else if( OP==apply )
   {
      RESULT_B = arg(0)->to_bool_vector( ct ) ;
      MAC_Variable const* v = MAC_Variable::object( arg(2)->to_string() ) ;
      MAC_Variable const* vic =
         ( nb_arguments()==4 ?
                        MAC_Variable::object( arg(3)->to_string() ) : 0 ) ;
      MAC_Context const* ctx = create_apply_context( 0, ct, v, vic ) ;
      MAC_Bool* d = static_cast<MAC_Bool*>( ctx->value( v ) ) ;
      MAC_Int* ic = 
         ( vic!=0 ? static_cast<MAC_Int*>( ctx->value( vic ) ) : 0 ) ;
      for( size_t i=0 ; i<RESULT_B.size() ; ++i )
      {
         d->set( RESULT_B(i) ) ;
         if( ic != 0 ) ic->set( i ) ;
         RESULT_B(i) = arg(1)->to_bool( ctx ) ;
      }
      ctx->destroy() ;
   }
   else if( OP==reverse )
   {
      boolVector const& v = arg(0)->to_bool_vector( ct ) ;
      size_t n = v.size() ;
      
      RESULT_B.resize(n) ;
      for( size_t i=0 ; i<n ; i++ )
         RESULT_B(i) = v(n-i-1) ;
   }
   
   return( RESULT_B ) ;
}

//----------------------------------------------------------------------
stringVector const&
MAC_VectorExp:: to_string_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: to_string_vector" ) ;
   MAC_CHECK_PRE( to_string_vector_PRE( ct ) ) ;

   if( OP==vect )
   {
      RESULT_S.re_initialize( nb_arguments() ) ;
      size_t idx = 0 ;
      for( MY_IT->start(); MY_IT->is_valid(); MY_IT->go_next() )
      {
         RESULT_S(idx++) =
            static_cast<MAC_Data const*>( MY_IT->item() )->to_string(ct) ;
      }
   }
   else if( OP==nvector )
   {
      RESULT_S.re_initialize( arg(0)->to_int(ct) ) ;
      RESULT_S.set( arg(1)->to_string(ct) ) ;
   }
   else if( OP==cond_vect )
   {
      RESULT_S.re_initialize( 0 ) ;
      for( MY_IT->start(); MY_IT->is_valid(); MY_IT->go_next() )
      {
         bool test = static_cast<MAC_Data const*>( MY_IT->item() )->to_bool(ct) ;
         MY_IT->go_next() ;
         if( test )
         {   
            RESULT_S.append(
               static_cast<MAC_Data const*>( MY_IT->item() )->to_string(ct) ) ;
         }
               
      }
      
   }
   else if( OP==apply )
   {
      RESULT_S = arg(0)->to_string_vector( ct ) ;
      MAC_Variable const* v = MAC_Variable::object( arg(2)->to_string() ) ;
      MAC_Variable const* vic =
         ( nb_arguments()==4 ?
                        MAC_Variable::object( arg(3)->to_string() ) : 0 ) ;
      MAC_Context const* ctx = create_apply_context( 0, ct, v, vic ) ;
      MAC_String* d = static_cast<MAC_String*>( ctx->value( v ) ) ;
      MAC_Int* ic = 
         ( vic!=0 ? static_cast<MAC_Int*>( ctx->value( vic ) ) : 0 ) ;
      for( size_t i=0 ; i<RESULT_S.size() ; ++i )
      {
         d->set( RESULT_S(i) ) ;
         if( ic != 0 ) ic->set( i ) ;
         RESULT_S(i) = arg(1)->to_string( ctx ) ;
      }
      ctx->destroy() ;
   }
   else if( OP==reverse )
   {
      stringVector const& v = arg(0)->to_string_vector( ct ) ;
      size_t n = v.size() ;
      
      RESULT_S.resize(n) ;
      for( size_t i=0 ; i<n ; i++ )
         RESULT_S(i) = v(n-i-1) ;
   }
   
   return( RESULT_S ) ;
}

//----------------------------------------------------------------------
MAC_Data*
MAC_VectorExp:: create_derivative( MAC_Object* a_owner,
                                   MAC_Variable const* var,
                                   MAC_Context const* ct ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;

   MAC_Data* result = 0 ;
   switch( OP )
   {
      case vect :
         {
            MAC_List* lst = MAC_List::create( 0 ) ;
            for( MY_IT->start(); MY_IT->is_valid(); MY_IT->go_next() )
            {
               MAC_Data const* val = static_cast<MAC_Data const*>( MY_IT->item() ) ;
               lst->append( val->create_derivative(lst,var,ct) ) ;
            }

            result = MAC_Expression::create( a_owner, "vector", lst ) ;
            lst->set_owner( result ) ;
         }
         break ;
      case size :
         result = MAC_Double::create( a_owner, 0.0 ) ;
         break ;
      case component :
         {
            MAC_List* lst = MAC_List::create( 0 ) ;
            lst->append( arg(0)->create_derivative( lst, var, ct ) ) ;
            lst->append( arg(1)->create_clone( lst ) ) ;
            result = MAC_Expression::create( a_owner, "component", lst ) ;
            lst->set_owner( result ) ;
         }
         break ;
      default :
         MAC_Error::object()->raise_internal( "not implemented" ) ;
   }     
   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Context const*
MAC_VectorExp:: create_apply_context( MAC_Object* a_owner,
                                      MAC_Context const* ctx,
                                      MAC_Variable const* v,
                                      MAC_Variable const* vic ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorExp:: create_apply_context" ) ;
   MAC_CHECK( OP == apply ) ;
   MAC_CHECK( v != 0 ) ;
   MAC_CHECK( v->data_type() == Double ||
              v->data_type() == Int ||
              v->data_type() == Bool ||
              v->data_type() == String ) ;
   MAC_CHECK( IMPLIES( vic != 0, vic->data_type() == Int ) ) ;

   MAC_ContextSimple* a_ctx = MAC_ContextSimple::create( 0 ) ;
   if( v->data_type() == Double )
   {
      MAC_Double* d = MAC_Double::create( a_ctx, 0. ) ;
      a_ctx->extend( v, d ) ;
   }
   else if( v->data_type() == Int )
   {
      MAC_Int* d = MAC_Int::create( a_ctx, 0 ) ;
      a_ctx->extend( v, d ) ;
   }
   else if( v->data_type() == Bool )
   {
      MAC_Bool* d = MAC_Bool::create( a_ctx, false ) ;
      a_ctx->extend( v, d ) ;
   }
   else if( v->data_type() == String )
   {
      MAC_String* d = MAC_String::create( a_ctx, "" ) ;
      a_ctx->extend( v, d ) ;
   }
   if( vic != 0 )
   {
      MAC_Int* d = MAC_Int::create( a_ctx, 0 ) ;
      a_ctx->extend( vic, d ) ;
   }
   
   MAC_Context* result = 0 ;
   if( ctx != 0 )
   {
      result = MAC_ContextPair::create( a_owner, ctx, a_ctx ) ;
      a_ctx->set_owner( result ) ;
   }
   else
   {
      a_ctx->set_owner( a_owner ) ;
      result = a_ctx ;
   }
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->has_variable( v ) ) ;
   MAC_CHECK_POST( IMPLIES( vic!=0, result->has_variable( vic ) ) ) ;
   return( result ) ;
}
