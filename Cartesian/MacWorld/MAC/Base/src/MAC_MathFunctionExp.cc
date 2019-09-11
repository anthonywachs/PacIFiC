#include <MAC_MathFunctionExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Double.hh>
#include <MAC_Error.hh>
#include <MAC_Int.hh>
#include <MAC_List.hh>
#include <MAC_Sequence.hh>
#include <MAC.hh>

#include <iostream>

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Sqr = new MAC_MathFunctionExp( "sqr", Sqr ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Sqrt = new MAC_MathFunctionExp( "sqrt", Sqrt ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Pow = new MAC_MathFunctionExp( "pow", Pow ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Exp = new MAC_MathFunctionExp( "exp", Exp ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Log = new MAC_MathFunctionExp( "log", Log ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Log10 = new MAC_MathFunctionExp( "log10", Log10 ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Sin = new MAC_MathFunctionExp( "sin", Sin ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Cos = new MAC_MathFunctionExp( "cos", Cos ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Tan = new MAC_MathFunctionExp( "tan", Tan ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_ASin = new MAC_MathFunctionExp( "asin", ASin ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_ACos = new MAC_MathFunctionExp( "acos", ACos ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_ATan = new MAC_MathFunctionExp( "atan", ATan ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_ATan2 = new MAC_MathFunctionExp( "atan2", ATan2 ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Sinh = new MAC_MathFunctionExp( "sinh", Sinh ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Cosh = new MAC_MathFunctionExp( "cosh", Cosh ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Tanh = new MAC_MathFunctionExp( "tanh", Tanh ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_ASinh = new MAC_MathFunctionExp( "asinh", ASinh ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_ACosh = new MAC_MathFunctionExp( "acosh", ACosh ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_ATanh = new MAC_MathFunctionExp( "atanh", ATanh ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_J0 = new MAC_MathFunctionExp( "j0", J0 ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_J1 = new MAC_MathFunctionExp( "j1", J1 ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Jn = new MAC_MathFunctionExp( "jn", Jn ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Y0 = new MAC_MathFunctionExp( "y0", Y0 ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Y1 = new MAC_MathFunctionExp( "y1", Y1 ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Yn = new MAC_MathFunctionExp( "yn", Yn ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Gamma = new MAC_MathFunctionExp( "gamma", Gamma ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_LGamma = new MAC_MathFunctionExp( "lgamma", LGamma ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Erf = new MAC_MathFunctionExp( "erf", Erf ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Erfc = new MAC_MathFunctionExp( "erfc", Erfc ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_IGamma =
                        new MAC_MathFunctionExp( "incomplete_gamma", IGamma ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_En = new MAC_MathFunctionExp( "En", En ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Ei = new MAC_MathFunctionExp( "Ei", Ei ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_ABS = new MAC_MathFunctionExp( "abs", Abs ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Floor = new MAC_MathFunctionExp( "floor", Floor ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_Ceil = new MAC_MathFunctionExp( "ceil", Ceil ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_DblEq = new MAC_MathFunctionExp( "double_equality", DblEq ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_RandDbl = new MAC_MathFunctionExp( "random_double", RandDbl ) ;

MAC_MathFunctionExp const* 
MAC_MathFunctionExp::PROTOTYPE_RandInt = new MAC_MathFunctionExp( "rand", RandInt ) ;

struct MAC_MathFunctionExp_ERROR
{
   static void n0( std::string const& f_name, std::string const& op_name ) ;
} ;

//----------------------------------------------------------------------
MAC_MathFunctionExp:: MAC_MathFunctionExp( std::string const& a_name,
                                           MathOp a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( a_op )
   , ARG1( 0 )
   , ARG2( 0 )
   , ARG3( 0 )
   , ARG4( 0 )
   , ARG5( 0 )
{
   MAC_LABEL( "MAC_MathFunctionExp:: MAC_MathFunctionExp" ) ;
}

//----------------------------------------------------------------------
MAC_MathFunctionExp*
MAC_MathFunctionExp:: create_replica(
                             MAC_Object* a_owner,
                             MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MathFunctionExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_MathFunctionExp* result = new MAC_MathFunctionExp( a_owner, 
                                                          name(), 
                                                          argument_list, 
                                                          OP ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_MathFunctionExp:: MAC_MathFunctionExp( MAC_Object* a_owner,
                                           std::string const& a_name,
                                           MAC_Sequence const* argument_list,
                                           MathOp a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , ARG1( nb_arguments()>0 ? arg(0) : 0 )
   , ARG2( nb_arguments()>1 ? arg(1) : 0 )
   , ARG3( nb_arguments()>2 ? arg(2) : 0 )
   , ARG4( nb_arguments()>3 ? arg(3) : 0 )
   , ARG5( nb_arguments()>4 ? arg(4) : 0 )
{
   MAC_LABEL( "MAC_MathFunctionExp:: MAC_MathFunctionExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_MathFunctionExp:: ~MAC_MathFunctionExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MathFunctionExp:: ~MAC_MathFunctionExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_MathFunctionExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MathFunctionExp:: data_type" ) ;

   MAC_Data::Type result = MAC_Data::Undefined ;

   switch( OP )
   {
      case Sqr :
      case Sqrt :
      case Pow :
      case Exp :
      case Log :
      case Log10 :
      case Sin :
      case Cos :
      case Tan :
      case ASin :
      case ACos :
      case ATan :
      case ATan2 :
      case Sinh :
      case Cosh :
      case Tanh :
      case ASinh :
      case ACosh :
      case ATanh :
      case J0 :
      case J1 :
      case Jn :
      case Y0 :
      case Y1 :
      case Yn :
      case Gamma :
      case LGamma :
      case IGamma :
      case Erf :
      case Erfc :
      case En :
      case Ei :
      case Floor :
      case Ceil :
      case RandDbl :
         result = MAC_Data::Double ;
         break ;
      case Abs :
         result = ARG1->data_type() ;
         break ;
      case DblEq :
         result = MAC_Data::Bool ;
         break ;
      case RandInt :
         result = MAC_Data::Int ;
         break ;
      default :
         MAC_MathFunctionExp_ERROR::n0( "data_type", name() ) ;
         break ;
   } 
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_MathFunctionExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "undefined" ;

   switch( OP )
   {
      case Sqr :
      case Sqrt :
      case Exp :
      case Log :
      case Log10 :
      case Sin :
      case Cos :
      case Tan :
      case ASin :
      case ACos :
      case ATan :
      case Sinh :
      case Cosh :
      case Tanh :
      case ASinh :
      case ACosh :
      case ATanh :
      case J0 :
      case J1 :
      case Y0 :
      case Y1 :
      case Gamma :
      case LGamma :
      case Erf :
      case Erfc :
      case Floor :
      case Ceil :
         result = name() + "(DS)" ;
         break ;
      case Jn :
      case Yn :
         result = name() + "(IS,DS)" ;
         break ;
      case Pow :
      case ATan2 :
         result =  name() + "(DS,DS)" ;
         break ;
      case IGamma :
         result =  name() + "(DS,DS[,DS,DS,IS])" ;
         break ;
      case En :
         result =  name() + "(IS,DS[,DS,DS,IS])" ;
         break ;
      case Ei :
         result =  name() + "(DS[,DS,DS,IS])" ;
         break ;
      case Abs :
         result =  name() + "(DS|IS)" ;
         break ;
      case DblEq :
         result = name() + "(DS,DS,DS,DS)" ;
         break ;
      case RandDbl :
      case RandInt :
         result = name() + "()" ;
         break ;
      default :
         MAC_MathFunctionExp_ERROR::n0( "usage", name() ) ;
         break ;
   }    
   
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_MathFunctionExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MathFunctionExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = false ;
   
   switch( OP )
   {
      case Sqr :
      case Sqrt :
      case Exp :
      case Log :
      case Log10 :
      case Sin :
      case Cos :
      case Tan :
      case ASin :
      case ACos :
      case ATan :
      case Sinh :
      case Cosh :
      case Tanh :
      case ASinh :
      case ACosh :
      case ATanh :
      case J0 :
      case J1 :
      case Y0 :
      case Y1 :
      case Gamma :
      case LGamma :
      case Erf :
      case Erfc :
      case Floor :
      case Ceil :
         result = ( some_arguments->count()==1 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = t0==Double ;
         }
         break ;
      case Abs :
         result = ( some_arguments->count()==1 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = t0==Double || t0==Int ;
         }
         break ;
      case Jn :
      case Yn :
         result = ( some_arguments->count()==2 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            result = ( t0==Int && t1==Double ) ;
         }
         break ;
      case Pow :
      case ATan2 :  
         result = ( some_arguments->count()==2 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            result = ( t0==Double && t1==Double ) ;
         }
         break ;
      case IGamma :
         result = ( some_arguments->count()==2 ||
                    some_arguments->count()==5 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            result = ( t0==Double && t1==Double ) ;
            if( result && some_arguments->count()==5 )
            {
               Type t2 = extract_arg( some_arguments, 2 )->data_type() ;
               Type t3 = extract_arg( some_arguments, 3 )->data_type() ;
               Type k4 = extract_arg( some_arguments, 4 )->data_type() ;
               result = ( t2==Double && t3==Double && k4==Int ) ;
            }
         }
         break ;
      case En :
         result = ( some_arguments->count()==2 ||
                    some_arguments->count()==5 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            result = ( t0==Int && t1==Double ) ;
            if( result && some_arguments->count()==5 )
            {
               Type t2 = extract_arg( some_arguments, 2 )->data_type() ;
               Type t3 = extract_arg( some_arguments, 3 )->data_type() ;
               Type k4 = extract_arg( some_arguments, 4 )->data_type() ;
               result = ( t2==Double && t3==Double && k4==Int ) ;
            }
         }
         break ;
      case Ei :
         result = ( some_arguments->count()==1 ||
                    some_arguments->count()==4 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = ( t0==Double ) ;
            if( result && some_arguments->count()==4 )
            {
               Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
               Type t2 = extract_arg( some_arguments, 2 )->data_type() ;
               Type t3 = extract_arg( some_arguments, 3 )->data_type() ;
               result = ( t1==Double && t2==Double && t3==Int ) ;
            }
         }
         break ;
      case DblEq :
         result = ( some_arguments->count()==4 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            Type t2 = extract_arg( some_arguments, 2 )->data_type() ;
            Type t3 = extract_arg( some_arguments, 3 )->data_type() ;
            result = ( t0==Double && t1==Double && t2==Double && t3==Double ) ;
         }
         break ;
      case RandDbl :
      case RandInt :
         result = ( some_arguments->count()==0 ) ;
         break ;
      default :
         MAC_MathFunctionExp_ERROR::n0( "matches_args", name() ) ;
         break ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_MathFunctionExp:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MathFunctionExp:: to_bool" ) ;
   MAC_CHECK_PRE( to_bool_PRE(ct) ) ;
   
   bool result = false ;

   switch( OP )
   {
      case DblEq :
         result = MAC::double_equality( ARG1->to_double( ct ),
                                        ARG2->to_double( ct ),
                                        ARG3->to_double( ct ),
                                        ARG4->to_double( ct ) ) ;
         break ;
      default :
         MAC_MathFunctionExp_ERROR::n0( "to_bool", name() ) ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
double
MAC_MathFunctionExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MathFunctionExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE(ct) ) ;
   
   double result = MAC::max_double() ;

   switch( OP )
   {
      case Sqr :
         result = MAC::sqr( ARG1->to_double(ct) ) ;
         break ;
      case Sqrt :
         result = MAC::sqrt( ARG1->to_double(ct) ) ;
         break ;
      case Pow :
         result = MAC::pow( ARG1->to_double(ct), ARG2->to_double(ct) ) ;
         break ;
      case Exp :
         result = MAC::exp( ARG1->to_double(ct) ) ;
         break ;
      case Log :
         result = MAC::log( ARG1->to_double(ct) ) ;
         break ;
      case Log10 :
         result = MAC::log10( ARG1->to_double(ct) ) ;
         break ;
      case Sin :
         result = MAC::sin( ARG1->to_double(ct) ) ;
         break ;
      case Cos :
         result = MAC::cos( ARG1->to_double(ct) ) ;
         break ;
      case Tan :
         result = MAC::tan( ARG1->to_double(ct) ) ;
         break ;
      case ASin :
         result = MAC::asin( ARG1->to_double(ct) ) ;
         break ;
      case ACos :
         result = MAC::acos( ARG1->to_double(ct) ) ;
         break ;
      case ATan :
         result = MAC::atan( ARG1->to_double(ct) ) ;
         break ;
      case ATan2 :
         result = MAC::atan2( ARG1->to_double(ct), ARG2->to_double(ct) ) ;
         break ;
      case Sinh :
         result = MAC::sinh( ARG1->to_double(ct) ) ;
         break ;
      case Cosh :
         result = MAC::cosh( ARG1->to_double(ct) ) ;
         break ;
      case Tanh :
         result = MAC::tanh( ARG1->to_double(ct) ) ;
         break ;
      case ASinh :
         result = MAC::asinh( ARG1->to_double(ct) ) ;
         break ;
      case ACosh :
         result = MAC::acosh( ARG1->to_double(ct) ) ;
         break ;
      case ATanh :
         result = MAC::atanh( ARG1->to_double(ct) ) ;
         break ;
      case J0 :
         result = MAC::j0( ARG1->to_double(ct) ) ;
         break ;
      case J1 :
         result = MAC::j1( ARG1->to_double(ct) ) ;
         break ;
      case Jn :
         result = MAC::jn( ARG1->to_int(ct), ARG2->to_double(ct) ) ;
         break ;
      case Y0 :
         result = MAC::y0( ARG1->to_double(ct) ) ;
         break ;
      case Y1 :
         result = MAC::y1( ARG1->to_double(ct) ) ;
         break ;
      case Yn :
         result = MAC::yn( ARG1->to_int(ct), ARG2->to_double(ct) ) ;
         break ;
      case Gamma :
         result = MAC::gamma( ARG1->to_double(ct) ) ;
         break ;
      case LGamma :
         result = MAC::lgamma( ARG1->to_double(ct) ) ;
         break ;
      case Erf :
         result = MAC::erf( ARG1->to_double(ct) ) ;
         break ;
      case Erfc :
         result = MAC::erfc( ARG1->to_double(ct) ) ;
         break ;
      case IGamma :
         {
            double const a = ARG1->to_double(ct) ;
            if( a <= 0 )
            {
               raise_error( "first argument should be stricktly positive" ) ;
            }
            double const x = ARG2->to_double(ct) ;
            if( x <= 0 )
            {
               raise_error( "second argument should be stricktly positive" ) ;
            }
            if( ARG3 == 0 )
            {
               result = MAC::incomplete_gamma( a, x ) ;
            }
            else
            {
               double const eps = ARG3->to_double(ct) ;
               if( eps <= 0 )
               {
                  raise_error( "third argument should be stricktly positive" ) ;
               }
               double const err = ARG4->to_double(ct) ;
               if( err <= 0 )
               {
                  raise_error( "fourth argument should be stricktly positive" ) ;
               }
               int const iter_max = ARG5->to_int(ct) ;
               if( iter_max <= 0 )
               {
                  raise_error( "fifth argument should be stricktly positive" ) ;
               }
               result = MAC::incomplete_gamma( a, x,
                                               eps, err, (size_t) iter_max ) ;
            }
         }
         break ;
      case En :
         {
            int const n = ARG1->to_int(ct) ;
            if( n < 0 )
            {
               raise_error( "first argument should be positive" ) ;
            }
            double const x = ARG2->to_double(ct) ;
            if( x < 0 )
            {
               raise_error( "second argument should be positive" ) ;
            }
            if( ARG3 == 0 )
            {
               if( x<1.E-30 && n<2 )
               {
                  raise_error( "first argument less than 2 and null second argument" ) ;
               }
               result = MAC::En( (size_t) n, x ) ;
            }
            else
            {
               double const eps = ARG3->to_double(ct) ;
               if( eps <= 0 )
               {
                  raise_error( "third argument should be stricktly positive" ) ;
               }
               double const err = ARG4->to_double(ct) ;
               if( err <= 0 )
               {
                  raise_error( "fourth argument should be stricktly positive" ) ;
               }
               int const iter_max = ARG5->to_int(ct) ;
               if( iter_max <= 0 )
               {
                  raise_error( "fifth argument should be stricktly positive" ) ;
               }
               if( x<eps && n<2 )
               {
                  raise_error( "first argument less than 2 and null second argument" ) ;
               }
               result = MAC::En( (size_t) n, x, eps, err, (size_t) iter_max ) ;
            }
         }
         break ;
      case Ei :
         {
            double const x = ARG1->to_double(ct) ;
            if( x < 0 )
            {
               raise_error( "first argument should be positive" ) ;
            }
            if( ARG2 == 0 )
            {
               result = MAC::Ei( x ) ;
            }
            else
            {
               double const eps = ARG2->to_double(ct) ;
               if( eps < 0 )
               {
                  raise_error( "second argument should be positive" ) ;
               }
               double const err = ARG3->to_double(ct) ;
               if( err <= 0 )
               {
                  raise_error( "third argument should be stricktly positive" ) ;
               }
               int const iter_max = ARG4->to_int(ct) ;
               if( iter_max <= 0 )
               {
                  raise_error( "fourth argument should be stricktly positive" ) ;
               }
               result = MAC::Ei( x, eps, err, (size_t) iter_max ) ;
            }
         }
         break ;
      case Abs :
         result = MAC::abs( ARG1->to_double(ct) ) ;
         break ;
      case Floor :
         result = MAC::floor( ARG1->to_double(ct) ) ;
         break ;
      case Ceil :
         result = MAC::ceil( ARG1->to_double(ct) ) ;
         break ;
      case  RandDbl :
         result = MAC::random_double() ;
         break ;
      default :
         MAC_MathFunctionExp_ERROR::n0( "to_double", name() ) ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
int
MAC_MathFunctionExp:: to_int( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MathFunctionExp:: to_int" ) ;
   MAC_CHECK_PRE( to_int_PRE(ct) ) ;
   
   int result = MAC::max_int() ;

   switch( OP )
   {
      case Abs :
         result = MAC::abs( ARG1->to_int(ct) ) ;
         break ;
      case RandInt :
         result = MAC::rand() ;
         break ;
      default :
         MAC_MathFunctionExp_ERROR::n0( "to_int", name() ) ;
         break ;
   }
   return result ;
}


//----------------------------------------------------------------------
MAC_Data*
MAC_MathFunctionExp:: create_derivative( MAC_Object* a_owner,
                                         MAC_Variable const* var,
                                         MAC_Context const* ct ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MathFunctionExp:: create_derivative" ) ;
   MAC_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   MAC_Data* result = 0 ;
   if( OP==Pow )
   {
      MAC_List * lst = MAC_List::create( 0 ) ;
      lst->append( ARG1->create_clone( lst ) ) ;
      MAC_Data* log = MAC_Expression::create( 0, "log", lst ) ;
      lst->set_owner( log ) ;
      
      lst = MAC_List::create( 0 ) ;
      lst->append( ARG2->create_derivative( lst, var, ct ) ) ;
      lst->append( log ) ; log->set_owner( lst ) ;
      MAC_Data* yprime_log = MAC_Expression::create( 0, "*", lst ) ;
      lst->set_owner( yprime_log ) ;
      
      lst = MAC_List::create( 0 ) ;
      lst->append( yprime_log ) ; yprime_log->set_owner( lst ) ;
      lst->append( create_clone( lst ) ) ;
      MAC_Data* yprime_log_x_pow_y = MAC_Expression::create( 0, "*", lst ) ;
      lst->set_owner( yprime_log_x_pow_y ) ;

      lst = MAC_List::create( 0 ) ;
      lst->append( ARG2->create_clone( lst ) ) ; 
      lst->append( ARG1->create_derivative( lst, var, ct ) ) ;
      MAC_Data* y_xprime = MAC_Expression::create( 0, "*", lst ) ;
      lst->set_owner( y_xprime ) ;
      
      lst = MAC_List::create( 0 ) ;
      lst->append( ARG2->create_clone( lst ) ) ; 
      lst->append( MAC_Double::create( lst, 1.0 ) ) ; 
      MAC_Data* ym1 = MAC_Expression::create( 0, "-", lst ) ;
      lst->set_owner( ym1 ) ;

      lst = MAC_List::create( 0 ) ;
      lst->append( ARG1->create_clone( lst ) ) ; 
      lst->append( ym1 ) ;  ym1->set_owner( lst ) ;
      MAC_Data* x_pow_ym1 = MAC_Expression::create( 0, "pow", lst ) ;
      lst->set_owner( x_pow_ym1 ) ;
      
      lst = MAC_List::create( 0 ) ;
      lst->append( y_xprime ) ; y_xprime->set_owner( lst ) ;
      lst->append( x_pow_ym1 ) ; x_pow_ym1->set_owner( lst ) ;
      MAC_Data* y_xprime_x_pow_ym1 = MAC_Expression::create( 0, "*", lst ) ;
      lst->set_owner( y_xprime_x_pow_ym1 ) ;
      
      lst = MAC_List::create( 0 ) ;
      lst->append( yprime_log_x_pow_y ) ;yprime_log_x_pow_y->set_owner( lst ) ;
      lst->append( y_xprime_x_pow_ym1 ) ;y_xprime_x_pow_ym1->set_owner( lst ) ;
      result = MAC_Expression::create( a_owner, "+", lst ) ;
      lst->set_owner( result ) ;
      
   }
   else
   {
      MAC_Data* composed_derivative = 0 ;
      MAC_Data const* der_variable = ARG1 ;
      MAC_List * lst = MAC_List::create( 0 ) ;
      switch( OP )
      {
         case Sqr :
            lst->append( MAC_Double::create( lst, 2.0 ) ) ;
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "*", lst ) ;
            lst->set_owner( composed_derivative ) ;
         
            break ;
         case Sqrt :
            lst->append( ARG1->create_clone( lst ) ) ;
            lst->append( MAC_Double::create( lst, -0.5 ) ) ;
            composed_derivative = MAC_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;
            lst = MAC_List::create( 0 ) ;
            lst->append( MAC_Double::create( lst, 0.5 ) ) ;
            lst->append( composed_derivative ) ; composed_derivative->set_owner( lst ) ;
            composed_derivative = MAC_Expression::create( 0, "*", lst ) ;
            lst->set_owner( composed_derivative ) ;
         
            break ;
         case Exp :
            composed_derivative = create_clone( 0 ) ;
            lst->destroy() ;
            break ;
         case Log :
            lst->append( MAC_Double::create( lst, 1.0 ) ) ;
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "/", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Sin :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "cos", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Cos :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "sin", lst ) ;
            lst->set_owner( composed_derivative ) ;
            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "unary_minus", lst ) ;
            lst->set_owner( composed_derivative ) ;         
            break ;
         case Tan :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "cos", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( MAC_Double::create( lst, -2.0 ) ) ;
            composed_derivative = MAC_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;         
            break ;
         case ASin :
         case ACos :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( MAC_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "-", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( MAC_Double::create( lst, -0.5 ) ) ;
            composed_derivative = MAC_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;
            if( OP==ASin ) break ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "unary_minus", lst ) ;
            lst->set_owner( composed_derivative ) ;         
            break ;
         case ATan :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( MAC_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "+", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( MAC_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "/", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Sinh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "cosh", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Cosh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "sinh", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Tanh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "cosh", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( MAC_Double::create( lst, -2.0 ) ) ;
            composed_derivative = MAC_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;         
            break ;
         case Erf :
         case Erfc :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "unary_minus", lst ) ;
            lst->set_owner( composed_derivative ) ;

            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "exp", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append(
               MAC_Double::create( lst,
                                   ( OP==Erf ? 1.0 : -1.0 )*2.0/MAC::sqrt( MAC::pi() ) ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "*", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case ASinh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( MAC_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "+", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( MAC_Double::create( lst, -0.5 ) ) ;
            composed_derivative = MAC_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case ACosh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( MAC_Double::create( lst, 1.0 ) ) ;
            composed_derivative = MAC_Expression::create( 0, "-", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( MAC_Double::create( lst, -0.5 ) ) ;
            composed_derivative = MAC_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case ATanh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( MAC_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "-", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = MAC_List::create( 0 ) ;
            lst->append( MAC_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = MAC_Expression::create( 0, "/", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Ei :
            // e(x)/x
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "exp", lst ) ;
            lst->set_owner( composed_derivative ) ;

            lst = MAC_List::create( 0 ) ;
            lst->append( composed_derivative ) ;
            composed_derivative->set_owner( lst ) ;
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = MAC_Expression::create( 0, "/", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case En :
            // n==0 ? -(x+1)*e(-x)/x^2 : -E(n-1,x)
            {
               der_variable = ARG2 ;
               int const n = ARG1->to_int( ct ) ;
               if( n == 0 )
               {
                  lst->append( ARG2->create_clone( lst ) ) ;
                  composed_derivative =
                     MAC_Expression::create( 0, "unary_minus", lst ) ;
                  lst->set_owner( composed_derivative ) ;

                  lst = MAC_List::create( 0 ) ;
                  lst->append( composed_derivative ) ;
                  composed_derivative->set_owner( lst ) ;
                  composed_derivative =
                     MAC_Expression::create( 0, "exp", lst ) ;
                  lst->set_owner( composed_derivative ) ;

                  lst = MAC_List::create( 0 ) ;
                  lst->append( MAC_Double::create( lst, -1. ) ) ;
                  lst->append( ARG2->create_clone( lst ) ) ;
                  MAC_Expression* exp = 
                     MAC_Expression::create( 0, "-", lst ) ;
                  lst->set_owner( exp ) ;

                  lst = MAC_List::create( 0 ) ;
                  lst->append( exp ) ;
                  exp->set_owner( lst ) ;
                  lst->append( composed_derivative ) ;
                  composed_derivative->set_owner( lst ) ;
                  composed_derivative =
                     MAC_Expression::create( 0, "*", lst ) ;
                  lst->set_owner( composed_derivative ) ;

                  lst = MAC_List::create( 0 ) ;
                  lst->append( composed_derivative ) ;
                  composed_derivative->set_owner( lst ) ;
                  lst->append( ARG2->create_clone( lst ) ) ;
                  composed_derivative =
                     MAC_Expression::create( 0, "/", lst ) ;
                  lst->set_owner( composed_derivative ) ;
                  
                  lst = MAC_List::create( 0 ) ;
                  lst->append( composed_derivative ) ;
                  composed_derivative->set_owner( lst ) ;
                  lst->append( ARG2->create_clone( lst ) ) ;
                  composed_derivative =
                     MAC_Expression::create( 0, "/", lst ) ;
                  lst->set_owner( composed_derivative ) ;
               }
               else
               {
                  lst->append( MAC_Int::create( lst, n-1 ) ) ;
                  lst->append( ARG2->create_clone( lst ) ) ;
                  if( ARG3 != 0 )
                  {
                     lst->append( ARG3->create_clone( lst ) ) ;
                     lst->append( ARG4->create_clone( lst ) ) ;
                     lst->append( ARG5->create_clone( lst ) ) ;
                  }
                  composed_derivative = MAC_Expression::create( 0, "En", lst ) ;
                  lst->set_owner( composed_derivative ) ;

                  lst = MAC_List::create( 0 ) ;
                  lst->append( composed_derivative ) ;
                  composed_derivative->set_owner( lst ) ;
                  composed_derivative =
                     MAC_Expression::create( 0, "unary_minus", lst ) ;
                  lst->set_owner( composed_derivative ) ;
               }
            }
            break ;
         case J0 :
         case J1 :
         case Jn :
         case Y0 :
         case Y1 :
         case Yn :
         case Gamma :
         case LGamma :
         case ATan2 :
         case Log10 :
         case IGamma :
         default :
            MAC_MathFunctionExp_ERROR::n0( "create_derivative", name() ) ;
            break ;
      }
      MAC_Data* der = der_variable->create_derivative( 0, var, ct ) ;
      lst = MAC_List::create( 0 ) ;
      lst->append( der ) ; der->set_owner( lst ) ;
      lst->append( composed_derivative ) ; composed_derivative->set_owner( lst ) ;
      result = MAC_Expression::create( a_owner, "*", lst ) ;
      lst->set_owner( result ) ;
   }   
   MAC_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
MAC_MathFunctionExp_ERROR:: n0( std::string const& f_name,
                                std::string const& op_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** MAC_MathFunctionExp::" + f_name +"\n" ;
   mesg += "    operation " + op_name + " not implemented." ;
   MAC_Error::object()->raise_internal( mesg ) ;
}
