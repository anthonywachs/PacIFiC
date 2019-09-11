#ifndef MAC_MATH_FUNCTION_HH
#define MAC_MATH_FUNCTION_HH

#include <MAC_Expression.hh>

/* 
Expressions defined by the usual mathematical functions :
      sqr, sqrt, pow, exp, log, log10,
      sin, cos, tan, asin, acos, atan, atan2,
      sinh, cosh, tanh, asinh, acosh, atanh
      j0, j1, jn, y0, y1, yn,
      gamma, lgamma, erf, erfc, incomplete_gamma, En, Ei,
      abs, floor, ceil
and by the comparison between double values :
      double_equality

PUBLISHED
*/

class MAC_MathFunctionExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value

      virtual bool to_bool( MAC_Context const* ct = 0 ) const ;
      
      virtual double to_double( MAC_Context const* ct = 0 ) const ;

      virtual int to_int( MAC_Context const* ct = 0 ) const ;
      
   //-- Formal calculus

      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct ) const ;

   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------

      MAC_MathFunctionExp( void ) ;
     ~MAC_MathFunctionExp( void ) ;
      MAC_MathFunctionExp( MAC_MathFunctionExp const& other ) ;
      MAC_MathFunctionExp& operator=( MAC_MathFunctionExp const& other ) ;

      enum MathOp { Sqr, Sqrt, Pow, Exp, Log, Log10,
                    Sin, Cos, Tan, ASin, ACos, ATan, ATan2,
                    Sinh, Cosh, Tanh, ASinh, ACosh, ATanh,
                    J0, J1, Jn, Y0, Y1, Yn,
                    Gamma, LGamma, Erf, Erfc,
                    IGamma, En, Ei, 
                    Abs, Floor, Ceil,
                    DblEq,
                    RandDbl, RandInt } ;
            
      MAC_MathFunctionExp( MAC_Object* a_owner,
                           std::string const& a_name,
                           MAC_Sequence const* argument_list,
                           MathOp a_op ) ;

      MAC_Data const* alternative_result( void ) const ;
      
   //-- Plug in

      MAC_MathFunctionExp( std::string const& a_name, MathOp a_op ) ;

      virtual MAC_MathFunctionExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;

      
   //-- Class attributes

      static MAC_MathFunctionExp const* PROTOTYPE_Sqr ;
      static MAC_MathFunctionExp const* PROTOTYPE_Sqrt ;
      static MAC_MathFunctionExp const* PROTOTYPE_Pow ;
      static MAC_MathFunctionExp const* PROTOTYPE_Exp ;
      static MAC_MathFunctionExp const* PROTOTYPE_Log ;
      static MAC_MathFunctionExp const* PROTOTYPE_Log10 ;
      static MAC_MathFunctionExp const* PROTOTYPE_Sin ;
      static MAC_MathFunctionExp const* PROTOTYPE_Cos ;
      static MAC_MathFunctionExp const* PROTOTYPE_Tan ;
      static MAC_MathFunctionExp const* PROTOTYPE_ASin ;
      static MAC_MathFunctionExp const* PROTOTYPE_ACos ;
      static MAC_MathFunctionExp const* PROTOTYPE_ATan ;
      static MAC_MathFunctionExp const* PROTOTYPE_ATan2 ;
      static MAC_MathFunctionExp const* PROTOTYPE_Sinh ;
      static MAC_MathFunctionExp const* PROTOTYPE_Cosh ;
      static MAC_MathFunctionExp const* PROTOTYPE_Tanh ;
      static MAC_MathFunctionExp const* PROTOTYPE_ASinh ;
      static MAC_MathFunctionExp const* PROTOTYPE_ACosh ;
      static MAC_MathFunctionExp const* PROTOTYPE_ATanh ;
      static MAC_MathFunctionExp const* PROTOTYPE_J0 ;
      static MAC_MathFunctionExp const* PROTOTYPE_J1 ;
      static MAC_MathFunctionExp const* PROTOTYPE_Jn ;
      static MAC_MathFunctionExp const* PROTOTYPE_Y0 ;
      static MAC_MathFunctionExp const* PROTOTYPE_Y1 ;
      static MAC_MathFunctionExp const* PROTOTYPE_Yn ;
      static MAC_MathFunctionExp const* PROTOTYPE_Gamma ;
      static MAC_MathFunctionExp const* PROTOTYPE_LGamma ;
      static MAC_MathFunctionExp const* PROTOTYPE_Erf ;
      static MAC_MathFunctionExp const* PROTOTYPE_Erfc ;
      static MAC_MathFunctionExp const* PROTOTYPE_IGamma ;
      static MAC_MathFunctionExp const* PROTOTYPE_En ;
      static MAC_MathFunctionExp const* PROTOTYPE_Ei ;
      static MAC_MathFunctionExp const* PROTOTYPE_ABS ;
      static MAC_MathFunctionExp const* PROTOTYPE_Floor ;
      static MAC_MathFunctionExp const* PROTOTYPE_Ceil ;
      static MAC_MathFunctionExp const* PROTOTYPE_DblEq ;
      static MAC_MathFunctionExp const* PROTOTYPE_RandDbl ;
      static MAC_MathFunctionExp const* PROTOTYPE_RandInt ;
      
   //-- Attributes
      
      MathOp const OP ;
      MAC_Data const* const ARG1 ;
      MAC_Data const* const ARG2 ;
      MAC_Data const* const ARG3 ;
      MAC_Data const* const ARG4 ;
      MAC_Data const* const ARG5 ;
} ;

#endif
