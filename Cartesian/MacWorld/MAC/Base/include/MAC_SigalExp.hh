#ifndef MAC_SIGAL_EXP_HH
#define MAC_SIGAL_EXP_HH

#include <MAC_Expression.hh>

#include <doubleVector.hh>
#include <intVector.hh>
#include <stringVector.hh>

/* Operator to form vector from sigal-like expressions.
   These operators are :
   
---
name : regular_vector

* Int case
   type      : IntVector
   arguments : Int, Int, Int

* Double case
   type      : DoubleVector
   arguments : Double, Double, Double

regular_vector( x, n, y ) returns a vector of (n+1) elements such that:
   - the first element is x ;
   - the last element is y ;
   - the intervals between two successive elements all have the same length
     (n being the number of these intervals).

---
name      : stretched_vector
arguments : Double, Double, Double, Double
type      : DoubleVector

stretched_vector( x, dx, dy, y ) returns a vector such that, as much as 
possible:
   - the first element is x ;
   - the last element is y ;
   - the interval between the first two elements has the length dx ;
   - the interval between the last two elements has the length dy.

---
name      : geometric_sequence
arguments : Double, Double, Int
type      : DoubleVector

geometric_sequence( x0, r, n ) returns the first (n+1) elements
of a geometric sequence of first term x0 and common ratio r, ie
the vector such that:
   - the first element is x0 ;
   - two successive elements x_{i} and x_{i+1} are related by 
        x_{i+1} = r * x_{i}
   
example :
   geometric_sequence( 1.0, 2.0, 4 ) = < 1.0 2.0 4.0 8.0 16.0 >

PUBLISHED
 */

class MAC_SigalExp : public MAC_Expression
{
   public: //-------------------------------------------------------
      
   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual stringVector const& to_string_vector( 
                                          MAC_Context const* ct = 0 ) const ;

      virtual intVector const& to_int_vector( 
                                          MAC_Context const* ct = 0 ) const ;

      virtual doubleVector const& to_double_vector( 
                                          MAC_Context const* ct = 0 ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------
      
      MAC_SigalExp( void ) ;
     ~MAC_SigalExp( void ) ;
      MAC_SigalExp( MAC_SigalExp const& other ) ;
      MAC_SigalExp& operator=( MAC_SigalExp const& other ) ;

      enum SigalOperator { Regu, Ratio, Concat, Geom } ;
      
      MAC_SigalExp( MAC_Object* a_owner,
                    SigalOperator a_op,
                    std::string const& a_name,
                    MAC_Sequence const* argument_list ) ;

      void realize( MAC_Context const* ct ) const ;

      //   The following relation gives `dx0' function of r and n:
      //
      //     ZL = | `x0'-`x1' |
      //
      //           n-1
      //     dx  * SUM r^i = ZL <=>  dx  (1-r^n)/(1-r) = ZL (r <> 1)
      //           i=0
      //
      //     dx  * r^(n-1) = dx-last
      //
      //   n is calculated using the optimized ratio :
      //
      //     r_opt = (ZL-`dx0')/(ZL-`dx1')
      //
      //     => n = 1 + log(`dx1'/`dx0') / log(r_opt)
      //          + 0.5 (to get the nearest integer value)
      //
      //   r is calculated in minimising the following function :
      //
      //     g(r) = (1.-`dx0'/dx(r))^2 + (r^(n-1)-`dx1'/dx(r))^2
      //          = (1.-A1*u(r))^2 + (r^(n-1)-An*u(r))^2
      //
      //     A1=`dx0'/ZL   An=`dx1'/ZL   u=ZL/dx=1/A=SUM r^i
      //
      //   The derivated function of u :
      //
      //     du(r)= SUM i r^(i-1)
      //
      //   The derivated function of G is :
      //
      //     F(r) = - 2 A1 (1 - A1 u) du
      //            + 2 (r^(n-1) - A2 u) ((n-1) r^(n-2) - A2 du)
      //
      //   We want F(r)=0. This relation is obtained through a Newton
      //   method, the derivative of F being calculated numerically.      
      void build_stretched_vector( double x0, double dx0,
                                   double x1, double dx1,
                                   doubleVector& result ) const ;
      

      void tageom(double dx1, double dxn, double zl,
                  double& dx, size_t& n, double& r ) const ;

      double rtaf( double a1, double an, size_t n, double r ) const ;
      
   //-- Plug in

      MAC_SigalExp( SigalOperator a_op, std::string const& a_name ) ;

      virtual MAC_SigalExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics
      
      virtual std::string const& usage( void ) const ;
      
      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes
            
      static MAC_SigalExp const* PROTOTYPE_Regu ;
      static MAC_SigalExp const* PROTOTYPE_Ratio ;
      static MAC_SigalExp const* PROTOTYPE_Concat ;
      static MAC_SigalExp const* PROTOTYPE_GeomSeq ;

   //-- Attributes      

      SigalOperator const OP ;
      mutable intVector RESULT_I ;
      mutable doubleVector RESULT_D ;
      mutable stringVector RESULT_S ;
} ;

#endif
