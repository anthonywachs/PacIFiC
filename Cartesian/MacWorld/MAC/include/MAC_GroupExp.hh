#ifndef MAC_GROUP_EXP_HH
#define MAC_GROUP_EXP_HH

#include <MAC_Expression.hh>
#include <doubleVector.hh>

/* 
Expressions grouping items according a given subdivision of a given set.

---
name      : unit_sort
arguments : Double, Double, Double, Int[, Bool]
type      : Int

unit_sort( x, x1, x2, N ) splits [x1,x2] into N intervals of equal
diameters and returns the index of the interval containing x. The lower
bound of indices is zero (included) and the upper bound is the
number of intervals (excluded).

example :  
   unit_sort( 0.0, 0.0, 1.0, 2 ) : value is 0 (first interval)
   unit_sort( 0.9, 0.0, 1.0, 2 ) : value is 1 (second interval)

The last optional argument (default is false) shifts the intervals of
half an interval.

example:
   unit_sort( 0.0000, 0.0, 1.0, 2, true ) : value is 0 (first interval)
   unit_sort( 0.2499, 0.0, 1.0, 2, true ) : value is 0 (first interval)
   unit_sort( 0.2501, 0.0, 1.0, 2, true ) : value is 1 (second interval)
   unit_sort( 0.5000, 0.0, 1.0, 2, true ) : value is 1 (second interval)
   unit_sort( 0.7499, 0.0, 1.0, 2, true ) : value is 1 (second interval)
   unit_sort( 0.7501, 0.0, 1.0, 2, true ) : value is 0 (first interval)
   unit_sort( 1.0000, 0.0, 1.0, 2, true ) : value is 0 (first interval)
   
---
name      : segm_sort
arguments : Double, doubleVector, Int[, Bool]
type      : Int

segm_sort( x, vec, N ) splits vec into N intervals of equal
number of values and returns the index of the interval containing x. The lower
bound of indices is zero (included) and the upper bound is the
number of intervals (excluded).

example :  
   segm_sort( 0.0, <0. 0.5 0.75 1.>, 2 ) : value is 0 (first interval)
   segm_sort( 0.9, <0. 0.5 0.75 1.>, 2 ) : value is 1 (second interval)
 
The last optional argument (default is false) shifts the intervals of
half an interval.
   
---
name      : segm2D_sort
arguments : 2D doubleVector, doubleVector1, Int1, doubleVector2, Int2[, Bool]
type      : Int

segm2D_sort( x, vec1, N1, vec2, N2 ) splits vec1 into N1 intervals of equal
number of values, vec2 into N2 intervals of equal number of values and returns
the index of the interval containing x.
The lower bound of indices is zero (included) and the upper bound is N1*N2
(excluded).
    with i1 = segm_sort(x(0),vec1, N1) (in [0...N1-1])
         i2 = segm_sort(x(1),vec2, N2) (in [0...N2-1])
    returns : i2*N1+i1 (in [0...N1*N2-1])
  
The last optional argument (default is false) shifts the intervals of
half an interval.

---
name      : segm3D_sort
arguments : 3D doubleVector, doubleVector1, Int1, doubleVector2, Int2, doubleVector3, Int3[, Bool]
type      : Int

segm3D_sort( x, vec1, N1, vec2, N2, vec3, N3 ) splits vec1 into N1 intervals
of equal number of values, vec2 into N2 intervals of equal number of values,
vec3 into N3 intervals of equal number of values, and returns the index of the
interval containing x.
The lower bound of indices is zero (included) and the upper bound is N1*N2*N3
(excluded).
    with i1 = segm_sort(x(0),vec1, N1) (in [0...N1-1])
         i2 = segm_sort(x(1),vec2, N2) (in [0...N2-1])
         i3 = segm_sort(x(2),vec3, N3) (in [0...N3-1])
    returns : i3*N1*N2+i2*N1+i1 (in [0...N1*N2*N3-1])

The last optional argument (default is false) shifts the intervals of
half an interval.
    
---

PUBLISHED
*/

class MAC_GroupExp : public MAC_Expression
{
   public: //--------------------------------------------------------------

   //-- Instance delivery and initialization

      /* Enable optimized evaluation:
          - some inner tables are set at the first called to `::to_int' function
          - WARNING: tables vec, vec1, vec2, vec3, and sizes N, N1, N2, N3
                     are supposed to be fixed between two calls of the function
                     (only x is varying).
          - then:
               initialization is in O( size(vec) ), performed only onces
               `::to_int' function is then in O(log(N)).
      */
      static void set_optimized_evaluation( void ) ;

      /* Disable optimized evaluation:
          - `::to_int' function is then in O( size(vec)+log(N) ).
      */
      static void unset_optimized_evaluation( void ) ;
      
   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual int to_int( MAC_Context const* ct ) const ;
      
         
   protected: //----------------------------------------------------------
            
   private: //------------------------------------------------------------

      MAC_GroupExp( void ) ;
     ~MAC_GroupExp( void ) ;
      MAC_GroupExp( MAC_GroupExp const& other ) ;
      MAC_GroupExp& operator=( MAC_GroupExp const& other ) ;

      enum GroupOp { unit_sort, segm_sort, segm2D_sort, segm3D_sort } ;
      
      MAC_GroupExp( MAC_Object* a_owner,
	            std::string const& a_name,
		    MAC_Sequence const* argument_list,
                    GroupOp a_op ) ;

      void initialize( std::string const& v_arg, doubleVector const& v,
                       std::string const& n_arg, int n,
                       doubleVector& x_table ) const ;

      int index( double x, doubleVector const& x_table, bool shift ) const ;
      
   //-- Plug in

      MAC_GroupExp( std::string const& a_name, GroupOp a_op ) ;

      virtual MAC_GroupExp* create_replica( 
                                  MAC_Object * a_owner,
				  MAC_Sequence const* argument_list ) const ;

   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static bool OPT_EVAL ;
      
      static MAC_GroupExp const* PROTOTYPE_UNIT_SORT ;
      static MAC_GroupExp const* PROTOTYPE_SEGM_SORT ;
      static MAC_GroupExp const* PROTOTYPE_SEGM2D_SORT ;
      static MAC_GroupExp const* PROTOTYPE_SEGM3D_SORT ;
            
   //-- Attributes
      
      GroupOp const OP ;

      mutable bool INITIALIZED ;
      mutable doubleVector X ;
      mutable doubleVector Y ;
      mutable doubleVector Z ;
} ;

#endif
