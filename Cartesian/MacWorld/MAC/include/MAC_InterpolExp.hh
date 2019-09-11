#ifndef MAC_INTERPOL_EXP_HH
#define MAC_INTERPOL_EXP_HH

#include <MAC_Expression.hh>
#include <doubleVector.hh>

/*
Expressions interpolating between given values.

---
name      : interpol
arguments : DoubleVector, DoubleVector, Double 
type      : double

interpol( x_values, fx_values, x ) computes at x the linear interpolation 
between the data points (x_values, fx_values).
More precisely:
   * x_values and fx_values are two vectors of the same size, say N
   * the elements of x_values are ordered by increasing values
        x_values(0) < x_values(1) < ... < x_values(N-1)
   * if x is smaller than x_values(0) then return fx_values(0)
     else if x is greater than x_values(N-1) then return fx_values(N-1)
     else find j such that x is between x_values(j) and x_values(j+1)
          return the linear interpolation between 
                 fx_values(j) and fx_values(j+1)

example :
   $DV_Xval = < 1. 2. 5. >
   $DV_Yval = < 1. 4. 2. >
   interpol( $DV_Xval, $DV_Yval, 3. ) : value is y given by 
                                        (y-4.0)/(3.0-2.0)=(2.0-4.0)/(5.0-2.0)

---
name      : interpol
arguments : String, Double 
type      : double

interpol( filename, x ) computes at x the linear interpolation between the
data points stored in file of name filename.
the value is identical as that of interpol( x_values, fx_values, x )
(see above) where filename stores in sequence x_values(0), fx_values(0), 
x_values(1), fx_values(1), ... , x_values(N-1), fx_values(N-1)

example :
   file "values.txt":
      1. 1.
      2. 4.
      5. 2.
   interpol( join( this_file_dir(), "values.txt" ), 3. ) : same value as above
      ie y given by (y-4.0)/(3.0-2.0)=(2.0-4.0)/(5.0-2.0)
   
PUBLISHED
*/

class MAC_InterpolExp : public MAC_Expression
{
   public: //-------------------------------------------------------
      
   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( MAC_Context const* ct = 0 ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------

      MAC_InterpolExp( void ) ;
     ~MAC_InterpolExp( void ) ;
      MAC_InterpolExp( MAC_InterpolExp const& other ) ;
      MAC_InterpolExp& operator=( MAC_InterpolExp const& other ) ;

      enum MAC_InterpolOp { lin_inter_1D } ;

      MAC_InterpolExp( MAC_Object* a_owner,
                       std::string const& a_name,
                       MAC_Sequence const* argument_list,
                       MAC_InterpolOp a_op ) ;

      void read_tables_1D( std::string const& filename,
                           doubleVector& X_table,
                           doubleVector& FX_table ) const ;
      
      double linear_interpol_1D( doubleVector const& X_table,
                                 doubleVector const& FX_table,
                                 double x ) const ;

      void check_tables_1D( doubleVector const& X_table,
                            doubleVector const& FX_table ) const ;

   //-- Plug in

      MAC_InterpolExp( std::string const& a_name,
                       MAC_InterpolOp a_op ) ;

      virtual MAC_Expression* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static MAC_InterpolExp const* PROTOTYPE_LIN_1D ;
      
   //-- Attributes

      MAC_InterpolOp const OP ;
      bool FROM_FILE ;
      mutable bool X1_IS_SET ;
      mutable doubleVector X1 ;
      mutable bool F_IS_SET ;
      mutable doubleVector FX1 ;
      mutable bool CHECK ;
} ;

#endif
