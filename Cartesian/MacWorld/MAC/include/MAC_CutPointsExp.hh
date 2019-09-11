#ifndef MAC_CUT_POINTS_EXP_HH
#define MAC_CUT_POINTS_EXP_HH

#include <MAC_Expression.hh>

/*
---
name     : middle_point
argument : Double, DoubleVector
type     : Double

middle_point(x,v) return the middle point of the interval defined by two
successive elements of v containing x.

Example:
   middle_point( 2.1, (1.,2.,3.,4.) ): value is 2.5

---
name     : middle_points
argument : DoubleVector[,DoubleVector]
type     : DoubleVector

middle_points(v) builds a vector containaing the middle points of
the intervals defined by two successive elements of v.

Example:
   middle_points( (1.,2.,3.,4.) ): value is (1.5,2.5,3.5)

middle_points(x_table,v) builds a vector containaing the middle points of
the intervals defined by two successive elements of v containing the elements
of x_table.

Example:
   middle_points( (1.1, 3.9), (1.,2.,3.,4.) ): value is (1.5,3.5)
   
---
name     : x_cut_points
argument : < x_table, y, [z] >
type     : doubleArray2D

x_cut_points(x_table,y) builds size(x_table)-1 points of 2 dimensions,
the i-th of which is (0.5*(x_table(i)+xtable(i+1)),y).
x_cut_points(x_table,y,z) builds size(x_table)-1 points of 3 dimensions,
the i-th of which is (0.5*(x_table(i)+xtable(i+1)),y,z).

example:
   x_cut_points((1.,2.,3.),3.9):
                 value is ( (1.5,3.9),(2.5,3.9) )
   x_cut_points((1.,2.,3.,5.),3.9,6.):
                 value is ( (1.5,3.9,6.),(2.5,3.9,6.),(4.,3.9,6.) )

---
name     : y_cut_points
argument : < x, y_table, [z] >
type     : doubleArray2D

y_cut_points(x,y_table) builds a set of size(y_table)-1 points of 2 dimensions,
the i-th of which is (x,0.5*(y_table(i)+ytable(i+1))).
y_cut_points(x,y_table,z) builds a set of size(y_table)-1 points of 3 dimensions,
the i-th of which is (x,0.5*(y_table(i)+ytable(i+1)),z).

example:
   y_cut_points(3.9,(1.,2.,3.),):
                 value is ( (3.9,1.5),(3.9,2.5) )
   y_cut_points(3.9,(1.,2.,3.,5.),6.):
                 value is ( (3.9,1.5,6.),(3.9,2.5,6.),(3.9,4.,6.) )

---
name     : z_cut_points
argument : < x, y, z_table >
type     : doubleArray2D

z_cut_points(x,y,z_table) builds a set of size(z_table)-1 points of 3 dimensions,
the i-th of which is (x,y,0.5*(z_table(i)+ztable(i+1))).

example:
   z_cut_points(3.9,6.,(1.,2.,3.,5.)):
                 value is ( (3.9,6.,1.5),(3.9,6.,2.5),(3.9,6.,4.) )
                 
PUBLISHED
*/

class MAC_CutPointsExp : public MAC_Expression
{
   public: //----------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value

      virtual double to_double( MAC_Context const* ct ) const ;
      
      virtual doubleVector const& to_double_vector(
                                        MAC_Context const* ct ) const ;
      
      virtual doubleArray2D const& to_double_array2D(
                                        MAC_Context const* ct ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      MAC_CutPointsExp( void ) ;
     ~MAC_CutPointsExp( void ) ;
      MAC_CutPointsExp( MAC_CutPointsExp const& other ) ;
      MAC_CutPointsExp& operator=( MAC_CutPointsExp const& other ) ;

      enum IS_CutOp { xcut, ycut, zcut, mid_point, mid_points } ;
      
      MAC_CutPointsExp( MAC_Object* a_owner,
                        std::string const& a_name,
                        MAC_Sequence const* argument_list,
                        IS_CutOp a_op ) ;

      void check_table( doubleVector const& verts_table ) const ;
      
      double m_pt( double x, doubleVector const& verts_table ) const ;
      
      void build_coords( doubleVector const& verts_table,
                         doubleVector& coords_table ) const ;
      
   //-- Plug in
      
      MAC_CutPointsExp( std::string const& a_name, IS_CutOp a_op ) ;

      virtual MAC_CutPointsExp* create_replica( 
                       MAC_Object* a_owner,
                       MAC_Sequence const* argument_list ) const ;
      
  //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments(
                      MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static MAC_CutPointsExp const* PROTOTYPE_MP ;
      static MAC_CutPointsExp const* PROTOTYPE_MPS ;
      static MAC_CutPointsExp const* PROTOTYPE_X ;
      static MAC_CutPointsExp const* PROTOTYPE_Y ;
      static MAC_CutPointsExp const* PROTOTYPE_Z ;
          
   //-- Attributes
      
      IS_CutOp const OP ;
} ;

#endif
