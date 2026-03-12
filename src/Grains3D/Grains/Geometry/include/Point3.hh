#ifndef _POINT3_HH_
#define _POINT3_HH_

#include "Group3.hh"
#include "Vector3.hh"


namespace solid
{
  /** @brief The class Point3.
  
  Point in a 3D space. From GJK Engine - A Fast and Robust GJK
  Implementation, Copyright (C) 1998  Gino van den Bergen. 
      
  @author Institut Francais du Petrole - 1999 - Creation 
  @author Institut Francais du Petrole - 2000 - Modification
  @author A.WACHS - 2009 - Modification */
  // ==========================================================================
  class Point3 : public Group3
  {
    public:
      /**@name Constructors */
      //@{
      /** @brief Default constructor 
      @param def value of all 3 components */
      inline Point3( double def = 0. ) : Group3( def ) {}

      /** @brief Constructor with 3 components as inputs 
      @param x 1st component
      @param y 2nd component
      @param z 3rd component*/
      inline Point3( double x, double y, double z ) : Group3( x, y, z ) {}

      /** @brief Copy constructor
      @param g copied Point3 object */
      inline Point3( Point3 const& g ) : Group3( g ) {}
      
      /** @brief Copy constructor
      @param g copied Group3 object */      
      inline Point3( Group3 const& g ) : Group3( g ) {}

      /** @brief Destructor */
      inline ~Point3() {}
      //@}


      /**@name Methods */
      //@{
      /** @brief Adds the same value to all 3 components
      @param dist motion value */
      void Move( double dist );
    
      /** @brief Adds values to the 3 components using a 3-component array
      @param dist 3-component array with motion in each direction */
      void Move( double const* dist );
    
      /** @brief Adds values to the 3 components using 3 scalars
      @param distX motion in x
      @param distY motion in y
      @param distZ motion in z */
      void Move( double distX, double distY, double distZ );

      /** @brief Distance between 2 points of type Point3
      @param point 2nd point */    
      inline double DistanceTo( Point3 const& point ) const
      {
        double a = m_comp[0] - point[0];
        double b = m_comp[1] - point[1];
        double c = m_comp[2] - point[2];
        return ( sqrt( a*a + b*b + c*c ) );
      }
    
      /** @brief Distance between the point of type Point3 and a point defined 
      by a 3-component array
      @param pp 2nd point defined by a 3-component array */
      inline double DistanceTo( double const* pp ) const
      {
        double a = m_comp[0] - pp[0];
        double b = m_comp[1] - pp[1];
        double c = m_comp[2] - pp[2];
        return ( sqrt( a*a + b*b + c*c ) );
      }
    
      /** @brief Distance between the point of type Point3 and a point defined 
      by 3 scalars
      @param x Position en X
      @param y Position en Y
      @param z Position en Z */
      inline double DistanceTo( double x, double y, double z ) const
      {
        double a = m_comp[0] - x;
        double b = m_comp[1] - y;
        double c = m_comp[2] - z;
        return ( sqrt( a*a + b*b + c*c ) );
      }
      //@}


      /**@name Operators */
      //@{      
      /** @brief Equal operator to another Group3 object
      @param g2 the other Group3 object */
      inline Point3& operator = ( Point3 const& g2 )
      {
        if ( &g2 != this )
        { m_comp[X] = g2.m_comp[X]; m_comp[Y] = g2.m_comp[Y]; m_comp[Z] = g2.m_comp[Z]; }
        return ( *this );
      }
      //@}      
  };
  
  static Point3 Point3Null; /**< (0.,0.,0.)  */
} // namespace solid

#endif

