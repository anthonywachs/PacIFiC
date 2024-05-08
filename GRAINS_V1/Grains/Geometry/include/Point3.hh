#ifndef _POINT3_HH_
#define _POINT3_HH_

#include "Group3.hh"
#include "Vector3.hh"


namespace solid
{
  /** @brief The class Point3.
  
  Point in a 3D space. From GJK Engine - A Fast and Robust GJK
  Implementation, Copyright (C) 1998  Gino van den Bergen. 
      
  @author G.FERRER - Institut Francais du Petrole - 1999 - Creation 
  @author F.PRADEL - Institut Francais du Petrole - 2000 - Modification
  @author D.PETIT  - 2000 - Modification 
  @author A.WACHS  - 2009 - Modification 
  @author A.WACHS  - 2019 - Modification */
  // ==========================================================================
  class Point3 : public Group3
  {
    public:
      /**@name Constructors */
      //@{
      /** @brief Default constructor 
      @param def value of all 3 components */
      Point3( double def = 0. );

      /** @brief Constructor with 3 components as inputs 
      @param x 1st component
      @param y 2nd component
      @param z 3rd component*/
      Point3( double x, double y, double z );

      /** @brief Copy constructor
      @param g copied Point3 object */
      Point3( Point3 const& g );
      
      /** @brief Copy constructor
      @param g copied Group3 object */      
      Point3( Group3 const& g );      

      /** @brief Destructor */
      ~Point3();
      //@}


      /**@name Methods */
      //@{
      /** @brief Adds the same value to all 3 components
      @param dist displacement value */
      void Move( double dist );
    
      /** @brief Adds values to the 3 components using a 3-component array
      @param dist 3-component array with displacement in each direction */
      void Move( double const* dist );
    
      /** @brief Adds values to the 3 components using 3 scalars
      @param distX displacement in x
      @param distY displacement in y
      @param distZ displacement in z */
      void Move( double distX, double distY, double distZ );

      /** @brief Distance between 2 points of type Point3
      @param point 2nd point */    
      double DistanceTo( Point3 const& point ) const;
    
      /** @brief Distance between the point of type Point3 and a point defined 
      by a 3-component array
      @param pp 2nd point defined by a 3-component array */
      double DistanceTo( double const* pp ) const;
    
      /** @brief Distance between the point of type Point3 and a point defined 
      by 3 scalars
      @param x Position en X
      @param y Position en Y
      @param z Position en Z */
      double DistanceTo( double x, double y, double z ) const;
      //@}


      /**@name Operators */
      //@{      
      /** @brief Equal operator to another Group3 object
      @param g2 the other Group3 object */
      Point3& operator = ( Point3 const& g2 );  
      //@}      
  };
  
  static Point3 Point3Null; /**< (0.,0.,0.)  */
} // namespace solid

#endif

