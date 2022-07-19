#ifndef _BBOX_HH_
#define _BBOX_HH_

#include "Vector3.hh"
#include "Point3.hh"
using namespace solid;

class BBox;
ostream& operator << ( ostream& f, BBox const& B ); 


/** @brief The class BBox.

    Bounding box oriented along the axis of the world reference frame (no
    orientation). From GJK Engine - A Fast and 
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.
    
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class BBox 
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    BBox(); 
  
    /** @brief Constructor with 2 corners as inputs, the 1st point with the 
    lowest coordinates and the 2nd point with the largest coordinates 
    @param min point with the lowest coordinates
    @param max point with the largest coordinates */
    BBox( Point3 const& min, Point3 const& max );

    /** @brief Copy constructor 
    @param bbox_ La boite de reference */
    BBox( BBox const& bbox_ );
  
    /** @brief Destructeur */
    ~BBox();
    //@}

    /**@name Methods */
    //@{
    /** @brief Sets the bounding box to the intersection of the 2 bounding boxes
    a and b. If there is no intersection, leaves the bounding box unchanged
    @param a 1st bounding box
    @param b 2nd bounding box */
    void closest( BBox const& a, BBox const& b );

    /** @brief Sets the bounding box to the union of the 2 bounding boxes
    a and b
    @param a 1st bounding box
    @param b 2nd bounding box */
    void enclose( BBox const& a, BBox const& b ); 

    /** @brief Returns the bounding box center */
    Point3 const& getCenter() const;

    /** @brief Returns the box extensions wrt its center, i.e., its edge half
    length vector */
    Vector3 const& getExtent() const;

    /** @brief Returns the ith minimum coordinate
    @param i coordinate index */
    double getLower( int i ) const;

    /** @brief Returns the ith maximum coordinate
    @param i coordinate index */
    double getUpper( int i ) const;

    /** @brief Extends the bounding box to a point p if p is outside the box
    @param p point */
    void include( Point3 const& p );
  
    /** @brief Sets the bounding box to the union of itself and another bounding
    box
    @param b the other bounding box */
    void include ( BBox const& b );

    /** @brief Returns whether a cubic box defined by its center and its half 
    edge length intersects the bounding box
    @param p center of the cubic box
    @param halfEdgeLength half-edge length of the cubic box */
    bool InZone( Point3 const* p, double halfEdgeLength ) const;

    /** @brief Returns whether a box defined by its center and its half 
    edge lengths intersects the bounding box 
    @param p center of the cubic box
    @param halfEdgeLength_X half-edge length of the box in the X direction 
    @param halfEdgeLength_Y half-edge length of the box in the Y direction     
    @param halfEdgeLength_Z half-edge length of the box in the Z direction */
    bool InZone( Point3 const* p, double halfEdgeLength_X, 
    	double halfEdgeLength_Y, double halfEdgeLength_Z ) const;

    /** @brief Returns the direction of longest edge */
    int longestAxis() const;

    /** @brief Sets the bounding box center 
    @param p new center */
    void setCenter( Point3 const& p );

    /** @brief Sets the bounding box to an empty bounding box. This is done by
    assigning minus infinity extensions */
    void setEmpty(); 

    /** @brief Sets the bounding boxes half lengths 
    @param v new extension vector */
    void setExtent( Vector3 const& v );

    /** @brief Sets the box dimensions using the point with the 
    lowest coordinates and the point with the largest coordinates 
    @param min point with the lowest coordinates
    @param max point with the largest coordinates */
    void setValue( Point3 const& min, Point3 const& max );

    /** @brief Returns the largest half length of the bounding box */
    double size() const ;
  
    /** @brief Debugging method
    @param s debugging message to be printed on the default error output cerr */
    void debug( char const* s) const;
    //@}
  

    /** @name Friend methods */
    //@{
    /** @brief Returns whether the 2 bounding boxes a and b intersect
    @param a 1st bounding box
    @param b 2nd bounding box */
    friend bool intersect( BBox const& a, BBox const& b );
  
    /** @brief Output operator
    @param f output stream
    @param B BBox object */
    friend ostream& operator << ( ostream& f, BBox const& B );   
    //@}


    /**@name Operators */
    //@{
    /** @brief Equal operator to another BBox object
    @param rhs the other BBox object */
    BBox& operator = ( BBox const& rhs );
    //@}


  private:
    /** @name Parameters */
    //@{  
    Point3 m_center; /**< bounding box center */
    Vector3 m_extent; /**< bounding box half lengths */
    //@}
};

#endif


