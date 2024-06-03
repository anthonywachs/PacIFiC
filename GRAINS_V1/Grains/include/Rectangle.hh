#ifndef _RECTANGLE_HH_
#define _RECTANGLE_HH_

#include "Convex.hh"
#include "ReaderXML.hh"

class Transform;
/** @brief The class Rectangle.

    Convex with a rectangular shape.

    @author A.YAZDANI - 2019 - Creation */
// ============================================================================
class Rectangle : public Convex
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with edge length as input parameters
    @param LX edge length in x
    @param LY edge length in y */
    Rectangle( double LX = 0., double LY = 0. );

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    Rectangle( istream& fileIn );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    Rectangle( DOMNode* root );

    /** @brief Destructor */
    ~Rectangle();
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the box */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns a vector of points describing the envelope of the box */
    vector<Point3> getEnvelope() const;

    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the point indices */
    vector<vector<int> > const* getFaces() const;

    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape. Here returns the code 444 */
    int getNbCorners() const;

    /** @brief Returns the area */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Rectangle support function, returns the support point P, i.e. the
    point on the surface of the box that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;

    /** @brief Returns the number of points to write the rectangle in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the rectangle
    in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the rectangle in a
    Paraview format
    @param f output stream
    @param transform geometric transformation
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns a list of points describing the rectangle in a
    Paraview format
    @param transform geometric transformation
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the rectangle in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	    int& last_offset ) const;

    /** Returns whether the convex shape is a rectangle */
    bool isRectangle() const;

    /** @ brief Returns whether a point lies inside the rectangle
    @param pt point */
    bool isIn( Point3 const& pt ) const;
    
    /** @brief Performs advanced comparison of the two rectangles and returns
    whether they match
    @param other the other rectangle */
    bool equalType_level2( Convex const* other ) const;     
    //@}


  protected:
    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference rectangle,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
    //@}


    /** @name Parameters*/
    //@{
    double m_LX; /**< legnth of the X edge */
    double m_LY; /**< legnth of the Y edge */
    vector<Point3> m_corners; /**< vector of the 4 corners/vertices */
    //@}


  private:
    /**@name Methods */
    //@{
    /** @brief Sets the corner/vertex coordinates */
    void setCorners();
};

#endif
