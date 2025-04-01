#ifndef _TRAPEZOIDALPRISM_HH_
#define _TRAPEZOIDALPRISM_HH_

#include "Convex.hh"
#include "ReaderXML.hh"
#include "Vector3.hh"
#include "Transform.hh"
#include "Error.hh"
using namespace solid;

class Transform;


/** @brief The class TrapezoidalPrism.

    Convex with a trapezoidal prism shape. 

    @author A.WACHS - 2024 - Creation */
// ============================================================================
class TrapezoidalPrism : public Convex
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with edge length as input parameters
    @param width edge length in x
    @param a lower edge length in y 
    @param b upper edge length in y 
    @param h height in z */
    TrapezoidalPrism( double width, double a, double b, double h );

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    TrapezoidalPrism( istream& fileIn );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    TrapezoidalPrism( DOMNode* root );

    /** @brief Destructor */
    ~TrapezoidalPrism();
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the trapezoidal prism */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns the edge half-lengths of the trapezoidal prism */
    double const* getExtent() const;

    /** @brief Returns a vector of points describing the envelope of the 
    trapezoidal prism */
    vector<Point3> getEnvelope() const;

    /** @brief Returns point[0] in face number i
    @param i face number */
    Point3 getFirstPointFace( int i ) const;

    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the point indices */
    vector<vector<int> > const* getFaces() const;

    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape. Here returns the code 1111 */
    int getNbCorners() const;

    /** @brief Returns the trapezoidal prism volume */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Trapezoidal prism support function, returns the support point P,
    i.e. the point on the surface of the trapezoidal prism that satisfies 
    max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;

    /** @brief Returns the number of points to write the trapezoidal prism in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the 
    trapezoidal prism in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the trapezoidal prism in a
    Paraview format
    @param f output stream
    @param transform geometric transformation
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns a list of points describing the trapezoidal prism in a
    Paraview format
    @param transform geometric transformation
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the trapezoidal prism in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;

    /** @ brief Returns whether a point lies inside the trapezoidal prism
    @param pt point */
    bool isIn( Point3 const& pt ) const;

    /** @brief Writes the trapezoidal prism in an OBJ format
    @param f output stream
    @param transform geometric transformation 
    @param firstpoint_number number of the 1st point */
    void write_convex_OBJ( ostream& f, Transform const& transform,
    	size_t& firstpoint_number ) const;

    /** @ Returns the bounding volume to trapezoidal prism */
    BVolume* computeBVolume( unsigned int type ) const;
    
    /** @brief Performs advanced comparison of the trapezoidal prisms and 
    returns whether they match
    @param other the other trapezoidal prism */
    bool equalType_level2( Convex const* other ) const;    
    //@}


  protected:
    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference trapezoidal 
    prism, i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
    //@}


    /** @name Parameters*/
    //@{
    double* m_extent; /**< vector containing the half-length of the edges */
    vector<Point3> m_corners; /**< vector of the 8 corners/vertices of the
    	trapezoidal prism */
    double m_sinAngle; /**< sine of the angle of the inclined face and the
    	z-axis at the shortest y-edge */
    static vector< vector<int> > m_allFaces; /**< Face/vertices numbering */
    //@}


  private:
    /**@name Methods */
    //@{
    /** @brief Sets the corner/vertex coordinates and the face/vertices
    numbering */
    void setCornersFaces();
    //@}
};

#endif
