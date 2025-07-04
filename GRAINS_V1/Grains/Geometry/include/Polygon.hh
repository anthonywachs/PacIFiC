#ifndef _POLYGON_HH_
#define _POLYGON_HH_

#include "Polytope.hh"
#include "Basic.hh"
#include "IndexArray.hh"
#include "ReaderXML.hh"
#include <string>
using namespace std;


/** @brief The class Polygon.

    Convex with a polygonal shape. From GJK Engine - A Fast and 
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.
    
    @author Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Polygon : public Polytope 
{
  public:  
    /**@name Constructors */
    //@{
    /** @brief Constructor with an input stream, a number of vertices, a
    reference to the array of vertices and a reference to the array of indices
    as input parameters
    @param fileIn input stream
    @param nb_point number of vertices
    @param ref reference to the array of vertices
    @param ia reference to the array of indices */
    Polygon( istream& fileIn, int nb_point, VertexBase& ref, IndexArray& ia );

    /** @brief Copy constructor
    @param copy copied Polytope object  */
    Polygon( Polygon const& copy );
 
    /** @brief Destructor */
    virtual ~Polygon();
    //@}
  

    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the polygon */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns the polygon surface area */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Polygon support function, returns the support point P, i.e. the
    point on the surface of the polygon that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;
  
    /** @brief Returns the number of elementary polytopes to write the polygon 
    in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes the polygon in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;
	
    /** @ brief Returns whether a point lies inside the polygon
    @param pt point */
    bool isIn( Point3 const& pt ) const;
    
    /** @brief Performs advanced comparison of the two polygons and returns
    whether they match
    @param other the other convex */
    bool equalType_level2( Convex const* other ) const;    
    //@}


    /** @name Methods static */
    //@{
    /** @brief Creates a polygon from an input stream
    @param fileIn inout stream */
    static Polygon* create( istream& fileIn );

    /** @brief Creates a polygon from an XML node
    @param root XML node */
    static Polygon* create( DOMNode* root );
    //@}
  

  private:
    /** @name Parameters */
    //@{
    IndexArray* m_cobound; /**< Indices of the vertices related to a given
    	vertex */
    mutable unsigned int m_curr_vertex; /**< index of the last vertex returned
    	by the support function */
    double* m_InertiaPoly; /**< polygon inertia tensor */
    double m_surfaceArea; /**< polygon surface area */
    //@}


    /**@name Constructor */
    //@{
    /** @brief Default constructor (forbidden) */
    Polygon();
    //@}


    /**@name Methods */
    //@{
    /** @brief Constructs the polygon knowing its description, i.e., vertices
    and indices
    @param nbedge number of edges of the polygon
    @param edge indices of the vertices in the edge, actually edge has a single
    element with all indices as the polygon perimeter is considered as a single
    edge */
    void BuildPolygon( int nbedge, IndexArray const* edge );
  
    /** @brief Computes the contribution to inertia and surface area of an edge
    defined by 2 consecutive vertices relative to the center of mass position
    @param P 1st vertex
    @param Q 2nd vertex */
    void computeSurfaceInertiaContrib( Point3 const& P, Point3 const& Q );
    
    /** @brief Reads the face indexing and builds the polygon by calling
    BuildPolygon 
    @param fileIn input stream */ 
    void readFaces( istream &fileIn );

    /** @brief Returns the circumscribed radius of the reference polygon,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
  
    /** @brief Allocates the inertia tensor array and sets its component and the
    surface area to 0 */
    void Initialisation();
  //@}   
};

#endif
