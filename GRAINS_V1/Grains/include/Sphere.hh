#ifndef _SPHERE_HH_
#define _SPHERE_HH_

#include "Convex.hh"
#include "Point3.hh"
#include "ReaderXML.hh"
using namespace solid;


/** @brief The class Sphere.

    Convex with a spherical shape. From GJK Engine - A Fast and
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.

    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Sphere : public Convex
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Constructor with radius as input parameter
    @param r sphere radius */
    Sphere( double r = 0 );

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    Sphere( istream& fileIn );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    Sphere( DOMNode* root );

    /** @brief Destructor */
    ~Sphere();
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the sphere */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns a vector of points describing the surface of the sphere.
    Here simply returns the point (0,0,0) as a convention */
    vector<Point3> getEnvelope() const;

    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the point indices. Returns a null pointer as a
    convention */
    vector<vector<int> > const* getFaces() const;

    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape. Here 1 as a convention */
    int getNbCorners() const;

    /** @brief Returns the sphere volume */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Sphere support function, returns the support point P, i.e. the
    point on the surface of the sphere that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;

    /** @brief Returns the number of points to write the sphere in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the sphere
    in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the sphere in a
    Paraview format
    @param f output stream
    @param transform geometric transformation
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns a list of points describing the sphere in a
    Paraview format
    @param transform geometric transformation
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the sphere in a STL format
    @param f output stream
    @param transform geometric transformation */
    void write_convex_STL( ostream& f, Transform const& transform ) const;

    /** @brief Writes the sphere in an OBJ format
    @param f output stream
    @param transform geometric transformation 
    @param firstpoint_number number of the 1st point */
    void write_convex_OBJ( ostream& f, Transform const& transform,
    	size_t& firstpoint_number ) const;	

    /** @brief Writes the sphere in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;

    /** @brief Sets the number of points per quarter of the equator line for
    Paraview post-processing, i.e., controls the number of facets in the sphere
    reconstruction in Paraview
    @param nbpts number of point per quarter of the equator line */
    static void SetvisuNodeNbPerQar( int nbpts );

    /** @ brief Returns whether a point lies inside the sphere
    @param pt point */
    bool isIn( Point3 const& pt ) const;
    
    /** @brief Performs advanced comparison of the two spheres and returns
    whether they match
    @param other the other sphere */
    bool equalType_level2( Convex const* other ) const;
    
    /** @brief Returns the sphere bounding volume
    @param type 1 = OBB, 2 = OBC */
    BVolume* computeBVolume( unsigned int type ) const;        
    //@}


  protected:
    /**@name Parameters */
    //@{
    double m_radius; /**< sphere radius */
    static int m_visuNodeNbPerQar; /**< number of points per quarter of the
    	equator line for Paraview post-processing */
    //@}


    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference sphere,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
    //@}


  private:
    /**@name Methods */
    //@{
    /** @brief Writes a triangular facet in the STL format based on its 3
    vertices and the center of mass coordinates
    @param f output stream
    @param GC center of mass coordinates
    @param pp1 point 1
    @param pp2 point 2
    @param pp3 point 3 */
    void write_STLfacet_sphere( ostream& f, Point3 const& GC,
  	Point3 const& pp1,
  	Point3 const& pp2,
  	Point3 const& pp3 ) const;
    //@}
};

#endif
