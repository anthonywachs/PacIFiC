#ifndef _POLYHEDRON_HH_
#define _POLYHEDRON_HH_

#include "Polytope.hh"
#include "IndexArray.hh"
#include "Basic.hh"
#include "ReaderXML.hh"



/** @brief The class Polyhedron.

    Convex with a polyhedral shape. From GJK Engine - A Fast and 
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Polyhedron : public Polytope 
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
    Polyhedron( istream& fileIn, int nb_point, VertexBase& ref, 
    	IndexArray& ia );

    /** @brief Copy constructor
    @param copy copied Polyhedron object  */
    Polyhedron( Polyhedron const& copy );
 
    /** @brief Destructor */
    virtual ~Polyhedron();
    //@}
  
  
   /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the polyhedron */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;
  
    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the vertex indices */
    vector<vector<int> > const* getFaces() const;

    /** @brief Returns the polyhedron volume */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Polyhedron support function, returns the support point P, i.e. 
    the point on the surface of the sphere that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;
  
    /** @brief Returns the number of elementary polytopes to write the
    polyhedron in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes the polyhedron in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;

    /** @brief Writes the polyhedron in a STL format
    @param f output stream
    @param transform geometric transformation */
    void write_convex_STL( ostream& f, Transform const& transform )
  	const;	
	
    /** @ brief Returns whether a point lies inside the polyhedron
    @param pt point */
    bool isIn( Point3 const& pt ) const;		
    //@}
  

    /** @name Static methods */
    //@{
    /** @brief Creates a polyhedron from an input stream
    @param fileIn inout stream */
    static Polyhedron* create( istream& fileIn );

    /** @brief Creates a polyhedron from an XML node
    @param root XML node */
    static Polyhedron* create( DOMNode* root );
    //@}


  private:
    /** @name Parameters */
    //@{
    IndexArray* m_cobound; /**< Indices of the vertices related to a given
    	vertex */
    mutable unsigned int m_curr_vertex; /**< index of the last vertex returned
    	by the support function */
    double* m_InertiaPoly; /**< polygon inertia tensor */
    double m_VolumePoly; /**< polyhedron volume */
    vector< vector<int> >* m_allFaces; /**< vertex indices in each face */ 
    //@}


    /**@name Constructor */
    //@{
    /** @brief Default constructor (forbidden) */
    Polyhedron();
    //@}
    

    /** @name Methods */
    //@{
    /** @brief Constructs the polyhedron knowing its description, i.e.,
    vertices, faces and indices
    @param nbface number of faces of the polyhedron
    @param face vertex indices in each face */
    void BuildPolyhedron( int nbface, IndexArray const* face );    
    
    /** @brief Reads the face indexing and builds the polygon by calling
    BuildPolygon 
    @param fileIn input stream */ 
    void readFaces( istream &fileIn );        

    /** @brief Returns the circumscribed radius of the reference polygon,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
  
    /** @brief Computes the contribution to inertia and volume of a tetrahedron
    defined by the center of mass (assuming that the center of mass is located 
    at (0,0,0)), the center of mass on a face and 2 consecutives vertices on 
    this face
    @param A2 center of mass of the face
    @param A3 a point of the face that H belongs to
    @param A4 the next point neighbor of A3 of the face that H belongs to */
    void computeVolumeInertiaContrib( Point3 const& A2, Point3 const& A3, 
	Point3 const& A4 );

    /** @brief Allocates the inertia tensor array and sets its component and the
    volume to 0 */
    void Initialisation();

    /** @brief Polyhedron support function when the vertex neighbor indexing is
    not known, returns the support point P, i.e. the point on the surface of 
    the polyhedron that satisfies max(P.v). Note: Never used so far, not even 
    clear when to use it, was used in orginal SOLID-2.0 software
    @param v direction vector */
    Point3 support20( Vector3 const& v ) const;
    //@}
};

#endif
