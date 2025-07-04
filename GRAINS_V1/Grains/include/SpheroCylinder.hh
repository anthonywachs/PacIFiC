#ifndef _SPHEROCYLINDER_HH_
#define _SPHEROCYLINDER_HH_

#include "Convex.hh"
#include "ReaderXML.hh"
#include "Vector3.hh"
#include "Transform.hh"
#include "Error.hh"
using namespace solid;

class Transform;


/** @brief The class SpheroCylinder.

    Convex with a SpheroCylinder shape. 

    @author A.WACHS - 2024 - Creation */
// ============================================================================
class SpheroCylinder : public Convex
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with radius and height as input parameters
    @param r radius of the elementary cylinder and the two spherical caps
    @param h height of the elementary cylinder */
    SpheroCylinder( double r, double h );

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    SpheroCylinder( istream& fileIn );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    SpheroCylinder( DOMNode* root );

    /** @brief Destructor */
    ~SpheroCylinder();
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the SpheroCylinder */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns a vector of points describing the surface of the 
    SpheroCylinder */
    vector<Point3> getEnvelope() const;

    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the point indices. Returns a null pointer as a
    convention */
    vector<vector<int> > const* getFaces() const;

    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape. Here returns the code 3333 */
    int getNbCorners() const;

    /** @brief Returns the SpheroCylinder volume */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief SpheroCylinder support function, returns the support point P,
    i.e. the point on the surface of the SpheroCylinder that satisfies 
    max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;

    /** @brief Returns the number of points to write the SpheroCylinder in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the 
    SpheroCylinder in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the SpheroCylinder in a
    Paraview format
    @param f output stream
    @param transform geometric transformation
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns a list of points describing the SpheroCylinder in a
    Paraview format
    @param transform geometric transformation
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the SpheroCylinder in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;

    /** @ brief Returns whether a point lies inside the SpheroCylinder
    @param pt point */
    bool isIn( Point3 const& pt ) const;
    
    /** @brief Writes the SpheroCylinder in an OBJ format
    @param f output stream
    @param transform geometric transformation 
    @param firstpoint_number number of the 1st point */
    void write_convex_OBJ( ostream& f, Transform const& transform,
    	size_t& firstpoint_number ) const;    

    /** @ Returns the bounding volume to SpheroCylinder */
    BVolume* computeBVolume( unsigned int type ) const;
    
    /** @brief Performs advanced comparison of the SpheroCylinders and 
    returns whether they match
    @param other the other SpheroCylinder */
    bool equalType_level2( Convex const* other ) const;
    
    /** @brief Sets the number of points per quarter of the elementary 
    spherocylinder perimeter for Paraview post-processing, i.e., controls the 
    number of facets in the spherocylinder reconstruction in Paraview
    @param nbpts number of point per quarter of the elementary cylinder 
    perimeter */
    static void SetvisuNodeNbPerQar( int nbpts );        
    //@}


  protected:
    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference spherocylinder,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
    //@}


    /** @name Parameters */
    //@{
    double m_radius; /**< Radius of the elementary cylinder and the two
    	spherical caps */
    double m_height; /**< Height of the elementary cylinder */
    static int m_visuNodeNbPerQar; /**< number of points over a quarter of 
    	the circular perimeter of the spherocylinder for Paraview 
	post-processing */	
    //@}
};

#endif
