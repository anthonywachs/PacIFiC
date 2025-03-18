#ifndef _SPHEROCYL_HH_
#define _SPHEROCYL_HH_

#include "Convex.hh"
#include "ReaderXML.hh"
#include "Vector3.hh"
#include "Transform.hh"
#include "Error.hh"
using namespace solid;

class Transform;


/** @brief The class SpheroCyl.

    Convex with a spherocylinder shape. 

    @author A.WACHS - 2024 - Creation */
// ============================================================================
class SpheroCyl : public Convex
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with radius and height as input parameters
    @param r radius of the elementary cylinder and the two spherical caps
    @param h height of the elementary cylinder */
    SpheroCyl( double r, double h );

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    SpheroCyl( istream& fileIn );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    SpheroCyl( DOMNode* root );

    /** @brief Destructor */
    ~SpheroCyl();
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the spherocylinder */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns a vector of points describing the envelope of the 
    spherocylinder */
    vector<Point3> getEnvelope() const;

    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the point indices. Returns a null pointer as a
    convention */
    vector<vector<int> > const* getFaces() const;

    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape. Here returns the code 3333 */
    int getNbCorners() const;

    /** @brief Returns the spherocylinder volume */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief spherocylinder support function, returns the support point P,
    i.e. the point on the surface of the spherocylinder that satisfies 
    max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;

    /** @brief Returns the number of points to write the spherocylinder in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the 
    spherocylinder in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the spherocylinder in a
    Paraview format
    @param f output stream
    @param transform geometric transformation
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns a list of points describing the spherocylinder in a
    Paraview format
    @param transform geometric transformation
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the spherocylinder in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;

    /** @ brief Returns whether a point lies inside the spherocylinder
    @param pt point */
    bool isIn( Point3 const& pt ) const;

    /** @ Returns the bounding volume to spherocylinder */
    BVolume* computeBVolume( unsigned int type ) const;
    
    /** @brief Performs advanced comparison of the spherocylinders and 
    returns whether they match
    @param other the other spherocylinder */
    bool equalType_level2( Convex const* other ) const;
    
    /** @brief Sets the number of points per quarter of the elementary cylinder 
    perimeter for Paraview post-processing, i.e., controls the number of facets
    in the sphero-cylinder reconstruction in Paraview
    @param nbpts number of point per quarter of the elementary cylinder 
    perimeter */
    static void SetvisuNodeNbPerQar( int nbpts );        
    //@}


  protected:
    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference trapezoidal 
    prism, i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
    //@}


    /** @name Parameters */
    //@{
    double m_radius; /**< Radius of the elementary cylinder and the two
    	spherical caps */
    double m_height; /**< Height of the elementary cylinder */
    static int m_visuNodeNbPerQar; /**< number of points over a quarter of 
    	the circular perimeter of the sphero-cylinder for Paraview 
	post-processing */	
    //@}
};

#endif
