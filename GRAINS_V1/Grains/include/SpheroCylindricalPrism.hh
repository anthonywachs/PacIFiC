#ifndef _SPHEROCYLINDRICALPRISM_HH_
#define _SPHEROCYLINDRICALPRISM_HH_

#include "Convex.hh"
#include "ReaderXML.hh"
#include "Vector3.hh"
#include "Transform.hh"
#include "Error.hh"
using namespace solid;

class Transform;


/** @brief The class SpheroCylindricalPrism.

    Convex with a spherocylindrical prism shape. 

    @author A.WACHS - 2024 - Creation */
// ============================================================================
class SpheroCylindricalPrism : public Convex
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with edge length as input parameters
    @param width edge length in x
    @param r eadius of the two half cylinders
    @param l length 
    @param h height */
    SpheroCylindricalPrism( double r, double l, double h );

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    SpheroCylindricalPrism( istream& fileIn );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    SpheroCylindricalPrism( DOMNode* root );

    /** @brief Destructor */
    ~SpheroCylindricalPrism();
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the spherocylindrical prism */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns a vector of points describing the envelope of the 
    spherocylindrical prism */
    vector<Point3> getEnvelope() const;

    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the point indices. Returns a null pointer as a
    convention */
    vector<vector<int> > const* getFaces() const;

    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape. Here returns the code 4444 */
    int getNbCorners() const;

    /** @brief Returns the spherocylindrical prism volume */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Spherocylindrical prism support function, returns the support 
    point P, i.e. the point on the surface of the spherocylindrical prism that 
    satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;

    /** @brief Returns the number of points to write the spherocylindrical 
    prism in a Paraview format */
    int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the 
    spherocylindrical prism in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the spherocylindrical prism 
    in a Paraview format
    @param f output stream
    @param transform geometric transformation
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns a list of points describing the spherocylindrical prism 
    in a Paraview format
    @param transform geometric transformation
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the spherocylindrical prism in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;

    /** @ brief Returns whether a point lies inside the spherocylindrical prism
    @param pt point */
    bool isIn( Point3 const& pt ) const;

    /** @ Returns the bounding volume to spherocylindrical prism */
    BVolume* computeBVolume( unsigned int type ) const;
    
    /** @brief Performs advanced comparison of the spherocylindrical prisms and 
    returns whether they match
    @param other the other spherocylindrical prism */
    bool equalType_level2( Convex const* other ) const;    
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
    double m_radius; /**< Radius of the two half cylinders */
    double m_length; /**< Length */ 
    double m_height; /**< Height */       
    static int m_visuNodeNbPerHalf; /**< number of points over 
    	the circular perimeter of each half cylinder (half of a cricle ) 
	for Paraview post-processing */	
    //@}
};

#endif
