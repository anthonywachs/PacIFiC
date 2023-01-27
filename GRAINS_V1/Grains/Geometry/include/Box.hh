#ifndef _BOX_HH_
#define _BOX_HH_

#include "Convex.hh"
#include "ReaderXML.hh"
#include "Vector3.hh"
#include "Transform.hh"
#include "Error.hh"
using namespace solid;

class Transform;
/** @brief The class Box.

    Convex with a box shape. From GJK Engine - A Fast and
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.

    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Box : public Convex
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with edge length as input parameters
    @param x edge length in x
    @param y edge length in y
    @param z edge length in z */
    Box( double x = 0., double y = 0., double z = 0. );

    /**@brief Constructor with a vector containing the edge half-lengths as
    input parameters
    @param extent_ vector containing the edge half-lengths */
    Box( Vector3 const& extent_ );

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    Box( istream& fileIn );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    Box( DOMNode* root );

    /** @brief Destructor */
    ~Box();
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

    /** @brief Returns the edge half-lengths of the box */
    Vector3 const& getExtent() const;

    /** @brief Returns a vector of points describing the envelope of the box */
    vector<Point3> getEnvelope() const;

    /** @brief Returns point[0] in face number i
    @param i face number */
    Point3 getFirstPointFace( int i ) const;

    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the point indices */
    vector<vector<int> > const* getFaces() const;

    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape. Here returns the code 444 for a 2D box and 666 for
    a 3D box */
    int getNbCorners() const;

    /** @brief Returns the box volume */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Box support function, returns the support point P, i.e. the
    point on the surface of the box that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;

    /** @brief Returns the number of points to write the box in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the box
    in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the box in a
    Paraview format
    @param f output stream
    @param transform geometric transformation
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns a list of points describing the box in a
    Paraview format
    @param transform geometric transformation
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the box in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;

    /** @brief Returns the contact point in the box reference frame and the
    overlapping distance between the box and a sphere. If no contact, returned
    contact point is world reference frame origin (0,0,0)
    @param SphereCenter sphere center of mass coordinates in the box reference
    frame
    @param SphereRadius sphere radius
    @param overlap distance between the box and the sphere (negative if contact)
    @param warningSphereCenterInBox if true, throw a contact error, otherwise
    set the overlap to a negative value of -10000 by default to tell that
    there is contact and returns the origin as contact point */
    Point3 IntersectionPointSPHERE( Point3 const& SphereCenter,
  	double const& SphereRadius, double& overlap,
<<<<<<< HEAD
	bool warningSphereCenterInBox = true ) const
    	throw(ContactError);
=======
	  bool warningSphereCenterInBox = true ) const;    
>>>>>>> NewGrains

    /** @brief Same as IntersectionPointSPHERE except that it returns a non zero
    normal distance only for a configuration sphere-face (i.e. returns zero with
    the minimal distance is wrt an edge or a corner)
    @param SphereCenter sphere center of mass coordinates in the box reference
    frame
    @param SphereRadius sphere radius
    @param gap distance between the box and the sphere (negative if contact)
    */
    Point3 ProjectedPointSPHERE( Point3 const& SphereCenter,
  	double const& SphereRadius,
	double& gap ) const;

    /** @ brief Returns whether a point lies inside the box
    @param pt point */
    bool isIn( Point3 const& pt ) const;

    /** @ Returns the bounding cylinder to box */
    BCylinder getBCylinder() const;
    //@}


  protected:
    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference box,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
    //@}


    /** @name Parameters*/
    //@{
    Vector3 m_extent; /**< vector containing the half-legnth of the edges */
    vector<Point3> m_corners; /**< vector of the 8 corners/vertices of the
    	box */
    vector<Point3>* m_corners2D_XY; /**< vector of the 4 corners/vertices in the
    	XY plan */
    static vector< vector<int> > m_allFaces; /**< Face/vertices numbering */
    //@}


  private:
    /**@name Methods */
    //@{
    /** @brief Sets the corner/vertex coordinates and the face/vertices
    numbering */
    void setCornersFaces();

    /** @brief Returns the contact point in the box reference frame and the
    overlapping distance between a sphere and a box corner. If no contact,
    returned contact point is world reference frame origin (0,0,0)
    @param SphereCenter sphere center of mass coordinates in the box reference
    frame
    @param SphereRadius sphere radius
    @param cornerNumber corner number
    @param overlap distance between the box and the sphere (negative if contact)
    */
    Point3 ContactCornerSPHERE( Point3 const& SphereCenter,
  	double const& SphereRadius,
	int cornerNumber,
	double& overlap ) const;

    /** @brief Returns the contact point in the box reference frame and the
    overlapping distance between a sphere and a box edge. If no contact,
    returned contact point is world reference frame origin (0,0,0)
    @param SphereCenter sphere center of mass coordinates in the box reference
    frame
    @param SphereRadius sphere radius
    @param cornerNumber corner number of one of the 2 tips of the edge
    @param projectionDirection plane in which the sphere center is projected
    @param overlap distance between the box and the sphere (negative if contact)
    */
    Point3 ContactEdgeSPHERE( Point3 const& SphereCenter,
  	double const& SphereRadius,
	int cornerNumber,
	int projectionDirection,
	double& overlap ) const;
    //@}
};

#endif
