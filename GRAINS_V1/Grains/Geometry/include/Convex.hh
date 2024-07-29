#ifndef _CONVEX_HH_
#define _CONVEX_HH_

#include "Shape.hh"
#include "BVolume.hh"
#include "Transform.hh"
#include "Point3.hh"
#include "ReaderXML.hh"
#include "WriterXML.hh"
#include <vector>
#include <list>
using namespace std;
using namespace solid;


class Convex;
ostream& operator << ( ostream& fileOut, Convex const& convex);
istream& operator >> ( istream& fileIn, Convex& convex);


// List of types of convex
enum ConvexType {
  SPHERE,
  DISC2D,
  POLYHEDRON,
  POLYGON,
  BOX,
  CONE,
  CYLINDER,
  POINT,
  SEGMENT,
  SUPERQUADRIC,
  RECTANGLE2D
};


/** @brief The class Convex.

    High level class describing convex shapes. From GJK Engine - A Fast and
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.

    @author D.PETIT - Institut Francais du Petrole - 2000 - Modification
    @author D.RAKOTONIRINA - IFP Energies Nouvelles - 2014 - Modification
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Convex : public Shape
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Destructor */
    virtual ~Convex();
    //@}


    /**@name Methods */
    //@{
    /** @brief Returns shape type */
    ShapeType getType() const;
    
    /** @brief Returns whether two convexes are of the same type 
    @param other the other convex
    @param level2 if true, performs advanced comparison */
    bool equalType( Convex const* other, bool const& level2 ) const;    
    //@}    


    /**@name Virtual methods */
    //@{
    /** @brief Returns the convex shape bounding box
    @param t geometric transformation */
    virtual BBox bbox( Transform const& t ) const;

    /** @brief Returns the convex shape bounding volume
    @param type 1 = OBB, 2 = OBC */
    virtual BVolume* computeBVolume( unsigned int type ) const;

    /** @brief Convex support function, returns the support point P, i.e. the
    point on the surface of the convex shape that satisfies max(P.v)
    @param v direction vector */
    virtual Point3 support( Vector3 const& v ) const = 0;

    /** @brief Computes and returns the circumscribed radius of the reference
    convex shape, i.e., without applying any transformation.
    Note: the circumscribed radius of the convex is unchanged
    by a translation (TRANSLATION) or a rotation (ROTATION) or a combination of
    these two, but is changed by a transformation of SCALING type. However, the
    use of scaling/stretching is not allowed in the definition of components in
    the input file */
    virtual double computeCircumscribedRadius() const = 0;

    /** @brief Returns the convex type */
    virtual ConvexType getConvexType() const = 0;

    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    virtual bool BuildInertia( double* inertia, double* inertia_1 ) const = 0;

    /** @brief Returns a clone copy of the convex */
    virtual Convex* clone() const = 0;

    /** @brief Returns a vector of points describing the envelope of the convex
    */
    virtual vector<Point3> getEnvelope() const;

    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the vertex indices */
    virtual vector<vector<int> > const* getFaces() const;

    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape having to vertices/corners as, e.g., a cylinder */
    virtual int getNbCorners() const;

    /** @brief Returns the volume of the convex shape */
    virtual double getVolume() const = 0;

    /** @brief Returns the number of points to write the convex shape in a
    Paraview format */
    virtual int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the convex
    shape in a Paraview format */
    virtual int numberOfCells_PARAVIEW() const;

    /** @brief Writes the points describing the convex shape in a
    Paraview format
    @param f output stream
    @param transform geometric transformation
    @param translation additional center of mass translation */
    virtual void write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the convex shape in a STL format
    @param f output stream
    @param transform geometric transformation */
    virtual void write_convex_STL( ostream& f, Transform const& transform )
  	const;

    /** @brief Returns a list of points describing the convex shape in a
    Paraview format
    @param transform geometric transformation
    @param translation additional center of mass translation */
    virtual list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL )
  	const;

    /** @brief Writes the convex shape in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    virtual void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;

    /** @brief Returns whether the convex shape is a sphere */
    virtual bool isSphere() const;

    /** @brief Returns an orientation vector describing the convex shape angular
    position
    @param transform geometric transformation */
    virtual Vector3 computeOrientationVector( Transform const* transform )
    	const;

    /** @ brief Returns whether a point lies inside the convex shape
    @param pt point */
    virtual bool isIn( Point3 const& pt ) const = 0;
    
    /** @brief Performs advanced comparison of the two convexes and returns
    whether they match
    @param other the other convex */
    virtual bool equalType_level2( Convex const* other ) const = 0;       
    //@}


    /** @name Friend methods */
    //@{
    /** @brief Output operator
    @param fileOut output stream
    @param convex Convex object */
    friend ostream& operator << ( ostream& fileOut, Convex const& convex );

    /** @brief Input operator
    @param fileIn input stream
    @param convex Convex object*/
    friend istream& operator >> ( istream& fileIn, Convex& convex );
    //@}


  protected:
    /**@name Contructors */
    //@{
    /** @brief Default constructor (forbidden except in derived classes) */
    Convex();
    //@}


    /**@name Virtual methods */
    //@{
    /** @brief Output operator (is called by <<)
    @param fileOut output stream */
    virtual void writeShape( ostream &fileOut ) const = 0;

    /** @brief Input operator (is called by >>)
    @param fileIn input stream */
    virtual void readShape( istream &fileIn ) = 0;
    //@}
};



/**@name Convex : External methods from the SOLID Software */
//@{
/** @brief Returns whether 2 convex shapes intersect
@param a convex shape A
@param b convex shape B
@param a2w geometric tramsformation describing convex A in the world reference
frame
@param b2w geometric tramsformation describing convex B in the world reference
frame
@param v initial direction of GJK */
bool intersect( Convex const& a, Convex const& b, Transform const& a2w,
	Transform const& b2w, Vector3& v );

/** @brief Returns whether 2 convex shapes intersect
@param a convex shape A
@param b convex shape B
@param b2a geometric tramsformation describing convex B in the reference frame
of A
@param v initial direction of GJK */
bool intersect( Convex const& a, Convex const& b, Transform const& b2a,
	Vector3& v );

/** @brief Returns whether 2 convex shapes intersect and if they intersect
returns an intersection point per convex shape in each convex reference frame
@param a convex shape A
@param b convex shape B
@param a2w geometric tramsformation describing convex A in the world reference
frame
@param b2w geometric tramsformation describing convex B in the world reference
frame
@param v initial direction of GJK
@param pa intersection point of A in the reference frame of A
@param pb intersection point of B in the reference frame of B */
bool common_point( Convex const& a, Convex const& b, Transform const& a2w,
	Transform const& b2w, Vector3& v, Point3& pa, Point3& pb );

/** @brief Returns whether 2 convex shapes intersect and if they intersect
returns an intersection point per convex shape in each convex reference frame
@param a convex shape A
@param b convex shape B
@param b2a geometric tramsformation describing convex B in the reference frame
of A
@param v initial direction of GJK
@param pa intersection point of A in the reference frame of A
@param pb intersection point of B in the reference frame of B */
bool common_point( Convex const& a, Convex const& b, Transform const& b2a,
	Vector3& v, Point3& pa, Point3& pb );

/** @brief Returns the minimal distance between 2 convex shapes and a point per
convex shape that represents the tips of the minimal distance segment
@param a convex shape A
@param b convex shape B
@param a2w geometric tramsformation describing convex A in the world reference
frame
@param b2w geometric tramsformation describing convex B in the world reference
frame
@param pa point representing one tip of the minimal distance segment on A
@param pb point representing the other tip of the minimal distance segment on
B
@param nbIter number of iterations of GJK for convergence */
double closest_points( Convex const& a, Convex const& b, Transform const& a2w,
	Transform const& b2w, Point3& pa, Point3& pb, int& nbIter );
//@}

#endif
