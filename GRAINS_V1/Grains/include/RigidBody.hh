#ifndef _RIGIDBODY_HH_
#define _RIGIDBODY_HH_

#include "Error.hh"
#include "Basic.hh"
#include "Transform.hh"
#include "Convex.hh"
#include "Vector3.hh"
#include "Point3.hh"
#include "BVolume.hh"
using namespace solid;

class Cell;
class Torsor;


/** @brief The class RigidBody.

    A combination of a convex shape and a transformation.

    @author F.PRADEL - Institut Francais du Petrole - 1999 - Creation
    @author G.FERRER - Institut Francais du Petrole - 2000 - Polymorphisme
    @author D.PETIT  - Institut Francais du Petrole - 2000 - Modification
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class RigidBody
{
  public:
    /**@name Constructeurs */
    //@{
    /** @brief Default constructor */
    RigidBody();

    /** @brief Constructor with a convex and a transformation as input
    parameters, used exclusively by compObstacle whose own shape is not defined
    @param convex_ convex
    @param position_ transformation */
    RigidBody( Convex* convex_, Transform const& position_ );

    /** @brief Copy constructor
    @param form copied RigidBody object */
    RigidBody( RigidBody const& form );

    /** @brief Destructor */
    virtual ~RigidBody() ;
    //@}


    /**@name Methods */
    //@{
    /** @brief Returns the bounding box of the rigid body in its current
    configuration */
    BBox BoxRigidBody() const;

    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns the distance to another rigid body
    @param neighbor the other rigid body */
    double DistanceTo( RigidBody const& neighbor ) const;

    /** @brief Returns whether there is geometric contact with another rigid
    body
    @param neighbor the other rigid body */
    bool isContact( RigidBody const& neighbor ) const;

    /** @brief Returns whether the rigid body is close to another rigid body in
    the sense of whether their respective bounding boxes overlap
    @param neighbor the other rigid body */
    bool isClose( RigidBody const& neighbor ) const;

    /** @brief Applies a rotation defined by a quaternion to the rigid body
    @param q quaternion representing the rotation */
    void Rotate( Quaternion const& q );

    /** @brief Applies a transformation trot to the right, i.e., this = this o
    trot, which means first trot then this
    @param trot transformation */
    void composeRightByTransform( Transform const& trot );

    /** @brief Applies a transformation trot to the left, i.e., this = trot o
    this, which means first this then trot
    @param trot transformation */
    void composeLeftByTransform( Transform const& trot );

    /** @brief Applies a rotation defined by a transformation trot to the left,
    i.e., this = trot o this, which means first this then trot. This composition
    leaves the origin unchanged but does not check that trot is indeed a
    rotation
    @param trot transformation � appliquer */
    void composeLeftByRotation( Transform const& trot );

    /** @brief Applies a translation to the left, i.e., this = translation
    o this, which means first this then translation.
    @param v translation vector */
    void composeLeftByTranslation( Vector3 const& v );

    /** @brief Copies the rigid body's transformation in a 1D array
    @param vit 1D array where transformation is copied
    @param i start index to copy in the 1D array */
    void copyTransform( double* vit, int i ) const;

    /** @brief Copies the rigid body's transformation in a 1D array composed on
    the left by a translation (useful for periodic particles in parallel)
    @param vit 1D array where transformation is copied
    @param i start index to copy in the 1D array
    @param vec translation vecteur */
    void copyTransform( double* vit, int i, Vector3 const& vec ) const;

    /** @ brief Returns whether a point lies inside the rigid body
    @param pt point */
    bool isIn( Point3 const& pt ) const;
    //@}


    /**@name Set methods */
    //@{
    /** @brief Sets the origin of the rigid body's transformation
    @param pos origin coordinates as a 3-element 1D array  */
    void setOrigin( double const* pos );

    /** @brief Sets the origin of the rigid body's transformation
    @param gx x-coordinate of the origin
    @param gy y-coordinate of the origin
    @param gz z-coordinate of the origin */
    void setOrigin( double gx, double gy, double gz );

    /** @brief Sets the origin of the rigid body's transformation
    @param pos origin coordinates as a Point3 */
    void setOrigin( Point3 const& pos );

    /** @brief Sets the rigid body's transformation with an 1D array of 12
    values (see class Transform for details)
    @param pos 1D array of values containing the tranformation coefficients */
    void setTransform( double const* pos );

    /** @brief Sets the rigid body's circumscribed radius
    @param r circumscribed radius */
    void setCircumscribedRadius( double r );

    /** @brief Sets the rigid body's transformation with a transformation
    @param transform_ transformation */
    void setTransform( Transform const& transform_ );
    //@}


    /**@name Get methods */
    //@{
    /** @brief Returns a pointer to the rigid body' center of mass position */
    Point3 const* getCentre() const;

    /** @brief Copies the rigid body' center of mass position in a 3-element 1D
    array
    @param pos 3-element 1D array where the rigid body' center of mass
    cooridnates are copied */
    void getCentre( double *pos ) const;

    /** @brief Returns a pointer to the convex shape */
    Convex const* getConvex() const;

    /** @brief Returns a pointer to the rigid body's transformation */
    Transform const* getTransform() const;

    /** @brief Returns a pointer to the rigid body's transformation */
    Transform* getTransform();

    /** @brief Returns the rigid body's circumscribed radius */
    double getCircumscribedRadius() const;

    /** @brief Returns the rigid body volume */
    double getVolume() const;

    /** @brief Returns the rigid body bounding cylinder */
    BVolume const& getBVolume() const;
    //@}


    /**@name I/O methods */
    //@{
    /** @brief Reads the rigid body's transformation from an input stream
    @param fileIn input stream */
    void readPosition( istream& fileIn );

    /** @brief Reads the rigid body's transformation from an input stream with
    the 2014 reload format
    @param fileIn input stream */
    void readPosition2014( istream& fileIn );

    /** @brief Reads the rigid body's transformation in binary format from an
    input stream with the 2014 reload format
    @param fileIn input stream */
    void readPosition2014_binary( istream& fileIn );

    /** @brief Writes the rigid body's transformation in an output stream
    @param fileOut output stream */
    void writePosition( ostream& fileOut ) const;

    /** @brief Writes the geometric features of the rigid body in its current
    position in an output stream in a format suitable to the coupling with a
    fluid solver. Note: this method works for discs, polygons, polyhedrons,
    spheres and 3D cylinders
    @param fluid output stream */
    void writePositionInFluid( ostream& fluid );

    /** @brief Writes the rigid body's "static" data, i.e., the convex geometric
    description only (without any transformation)
    @param fileOut output stream */
    virtual void writeStatic( ostream& fileOut ) const;

    /** @brief Writes a list of points describing the rigid body's convex shape
    in a Paraview format
    @param f output stream
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f,
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the rigid body's convex shape in a STL format
    @param f output stream */
    void write_convex_STL( ostream& f ) const;

    /** @brief Returns a list of points describing the rigid body's convex shape
    in a Paraview format
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Vector3 const* translation = NULL )
    	const;
    //@}


    /**@name Friend methods */
    //@{
    /** @brief Returns whether 2 rigid bodies intersect
    @param a rigid body A
    @param b rigid body B */
    friend bool intersect( RigidBody const& a, RigidBody const& b );
    //@}


  protected:
    /**@name Parameters */
    //@{
    Transform m_transform; /**< rigid body's transformation: center of mass
  	position and body orientation */
    Convex *m_convex; /**< convex shape */
    BVolume* m_boundingVolume; /** The bounding volume of the convex body **/
    double m_circumscribedRadius; /**< circumscribed radius */
    double m_volume; /**< rigid body volume */
    //@}
};


/** @brief Intersection entre les deux RigidBody indiquees ? */
bool intersect( const RigidBody &a, const RigidBody &b );

#endif
