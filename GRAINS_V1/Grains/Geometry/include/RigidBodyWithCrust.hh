#ifndef _FORMEVDW_HH_
#define _FORMEVDW_HH_

#include "RigidBody.hh"
#include "ReaderXML.hh"
#include "WriterXML.hh"
#include "Error.hh"

class PointContact;


/** @brief The class RigidBodyWithCrust.

    A combination of a convex shape, a transformation and a surface crust.

    @author G. FERRER - Institut Francais du Petrole - 2002 - Creation
    @author D. RAKOTONIRINA - IFP Energies Nouvelles - Sept. 2014
    - Modification
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class RigidBodyWithCrust : public RigidBody
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    RigidBodyWithCrust();

    /** @brief Copy constructor
    @param rbwc copied RigidBody object */
    RigidBodyWithCrust( RigidBodyWithCrust const& rbwc );

    /** @brief Constructor with a convex and a transformation as input
    parameters, used exclusively by compObstacle whose own shape is not defined
    @param convex_ convex
    @param position_ transformation */
    RigidBodyWithCrust( Convex* convex_, Transform const& position_ );

    /** @brief Constructor from an input stream and a convex type
    @param fileIn input stream
    @param type convex type, used in case of elementary particle only */
    RigidBodyWithCrust( istream &fileIn, string type = ""  );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    RigidBodyWithCrust( DOMNode* root );

    /** Destructor */
    virtual ~RigidBodyWithCrust() ;
    //@}


    /**@name Methods */
    //@{
    /** @brief Returns the bounding box extended by the crust thickness */
    BBox BoxRigidBody() const;

    /** @brief Returns the features of the contact: contact point location,
    overlap vector (vector joining the points on each rigid body surface that
    realize the minimal distance between the shrunk rigid bodies, divided by the
    minimal distance between the shrunk rigid bodies and multiplied by the
    sum of the crust thicknesses minus the minimal distance between the shrunk
    rigid bodies, i.e., minus the overlap), overlap distance = minimal distance
    between the shrunk rigid bodies minus the sum of the crust thicknesses.
    Note: contact exists if overlap distance is negative, i.e., minimal distance
    between the shrunk rigid bodies < sum of the crust thicknesses
    @param neighbor the other rigid body
    @exception if the minimal distance between the shrunk rigid bodies is 0, in
    practice less than EPSILON defined in Basic.H, which means that the shrunk
    rigid bodies already touch or overlap */
    PointContact ClosestPoint( RigidBodyWithCrust &neighbor )
      throw(ContactError);

    /** @brief Returns the features of the contact when the overlap computed by
    ClosestPoint is too large, the method artificially increases the size of the
    crust thickness for this particular contact detection by a factor > 1 and
    imposes the overlap distance to minus the sum of the crust thicknesses. This
    is useful when a few contacts involve a slightly large overlap due for
    instance (i) to the user's contact parameters not being properly set, (ii)
    an unexpectedly high colliding velocity between 2 rigid bodies that
    constitutes a rare by physically meaningful event or (iii) when a contact
    has not been detected by GJK for a few time steps preceding the call to
    this method that constitutes an extremely rare event.
    @param neighbor the other rigid body
    @param factor factor by which the crust thickness is increased for this
    particular rigid body collision
    @param id ID number of 1st component
    @param id_neighbor ID number of 2nd component
    @exception if the minimal distance between the shrunk rigid bodies with
    crust thickness artificially increased is 0, in practice less than EPSILON
    defined in Basic.H, which means that the shrunk rigid bodies with
    crust thickness artificially already touch or overlap */
    PointContact ClosestPoint_ErreurHandling(
    	RigidBodyWithCrust const& neighbor, double const& factor, int const& id,
	int const& id_neighbor )
      throw(ContactError);

    /** @brief Returns whether the rigid body is close to another rigid body in
    the sense of whether their respective bounding boxes minus their crust
    thickness overlap
    @param neighbor the other rigid body */
    bool isClose( RigidBodyWithCrust const& neighbor ) const;

    /** @brief Returns whether there is geometric contact with another rigid
    body in the sense of ClosestPoint, i.e., if minimal distance
    between the shrunk rigid bodies < sum of the crust thicknesses
    @param neighbor the other rigid body */
    bool isContact( RigidBodyWithCrust& neighbor );

    /** @brief Returns whether the rigid body is close to another rigid body in
    the sense of whether their respective bounding boxes plus their crust
    thickness overlap. Slightly different from isClose, needs further
    clarification.
    @param a 1st rigid body
    @param b 2nd rigid body */
    friend bool intersect( RigidBodyWithCrust const& a,
    	RigidBodyWithCrust const& b );
    //@}


    /**@name Get methods */
    //@{
    /** @brief Returns the crust thickness */
    double getCrustThickness() const;

    /** @brief Returns a pointer to the rigid body's transformation with the
    scaling by the crust thickness to shrink the rigid body */
    Transform const* getTransformWithCrust();

    /** @brief Returns the rigid body's transformation with the
    scaling by the crust thickness multiplied by a factor > 1 to shrink the
    rigid body. The shrinkage is capped by a minimum scaling usually set to a
    value between 0.5 and 1 when this method is called
    @param factor factor by which the crust thickness is increased
    @param min_scaling the rigid body cannot be shrunk by more than this value
    usually set between 0.5 and 1 */
    Transform getTransformWithCrust( double const& factor,
	double const& min_scaling = 0.5 ) const;
    //@}


    /**@name Set methods */
    //@{
    /** @brief Sets the boolean that tells that the rigid body's transformation
    with the scaling by the crust thickness to shrink the rigid bodies has
    already been computed to false */
    void initialize_transformWithCrust_to_notComputed();

    /** @brief Sets the crust thickness
    @param cthickness_ crust thickness */
    void setCrustThickness( double cthickness_ );
    //@}


    /**@name I/O methods */
    //@{
    /** @brief Writes the rigid body's "static" data, i.e., the convex geometric
    description only (without any transformation)
    @param fileOut output stream */
    void writeStatic( ostream& fileOut )
    	const;
    //@}


protected:
  //@name Parameters */
  //@{
  double m_crustThickness; /**< crust thickness */
  Vector3* m_scaling; /**< Diagonal coefficients of the scaling matrix related
  	to the crust thickness */
  Transform* m_transformWithCrust; /** Transformation corresponding to a
  	composition of (i) the scaling to shrink the rigid body by the crust
	thickness and (ii) the position of the center of mass and the angular
	position of the rigid body */
  bool m_transformWithCrust_computed; /** whether m_positionWithCrust has been
  	computed or not */
  //@}


private:

};

/** @brief Returns the features of the contact when the 2 rigid bodies are
spheres, i.e., a SPHERE-SPHERE contact
@param rbA 1st rigid body
@param rbB 2nd rigid body */
PointContact ClosestPointSPHERE( RigidBodyWithCrust const& rbA,
	RigidBodyWithCrust const& rbB )
  throw(ContactError);

/** @brief Returns the features of the contact when the 1 rigid body is a sphere
and the other rigid body is a box, i.e., a SPHERE-BOX contact
@param rbA 1st rigid body
@param rbB 2nd rigid body */
PointContact ClosestPointSPHEREBOX( RigidBodyWithCrust const& rbA,
	RigidBodyWithCrust const& rbB )
  throw(ContactError);

/** @brief Returns whether there is geometric contact with another rigid body
in the sense of ClosestPoint when the 2 rigid bodies are spheres, i.e., a
SPHERE-SPHERE contact
@param rbA 1st rigid body
@param rbB 2nd rigid body */
bool isContactSPHERE( RigidBodyWithCrust const& rbA,
	RigidBodyWithCrust const& rbB );

/** @brief Returns whether there is geometric contact with another rigid body
in the sense of ClosestPoint when 1 rigid body is a sphere and the other rigid
body is a box, i.e., a SPHERE-BOX contact
@param rbA 1st rigid body
@param rbB 2nd rigid body */
bool isContactSPHEREBOX( RigidBodyWithCrust const& rbA,
	RigidBodyWithCrust const& rbB );

/** @brief Returns the features of the contact when the 2 rigid bodies
are cylinders
@param rbA 1st rigid body
@param rbB 2nd rigid body */
PointContact ClosestPointCYLINDERS( RigidBodyWithCrust const& rbA,
  RigidBodyWithCrust const& rbB )
  throw(ContactError);

/** @brief Returns whether there is a contact between the circumscribed
cylinders of two rigid bodies
@param rbA 1st rigid body
@param rbB 2nd rigid body */
bool isContactCYLINDERS( RigidBodyWithCrust const& rbA,
  RigidBodyWithCrust const& rbB );

/** @brief Returns the features of the contact when the 1 rigid body is
a rectangle
@param rbA 1st rigid body
@param rbB 2nd rigid body */
PointContact ClosestPointRECTANGLE( RigidBodyWithCrust const& rbA ,
  RigidBodyWithCrust const& rbB )
  throw(ContactError);

#endif
