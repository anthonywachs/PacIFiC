#ifndef _BCYLINDER_HH_
#define _BCYLINDER_HH_

#include "Point3.hh"
#include "Vector3.hh"
#include "Matrix.hh"
#include "Transform.hh"
#include "PointContact.hh"
using namespace solid;

class BCylinder;
// ostream& operator << ( ostream& f, BCylinder const& B );

/** @brief The class BCylinder.

    Bounding cylinder oriented along the axis of the world reference frame

    @author A.YAZDANI - 2022 - Creation */
// ============================================================================
class BCylinder
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    BCylinder();

    /** @brief Constructor with radius, height, and axis
    @param r radius
    @param h height
    @param e axis */
    BCylinder( double r, double h, Vector3 const& v );

    /** @brief Copy constructor
    @param bcylinder_ reference bounding cylinder */
    BCylinder( BCylinder const& bcylinder_ );

    /** @brief Destructeur */
    ~BCylinder();
    //@}

    /**@name Methods */
    //@{
    /** @brief Returns the bounding cylinder radius */
    double getRadius() const;

    /** @brief Returns the bounding cylinder height */
    double getHeight() const;

    /** @brief Returns the bounding cylinder axis */
    Vector3 const& getAxis() const;

    /** @brief Sets the bounding cylinder radius
    @param r new radius */
    void setRadius( double r );

    /** @brief Sets the bounding cylinder height
    @param h new height */
    void setHeight( double h );

    /** @brief Sets the bounding cylinder axis
    @param v new axis */
    void setAxis( Vector3 const& v );
    //@}


    /** @name Friend methods */
    //@{
    /** @brief Returns the contact point of two cylinders
    @param a 1st bounding cylinder
    @param b 2nd bounding cylinder
    @param a2w transformation of the first cylinder
    @param b2w transformation of the second cylinder */
    friend PointContact intersect( BCylinder const& a, BCylinder const& b,
                                   Transform const& a2w, Transform const& b2w);

    /** @brief Returns whether 2 cylinders are in contact
    @param a 1st bounding cylinder
    @param b 2nd bounding cylinder
    @param a2w transformation of the first cylinder
    @param b2w transformation of the second cylinder */
    friend bool isContact( BCylinder const& a, BCylinder const& b,
                           Transform const& a2w, Transform const& b2w );

    /** @brief Output operator
    @param f output stream
    @param B BCylinder object */
    friend ostream& operator << ( ostream& f, BCylinder const& B );
    //@}


  private:
    /** @name Parameters */
    //@{
    double m_radius; /**< bounding cylinder radius */
    double m_height; /**< bounding cylinder height */
    Vector3 m_axis; /**< bounding cylinder axis */
    //@}
};


/** @name BCylinder : External methods */
//@{
/** @brief Wrapper function for analytical contacts between cylinders.
@param rA cylinder A radius
@param hA cylinder A height
@param rB cylinder B radius
@param hB cylinder B height
@param e_B2A orientation of cylinder B w.r.t. to cylinder A - A is along Z axis
@param x_B2A position of cylinder B center w.r.t. to cylinder A center (origin)
@param e_A2B orientation of cylinder A w.r.t. to cylinder B - B is along Z axis
@param x_A2B position of cylinder A center w.r.t. to cylinder B center (origin)
@param method type of contact: F2F, F2B, ...
@param ptCont contact variables */
void BCylinderContactWrapper( double rA, double hA, double rB, double hB,
                              Vector3 const& e_B2A, Point3 const& x_B2A,
                              Vector3 const& e_A2B, Point3 const& x_A2B,
                              int method, PointContact& ptCont );

/** @brief Returns the contact point of two cylinders in the world of cylidner A
if the contact is either Face-Face or Band-Band (Parallel). Cylinder A is at
origin oriented along Z-axis.
@param rA cylinder A radius
@param hA cylinder A height
@param rB cylinder B radius
@param hB cylinder B height
@param e orientation of cylinder B w.r.t. to cylinder A - A is along Z axis
@param x position of cylinder B center w.r.t. to cylinder A center (origin)
@param ptCont contact variables */
void F2FB2BParContact( double rA, double hA, double rB, double hB,
                       Vector3 const& e, Point3 const& x, PointContact& ptCont );

/** @brief Returns the contact point of two cylinders in the world of cylidner A
if the contact is Face-Band. Cylinder A is at origin oriented along Z-axis.
@param rA cylinder A radius
@param hA cylinder A height
@param rB cylinder B radius
@param hB cylinder B height
@param e orientation of cylinder B w.r.t. to cylinder A - A is along Z axis
@param x position of cylinder B center w.r.t. to cylinder A center (origin)
@param ptCont contact variables */
void F2BContact( double rA, double hA, double rB, double hB,
                 Vector3 const& e, Point3 const& x, PointContact& ptCont );


/** @brief Returns the contact point of two cylinders in the world of cylidner A
if the contact is Face-Edge. Cylinder A is at origin oriented along Z-axis.
@param rA cylinder A radius
@param hA cylinder A height
@param rB cylinder B radius
@param hB cylinder B height
@param e orientation of cylinder B w.r.t. to cylinder A - A is along Z axis
@param x position of cylinder B center w.r.t. to cylinder A center (origin)
@param ptCont contact variables */
void F2EContact( double rA, double hA, double rB, double hB,
                 Vector3 const& e, Point3 const& x, PointContact& ptCont );

/** @brief Returns the contact point of two cylinders in the world of cylidner A
if the contact is Band-Band (Skewed). Cylinder A is at origin oriented along
Z-axis.
@param rA cylinder A radius
@param hA cylinder A height
@param rB cylinder B radius
@param hB cylinder B height
@param e orientation of cylinder B w.r.t. to cylinder A - A is along Z axis
@param x position of cylinder B center w.r.t. to cylinder A center (origin)
@param ptCont contact variables */
void B2BSkewContact( double rA, double hA, double rB, double hB,
                     Vector3 const& e, Point3 const& x, PointContact& ptCont );

/** @brief Returns the contact point of two cylinders in the world of cylidner A
if the contact is Band-Edge. Cylinder A is at origin oriented along Z-axis.
@param rA cylinder A radius
@param hA cylinder A height
@param rB cylinder B radius
@param hB cylinder B height
@param e orientation of cylinder B w.r.t. to cylinder A - A is along Z axis
@param x position of cylinder B center w.r.t. to cylinder A center (origin)
@param ptCont contact variables */
void B2EContact( double rA, double hA, double rB, double hB,
                 Vector3 const& e, Point3 const& x, PointContact& ptCont );

/** @brief Returns the contact point of two cylinders in the world of cylidner A
if the contact is Edge-Edge. Cylinder A is at origin oriented along Z-axis.
@param rA cylinder A radius
@param hA cylinder A height
@param rB cylinder B radius
@param hB cylinder B height
@param e orientation of cylinder B w.r.t. to cylinder A - A is along Z axis
@param x position of cylinder B center w.r.t. to cylinder A center (origin)
@param ptCont contact variables */
void E2EContact( double rA, double hA, double rB, double hB,
                 Vector3 const& e, Point3 const& x, PointContact& ptCont );

/** @brief Returns the band (ptA) and Face (ptB) points of the given circle
// intersecting with a cylinder oriented along Z axis and centered at the
// origin with radius rA and height hA
@param */
void edgePointClose2Z( double rA, double r, Vector3 const& e, Point3 const& c,
                       Point3& ptA );

/** @brief Returns the band (ptA) and Face (ptB) points of the given circle
intersecting with a cylinder oriented along Z axis and centered at the  origin
with radius R and height h
@param rA radius of the master circle
@param hA height of the master circle
@param rB radius of the slave circle
@param hB height of the slave circle
@param e orientation of the slave circle
@param c center of the slave circle
@param ptA contact point with Band
@param ptB contact point with Face */
void edgePointsOnCyl( double rA, double hA, double rB, Vector3 const& e,
                      Point3 const& c, Point3& ptA, Point3& ptB );

/** @brief // Returns the solutions to a quartic equation
apearing in this class
@param a coefficient of x^4
@param b coefficient of x^3
@param c coefficient of x^2
@param d coefficient of x^1
@param e coefficient of 1
@param sol array of solutions */
void solveQuartic( double a, double b, double c, double d, double e,
                   double sol[4] );

/** @brief Returns a REAL solution to a cubic equation
apearing in this class
@param a coefficient of x^3
@param b coefficient of x^2
@param c coefficient of x^1
@param d coefficient of 1 */
double solveCubicReal( double a, double b, double c, double d );

/** @brief Returns a SPECIFIC solution for the quadrature ax^2 + bx + c = 0
apearing in this class
@param a coefficient of x^2
@param b coefficient of x
@param c coefficient of 1
@param sol array of solutions */
void solveQuadratic( double a, double b, double c, double sol[2] );

/** @brief Type-proof sign function
@param val parameter whose sign is needed */
template < typename T > int sgn( T val );
//@}
#endif
