#ifndef _BCYLINDER_HH_
#define _BCYLINDER_HH_

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
                       Vector3 const& e, Point3 const& x, PointContact& ptCont);

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

/** @brief Returns the contact point of two cylinders in the global world if
the contact is Edge-Edge.
@param rA cylinder A radius
@param hA cylinder A height
@param rB cylinder B radius
@param hB cylinder B height
@param e1 orientation of cylinder B w.r.t. to cylinder A - A is along Z axis
@param x1 position of cylinder B center w.r.t. to cylinder A center (origin)
@param ptCont contact variables */
void E2EContact( double rA, double hA, double rB, double hB,
                 Vector3 const& e1, Point3 const& x1, PointContact& ptCont );

// /** @brief Returns the real solutions to the quartic equation
// x^4 + bx^3 + cx^2 + dx + e = 0
// @param b coefficient of x^3
// @param c coefficient of x^2
// @param d coefficient of x^1
// @param e coefficient of 1
// @param sol array of solutions
// @param nbRoots number of real roots */
// void solveQuartic( double const b, double const c, double const d,
//                    double const e, double sol[4], int& nbRoots );

// /** @brief Returns the solutions to the quadrature ax^2 + bx + c = 0
// @param a coefficient of x^2
// @param b coefficient of x
// @param c coefficient of 1
// @param sol array of solutions */
// void solveQuadratic( double const a, double const b, double const c,
//                      double sol[2] );
//
// /** @brief Returns the rotation matrix, trasforming v to z-axis
// @param v the source vector whose transformation to the z-axis is needed
// @param rotMat the rotation matrix */
// void rotateVec2VecZ( Vector3 const& v, Matrix& rotMat );
//
// /** @brief Type-proof sign function
// @param val parameter whose sign is needed */
// template < typename T > int sgn( T val );
// //@}
#endif
