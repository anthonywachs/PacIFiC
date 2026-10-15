#ifndef _GJK_SHAPEDATA_HH_
#define _GJK_SHAPEDATA_HH_

#include "Basic.hh"
#include "Box.hh"
#include "Cone.hh"
#include "Convex.hh"
#include "Cylinder.hh"
#include "MiscMath.hh"
#include "Rectangle.hh"
#include "RigidBody.hh"
#include "Sphere.hh"
#include "Superquadric.hh"
#include "Vector3.hh"
#include "VectorMath.hh"

/** @name Structs */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief POD struct holding all data needed to evaluate a convex support function on the GPU
    without virtual dispatch. Layout of params[] by shape type:
      SPHERE:         params[0] = radius
      BOX:            params[0] = ex, [1] = ey, [2] = ez
      CYLINDER:       params[0] = radius, [1] = halfHeight
      CONE:           params[0] = bottomRadius, [1] = quarterHeight, [2] = sinAngle
      SUPERQUADRIC:   params[0] = a, [1] = b, [2] = c, [3] = n1, [4] = n2
      RECTANGLE:      params[0] = LX, [1] = LY
    crust = crustThickness  (per-body full crust, matching GJK erosion convention). */
template <typename T>
struct ShapeData
{
    ConvexType type;
    T          crust;                // crustThickness (per-body full crust)
    T          params[5];            // shape-specific support params
    T          circumscribedRadius;  // circumscribed sphere radius
    T          invMass;              // 1/mass (0 for fixed/obstacle bodies)
    uint       material;             // material ID for contact hash
};

// -------------------------------------------------------------------------------------------------
/** @brief POD struct holding bounding-volume data for the GPU BV pre-filter, stored separately
    from ShapeData so that the GJK kernel does not carry unused BV fields in registers.
    boundingBox      = half-extents {dx, dy, dz}  (from computeBoundingBox())
    boundingCylinder = {radius, halfHeight, axisIndex}  (from computeBoundingCylinder()) */
template <typename T>
struct BVData
{
    Vector3<T> boundingBox;          // half-extents {dx, dy, dz}
    Vector3<T> boundingCylinder;     // {radius, halfHeight, axisIndex(0=X,1=Y,2=Z)}
    T          circumscribedRadius;  // circumscribed sphere radius (for fast reject)
};
//@}

// =================================================================================================
/** @brief Utility functions to fill ShapeData and BVData from RigidBody/Convex using virtual
    dispatch once at init time, and to evaluate the raw support function for a ShapeData on the
    device without virtual dispatch.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
/** @name Structs */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Fills a ShapeData struct from a RigidBody, reading the concrete shape geometry.
    Uses virtual dispatch once on the init path (not in the GJK hot loop).
    @param sd output ShapeData to fill
    @param rb source RigidBody (must not be null) */
template <typename T>
__HOSTDEVICE__ void fillShapeData(ShapeData<T>& sd, const RigidBody<T>* rb)
{
    const Convex<T>* conv = rb->getConvex();
    sd.type               = conv->getConvexType();
    sd.crust              = rb->getCrustThickness();
    switch(sd.type)
    {
    case SPHERE:
        sd.params[0] = static_cast<const Sphere<T>*>(conv)->getRadius();
        break;
    case BOX:
    {
        const Vector3<T> e = static_cast<const Box<T>*>(conv)->getExtent();
        sd.params[0]       = e[X];
        sd.params[1]       = e[Y];
        sd.params[2]       = e[Z];
        break;
    }
    case CYLINDER:
    {
        const Cylinder<T>* c = static_cast<const Cylinder<T>*>(conv);
        sd.params[0]         = c->getRadius();
        sd.params[1]         = c->getHeight() / T(2);
        break;
    }
    case CONE:
    {
        const Cone<T>* c = static_cast<const Cone<T>*>(conv);
        sd.params[0]     = c->getRadius();
        sd.params[1]     = c->getHeight() / T(4);
        sd.params[2]     = c->getSinAngle();
        break;
    }
    case SUPERQUADRIC:
    {
        const Superquadric<T>* sq = static_cast<const Superquadric<T>*>(conv);
        const Vector3<T>       ex = sq->getExtent();
        const Vector3<T>       en = sq->getExponent();
        sd.params[0]              = ex[X];
        sd.params[1]              = ex[Y];
        sd.params[2]              = ex[Z];
        sd.params[3]              = en[X];
        sd.params[4]              = en[Y];
        break;
    }
    case RECTANGLE:
    {
        const Vector3<T> e = static_cast<const Rectangle<T>*>(conv)->getExtent();
        sd.params[0]       = e[X];
        sd.params[1]       = e[Y];
        break;
    }
    default:
        break;
    }
    sd.circumscribedRadius = rb->getCircumscribedRadius();
    const T mass           = rb->getMass();
    sd.invMass             = (mass == T(0)) ? T(0) : T(1) / mass;
    sd.material            = rb->getMaterial();
}

// -------------------------------------------------------------------------------------------------
/** @brief Fills a BVData struct from a RigidBody/Convex using virtual dispatch once at init time.
    @param bv  output BVData to fill
    @param rb  source RigidBody (must not be null) */
template <typename T>
__HOSTDEVICE__ void fillBVData(BVData<T>& bv, const RigidBody<T>* rb)
{
    const Convex<T>* conv  = rb->getConvex();
    bv.boundingBox         = conv->computeBoundingBox();
    bv.boundingCylinder    = conv->computeBoundingCylinder();
    bv.circumscribedRadius = rb->getCircumscribedRadius();
}

// -------------------------------------------------------------------------------------------------
/** @brief Evaluates the raw support function for a ShapeData (without crust erosion).
    Exactly mirrors the virtual support(v) in each concrete Convex subclass.
    @param sd ShapeData describing the shape
    @param v  search direction (need not be normalised)
    @return   support point in shape-local frame */
template <typename T>
__HOSTDEVICE__ Vector3<T> device_support_raw(const ShapeData<T>& sd, const Vector3<T>& v)
{
    switch(sd.type)
    {
    case SPHERE:
    {
        const T len = norm(v);
        if(len < EPS<T>)
            return Vector3<T>(T(0), T(0), T(0));
        return (sd.params[0] / len) * v;
    }
    case BOX:
    {
        return Vector3<T>(v[X] < T(0) ? -sd.params[0] : sd.params[0],
                          v[Y] < T(0) ? -sd.params[1] : sd.params[1],
                          v[Z] < T(0) ? -sd.params[2] : sd.params[2]);
    }
    case CYLINDER:
    {
        const T radius     = sd.params[0];
        const T halfHeight = sd.params[1];
        const T s          = sqrt(v[X] * v[X] + v[Z] * v[Z]);
        if(s > EPS<T>)
        {
            const T d = radius / s;
            if(fabs(v[Y]) < EPS<T>)
                return Vector3<T>(v[X] * d, T(0), v[Z] * d);
            else
                return Vector3<T>(v[X] * d, v[Y] < T(0) ? -halfHeight : halfHeight, v[Z] * d);
        }
        else
            return Vector3<T>(T(0), v[Y] < T(0) ? -halfHeight : halfHeight, T(0));
    }
    case CONE:
    {
        const T bottomRadius  = sd.params[0];
        const T quarterHeight = sd.params[1];
        const T sinAngle      = sd.params[2];
        if(v[Y] > norm(v) * sinAngle)
            return Vector3<T>(T(0), T(3) * quarterHeight, T(0));
        const T s = sqrt(v[X] * v[X] + v[Z] * v[Z]);
        if(s > EPS<T>)
        {
            const T d = bottomRadius / s;
            return Vector3<T>(v[X] * d, -quarterHeight, v[Z] * d);
        }
        return Vector3<T>(T(0), -quarterHeight, T(0));
    }
    case SUPERQUADRIC:
    {
        const T a     = sd.params[0];
        const T b     = sd.params[1];
        const T c     = sd.params[2];
        const T n1    = sd.params[3];
        const T n2    = sd.params[4];
        const T abvx  = fabs(v[X]);
        const T abvy  = fabs(v[Y]);
        const T abvz  = fabs(v[Z]);
        const T signx = T(sgn(v[X]));
        const T signy = T(sgn(v[Y]));
        const T signz = T(sgn(v[Z]));
        if(abvx == T(0))
        {
            if(abvy == T(0))
                return Vector3<T>(T(0), T(0), signz * c);
            const T alpha = pow(c / b * abvz / abvy, T(1) / (n1 - T(1)));
            const T yt    = T(1) / pow(T(1) + pow(alpha, n1), T(1) / n1);
            return Vector3<T>(T(0), signy * b * yt, signz * alpha * c * yt);
        }
        const T alpha = pow(b / a * abvy / abvx, T(1) / (n2 - T(1)));
        const T temp  = T(1) + pow(alpha, n2);
        const T gamma = pow(temp, (n1 - n2) / (n2 * (n1 - T(1))));
        const T beta  = gamma * pow(c / a * abvz / abvx, T(1) / (n1 - T(1)));
        const T xt    = T(1) / pow(pow(temp, n1 / n2) + pow(beta, n1), T(1) / n1);
        return Vector3<T>(signx * a * xt, signy * alpha * b * xt, signz * beta * c * xt);
    }
    case RECTANGLE:
    {
        return Vector3<T>(v[X] < T(0) ? -sd.params[0] : sd.params[0],
                          v[Y] < T(0) ? -sd.params[1] : sd.params[1],
                          T(0));
    }
    default:
        return Vector3<T>(T(0), T(0), T(0));
    }
}

// -------------------------------------------------------------------------------------------------
/** @brief Evaluates the crust-eroded support function using ShapeData (vtable-free).
    Equivalent to Convex<T>::support(v, crustArg, invNorm) but without virtual dispatch.
    @param sd      ShapeData describing the shape
    @param v       search direction (need not be normalised)
    @param crustArg crust thickness (per-body full crustThickness)
    @param invNorm 1 / norm(v), pre-computed by the caller
    @return        crust-eroded support point in shape-local frame */
template <typename T>
__HOSTDEVICE__ Vector3<T>
               device_support(const ShapeData<T>& sd, const Vector3<T>& v, T crustArg, T invNorm)
{
    return device_support_raw(sd, v) - crustArg * invNorm * v;
}

// -------------------------------------------------------------------------------------------------
/** @brief Fills compact ShapeData and BVData tables indexed by shapeId.
    One entry per unique shape prototype (size = nUniqueShapes, NOT nComponents).
    repSlots[k] is any component slot whose shapeId == k -- all slots with the same shapeId
    share identical RigidBody* and produce identical entries.
    Must be called on HOST with HOST-resident rigidBodies and repSlots.
    @param sdOut         output ShapeData array (size >= nUniqueShapes)
    @param bvOut         output BVData array (size >= nUniqueShapes)
    @param rbs           per-component RigidBody* array (size >= nComponents)
    @param repSlots      representative slot index for each unique shapeId (size >= nUniqueShapes)
    @param nUniqueShapes number of unique shape prototypes (max shapeId + 1) */
template <typename T>
__HOST__ void buildShapeAndBVData(ShapeData<T>*              sdOut,
                                  BVData<T>*                 bvOut,
                                  const RigidBody<T>* const* rbs,
                                  const uint*                repSlots,
                                  uint                       nUniqueShapes)
{
    for(uint k = 0; k < nUniqueShapes; ++k)
    {
        fillShapeData(sdOut[k], rbs[repSlots[k]]);
        fillBVData(bvOut[k], rbs[repSlots[k]]);
    }
}
//@}

#endif
