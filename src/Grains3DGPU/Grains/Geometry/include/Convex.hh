#ifndef _CONVEX_HH_
#define _CONVEX_HH_

#include "Transform3.hh"

// Convex types
enum ConvexType
{
    SPHERE,
    BOX,
    CYLINDER,
    CONE,
    SUPERQUADRIC,
    POLYHEDRON,
    RECTANGLE
};

// =================================================================================================
/** @brief The class Convex. Convex bodies - The base class for various particle shapes.

    @author A.Yazdani - 2023 - Construction
    @author A.Yazdani - 2024 - Modificiation */
// =================================================================================================
template <typename T>
class Convex
{
protected:
    /**@name Contructors */
    //@{
    /** @brief Default constructor (forbidden except in derived classes) */
    __HOSTDEVICE__
    Convex();
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Destructor */
    __HOSTDEVICE__
    virtual ~Convex();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Returns the convex type */
    __HOSTDEVICE__
    virtual ConvexType getConvexType() const = 0;
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a clone of the convex */
    __HOSTDEVICE__
    virtual Convex<T>* clone() const = 0;

    /** @brief Returns the volume of the convex shape */
    __HOSTDEVICE__
    virtual T computeVolume() const = 0;

    /** @brief Computes the diagonal inertia tensor
        @param inertia diagonal inertia tensor (3 components: Ixx, Iyy, Izz) */
    __HOSTDEVICE__
    virtual void computeInertia(T (&inertia)[3]) const = 0;

    /** @brief Computes and returns the circumscribed radius of the reference convex shape */
    __HOSTDEVICE__
    virtual T computeCircumscribedRadius() const = 0;

    /** @brief Returns the half-length of the bounding box fitted to the convex without considering
        the transformation */
    __HOSTDEVICE__
    virtual Vector3<T> computeBoundingBox() const = 0;

    /** @brief Returns the tightest bounding cylinder fitted to the convex in body-local frame.
        The returned Vector3 encodes: [0] radius, [1] half-height, [2] axis index (0=X, 1=Y, 2=Z).
     */
    __HOSTDEVICE__
    virtual Vector3<T> computeBoundingCylinder() const = 0;

    /** @brief Convex support function, returns the support point P, i.e. the point on the surface
        of the convex shape that satisfies max(P.v)
        @param v direction vector */
    __HOSTDEVICE__
    virtual Vector3<T> support(const Vector3<T>& v) const = 0;

    /** @brief Eroded convex support function, returns the support point of the
        shape shrunk inward by crust along v: support(v) - crust/norm(v) * v
        @param v direction vector
        @param crust erosion thickness
        @param invNorm precomputed 1/norm(v) to avoid redundant sqrt */
    __HOSTDEVICE__
    Vector3<T> support(const Vector3<T>& v, T crust, T invNorm) const;

    /** @brief Returns whether point p lies in the convex shape
        @param p point */
    __HOSTDEVICE__
    virtual bool isInside(const Vector3<T>& p) const;
    //@}

    /** @name I/O methods */
    //@{
    /** @brief Input operator
        @param fileIn input stream */
    __HOST__
    virtual void readConvex(std::istream& fileIn) = 0;

    /** @brief Output operator
        @param fileOut output stream */
    __HOST__
    virtual void writeConvex(std::ostream& fileOut) const = 0;

    /** @brief Returns the number of points to write the convex in a Paraview format */
    __HOST__
    virtual int numberOfPoints_PARAVIEW() const = 0;

    /** @brief Returns the number of elementary polytopes to write the convex in a Paraview format
     */
    __HOST__
    virtual int numberOfCells_PARAVIEW() const = 0;

    /** @brief Writes the list of points describing the convex to an stream
        @param f output stream
        @param transform geometric transformation
        @param translation additional center of mass translation */
    __HOST__
    void writePoints_PARAVIEW(std::ostream&        f,
                              const Transform3<T>& transform,
                              const Vector3<T>*    translation = NULL) const;

    /** @brief Returns a list of points describing the convex in a Paraview format
        @param transform geometric transformation
        @param translation additional center of mass translation */
    __HOST__
    virtual std::list<Vector3<T>> writePoints_PARAVIEW(const Transform3<T>& transform,
                                                       Vector3<T> const*    translation) const
        = 0;

    /** @brief Writes the connectivity of the convex in a Paraview format
        @param connectivity connectivity of Paraview polytopes
        @param offsets connectivity offsets
        @param cellstype Paraview polytopes type
        @param firstpoint_globalnumber global number of the 1st point
        @param last_offset last offset used for the previous convex shape */
    __HOST__
    virtual void writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                          std::list<uint>& offsets,
                                          std::list<uint>& cellstype,
                                          uint&            firstpoint_globalnumber,
                                          uint&            last_offset) const
        = 0;
    //@}
};

/** @name External Methods - I/O methods */
//@{
/** @brief Output operator for Convex: delegates to virtual writeConvex
    @param fileOut output stream
    @param convex convex object */
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Convex<T>& convex);

/** @brief Input operator for Convex: delegates to virtual readConvex
    @param fileIn input stream
    @param convex convex object */
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Convex<T>& convex);
//@}

#endif