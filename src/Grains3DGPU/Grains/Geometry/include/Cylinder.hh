#ifndef _CYLINDER_HH_
#define _CYLINDER_HH_

#include "Convex.hh"
#include "ReaderXML.hh"

// =================================================================================================
/** @brief The class Cylinder.

    Convex with the shape of a cylinder.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class Cylinder : public Convex<T>
{
protected:
    /** @name Parameters */
    //@{
    T m_radius;     /**< cylinder radius */
    T m_halfHeight; /**< cylinder half height */
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with inputs
        @param r cylinder radius
        @param h cylinder height */
    __HOSTDEVICE__
    Cylinder(T r = T(0), T h = T(0));

    /** @brief Constructor with an input stream
        @param fileIn input stream */
    __HOST__
    Cylinder(std::istream& fileIn);

    /** @brief Constructor with an XML node as an input parameter
        @param root XML node */
    __HOST__
    Cylinder(DOMNode* root);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Cylinder();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the convex type */
    __HOSTDEVICE__
    ConvexType getConvexType() const final;

    /** @brief Gets the radius */
    __HOSTDEVICE__
    T getRadius() const;

    /** @brief Gets the height */
    __HOSTDEVICE__
    T getHeight() const;
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a clone of the cylinder */
    __HOSTDEVICE__
    Convex<T>* clone() const final;

    /** @brief Returns the cylinder volume */
    __HOSTDEVICE__
    T computeVolume() const final;

    /** @brief Computes the inertia tensor
        @param inertia inertia tensor */
    __HOSTDEVICE__
    void computeInertia(T (&inertia)[3]) const final;

    /** @brief Computes and returns the circumscribed radius of the cylinder */
    __HOSTDEVICE__
    T computeCircumscribedRadius() const final;

    /** @brief Returns the half-length of the bounding box fitted to the cylinder without
        considering the transformation */
    __HOSTDEVICE__
    Vector3<T> computeBoundingBox() const final;

    /** @brief Returns the tightest bounding cylinder fitted to the cylinder in body-local frame. */
    __HOSTDEVICE__
    Vector3<T> computeBoundingCylinder() const final;

    /** @brief Cylinder support function, returns the support point P, i.e. point on the surface of
        the cylinder that satisfies max(P.v)
        @param v direction */
    __HOSTDEVICE__
    Vector3<T> support(const Vector3<T>& v) const final;
    //@}

    /** @name I/O methods */
    //@{
    /** @brief Input operator
        @param fileIn input stream */
    __HOST__
    void readConvex(std::istream& fileIn) final;

    /** @brief Output operator
        @param fileOut output stream */
    __HOST__
    void writeConvex(std::ostream& fileOut) const final;

    /** @brief Returns the number of points to write the cylinder in a Paraview format. */
    __HOST__
    int numberOfPoints_PARAVIEW() const final;

    /** @brief Returns the number of elementary polytopes to write the cylinder in a Paraview
        format. */
    __HOST__
    int numberOfCells_PARAVIEW() const final;

    /** @brief Returns a list of points describing the cylinder in a Paraview format.
        @param transform geometric transformation
        @param translation additional center of mass translation */
    __HOST__
    std::list<Vector3<T>> writePoints_PARAVIEW(const Transform3<T>& transform,
                                               Vector3<T> const*    translation) const final;

    /** @brief Writes the connectivity of the cylinder in a Paraview format
        @param connectivity connectivity of Paraview polytopes
        @param offsets connectivity offsets
        @param cellstype Paraview polytopes type
        @param firstpoint_globalnumber global number of the 1st point
        @param last_offset last offset used for the previous convex shape */
    __HOST__
    void writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                  std::list<uint>& offsets,
                                  std::list<uint>& cellstype,
                                  uint&            firstpoint_globalnumber,
                                  uint&            last_offset) const final;
    //@}
};

#endif
