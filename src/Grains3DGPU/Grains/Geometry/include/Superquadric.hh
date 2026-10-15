#ifndef _SUPERQUADRIC_HH_
#define _SUPERQUADRIC_HH_

#include "Convex.hh"
#include "ReaderXML.hh"

// =================================================================================================
/** @brief The class Superquadric.

    Convex with the shape of a superquadric.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class Superquadric : public Convex<T>
{
protected:
    /** @name Parameters */
    //@{
    T m_a;  /**< scale parameter along the x-axis */
    T m_b;  /**< scale parameter along the y-axis */
    T m_c;  /**< scale parameter along the z-axis */
    T m_n1; /**< first exponent */
    T m_n2; /**< second exponent */
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with inputs
        @param a scale parameter along the x-axis
        @param b scale parameter along the y-axis
        @param c scale parameter along the z-axis
        @param n1 first exponent
        @param n2 second exponent */
    __HOSTDEVICE__
    Superquadric(T a = T(0), T b = T(0), T c = T(0), T n1 = T(2), T n2 = T(2));

    /** @brief Constructor with an input stream
        @param fileIn input stream */
    __HOST__
    Superquadric(std::istream& fileIn);

    /** @brief Constructor with an XML node as an input parameter
        @param root XML node */
    __HOST__
    Superquadric(DOMNode* root);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Superquadric();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the convex type */
    __HOSTDEVICE__
    ConvexType getConvexType() const final;

    /** @brief Gets the extent in a Vector3 format */
    __HOSTDEVICE__
    Vector3<T> getExtent() const;

    /** @brief Gets the exponents (blockiness) in a Vector3 format (Z=0) */
    __HOSTDEVICE__
    Vector3<T> getExponent() const;
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a clone of the superquadric */
    __HOSTDEVICE__
    Convex<T>* clone() const final;

    /** @brief Returns the Superquadric volume */
    __HOSTDEVICE__
    T computeVolume() const final;

    /** @brief Computes the inertia tensor
        @param inertia inertia tensor */
    __HOSTDEVICE__
    void computeInertia(T (&inertia)[3]) const final;

    /** @brief Computes and returns the circumscribed radius of the Superquadric */
    __HOSTDEVICE__
    T computeCircumscribedRadius() const final;

    /** @brief Returns the half-length of the bounding box fitted to the superquadric without
        considering the transformation */
    __HOSTDEVICE__
    Vector3<T> computeBoundingBox() const final;

    /** @brief Returns the tightest bounding cylinder fitted to the superquadric in body-local
     * frame. */
    __HOSTDEVICE__
    Vector3<T> computeBoundingCylinder() const final;

    /** @brief Superquadric support function, returns the support point P, i.e. point on the surface
        of the Superquadric that satisfies max(P.v)
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

    /** @brief Returns the number of points to write the superquadric in a Paraview format */
    __HOST__
    int numberOfPoints_PARAVIEW() const final;

    /** @brief Returns the number of elementary polytopes to write the superquadric in a Paraview
        format */
    __HOST__
    int numberOfCells_PARAVIEW() const final;

    /** @brief Returns a list of points describing the superquadric in a Paraview format
        @param transform geometric transformation
        @param translation additional center of mass translation */
    __HOST__
    std::list<Vector3<T>> writePoints_PARAVIEW(const Transform3<T>& transform,
                                               Vector3<T> const*    translation) const final;

    /** @brief Writes the connectivity of the superquadric in a Paraview format
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
