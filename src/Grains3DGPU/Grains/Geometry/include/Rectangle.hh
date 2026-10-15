#ifndef _RECTANGLE_HH_
#define _RECTANGLE_HH_

#include "Convex.hh"
#include "ReaderXML.hh"

// =================================================================================================
/** @brief The class Rectangle.

    Convex with the shape of a rectangle.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class Rectangle : public Convex<T>
{
protected:
    /** @name Parameters */
    //@{
    T m_LX; /**< half-legnth of the X edge */
    T m_LY; /**< half-legnth of the Y edge */
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with the dimensions of the rectangle
        @param x 1st component
        @param y 2nd component */
    __HOSTDEVICE__
    Rectangle(T x = T(0), T y = T(0));

    /** @brief Constructor with an input stream
        @param fileIn input stream */
    __HOST__
    Rectangle(std::istream& fileIn);

    /** @brief Constructor with an XML node as an input parameter
        @param root XML node */
    __HOST__
    Rectangle(DOMNode* root);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Rectangle();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Returns the convex type */
    __HOSTDEVICE__
    ConvexType getConvexType() const final;

    /** @brief Returns the edge lengths in a Vector3 format with Z = 0 */
    __HOSTDEVICE__
    Vector3<T> getExtent() const;
    //@}

    /** @name Set methods */
    //@{
    /** @brief Sets values of the edge length
        @param x 1st component
        @param y 2nd component */
    __HOSTDEVICE__
    void setExtent(T x, T y);
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a clone of the rectangle */
    __HOSTDEVICE__
    Convex<T>* clone() const final;

    /** @brief Returns the rectangle volume (area) */
    __HOSTDEVICE__
    T computeVolume() const final;

    /** @brief Computes the inertia tensor
        @param inertia inertia tensor */
    __HOSTDEVICE__
    void computeInertia(T (&inertia)[3]) const final;

    /** @brief Returns the circumscribed radius of the box */
    __HOSTDEVICE__
    T computeCircumscribedRadius() const final;

    /** @brief Returns the half-length of the bounding box fitted to the rectangle without
        considering the transformation */
    __HOSTDEVICE__
    Vector3<T> computeBoundingBox() const final;

    /** @brief Returns the tightest bounding cylinder fitted to the rectangle in body-local frame.
     */
    __HOSTDEVICE__ Vector3<T> computeBoundingCylinder() const final;

    /** @brief Rectnagle support function, returns the support point P, i.e.
        the point on the surface of the rectangle that satisfies max(P.v)
        @param v direction */
    __HOSTDEVICE__
    Vector3<T> support(const Vector3<T>& v) const final;

    /** @brief Returns whether point p lies in the rectangle
        @param p point */
    __HOSTDEVICE__
    bool isInside(const Vector3<T>& p) const final;
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

    /** @brief Returns the number of points to write the rectangle in a Paraview format. */
    __HOST__
    int numberOfPoints_PARAVIEW() const final;

    /** @brief Returns the number of elementary polytopes to write the rectangle in a Paraview
        format. */
    __HOST__
    int numberOfCells_PARAVIEW() const final;

    /** @brief Returns a list of points describing the rectangle in a Paraview format.
        @param transform geometric transformation
        @param translation additional center of mass translation */
    __HOST__
    std::list<Vector3<T>> writePoints_PARAVIEW(const Transform3<T>& transform,
                                               Vector3<T> const*    translation) const final;

    /** @brief Writes the connectivity of the rectangle in a Paraview format
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
