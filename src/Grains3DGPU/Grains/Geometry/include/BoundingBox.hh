#ifndef _BOUNDINGBOX_HH_
#define _BOUNDINGBOX_HH_

#include "Convex.hh"

// Primary template
template <typename>
class Convex;

// =================================================================================================
/** @brief The class BoundingBox.

    Bounding Boxes (AABB/BoundingBox) collision detection for convex bodies.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class BoundingBox
{
protected:
    /** @name Parameters */
    //@{
    Vector3<T> m_extent; /**< vector containing the half-legnth of edges */
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    __HOSTDEVICE__
    BoundingBox();

    /** @brief Constructor with 3 components as inputs
        @param x 1st component
        @param y 2nd component
        @param z 3rd component */
    __HOSTDEVICE__
    BoundingBox(T x, T y, T z);

    /** @brief Constructor with a vector containing the edge half-lengths
        @param extent_ vector of half-lengths */
    __HOSTDEVICE__
    BoundingBox(const Vector3<T>& extent_);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~BoundingBox();
    //@}

    /** @name Methods */
    //@{
    /** @brief Sets values of the edge length
        @param x 1st component
        @param y 2nd component
        @param z 3rd component */
    __HOSTDEVICE__
    void setExtent(T x, T y, T z);

    /** @brief Sets the extent values according to the convex which we want to fit the BoundingBox
        around.
        @param convex convex object */
    __HOSTDEVICE__
    void setExtent(Convex<T> const& convex);

    /** @brief Gets values of the edge length
        @param x 1st component
        @param y 2nd component
        @param z 3rd component */
    __HOSTDEVICE__
    Vector3<T> const getExtent() const;

    /** @brief Returns a clone of the bounding box */
    __HOSTDEVICE__
    BoundingBox<T>* clone() const;
    //@}
};

typedef BoundingBox<float>  BoundingBoxF;
typedef BoundingBox<double> BoundingBoxD;

#endif