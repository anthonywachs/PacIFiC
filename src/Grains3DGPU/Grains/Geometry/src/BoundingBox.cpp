#include "BoundingBox.hh"

// -------------------------------------------------------------------------------------------------
// Constructor
template <typename T>
__HOSTDEVICE__ BoundingBox<T>::BoundingBox()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with half edge length as input parameters
template <typename T>
__HOSTDEVICE__ BoundingBox<T>::BoundingBox(T x, T y, T z)
    : m_extent(Vector3<T>(x, y, z))
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with a vector containing the edge half-lengths
template <typename T>
__HOSTDEVICE__ BoundingBox<T>::BoundingBox(const Vector3<T>& extent_)
    : m_extent(extent_)
{
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ BoundingBox<T>::~BoundingBox()
{
}

// -------------------------------------------------------------------------------------------------
// Sets values of the edge length
template <typename T>
__HOSTDEVICE__ void BoundingBox<T>::setExtent(T x, T y, T z)
{
    m_extent = Vector3<T>(x, y, z);
}

// -------------------------------------------------------------------------------------------------
// Sets extent values using ccording to the convex which we want to fit the
// bounding box around.
template <typename T>
__HOSTDEVICE__ void BoundingBox<T>::setExtent(Convex<T> const& convex)
{
    m_extent = convex.computeBoundingBox();
}

// -------------------------------------------------------------------------------------------------
// Gets values of the edge length
template <typename T>
__HOSTDEVICE__ Vector3<T> const BoundingBox<T>::getExtent() const
{
    return (m_extent);
}

// -------------------------------------------------------------------------------------------------
// Returns a clone of the box
template <typename T>
__HOSTDEVICE__ BoundingBox<T>* BoundingBox<T>::clone() const
{
    return (new BoundingBox<T>(m_extent[X], m_extent[Y], m_extent[Z]));
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class BoundingBox<float>;
template class BoundingBox<double>;