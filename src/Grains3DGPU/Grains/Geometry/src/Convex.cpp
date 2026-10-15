#include "Convex.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ Convex<T>::Convex()
{
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Convex<T>::~Convex()
{
}

// -------------------------------------------------------------------------------------------------
// Eroded convex support function
template <typename T>
__HOSTDEVICE__ Vector3<T> Convex<T>::support(const Vector3<T>& v, T crust, T invNorm) const
{
    return support(v) - (crust * invNorm) * v;
}

// -------------------------------------------------------------------------------------------------
// Returns whether point p lies in the convex shape
// @param p point
template <typename T>
__HOSTDEVICE__ bool Convex<T>::isInside(const Vector3<T>& p) const
{
    // Default implementation (for convex shapes that are not defined)
    return true;
}

// -------------------------------------------------------------------------------------------------
//
template <typename T>
__HOST__ void Convex<T>::writePoints_PARAVIEW(std::ostream&        f,
                                              const Transform3<T>& transform,
                                              const Vector3<T>*    translation) const
{
    std::list<Vector3<T>> points = writePoints_PARAVIEW(transform, translation);

    for(auto& point : points)
    {
        f << point[X] << " " << point[Y] << " " << point[Z] << std::endl;
    }
}

// -------------------------------------------------------------------------------------------------
// Output operator for Convex
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Convex<T>& convex)
{
    convex.writeConvex(fileOut);
    return fileOut;
}

// -------------------------------------------------------------------------------------------------
// Input operator for Convex
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Convex<T>& convex)
{
    convex.readConvex(fileIn);
    return fileIn;
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Convex<float>;
template class Convex<double>;

#define X(T)                                                                \
    template std::ostream& operator<< <T>(std::ostream&, const Convex<T>&); \
    template std::istream& operator>> <T>(std::istream&, Convex<T>&);
X(float)
X(double)
#undef X