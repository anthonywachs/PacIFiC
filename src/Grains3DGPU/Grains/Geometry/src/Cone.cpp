#include "Cone.hh"
#include "VectorMath.hh"

// multiple of 4
#define visuNodeNbOnPer 32

// -------------------------------------------------------------------------------------------------
// Constructor with radius and height as input parameters
template <typename T>
__HOSTDEVICE__ Cone<T>::Cone(T r, T h)
    : m_bottomRadius(r)
    , m_quarterHeight(h / T(4))
    , m_sinAngle(r / sqrt(r * r + h * h))
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with an input stream
template <typename T>
__HOST__ Cone<T>::Cone(std::istream& fileIn)
{
    readConvex(fileIn);
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
template <typename T>
__HOST__ Cone<T>::Cone(DOMNode* root)
{
    m_bottomRadius  = T(ReaderXML::getNodeAttr_Double(root, "Radius"));
    m_quarterHeight = T(ReaderXML::getNodeAttr_Double(root, "Height")) / T(4);
    m_sinAngle
        = m_bottomRadius
          / sqrt(m_bottomRadius * m_bottomRadius + T(16) * m_quarterHeight * m_quarterHeight);
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ Cone<T>::~Cone()
{
}

// -------------------------------------------------------------------------------------------------
// Returns the convex type
template <typename T>
__HOSTDEVICE__ ConvexType Cone<T>::getConvexType() const
{
    return (CONE);
}

// -------------------------------------------------------------------------------------------------
// Returns the radius
template <typename T>
__HOSTDEVICE__ T Cone<T>::getRadius() const
{
    return (m_bottomRadius);
}

// -------------------------------------------------------------------------------------------------
// Returns the height
template <typename T>
__HOSTDEVICE__ T Cone<T>::getHeight() const
{
    return (T(4) * m_quarterHeight);
}

// -------------------------------------------------------------------------------------------------
// Returns the sine of the half-angle
template <typename T>
__HOSTDEVICE__ T Cone<T>::getSinAngle() const
{
    return (m_sinAngle);
}

// -------------------------------------------------------------------------------------------------
// Returns a clone of the cone
template <typename T>
__HOSTDEVICE__ Convex<T>* Cone<T>::clone() const
{
    return (new Cone<T>(m_bottomRadius, T(4) * m_quarterHeight));
}

// -------------------------------------------------------------------------------------------------
// Returns the volume of the cone
template <typename T>
__HOSTDEVICE__ T Cone<T>::computeVolume() const
{
    return (T(4) * m_quarterHeight * PI<T> * m_bottomRadius * m_bottomRadius / T(3));
}

// -------------------------------------------------------------------------------------------------
// Computes the diagonal inertia tensor
template <typename T>
__HOSTDEVICE__ void Cone<T>::computeInertia(T (&inertia)[3]) const
{
    T r2 = m_bottomRadius * m_bottomRadius;
    T c  = T(.2) * PI<T> * m_quarterHeight * r2;
    // Diagonal components: Ixx = Izz (perpendicular to axis), Iyy (along axis)
    inertia[0] = inertia[2] = c * (T(4) * m_quarterHeight * m_quarterHeight + r2);
    inertia[1]              = T(2) * c * r2;
}

// -------------------------------------------------------------------------------------------------
// Returns the circumscribed radius of the cone
template <typename T>
__HOSTDEVICE__ T Cone<T>::computeCircumscribedRadius() const
{
    return (max(sqrt(m_bottomRadius * m_bottomRadius + m_quarterHeight * m_quarterHeight),
                T(3) * m_quarterHeight));
}

// -------------------------------------------------------------------------------------------------
// Returns the bounding box to the cone
// TODO: OBB should be shifted one unit to the right because of center of mass
template <typename T>
__HOSTDEVICE__ Vector3<T> Cone<T>::computeBoundingBox() const
{
    return (Vector3<T>(m_bottomRadius, T(2) * m_quarterHeight, m_bottomRadius));
}

// -------------------------------------------------------------------------------------------------
// Returns the bounding cylinder to Cone
template <typename T>
__HOSTDEVICE__ Vector3<T> Cone<T>::computeBoundingCylinder() const
{
    // [radius, halfHeight, axisIndex=Y(1)]
    return Vector3<T>(m_bottomRadius, T(2) * m_quarterHeight, T(1));
}

// -------------------------------------------------------------------------------------------------
// Cone support function, returns the support point P, i.e. the point on the
// surface of the cone that satisfies max(P.v)
template <typename T>
__HOSTDEVICE__ Vector3<T> Cone<T>::support(const Vector3<T>& v) const
{
    if(v[Y] > norm(v) * m_sinAngle)
        return (Vector3<T>(T(0), T(3) * m_quarterHeight, T(0)));
    else
    {
        T s = sqrt(v[X] * v[X] + v[Z] * v[Z]);
        if(s > EPS<T>)
        {
            T d = m_bottomRadius / s;
            return (Vector3<T>(v[X] * d, -m_quarterHeight, v[Z] * d));
        }
        else
            return (Vector3<T>(T(0), -m_quarterHeight, T(0)));
    }
}

// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ void Cone<T>::readConvex(std::istream& fileIn)
{
    fileIn >> m_bottomRadius >> m_quarterHeight;
    m_quarterHeight /= T(4);
    m_sinAngle
        = m_bottomRadius
          / sqrt(m_bottomRadius * m_bottomRadius + T(16) * m_quarterHeight * m_quarterHeight);
}

// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ void Cone<T>::writeConvex(std::ostream& fileOut) const
{
    fileOut << "Cone: " << m_bottomRadius << ", " << T(4) * m_quarterHeight << ".\n";
}

// -------------------------------------------------------------------------------------------------
// Returns the number of points to write the cone in a Paraview format
template <typename T>
__HOST__ int Cone<T>::numberOfPoints_PARAVIEW() const
{
    return (visuNodeNbOnPer + 2);
}

// -------------------------------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the cone in a
// Paraview format
template <typename T>
__HOST__ int Cone<T>::numberOfCells_PARAVIEW() const
{
    return (visuNodeNbOnPer);
}

// -------------------------------------------------------------------------------------------------
// Returns a list of points describing the cone in a Paraview format
template <typename T>
__HOST__ std::list<Vector3<T>> Cone<T>::writePoints_PARAVIEW(const Transform3<T>& transform,
                                                             Vector3<T> const*    translation) const
{
    list<Vector3<T>> ParaviewPoints;
    Vector3<T>       pp, p;
    T                dtheta = TWO_PI<T> / visuNodeNbOnPer;

    // Disk rim
    p[Y] = -m_quarterHeight;
    for(int i = 0; i < visuNodeNbOnPer; ++i)
    {
        p[X] = m_bottomRadius * cos(i * dtheta);
        p[Z] = m_bottomRadius * sin(i * dtheta);
        pp   = transform(p);
        if(translation)
            pp += *translation;
        ParaviewPoints.push_back(pp);
    }

    // Disk center
    p[X] = T(0);
    p[Y] = -m_quarterHeight;
    p[Z] = T(0);
    pp   = transform(p);
    if(translation)
        pp += *translation;
    ParaviewPoints.push_back(pp);

    // Upper tip
    p[X] = T(0);
    p[Y] = T(3) * m_quarterHeight;
    p[Z] = T(0);
    pp   = transform(p);
    if(translation)
        pp += *translation;
    ParaviewPoints.push_back(pp);

    return (ParaviewPoints);
}

// -------------------------------------------------------------------------------------------------
// Writes the connectivity of the cone in a Paraview format
template <typename T>
__HOST__ void Cone<T>::writeConnection_PARAVIEW(std::list<uint>& connectivity,
                                                std::list<uint>& offsets,
                                                std::list<uint>& cellstype,
                                                uint&            firstpoint_globalnumber,
                                                uint&            last_offset) const
{
    for(uint i = 0; i < visuNodeNbOnPer - 1; ++i)
    {
        connectivity.push_back(firstpoint_globalnumber + i);
        connectivity.push_back(firstpoint_globalnumber + i + 1);
        connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer);
        connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer + 1);
        last_offset += 4;
        offsets.push_back(last_offset);
        cellstype.push_back(10);
    }
    connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer - 1);
    connectivity.push_back(firstpoint_globalnumber);
    connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer);
    connectivity.push_back(firstpoint_globalnumber + visuNodeNbOnPer + 1);
    last_offset += 4;
    offsets.push_back(last_offset);
    cellstype.push_back(10);

    firstpoint_globalnumber += visuNodeNbOnPer + 2;
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Cone<float>;
template class Cone<double>;

#undef visuNodeNbOnPer