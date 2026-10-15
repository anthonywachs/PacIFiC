#include "InsertionWindow.hh"
#include "GrainsUtils.hh"
#include "VectorMath.hh"
#include <ctime>

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOST__ InsertionWindow<T>::InsertionWindow()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with XML node and the type of the seed
template <typename T>
__HOST__ InsertionWindow<T>::InsertionWindow(DOMNode* dn, RandomGeneratorSeed seed)
{
    GoutWI(12, "Reading insertion window ...");

    // Setting up the random generator engine with the given seed
    if(seed == RGS_DEFAULT)
        m_randGenerator.seed(0);
    else if(seed == RGS_RANDOM)
        m_randGenerator.seed(time(NULL));

    std::string nType = ReaderXML::getNodeAttr_String(dn, "Type");
    if(nType == "Box")
    {
        m_type         = BOXWINDOW;
        DOMNode* nP1   = ReaderXML::getNode(dn, "MinPoint");
        T        xVal1 = T(ReaderXML::getNodeAttr_Double(nP1, "X"));
        T        yVal1 = T(ReaderXML::getNodeAttr_Double(nP1, "Y"));
        T        zVal1 = T(ReaderXML::getNodeAttr_Double(nP1, "Z"));
        m_v1           = Vector3<T>(xVal1, yVal1, zVal1);
        DOMNode* nP2   = ReaderXML::getNode(dn, "MaxPoint");
        T        xVal2 = T(ReaderXML::getNodeAttr_Double(nP2, "X"));
        T        yVal2 = T(ReaderXML::getNodeAttr_Double(nP2, "Y"));
        T        zVal2 = T(ReaderXML::getNodeAttr_Double(nP2, "Z"));
        m_v2           = Vector3<T>(xVal2, yVal2, zVal2);
        GoutWI(15,
               "Box insertion window with min and max points",
               Vector3ToString(m_v1),
               "and",
               Vector3ToString(m_v2) + ".");
        GoutWI(12, "Reading insertion window completed!");
    }
    else if(nType == "Annulus")
    {
        m_type         = ANNULUSWINDOW;
        DOMNode* nP1   = ReaderXML::getNode(dn, "BottomPoint");
        T        xVal1 = T(ReaderXML::getNodeAttr_Double(nP1, "X"));
        T        yVal1 = T(ReaderXML::getNodeAttr_Double(nP1, "Y"));
        T        zVal1 = T(ReaderXML::getNodeAttr_Double(nP1, "Z"));
        m_v1           = Vector3<T>(xVal1, yVal1, zVal1);
        DOMNode* nP2   = ReaderXML::getNode(dn, "TopPoint");
        T        xVal2 = T(ReaderXML::getNodeAttr_Double(nP2, "X"));
        T        yVal2 = T(ReaderXML::getNodeAttr_Double(nP2, "Y"));
        T        zVal2 = T(ReaderXML::getNodeAttr_Double(nP2, "Z"));
        m_v2           = Vector3<T>(xVal2, yVal2, zVal2);
        DOMNode* nR    = ReaderXML::getNode(dn, "Radius");
        m_iRad         = T(ReaderXML::getNodeAttr_Double(nR, "Inner"));
        m_oRad         = T(ReaderXML::getNodeAttr_Double(nR, "Outter"));
        GoutWI(15,
               "Annulus insertion window with bottom point",
               Vector3ToString(m_v1),
               ", direction",
               Vector3ToString(m_v2),
               ", outter radius",
               m_oRad,
               ", and inner radius",
               m_iRad,
               ".");
        GoutWI(12, "Reading insertion window completed!");
    }
    else
        GAbort("Insertion window type is not supported");
}

// -------------------------------------------------------------------------------------------------
// Constructor for box window with min and max points
template <typename T>
__HOST__ InsertionWindow<T>::InsertionWindow(const Vector3<T>& minPoint,
                                             const Vector3<T>& maxPoint,
                                             unsigned          seed)
    : m_v1(minPoint)
    , m_v2(maxPoint)
    , m_iRad(T(0))
    , m_oRad(T(0))
    , m_dist(T(0), T(1))
    , m_type(BOXWINDOW)
{
    m_randGenerator.seed(seed);
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOST__ InsertionWindow<T>::~InsertionWindow()
{
}

// -------------------------------------------------------------------------------------------------
// Generates a random number with uniform distribution in window
template <typename T>
__HOST__ Vector3<T> InsertionWindow<T>::generateRandomPoint()
{
    Vector3<T> out;
    if(m_type == BOXWINDOW)
    {
        out = Vector3<T>(m_v1[X] + m_dist(m_randGenerator) * (m_v2[X] - m_v1[X]),
                         m_v1[Y] + m_dist(m_randGenerator) * (m_v2[Y] - m_v1[Y]),
                         m_v1[Z] + m_dist(m_randGenerator) * (m_v2[Z] - m_v1[Z]));
    }
    else if(m_type == ANNULUSWINDOW)
    {
        // Axis from bottom to top and its length
        Vector3<T> axis = m_v2 - m_v1;
        T          len  = sqrt(axis * axis);
        GAssert(len > HIGHEPS<T>, "Annulus insertion window has zero height!");
        Vector3<T> k = (T(1) / len) * axis;  // unit axis

        // Build an orthonormal basis (u, v) spanning the disk plane
        Vector3<T> ref
            = (fabs(k[Z]) < T(0.999)) ? Vector3<T>(T(0), T(0), T(1)) : Vector3<T>(T(0), T(1), T(0));
        Vector3<T> uvec = k ^ ref;  // perpendicular to axis
        T          un   = sqrt(uvec * uvec);
        if(un < HIGHEPS<T>)
        {
            // Fallback if ref nearly parallel to axis
            ref  = Vector3<T>(T(1), T(0), T(0));
            uvec = k ^ ref;
            un   = sqrt(uvec * uvec);
            GAssert(un > HIGHEPS<T>, "Failed to construct orthonormal basis for annulus window!");
        }
        uvec            = (T(1) / un) * uvec;
        Vector3<T> vvec = k ^ uvec;

        // Sample random angle theta in [0, 2*pi)
        T theta = m_dist(m_randGenerator) * TWO_PI<T>;

        // Sample random radius r with area-uniform distribution in [iRad, oRad]
        T ri = m_iRad;
        T ro = m_oRad;
        if(ro < ri)
            std::swap(ri, ro);
        GAssert(ro >= T(0), "Annulus outer radius must be non-negative!");
        T u = m_dist(m_randGenerator);
        T r = sqrt((T(1) - u) * ri * ri + u * ro * ro);

        // Sample random height along the axis segment [m_v1, m_v2]
        T          t = m_dist(m_randGenerator);
        Vector3<T> h = m_v1 + t * axis;

        // Assemble point in the local annulus plane perpendicular to axis
        out = h + r * (cos(theta) * uvec + sin(theta) * vvec);
    }
    return (out);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class InsertionWindow<float>;
template class InsertionWindow<double>;