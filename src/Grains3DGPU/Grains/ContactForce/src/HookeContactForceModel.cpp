#include "HookeContactForceModel.hh"
#include "GrainsUtils.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ HookeContactForceModel<T>::HookeContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node
template <typename T>
__HOST__ HookeContactForceModel<T>::HookeContactForceModel(DOMNode* root)
{
    GAssert(ReaderXML::hasNodeAttr(root, "kn"), "kn not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "en"), "en not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "etat"), "etat not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "muc"), "muc not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "kr"), "kr not defined! Aborting Grains!");

    m_kn   = T(ReaderXML::getNodeAttr_Double(root, "kn"));
    m_en   = T(ReaderXML::getNodeAttr_Double(root, "en"));
    m_etat = T(ReaderXML::getNodeAttr_Double(root, "etat"));
    m_muc  = T(ReaderXML::getNodeAttr_Double(root, "muc"));
    m_kr   = T(ReaderXML::getNodeAttr_Double(root, "kr"));

    m_muen = log(m_en) / sqrt(PI<T> * PI<T> + log(m_en) * log(m_en));
}

// -------------------------------------------------------------------------------------------------
// Constructor with five values as contact parameters
template <typename T>
__HOSTDEVICE__ HookeContactForceModel<T>::HookeContactForceModel(T kn, T en, T etat, T muc, T kr)
    : m_kn(kn)
    , m_en(en)
    , m_etat(etat)
    , m_muc(muc)
    , m_kr(kr)
{
    m_muen = log(m_en) / sqrt(PI<T> * PI<T> + log(m_en) * log(m_en));
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ HookeContactForceModel<T>::~HookeContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Gets the ContactForceModel type
template <typename T>
__HOSTDEVICE__ ContactForceModelType HookeContactForceModel<T>::getContactForceModelType() const
{
    return (HOOKE);
}

// -------------------------------------------------------------------------------------------------
// Gets the parameters of the Hooke contact force model
template <typename T>
__HOSTDEVICE__ void HookeContactForceModel<T>::getContactForceModelParameters(
    T& kn, T& en, T& etat, T& muc, T& kr) const
{
    kn   = m_kn;
    en   = m_en;
    etat = m_etat;
    muc  = m_muc;
    kr   = m_kr;
}

// -------------------------------------------------------------------------------------------------
// Performs forces & torques computation
template <typename T>
__HOSTDEVICE__ void
    HookeContactForceModel<T>::performForcesCalculus(const Vector3<T>& contactVector,
                                                     const Vector3<T>& relVelocityAtContact,
                                                     const Vector3<T>& relAngVelocity,
                                                     const T           overlapDistance,
                                                     const T           averageMass,
                                                     Vector3<T>&       delFN,
                                                     Vector3<T>&       delFT,
                                                     Vector3<T>&       delM) const
{
    // Notes:
    // - contactVector is a unit vector pointing from B to A
    // - overlapDistance is negative when there is penetration

    // Normal linear elastic force
    // We do this here as we want to modify the penetration vector later
    delFN = m_kn * overlapDistance * contactVector;

    // Unit normal vector at contact point
    Vector3<T> v_n = (relVelocityAtContact * contactVector) * contactVector;
    Vector3<T> v_t = relVelocityAtContact - v_n;

    // Unit tangential vector along relative velocity at contact point
    T          normv_t = norm(v_t);
    Vector3<T> tangent(0, 0, 0);
    if(normv_t > EPS<T>)
        tangent = v_t / normv_t;

    // Normal dissipative force
    T gamman = -T(2) * m_muen * sqrt(averageMass * m_kn);
    delFN -= gamman * v_n;
    T normFN = norm(delFN);

    // Tangential dissipative force
    // If m_etat = -1, we compute its value such that gamma_n = gamma_t, i.e., same damping in the
    // normal and tangential directions
    T etat   = (m_etat == T(-1)) ? (-m_muen * sqrt(m_kn / averageMass)) : m_etat;
    delFT    = (-T(2) * etat * averageMass) * v_t;
    T normFT = norm(delFT);

    // Tangential Coulomb saturation
    T fn = m_muc * normFN;
    if(fn < normFT)
        delFT = (-fn) * tangent;

    // Rolling resistance moment
    delM = Vector3<T>(0, 0, 0);
    if(m_kr)
    {
        // Relative angular velocity at contact point
        Vector3<T> wn     = (relAngVelocity * contactVector) * contactVector;
        Vector3<T> wt     = relAngVelocity - wn;
        T          normwt = norm(wt);

        // Anti-spinning effect along the normal wn
        delM = -m_kr * normFN * T(0.001) * wn;

        // Classical rolling resistance moment
        if(normwt > EPS<T>)
            delM -= m_kr * normFN * wt;
    }
}

// -------------------------------------------------------------------------------------------------
// Returns a torce based on the contact information
template <typename T>
__HOSTDEVICE__ void HookeContactForceModel<T>::computeForces(const ContactInfo<T>& contactInfos,
                                                             const Vector3<T>& relVelocityAtContact,
                                                             const Vector3<T>& relAngVelocity,
                                                             const Vector3<T>& vA,
                                                             const Vector3<T>& vB,
                                                             ContactHistory<T>* contactHistory,
                                                             Torce<T>&          torceA,
                                                             Torce<T>&          torceB) const
{
    // Note: contactHistory is unused for non-memory models
    (void)contactHistory;

    // Get snapshot with all contact information
    auto snapshot = contactInfos.getSnapshot();

    // Compute contact force and torque
    Vector3<T> delFN, delFT, delM;
    performForcesCalculus(snapshot.contactVector,
                          relVelocityAtContact,
                          relAngVelocity,
                          snapshot.overlapDistance,
                          snapshot.averageMass,
                          delFN,
                          delFT,
                          delM);

    delFN += delFT;
    torceA.addForce(delFN, snapshot.contactPoint - vA);
    torceB.addForce(-delFN, snapshot.contactPoint - vB);
    if(m_kr)
    {
        torceA.addTorque(delM);
        torceB.addTorque(-delM);
    }
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class HookeContactForceModel<float>;
template class HookeContactForceModel<double>;