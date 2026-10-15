#include "HookeMemoryContactForceModel.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "QuaternionMath.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ HookeMemoryContactForceModel<T>::HookeMemoryContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node
template <typename T>
__HOST__ HookeMemoryContactForceModel<T>::HookeMemoryContactForceModel(DOMNode* root)
{
    GAssert(ReaderXML::hasNodeAttr(root, "kn"), "kn not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "en"), "en not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "kt"), "kt not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "etat"), "etat not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "muc"), "muc not defined! Aborting Grains!");
    GAssert(GrainsParameters<T>::m_dt > T(0), "Time step not defined! Aborting Grains!");

    m_kn   = T(ReaderXML::getNodeAttr_Double(root, "kn"));
    m_en   = T(ReaderXML::getNodeAttr_Double(root, "en"));
    m_kt   = T(ReaderXML::getNodeAttr_Double(root, "kt"));
    m_etat = T(ReaderXML::getNodeAttr_Double(root, "etat"));
    m_muc  = T(ReaderXML::getNodeAttr_Double(root, "muc"));
    if(ReaderXML::hasNodeAttr(root, "mur"))
        m_mur = T(ReaderXML::getNodeAttr_Double(root, "mur"));
    if(ReaderXML::hasNodeAttr(root, "etarpf"))
        m_etarpf = T(ReaderXML::getNodeAttr_Double(root, "etarpf"));

    m_muen = log(m_en) / sqrt(PI<T> * PI<T> + log(m_en) * log(m_en));
    m_dt   = GrainsParameters<T>::m_dt;
}

// -------------------------------------------------------------------------------------------------
// Constructor with eight values as contact parameters
template <typename T>
__HOSTDEVICE__ HookeMemoryContactForceModel<T>::HookeMemoryContactForceModel(
    T kn, T en, T kt, T etat, T muc, T mur, T etarpf, T dt)
    : m_kn(kn)
    , m_en(en)
    , m_kt(kt)
    , m_etat(etat)
    , m_muc(muc)
    , m_mur(mur)
    , m_etarpf(etarpf)
    , m_dt(dt)
{
    m_muen = log(m_en) / sqrt(PI<T> * PI<T> + log(m_en) * log(m_en));
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ HookeMemoryContactForceModel<T>::~HookeMemoryContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Gets the ContactForceModel type
template <typename T>
__HOSTDEVICE__ ContactForceModelType
    HookeMemoryContactForceModel<T>::getContactForceModelType() const
{
    return (HOOKEMEMORY);
}

// -------------------------------------------------------------------------------------------------
// Gets the parameters of the HookeMemory contact force model
template <typename T>
__HOSTDEVICE__ void HookeMemoryContactForceModel<T>::getContactForceModelParameters(
    T& kn, T& en, T& kt, T& etat, T& muc, T& mur, T& etarpf) const
{
    kn     = m_kn;
    en     = m_en;
    kt     = m_kt;
    etat   = m_etat;
    muc    = m_muc;
    mur    = m_mur;
    etarpf = m_etarpf;
}

// -------------------------------------------------------------------------------------------------
// Performs forces & torques computation with optional memory tracking
template <typename T>
__HOSTDEVICE__ void
    HookeMemoryContactForceModel<T>::performForcesCalculus(const Vector3<T>&  contactVector,
                                                           const Vector3<T>&  relVelocityAtContact,
                                                           const Vector3<T>&  relAngVelocity,
                                                           const T            overlapDistance,
                                                           const T            averageMass,
                                                           ContactHistory<T>* contactHistory,
                                                           Vector3<T>&        delFN,
                                                           Vector3<T>&        delFT,
                                                           Vector3<T>&        delM) const
{
    // Notes:
    // - contactVector is a unit vector pointing from B to A (the normal)
    // - overlapDistance is negative when there is penetration

    // Penetration vector
    Vector3<T> normal = contactVector;

    // Relative velocity components
    Vector3<T> v_n = (relVelocityAtContact * normal) * normal;
    Vector3<T> v_t = relVelocityAtContact - v_n;

    // 1) Compute normal force
    // Normal linear elastic force
    delFN = m_kn * overlapDistance * normal;

    // Normal dissipative force
    T gamman = -T(2) * m_muen * sqrt(averageMass * m_kn);
    delFN -= gamman * v_n;
    T normFN = norm(delFN);

    // 2) Compute tangential force with memory
    // Check if this is a new contact (previousNormal is zero)
    bool contactExisted = (norm(contactHistory->m_previousNormal) > EPS<T>);

    // Rotate previous cumulative displacement to current contact plane
    if(contactExisted)
    {
        Quaternion<T> qrot;
        qrot.setRotFromTwoVectors(contactHistory->m_previousNormal, normal);
        contactHistory->m_tangentialDisplacement = qrot >> contactHistory->m_tangentialDisplacement;
    }

    // Add contribution of current timestep
    contactHistory->m_tangentialDisplacement += m_dt * m_kt * v_t;

    // Update the normal vector in history
    contactHistory->m_previousNormal = normal;

    // Compute tangential force direction
    // If m_etat = -1, we compute its value such that gamma_n = gamma_t, i.e., same damping in the
    // normal and tangential directions
    T          etat        = (m_etat == T(-1)) ? (-m_muen * sqrt(m_kn / averageMass)) : m_etat;
    Vector3<T> viscousFT   = (-T(2) * etat * averageMass) * v_t;
    Vector3<T> tentativeFT = -contactHistory->m_tangentialDisplacement + viscousFT;
    T          normFT      = norm(tentativeFT);
    Vector3<T> tangentDir = (normFT > EPS<T>) ? tentativeFT / normFT : Vector3<T>(T(0), T(0), T(0));

    // Compute tangential force with Coulomb limit
    if(normFT <= m_muc * normFN)
    {
        // Below Coulomb limit
        delFT = normFT * tangentDir;
    }
    else
    {
        // Above Coulomb limit - apply saturation and adjust history
        delFT = m_muc * normFN * tangentDir;
        if(m_kt > EPS<T>)
        {
            contactHistory->m_tangentialDisplacement
                = (-m_muc * normFN * tangentDir + viscousFT) / m_kt;
        }
    }

    // 3) Compute rolling resistance torque with memory (if applicable)
    delM = Vector3<T>(T(0), T(0), T(0));
    if(m_mur > EPS<T>)
    {
        // Using Jiang et al (2005, 2015) formulation
        // Note: Req would require particle radii - using a simplified approach here
        // For now, use a characteristic length scale based on penetration
        T Req  = T(1);  // Placeholder - should be computed from particle radii if available
        T kr   = T(3) * m_kn * m_mur * m_mur * Req * Req;       // torque spring stiffness
        T etar = T(3) * (-gamman) * m_mur * m_mur * Req * Req;  // torque dissipative coefficient
        T maxNormMk = m_mur * Req * normFN;                     // saturation torque

        // Rotate previous rolling friction spring to current plane
        if(contactExisted)
        {
            Quaternion<T> qrot;
            qrot.setRotFromTwoVectors(contactHistory->m_previousNormal, normal);
            contactHistory->m_rollingDisplacement = qrot >> contactHistory->m_rollingDisplacement;
        }

        // Update rolling friction spring
        contactHistory->m_rollingDisplacement -= kr * m_dt * relAngVelocity;
        T normMk = norm(contactHistory->m_rollingDisplacement);

        // Apply saturation
        if(normMk > maxNormMk)
        {
            contactHistory->m_rollingDisplacement *= maxNormMk / normMk;
        }

        delM = contactHistory->m_rollingDisplacement - m_etarpf * etar * relAngVelocity;
    }
}

// -------------------------------------------------------------------------------------------------
// Returns a force based on the contact information with memory tracking
template <typename T>
__HOSTDEVICE__ void
    HookeMemoryContactForceModel<T>::computeForces(const ContactInfo<T>& contactInfos,
                                                   const Vector3<T>&     relVelocityAtContact,
                                                   const Vector3<T>&     relAngVelocity,
                                                   const Vector3<T>&     vA,
                                                   const Vector3<T>&     vB,
                                                   ContactHistory<T>*    contactHistory,
                                                   Torce<T>&             torceA,
                                                   Torce<T>&             torceB) const
{
    // Get snapshot with all contact information
    auto snapshot = contactInfos.getSnapshot();

    // Compute contact forces and torques
    Vector3<T> delFN, delFT, delM;
    performForcesCalculus(snapshot.contactVector,
                          relVelocityAtContact,
                          relAngVelocity,
                          snapshot.overlapDistance,
                          snapshot.averageMass,
                          contactHistory,
                          delFN,
                          delFT,
                          delM);

    // Apply forces and torques
    Vector3<T> totalForce = delFN + delFT;
    torceA.addForce(totalForce, snapshot.contactPoint - vA);
    torceB.addForce(-totalForce, snapshot.contactPoint - vB);

    if(m_mur > EPS<T>)
    {
        torceA.addTorque(delM);
        torceB.addTorque(-delM);
    }
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class HookeMemoryContactForceModel<float>;
template class HookeMemoryContactForceModel<double>;
