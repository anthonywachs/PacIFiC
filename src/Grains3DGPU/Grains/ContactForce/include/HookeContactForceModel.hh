#ifndef _HOOKECONTACTFORCEMODEL_HH_
#define _HOOKECONTACTFORCEMODEL_HH_

#include "ContactForceModel.hh"
#include "ContactTable.hh"
#include "ReaderXML.hh"

// =================================================================================================
/** @brief The class HookeContactForceModel.

    Contact force model involving a normal Hookean spring, a normal Dashpot and
    a tangential Coulomb friction (HO-D-C) to compute the force and torque
    induced by the contact between two rigid components.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class HookeContactForceModel : public ContactForceModel<T>
{
private:
    /** @name Parameters */
    //@{
    /** \brief Normal stiffness coefficient */
    T m_kn;
    /** \brief Normal restitution coefficient */
    T m_en;
    /** \brief log(m_en) / sqrt( PI * PI + log(m_en) * log(m_en) ) */
    T m_muen;
    /** \brief Tangential damping coefficient */
    T m_etat;
    /** \brief Tangential Coulomb friction coefficient */
    T m_muc;
    /** \brief Rolling resistance coefficient */
    T m_kr;
    //@}

public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    __HOSTDEVICE__
    HookeContactForceModel();

    /** @brief Constructor with an XML node
        @param root XML node */
    __HOST__
    HookeContactForceModel(DOMNode* root);

    /** @brief Constructor with five values as contact parameters
        @param kn normal stiffness coefficient
        @param en normal restitution coefficient
        @param etat tangential damping coefficient
        @param muc tangential Coulomb friction coefficient
        @param kr rolling resistance coefficient */
    __HOSTDEVICE__
    HookeContactForceModel(T kn, T en, T etat, T muc, T kr);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~HookeContactForceModel();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the ContactForceModel type */
    __HOSTDEVICE__
    ContactForceModelType getContactForceModelType() const final;

    /** @brief Gets the parameters of the Hooke contact force model
        @param kn normal stiffness coefficient
        @param en normal restitution coefficient
        @param etat tangential damping coefficient
        @param muc tangential Coulomb friction coefficient
        @param kr rolling resistance coefficient */
    __HOSTDEVICE__
    void getContactForceModelParameters(T& kn, T& en, T& etat, T& muc, T& kr) const;
    //@}

    /**@name Methods */
    //@{
    /** @brief Returns a torce based on the contact information
        @param contactVector geometric contact vector
        @param relVelocityAtContact relative velocity at the contact point
        @param relAngVelocity relative angular velocity
        @param overlapDistance overlap distance
        @param averageMass average mass of the two components
        @param delFN normal force
        @param delFT tangential force
        @param delM torque */
    __HOSTDEVICE__
    void performForcesCalculus(const Vector3<T>& contactVector,
                               const Vector3<T>& relVelocityAtContact,
                               const Vector3<T>& relAngVelocity,
                               const T           overlapDistance,
                               const T           averageMass,
                               Vector3<T>&       delFN,
                               Vector3<T>&       delFT,
                               Vector3<T>&       delM) const;

    /** @brief Returns a torce based on the contact information
        @param contactInfos geometric contact features
        @param relVelocityAtContact relative velocity at the contact point
        @param relAngVelocity relative angular velocity
        @param vA position of the first component
        @param vB position of the second component
        @param contactHistory pointer to contact history (unused, always nullptr)
        @param torceA computed force and torque for the first component
        @param torceB computed force and torque for the second component */
    __HOSTDEVICE__
    void computeForces(const ContactInfo<T>& contactInfos,
                       const Vector3<T>&     relVelocityAtContact,
                       const Vector3<T>&     relAngVelocity,
                       const Vector3<T>&     vA,
                       const Vector3<T>&     vB,
                       ContactHistory<T>*    contactHistory,
                       Torce<T>&             torceA,
                       Torce<T>&             torceB) const final;
    //@}
};

#endif
