#ifndef _CONTACTFORCEMODEL_HH_
#define _CONTACTFORCEMODEL_HH_

#include "ContactInfo.hh"
#include "ContactTable.hh"
#include "Torce.hh"
#include "Vector3.hh"

// ContactForceModel types
enum ContactForceModelType
{
    HOOKE,
    HOOKEMEMORY
};

// =================================================================================================
/** @brief The class ContactForceModel.

    Defines the contact forces between two colliding components and computes
    these contact forces.

    @author A.YAZDANI - 2024 - Construction */
// =================================================================================================
template <typename T>
class ContactForceModel
{
protected:
    /**@name Contructors */
    //@{
    /** @brief Default constructor (forbidden except in derived classes) */
    __HOSTDEVICE__
    ContactForceModel();

    /** @brief Copy constructor
        @param cf ContactForceModel object to be copied */
    __HOSTDEVICE__
    ContactForceModel(ContactForceModel<T> const& cf);
    //@}

public:
    /**@name Contructors */
    //@{
    /** @brief Destructor */
    __HOSTDEVICE__
    virtual ~ContactForceModel();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Returns the ContactForceModel type */
    __HOSTDEVICE__
    virtual ContactForceModelType getContactForceModelType() const = 0;
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a torce based on the contact information
        @param contactInfos geometric contact features
        @param relVelocityAtContact relative velocity at the contact point
        @param relAngVelocity relative angular velocity
        @param vA position of the first component
        @param vB position of the second component
        @param contactHistory pointer to contact history
        @param torceA computed force and torque for the first component
        @param torceB computed force and torque for the second component */
    __HOSTDEVICE__
    virtual void computeForces(const ContactInfo<T>& contactInfos,
                               const Vector3<T>&     relVelocityAtContact,
                               const Vector3<T>&     relAngVelocity,
                               const Vector3<T>&     vA,
                               const Vector3<T>&     vB,
                               ContactHistory<T>*    contactHistory,
                               Torce<T>&             torceA,
                               Torce<T>&             torceB) const
        = 0;
    //@}
};

#endif
