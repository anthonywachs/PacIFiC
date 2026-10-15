#ifndef _CONTACTINFO_HH_
#define _CONTACTINFO_HH_

#include <limits>
#include <type_traits>

#include "Vector3.hh"

// =================================================================================================
/** @brief The class ContactInfo.

    Contains all the features of a contact point.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class ContactInfo
{
protected:
    /** @name Parameters */
    //@{
    /** \brief Contact point */
    Vector3<T> m_contactPoint;
    /** \brief Contact vector */
    Vector3<T> m_contactVector;
    /** \brief Overlap distance */
    T m_overlapDistance;
    /** \brief Average mass */
    T m_averageMass;
    /** \brief Average radius */
    T m_averageRadius;
    /** \brief Contact hash */
    uint m_contactHash;
    //@}

public:
    /** @name Structures */
    //@{
    /** @brief Snapshot of commonly used properties for fast bulk access */
    struct Snapshot
    {
        Vector3<T> contactPoint    = Vector3<T>(T(0), T(0), T(0));
        Vector3<T> contactVector   = Vector3<T>(T(0), T(0), T(0));
        T          overlapDistance = std::numeric_limits<T>::max();
        T          averageMass     = T(0);
        T          averageRadius   = T(0);
        uint       contactHash     = std::numeric_limits<uint>::max();
    };
    //@}

    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    __HOSTDEVICE__
    ContactInfo();

    /** @brief Constructor with contact point location in the world reference frame, overlap vector,
        overlap distance and number of iterations of GJK as input parameters.
        @param pt contact point
        @param vec contact vector
        @param distance_ overlap distance */
    __HOSTDEVICE__
    ContactInfo(const Vector3<T>& pt, const Vector3<T>& vec, T overlap);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~ContactInfo();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the contact point */
    __HOSTDEVICE__
    Vector3<T> getContactPoint() const;

    /** @brief Gets the contact vector */
    __HOSTDEVICE__
    Vector3<T> getContactVector() const;

    /** @brief Gets the overlap distance */
    __HOSTDEVICE__
    T getOverlapDistance() const;

    /** @brief Gets the average mass */
    __HOSTDEVICE__
    T getAverageMass() const;

    /** @brief Gets the average radius */
    __HOSTDEVICE__
    T getAverageRadius() const;

    /** @brief Gets the contact hash */
    __HOSTDEVICE__
    uint getContactHash() const;

    /** @brief Gets a POD snapshot containing all contact fields */
    __HOSTDEVICE__
    Snapshot getSnapshot() const;
    //@}

    /** @name Set methods */
    //@{
    /** @brief Sets the contact point
        @param p contact point */
    __HOSTDEVICE__
    void setContactPoint(const Vector3<T>& p);

    /** @brief Sets the contact vector
        @param v overlap vector */
    __HOSTDEVICE__
    void setContactVector(const Vector3<T>& v);

    /** @brief Sets the overlap distance
        @param d overlap distance  */
    __HOSTDEVICE__
    void setOverlapDistance(T d);

    /** @brief Sets the average mass
        @param avgMass average mass of contacting particles */
    __HOSTDEVICE__
    void setAverageMass(T avgMass);

    /** @brief Sets the average radius
        @param avgRadius average radius of contacting particles */
    __HOSTDEVICE__
    void setAverageRadius(T avgRadius);

    /** @brief Sets the contact hash from two material IDs */
    __HOSTDEVICE__
    void setContactHash(uint hash);

    /** @brief Sets all contact fields from a POD snapshot */
    __HOSTDEVICE__
    void setSnapshot(const Snapshot& s);
    //@}

    /** @name Operators */
    //@{
    /** @brief Equality operator
        @param other the ContactInfo to compare with
        @return true if all members are equal */
    __HOSTDEVICE__
    bool operator==(const ContactInfo<T>& other) const;

    /** @brief Inequality operator
        @param other the ContactInfo to compare with
        @return true if any member is different */
    __HOSTDEVICE__
    bool operator!=(const ContactInfo<T>& other) const;
    //@}
};

/** @name External Methods - I/O methods */
//@{
/** @brief Output operator
    @param fileIn input stream
    @param c contact point object */
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const ContactInfo<T>& c);

/** @brief Input operator
    @param fileIn input stream
    @param c contact point object */
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, ContactInfo<T>& c);
//@}

#endif