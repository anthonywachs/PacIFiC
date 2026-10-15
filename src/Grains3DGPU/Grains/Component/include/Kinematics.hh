#ifndef _KINEMATICS_HH_
#define _KINEMATICS_HH_

#include "Vector3.hh"

// =================================================================================================
/** @brief The class Kinematics.

    Manages the kinematics (velocities or accelerations) of components.

    @author A.YAZDANI - 2024 - Construction */
// =================================================================================================
template <typename T>
class Kinematics
{
protected:
    /** @name Parameters */
    //@{
    Vector3<T> m_translational; /**< Translational component */
    Vector3<T> m_angular;       /**< Angular component */
    //@}

public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    __HOSTDEVICE__
    Kinematics();

    /** @brief Constructor with two vectors as input parameters
        @param translational translational components
        @param angular angular component */
    __HOSTDEVICE__
    Kinematics(const Vector3<T>& translational, const Vector3<T>& angular);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Kinematics();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the translational component of the kinematics */
    __HOSTDEVICE__
    Vector3<T> getTranslationalComponent() const;

    /** @brief Gets the angular component of the kinematics */
    __HOSTDEVICE__
    Vector3<T> getAngularComponent() const;
    //@}

    /** @name Set methods */
    //@{
    /** @brief Sets the translational component of the kinematics
        @param translational translational component */
    __HOSTDEVICE__
    void setTranslationalComponent(const Vector3<T>& translational);

    /** @brief Sets the angular component of the kinematics
        @param angular angular component */
    __HOSTDEVICE__
    void setAngularComponent(const Vector3<T>& angular);
    //@}

    /** @name Methods */
    //@{
    /** @brief Adds a vector to the translational component
        @param translational vector to add to the translational component */
    __HOSTDEVICE__
    void addToTranslationalComponent(const Vector3<T>& translational);

    /** @brief Adds a angular velocity to the angular velocity
        @param angular vector to add to the angular component */
    __HOSTDEVICE__
    void addToAngularComponent(const Vector3<T>& angular);

    /** @brief Returns the total velocity U + om x R given R
        @param R arm vector */
    __HOSTDEVICE__
    Vector3<T> kinematicsAtPoint(const Vector3<T>& R) const;
    //@}
};

/** @name External Methods - I/O methods */
//@{
/** @brief Input operator
    @param fileIn input stream
    @param k kinematics */
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Kinematics<T>& k);

/** @brief Output operator
    @param fileOut output stream
    @param k kinematics */
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Kinematics<T>& k);
//@}

#endif
