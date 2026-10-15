#ifndef _TORCE_HH_
#define _TORCE_HH_

#include "Vector3.hh"

// =================================================================================================
/** @brief The class Torce.

    A torque and a force (Torque + Force = Torce) imposed on the center of mass
    of a component.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class Torce
{
private:
    /**@name Parameter */
    //@{
    Vector3<T> m_torque; /**< torques exerted at the center of mass */
    Vector3<T> m_force;  /**< force exerted on the component */
    //@}

public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    __HOSTDEVICE__
    Torce();

    /** @brief Constructor with a torque and a force as input parameters
        @param t torque
        @param f force */
    __HOSTDEVICE__
    Torce(const Vector3<T>& t, const Vector3<T>& f);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~Torce();
    //@}

    /**@name Get methods */
    //@{
    /** @brief Gets the total torque of the torce */
    __HOSTDEVICE__
    Vector3<T> getTorque() const;

    /** @brief Gets the total force of the torce */
    __HOSTDEVICE__
    Vector3<T> getForce() const;
    //@}

    /**@name Set methods */
    //@{
    /** @brief Sets the torque of the torce
        @param t torque */
    __HOSTDEVICE__
    void setTorque(const Vector3<T>& t);

    /** @brief Sets the force of the torce
        @param f force */
    __HOSTDEVICE__
    void setForce(const Vector3<T>& f);
    //@}

    /**@name Methods */
    //@{
    /** @brief Resets the torce */
    __HOSTDEVICE__
    void reset();

    /** @brief Adds a torque to the torce
        @param t added torque */
    __HOSTDEVICE__
    void addTorque(const Vector3<T>& t);

    /** @brief Adds a force to the torce
        @param f added force */
    __HOSTDEVICE__
    void addForce(const Vector3<T>& f);

    /** @brief Adds a force to the torce with accounting for the additional torque.
        @param f added force
        @param p point of application */
    __HOSTDEVICE__
    void addForce(const Vector3<T>& f, const Vector3<T>& p);

    /** @brief Atomically adds a torque to the torce (GPU safe)
        @param t added torque */
    __device__ void addTorqueAtomic(const Vector3<T>& t);

    /** @brief Atomically adds a force to the torce (GPU safe)
        @param f added force */
    __device__ void addForceAtomic(const Vector3<T>& f);

    /** @brief Atomically adds both force and torque from another Torce (GPU safe)
        @param other Torce containing force and torque to add */
    __device__ void addTorceAtomic(const Torce<T>& other);
    //@}
};

/** @name External Methods - I/O methods */
//@{
/** @brief Input operator
@param fileIn input stream
@param t torce */
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Torce<T>& t);

/** @brief Output operator
@param fileOut output stream
@param t torce */
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Torce<T>& t);
//@}

typedef Torce<float>  TorceF;
typedef Torce<double> TorceD;

#endif
