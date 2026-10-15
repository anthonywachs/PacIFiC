#ifndef _RIGIDBODY_HH_
#define _RIGIDBODY_HH_

#include <limits>

#include "Convex.hh"
#include "Kinematics.hh"
#include "Quaternion.hh"
#include "ReaderXML.hh"
#include "Torce.hh"

// =================================================================================================
/** @brief The class RigidBody.

    Rigid bodies comprising their shapes and physical attributes.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class RigidBody
{
protected:
    /** @name Parameters */
    //@{
    /** \brief Convex shape */
    Convex<T>* m_convex;
    /** \brief Crust thickness */
    T m_crustThickness;
    /** \brief Circumscribed radius */
    T m_circumscribedRadius;
    /** \brief Mass of the rigid body */
    T m_mass;
    /** \brief Equivalent radius */
    T m_equivalentRadius;
    /** \brief Material ID */
    uint m_material;
    /** \brief Diagonal inertia tensor (3 components: Ixx, Iyy, Izz) */
    T m_inertia[3];
    //@}

public:
    /** @name Structures */
    //@{
    /** @brief Snapshot of commonly used properties for fast bulk access */
    struct PropertiesSnapshot
    {
        Convex<T>* convex              = nullptr;
        T          crustThickness      = T(0);
        T          circumscribedRadius = T(0);
        T          mass                = T(0);
        T          equivalentRadius    = T(0);
        uint       material            = std::numeric_limits<uint>::max();
    };
    //@}

    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    __HOSTDEVICE__
    RigidBody();

    /** @brief Constructor with a convex, crust thickness, material, and density.
        @param convex convex
        @param ct crust thickness of the rigid body
        @param density density
        @param material material ID */
    __HOSTDEVICE__
    RigidBody(Convex<T>* convex, T ct, T density, uint material);

    /** @brief Copy constructor
        @param rb RigidBody object to be copied */
    __HOSTDEVICE__
    RigidBody(RigidBody<T> const& rb);

    /** @brief Copy assignment operator
        @param other RigidBody object to be assigned */
    __HOSTDEVICE__
    RigidBody<T>& operator=(const RigidBody<T>& other);

    /** @brief Move constructor
          @param other RigidBody object to be moved */
    __HOSTDEVICE__
    RigidBody(RigidBody<T>&& other);

    /** @brief Move assignment operator
        @param other RigidBody object to be moved */
    __HOSTDEVICE__
    RigidBody<T>& operator=(RigidBody<T>&& other);

    /** @brief Constructor with an XML input
        @param root XML input */
    __HOST__
    RigidBody(DOMNode* root);

    /** @brief Destructor */
    __HOSTDEVICE__
    ~RigidBody();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Gets the rigid body's convex */
    __HOSTDEVICE__
    Convex<T>* getConvex() const;

    /** @brief Gets the rigid body's crust thickness */
    __HOSTDEVICE__
    T getCrustThickness() const;

    /** @brief Gets the circumcribed radius of the rigid body */
    __HOSTDEVICE__
    T getCircumscribedRadius() const;

    /** @brief Gets the rigid body's mass */
    __HOSTDEVICE__
    T getMass() const;

    /** @brief Gets the equivalent radius of the rigid body */
    __HOSTDEVICE__
    T getEquivalentRadius() const;

    /** @brief Gets the rigid body's material ID */
    __HOSTDEVICE__
    uint getMaterial() const;

    /** @brief Gets the rigid body's diagonal inertia
        @param inertia the destination for diagonal inertia (3 components) */
    __HOSTDEVICE__
    void getInertia(T (&inertia)[3]) const;

    /** @brief Gets the rigid body's volume */
    __HOSTDEVICE__
    T getVolume() const;

    /** @brief Gets a snapshot of commonly used properties (excluding inertia) */
    __HOSTDEVICE__
    PropertiesSnapshot getPropertiesSnapshot() const;
    //@}

    /**@name Set methods */
    //@{
    /** @brief Sets the rigid body's inertia and circumscribed radius */
    __HOSTDEVICE__
    void setInertia();

    /** @brief Overrides mass and principal inertia values for a composite master sub-body.
        Unlike setInertia(), this bypasses the shape-derived computation and directly injects
        the composite totals (computed via parallel-axis theorem + diagonalization at init).
        @param mass total composite mass
        @param ixx principal moment Ixx
        @param iyy principal moment Iyy
        @param izz principal moment Izz */
    __HOSTDEVICE__
    void setCompositeProperties(T mass, T ixx, T iyy, T izz);
    //@}

    /**@name Methods */
    //@{
    /** @brief Computes the acceleration of the rigid body as a kinematics object after imposing a
        torce (Torque + Force). The assumption is that the torce is given in the body-fixed
        coordinate system, hence, there is no need to have the quaternion.
        @param omega angular velocity in the body-fixed coordinate system
        @param t imposed torce in the body-fixed coordinate system */
    __HOSTDEVICE__
    Kinematics<T> computeMomentum(const Vector3<T>& omega, const Torce<T>& t) const;

    /** @brief Computes the acceleration of the rigid body as a kinematics object after imposing a
        torce (Torque + Force). The assumption is that the torce is given in the space-fixed
        coordinate system.
        @param omega angular velocity in the space-fixed coordinate system
        @param t imposed torce in the space-fixed coordinate system
        @param q quaternion of rotation from space to body coordinate systems */
    __HOSTDEVICE__
    Kinematics<T>
        computeMomentum(const Vector3<T>& omega, const Torce<T>& t, const Quaternion<T>& q) const;
    //@}
};

#endif
