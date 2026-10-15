#ifndef _GRAINSGPU_HH_
#define _GRAINSGPU_HH_

#include "Grains.hh"
#include "ReaderXML.hh"

// =================================================================================================
/** @brief The class GrainsGPU.

    Standard Grains3D application on GPU.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class GrainsGPU : public Grains<T>
{
protected:
    /** \brief Memory buffer for rigid bodies on the device. */
    GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE> m_d_rigidBodyList;
    /** \brief GPU manager of the components in the simulation. */
    std::unique_ptr<ComponentManager<T, MemType::DEVICE>> m_d_components;
    /** \brief Buffer of contact forces. */
    GrainsMemBuffer<ContactForceModel<T>*, MemType::DEVICE> m_d_contactForce;
    /** \brief Buffer of time integrators. */
    GrainsMemBuffer<TimeIntegrator<T>*, MemType::DEVICE> m_d_timeIntegrator;
    //@}

public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor */

    GrainsGPU();

    /** @brief Destructor */
    ~GrainsGPU();
    //@}

    /** @name High-level methods */
    //@{
    /** @brief Sets up the GPU and its parameters */
    void setupGPUDevice();

    /** @brief Tasks to perform before time-stepping.
        @param rootElement XML root */
    void initialize(DOMElement* rootElement) final;

    /** @brief Runs the simulation over the prescribed time interval */
    void simulate() final;

    // /** @brief Tasks to perform after time-stepping */
    // virtual void finalize();
    //@}

    /**@name Low-level methods */
    //@{
    /** @brief Construction of the simulation: linked cell, particles & obstacles, domain
        decomposition.
        @param rootElement XML root */
    void Construction(DOMElement* rootElement);

    /** @brief External force definition
        @param rootElement XML root */
    void Forces(DOMElement* rootElement);

    /** @brief Additional features of the simulation: insertion, post-processing.
        @param rootElement XML root */
    void AdditionalFeatures(DOMElement* rootElement);
    //@}
};

#endif
