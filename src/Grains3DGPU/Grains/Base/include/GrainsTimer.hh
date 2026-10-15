#ifndef _GRAINSTIMER_HH_
#define _GRAINSTIMER_HH_

#include <chrono>
#include <iomanip>
#include <iostream>

/** @name Enumerations */
//@{
/** @brief Timed stages inside the main Grains simulation loop.
    Values are contiguous from 0 and used directly as array indices. */
enum class SimStage : uint8_t
{
    Total = 0,             ///< Full simulation wall-clock span
    Insertion,             ///< insertParticles() [+ copyTo(device) for GPU]
    DetectCollisions,      ///< ComponentManager::detectCollisions()
    ComputeContactForces,  ///< ComponentManager::computeContactForces()
    MoveParticles,         ///< ComponentManager::moveParticles()
    AdvanceVelocity,       ///< ComponentManager::advanceVelocity() (KDK only)
    PostProcess,           ///< Grains::postProcess()
    COUNT
};

/** @brief Timed sub-stages inside the CollisionDetectionModule pipeline.
    Values are contiguous from 0 and used directly as array indices. */
enum class CDMStage : uint8_t
{
    Sort = 0,           ///< CollisionDetectionModule::sortParticles()
    NeighborList,       ///< CollisionDetectionModule::updateNeighborList()
    RelativeTransform,  ///< CollisionDetectionModule::computeRelativeTransformations()
    BVFilter,           ///< CollisionDetectionModule::filterPairsBV()
    NarrowPhase,        ///< CollisionDetectionModule::detectCollisionsComponents*()
    Transform,          ///< CollisionDetectionModule::transformContactInfo()
    Total,              ///< Full CollisionDetectionModule::run() wall-clock span
    COUNT
};

/** @brief Timed sub-stages inside the ForceModule (contact force computation) pipeline.
        On the HOST path FlagAndCompact and ReduceTorces contain no work and accumulate ~0 s. */
enum class FMStage : uint8_t
{
    FlagAndCompact = 0,  ///< [GPU] flagActivePairs_Kernel + CUB compaction; no-op on CPU
    ComputeForces,       ///< computeContactForces_Kernel / sequential per-pair loop
    ReduceTorces,    ///< [GPU] reduceTorces_Kernel atomic per-particle accumulation; no-op on CPU
    ExternalForces,  ///< addExternalForces gravity application
    AssembleComposites,  ///< assembleCompositeTorces sub-body -> master accumulation
    Total,               ///< Full ForceModule::run() wall-clock span
    COUNT
};
//@}

// =================================================================================================
/** @brief Wall-clock timer for the main simulation-loop stages.

        Call enable(isGPU) once at setup -- the timer records whether it is running on the GPU and
        issues cudaDeviceSynchronize() inside stop() only when needed.  When disabled,
        start()/stop() are single always-false predicted branches with no clock calls.

        Usage:
        @code
            // in Grains.cpp after reading XML:
            GP::m_simTimer.enable(GP::m_isGPU);

            // in simulate():
            auto& t = GP::m_simTimer;
            t.start(SimStage::Total);
            t.start(SimStage::Insertion);
            components->insertParticles(...);
            t.stop(SimStage::Insertion);   // syncs device automatically if running on GPU
            ...
            t.stop(SimStage::Total);
            if(t.isEnabled()) t.printSummary();
        @endcode

        @author A.Yazdani - 2026 - Construction */
// =================================================================================================
class GrainsSimTimer
{
    using Clock                     = std::chrono::high_resolution_clock;
    static constexpr uint8_t kCount = static_cast<uint8_t>(SimStage::COUNT);

    bool              m_enabled = false;
    bool              m_isGPU   = false;
    double            m_accum[kCount]{};
    Clock::time_point m_t0[kCount]{};

public:
    /** @brief Enables the timer.  @p isGPU controls whether stop() syncs the device. */
    void enable(bool isGPU = false) noexcept
    {
        m_enabled = true;
        m_isGPU   = isGPU;
    }
    void disable() noexcept
    {
        m_enabled = false;
    }
    bool isEnabled() const noexcept
    {
        return m_enabled;
    }
    void reset() noexcept
    {
        for(uint8_t i = 0; i < kCount; ++i)
            m_accum[i] = 0.0;
    }

    void start(SimStage s) noexcept
    {
        if(m_enabled)
            m_t0[static_cast<uint8_t>(s)] = Clock::now();
    }

    void stop(SimStage s) noexcept
    {
        if(m_enabled)
        {
            if(m_isGPU)
                cudaDeviceSynchronize();
            const uint8_t i = static_cast<uint8_t>(s);
            m_accum[i] += std::chrono::duration<double>(Clock::now() - m_t0[i]).count();
        }
    }

    double operator[](SimStage s) const noexcept
    {
        return m_accum[static_cast<uint8_t>(s)];
    }

    void printSummary() const
    {
        const double tot = (*this)[SimStage::Total] > 0.0 ? (*this)[SimStage::Total] : 1e-12;
        auto         row = [&](const char* label, double t) {
            std::cout << "  " << std::left << std::setw(30) << label << std::right << std::setw(10)
                      << std::fixed << std::setprecision(4) << t << " s   (" << std::setw(6)
                      << std::setprecision(2) << 100.0 * t / tot << " %)\n";
        };
        const double staged = (*this)[SimStage::Insertion] + (*this)[SimStage::DetectCollisions]
                              + (*this)[SimStage::ComputeContactForces]
                              + (*this)[SimStage::MoveParticles]
                              + (*this)[SimStage::AdvanceVelocity] + (*this)[SimStage::PostProcess];
        const double other = (tot > staged) ? (tot - staged) : 0.0;

        std::cout << '\n' << std::string(80, '=') << '\n';
        std::cout << "  Simulation Loop Performance Summary\n";
        std::cout << std::string(80, '-') << '\n';
        row("Insertion:", (*this)[SimStage::Insertion]);
        row("Collision Detection:", (*this)[SimStage::DetectCollisions]);
        row("Contact Forces:", (*this)[SimStage::ComputeContactForces]);
        row("Move Particles:", (*this)[SimStage::MoveParticles]);
        row("Advance Velocity:", (*this)[SimStage::AdvanceVelocity]);
        row("Post-Processing:", (*this)[SimStage::PostProcess]);
        row("Other:", other);
        std::cout << std::string(80, '-') << '\n';
        row("Total:", tot);
        std::cout << std::string(80, '=') << "\n\n";
    }
};

// =================================================================================================
/** @brief Wall-clock timer for CollisionDetectionModule sub-stage timing.

    Call enable(isGPU) once at setup.  stop() issues cudaDeviceSynchronize() when the timer was
    enabled for GPU execution.  When disabled all calls are branch-only no-ops.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
class GrainsCDMTimer
{
    using Clock                     = std::chrono::high_resolution_clock;
    static constexpr uint8_t kCount = static_cast<uint8_t>(CDMStage::COUNT);

    bool              m_enabled = false;
    bool              m_isGPU   = false;
    double            m_accum[kCount]{};
    Clock::time_point m_t0[kCount]{};

public:
    void enable(bool isGPU = false) noexcept
    {
        m_enabled = true;
        m_isGPU   = isGPU;
    }
    void disable() noexcept
    {
        m_enabled = false;
    }
    bool isEnabled() const noexcept
    {
        return m_enabled;
    }
    void reset() noexcept
    {
        for(uint8_t i = 0; i < kCount; ++i)
            m_accum[i] = 0.0;
    }

    void start(CDMStage s) noexcept
    {
        if(m_enabled)
            m_t0[static_cast<uint8_t>(s)] = Clock::now();
    }

    void stop(CDMStage s) noexcept
    {
        if(m_enabled)
        {
            if(m_isGPU)
                cudaDeviceSynchronize();
            const uint8_t i = static_cast<uint8_t>(s);
            m_accum[i] += std::chrono::duration<double>(Clock::now() - m_t0[i]).count();
        }
    }

    double operator[](CDMStage s) const noexcept
    {
        return m_accum[static_cast<uint8_t>(s)];
    }

    void printSummary() const
    {
        const double total       = (*this)[CDMStage::Total];
        const double subStageSum = (*this)[CDMStage::Sort] + (*this)[CDMStage::NeighborList]
                                   + (*this)[CDMStage::RelativeTransform]
                                   + (*this)[CDMStage::BVFilter] + (*this)[CDMStage::NarrowPhase]
                                   + (*this)[CDMStage::Transform];
        const double denom = total > 0.0 ? total : (subStageSum > 0.0 ? subStageSum : 1e-12);
        const double other = (total > subStageSum) ? (total - subStageSum) : 0.0;
        auto         row   = [&](const char* label, double t) {
            std::cout << "  " << std::left << std::setw(28) << label << std::right << std::setw(10)
                      << std::fixed << std::setprecision(4) << t << " s   (" << std::setw(6)
                      << std::setprecision(2) << 100.0 * t / denom << " %)\n";
        };

        std::cout << '\n' << std::string(80, '=') << '\n';
        std::cout << "  CDM Sub-stage Breakdown\n";
        std::cout << std::string(80, '-') << '\n';
        row("Sort:", (*this)[CDMStage::Sort]);
        row("Neighbor List:", (*this)[CDMStage::NeighborList]);
        row("Relative Transforms:", (*this)[CDMStage::RelativeTransform]);
        row("BV Filter:", (*this)[CDMStage::BVFilter]);
        row("GJK Narrow Phase:", (*this)[CDMStage::NarrowPhase]);
        row("Contact Transform:", (*this)[CDMStage::Transform]);
        row("Other:", other);
        std::cout << std::string(80, '-') << '\n';
        row("Total:", total > 0.0 ? total : subStageSum);
        std::cout << std::string(80, '=') << "\n";
    }
};

// =================================================================================================
/** @brief Wall-clock timer for ForceModule sub-stage timing.

    Call enable(isGPU) once at setup.  stop() issues cudaDeviceSynchronize() when the timer was
    enabled for GPU execution.  When disabled all calls are branch-only no-ops.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
class GrainsForceTimer
{
    using Clock                     = std::chrono::high_resolution_clock;
    static constexpr uint8_t kCount = static_cast<uint8_t>(FMStage::COUNT);

    bool              m_enabled = false;
    bool              m_isGPU   = false;
    double            m_accum[kCount]{};
    Clock::time_point m_t0[kCount]{};

public:
    void enable(bool isGPU = false) noexcept
    {
        m_enabled = true;
        m_isGPU   = isGPU;
    }
    void disable() noexcept
    {
        m_enabled = false;
    }
    bool isEnabled() const noexcept
    {
        return m_enabled;
    }
    void reset() noexcept
    {
        for(uint8_t i = 0; i < kCount; ++i)
            m_accum[i] = 0.0;
    }

    void start(FMStage s) noexcept
    {
        if(m_enabled)
            m_t0[static_cast<uint8_t>(s)] = Clock::now();
    }

    void stop(FMStage s) noexcept
    {
        if(m_enabled)
        {
            if(m_isGPU)
                cudaDeviceSynchronize();
            const uint8_t i = static_cast<uint8_t>(s);
            m_accum[i] += std::chrono::duration<double>(Clock::now() - m_t0[i]).count();
        }
    }

    double operator[](FMStage s) const noexcept
    {
        return m_accum[static_cast<uint8_t>(s)];
    }

    void printSummary() const
    {
        const double total = (*this)[FMStage::Total];
        const double subStageSum
            = (*this)[FMStage::FlagAndCompact] + (*this)[FMStage::ComputeForces]
              + (*this)[FMStage::ReduceTorces] + (*this)[FMStage::ExternalForces]
              + (*this)[FMStage::AssembleComposites];
        const double denom = total > 0.0 ? total : (subStageSum > 0.0 ? subStageSum : 1e-12);
        const double other = (total > subStageSum) ? (total - subStageSum) : 0.0;
        auto         row   = [&](const char* label, double t) {
            std::cout << "  " << std::left << std::setw(28) << label << std::right << std::setw(10)
                      << std::fixed << std::setprecision(4) << t << " s   (" << std::setw(6)
                      << std::setprecision(2) << 100.0 * t / denom << " %)\n";
        };

        std::cout << '\n' << std::string(80, '=') << '\n';
        std::cout << "  FM Sub-stage Breakdown\n";
        std::cout << std::string(80, '-') << '\n';
        row("Flag & Compact:", (*this)[FMStage::FlagAndCompact]);
        row("Compute Forces:", (*this)[FMStage::ComputeForces]);
        row("Reduce Torces:", (*this)[FMStage::ReduceTorces]);
        row("External Forces:", (*this)[FMStage::ExternalForces]);
        row("Assemble Composites:", (*this)[FMStage::AssembleComposites]);
        row("Other:", other);
        std::cout << std::string(80, '-') << '\n';
        row("Total:", total > 0.0 ? total : subStageSum);
        std::cout << std::string(80, '=') << "\n";
    }
};

#endif
