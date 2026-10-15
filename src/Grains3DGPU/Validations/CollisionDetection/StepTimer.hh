#ifndef _STEPTIMER_HH_
#define _STEPTIMER_HH_

#include <chrono>

// =================================================================================================
/** @brief Timer utility for measuring execution time */
class StepTimer
{
private:
    std::chrono::high_resolution_clock::time_point m_start;
    std::chrono::high_resolution_clock::time_point m_end;
    bool                                           m_running;

public:
    StepTimer()
        : m_running(false)
    {
    }

    void start()
    {
        m_start   = std::chrono::high_resolution_clock::now();
        m_running = true;
    }

    void stop()
    {
        m_end     = std::chrono::high_resolution_clock::now();
        m_running = false;
    }

    double getElapsedMilliseconds() const
    {
        if(m_running)
            return 0.0;
        return std::chrono::duration<double, std::milli>(m_end - m_start).count();
    }

    double getElapsedSeconds() const
    {
        return getElapsedMilliseconds() / 1000.0;
    }
};

#endif