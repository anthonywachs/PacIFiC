#include "TimeIntegratorFactory.hh"
#include "FirstOrderExplicit.hh"
#include "SecondOrderLeapFrog.hh"

/* ============================================================================================== */
/* Low-Level Methods                                                                              */
/* ============================================================================================== */
// GPU kernel to construct the timeIntegrator on device.
// This is mandatory as we cannot access device memory addresses on the host
// So, we pass a device memory address to a kernel.
// Memory address is then populated within the kernel.
template <typename T>
__GLOBAL__ void createTimeIntegratorKernel(TimeIntegrator<T>**      TI,
                                           const uint               index,
                                           const TimeIntegratorType timeIntegratorType,
                                           const T                  dt)
{
    uint tID = blockIdx.x * blockDim.x + threadIdx.x;
    if(tID > 0)
        return;

    if(timeIntegratorType == FIRSTORDEREXPLICIT)
        TI[index] = new FirstOrderExplicit<T>(dt);
    else if(timeIntegratorType == SECONDORDERLEAPFROG)
        TI[index] = new SecondOrderLeapFrog<T>(dt);
    else
        GAbort("Time integrator is not implemented for GPU!", "Aborting Grains!");
}

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
// Creates and stores a TimeIntegrator object in the host memory.
template <typename T>
__HOST__ void TimeIntegratorFactory<T>::create(
    DOMNode* root, T dt, GrainsMemBuffer<TimeIntegrator<T>*, MemType::HOST>& TI)
{
    TI.initialize(1);  // Initialize memory for one time integrator

    std::string type = ReaderXML::getNodeAttr_String(root, "Type");
    if(type == "FirstOrderExplicit")
        TI[0] = new FirstOrderExplicit<T>(dt);
    else if(type == "SecondOrderLeapFrog")
        TI[0] = new SecondOrderLeapFrog<T>(dt);
    else
        GAbort("Unknown time integration! Aborting Grains!");
}

// -------------------------------------------------------------------------------------------------
// Constructs a TimeIntegrator object on device.
template <typename T>
__HOST__ void TimeIntegratorFactory<T>::copyHostToDevice(
    GrainsMemBuffer<TimeIntegrator<T>*, MemType::HOST>&   h_TI,
    GrainsMemBuffer<TimeIntegrator<T>*, MemType::DEVICE>& d_TI)
{
    // Allocate the device memory for the time integrator
    d_TI.initialize(h_TI.getSize());
    for(uint i = 0; i < h_TI.getSize(); ++i)
    {
        if(h_TI[i] == nullptr)
            continue;

        // Extracting info from the host side object
        const TimeIntegratorType timeIntegratorType = h_TI[i]->getTimeIntegratorType();
        const T                  dt                 = h_TI[i]->getTimeStep();
        createTimeIntegratorKernel<<<1, 1>>>(d_TI.getData(), i, timeIntegratorType, dt);
    }
    cudaDeviceSynchronize();
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class TimeIntegratorFactory<float>;
template class TimeIntegratorFactory<double>;
