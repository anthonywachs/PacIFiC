#include "GrainsParameters.hh"

// -------------------------------------------------------------------------------------------------
// Static variables
/* GPU */
template <typename T>
bool GrainsParameters<T>::m_isGPU = false;
template <typename T>
cudaDeviceProp GrainsParameters<T>::m_GPU = {};

/* Dynamic Settings */
template <typename T>
SimulationState<T> GrainsParameters<T>::m_simulationState;

/* Spatial */
template <typename T>
Vector3<T> GrainsParameters<T>::m_origin = Vector3<T>(0, 0, 0);
template <typename T>
Vector3<T> GrainsParameters<T>::m_maxCoordinate = Vector3<T>(0, 0, 0);
template <typename T>
bool GrainsParameters<T>::m_isPeriodic = false;

/* Collision Detection */
template <typename T>
CollisionDetectionParameters<T> GrainsParameters<T>::m_collisionDetection;

/* Material */
template <typename T>
std::unordered_map<std::string, uint> GrainsParameters<T>::m_materialMap;
template <typename T>
uint GrainsParameters<T>::m_numContactPairs = 0;
template <typename T>
bool GrainsParameters<T>::m_isContactWithMemory = false;
template <typename T>
bool GrainsParameters<T>::m_useCompaction = true;

/* Temporal */
template <typename T>
T GrainsParameters<T>::m_tStart;
template <typename T>
T GrainsParameters<T>::m_tEnd;
template <typename T>
T GrainsParameters<T>::m_dt;
template <typename T>
bool GrainsParameters<T>::m_isLeapFrog = false;

/* Physical */
template <typename T>
Vector3<T> GrainsParameters<T>::m_gravity = Vector3<T>(0, 0, 0);

/* Post-Processing */
template <typename T>
std::queue<T> GrainsParameters<T>::m_tSave;
template <typename T>
uint GrainsParameters<T>::m_verbosityFrequency = 1;
template <typename T>
GrainsSimTimer GrainsParameters<T>::m_simTimer;
template <typename T>
GrainsCDMTimer GrainsParameters<T>::m_cdmTimer;
template <typename T>
GrainsForceTimer GrainsParameters<T>::m_fmTimer;

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class GrainsParameters<float>;
template class GrainsParameters<double>;