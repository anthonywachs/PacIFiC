#include "ContactForceModel.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ ContactForceModel<T>::ContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Copy constructor
template <typename T>
__HOSTDEVICE__ ContactForceModel<T>::ContactForceModel(ContactForceModel<T> const& cf)
{
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ ContactForceModel<T>::~ContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class ContactForceModel<float>;
template class ContactForceModel<double>;