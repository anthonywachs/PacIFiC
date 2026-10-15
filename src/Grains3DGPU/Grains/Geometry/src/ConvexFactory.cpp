#include "ConvexFactory.hh"
#include "Box.hh"
#include "Cone.hh"
#include "Convex.hh"
#include "Cylinder.hh"
#include "GrainsUtils.hh"
#include "Rectangle.hh"
#include "Sphere.hh"
#include "Superquadric.hh"

// -------------------------------------------------------------------------------------------------
// Construct a convex with an XML node as an input parameter
template <typename T>
Convex<T>* ConvexFactory<T>::create(DOMNode* root)
{
    Convex<T>* convex;
    DOMNode*   element = ReaderXML::getNodeNext(root);
    string     type    = ReaderXML::getNodeName(element);

    if(type == "Sphere")
        convex = new Sphere<T>(element);
    else if(type == "Box")
        convex = new Box<T>(element);
    else if(type == "Cylinder")
        convex = new Cylinder<T>(element);
    else if(type == "Cone")
        convex = new Cone<T>(element);
    else if(type == "Superquadric")
        convex = new Superquadric<T>(element);
    else if(type == "Rectangle")
        convex = new Rectangle<T>(element);
    else
        GAbort("Invalid convex type:", type.c_str(), "Aborting Grains!");

    return (convex);  // returns the host-side pointer to convex object
}

// -------------------------------------------------------------------------------------------------
// Construct a convex with a type and a input stream as input parameters
template <typename T>
Convex<T>* ConvexFactory<T>::create(string& type, istream& fileIn)
{
    Convex<T>* convex = NULL;
    if(type == "Sphere")
        convex = new Sphere<T>(fileIn);
    else if(type == "Box")
        convex = new Box<T>(fileIn);
    else if(type == "Cylinder")
        convex = new Cylinder<T>(fileIn);
    else if(type == "Cone")
        convex = new Cone<T>(fileIn);
    else if(type == "Superquadric")
        convex = new Superquadric<T>(fileIn);
    else if(type == "Rectangle")
        convex = new Rectangle<T>(fileIn);
    else
        GAbort("Invalid convex type:", type.c_str(), "Aborting Grains!");

    return (convex);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class ConvexFactory<float>;
template class ConvexFactory<double>;