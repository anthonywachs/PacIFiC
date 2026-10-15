#ifndef _CONVEXFACTORY_HH_
#define _CONVEXFACTORY_HH_

#include "Convex.hh"
#include "ReaderXML.hh"

// =================================================================================================
/** @brief The class ConvexFactory.

    Static class that constructs a convex using input data from an XML node or a stream. Create a
    convex of the derived type and returns a pointer to the mother class Convex.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation
    @author D. RAKOTONIRINA - IFP Energies nouvelles - Oct. 2014 - Modification
    @author A.WACHS - 2019 - Major cleaning & refactoring
    @author A.Yazdani - 2024 - Modification for GPUs */
// =================================================================================================
template <typename T>
class ConvexFactory
{
private:
    /** @name Constructors */
    //@{
    /** @brief Constructor */
    ConvexFactory();

    /** @brief Copy constructor.
    @param cb the other ConvexFactory object */
    ConvexFactory(ConvexFactory<T> const& cb);

    /** @brief Equal operator to another ConvexFactory object.
        @param cb the other ConvexFactory object */
    ConvexFactory<T>& operator=(ConvexFactory<T> const& cb);

    /** @brief Destructor */
    ~ConvexFactory();
    //@}

public:
    /** @name Methods */
    //@{
    /** @brief Construct a convex with an XML node as an input parameter.
        @param root XML node */
    static Convex<T>* create(DOMNode* root);

    /** @brief Construct a convex with a type and a input stream as input parameters.
        @param type convex type
        @param fileIn input stream */
    static Convex<T>* create(string& type, std::istream& fileIn);
    //@}
};

#endif
