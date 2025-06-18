#ifndef _CONVEXBUILDERFACTORY_HH_
#define _CONVEXBUILDERFACTORY_HH_

class Convex;

#include <iostream>
#include <string>
using namespace std;

#include "ReaderXML.hh"


/** @brief The class ConvexBuilderFactory.

    Static class that constructs a convex using input data from an XML node or a
    stream. Create a vonvex of the derived type and returns a pointer to the
    mother class Convex. Constructors, destructor and equal operator are 
    forbidden.

    @author Institut Francais du Petrole - 2003 - Creation
    @author D.RAKOTONIRINA - IFP Energies nouvelles - Oct. 2014
    - Modification 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ConvexBuilderFactory
{
  public:
    /** @name Methods Static */
    //@{
    /** @brief Construct a convex with an XML node as an input parameter
    @param root XML node */
    static Convex* create( DOMNode *root );

    /** @brief Construct a convex with a type and a input stream as input 
    parameters
    @param type convex type
    @param fileIn input stream */
    static Convex* create( string &type, istream& fileIn );
    //@}


  private:
    /** @name Constructors */
    //@{
    /** @brief Constructor */
    ConvexBuilderFactory();
    
    /** @brief Copy constructor */
    ConvexBuilderFactory( ConvexBuilderFactory const& cb );    

    /** @brief Destructor */
    ~ConvexBuilderFactory();
    
    /** @brief Equal operator to another ConvexBuilderFactory object
    @param cb the other ConvexBuilderFactory object */
    ConvexBuilderFactory& operator = ( ConvexBuilderFactory const& cb );
    //@}
};

#endif
