#ifndef _OBSTACLEBUILDERFACTORY_HH_
#define _OBSTACLEBUILDERFACTORY_HH_

#include "ReaderXML.hh"
#include <iostream>
using namespace std;

class Obstacle;


/** @brief The class ObstacleBuilderFactory.

    Creates the appropriate obstacle depending on options.

    @author Institut Francais du Petrole - 2003 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ObstacleBuilderFactory
{
  public:
    /**@name Static methods */
    //@{
    /** @brief Creates an obstacle from an XML node
    @param root XML node */
    static Obstacle* create( DOMNode* root );

    /** @brief Creates and reloads an obstacle from a stream
    @param tag obstacle type
    @param mother higher level obstacle 
    @param file input stream */
    static void reload( string const& tag, Obstacle& mother, 
	istream& file );
    //@}


  private:
    /**@name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    ObstacleBuilderFactory() {}

    /** @brief Destructor (forbidden) */
    ~ObstacleBuilderFactory() {}
    //@}
};

#endif
