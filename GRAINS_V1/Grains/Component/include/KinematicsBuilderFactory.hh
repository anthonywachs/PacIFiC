#ifndef _KINEMATICSBUILDERFACTORY_HH_
#define _KINEMATICSBUILDERFACTORY_HH_

#include <string>
#include <iostream>
using namespace std;

class ParticleKinematics;
class Convex;


/** @brief The class KinematicsBuilderFactory.

    Creates the appropriate particle kinematics depending on the space 
    dimension and the rigid body shape.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class KinematicsBuilderFactory
{
  public:
    /**@name Methods Static */
    //@{
    /** @brief Creates and returns the particle kinematics depending on the 
    space dimension and the particle rigid body shape. Note about convention: 
    the pointer to the convex is NULL for a composite particle
    @param convex_ particle rigid body convex type */
    static ParticleKinematics* create( Convex const* convex_ );

    /** @brief Creates and returns the particle kinematics from an input stream
    and the particle rigid body shape. Note about convention: 
    the pointer to the convex is NULL for a composite particle
    @param fileIn input stream
    @param convex_ particle rigid body convex type */
    static ParticleKinematics* read( istream& fileIn, Convex const* convex_ );
    //@}


  private:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor (forbidden) */
    KinematicsBuilderFactory() {}

    /** @brief Destructor (forbidden) */
    ~KinematicsBuilderFactory() {}
    //@}
};

#endif
    
