#ifndef _KINEMATICS_HH_
#define _KINEMATICS_HH_

/** @brief The class Kinematics.

    Manages the kinematics of mobile rigid bodies.

    @author G.FERRER - Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Kinematics
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Destructor */
    virtual ~Kinematics() {};
    //@}

  protected:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    Kinematics() {};
    //@}
};

#endif
