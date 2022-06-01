#ifndef _TIMEINTEGRATORBUILDERFACTORY_HH_
#define _TIMEINTEGRATORBUILDERFACTORY_HH_

#include <string>
using std::string;
class TimeIntegrator;


/** @brief The class TimeIntegratorBuilderFactory.

    Creates the numerical scheme for the time integration of the Newton's law 
    and the kinematic equations. 

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation 
    @author A.WACHS - 2021 - Major cleaning & refactoring */
// ============================================================================
class TimeIntegratorBuilderFactory
{
  public:
    /**@name Methods Static */
    //@{
    /** @brief Creates and returns the time integration scheme
    @return L'integrateur en time */
    static TimeIntegrator* create();
    //@}


  private:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor (forbidden) */
    TimeIntegratorBuilderFactory() {}

    /** @brief Destructor (forbidden) */
    ~TimeIntegratorBuilderFactory() {}
    //@}
};

#endif
