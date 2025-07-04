#ifndef _APPBROWNIAN_HH_
#define _APPBROWNIAN_HH_

#include "App.hh"
#include "ReaderXML.hh"
#include <vector>
using namespace std;


/** @brief The class AppBrownian.

    Computes a random Brownian force.

    @author A.WACHS - 2024 - Creation */
// ============================================================================
class AppBrownian : public App
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Constructor with an XML node as an input parameter 
    @param root XML node 
    @param rank process rank 
    @param error returned value, 0 is OK, non zero is error */
    AppBrownian( DOMNode* root, int rank, size_t& error );

    /** @brief Destructor */
    ~AppBrownian();
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes AppBrownian forces exerted on rigid bodies
    @param time physical time
    @param dt time step magnitude
    @param particles active particles */
    void ComputeForces( double time, double dt,
  	list<Particle*> const* particles );
    //@}


  private:
    /** @name Parameters */
    //@{
    double m_mag; /**< maximum magnitude of the AppBrownian force */
    //@}

    
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    AppBrownian();
    //@}    
};

#endif
