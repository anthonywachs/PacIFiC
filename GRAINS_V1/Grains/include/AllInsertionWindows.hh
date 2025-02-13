#ifndef _ALLINSERTIONWINDOWS_HH_
#define _ALLINSERTIONWINDOWS_HH_

#include "Basic.hh"
#include "Error.hh"
#include <vector>
using namespace std;
#include "Window.hh"


/** @brief The class AllInsertionWindows.

    Manages all particle insertion windows.

    @author A.WACHS - 2025 - Creation */
// ============================================================================
class AllInsertionWindows
{
  public:
    /** @name Constructors & Destructor */
    //@{
    /** @brief Default constructor */
    AllInsertionWindows();

    /** @brief Destructor */
    ~AllInsertionWindows();
    //@}


    /** @name Methods */
    //@{
    /** @brief Adds a new insertion window
    @param iwindow new insertion window */
    void addWindow( Window const& iwindow );
    
    /** @brief Converts the vector of insertion windows to a list of rigid
    bodies */
    void convertToRigidBodies( list<RigidBody*>& iwlist ) const;
    
    /** @brief Returns a pointer to a randomly selected window */
    Window const* getRandomWindow() const;
    
    /** @brief Associates the imposed velocity to an insertion window
    @param impvel imposed velocity */
    void LinkImposedMotion( ObstacleImposedVelocity* impvel );
    
    /** @brief Moves the insertion windows if any has an imposed motion
    @param time physical time 
    @param dt time step magnitude */
    void Move( double time, double dt );                    
    //@}


  private:
    /** @name Parameters */
    //@{
    vector<Window> m_insertion_windows; /**< Insertion windows */
    list<ObstacleImposedVelocity*> m_AllImposedVelocitiesOnWindows; /**< list
  	of all imposed velocities on insertion windows */      
    //@}
};

#endif
