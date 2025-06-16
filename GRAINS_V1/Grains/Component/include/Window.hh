#ifndef _WINDOW_HH_
#define _WINDOW_HH_

#include "Basic.hh"
#include "Point3.hh"
#include "ObstacleKinematicsVelocity.hh"
using namespace std;

class RigidBody;


/** @brief Insertion window type */
enum WindowType
{
  WINDOW_BOX, /**< Box */
  WINDOW_CYLINDER, /**< Cylinder */
  WINDOW_ANNULUS, /**< Annulus */
  WINDOW_LINE, /**< Line */
  WINDOW_NONE /**< unknown */
};


/** @brief The class Window.

    Geometric subspace of simple shape.

    @author A.WACHS - 2025 - Creation */
// ============================================================================
class Window
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    Window();
    
    /** @brief Copy constructor 
    @param copy copied Component object */
    Window( Window const& copy );    

    /** @brief Destructor */
    ~Window();
    //@}


    /** @name Methods */
    //@{
    /** @brief Reads a window from an XML node
    @param nWindow XML node
    @param oshift empty string to shift the output 
    @param rank process rank */
    bool readWindow( DOMNode* nWindow, string const& oshift,
    	int const& rank );
    
    /** @brief Returns a point randomly selected within the window */
    Point3 getInsertionPoint() const;
    
    /** @brief Sets the window as a box
    @param ptA lower-left-behind corner 
    @param ptB upper-right-front corner */
    void setAsBox( Point3 const& ptA, Point3 const& ptB );
    
    /** @brief Sets the window as a cylinder
    @param bottomC bottom centre of the cylinder 
    @param radius_ cylinder radius
    @param height_ cylinder height
    @param axisdir_str axis direction as a string */
    void setAsCylinder( Point3 const& bottomC, double const& radius_,
    	 double const& height_, string const& axisdir_str );    
    
    /** @brief Adds the window as a rigid body to a list of rigid bodies
    @param iwlist list of rigid bodies */
    void addAsRigidBody( list<RigidBody*>& iwlist ) const; 
    
    /** @brief Shifts the 2 points defined the window by a specified amount in a
    specified direction
    @param geoshift magnitude of the shift
    @param dir direction of the shift (x, y or z) */
    void shiftWindow( double const& geoshift, Direction dir ); 
    
    /** @brief Associates the imposed velocity to the window
    @param impvel imposed velocity */
    bool LinkImposedMotion( ObstacleImposedVelocity* impvel ); 
    
    /** @brief Moves the window if it has an imposed motion
    @param time physical time 
    @param dt time step magnitude */
    void Move( double time, double dt );
    
    /** @brief Resets kinematics to 0 */
    void resetKinematics();                  
    //@}


    /** @name Accessors */
    //@{
    /** @brief Returns a pointer to point A */
    Point3 const* getPointA() const;
    
    /** @brief Returns a pointer to point B */
    Point3 const* getPointB() const;
    
    /** @brief Returns the cylinder radius */
    double getRadius() const; 
    
    /** @brief Returns the cylinder height */
    double getHeight() const;
    
    /** @brief Returns the window type */
    WindowType getType() const;
    
    /** @brief Returns the cylinder axis direction */
    Direction getAxisDirection() const;                    
    //@}
    
    
  private:
    /** @name Parameters */
    //@{
    WindowType m_ftype; /**< Window type */
    Point3 m_ptA; /**< Box 1st corner or center of lower disk of the cylinder */
    Point3 m_ptB; /**< Box 2nd corner */
    double m_radius; /**< Cylinder radius */
    double m_radius_int; /**< Inner cylinder radius in case of annulus */
    double m_height; /**< Cylinder height */
    Direction m_axisdir; /**< Cylinder axis direction */
    string m_name; /**< Window name */
    ObstacleKinematicsVelocity* m_kinematics; /**< obstacle kinematics with
    	imposed velocity */    
    //@}   
};

#endif
