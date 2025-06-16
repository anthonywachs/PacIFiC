#ifndef _RoughWall_HH_
#define _RoughWall_HH_

#include "CompositeObstacle.hh"
#include "SimpleObstacle.hh"

#include <list>
#include <map>
using namespace std;

#include "ReaderXML.hh"


/** @brief The class RoughWall

    Cylindrical shell made of boxes.

    @author A.WACHS - 2024 - Creation */
// ============================================================================
class RoughWall : public CompositeObstacle
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    RoughWall( DOMNode* root );

    /** @brief Constructor with name as input parameter
    @param s obstacle name */
    RoughWall( string const& s );

    /** @brief Destructor */
    ~RoughWall();
    //@}


    /** @name Methods */
    //@{    
    /** @brief Checks if there is anything special to do about periodicity and
    if there is applies periodicity 
    @param LC linked-cell grid */
    void periodicity( LinkedCell* LC );         
    //@}


    /** @name I/O methods */
    //@{
    /** @brief Reloads the rough wall and links it to the higher level 
    obstacle in the obstacle tree
    @param mother higher level obstacle
    @param file input stream */
    virtual void reload( Obstacle& mother, istream& file ) ; 
    
    /** @brief Outputs the rough wall for reload
    @param fileSave output stream */
    virtual void write( ostream& fileSave ) const;       
    //@}
    

  private:
    /** @name Parameters */
    //@{  
    double m_shift; /**< position shift wrt z=0  */
    size_t* m_nb_spheres; /**< number of spheres in each direction */	
    double m_radius; /**< sphere radius */
    double m_random_mag; /**< position randomness magnitude */
    double m_periodic_ext; /**< roughness periodic extension length */
    bool m_periodic; /**< whether to handle periodicity */
    list<size_t> m_periodic_directions; /**< list of periodic directions, cannot
    	contain more than 2 directions */
    bool m_restrict_box_geommotion; /**< the box does not move geometrically in
    	the periodic directions if set to true */
    Point3 m_domain_origin; /**< global domain origin */
    Vector3 m_domain_size; /**< global domain size */
    Vector3 m_lper; /**< actual periodic domain extension per direction */ 
    map<int,Point3> m_pos_nm1;	       
    //@}


    /** @name Constructors */
    //@{
    /** @brief Copy constructor
    @param copy copied RoughWall */
    RoughWall( RoughWall const& copy );   
    //@} 
    
    
    /** @name Methods */
    //@{    
    /** @brief Returns whether a position is in the global domain in the 
    directions specified by m_periodic_directions */
    bool isInDomain( Point3 const* pos );         
    //@}      
};

#endif
