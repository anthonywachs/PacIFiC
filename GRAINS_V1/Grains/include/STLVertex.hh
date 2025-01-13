#ifndef _STLVERTEX_HH_
#define _STLVERTEX_HH_

#include "Point3.hh"
#include "Vector3.hh"
using namespace std;
using namespace solid;


/** @brief The class STLVertex.

    Definition of a vertex in a STL triangulation.

    @author A.MORENTE - 2023 - Creation */
// ============================================================================
class STLVertex
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    STLVertex();
    
    /** @brief Contructor with parameters */
    STLVertex( double x, double y, double z, Vector3 const& ne, 
    	size_t const& ide );    

    /** @brief Destructor */
    ~STLVertex();
    //@}


    /** @name Core methods */
    //@{  	
    //@}


    /**@name Class Friend */
    //@{
    friend class STLObstacle; /**< STL obstacle the STLVertex is part of */
    friend class STLTriangle; /**< STL triangle the STLVertex is part of */    
    //@}

  
  private:        
    /** @name Parameters */
    //@{
    Point3 m_p; /**< vertex coordinates */
    Vector3 m_n; /**< normal vector */  
    size_t m_id; /**< ID number */ 
    //@}    
};

#endif
  
