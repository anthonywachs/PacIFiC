#ifndef _STLTRIANGLE_HH_
#define _STLTRIANGLE_HH_


#include "STLVertex.hh"
#include "Vector3.hh"
#include <tuple>
using namespace std;
using namespace solid;


/** @brief The class STLTriangle.

    Definition of a triangle in a STL triangulation.

    @author A.MORENTE - 2023 - Creation */
// ============================================================================
class STLTriangle
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    STLTriangle();
    
    /** @brief Contructor with parameters */
    STLTriangle( tuple<STLVertex*,STLVertex*,STLVertex*> ve, 
    	Vector3 const& ne, size_t const& ide );    

    /** @brief Destructor */
    ~STLTriangle();
    //@}


    /** @name Get methods */
    //@{  
    /** @brief Returns the surface area of the triangle */
    double getSurfaceArea() const;
    //@}
    
    
    /** @name Core methods */
    //@{  
    /** @brief Computes the surface area of the triangle */
    void computeSurfaceArea();
    //@}    


    /**@name Class Friend */
    //@{
    friend class STLObstacle; /**< STL obstacle the STLTriangle is part of */   
    //@}

  
  private:        
    /** @name Parameters */
    //@{
    tuple<STLVertex*,STLVertex*,STLVertex*> m_v; /**< pointers to the three
    	vertices of the triangle */
    Vector3 m_n; /**< normal vector */  
    size_t m_id; /**< ID number */
    double m_surfacearea; /**< surface area */
    //@}    
};

#endif
