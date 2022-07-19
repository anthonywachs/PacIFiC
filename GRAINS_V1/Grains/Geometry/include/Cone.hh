#ifndef _CONE_HH_
#define _CONE_HH_

#include "Convex.hh"


/** @brief The class Cone.

    Convex with a conical shape. From GJK Engine - A Fast and 
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Cone : public Convex {
  public:
    /**@name Constructors */
    //@{
    /** @brief Constructor with flat base radius and height as input parameters
    @param r flat base radius
    @param h height */
    Cone( double r = 0, double h = 0 ); 
  
    /** @brief Constructor with an input stream
    @param fileIn input stream */
    Cone( istream& fileIn );
  
    /** @brief Destructor */
    ~Cone();
    //@}
  

    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Cone support function, returns the support point P, i.e. the
    point on the surface of the sphere that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;

    /** @brief Returns a clone of the cone */
    Convex* clone() const;

    /** @brief Returns the cone volume */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );
    
    /** @ brief Returns whether a point lies inside the cone
    @param pt point */
    bool isIn( Point3 const& pt ) const;    
    //@}
  

  protected:
    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference box,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
    //@} 


    /** @name Parameters */
    //@{
    double m_bottomRadius; /**< radius of the flat base */
    double m_quarterHeight; /**< quarter-height */
    double m_sinAngle; /**< sine of the half-angle at the center */
    //@}  
};

#endif
