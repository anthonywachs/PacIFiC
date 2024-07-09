#ifndef _POINTC_HH_
#define _POINTC_HH_

#include "Convex.hh"


/** @brief The class PointC.

    Convex with a point shape. 
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class PointC : public Convex 
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    PointC();

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    PointC( istream& fileIn );
  
    /** @brief Destructor */
    ~PointC();
    //@}

  
    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the point */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns the point volume, 0 by convention */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Point support function, returns the support point P, i.e. the
    point on the surface of the Point that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;
    
    /** @ brief Returns whether a point lies inside the point (returns false
    by convention)
    @param pt point */
    bool isIn( Point3 const& pt ) const;    

    /** @ Returns the bounding volume to point */
    BVolume* computeBVolume( unsigned int type ) const;
    //@}


  private:
    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference disc,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
  //@}    
};

#endif
