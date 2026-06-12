#ifndef _POINTC_HH_
#define _POINTC_HH_

#include "Convex.hh"


/** @brief The class PointC.

    Convex with a point shape. 
    
    @author Institut Francais du Petrole - 2000 - Creation
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

    /** @brief Returns the number of points to write the point in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;
  
    /** @brief Returns the number of elementary polytopes to write the point 
    in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the point in a
    Paraview format 
    @param f output stream
    @param transform geometric transformation 
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, 
  	Vector3 const* translation = NULL ) const;
	
    /** @brief Returns a list of points describing the point in a
    Paraview format 
    @param transform geometric transformation 
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const; 
  
    /** @brief Writes the point in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const; 
    
    /** @ brief Returns whether a point lies inside the point (returns false
    by convention)
    @param pt point */
    bool isIn( Point3 const& pt ) const;

    /** @brief Returns the bounding volume to point */
    BVolume* computeBVolume( unsigned int type ) const;
    
    /** @brief Performs advanced comparison of the two points and returns
    whether they match
    @param other the other point */
    bool equalType_level2( Convex const* other ) const;         
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
