#ifndef _SEGMENT_HH_
#define _SEGMENT_HH_

#include "Convex.hh"


/** @brief The class Sphere.

    Convex with a segment shape. 
    
    @author Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Segment : public Convex 
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Constructor with the length as an input parameter
    @param x segment length */
    Segment( double x = 0 );

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    Segment( istream& fileIn );

    /** @brief Destructor */
    ~Segment();
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia
    tensor. Inertia and inverse of inertia are 0 by convention
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;
  
    /** @brief Returns a clone of the segment */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns the segment length */
    double getLength() const;

    /** @brief Returns the segment volume, here 0 by convention */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Segment support function, returns the support point P, i.e. the
    point on the surface of the sphere that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;
  
    /** @brief Return the transformation associated to a direction vector and a
    point correponding to the center of mass of the segment 
    @param v direction vector
    @param gc center of mass coordinates */
    static Transform computeTransform( Vector3 const& v, Point3 const& gc );

    /** @brief Returns the number of points to write the segment in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;
  
    /** @brief Returns the number of elementary polytopes to write the segment 
    in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the segment in a
    Paraview format 
    @param f output stream
    @param transform geometric transformation 
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, 
  	Vector3 const* translation = NULL ) const;
	
    /** @brief Returns a list of points describing the segment in a
    Paraview format 
    @param transform geometric transformation 
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const; 
  
    /** @brief Writes the segment in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const; 

    /** @ brief Returns whether a point lies inside the segment
    @param pt point */
    bool isIn( Point3 const& pt ) const;

    /** @ Returns the bounding volume to segment */
    BVolume* computeBVolume( unsigned int type ) const;
    
    /** @brief Performs advanced comparison of the two segments and returns
    whether they match
    @param other the other segment */
    bool equalType_level2( Convex const* other ) const;    
    //@}


  protected:
    /** @name Parameters */
    //@{
    double m_halflength; /**< half length of the segment */
    //@}


    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference sphere,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
    //@}   
};

#endif
