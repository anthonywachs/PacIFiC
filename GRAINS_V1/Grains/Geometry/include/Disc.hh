#ifndef _DISC_HH_
#define _DISC_HH_

#include "Convex.hh"
#include "ReaderXML.hh"


/** @brief The class Disc.

    Convex with a 2D disc shape. From GJK Engine - A Fast and 
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.
    
    @author Institut Francais du Petrole - 2003 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Disc : public Convex 
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Constructor with radius as input parameter
    @param r disc radius */
    Disc( double r = 0 );

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    Disc( istream& fileIn );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    Disc( DOMNode* root );
  
    /** @brief Destructor */
    ~Disc();
    //@}
  

    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the disc */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns a vector of points describing the surface of the disc.
    Here simply returns the point (0,0,0) as a convention */
    vector<Point3> getEnvelope() const;

    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape. Here 1 as a convention */
    int getNbCorners() const;  

    /** @brief Returns the disc surface area */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Disc support function, returns the support point P, i.e. the
    point on the surface of the disc that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;
  
    /** @brief Returns the number of points to write the disc in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;
  
    /** @brief Returns the number of elementary polytopes to write the sphere 
    in a Paraview format */
    int numberOfCells_PARAVIEW() const;
    
    /** @brief Writes a list of points describing the disc in a
    Paraview format 
    @param f output stream
    @param transform geometric transformation 
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, 
  	Vector3 const* translation = NULL ) const;
	
    /** @brief Returns a list of points describing the disc in a
    Paraview format 
    @param transform geometric transformation 
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const; 
  
    /** @brief Writes the disc in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const; 

    /** @brief Sets the number of points over the disc perimeter for
    Paraview post-processing, i.e., controls the number of facets in the
    disc reconstruction in Paraview
    @param nbpts number of point over the cylinder perimeter */
    static void SetvisuNodeNbOverPer( int nbpts ); 
    
    /** @ brief Returns whether a point lies inside the disc
    @param pt point */
    bool isIn( Point3 const& pt ) const;

    /** @ Returns the bounding volume to disc */
    BVolume* computeBVolume( unsigned int type ) const;
    
    /** @brief Performs advanced comparison of the two discs and returns
    whether they match
    @param other the other disc */
    bool equalType_level2( Convex const* other ) const;     
    //@}
  

  protected:
    /**@name Parameters */
    //@{
    double m_radius; /**< disc radius */
    static int m_visuNodeNb; /**< number of points over the perimeter of the 
    	disc for Paraview post-processing */
    //@}  
  
  
    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference disc,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
    //@} 
};

#endif
