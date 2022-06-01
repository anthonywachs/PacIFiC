#ifndef _CYLINDER_HH_
#define _CYLINDER_HH_

#include "Convex.hh"
#include "ReaderXML.hh"


/** @brief The class Cylinder.

    Convex with a circular cylindrical shape. From GJK Engine - A Fast and 
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Cylinder : public Convex 
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Constructor with radius and height as input parameters
    @param r radius
    @param h height */
    Cylinder( double r = 0, double h = 0 );

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    Cylinder( istream& fileIn );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    Cylinder( DOMNode* root );

    /** @brief Destructor */
    ~Cylinder();
    //@}
  

    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns a clone of the cylinder */
    Convex* clone() const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns a vector of points describing the envelope of the
    cylinder. Here simply returns 3 points as follows: center of bottom circular
    face, center of top circular face and an arbitrary point on the lateral
    surface of the cylinder */
    vector<Point3> getEnvelope() const;

    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the point indices. Returns a null pointer as a
    convention */
    vector<vector<int> > const* getFaces() const;

    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape. Here returns the code 777 */
    int getNbCorners() const;
     
    /** @brief Cylinder support function, returns the support point P, i.e. the
    point on the surface of the cylinder that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;

    /** @brief Returns the cylinder volume */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );
  
    /** @brief Returns the number of points to write the cylinder in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;
  
    /** @brief Returns the number of elementary polytopes to write the cylinder 
    in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the cylinder in a
    Paraview format 
    @param f output stream
    @param transform geometric transformation 
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, 
  	Vector3 const* translation = NULL ) const;
	
    /** @brief Returns a list of points describing the cylinder in a
    Paraview format 
    @param transform geometric transformation 
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const; 
  
    /** @brief Writes the cylinder in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const; 
	
    /** @ brief Returns whether a point lies inside the cylinder
    @param pt point */
    bool isIn( Point3 const& pt ) const;	  
    //@}
  

  protected:
    /** @name Parameters */
    //@{  
    double m_radius; /**< cylinder radius */  
    double m_halfHeight; /**< cylinder half height */
    static int m_visuNodeNbOnPer; /**< number of points over the circular edges
    	for Paraview post-processing */  
    //@}


    /**@name Methods */
    //@{
    /** @brief Returns the circumscribed radius of the reference sphere,
    i.e., without applying any transformation */
    double computeCircumscribedRadius() const;
    //@}   
};

#endif
