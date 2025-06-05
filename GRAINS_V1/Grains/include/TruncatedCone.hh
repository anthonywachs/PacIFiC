#ifndef _TRUNCATEDCONE_HH_
#define _TRUNCATEDCONE_HH_

#include "Convex.hh"
#include "ReaderXML.hh"


/** @brief The class TruncatedCone.

    Convex with a truncated conical shape.

    @author M.BARCET - 2023 - Creation 
    @author A.WACHS - 2025 - Clean up    */
// ============================================================================
class TruncatedCone : public Convex
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Constructor with radius and height as input parameters
    @param r_bot bottom base radius
    @param r_top top base radius
    @param h height */
    TruncatedCone( double r_bot = 0., double r_top = 0., double h = 0. ); 

    /** @brief Constructor with an input stream
    @param fileIn input stream */
    TruncatedCone( istream& fileIn );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    TruncatedCone( DOMNode* root );

    /** @brief Destructor */
    ~TruncatedCone();
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    bool BuildInertia( double* inertia, double* inertia_1 ) const;

    /** @brief Returns the convex type */
    ConvexType getConvexType() const;

    /** @brief Returns a vector of points describing the envelope of the
    truncated cone. Here simply returns 3 points as follows: center of bottom 
    circular face, center of top circular face and an arbitrary point on the 
    lateral surface of the truncated cone */
    vector<Point3> getEnvelope() const;
    
    /** @brief Returns a pointer to a 2D array describing the relationship
    between the face indices and the point indices. Returns a null pointer as a
    convention */
    vector<vector<int> > const* getFaces() const; 
    
    /** @brief Returns the number of vertices/corners or a code corresponding to
    a specific convex shape. Here returns the code 777 */
    int getNbCorners() const;       

    /** @brief Truncated cone support function, returns the support point P, 
    i.e. the point on the surface of the truncated cone that satisfies max(P.v)
    @param v direction vector */
    Point3 support( Vector3 const& v ) const;

    /** @brief Returns a clone of the truncated cone */
    Convex* clone() const;

    /** @brief Returns the truncated cone volume */
    double getVolume() const;

    /** @brief Output operator
    @param fileOut output stream */
    void writeShape( ostream &fileOut ) const;

    /** @brief Input operator
    @param fileIn input stream */
    void readShape( istream &fileIn );

    /** @brief Returns the number of points to write the truncated cone in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the 
    truncated cone in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes a list of points describing the truncated cone in a
    Paraview format
    @param f output stream
    @param transform geometric transformation
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns a list of points describing the truncated cone in a
    Paraview format
    @param transform geometric transformation
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the truncated cone in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const;

    /** @ brief Returns whether a point lies inside the truncated cone
    @param pt point */
    bool isIn( Point3 const& pt ) const;

    /** @brief Writes the truncated cone in an OBJ format
    @param f output stream
    @param transform geometric transformation 
    @param firstpoint_number number of the 1st point */
    void write_convex_OBJ( ostream& f, Transform const& transform,
    	size_t& firstpoint_number ) const;     

    /** @ Returns the bounding volume to truncated cone */
    BVolume* computeBVolume( unsigned int type ) const;
    
    /** @brief Performs advanced comparison of the two truncated cones and 
    returns whether they match
    @param other the other cone */
    bool equalType_level2( Convex const* other ) const;
    
    /** @brief Sets the number of points over the truncated cone perimeter for
    Paraview post-processing, i.e., controls the number of facets in the
    truncated cone reconstruction in Paraview
    @param nbpts number of point over the cone perimeter */
    static void SetvisuNodeNbOverPer( int nbpts ); 
    //@}


  protected:
    /** @name Parameters */
    //@{
    double m_bottomRadius; /**< radius of the flat bottom base */
    double m_topRadius; /**< radius of the flat top base */
    double m_halfHeight; /**< truncated cone half height */
    double m_bottomHeight; /**< from center of mass to bottom base */
    double m_topHeight; /**< from center of mass to top base */
    double m_Hprim; /**< from bottom to top of the full cone */
    double m_sinAngle; /**< sine of the angle of the side wall with bottom to
    	top line */
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
