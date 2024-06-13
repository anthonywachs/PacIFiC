#ifndef _STLOBSTACLE_HH_
#define _STLOBSTACLE_HH_

#include "SimpleObstacle.hh"
#include "STLVertex.hh"
#include "STLTriangle.hh"
#include <list>
#include <set>
using namespace std;

#include "ReaderXML.hh"


/** @brief The class STLObstacle.

    A simple obstacle defined by an STL triangulation.

    @author A.MORENTE - 2023 - Creation */
// ============================================================================
class STLObstacle : public SimpleObstacle
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Constructor with name as input parameter
    @param s obstacle name */
    STLObstacle( string const& s );
    
    /** @brief Constructor with name and file name as input parameters
    @param s obstacle name 
    @param filename file name */
    STLObstacle( string const& s, string const& filename );    

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    STLObstacle( DOMNode* root );

    /** @brief Destructor */
    ~STLObstacle();
    //@}


    /** @name Accessors */
    //@{
    //@}


    /** @name Methods */
    //@{    
    /** @brief Reads the STL from a file
    @param filename file name */
    void readSTL( string const& filename );

    /** @brief Is a point already a vertex in the list of all vertices ?
    @param x x-coordinate of the point
    @param y y-coordinate of the point
    @param z z-coordinate of the point */
    int isPInV( double x, double y, double z ); 
    
    /** @brief Searches a vertex by its id number, returns a pointer to the
    vertex if found, returns NULL otherwise
    @param ide id number */
    STLVertex* VidSearch( size_t ide );                         

    //int isVinT(STLVertex &vx, STLTriangle &T);
    //int isVinV(STLVertex &v1, STLVertex &v2);

    /** @brief Moves the STL obstacle and returns a list of moved
    obstacles (here itself)
    @param time physical time
    @param dt time step magnitude
    @param motherCompositeHasImposedVelocity whether the composite that the 
    obstacle belongs to has a non-zero imposed velocity
    @param motherCompositeHasImposedForce whether the composite that the 
    obstacle belongs to has a non-zero imposed force */
    virtual list<SimpleObstacle*> Move( double time,
	double dt, bool const& motherCompositeHasImposedVelocity,
        bool const& motherCompositeHasImposedForce ) ;

    /** @brief Contact between a simple obstacle and a component. If contact
    exists, computes the contact force and torque and adds to each component
    @exception ContactError if overlapping distance is larger than the sum of
    the crust thicknesses of the components
    @param voisin the other component
    @param dt time step magnitude
    @param time physical time
    @param LC linked-cell grid */
    virtual void InterAction( Component* voisin,
	double dt, double const& time, LinkedCell* LC );

    /** @brief Searches and stores all contact points between two components
    @exception ContactError if overlapping distance is larger than the sum of
    the crust thicknesses of the components
    @param voisin the other component
    @param dt time step magnitude
    @param time physical time
    @param LC linked-cell grid
    @param listContact list of information about contacts */
    void SearchContact( Component* voisin, double dt,
      double const& time, LinkedCell *LC,
      list<ContactInfos*>& listContact );

    /** @brief Returns whether there is geometric contact with another
    component 
    @param voisin the other component */
    virtual bool isContact( Component const* voisin ) const;

    /** @brief Returns whether there is geometric contact with another
    component accounting for crust thickness 
    @param voisin the other component */
    virtual bool isContactWithCrust( Component const* voisin ) const;

    /** @brief Returns whether there is geometric proximity with another
    component in the sense of whether their respective bounding boxes overlap 
    @param voisin the other component */
    virtual bool isClose( Component const* voisin ) const;

    /** @brief Returns whether there is geometric proximity with another
    component in the sense of whether their respective bounding boxes minus 
    their crust thickness overlap 
    @param voisin the other component */
    virtual bool isCloseWithCrust( Component const* voisin ) const; 

    /** @brief Rotates the obstacle with a quaternion
    @param rotation the quaternion defining the rotation */
    virtual void Rotate( Quaternion const& rotation );

    /** @brief Translates the obstacle
    @param translation translation vector */
    virtual void Translate( Vector3 const& translation );

    /** @brief Returns the maximum of the absolute value of the obstacle
    velocity in each direction */
    virtual Vector3 vitesseMaxPerDirection() const;
  
    /** @brief Returns whether the component is an STL obstacle */
    virtual bool isSTLObstacle() const;  
    //@}


    /** @name Set methods */
    //@{
    /** @brief Sets the STL normal vectors at all vertices */    
    void STLnormalsAtVertices();
    //@}


    /** @name I/O methods */
    //@{
    /** @brief Reloads the STL obstacle and links it to the higher level
    obstacle in the obstacle tree
    @param mother higher level obstacle
    @param file input stream */
    virtual void reload( Obstacle& mother, istream& file ) ;

    /** @brief Outputs the STL obstacle for reload
    @param fileSave output stream */
    virtual void write( ostream& fileSave ) const;

    /** @brief Returns the number of points to write the STL obstacle in a
    Paraview format */
    virtual int numberOfPoints_PARAVIEW() const ;

    /** @brief Returns the number of elementary polytopes to write the
    STL obstacle shape in a Paraview format */
    virtual int numberOfCells_PARAVIEW() const;

    /** @brief Returns a list of points describing the STL obstacle in a
    Paraview format
    @param translation additional center of mass translation */
    virtual list<Point3> get_polygonsPts_PARAVIEW(
  	Vector3 const* translation = NULL ) const ;

    /** @brief Writes the points describing the STL obstacle in a
    Paraview format
    @param f output stream
    @param translation additional center of mass translation */
    virtual void write_polygonsPts_PARAVIEW( ostream& f,
  	Vector3 const* translation = NULL ) const ;

    /** @brief Writes the STL obstacle in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    virtual void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const ;

    /**  @brief Outputs information to be transferred to the fluid
    @param fluid output stream */
    virtual void writePositionInFluid( ostream& fluid );
    
    /** @brief Displays all vertices */
    void displayAllSTLVertices(); 
    
    /** @brief Displays all triangles */           
    void displayAllSTLTriangles();    
    //@}


    /** @name STL geometric utilities */
    //@{
    /** @brief ???
    @param P ???
    @param P1 ???
    @param P2 ??? 
    @param P3 ??? */
    static int intersect( Point3 P, Point3 P1, Point3 P2, Point3 P3 );
    
    /** @brief Returns the dot product of two vectors
    @param vect_A[] first vector
    @param vect_B[] second vector */
    static double dotProduct( double vect_A[], double vect_B[] );

    /** @brief Computes the cross product of two vectors
    @param vect_A[]  first vector
    @param vect_B[]  second vector
    @param cross_P[] result */
    static void crossProduct( double vect_A[],
	double vect_B[],
	double cross_P[] );

    /** @brief Computes the difference between two vectors
    @param vect_A[] first vector
    @param vect_B[] second vector
    @param diff_P[] result */
    static void diffProduct( double vect_A[],
	double vect_B[],
	double diff_P[] );

    /** @brief Computes the sign of the signed volume of the tetrahedron 
    (A,B,C,D): 1 -> positive or 0 -> negative
    @param vect_A[] first vector
    @param vect_B[] second vector
    @param vect_C[] third vector
    @param vect_D[] fourth vector */
    static int orient3d( double vect_A[],
	double vect_B[],
	double vect_C[],
	double vect_D[] );

    /** @brief Returns 1 if a segment [q1 q2] intersects a triangle 
    (tri1, tri2, tri3)
    @param q1[] first point
    @param q2[] second point
    @param tri1[] first vertex
    @param tri2[] second vertex
    @param tri3[] third vertex */
    static int intersect3d(double q1[], 
	double q2[], 
	double tri1[], 
	double tri2[], 
	double tri3[] );    
    //@}


  protected:    
    /**@name Parameters */
    //@{
    int m_npls; /**< Number of triangles */
    vector<tuple<double,double,double>> m_llvls; /**< vector containing the
    	vertices of the triangulation */
    vector<tuple<double,double,double>> m_llvns; /**< vector containing the
    	normals associated to the triangles i.e. one per 3 vertices */
    list<STLVertex*> m_allSTLVertices; /**< list containing the
    	vertices of the triangulation */
    vector<STLTriangle> m_allSTLTriangles; /**< vector containing the
    	vertices of the triangulation */    
    //@}
    
    
    /** @name Methods */
    //@{  
    /** @brief Computes center of mass position */
    pair<Point3,double> computeCenterOfMass();
    //@}    

};

#endif
