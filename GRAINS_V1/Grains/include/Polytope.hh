#ifndef _POLYTOPE_HH_
#define _POLYTOPE_HH_

#include "Convex.hh"
#include "IndexArray.hh"
#include "VertexBase.hh"
#include <string>
using namespace std;


/** @brief The class Polytope.

    Convex with a polygonal/polyedral shape. From GJK Engine - A Fast and 
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.
    
    @author D.PETIT - Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Polytope : public Convex 
{
  public:
    /** @name Virtual methods */
    //@{
    /** @brief Computes the inertia tensor and the inverse of the inertia tensor
    @param inertia inertia tensor
    @param inertia_1 inverse of the inertia tensor */
    virtual bool BuildInertia( double* inertia, double* inertia_1 ) const = 0;

    /** @brief Returns a clone of the polytope */
    virtual Convex* clone() const = 0;

    /** @brief Returns the polytope volume */
    virtual double getVolume() const = 0;

    /** @brief Sphere support function, returns the support point P, i.e. the
    point on the surface of the sphere that satisfies max(P.v)
    @param v direction vector */
    virtual Point3 support( Vector3 const& v ) const = 0;
    //@}
  

    /**@name Methods */
    //@{
    /** @brief Returns the number of vertices of the polytope */
    int numVerts() const;

    /** @brief Returns a vector of points describing the envelope of the
    polytope */
    vector<Point3> getEnvelope() const;
  
    /** @brief Returns the number of vertices/corners */
    int getNbCorners() const; 
  
    /** @brief Returns the number of points to write the polytope in a
    Paraview format */
    virtual int numberOfPoints_PARAVIEW() const;
  
    /** @brief Writes a list of points describing the polytope in a
    Paraview format 
    @param f output stream
    @param transform geometric transformation 
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, 
  	Vector3 const* translation = NULL ) const;
	
    /** @brief Returns a list of points describing the polytope in a
    Paraview format 
    @param transform geometric transformation 
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation = NULL ) const; 
	
    /** @ brief Returns whether a point lies inside the polytope
    @param pt point */
    virtual bool isIn( Point3 const& pt ) const;

    /** @ Returns the bounding volume to polytope */
    virtual BVolume* computeBVolume( unsigned int type ) const;
    
    /** @brief Performs advanced comparison of the two convexes and returns
    whether they match
    @param other the other convex */
    virtual bool equalType_level2( Convex const* other ) const;     	
    //@}

  
    /**@name Operators */
    //@{
    /** @brief ith vertex accessor
    @param i vextex index */
    Point3 const& operator [] ( int i ) const;

    /** @brief ith vertex accessor
    @param i vextex index */
    Point3& operator [] ( int i );
    //@}


  protected:
    /**@name Parameters */
    //@{
    VertexBase const& m_base; /**< reference to the array containing the
    	vertices of the polytope */
    IndexArray const& m_index; /**< reference to the index array giving the
    	index of a vertex in the array of vertices m_base */
    string m_fichPoly; /**< file containing the vertices and faces describing
    	the polytope */
    //@}


    /**@name Constructors */
    //@{
    /** @brief Constructor with an input stream, a number of vertices, a
    reference to the array of vertices and a reference to the array of indices
    as input parameters
    @param fileIn input stream
    @param nb_point number of vertices
    @param ref reference to the array of vertices
    @param ia reference to the array of indices */
    Polytope( istream& fileIn, int nb_point, VertexBase& ref, IndexArray& ia );

    /** @brief Copy constructor
    @param copy copied Polytope object  */
    Polytope( Polytope const& copy );
 
    /** @brief Destructor */
    virtual ~Polytope();
    //@}

  
  private:
    /**@name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    Polytope();
    //@}
  
};

#endif
