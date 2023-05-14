#ifndef _SUPERQUADRIC_H_
#define _SUPERQUADRIC_H_

#include "Convex.hh"
#include "ReaderXML.hh"


/** @brief The class Superquadric.

    Convex with a superquadrical shape. From GJK Engine - A Fast and
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.

    @author M.ROULE - Ecole polytechnique - 2020 - Creation
    @author A.Yazdani - 2022 - Major cleaning & refactoring */
// ============================================================================
class Superquadric : public Convex
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur avec initialisation
  @param a scale parameter along the x-axis
  @param b scale parameter along the y-axis
  @param c scale parameter along the z-axis
  @param n1 first exponent
  @param n2 second exponent */
  Superquadric( double a0 = 0., double b0 = 0., double c0 = 0.,
    double n1 = 2., double n2 = 2.);

    /** @brief Constructor with an input stream
    @param fileIn input stream */
  Superquadric( istream& fileIn );

  /** @brief Constructor with an XML node as an input parameter
  @param root XML node */
  Superquadric( DOMNode* root );

  /** @brief Destructor */
  ~Superquadric();
  //@}


  /** @name Methods */
  //@{
  /** @brief Computes the inertia tensor and the inverse of the inertia tensor
  @param inertia inertia tensor
  @param inertia_1 inverse of the inertia tensor */
  bool BuildInertia( double* inertia, double* inertia_1 ) const;

  /** @brief Returns a clone of the superquadric */
  Convex* clone() const;

  /** @brief Returns the convex type */
  ConvexType getConvexType() const;

  /** @brief Returns a vector of points describing the envelope of superquadric.
  Here returns the peak(sharp/corner?) points */
  vector<Point3> getEnvelope() const;

  /** @brief Returns a pointer to a 2D array describing the relationship
  between the face indices and the point indices. Returns a null pointer as a
  convention */
  vector< vector < int > > const* getFaces() const;

  /** @brief Returns the number of vertices/corners or a code corresponding to
  a specific convex shape. Here 777(?) as a convention */
  int getNbCorners() const;

  /** @brief Returns the superquadric volume */
  double getVolume() const;

  /** @brief Output operator
  @param fileOut output stream */
  void writeShape( ostream &fileOut ) const;

  /** @brief Input operator
  @param fileIn input stream */
  void readShape( istream &fileIn );

  /** @brief Superquadric support function, returns the support point P, i.e. the
  point on the surface of the sphere that satisfies max(P.v)
  @param v direction vector */
  Point3 support( Vector3 const& v ) const;

  /** @brief Returns the number of points to write the superquadric in a
  Paraview format */
  int numberOfPoints_PARAVIEW() const;

  /** @brief Returns the number of elementary polytopes to write the superquadric
  in a Paraview format */
  int numberOfCells_PARAVIEW() const;

  /** @brief Writes a list of points describing the superquadric in a
  Paraview format
  @param f output stream
  @param transform geometric transformation
  @param translation additional center of mass translation */
  void write_polygonsPts_PARAVIEW( ostream& f,
    Transform const& transform,
    Vector3 const* translation = NULL ) const;

  /** @brief Returns a list of points describing the superquadric in a
  Paraview format
  @param transform geometric transformation
  @param translation additional center of mass translation */
  list<Point3> get_polygonsPts_PARAVIEW( Transform const& transform,
    Vector3 const* translation = NULL ) const;

  /** @brief Writes the superquadric in a Paraview format
  @param connectivity connectivity of Paraview polytopes
  @param offsets connectivity offsets
  @param cellstype Paraview polytopes type
  @param firstpoint_globalnumber global number of the 1st point
  @param last_offset last offset used for the previous convex shape */
  void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
    int& last_offset ) const;

  /** @brief Returns whether the convex shape is a superquadric */
  bool isSuperquadric() const;

  /** @ brief Returns whether a point lies inside the superquadric
  @param pt point */
  bool isIn( Point3 const& pt ) const;

  /** @ Returns the bounding volume to superquadric */
  BVolume* computeBVolume( unsigned int type ) const;
  //@}

protected:
  /** @name Parameters */
  //@{
  double m_a; /**< scale parameter along the x-axis */
  double m_b; /**< scale parameter along the y-axis */
  double m_c; /**< scale parameter along the z-axis */
  double m_n1; /**< first exponent */
  double m_n2; /**< second exponent */
  static int visuNodeNbOnPer; /**< number of latitudes/longitudes
    (spherical coordinates) for Paraview post-processing */
  //@}


  /**@name Methods */
  //@{
  /** @brief Returns the circumscribed radius of the reference superquadric,
  i.e., without applying any transformation */
  double computeCircumscribedRadius() const;
  //@}
};

#endif
