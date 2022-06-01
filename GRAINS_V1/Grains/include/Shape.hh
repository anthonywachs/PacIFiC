#ifndef _SHAPE_HH_
#define _SHAPE_HH_

#include "BBox.hh"
#include "Point3.hh"
using namespace solid;

class Transform;

enum ShapeType {
  CONVEX,  
  COMPLEX
};


/** @brief The class Shape.

    High level class describing the type of shape. Only the type convex is
    implemented in Grains3D. From GJK Engine - A Fast and Robust GJK
    Implementation, Copyright (C) 1998  Gino van den Bergen.

    @author F.PRADEL - Institut Francais du Petrole - 2000 - Modification 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Shape 
{
  public:
    /**@name Constructor & Destructor */
    //@{
    /** @brief Default constructor */
    Shape();
  
    /** @brief Destructor */
    virtual ~Shape();
    //@}


    /**@name Purely virtual methods */
    //@{
    /** @brief Returns the BBox 
    @param t transformation associated to the shape */
    virtual BBox bbox( Transform const& t ) const = 0;

    /** @brief Returns the type of shape  */
    virtual ShapeType getType() const = 0;
    //@}
};

typedef Shape const* ShapePtr;

typedef bool (*Intersect)( Shape const&, Shape const&, 
	Transform const&, Transform const&, Vector3& );

typedef bool (*Common_point)( Shape const&, Shape const&, 
	Transform const&, Transform const&,
	Vector3& , Point3& , Point3& );

#endif
