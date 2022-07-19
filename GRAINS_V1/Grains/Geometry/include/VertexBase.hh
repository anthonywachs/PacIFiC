#ifndef _VERTEXBASE_HH_
#define _VERTEXBASE_HH_

#include "Point3.hh"
using namespace solid;


/** @brief The class VertexBase.

    To manage an array of vertices. From GJK Engine - A Fast and 
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.
    
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class VertexBase 
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    VertexBase();
  
    /** @brief Constructor with a pointer to the array of vertices as input
    parameter
    @param ptr pointer to the array of vertices */
    VertexBase( void const* ptr );

    /** @brief Destructor */
    ~VertexBase();
    //@}

  
    /** @name Methods */
    //@{
    /** @brief Returns a pointer to the array of vertices */
    void const* getPointer() const;

    /** @brief ith vertex accessor
    @param i vertex index */
    Point3& operator [] ( int i ) const ;
    //@}
  
  private:
    /** @name Parameters */
    //@{
    void const* m_base; /**< pointer to the memory space containing the 
    	vertices */
    //@}
};


// --------------------------------------------------------------------
// Returns a pointer to the array of vertices
inline const void *VertexBase::getPointer() const 
{ 
  return ( m_base ); 
}

// --------------------------------------------------------------------
// ith vertex accessor
inline Point3& VertexBase::operator [] ( int i ) const 
{
  return ( ((Point3 *)m_base)[i] ); 
}

#endif
