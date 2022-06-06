#ifndef _INDEXARRAY_HH_
#define _INDEXARRAY_HH_

#include "Basic.hh"

/** @brief The class IndexArray.

    To manage an array of indices. From GJK Engine - A Fast and 
    Robust GJK Implementation, Copyright (C) 1998  Gino van den Bergen.
    
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class IndexArray 
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    IndexArray();

    /** @brief Constructor with the number of elements as an input parameter
    @param n number of elements */
    IndexArray( int n );

    /** @brief Constructor with the number of elements and an array as input 
    parameters
    @param n number of elements 
    @param v array to be copied from 0 to n */
    IndexArray( int n, unsigned int const v[] );

    /** @brief Destructor */
    ~IndexArray();
    //@}
  

    /**@name Methods */
    //@{
    /** @brief Returns a pointer to the array of indices */
    unsigned int* getAdress() const;

    /** @brief Returns the number of elements in the array of indices */
    int size() const;
    //@}


    /**@name Operators */
    //@{
    /** @brief ith element accessor
    @param i element index */
    int operator [] ( int i ) const ;
    //@}
  

  private:
    /**@name Constructors */
    //@{
    /** @brief Copy constructor (forbidden)
    @param Ia object to be copied */
    IndexArray( IndexArray const& Ia );
    //@}


    /**@name Operators */
    //@{
    /** @brief Equal operator (forbidden) 
    @param Ia object to compare to */
    IndexArray& operator = ( IndexArray const& Ia );
    //@}

  
    /** @name Parameters */
    //@{
    unsigned int *m_indices; /**< pointer to the array of indices */
    int m_count; /**< number of elements in the array of indices */
    //@}
};
  
// --------------------------------------------------------------------
// ith element accessor
inline int IndexArray::operator[](int i) const { return m_indices[i]; }

// --------------------------------------------------------------------
// Returns a pointer to the array of indices
inline unsigned int* IndexArray::getAdress() const { return m_indices; }

// --------------------------------------------------------------------
// Returns the number of elements in the array of indices
inline int IndexArray::size() const { return m_count; }

#endif

