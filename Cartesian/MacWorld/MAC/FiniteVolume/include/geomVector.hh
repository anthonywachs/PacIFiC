#ifndef geomVector_HH
#define geomVector_HH

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <iostream>
#include <vector>

//class GE_Point;


/** @brief The Class geomVector.

For the definition of a vector the size of which is in [1:3].

@author A.Wachs - Particulate flow project 2007-2009 */

class geomVector 
{
  private :
     size_t vecSize; /**< number of components */
     double* xx; /**< elements */
    
  public :  
    /** @name Constructors & Destructor */
    //@{    
    /** @brief Constructor without argument */
    geomVector();    
    
    /** @brief Constructor with arguments
    @param vecSize_ number of components */
    geomVector( const size_t &vecSize_ );
    
    /** @brief Constructor with arguments
    @param x x coordinate
    @param y y coordinate */
    geomVector( const double &x, const double &y );
    
    /** @brief Constructor with arguments
    @param x x coordinate
    @param y y coordinate 
    @param z z coordinate */
    geomVector( const double &x, const double &y, const double &z );   
    
    /** @brief Copy constructor */
    geomVector( const geomVector &other );
    
    /** @brief Destructor */
    ~geomVector();
    //@}

    
    /** @name GET methods */
    //@{      
    /** @brief Get a pointer to internal data */
    double const* data() const ;

    /** @brief Get vector size */
    size_t getVecSize() const;        
    //@}

    
    /** @name SET methods */
    //@{ 
    /** @brief Set the vector size and nullify its content
    @param vecSize_ vector size */
    void resize( const size_t &vecSize_ );

    /** @brief Set all components to the same value
    @param val value */
    void set( const double &val );
    
    /** @brief Set a 2D vector components
    @param x x coordinate
    @param y y coordinate */
    void set( const double &x, const double &y );
    
    /** @brief Set a 3D vector components
    @param x x coordinate
    @param y y coordinate 
    @param z z coordinate */
    void set( const double &x, const double &y, const double &z );
    
    /** @brief Add one component of the vector 
    @param index component number 
    @param value value */
    void addOneComp( const size_t index, const double &value );        

    /** @brief Set the vector 
    @param newVec vector to be copied */
    void setVec( const std::vector<double> &newVec );       

    /** @brief Set all components to zero */
    void setVecZero(); 
    
//     /** @brief Set vector components from a GE_Point
//     @param gep GE_Point */
//     void set( const GE_Point *gep );              
    //@}             


    /** @name Operators */
    //@{ 
    /** @brief Operator = 
    @param other right hand side vector */
    geomVector& operator=( geomVector const& other );
       
    /** @brief Operator ==
    @param other right hand side vector */
    bool operator==( geomVector const &other ) const;    
    
    /** @brief Operator * by a scalar 
    @param c the scalar */
    geomVector operator*( const double &c ) const;             

    /** @brief Operator + 
    @param secondVec the other vector to be added */
    geomVector operator+( const geomVector &secondVec ) const; 

    /** @brief Operator - 
    @param secondVec the other vector to be subtracted */
    geomVector operator-( const geomVector &secondVec ) const; 

    /** @brief Operator -= 
    @param secondVec the other vector to be subtracted */
    geomVector operator-=( const geomVector &secondVec ); 

    /** @brief Operator += 
    @param secondVec the other vector to be added */
    geomVector operator+=( const geomVector &secondVec ); 

    /** @brief Operator *= by a scalar 
    @param c the scalar */
    geomVector operator*=( const double &c ); 
    
    /** @brief Operator /= by a scalar 
    @param c the scalar */
    geomVector operator/=( const double &c );

    /** @brief Scalar product (*,*) 
    @param secondVec the other vector to be multiplied with */
    double operator,( const geomVector &secondVec ) const;  
      
    /** @brief Cross product (*^*) 
    @param secondVec the other vector to be multiplied with */
    geomVector operator^( const geomVector &secondVec ) const; 
    
    /** @brief Access operator 
    @param i index */
    double& operator()( size_t i );
    
    /** @brief Access operator 
    @param i index */
    const double& operator()( size_t i ) const;         

    /** @brief Compute the norm of the vector */
    double calcNorm() const;
    
    /** @brief Compute the norm square of the vector */
    double calcNormSquare() const;
        
    /** @brief Compute distance from another vector i.e. norm(*this-secondVec)
    @param secondVec the other vector */
    double calcDist( const geomVector &secondVec ) const;

    /** @brief Compute distance from another vector i.e. norm(*this-secondVec)
    @param secondVec the other vector */
    double calcDist( geomVector const* secondVec ) const;
    
//     /** @brief Compute distance from another GE_Point 
//     @param secondPoint the other point */    
//     double calcDist( const GE_Point* secondPoint ) const;
    
    /** @brief Compute distance from another 1D point defined by its coordinate
    @param x coordinate of the other point */    
    double calcDist( double x ) const;
    
    /** @brief Compute distance from another 2D point defined by its coordinates
    @param x x-coordinate of the other point 
    @param y y-coordinate of the other point */       
    double calcDist( double x, double y ) const;

    /** @brief Compute distance from another 3D point defined by its coordinates
    @param x x-coordinate of the other point 
    @param y y-coordinate of the other point 
    @param z z-coordinate of the other point */
    double calcDist( double x, double y, double z ) const;    
    
    /** @brief Operator translate 
    @param translation_vector translation vector */
    void translate( const geomVector &translation_vector );
    
    /** @brief Operator retract/extend with respect to a point 
    @param gravity_center reference point
    @param amount retraction/extension magnitude */
    void retract( const geomVector &gravity_center, const double &amount );
    //@}    

    
    /** @brief Display a vector 
    @param f output stream 
    @param P the vector to be displayed */
    friend std::ostream& operator <<( std::ostream& f, geomVector const &P );
    
    /** @brief Read a vector 
    @param f input stream 
    @param P the vector to be read */    
    friend std::istream& operator >>( std::istream& f, geomVector &P );  

    /** @brief Return whether the vector is a zero vector i.e. all components
    are zero */ 
    bool isZeroVal() const;
      
    /** @brief Check properties of *this */
    bool check_invariant() const;
    
//     /** @brief Copy coordinates to a GE_Point 
//     @param gep GE_Point */
//     void copy_coordinates( GE_Point *gep ) const;    
           
};

/** @brief Operator * by a scalar 
@param c scalar 
@param secondVec the other vector to be multiplied by the scalar */ 
geomVector operator*( const double &c, const geomVector &secondVec );


//----------------------------------------------------------------------
inline double& 
geomVector::operator()( size_t i )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator()" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK( int(i) >= 0 );
   MAC_CHECK( i < vecSize );
   	
   return xx[i];

}

//----------------------------------------------------------------------
inline const double& 
geomVector::operator()( size_t i ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator()" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK( int(i) >= 0 );
   MAC_CHECK( i < vecSize );
   	
   return xx[i];

}

//----------------------------------------------------------------------
inline size_t 
geomVector::getVecSize() const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: getVecSize" ) ;
   MAC_CHECK_INV( check_invariant() );
	
   return vecSize;
   
} 

//----------------------------------------------------------------------
inline double const*
geomVector::data() const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: getVec" ) ;
   MAC_CHECK_INV( check_invariant() );
	
   return xx;
   
}   

#endif
