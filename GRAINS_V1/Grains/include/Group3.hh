#ifndef _GROUP3_HH_
#define _GROUP3_HH_

#include "Basic.hh"

#include <iostream>
using namespace std;


class Group3;
ostream &operator << (ostream &fileOut, const Group3 &objet);
istream &operator >> (istream &fileIn, Group3 &objet);


namespace solid
{
  /** @brief The class Group3.
  
  Geometric object with 3 components. From GJK Engine - A Fast and Robust GJK
  Implementation, Copyright (C) 1998  Gino van den Bergen.  
      
  @author D.PETIT - 2000 - Creation 
  @author A.WACHS - 2019 - Update */
  // ==========================================================================
  class Group3
  {
    public:
      /**@name Constructors */
      //@{
      /** @brief Default constructor 
      @param def value of all 3 components */
      Group3( double def = 0. );

      /** @brief Constructor with 3 components as inputs 
      @param x 1st component
      @param y 2nd component
      @param z 3rd component */
      Group3( double x, double y, double z );

      /** @brief Copy constructor
      @param g copied Group3 object */
      Group3( Group3 const& g );

      /** @brief Destructor */
      ~Group3();
      //@}
  

      /**@name Get methods */
      //@{ 
      /** @brief Returns the const pointer to the array */
      double const* getValue() const;
    
      /** @brief Returns the pointer to the array */
      double *getValue();    

      /** @brief Returns the number of components, always return 3 */
      int size() const;    
      //@}
    
    
      /**@name Set methods */
      //@{ 
      /** @brief Modifies the value of the 3 components 
      @param x 1st component
      @param y 2nd component
      @param z 3rd component */
      void setValue( double x, double y, double z );
    
      /** @brief Modifies the value of the 3 components 
      @param g 3-element array */
      void setValue( double const* g ); 
    
      /** @brief Nullifies all components */
      void reset();            
      //@}


      /**@name Methods */
      //@{ 
      /** @brief Writes the object with a high precision format given by
      POSITIONFORMAT defined in GrainsExec.hh
      @param fileOut output stream */
      void writeGroup3( ostream& fileOut ) const;
    
      /** @brief Writes the object in binary format
      @param fileOut output stream */
      void writeGroup3_binary( ostream& fileOut ); 
    
      /** @brief Reads the object in binary format
      @param StreamIN input stream */
      void readGroup3_binary( istream& StreamIN );            
      //@}
                 

      /**@name Operators */
      //@{
      /** @brief ith component accessor
      @param i component index */
      double& operator [] ( size_t i );
    
      /** @brief ith component accessor
      @param i component index */
      double const& operator [] ( size_t i ) const;

      /** @brief Unitary operator -. Returns an object with negative 
      components */
      Group3 operator - () const;

      /** @brief double product
      @param g 2nd Group3 object */
      double operator * ( Group3 const& g ) const;
    
      /** @brief Multiplication by a scalar of the form Group3 * scalar
      @param d multiplication factor */
      Group3 operator * ( double d ) const;

      /** @brief Division by a scalar
      @param d division factor */
      Group3 operator / ( double d ) const;

      /** @brief Addition
      @param g2 2nd Group3 object */
      Group3 operator + ( Group3 const& g2 ) const;

      /** @brief Subtraction
      @param g2 2nd object */
      Group3 operator - ( Group3 const& g2 ) const;

      /** @brief Comparaison operator
      @param g2 2nd Group3 object */
      bool operator == ( Group3 const& g2) const;
    
      /** @brief Difference operator
      @param g2 2nd object */
      bool operator != ( Group3 const& g2 );

      /** @brief Equal operator to another Group3 object
      @param g2 the other Group3 object */
      Group3& operator = ( Group3 const& g2 );
    
      /** @brief Equal operator to a scalar, all components are equal to the
      same scalar
      @param v value set to all 3 components */
      void operator = ( double v );

      /** @brief Unitary operator *= by a scalar
      @param d multiplication factor */
      Group3& operator *= ( double d );

      /** @brief Unitary operator /= by a scalar
      @param d multiplication factor */
      Group3& operator /= ( double d );

      /** @brief Operator += 
      @param g2 2nd Group3 object */
      Group3& operator += ( Group3 const& g2 );

      /** @brief Operator -= 
      @param g2 2nd Group3 object */
      Group3& operator -= ( Group3 const& g2 );
      //@}


      /**@name Friend methods */
      //@{
      /** @brief Mixed product of 3 Group3 objects
      @param g1 1st Group3 object
      @param g2 2nd Group3 object
      @param g3 3rd Group3 object */
      friend double triple( Group3 const& g1, 
	Group3 const& g2, Group3 const& g3 );

      /** @brief Multiplication by a scalar of the form scalar * Group3
      @param d multiplication factor
      @param g Group3 object */
      friend Group3 operator * ( double d, Group3 const& g );

      /** @brief Output operator
      @param fileOut output stream
      @param g Group3 object */
      friend ostream& operator << ( ostream& fileOut, Group3 const& g );

      /** @brief Input operator
      @param fileIn input stream
      @param g Group3 object */
      friend istream& operator >> ( istream& fileIn, Group3& g );
      //@}


    protected:
      /**@name Parameters */
      //@{
      double m_comp[3]; /**< array of 3 components */
      //@}
  };

  /**@name Group3 : External methods */
  //@{
  /** @brief Mixed product of 3 Group3 objects
  @param g1 1st Group3 object
  @param g2 2nd Group3 object
  @param g3 3rd Group3 object */ 
  double triple( Group3 const& g1, Group3 const& g2, Group3 const& g3 );

  /** @brief Multiplication by a scalar of the form scalar * Group3
  @param d multiplication factor
  @param g Group3 object */  
  Group3 operator * ( double d, Group3 const& g );
  //@}
  
  static Group3 OrigineGroup3; /**< Origine (0.,0.,0.)  */

} // end namespace solid

#endif
