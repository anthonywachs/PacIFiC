#ifndef _TORSOR_HH_
#define _TORSOR_HH_

#include "Vector3.hh"
#include "Point3.hh"
using namespace solid;


class Torsor;
ostream& operator << ( ostream& fileOut, Torsor const& object );
istream& operator >> ( istream& fileIn, Torsor& object );


/** @brief The class Torsor.

    A force, a point of reference of the torque (generally the center of mass
    of the corresponding rigid body) and the associated torque

    @author F.PRADEL - Institut Francais du Petrole - 1999 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Torsor
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    Torsor();

    /** @brief Constructor with the force, the reference point and the torque as
    input parameters
    @param pt reference point
    @param f force
    @param m torque */
    Torsor( Point3 const& pt, Vector3 const& f, Vector3 const& m );

    /** @brief Constructor with the force and the reference point as
    input parameters. The torque is initialized to (0,0,0)
    @param pt reference point
    @param f force */
    Torsor( Point3 const& pt, Vector3 const& f );

    /** @brief Constructor with the force and the reference point as
    input parameters. The force and the torque are initialized to (0,0,0)
    @param pt reference point */
    Torsor( Point3 const& pt );

    /** @brief Constructor with the force, the reference point and the torque as
    input parameters
    @param point reference point
    @param load force
    @param m torque */
    Torsor( double* point, double* load, double* m );

    /** @brief Copy constructor
    @param tau copied Torsor object */
    Torsor( Torsor const& tau );

    /** @brief Destructor */
    ~Torsor();
    //@}


    /**@name Methods */
    //@{
    /** @brief Modifies the reference point of the torsor and recalculates the
    torque at the new reference point
    @param newpoint new reference point */
    void ChangeReferencePoint( Point3 const& newpoint );
  
    /** @brief Adds a force whose point of application is the reference point of
    the torsor (no torque contribution)
    @param f_ the added force */
    void addForce( Vector3 const& f_ ); 
  
    /** @brief Adds a force whose point of application is different from the 
    reference point of the torsor (additional torque contribution)
    @param f_ the added force 
    @param point point of application of the force */
    void addForce( Point3 const& point, Vector3 const& f_ );   
  
    /** @brief Adds a force whose point of application is the reference point of
    the torsor (no torque contribution)
    @param fx x-component of the force
    @param fy y-component of the force   
    @param fz z-component of the force */
    void addForce( double fx, double fy, double fz );   
  
    /** @brief Adds a torque
    @param torque_ the added torque */
    void addTorque( Vector3 const& torque_ );
  
    /** @brief Adds a torque
    @param mx x-component of the torque
    @param my y-component of the torque
    @param mz z-component of the torque */
    void addTorque( double mx, double my, double mz );   
  
    /** @brief Copies the torsor in a 1D array
    @param fm 1D array where torsor is copied
    @param i start index to copy in the 1D array */
    void copyForceTorque( double* fm, int i ) const;         
    //@}


    /**@name Get methods */
    //@{
    /** @brief Returns a pointer to the total force of the torsor */
    Vector3 const* getForce() const;

    /** @brief Returns a pointer to the total torque of the torsor */
    Vector3 const* getTorque() const;

    /** @brief Returns a pointer to the reference point of the torsor */
    Point3 const* getPoint() const;
    //@}  


    /**@name Set methods */
    //@{
    /** @brief Sets the total force of the torsor
    @param f_ new total force */
    void setForce( Vector3 const& f_ );
    
    /** @brief Sets the reference point of the torsor
    @param point new reference point */
    void setPoint( Point3 const& point );
    
    /** @brief Sets the total torque of the torsor
    @param t_ new total torque */
    void setTorque( Vector3 const& t_ );    
  
    /** @brief Sets the total force of the torsor, the reference point of the
    torsor and initializes the torque to (0,0,0)
    @param point new reference point 
    @param f_ new total force */
    void setToBodyForce( Point3 const& point, Vector3 const& f_ );    
    //@}  


    /**@name Operators */
    //@{
    /** @brief Unitary operator -. Returns a torsor with negative force and
    torque but keeps the reference point unchanged */
    Torsor operator - ();

    /** @brief Multiplication by a scalar of the form Torsor * scalar. Leaves
    the reference point unchanged 
    @param d multiplication factor */
    Torsor operator * ( double d );

    /** @brief Addition of 2 torsors at the reference point of the 1st torsor
    @param k2 torsor to be added */
    Torsor operator + ( Torsor& k2 );

    /** @brief Subtraction of 2 torsors at the reference ponit of the 1st torsor
    @param k2 torsor to be subtracted */
    Torsor operator - ( Torsor &k2 );
  
    /** @brief Comparaison operator
    @param top2 2nd Torsor object */
    bool operator == ( Torsor const& top2 );

    /** @brief Difference operator
    @param top2 2nd Torsor object */
    bool operator != ( Torsor const& top2 );

    /** @brief Equal operator to another Torsor object
    @param rhs the other Torsor object */
    Torsor& operator = ( Torsor const& rhs );

    /** @brief Operator +=. Leaves the reference point unchanged
    @param k2 torsor to be added */
    Torsor& operator += ( Torsor const& k2 ); 

    /** @brief Operator -=. Leaves the reference point unchanged
    @param k2 torsor to be subtracted */
    Torsor& operator -= ( Torsor const& k2 );
    //@}


    /**@name I/O methods */
    //@{
    /** @brief Write the object in a output stream
    @param fileOut output stream */
    void write( ostream& fileOut );
    //@}


    /**@name Friend methods */
    //@{
    /** @brief Output operator
    @param fileOut output stream
    @param object Torsor object */  
    friend ostream& operator << ( ostream& fileOut, Torsor const& object );

    /** @brief Input operator
    @param fileIn input stream
    @param object Torsor object */
    friend istream& operator >> ( istream& fileIn, Torsor& object );
    //@}


  private:
    /**@name Parameter */
    //@{
    Point3 m_refPoint; /**< Reference point of the torsor */
    Vector3 m_totalForce; /**< Sum of all forces */
    Vector3 m_totalTorque; /**< Sum of all torques exerted at the reference 
    	point */
    //@}
};

#endif
