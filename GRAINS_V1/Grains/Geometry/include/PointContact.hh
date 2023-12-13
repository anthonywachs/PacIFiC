#ifndef _POINTCONTACT_HH_
#define _POINTCONTACT_HH_

#include "Point3.hh"
#include "Vector3.hh"
using namespace solid;

class Component;
class PointContact;
ostream& operator << ( ostream& Out, PointContact const& pc );

/** @brief The class PointContact.

    Contains all the features of a contact point between 2 rigid bodies.

    @author GRAINS Project - IFP - 2007 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class PointContact
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    PointContact();

    /** @brief Constructor with contact point location in the world reference
    frame, overlap vector, overlap distance and number of iterations of GJK as 
    input parameters
    @param point_  contact point location in the world reference frame
    @param ov_ overlap vector
    @param distance_ overlap distance 
    @param num_iter_ number of iterations of GJK */
    PointContact( Point3 const& point_, Vector3 const& ov_, 
	double distance_, int num_iter_ );

    /** @brief Constructor with contact point location in the world reference
    frame, the point in the 1st rigid body that realizes the minimal distance in
    the rigid body reference frame, the point in the 2nd rigid body that 
    realizes the minimal distance in the rigid body reference frame, overlap 
    vector, overlap distance and number of iterations of GJK as input parameters
    @param point_  contact point location in the world reference frame
    @param pointA_ point in the 1st rigid body that realizes the minimal 
    distance in the rigid body reference frame
    @param pointB_ point in the 2nd rigid body that realizes the minimal 
    distance in the rigid body reference frame
    @param ov_ overlap vector
    @param distance_ overlap distance 
    @param num_iter_ number of iterations of GJK */
    PointContact( Point3 const& point_, Point3 const& pointA_, 
    	Point3 const& pointB_, Vector3 const& ov_, 
	double distance_, int num_iter_ );

    /** @brief Copy constructor
    @param pc_ copied PointContact object */
    PointContact( PointContact const& pc_ );

    /** @brief Destructor */
    ~PointContact();
    //@}


    /** @name Get methods */
    //@{
    /** @brief Returns the contact point in the world reference frame */
    Point3 getContact() const;
  
    /** @brief Returns the point in the 1st rigid body that realizes the minimal
    distance in the rigid body reference frame */
    Point3 getContactA() const;
  
    /** @brief Returns the point in the 2nd rigid body that realizes the minimal
    distance in the rigid body reference frame */
    Point3 getContactB() const;
  
    /** @brief Returns the overlap distance */
    double getOverlapDistance() const;
  
    /** @brief Returns the overlap vector */
    Vector3 getOverlapVector() const;
  
    /** @brief Returns the number of iterations of GJK for convergence */
    int getNbIterGJK() const;
    //@}


    /** @name Set methods */
    //@{
    /** @brief Sets the contact point in the world reference frame
    @param point contact point */
    void setContact( Point3 const& point );
  
    /** @brief Sets the overlap distance
    @param dist overlap distance  */
    void setOverlapDistance( double dist );
  
    /** @brief Sets the overlap vector
    @param vec overlap vector */
    void setOverlapVector( Vector3 const& vec );
  //@}


    /**@name Operators */
    //@{
    /** @brief Equal operator to another PointContact object
    @param rhs the other PointContact object */
    PointContact& operator = ( PointContact const& rhs );
    
    /** @brief Output operator
    @param Out output stream
    @param pc PointContact object */
    friend ostream& operator << ( ostream& Out, PointContact const& pc );
    //@}


  private:
    /** @name Parameters */
    //@{
    Point3 m_contact; /**< contact point location in the world reference 
    	frame */
    Point3 m_contactA; /**< point in the 1st rigid body that realizes the 
    	minimal distance in the rigid body reference frame */
    Point3 m_contactB; /**< point in the 2nd rigid body that realizes the 
    	minimal distance in the rigid body reference frame */
    Vector3 m_overlapvector; /**< overlap vector */
    double m_overlapdistance; /**< overlap distance */
    int m_nbIterGJK; /**< number of iterations of GJK for convergence */
    //@}
};

static PointContact PointNoContact( OriginePoint, OriginePoint, OriginePoint,
	Vector3Nul, 1.e20, 0 );	

#endif
