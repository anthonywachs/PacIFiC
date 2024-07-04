#include "Torsor.hh"
#include "Point3.hh"
#include "Vector3.hh"


//-----------------------------------------------------------------------------
// Default constructor         
Torsor::Torsor()
{}




//-----------------------------------------------------------------------------
// Constructor with the force, the reference point and the torque as input 
// parameters
Torsor::Torsor( Point3 const& pt, Vector3 const& f, Vector3 const& m )
{
  m_refPoint = pt;
  m_totalForce = f;
  m_totalTorque = m;
}




//-----------------------------------------------------------------------------
// Constructor with the force and the reference point as
// input parameters. The torque is initialized to (0,0,0)
Torsor::Torsor( Point3 const& pt, Vector3 const& f )
{
  m_refPoint = pt;
  m_totalForce = f;
}




//-----------------------------------------------------------------------------
// Constructor with the force and the reference point as
// input parameters. The force and the torque are initialized to (0,0,0)
Torsor::Torsor( Point3 const& pt )
{
  m_refPoint = pt;
}




//-----------------------------------------------------------------------------
// Constructor with the force, the reference point and the torque as
// input parameters
Torsor::Torsor( double* point, double* load, double* m )
{
  m_refPoint = Point3( point[X], point[Y], point[Z] );
  m_totalForce = Vector3( load[X], load[Y], load[Z] );
  m_totalTorque = Vector3( m[X], m[Y], m[Z] );
}




//-----------------------------------------------------------------------------
// Copy constructor
Torsor::Torsor( Torsor const& tau )
{
  m_refPoint = tau.m_refPoint;
  m_totalForce = tau.m_totalForce;
  m_totalTorque = tau.m_totalTorque;
}




//-----------------------------------------------------------------------------
// Destructor
Torsor::~Torsor()
{}




// ----------------------------------------------------------------------------
// Modify the reference point of the torsor and recalculate the
// torque at the new reference point
void Torsor::ChangeReferencePoint( Point3 const& newpoint )
{
  Vector3 vecteurA = m_refPoint - newpoint;
  m_totalTorque += vecteurA ^ m_totalForce;
  m_refPoint = newpoint;
}




// ----------------------------------------------------------------------------
// Returns a pointer to the reference point of the torsor
Point3 const* Torsor::getPoint() const
{
  return ( &m_refPoint );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the total force of the torsor
Vector3 const* Torsor::getForce() const
{
  return ( &m_totalForce );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the total torque of the torsor
Vector3 const* Torsor::getTorque() const
{
  return ( &m_totalTorque );
}




// ----------------------------------------------------------------------------
// Sets the total force of the torsor
void Torsor::setForce( Vector3 const& f_ )
{
  m_totalForce = f_;
}





// ----------------------------------------------------------------------------
// Sets the reference point of the torsor
void Torsor::setPoint( Point3 const& point )
{
  m_refPoint = point;
}




// ----------------------------------------------------------------------------
// Sets the total torque of the torsor
void Torsor::setTorque( Vector3 const& t_ )
{
  m_totalTorque = t_;
}




// ----------------------------------------------------------------------------
// Sets the total force of the torsor, the reference point of the
// torsor and initializes the torque to (0,0,0)
void Torsor::setToBodyForce( Point3 const& point, Vector3 const& f_ )
{
  m_refPoint = point;  
  m_totalForce = f_;
  m_totalTorque[X] = m_totalTorque[Y] = m_totalTorque[Z] = 0.;
}




// ----------------------------------------------------------------------------
// Adds a force whose point of application is the reference point of
// the torsor (no torque contribution)
void Torsor::addForce( Vector3 const& f_ )
{
  m_totalForce += f_;
}




// ----------------------------------------------------------------------------
// Adds a force whose point of application is different from the 
// reference point of the torsor (additional torque contribution)
void Torsor::addForce( Point3 const& point, Vector3 const& f_ )
{
  m_totalForce += f_;
  Vector3 vecteurA = point - m_refPoint;
  m_totalTorque += (vecteurA ^ f_);
} 




// ----------------------------------------------------------------------------
// Adds a force whose point of application is the reference point of
// the torsor (no torque contribution)
void Torsor::addForce( double fx, double fy, double fz )
{
  m_totalForce[X] += fx;
  m_totalForce[Y] += fy;  
  m_totalForce[Z] += fz;  
} 




// ----------------------------------------------------------------------------
// Adds a torque
void Torsor::addTorque( Vector3 const& torque_ )
{
  m_totalTorque += torque_;
} 




// ----------------------------------------------------------------------------
// Adds a torque
void Torsor::addTorque( double mx, double my, double mz )
{
  m_totalTorque[X] += mx;
  m_totalTorque[Y] += my;  
  m_totalTorque[Z] += mz;  
} 

  


// ----------------------------------------------------------------------------
// Equal operator to another Torsor object
Torsor& Torsor::operator = ( Torsor const& rhs )
{
  if ( &rhs != this )
  {  
    m_refPoint = rhs.m_refPoint;
    m_totalForce = rhs.m_totalForce;
    m_totalTorque = rhs.m_totalTorque;
  }
  return ( *this );
}




// ----------------------------------------------------------------------------
// Addition of 2 torsors at the reference ponit of the 1st torsor
Torsor Torsor::operator + ( Torsor& k2 )
{
  Torsor somme;
  somme.m_refPoint = m_refPoint;
  somme.m_totalForce = m_totalForce + k2.m_totalForce;
  if ( m_refPoint == k2.m_refPoint ) 
  {
    somme.m_totalTorque = m_totalTorque + k2.m_totalTorque;
  } 
  else 
  {
    Vector3 vecteurA = k2.m_refPoint - m_refPoint;
    somme.m_totalTorque = m_totalTorque + k2.m_totalTorque 
    	+ ( vecteurA ^ k2.m_totalForce );
  }   
  return ( somme );
}




// ----------------------------------------------------------------------------
// Subtraction of 2 torsors at the reference ponit of the 1st torsor
Torsor Torsor::operator - ( Torsor& k2 )
{
  Torsor somme;
  somme.m_refPoint = m_refPoint;
  somme.m_totalForce = m_totalForce - k2.m_totalForce;
  if ( m_refPoint == k2.m_refPoint ) 
  {
    somme.m_totalTorque = m_totalTorque - k2.m_totalTorque;
  } 
  else 
  {
    Vector3 vecteurA = k2.m_refPoint - m_refPoint;
    somme.m_totalTorque = m_totalTorque - k2.m_totalTorque 
    	- ( vecteurA ^ k2.m_totalForce );
  }   
  return ( somme );
}




// ----------------------------------------------------------------------------
// Unitary operator -. Returns a torsor with negative force and
// torque but keeps the reference point unchanged
Torsor Torsor::operator - ()
{
  Vector3 tmpforce, tmpmom;
  tmpforce = -1.0 * m_totalForce;
  tmpmom = -1.0 * m_totalTorque;
  
  return ( Torsor( m_refPoint, tmpforce, tmpmom ) );
}




// ----------------------------------------------------------------------------
// Operator -=. Leaves the reference point unchanged
Torsor& Torsor::operator -= ( Torsor const& k2 )
{
  m_totalForce -= k2.m_totalForce;
  if ( m_refPoint == k2.m_refPoint ) 
  {
    m_totalTorque -= k2.m_totalTorque;
  } 
  else 
  {
    Vector3 vecteurA = k2.m_refPoint - m_refPoint;
    m_totalTorque -= k2.m_totalTorque - ( vecteurA ^ k2.m_totalForce );
  }   
  return ( *this );
}




// ----------------------------------------------------------------------------
// Operator +=. Leaves the reference point unchanged
Torsor& Torsor::operator += ( Torsor const& k2 )
{
  m_totalForce += k2.m_totalForce;
  if ( m_refPoint == k2.m_refPoint ) 
  {
    m_totalTorque += k2.m_totalTorque;
  } 
  else 
  {
    Vector3 vecteurA = k2.m_refPoint - m_refPoint;
    vecteurA.round();
    m_totalTorque += k2.m_totalTorque + ( vecteurA ^ k2.m_totalForce );
  }   
  return ( *this );
}




// ----------------------------------------------------------------------------
// Comparaison operator
bool Torsor::operator == ( Torsor const& top2 )
{
  Vector3 vecteurA = top2.m_refPoint - m_refPoint;  
  return ( m_totalTorque == 
  		( top2.m_totalTorque + ( vecteurA ^ top2.m_totalForce ) ) 
  	&& m_totalForce == top2.m_totalForce );
}




// ----------------------------------------------------------------------------
// Difference operator
bool Torsor::operator != ( Torsor const& top2 )
{
 return ( !(*this == top2) );
}




// ----------------------------------------------------------------------------
// Multiplication by a scalar of the form Torsor * scalar. Leaves
// the reference point unchanged
Torsor Torsor::operator * ( double d )
{
  Torsor multi;
  multi.m_refPoint = m_refPoint;
  multi.m_totalForce = m_totalForce * d;
  multi.m_totalTorque = m_totalTorque * d;
  return ( multi );
}




// ----------------------------------------------------------------------------
// Write the object in a output stream
void Torsor::write( ostream& fileOut )
{
  fileOut << "*PtContact\t" << m_refPoint;
  fileOut << "*Fn+Ft\t"     << m_totalForce;
  if ( Norm( m_totalTorque ) < 1.e20 )
    fileOut << "*Moment\t"    << m_totalTorque;
  else
    fileOut << "*Moment   0.   0.   0.\n";
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, Torsor const& object )
{
  fileOut << "Point = " << object.m_refPoint
	  << "Force = " << object.m_totalForce
	  << "Moment = " << object.m_totalTorque;
  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, Torsor& object )
{
  fileIn >> object.m_refPoint >> object.m_totalForce >> object.m_totalTorque;
  return ( fileIn );
}




// ----------------------------------------------------------------------------
// Copies the torsor in a 1D array
void Torsor::copyForceTorque( double* fm, int i ) const
{
  for (int j=0 ;j<3; j++) fm[i+j] = m_totalForce[j];
  for (int j=0 ;j<3; j++) fm[i+3+j] = m_totalTorque[j];    
} 
