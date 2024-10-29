#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "Transform.hh"
#include <string>
#include <iostream>
#include <sstream>
using namespace std;

size_t Transform::m_sizeofTransform = Matrix::m_sizeofMatrix + 
	solid::Group3::m_sizeofGroup3 + sizeof( unsigned int );


// --------------------------------------------------------------------------
// Default constructor. Origin is (0,0,0) and matrix is identity
Transform::Transform() 
  : m_type( IDENTITY )
{
  m_basis.setIdentity();
  m_origin.setValue( 0, 0, 0 );
}




// --------------------------------------------------------------------------
// Constructor with origin coordinates as inputs and matrix is identity
Transform::Transform( double gx, double gy, double gz )
{
  m_basis.setIdentity();
  m_origin.setValue( gx, gy, gz );
  if ( ( gx * gx + gy * gy + gz * gz ) != 0. ) m_type = TRANSLATION;
  else m_type = IDENTITY;
}




// --------------------------------------------------------------------------
// Constructor with a 1D array of 12 values as inputs containing the
// rotation matrix coefficients following by the origin coordinates
Transform::Transform( double const m[12] ) 
{ 
  setValue( m );    
}




// --------------------------------------------------------------------------
// Copy constructor
Transform::Transform( Transform const& other )
{
  m_basis = other.m_basis;
  m_origin = other.m_origin;
  m_type = other.m_type;  
}




// --------------------------------------------------------------------------
// Destructor
Transform::~Transform()
{}




// --------------------------------------------------------------------------
// Set the transformation to the inverse of another transformation t
void Transform::setToInverseTransform( Transform const& t ) 
{
  m_basis = t.m_type & SCALING ? inverse( t.m_basis ) : transpose( t.m_basis );
  m_origin.setValue( - m_basis[X] * t.m_origin, - m_basis[Y] * t.m_origin, 
	- m_basis[Z] * t.m_origin );  
  m_type = t.m_type;
}




// --------------------------------------------------------------------------
// Returns true if the transformation type is IDENTITY
bool Transform::isIdentity() const
{
  return ( m_type == IDENTITY );
}




// --------------------------------------------------------------------------
// Composition of affine transformations: this = t2 o t1 (t1 first
// followed by t2)
void Transform::setToCompositionOfTransform( Transform const& t1, 
	Transform const& t2 ) 
{  
  m_basis = t2.m_basis * t1.m_basis;
  m_origin = t2(t1.m_origin);
  m_type = t2.m_type | t1.m_type;  
}




// --------------------------------------------------------------------------
// Composition with a scaling transformation: this = this o scaling
void Transform::composeWithScaling( double x, double y, double z )
{
  m_basis.multiplyByScalingMatrix( x, y, z );
  m_type |= SCALING;  
} 




// --------------------------------------------------------------------------
// Output operator (calls <<)
void Transform::writeShape( ostream& fileOut ) const 
{
  fileOut << *this;
}




// --------------------------------------------------------------------------
// Input operator (calls >>)
void Transform::readShape( istream& fileIn ) 
{
  fileIn >> *this;
}




// --------------------------------------------------------------------------
// Composition on the left by a rotation described by a quaternion:
// this = rot(quaternion) o this (this first followed by rot(quaternion))
void Transform::composeLeftByRotation( const Quaternion& q ) 
{
  Mat3 const& mm = m_basis.getValue(); 
  double qx = q[X], qy = q[Y], qz = q[Z], qw = q[W];
  double px, py, pz, pw;
  
  // We compute below the matrix product M_rot(q).basis where
  // M_rot(q) is the rotation matrix corresponding to the quaternion q
  for (int i=0;i<3;++i)
  {
    px = qy *  mm[Z][i] - qz * mm[Y][i] + qw * mm[X][i];
    py = qz *  mm[X][i] - qx * mm[Z][i] + qw * mm[Y][i];    
    pz = qx *  mm[Y][i] - qy * mm[X][i] + qw * mm[Z][i]; 
    pw = - qx * mm[X][i] - qy * mm[Y][i] - qz * mm[Z][i];
    m_basis[X][i] = qy * pz - qz * py - pw * qx + qw * px;
    m_basis[Y][i] = qz * px - qx * pz - pw * qy + qw * py;          
    m_basis[Z][i] = qx * py - qy * px - pw * qz + qw * pz;    
  }
  m_type |= ROTATION;    
}




// --------------------------------------------------------------------------
// Sets the transformation to the identity
void Transform::setIdentity() 
{
  m_basis.setIdentity();
  m_origin.setValue( 0., 0., 0. );
  m_type = IDENTITY;
}




// --------------------------------------------------------------------------
// Sets the origin of the transformation
void Transform::setOrigin( double const* pos ) 
{
  m_origin[X] = pos[X];
  m_origin[Y] = pos[Y];
  m_origin[Z] = pos[Z];
  m_type |= TRANSLATION;
}




// --------------------------------------------------------------------------
// Sets the origin of the transformation
void Transform::setOrigin( double gx, double gy, double gz ) 
{
  m_origin[X] = gx;
  m_origin[Y] = gy;
  m_origin[Z] = gz;
  m_type |= TRANSLATION;  
}




// --------------------------------------------------------------------------
// Sets the origin of the transformation
void Transform::setOrigin( Point3 const& pos ) 
{
  m_origin = pos;
  m_type |= TRANSLATION;  
}




// --------------------------------------------------------------------------
// Sets the matrix part of the transformation
void Transform::setBasis( const Matrix &basis_ )
{
  m_basis = basis_;
  bool OriginIsZero = true ;
  if ( ( m_origin[X] * m_origin[X] + m_origin[Y] * m_origin[Y] 
  	+ m_origin[Z] * m_origin[Z] ) > EPSILON2 ) OriginIsZero = false ;
  if ( m_basis.isDiagonal() )
  {
    if ( m_basis.isIdentity() ) 
    {
      if ( OriginIsZero ) m_type = IDENTITY ;
      else m_type = TRANSLATION ; 
    }
    else
      m_type = SCALING ;
  }
  else
  {    
    if ( OriginIsZero )    
      if ( fabs( m_basis.determinant() - 1. ) < EPSILON2 )
        m_type = ROTATION ;
      else
        m_type = LINEAR ;
    else m_type = AFFINE ;
  }   
}




// --------------------------------------------------------------------------
// Sets the transformation with an 1D array of 12 values as inputs 
void Transform::setValue( const double m[12] ) 
{
  m_basis.setValue( m );
  m_origin.setValue( &m[9] );
  bool OriginIsZero = true ;
  if ( ( m[9] * m[9] + m[10] * m[10] + m[11] * m[11] ) > EPSILON2 ) 
    OriginIsZero = false ;
  if ( m_basis.isDiagonal() )
  {
    if ( m_basis.isIdentity() ) 
    {
      if ( OriginIsZero ) m_type = IDENTITY ;
      else m_type = TRANSLATION ; 
    }
    else
      m_type = SCALING ;
  }
  else
  {    
    if ( OriginIsZero )    
      if ( fabs( m_basis.determinant() - 1. ) < EPSILON2 )
        m_type = ROTATION ;
      else
        m_type = LINEAR ;
    else m_type = AFFINE ;
  } 
}




// --------------------------------------------------------------------------
// Composition on the left by a translation:
// this = trans(vector) o this (this first followed by trans(vector))
void Transform::composeLeftByTranslation( Vector3 const& v )
{ 
  m_origin += v;
  m_type |= TRANSLATION;
}




// --------------------------------------------------------------------------
// Returns a point corresponding to applying the transformation to the 
// input point
Point3 Transform::operator () ( Point3 const& p ) const 
{
  return ( Point3( m_basis[X] * p + m_origin[X], 
	m_basis[Y] * p + m_origin[Y], 
	m_basis[Z] * p + m_origin[Z] ) );
}




// --------------------------------------------------------------------------
// Equal operator to another Transform object
Transform& Transform::operator = ( Transform const& transform_ )
{
  if ( &transform_ != this )
  {  
    m_basis = transform_.m_basis;
    m_origin = transform_.m_origin;
    m_type = transform_.m_type;
  }
  
  return ( *this );
}




// --------------------------------------------------------------------------
// Composition on the right by another affine transformation: 
// this = this o t (t first followed by this)
void Transform::composeRightByTransform( Transform const& t ) 
{
  m_origin += m_basis * t.m_origin;
  m_basis *= t.m_basis;
  m_type |= t.m_type; 
}




// --------------------------------------------------------------------------
// Composition on the left by another affine transformation: 
// this = t o this (this first followed by t)
void Transform::composeLeftByTransform( Transform const& t ) 
{
   m_origin = t.m_origin + t.m_basis * m_origin;
   m_basis = t.m_basis * m_basis ;
   m_type = t.m_type | m_type ;      
}




// --------------------------------------------------------------------------
// Composition on the left by a rotation described by a transform:
// this = rot o this (this first followed by rot)
// This composition leaves the origin unchanged but does not check that rot 
// is indeed a rotation 
void Transform::composeLeftByRotation( Transform const& t ) 
{
   m_basis = t.m_basis * m_basis ;
   m_type = t.m_type | m_type ;      
}




// --------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, Transform const& t )
{  
  fileOut << "Type = " << t.m_type << endl;
  fileOut << "*Position\n";
  fileOut << t.m_origin << endl;
  fileOut << "*Orientation\n";
  fileOut << t.m_basis;

  return ( fileOut );
}




// --------------------------------------------------------------------------
// Writes the object with a high precision format given by
// FORMAT16DIGITS defined in GrainsExec.hh
void Transform::writeTransform( ostream& fileOut ) const
{
  fileOut << "*Position " << m_type << endl;
  fileOut << GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  	m_origin[X]) << " " << 
	GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  	m_origin[Y]) << " " << 
	GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  	m_origin[Z]) << endl;	
  m_basis.writeMatrix( fileOut );
}




// --------------------------------------------------------------------------
// Writes the object with a high precision format given by
// FORMAT16DIGITS defined in GrainsExec.hh and the 2014 reload format
void Transform::writeTransform2014( ostream& fileOut ) const
{
  fileOut << GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  	m_origin[X]) << " " << 
	GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  	m_origin[Y]) << " " << 
	GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  	m_origin[Z]) << " ";	
  m_basis.writeMatrix2014( fileOut );
  fileOut << " " << m_type;
}




// --------------------------------------------------------------------------
// Writes the object in binary format with the 2014 reload format
void Transform::writeTransform2014_binary( ostream& fileOut )
{
  fileOut.write( reinterpret_cast<char*>( &m_origin[X] ), sizeof( double ) );
  fileOut.write( reinterpret_cast<char*>( &m_origin[Y] ), sizeof( double ) );  
  fileOut.write( reinterpret_cast<char*>( &m_origin[Z] ), sizeof( double ) );
  m_basis.writeMatrix2014_binary( fileOut );
  fileOut.write( reinterpret_cast<char*>( &m_type ), sizeof( unsigned int ) );
}




// --------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, Transform& t )
{
  string cle;
  
  fileIn >> t.m_type;  
  fileIn >> t.m_origin;
  fileIn >> t.m_basis;

  return ( fileIn );
}




// --------------------------------------------------------------------------
// Reads the object with the 2014 reload format
void Transform::readTransform2014( istream& StreamIN )
{
  StreamIN >> m_origin >> m_basis >> m_type;
}




// --------------------------------------------------------------------------
// Reads the object in binary format with the 2014 reload format
void Transform::readTransform2014_binary( istream& StreamIN )
{
  StreamIN.read( reinterpret_cast<char*>( &m_origin[X] ), sizeof( double ) );
  StreamIN.read( reinterpret_cast<char*>( &m_origin[Y] ), sizeof( double ) );  
  StreamIN.read( reinterpret_cast<char*>( &m_origin[Z] ), sizeof( double ) );
  m_basis.readMatrix2014_binary( StreamIN );
  StreamIN.read( reinterpret_cast<char*>( &m_type ), sizeof( unsigned int ) );
}




// --------------------------------------------------------------------------
// Initialization of the transformation from an XML node
// The origin of the transformation is imposed at (0,0,0)
// The matrix must be a rotation matrix, we do not allow stretching/scaling
void Transform::load( DOMNode* root )
{
  // Initialization of the 1D array containing all the transformation
  // coefficients ordered as Mxx, Mxy, Mxz, Myx, Myy, Myz, Mzx, Mzy, Mzz,
  // Ox, Oy, Oz  
  double t[12] = {1., 0., 0., 
		  0., 1., 0., 
		  0., 0., 1., 
		  0., 0., 0.};

  // Origin
  DOMNode* cpos = ReaderXML::getNode( root, "Centre" );
  if ( cpos ) 
  {
    t[9] = ReaderXML::getNodeAttr_Double(cpos, "X");
    t[10] = ReaderXML::getNodeAttr_Double(cpos, "Y");
    t[11] = ReaderXML::getNodeAttr_Double(cpos, "Z");
  }

  // Rotation matrix
  // We only allow rotation matrices, stretching/scaling is prohibited
  DOMNode* angpos = ReaderXML::getNode( root, "AngularPosition" ); 
  assert( angpos != NULL );   
  string mode = ReaderXML::getNodeAttr_String( angpos, "Type" );
  if ( mode == "Matrix" ) 
  {
    string values = ReaderXML::getNodeValue_String( angpos );
    istringstream inValues( values.c_str() ); 
    inValues >> t[0] >> t[1] >> t[2] 
	>> t[3] >> t[4] >> t[5] 
	>> t[6] >> t[7] >> t[8];
	
    // Check that the matrix is a rotation matrix
    m_basis.setValue( t );
    if ( !m_basis.isRotation() )
    {
      cout << "A matrix in one of the AngularPosition XML nodes is"
      	<< " not a rotation matrix !!!" << endl;
      exit(1);
    }
  } 
  
  // Set the transformation
  setValue( t );
}




// --------------------------------------------------------------------------
// Copies the transformation in a 1D array
void Transform::copyTransform( double* vit, int i ) const
{
  m_basis.copyMatrix( vit, i );
  for (int j=0;j<3;++j) vit[i+9+j] = m_origin[j];
}




// --------------------------------------------------------------------------
// Copies the transformation in a 1D array composed on the left by a
// translation (useful for periodic particles in parallel)
void Transform::copyTransform( double* vit, int i, Vector3 const& vec ) const
{
  m_basis.copyMatrix( vit, i );
  for (int j=0;j<3;++j) vit[i+9+j] = m_origin[j] + vec[j];
}




// --------------------------------------------------------------------------
// Returns the matrix part of the transformation
Matrix const& Transform::getBasis() const
{
  return ( m_basis ) ;
}




// --------------------------------------------------------------------------
// Returns the transformation type
unsigned int Transform::getType() const
{
  return ( m_type ) ;
}




// --------------------------------------------------------------------------
// Returns a pointer to the origin of the transformation
Point3 const* Transform::getOrigin() const
{
  return ( &m_origin ) ;
}
