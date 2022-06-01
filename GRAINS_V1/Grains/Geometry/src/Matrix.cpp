#include "Matrix.hh"


// ----------------------------------------------------------------------------
// Default constructor. Matrix is initialized to the identity matrix
Matrix::Matrix()
{
  setValue( 1., 0., 0., 0., 1., 0., 0., 0., 1. ); 
}




// ----------------------------------------------------------------------------
// Constructor with an 1D array of values as inputs
Matrix::Matrix( double const* m )
{ 
  setValue( m ); 
}




// ----------------------------------------------------------------------------
// Constructor with a quaternion as input to initialized the matrix
// as a rotation matrix. The quaternion is not required to be unitary as it
// will be normalized by dividing by its norm.
Matrix::Matrix( Quaternion const& q ) 
{ 
  setRotation( q ); 
}




// ----------------------------------------------------------------------------
// Constructor of a diagonal matrix with the 3 diagonal coefficients as inputs
Matrix::Matrix( double x, double y, double z ) 
{ 
  setScaling( x, y, z ); 
}




// ----------------------------------------------------------------------------
// Constructor with all the 9 coefficients as inputs
Matrix::Matrix( double xx, double xy, double xz,
	       double yx, double yy, double yz,
	       double zx, double zy, double zz ) 
{ 
  setValue( xx, xy, xz, yx, yy, yz, zx, zy, zz );
}




// ----------------------------------------------------------------------------
// Copy constructor
Matrix::Matrix( Matrix const& other )
{
  for ( int i=0;i<3;++i )
    for ( int j=0;j<3;++j )
      m_elem[i][j] = other.m_elem[i][j];
}




// ----------------------------------------------------------------------------
// Destructor
Matrix::~Matrix()
{}




// ----------------------------------------------------------------------------
// Returns a matrix with all coefficients have the absolute value of
// the coefficients of this
Matrix Matrix::absolute() const 
{
  return ( Matrix( 
  	fabs( m_elem[X][X] ), fabs( m_elem[X][Y] ), fabs( m_elem[X][Z] ),
	fabs( m_elem[Y][X] ), fabs( m_elem[Y][Y] ), fabs( m_elem[Y][Z] ),
	fabs( m_elem[Z][X] ), fabs( m_elem[Z][Y] ), fabs( m_elem[Z][Z] ) ) );
}




// ----------------------------------------------------------------------------
// Returns the adjoint matrix of this. The adjoint matrix is the
// transpose of the cofactor matrix
Matrix Matrix::adjoint() const 
{
  return ( Matrix( 
  	m_elem[Y][Y] * m_elem[Z][Z] - m_elem[Y][Z] * m_elem[Z][Y],
	m_elem[X][Z] * m_elem[Z][Y] - m_elem[X][Y] * m_elem[Z][Z],
	m_elem[X][Y] * m_elem[Y][Z] - m_elem[X][Z] * m_elem[Y][Y],
	m_elem[Y][Z] * m_elem[Z][X] - m_elem[Y][X] * m_elem[Z][Z],
	m_elem[X][X] * m_elem[Z][Z] - m_elem[X][Z] * m_elem[Z][X],
	m_elem[X][Z] * m_elem[Y][X] - m_elem[X][X] * m_elem[Y][Z],
	m_elem[Y][X] * m_elem[Z][Y] - m_elem[Y][Y] * m_elem[Z][X],
	m_elem[X][Y] * m_elem[Z][X] - m_elem[X][X] * m_elem[Z][Y],
	m_elem[X][X] * m_elem[Y][Y] - m_elem[X][Y] * m_elem[Y][X] ) );
}




// ----------------------------------------------------------------------------
// Returns the determinant of the matrix
double Matrix::determinant() const 
{ 
  return ( triple( (*this)[X], (*this)[Y], (*this)[Z] ) );
}




// ----------------------------------------------------------------------------
// Returns a 3x3 array describing the matrix
Mat3& Matrix::getValue()       
{ 
  return ( m_elem ); 
}




// ----------------------------------------------------------------------------
// Returns a 3x3 array describing the matrix
Mat3 const& Matrix::getValue() const 
{ 
  return ( m_elem ); 
}




// ----------------------------------------------------------------------------
// Returns the inverse matrix
Matrix Matrix::inverse() const 
{
  Vector3 co( m_elem[Y][Y] * m_elem[Z][Z] - m_elem[Y][Z] * m_elem[Z][Y],
	m_elem[Y][Z] * m_elem[Z][X] - m_elem[Y][X] * m_elem[Z][Z], 
	m_elem[Y][X] * m_elem[Z][Y] - m_elem[Y][Y] * m_elem[Z][X] );
  double d = (*this)[X] * co;
//  assert( !eqz( d ) ); EPSILON = 1.e-10
  double s = 1 / d;
  return ( Matrix( co[X] * s,
	( m_elem[X][Z] * m_elem[Z][Y] - m_elem[X][Y] * m_elem[Z][Z] ) * s,
	( m_elem[X][Y] * m_elem[Y][Z] - m_elem[X][Z] * m_elem[Y][Y] ) * s,
	co[Y] * s,
	( m_elem[X][X] * m_elem[Z][Z] - m_elem[X][Z] * m_elem[Z][X] ) * s,
	( m_elem[X][Z] * m_elem[Y][X] - m_elem[X][X] * m_elem[Y][Z] ) * s,
	co[Z] * s,
	( m_elem[X][Y] * m_elem[Z][X] - m_elem[X][X] * m_elem[Z][Y] ) * s,
	( m_elem[X][X] * m_elem[Y][Y] - m_elem[X][Y] * m_elem[Y][X] ) * s ) );
}




// ----------------------------------------------------------------------------
// Sets the matrix to the identity matrix
void Matrix::setIdentity()
{ 
  setValue( 1., 0., 0., 0., 1., 0., 0., 0., 1. ); 
}




// ----------------------------------------------------------------------------
// Sets the matrix to a rotation matrix with a quaternion. The 
// quaternion is not required to be unitary as it will be normalized by 
// dividing by its norm.
void Matrix::setRotation( Quaternion const& q ) 
{
  double d = Norm2( q );
  assert( !eqz( d ) );
  double s = 2. / d;
  double xs = q[X] * s,   ys = q[Y] * s,   zs = q[Z] * s;
  double wx = q[W] * xs,  wy = q[W] * ys,  wz = q[W] * zs;
  double xx = q[X] * xs,  xy = q[X] * ys,  xz = q[X] * zs;
  double yy = q[Y] * ys,  yz = q[Y] * zs,  zz = q[Z] * zs;

  setValue( 1. - ( yy + zz ), xy - wz         , xz + wy,
	   xy + wz          , 1. - ( xx + zz ), yz - wx,
	   xz - wy          , yz + wx         , 1. - ( xx + yy ) );
}




// ----------------------------------------------------------------------------
// Sets the matrix to a diagonal matrix
void Matrix::setScaling( double x, double y, double z ) 
{
  setValue( x, 0, 0, 0, y, 0, 0, 0, z ); 
}




// ----------------------------------------------------------------------------
// Sets the matrix with an 1D array of 9 values as inputs
void Matrix::setValue( double const* m ) 
{
  m_elem[X][X] = *m++; m_elem[X][Y] = *m++; m_elem[X][Z] = *m++; 
  m_elem[Y][X] = *m++; m_elem[Y][Y] = *m++; m_elem[Y][Z] = *m++; 
  m_elem[Z][X] = *m++; m_elem[Z][Y] = *m++; m_elem[Z][Z] = *m;
}




// ----------------------------------------------------------------------------
// Sets the matrix with all the 9 coefficients as inputs
void Matrix::setValue( double xx, double xy, double xz, 
		double yx, double yy, double yz, 
		double zx, double zy, double zz ) 
{
  m_elem[X][X] = xx; m_elem[X][Y] = xy; m_elem[X][Z] = xz;
  m_elem[Y][X] = yx; m_elem[Y][Y] = yy; m_elem[Y][Z] = yz;
  m_elem[Z][X] = zx; m_elem[Z][Y] = zy; m_elem[Z][Z] = zz;
}




// ----------------------------------------------------------------------------
// Returns the scalar product of a column of the matrix and a vector
double Matrix::tdot( int i, Vector3 const& v ) const 
{
  return ( m_elem[X][i] * v[X] + m_elem[Y][i] * v[Y] + m_elem[Z][i] * v[Z] );
}




// ----------------------------------------------------------------------------
// Returns the transposed matrix
Matrix Matrix::transpose() const 
{
  return ( Matrix( m_elem[X][X], m_elem[Y][X], m_elem[Z][X],
	m_elem[X][Y], m_elem[Y][Y], m_elem[Z][Y],
	m_elem[X][Z], m_elem[Y][Z], m_elem[Z][Z] ) );
}




// ----------------------------------------------------------------------------
// Copies the matrix in a 1D array
void Matrix::copyMatrix( double *vit, int i ) const
{
  for ( int j=0;j<3;++j ) 
  {
    vit[i+j] = m_elem[X][j];
    vit[i+3+j] = m_elem[Y][j];
    vit[i+6+j] = m_elem[Z][j]; 
  }       
}




// ----------------------------------------------------------------------------
// i-th row accessor
Vector3& Matrix::operator [] ( int i ) 
{ 
  return ( *( Vector3* )m_elem[i] ); 
}




// ----------------------------------------------------------------------------
// i-th row accessor
Vector3 const& Matrix::operator [] ( int i ) const 
{ 
  return ( *( Vector3* )m_elem[i] ); 
}




// ----------------------------------------------------------------------------
// Matrix-vector product
Vector3 operator * ( Matrix const& m, Vector3 const& v ) 
{
  return ( Vector3( m[X] * v, m[Y] * v, m[Z] * v ) );
}




// ----------------------------------------------------------------------------
// Vector-matrix product
Vector3 operator*( Vector3 const& v, Matrix const& m ) 
{
  return ( Vector3( m.tdot( X, v ), m.tdot( Y, v ), m.tdot( Z, v ) ) );
}




// ----------------------------------------------------------------------------
// Equal operator to another matrix
Matrix& Matrix::operator=( Matrix const& m )
{
  if ( &m != this )
  {
    for ( int i=0;i<3;++i )
      for ( int j=0;j<3;++j )
        m_elem[i][j] = m.m_elem[i][j];
  }

  return ( *this );
} 




// ----------------------------------------------------------------------------
// Matrix-matrix product
Matrix operator * ( Matrix const& m1, Matrix const& m2 ) 
{
  return Matrix(
    m1[X][X] * m2[X][X] + m1[X][Y] * m2[Y][X] + m1[X][Z] * m2[Z][X],
    m1[X][X] * m2[X][Y] + m1[X][Y] * m2[Y][Y] + m1[X][Z] * m2[Z][Y],
    m1[X][X] * m2[X][Z] + m1[X][Y] * m2[Y][Z] + m1[X][Z] * m2[Z][Z],
    m1[Y][X] * m2[X][X] + m1[Y][Y] * m2[Y][X] + m1[Y][Z] * m2[Z][X],
    m1[Y][X] * m2[X][Y] + m1[Y][Y] * m2[Y][Y] + m1[Y][Z] * m2[Z][Y],
    m1[Y][X] * m2[X][Z] + m1[Y][Y] * m2[Y][Z] + m1[Y][Z] * m2[Z][Z],
    m1[Z][X] * m2[X][X] + m1[Z][Y] * m2[Y][X] + m1[Z][Z] * m2[Z][X],
    m1[Z][X] * m2[X][Y] + m1[Z][Y] * m2[Y][Y] + m1[Z][Z] * m2[Z][Y],
    m1[Z][X] * m2[X][Z] + m1[Z][Y] * m2[Y][Z] + m1[Z][Z] * m2[Z][Z] );
}




// ----------------------------------------------------------------------------
// Multiplies to the right by a scaling matrix defined by its 3
// diagonal coefficients
void Matrix::multiplyByScalingMatrix( double x, double y, double z )
{
  m_elem[X][X] *= x ;
  m_elem[X][Y] *= y ; 
  m_elem[X][Z] *= z ;
  m_elem[Y][X] *= x ;
  m_elem[Y][Y] *= y ; 
  m_elem[Y][Z] *= z ;
  m_elem[Z][X] *= x ;
  m_elem[Z][Y] *= y ; 
  m_elem[Z][Z] *= z ;         
}




// ----------------------------------------------------------------------------
// Operator +=
Matrix& Matrix::operator += ( const Matrix& m )
{
  setValue( 
  	m_elem[X][X] + m[X][X], m_elem[X][Y] + m[X][Y], m_elem[X][Z] + m[X][Z],
	m_elem[Y][X] + m[Y][X], m_elem[Y][Y] + m[Y][Y], m_elem[Y][Z] + m[Y][Z],
	m_elem[Z][X] + m[Z][X], m_elem[Z][Y] + m[Z][Y], m_elem[Z][Z] + m[Z][Z] 
	);
  return ( *this );
}




// ----------------------------------------------------------------------------
// Operator *=
Matrix& Matrix::operator *= ( Matrix const& m ) 
{
  setValue( 
    m_elem[X][X] * m[X][X] + m_elem[X][Y] * m[Y][X] + m_elem[X][Z] * m[Z][X],
    m_elem[X][X] * m[X][Y] + m_elem[X][Y] * m[Y][Y] + m_elem[X][Z] * m[Z][Y],
    m_elem[X][X] * m[X][Z] + m_elem[X][Y] * m[Y][Z] + m_elem[X][Z] * m[Z][Z],
    m_elem[Y][X] * m[X][X] + m_elem[Y][Y] * m[Y][X] + m_elem[Y][Z] * m[Z][X],
    m_elem[Y][X] * m[X][Y] + m_elem[Y][Y] * m[Y][Y] + m_elem[Y][Z] * m[Z][Y],
    m_elem[Y][X] * m[X][Z] + m_elem[Y][Y] * m[Y][Z] + m_elem[Y][Z] * m[Z][Z],
    m_elem[Z][X] * m[X][X] + m_elem[Z][Y] * m[Y][X] + m_elem[Z][Z] * m[Z][X],
    m_elem[Z][X] * m[X][Y] + m_elem[Z][Y] * m[Y][Y] + m_elem[Z][Z] * m[Z][Y],
    m_elem[Z][X] * m[X][Z] + m_elem[Z][Y] * m[Y][Z] + m_elem[Z][Z] * m[Z][Z] );
  return ( *this );
}




// ----------------------------------------------------------------------------
// Transposed matrix-matrix product 
Matrix multTransposeLeft( Matrix const& m1, Matrix const& m2 ) 
{
  return ( Matrix(
    m1[X][X] * m2[X][X] + m1[Y][X] * m2[Y][X] + m1[Z][X] * m2[Z][X],
    m1[X][X] * m2[X][Y] + m1[Y][X] * m2[Y][Y] + m1[Z][X] * m2[Z][Y],
    m1[X][X] * m2[X][Z] + m1[Y][X] * m2[Y][Z] + m1[Z][X] * m2[Z][Z],
    m1[X][Y] * m2[X][X] + m1[Y][Y] * m2[Y][X] + m1[Z][Y] * m2[Z][X],
    m1[X][Y] * m2[X][Y] + m1[Y][Y] * m2[Y][Y] + m1[Z][Y] * m2[Z][Y],
    m1[X][Y] * m2[X][Z] + m1[Y][Y] * m2[Y][Z] + m1[Z][Y] * m2[Z][Z],
    m1[X][Z] * m2[X][X] + m1[Y][Z] * m2[Y][X] + m1[Z][Z] * m2[Z][X],
    m1[X][Z] * m2[X][Y] + m1[Y][Z] * m2[Y][Y] + m1[Z][Z] * m2[Z][Y],
    m1[X][Z] * m2[X][Z] + m1[Y][Z] * m2[Y][Z] + m1[Z][Z] * m2[Z][Z] ) );
}




// ----------------------------------------------------------------------------
// Returns the matrix determinant
double determinant( Matrix const& m ) 
{ 
  return ( m.determinant() ); 
}




// ----------------------------------------------------------------------------
// Returns a matrix with all coefficients have the absolute value of 
// the coefficients of the matrix
Matrix absolute( Matrix const& m ) 
{ 
  return ( m.absolute() ); 
}




// ----------------------------------------------------------------------------
// Returns the transposed matrix
Matrix transpose( Matrix const& m ) 
{
  return ( m.transpose() ); 
}




// ----------------------------------------------------------------------------
// Returns the adjoint matrix 
Matrix adjoint( Matrix const& m ) 
{
  return ( m.adjoint() ); 
}




// ----------------------------------------------------------------------------
// Returns the inverse matrix
Matrix inverse( Matrix const& m ) 
{ 
  return ( m.inverse() ); 
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, Matrix const& m )
{
  return ( fileOut << m[X] << endl << m[Y] << endl << m[Z] );
}




// ----------------------------------------------------------------------
// Writes the object with a high precision format given by
// POSITIONFORMAT defined in GrainsExec.hh
void Matrix::writeMatrix( ostream& fileOut ) const
{
  (*this)[X].writeGroup3( fileOut );
  fileOut << endl;
  (*this)[Y].writeGroup3( fileOut );
  fileOut << endl;  
  (*this)[Z].writeGroup3( fileOut );
}




// ----------------------------------------------------------------------
// Writes the object with a high precision format given by
// POSITIONFORMAT defined in GrainsExec.hh and the 2014 reload format
void Matrix::writeMatrix2014( ostream& fileOut ) const
{
  (*this)[X].writeGroup3( fileOut );
  fileOut << " ";
  (*this)[Y].writeGroup3( fileOut );
  fileOut << " ";  
  (*this)[Z].writeGroup3( fileOut );
}




// ----------------------------------------------------------------------
// Writes the object in binary format with the 2014 reload format
void Matrix::writeMatrix2014_binary( ostream& fileOut )
{
  fileOut.write( reinterpret_cast<char*>( &m_elem[X][X] ), sizeof( double ) );
  fileOut.write( reinterpret_cast<char*>( &m_elem[X][Y] ), sizeof( double ) );  
  fileOut.write( reinterpret_cast<char*>( &m_elem[X][Z] ), sizeof( double ) );  
  fileOut.write( reinterpret_cast<char*>( &m_elem[Y][X] ), sizeof( double ) );
  fileOut.write( reinterpret_cast<char*>( &m_elem[Y][Y] ), sizeof( double ) );  
  fileOut.write( reinterpret_cast<char*>( &m_elem[Y][Z] ), sizeof( double ) );  
  fileOut.write( reinterpret_cast<char*>( &m_elem[Z][X] ), sizeof( double ) );
  fileOut.write( reinterpret_cast<char*>( &m_elem[Z][Y] ), sizeof( double ) );  
  fileOut.write( reinterpret_cast<char*>( &m_elem[Z][Z] ), sizeof( double ) );
}




// -----------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, Matrix& m )
{
  return ( fileIn >> m[X] >> m[Y] >> m[Z] );
}




// ----------------------------------------------------------------------
// Reads the object in binary format with the 2014 reload format
void Matrix::readMatrix2014_binary( istream& StreamIN )
{
  double* tab = new double[9] ;
  StreamIN.read( reinterpret_cast<char*>( tab ), 9 * sizeof( double ) );

  m_elem[X][X] = tab[0];
  m_elem[X][Y] = tab[1]; 
  m_elem[X][Z] = tab[2]; 
  m_elem[Y][X] = tab[3]; 
  m_elem[Y][Y] = tab[4]; 
  m_elem[Y][Z] = tab[5]; 
  m_elem[Z][X] = tab[6]; 
  m_elem[Z][Y] = tab[7]; 
  m_elem[Z][Z] = tab[8];
   
  delete[] tab;
} 




// -----------------------------------------------------------------------
// Returns whether the matrix is diagonal wrt EPSILON2 defined in Basic.hh
bool Matrix::isDiagonal() const
{
  return ( fabs( m_elem[X][Y] ) > EPSILON2 || fabs( m_elem[X][Z] ) > EPSILON2
	|| fabs( m_elem[Y][X] ) > EPSILON2 || fabs( m_elem[Y][Z] ) > EPSILON2  
	|| fabs( m_elem[Z][X] ) > EPSILON2 || fabs( m_elem[Z][Y] ) > EPSILON2 ?
	false : true );      
}




// -----------------------------------------------------------------------
// Returns whether the matrix is the identity matrix wrt tol. Default value 
// of tol is EPSILON2 defined in Basic.hh
bool Matrix::isIdentity( double tol ) const
{
  return ( fabs( m_elem[X][Y] ) < tol && fabs( m_elem[X][Z] ) < tol
	&& fabs( m_elem[Y][X] ) < tol && fabs( m_elem[Y][Z] ) < tol  
	&& fabs( m_elem[Z][X] ) < tol && fabs( m_elem[Z][Y] ) < tol 
	&& ( fabs( m_elem[X][X] ) - 1. ) < tol
	&& ( fabs( m_elem[Y][Y] ) - 1. ) < tol
	&& ( fabs( m_elem[Z][Z] ) - 1. ) < tol ? true : false );  
}




// -----------------------------------------------------------------------
// Returns whether the matrix is a rotation matrix wrt tol. Default value 
// of tol is EPSILON defined in Basic.hh
bool Matrix::isRotation( double tol ) const
{
  bool isRota = true;

  // Check that MMt = I
  Matrix MMt = multTransposeLeft( *this, *this );
  if ( !MMt.isIdentity( tol ) ) isRota = false;
  
  // Check that determinant is 1
  if ( isRota )
    if ( fabs( this->determinant() - 1. ) > tol ) isRota = false;
    
  return ( isRota );
}
