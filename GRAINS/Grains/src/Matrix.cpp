#include "Matrix.H"


// ----------------------------------------------------------------------------
// Constructeur par defaut
// D.PETIT - Octo.2000 - Creation
Matrix::Matrix()
{
  setValue( 1., 0., 0., 0., 1., 0., 0., 0., 1. ); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix::Matrix( const Scalar *m )
{ 
  setValue( m ); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix::Matrix( const Quaternion& q ) 
{ 
  setRotation( q ); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix::Matrix( Scalar x, Scalar y, Scalar z ) 
{ 
  setScaling( x, y, z ); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix::Matrix( Scalar xx, Scalar xy, Scalar xz,
	       Scalar yx, Scalar yy, Scalar yz,
	       Scalar zx, Scalar zy, Scalar zz ) 
{ 
  setValue( xx, xy, xz, yx, yy, yz, zx, zy, zz );
}




// ----------------------------------------------------------------------------
// Constructeur par copie
// A.WACHS - Aout 2009 - Creation
Matrix::Matrix( const Matrix &other )
{
  for ( int i=0;i<3;++i )
    for ( int j=0;j<3;++j )
      elem[i][j] = other.elem[i][j];
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix::~Matrix()
{
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix Matrix::absolute() const 
{
  return Matrix( fabs( elem[X][X] ), fabs( elem[X][Y] ), fabs( elem[X][Z] ),
		fabs( elem[Y][X] ), fabs( elem[Y][Y] ), fabs( elem[Y][Z] ),
		fabs( elem[Z][X] ), fabs( elem[Z][Y] ), fabs( elem[Z][Z] ) );
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix Matrix::adjoint() const 
{
  return Matrix( elem[Y][Y] * elem[Z][Z] - elem[Y][Z] * elem[Z][Y],
		elem[X][Z] * elem[Z][Y] - elem[X][Y] * elem[Z][Z],
		elem[X][Y] * elem[Y][Z] - elem[X][Z] * elem[Y][Y],
		elem[Y][Z] * elem[Z][X] - elem[Y][X] * elem[Z][Z],
		elem[X][X] * elem[Z][Z] - elem[X][Z] * elem[Z][X],
		elem[X][Z] * elem[Y][X] - elem[X][X] * elem[Y][Z],
		elem[Y][X] * elem[Z][Y] - elem[Y][Y] * elem[Z][X],
		elem[X][Y] * elem[Z][X] - elem[X][X] * elem[Z][Y],
		elem[X][X] * elem[Y][Y] - elem[X][Y] * elem[Y][X] );
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Scalar Matrix::determinant() const 
{ 
  return triple( (*this)[X], (*this)[Y], (*this)[Z] );
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Mat3& Matrix::getValue()       
{ 
return elem; 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
const Mat3& Matrix::getValue() const 
{ 
  return elem; 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix Matrix::inverse() const 
{
  Vecteur co( elem[Y][Y] * elem[Z][Z] - elem[Y][Z] * elem[Z][Y],
	    elem[Y][Z] * elem[Z][X] - elem[Y][X] * elem[Z][Z], 
	    elem[Y][X] * elem[Z][Y] - elem[Y][Y] * elem[Z][X] );
  Scalar d = (*this)[X] * co;
//  assert( !eqz( d ) ); EPSILON = 1.e-10
  Scalar s = 1 / d;
  return Matrix( co[X] * s,
		( elem[X][Z] * elem[Z][Y] - elem[X][Y] * elem[Z][Z] ) * s,
		( elem[X][Y] * elem[Y][Z] - elem[X][Z] * elem[Y][Y] ) * s,
		co[Y] * s,
		( elem[X][X] * elem[Z][Z] - elem[X][Z] * elem[Z][X] ) * s,
		( elem[X][Z] * elem[Y][X] - elem[X][X] * elem[Y][Z] ) * s,
		co[Z] * s,
		( elem[X][Y] * elem[Z][X] - elem[X][X] * elem[Z][Y] ) * s,
		( elem[X][X] * elem[Y][Y] - elem[X][Y] * elem[Y][X] ) * s );
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
void Matrix::setIdentity()
{ 
  setValue( 1., 0., 0., 0., 1., 0., 0., 0., 1. ); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
void Matrix::setRotation( const Quaternion& q ) 
{
  Scalar d = Norm2( q );
  assert( !eqz( d ) );
  Scalar s = 2. / d;
  Scalar xs = q[X] * s,   ys = q[Y] * s,   zs = q[Z] * s;
  Scalar wx = q[W] * xs,  wy = q[W] * ys,  wz = q[W] * zs;
  Scalar xx = q[X] * xs,  xy = q[X] * ys,  xz = q[X] * zs;
  Scalar yy = q[Y] * ys,  yz = q[Y] * zs,  zz = q[Z] * zs;

  setValue( 1 - ( yy + zz ),       xy - wz,       xz + wy,
	   xy + wz      , 1 - ( xx + zz ),       yz - wx,
	   xz - wy      , yz + wx      , 1 - ( xx + yy ) );
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
void Matrix::setScaling( Scalar x, Scalar y, Scalar z ) 
{
  setValue( x, 0, 0, 0, y, 0, 0, 0, z ); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
void Matrix::setValue( const Scalar *m ) 
{
  elem[X][X] = *m++; elem[Y][X] = *m++; elem[Z][X] = *m++; m++;
  elem[X][Y] = *m++; elem[Y][Y] = *m++; elem[Z][Y] = *m++; m++;
  elem[X][Z] = *m++; elem[Y][Z] = *m++; elem[Z][Z] = *m;
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
void Matrix::setValue( Scalar xx, Scalar xy, Scalar xz, 
		Scalar yx, Scalar yy, Scalar yz, 
		Scalar zx, Scalar zy, Scalar zz ) 
{
  elem[X][X] = xx; elem[X][Y] = xy; elem[X][Z] = xz;
  elem[Y][X] = yx; elem[Y][Y] = yy; elem[Y][Z] = yz;
  elem[Z][X] = zx; elem[Z][Y] = zy; elem[Z][Z] = zz;
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Scalar Matrix::tdot( int i, const Vecteur& v ) const 
{
  return elem[X][i] * v[X] + elem[Y][i] * v[Y] + elem[Z][i] * v[Z];
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix Matrix::transpose() const 
{
  return Matrix( elem[X][X], elem[Y][X], elem[Z][X],
		elem[X][Y], elem[Y][Y], elem[Z][Y],
		elem[X][Z], elem[Y][Z], elem[Z][Z] );
}




// ----------------------------------------------------------------------------
// Copie de la matrice de transformation
void Matrix::copyMatrix( double *vit, int i ) const
{
  for ( int j=0;j<3;++j ) vit[i+j]=elem[j][0];
  for ( int j=0;j<3;++j ) vit[i+4+j]=elem[j][1];
  for ( int j=0;j<3;++j ) vit[i+8+j]=elem[j][2];    
  vit[i+3]=0.; 
  vit[i+7]=0.;   
  vit[i+11]=0.;    
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Vecteur& Matrix::operator[] ( int i ) 
{ 
  return *( Vecteur * )elem[i]; 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
const Vecteur& Matrix::operator[] ( int i ) const 
{ 
  return *( Vecteur * )elem[i]; 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Vecteur operator*( const Matrix& m, const Vecteur& v ) 
{
  return Vecteur( m[X] * v, m[Y] * v, m[Z] * v );
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Vecteur operator*( const Vecteur& v, const Matrix& m ) 
{
  return Vecteur( m.tdot( X, v ), m.tdot( Y, v ), m.tdot( Z, v ) );
}




// ----------------------------------------------------------------------------
// Operation =
// A.WACHS - Aout 2009 - Creation
Matrix& Matrix::operator=( const Matrix& other )
{
  if ( &other != this )
  {
    for ( int i=0;i<3;++i )
      for ( int j=0;j<3;++j )
        elem[i][j] = other.elem[i][j];
  }

  return ( *this );
} 




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix operator*( const Matrix& m1, const Matrix& m2 ) 
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
// Multiplie à droite par une matrice de scaling
void Matrix::multiplyByScalingMatrix( Scalar x, Scalar y, Scalar z )
{
  elem[X][X] *= x ;
  elem[X][Y] *= y ; 
  elem[X][Z] *= z ;
  elem[Y][X] *= x ;
  elem[Y][Y] *= y ; 
  elem[Y][Z] *= z ;
  elem[Z][X] *= x ;
  elem[Z][Y] *= y ; 
  elem[Z][Z] *= z ;         
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix& Matrix::operator+=( const Matrix& m )
{
  setValue( elem[X][X] + m[X][X],elem[X][Y] + m[X][Y],elem[X][Z] + m[X][Z],
	   elem[Y][X] + m[Y][X],elem[Y][Y] + m[Y][Y],elem[Y][Z] + m[Y][Z],
	   elem[Z][X] + m[Z][X],elem[Z][Y] + m[Z][Y],elem[Z][Z] + m[Z][Z] );
  return *this;
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix& Matrix::operator*=( const Matrix& m ) 
{
  setValue( elem[X][X] * m[X][X] + elem[X][Y] * m[Y][X] + elem[X][Z] * m[Z][X],
	   elem[X][X] * m[X][Y] + elem[X][Y] * m[Y][Y] + elem[X][Z] * m[Z][Y],
	   elem[X][X] * m[X][Z] + elem[X][Y] * m[Y][Z] + elem[X][Z] * m[Z][Z],
	   elem[Y][X] * m[X][X] + elem[Y][Y] * m[Y][X] + elem[Y][Z] * m[Z][X],
	   elem[Y][X] * m[X][Y] + elem[Y][Y] * m[Y][Y] + elem[Y][Z] * m[Z][Y],
	   elem[Y][X] * m[X][Z] + elem[Y][Y] * m[Y][Z] + elem[Y][Z] * m[Z][Z],
	   elem[Z][X] * m[X][X] + elem[Z][Y] * m[Y][X] + elem[Z][Z] * m[Z][X],
	   elem[Z][X] * m[X][Y] + elem[Z][Y] * m[Y][Y] + elem[Z][Z] * m[Z][Y],
	   elem[Z][X] * m[X][Z] + elem[Z][Y] * m[Y][Z] + elem[Z][Z] * m[Z][Z] );
  return *this;
}


// ============================================================================


// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix multTransposeLeft( const Matrix& m1, const Matrix& m2 ) 
{
  return Matrix(
    m1[X][X] * m2[X][X] + m1[Y][X] * m2[Y][X] + m1[Z][X] * m2[Z][X],
    m1[X][X] * m2[X][Y] + m1[Y][X] * m2[Y][Y] + m1[Z][X] * m2[Z][Y],
    m1[X][X] * m2[X][Z] + m1[Y][X] * m2[Y][Z] + m1[Z][X] * m2[Z][Z],
    m1[X][Y] * m2[X][X] + m1[Y][Y] * m2[Y][X] + m1[Z][Y] * m2[Z][X],
    m1[X][Y] * m2[X][Y] + m1[Y][Y] * m2[Y][Y] + m1[Z][Y] * m2[Z][Y],
    m1[X][Y] * m2[X][Z] + m1[Y][Y] * m2[Y][Z] + m1[Z][Y] * m2[Z][Z],
    m1[X][Z] * m2[X][X] + m1[Y][Z] * m2[Y][X] + m1[Z][Z] * m2[Z][X],
    m1[X][Z] * m2[X][Y] + m1[Y][Z] * m2[Y][Y] + m1[Z][Z] * m2[Z][Y],
    m1[X][Z] * m2[X][Z] + m1[Y][Z] * m2[Y][Z] + m1[Z][Z] * m2[Z][Z] );
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Scalar determinant( const Matrix& m ) 
{ 
  return m.determinant(); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix absolute( const Matrix& m ) 
{ 
  return m.absolute(); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix transpose( const Matrix& m ) 
{
  return m.transpose(); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix adjoint( const Matrix& m ) 
{
  return m.adjoint(); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
Matrix inverse( const Matrix& m ) 
{ 
return m.inverse(); 
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
ostream& operator<<( ostream& fileOut, const Matrix& m ) 
{
  return fileOut << m[X] << m[Y] << m[Z];
}




// ----------------------------------------------------------------------
// Operateur d'ecriture
// A.WACHS - Janv.2011 - Creation.
void Matrix::writeMatrix( ostream &fileOut ) const
{
  (*this)[X].writeGroup3( fileOut );
  fileOut << endl;
  (*this)[Y].writeGroup3( fileOut );
  fileOut << endl;  
  (*this)[Z].writeGroup3( fileOut );
}




// ----------------------------------------------------------------------
// Operateur d'ecriture au format de reload 2014
// A.WACHS - Aout 2014 - Creation.
void Matrix::writeMatrix2014( ostream &fileOut ) const
{
  (*this)[X].writeGroup3(fileOut);
  fileOut << " ";
  (*this)[Y].writeGroup3(fileOut);
  fileOut << " ";  
  (*this)[Z].writeGroup3(fileOut);
}




// ----------------------------------------------------------------------
// Operateur d'ecriture en binaire au format de reload 2014
// A.WACHS - Aout 2014 - Creation.
void Matrix::writeMatrix2014_binary( ostream &fileOut )
{
  fileOut.write( reinterpret_cast<char*>( &elem[X][X] ), sizeof(double) );
  fileOut.write( reinterpret_cast<char*>( &elem[X][Y] ), sizeof(double) );  
  fileOut.write( reinterpret_cast<char*>( &elem[X][Z] ), sizeof(double) );  
  fileOut.write( reinterpret_cast<char*>( &elem[Y][X] ), sizeof(double) );
  fileOut.write( reinterpret_cast<char*>( &elem[Y][Y] ), sizeof(double) );  
  fileOut.write( reinterpret_cast<char*>( &elem[Y][Z] ), sizeof(double) );  
  fileOut.write( reinterpret_cast<char*>( &elem[Z][X] ), sizeof(double) );
  fileOut.write( reinterpret_cast<char*>( &elem[Z][Y] ), sizeof(double) );  
  fileOut.write( reinterpret_cast<char*>( &elem[Z][Z] ), sizeof(double) );       
}




// ----------------------------------------------------------------------------
//  Gino van den Bergen -  Eindhoven University of Technology - Creation 
istream &operator >> ( istream &fileIn, Matrix &m ) 
{
  return fileIn >> m[X] >> m[Y] >> m[Z];
}




// ----------------------------------------------------------------------
// Lecture en binaire de l'objet avec le format de reload 2014 en binaire
void Matrix::readMatrix2014_binary( istream &StreamIN )
{
  double *tab = new double[9] ;
  StreamIN.read( reinterpret_cast<char*>( tab ), 9 * sizeof(double) );

  elem[X][X] = tab[0];
  elem[X][Y] = tab[1]; 
  elem[X][Z] = tab[2]; 
  elem[Y][X] = tab[3]; 
  elem[Y][Y] = tab[4]; 
  elem[Y][Z] = tab[5]; 
  elem[Z][X] = tab[6]; 
  elem[Z][Y] = tab[7]; 
  elem[Z][Z] = tab[8];
  
  delete[] tab;
} 
