/*
  GJK Engine - A Fast and Robust GJK Implementation 
  Copyright (C) 1998  Gino van den Bergen

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Library General Public
  License as published by the Free Software Foundation; either
  version 2 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Library General Public License for more details.

  You should have received a copy of the GNU Library General Public
  License along with this library; if not, write to the Free
  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  Please send remarks, questions and bug reports to gino@win.tue.nl,
  or write to:
                  Gino van den Bergen
		  Department of Mathematics and Computing Science
		  Eindhoven University of Technology
		  P.O. Box 513, 5600 MB Eindhoven, The Netherlands
*/
// ============================================================================
#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "Transform.H"
#include <string>
#include <iostream>
#include <sstream>
using namespace std;


// ----------------------------------------------------------------------------
// Constructeur par defaut
Transform::Transform() :
  type(IDENTITY)
{
  basis.setIdentity();
  origin.setValue(0, 0, 0);
}




// ----------------------------------------------------------------------------
// Constructeur avec initialisation de l'origine */
Transform::Transform( const double &gx, const double &gy, const double &gz )
{
  basis.setIdentity();
  origin.setValue(gx, gy, gz);
}




// --------------------------------------------------------------------
// Gino van den Bergen - Eindhoven University of Technology - Creation 
Transform::Transform( const Scalar m[16] ) 
{ 
  setValue(m); 
}




// ----------------------------------------------------------------------------
// Constructeur par copie
// A.WACHS - Aout 2009 - Creation
Transform::Transform( const Transform &other )
{
  basis  = other.basis;
  origin = other.origin;
  type   = other.type;  
}




// ----------------------------------------------------------------------------
void Transform::invert( const Transform& t ) 
{
  basis = t.type & SCALING ? inverse(t.basis) : transpose(t.basis);
  origin.setValue(-basis[X] * t.origin, 
		  -basis[Y] * t.origin, 
		  -basis[Z] * t.origin);  
  type = t.type;
}




// ----------------------------------------------------------------------------
// La matrice est elle identite
bool Transform::isIdentity() const
{
  return (type == IDENTITY);
}




// ----------------------------------------------------------------------------
void Transform::mult( const Transform& t1, const Transform& t2 ) 
{
  basis = t1.basis * t2.basis;
  origin = t1(t2.origin);
  type = t1.type | t2.type;
}




// ----------------------------------------------------------------------------
// Composition avec une transformation de scaling 
void Transform::composeWithScaling( Scalar x, Scalar y, Scalar z )
{
  basis.multiplyByScalingMatrix( x, y, z );
} 



// ----------------------------------------------------------------------------
void Transform::multInverseLeft( const Transform& t1, const Transform& t2 ) 
{
  Vecteur v = t2.origin - t1.origin;
  if (t1.type & SCALING) {
    Matrix inv = inverse(t1.basis);
    basis = inv * t2.basis;
    origin = inv * v;
  }
  else {
    basis = multTransposeLeft(t1.basis, t2.basis);
    origin = v * t1.basis;
  }
  type = t1.type | t2.type;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// Ecriture de la Transformation
void Transform::printClass( ostream &fileOut ) const 
{
  fileOut << *this;
}




// ----------------------------------------------------------------------------
// Operateur d'ecriture.
void Transform::readClass( istream &fileIn ) 
{
  fileIn >> *this;
}




// ----------------------------------------------------------------------------
void Transform::rotate( const Quaternion& q ) 
{ 
  basis *= Matrix(q); 
  type |= ROTATION; 
}




// ----------------------------------------------------------------------------
// D.PETIT - aout 2000 - Creation
void Transform::rotateOrientation( const Quaternion& q ) 
{
//   Matrix m = basis.transpose();
//   m[X].Tourne(q);
//   m[Y].Tourne(q);
//   m[Z].Tourne(q);
//   basis = m.transpose();  

  // Nouvelle methode plus rapide en developpant l'expression q.basis.qt
  const Mat3& mm = basis.getValue(); 
  double qx = q[X], qy = q[Y], qz = q[Z], qw = q[W];
  double px, py, pz, pw;
  for (int i=0;i<3;++i)
  {
    px = qy *  mm[Z][i] - qz * mm[Y][i] + qw * mm[X][i];
    py = qz *  mm[X][i] - qx * mm[Z][i] + qw * mm[Y][i];    
    pz = qx *  mm[Y][i] - qy * mm[X][i] + qw * mm[Z][i]; 
    pw = - qx * mm[X][i] - qy * mm[Y][i] - qz * mm[Z][i];
    basis[X][i] = qy * pz - qz * py - pw * qx + qw * px;
    basis[Y][i] = qz * px - qx * pz - pw * qy + qw * py;          
    basis[Z][i] = qx * py - qy * px - pw * qz + qw * pz;    
  } 
}




// ----------------------------------------------------------------------------
void Transform::scale( Scalar x, Scalar y, Scalar z ) 
{ 
  basis *= Matrix(x, y, z);  
  type |= SCALING;
}




// ----------------------------------------------------------------------------
void Transform::setIdentity() 
{
  basis.setIdentity();
  origin.setValue(0, 0, 0);
  type = IDENTITY;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil 2000 - Creation
void Transform::setOrigin( const Scalar* pos ) 
{
  origin[X] = pos[X];
  origin[Y] = pos[Y];
  origin[Z] = pos[Z];
}




// ----------------------------------------------------------------------------
// A.WACHS - Mai 2009 - Creation
void Transform::setOrigin( const double &gx, const double &gy, 
	const double &gz ) 
{
  origin[X] = gx;
  origin[Y] = gy;
  origin[Z] = gz;
}




// ----------------------------------------------------------------------------
// G.FERRER - Juin.2003 - Creation
void Transform::setOrigin( const Point &pos ) 
{
  origin = pos;
}




// ----------------------------------------------------------------------------
// Modifie la matrice de transformation
void Transform::setBasis( const Matrix &basis_ )
{
  basis = basis_;
  type = AFFINE;  
}




// ----------------------------------------------------------------------------
void Transform::setValue( const Scalar m[16] ) 
{
  basis.setValue(m);
  origin.setValue(&m[12]);
  type = AFFINE;
}




// ----------------------------------------------------------------------------
void Transform::translate( const Vecteur& v )
{ 
  origin += v;
  type |= TRANSLATION;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// Passage du point du repere local au repre global
Point Transform::operator() ( const Point& p ) const 
{
  return Point(basis[X] * p + origin[X], 
	       basis[Y] * p + origin[Y], 
	       basis[Z] * p + origin[Z]);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Operateur d'egalite : affecte la transformation donnee a *this
Transform& Transform::operator= ( const Transform& transform )
{
  if ( &transform != this )
  {  
    basis  = transform.basis;
    origin = transform.origin;
    type   = transform.type;
  }
  
  return *this;
}




// ----------------------------------------------------------------------------
// D.RAKOTONIRINA - Fev.2014 - Creation
// Composition a gauche par une autre transformation affine
// !!! Resultat: d'abord t puis this !!! 
void Transform::composeTransformLeft( const Transform& t ) 
{
  origin += basis * t.origin;
  basis *= t.basis;
  type |= t.type; 
}




// ----------------------------------------------------------------------------
// D.RAKOTONIRINA - Fev.2014 - Creation
// Composition a droite par une autre transformation affine
// !!! Resultat: d'abord this puis t !!! 
void Transform::composeTransformRight( const Transform& t ) 
{
   origin = t.origin + t.basis * origin;
   basis = t.basis * basis ;
   type = t.type | type ;      
}




// ----------------------------------------------------------------------------
// D.RAKOTONIRINA - Fev.2014 - Creation
// Composition a droite par une rotation rot par rapport a l'origine 
// de this
// !!! Resultat: d'abord this puis r !!! 
// !!! Cette composition ne change pas l'origine de this !!! 
void Transform::composeRotationRight( const Transform& t ) 
{
   basis = t.basis * basis ;
   type = t.type | type ;      
}




// ----------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// Operateur d'ecriture
// A.WACHS - Mars.2010 - Modif.
ostream &operator << ( ostream &fileOut, const  Transform &t ) 
{
  fileOut << "%Position&Orientation\n";
  fileOut << "*Position\n";
  fileOut << t.origin;
  fileOut << "*Orientation\n";
  fileOut << t.basis;
  fileOut << "*Type\n";
  fileOut << t.type << endl;
  return fileOut;
}




// ----------------------------------------------------------------------
// Operateur d'ecriture
// A.WACHS - Janv.2011 - Creation.
void Transform::writeTransform( ostream &fileOut ) const
{
  fileOut << "%Position&Orientation\n";
  fileOut << "*Position\n";
  fileOut << Grains_Exec::doubleToString(ios::scientific,POSITIONFORMAT,
  	origin[X]) << " " << 
	Grains_Exec::doubleToString(ios::scientific,POSITIONFORMAT,
  	origin[Y]) << " " << 
	Grains_Exec::doubleToString(ios::scientific,POSITIONFORMAT,
  	origin[Z]) << endl;	
  fileOut << "*Orientation\n";
  basis.writeMatrix( fileOut );
  fileOut <<  endl;
  fileOut << "*Type\n";
  fileOut << type;
}




// ----------------------------------------------------------------------
// Operateur d'ecriture avec le format de reload 2014
// A.WACHS - Aout 2014 - Creation.
void Transform::writeTransform2014( ostream &fileOut ) const
{
  fileOut << Grains_Exec::doubleToString(ios::scientific,POSITIONFORMAT,
  	origin[X]) << " " << 
	Grains_Exec::doubleToString(ios::scientific,POSITIONFORMAT,
  	origin[Y]) << " " << 
	Grains_Exec::doubleToString(ios::scientific,POSITIONFORMAT,
  	origin[Z]) << " ";	
  basis.writeMatrix2014( fileOut );
  fileOut << " " << type;
}




// ----------------------------------------------------------------------
// Operateur d'ecriture en binaire avec le format de reload 2014
// A.WACHS - Aout 2014 - Creation.
void Transform::writeTransform2014_binary( ostream &fileOut )
{
  fileOut.write( reinterpret_cast<char*>( &origin[X] ), sizeof(double) );
  fileOut.write( reinterpret_cast<char*>( &origin[Y] ), sizeof(double) );  
  fileOut.write( reinterpret_cast<char*>( &origin[Z] ), sizeof(double) );
  basis.writeMatrix2014_binary( fileOut );
  fileOut.write( reinterpret_cast<char*>( &type ), sizeof(unsigned int) );
}




// ----------------------------------------------------------------------
// D.PETIT - Juil. 2000 - Creation
// Operateur de lecture.
istream &operator >> ( istream &fileIn, Transform &t ) 
{
  string cle;
  fileIn >> cle >> t.origin;
  fileIn >> cle >> t.basis;
  fileIn >> cle >> t.type;
  return fileIn;
}




// ----------------------------------------------------------------------
// Lecture de l'objet sur le flux d'entre avec le format de reload 2014
// A.WACHS - Aout 2014 - Creation
void Transform::readTransform2014( istream &StreamIN )
{
  StreamIN >> origin >> basis >> type;
}




// ----------------------------------------------------------------------
// Lecture de l'objet sur le flux d'entre avec le format de reload 2014
// en binaire
// A.WACHS - Dec 2014 - Creation
void Transform::readTransform2014_binary( istream &StreamIN )
{
  StreamIN.read( reinterpret_cast<char*>( &origin[X] ), sizeof(double) );
  StreamIN.read( reinterpret_cast<char*>( &origin[Y] ), sizeof(double) );  
  StreamIN.read( reinterpret_cast<char*>( &origin[Z] ), sizeof(double) );
  basis.readMatrix2014_binary( StreamIN );
  StreamIN.read( reinterpret_cast<char*>( &type ), sizeof(unsigned int) );      
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Persistance du Transform
void Transform::load( DOMNode* root )
{
  DOMNode* orientation = ReaderXML::getNode(root, "Orientation");

  assert(orientation != NULL);

  // valeurs initiales de la transformation : bizarrement 16 dans Transform
  double t[16] = {1., 0., 0., 0.,
		  0., 1., 0., 0.,
		  0., 0., 1., 0.,
		  0., 0., 0., 0.};

  DOMNode* cpos = ReaderXML::getNode(root, "Centre");
  if (cpos) {
    t[12] = ReaderXML::getNodeAttr_Double(cpos, "X");
    t[13] = ReaderXML::getNodeAttr_Double(cpos, "Y");
    t[14] = ReaderXML::getNodeAttr_Double(cpos, "Z");
  }

  string mode = ReaderXML::getNodeAttr_String(orientation, "Type");
  if (mode == "Matrice") {
    string values = ReaderXML::getNodeValue_String(orientation);
    istringstream inValues(values.c_str()); // DGR: position.setValue surprenant
    inValues >> t[0] >> t[4] >> t[8] 
	     >> t[1] >> t[5] >> t[9] 
	     >> t[2] >> t[6] >> t[10];
    setValue(t);  

  } else {
    setValue(t);
    type = IDENTITY;
  }
}




// ----------------------------------------------------------------------------
// CopyTransform (poly?dres)
void Transform::copyTransform( double *vit, int i ) const
{
  basis.copyMatrix(vit,i);
  for (int j=0;j<3;++j) vit[i+12+j]=origin[j];
  vit[i+15]=0.;  
}




// ----------------------------------------------------------------------------
// CopyTransform (poly?dres)
void Transform::copyTransform( double *vit, int i, Vecteur const& vec ) const
{
  basis.copyMatrix(vit,i);
  for (int j=0;j<3;++j) vit[i+12+j]=origin[j]+vec[j];
  vit[i+15]=0.;  
}
