#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "Quaternion.H"
#include <math.h>


// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Constructeur par defaut
Quaternion::Quaternion() :
  w(1.), vqt(0.)
{}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Constructeur avec initialisation
Quaternion::Quaternion( const Scalar q, const Scalar d ) :
  w(d), vqt(q)
{}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Constructeur avec initialisation
Quaternion::Quaternion( const Vecteur &vecteur, const Scalar d ) :
  w(d), vqt(vecteur)
{}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Constructeur avec initialisation
Quaternion::Quaternion( Scalar x, Scalar y, Scalar z, Scalar d ) :
  w(d), vqt(Vecteur(x,y,z))
{}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Constructeur par copie
Quaternion::Quaternion( const Quaternion& q ) :
  w(q.w), vqt(q.vqt)
{}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Destructeur
Quaternion::~Quaternion()
{}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// creation du conjugue du quaternion
Quaternion Quaternion::Conjugate() const
{
  return (Quaternion(-vqt,w));
}




// ----------------------------------------------------------------------------
// Acces a la composante scalaire
// F.PRADEL - Janv.2000 - Creation
Scalar Quaternion::getScalaire() const
{
  return (w);
}




// ----------------------------------------------------------------------------
// Acces a la composante vecteur
// G.FERRER - Juil.2000 - Creation
Vecteur const* Quaternion::getVecteur() const
{
  return &vqt;
}




// ----------------------------------------------------------------------------
// Renvoie le vecteur rotation associe qui est renvoye
// F.PRADEL - Juin.2001 - Creation
Vecteur Quaternion::getVecteurRotation() const
{
  Scalar norm=Norm(vqt);
  Vecteur vTmp;
  if (norm>EPS) {
    vTmp=(2.0*atan(norm/w)/norm)*vqt;
  } else {
    vTmp=Vecteur(0.0,0.0,0.0);
  }
  return vTmp;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// creation de l'inverse du quaternion
Quaternion Quaternion::Inverse() const
{
  return (Conjugate() * (1./ABS(*this)));
}




// ----------------------------------------------------------------------------
// Modification du quaternion
// D.PETIT - Aout.2000 - Creation
void Quaternion::setQuaternion( const Vecteur &vecteur, Scalar scalaire )
{
  w = scalaire;
  vqt = vecteur;
}




// ----------------------------------------------------------------------------
// Modification du quaternion a partir d'un vecteur rotation
// F.PRADEL - Juin 2001 - Creation
void Quaternion::setQuaternion( const Vecteur &vecteur )
{
  Scalar phi   = Norm(vecteur);
  Scalar qstar = sin (phi/2.0);
  w = cos (phi/2.0);
  if (qstar > EPS) {
    qstar /= phi;
  } else {
    qstar = 0.5 * (1.0 - phi*phi/6.0 + phi*phi*phi*phi/120.0);
  }
  vqt = qstar * vecteur;
}




// ----------------------------------------------------------------------------
// Modification du quaternion
// D.PETIT - Sept.2000 - Creation
void Quaternion::setQuaternion( const Scalar &vecteur0,
	const Scalar &vecteur1,
	const Scalar &vecteur2,
	Scalar scalaire )
{
  vqt[X] = vecteur0;
  vqt[Y] = vecteur1;
  vqt[Z] = vecteur2;
  w = scalaire;
}




// ----------------------------------------------------------------------------
// Modification de la partie scalaire
// D.PETIT - Aout.2000 - Creation
void Quaternion::setScalaire( Scalar scalaire )
{
  w = scalaire;
}




// ----------------------------------------------------------------------------
// Modification de la partie vecteur
// D.PETIT - Aout.2000 - Creation
void Quaternion::setVecteur( const Vecteur &vecteur )
{
  vqt = vecteur;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Affectation
Quaternion& Quaternion::operator = ( const Quaternion& rhs )
{
  if ( &rhs != this )
  {
    w   = rhs.w;
    vqt = rhs.vqt;
  }
  return *this;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// Egalisation affectation a un scalaire
Quaternion Quaternion::operator = ( const Scalar rhs )
{
  w   = rhs;
  vqt = rhs;
  return *this;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// A.WACHS - Sept 2009 - Modif
// indexation des composantes des points
Scalar& Quaternion::operator [] ( int i )
{
  return (i==3 ? w : vqt[i]);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// A.WACHS - Sept 2009 - Modif
// indexation des composantes des points
Scalar Quaternion::operator [] ( int i ) const
{
  return (i==3 ? w :  vqt[i]);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
Quaternion Quaternion::operator - ()
{
  return Quaternion(-vqt,-w);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
Quaternion& Quaternion::operator += ( const Quaternion& rhs )
{
  w   += rhs.w;
  vqt += rhs.vqt;
  return *this;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
Quaternion Quaternion::operator + ( const Quaternion& rhs ) const
{
  return Quaternion(vqt+rhs.vqt, w+rhs.w);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
Quaternion& Quaternion::operator -= ( const Quaternion& rhs )
{
  w   -= rhs.w;
  vqt -= rhs.vqt;
  return *this;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
Quaternion Quaternion::operator - ( const Quaternion& rhs )
{
  return Quaternion(vqt-rhs.vqt, w-rhs.w);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
Quaternion& Quaternion::operator *= ( Scalar d )
{
  w   *= d;
  vqt *= d;
  return *this;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// multiplication par un scalaire
Quaternion Quaternion::operator * ( Scalar d )
{
  return Quaternion(d*vqt, d*w);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
Quaternion operator * ( Scalar d, const Quaternion& rhs )
{
  Quaternion result(rhs);
  result*=d;
  return result;
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
Quaternion Quaternion::operator*( const Quaternion& rhs ) const
{
  Scalar tmp = (w*rhs.w) - (vqt*rhs.vqt);
  Vecteur vtmp = (vqt^rhs.vqt) + (w*rhs.vqt) + (rhs.w*vqt);
  return Quaternion(vtmp, tmp);
}




// ----------------------------------------------------------------------------
// A.WACHS - Mai 2009 - Creation
Quaternion Quaternion::operator,( const Vecteur& rhs ) const
{
  Scalar tmp = - vqt*rhs;
  Vecteur vtmp = (vqt^rhs) + (w*rhs);
  return Quaternion(vtmp, tmp);
}




// ----------------------------------------------------------------------------
// A.WACHS - Mai 2009 - Creation
Quaternion Quaternion::multLeftVec( const Vecteur& lhs ) const
{
  Scalar tmp = -lhs*vqt;
  Vecteur vtmp = (lhs^vqt) + (lhs*w);
  return Quaternion(vtmp, tmp);
}




// ----------------------------------------------------------------------------
// A.WACHS - Mai 2009 - Creation
Quaternion operator,( const Vecteur& lhs, const Quaternion& q )
{
  return q.multLeftVec(lhs);
}




// ----------------------------------------------------------------------------
// A.WACHS - Aout 2009 - Creation
Vecteur Quaternion::multToVecteur( const Quaternion& rhs ) const
{
  Vecteur vtmp = (vqt^rhs.vqt) + (w*rhs.vqt) + (rhs.w*vqt);
  return vtmp;
}




// ----------------------------------------------------------------------------
// A.WACHS - Aout 2011 - Creation
Vecteur Quaternion::multConjugateToVecteur( const Quaternion& rhs ) const
{
  Vecteur vtmp = - (vqt^rhs.vqt) - (w*rhs.vqt) + (rhs.w*vqt);
  return vtmp;
}



/* Build a unit quaternion representing the rotation
 * from u to v. The input vectors need not be normalised. */
void Quaternion::setFromRotTwoVectors( const Vecteur& u, const Vecteur& v)
{
    double norm_u_norm_v = sqrt( (u*u) * (v*v) );
    double real_part = norm_u_norm_v + u*v;
    Vecteur vect;

    if (real_part < 1.e-6 * norm_u_norm_v)
    {
        /* If u and v are exactly opposite, rotate 180 degrees
         * around an arbitrary orthogonal axis. Axis normalisation
         * can happen later, when we normalise the quaternion. */
        real_part = 0.;
        vect = fabs(u[0]) > fabs(u[2]) ? Vecteur(-u[1], u[0], 0.)
                                : Vecteur(0., -u[2], u[1]);
    }
    else
    {
        /* Otherwise, build quaternion the standard way. */
        vect = u^v;
    }

    *this = Quaternion(vect[0],vect[1],vect[2],real_part) * (1/
            Norm( Quaternion(vect[0],vect[1],vect[2],real_part)) );
}


Vecteur Quaternion::rotateVector(const Vecteur v)
{
  Vecteur v_rotated;
  v_rotated = (w*w - Norm(vqt)*Norm(vqt))*v + 2*(v*vqt)*vqt
              + 2*w*(vqt^v);
  return v_rotated ;
}


// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
bool Quaternion::operator== ( const Quaternion& rhs )
{
  return (w==rhs.w &&
	  vqt[0]==rhs.vqt[0] && vqt[1]==rhs.vqt[1] && vqt[2]==rhs.vqt[2]);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
bool Quaternion::operator!= ( const Quaternion& rhs )
{
  return !(*this==rhs);
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// surcharge de l'operateur de sortie d'un quaternion : scalaire
// puis les composantes du vecteur.
ostream &operator << ( ostream &fileOut, const Quaternion &objet )
{
  fileOut << objet.w << '\t' << objet.vqt;
  return (fileOut);
}
istream &operator >> (istream &fileIn, Quaternion &objet)
{
  fileIn >> objet.w >> objet.vqt;
  return (fileIn);
}




// ----------------------------------------------------------------------
// Operateur d'ecriture
// A.WACHS - Janv.2011 - Creation.
void Quaternion::writeQuaternion( ostream &fileOut ) const
{
  fileOut << Grains_Exec::doubleToString(ios::scientific,POSITIONFORMAT,
  	w) << " ";
  vqt.writeGroup3(fileOut);
}




// ----------------------------------------------------------------------
// Operateur d'ecriture en binaire
// A.WACHS - Dec.2014 - Creation.
void Quaternion::writeQuaternion_binary( ostream &fileOut )
{
  fileOut.write( reinterpret_cast<char*>( &w ), sizeof(double) );
  vqt.writeGroup3_binary( fileOut );
}




// ----------------------------------------------------------------------
// Operateur de lecture en binaire
// A.WACHS - Dec.2014 - Creation.
void Quaternion::readQuaternion_binary( istream &StreamIN )
{
  StreamIN.read( reinterpret_cast<char*>( &w ), sizeof(double) );
  vqt.readGroup3_binary( StreamIN );
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// norme d'un quaternion
Scalar ABS( const Quaternion& v )
{
  return (sqrt(v.w*v.w +
	       v.vqt[0]*v.vqt[0] + v.vqt[1]*v.vqt[1] + v.vqt[2]*v.vqt[2]));
}




// ----------------------------------------------------------------------------
// F.PRADEL - Janvier 2000 - Creation
// norme d'un quaternion
Scalar Norm( const Quaternion& q )
{
  return (sqrt(q.vqt[X]*q.vqt[X]+q.vqt[Y]*q.vqt[Y]+q.vqt[Z]*q.vqt[Z]+q.w*q.w));
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin 2000 - Creation
// operateur de norme du quaternion au carre
Scalar Norm2( const Quaternion& qt )
{
  return (qt.vqt[X]*qt.vqt[X]+qt.vqt[Y]*qt.vqt[Y]+qt.vqt[Z]*qt.vqt[Z]
  	+qt.w*qt.w);
}
