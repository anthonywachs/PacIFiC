#include "Grains_BuilderFactory.H"
#include "Forme.H"
#include "Convex_BuilderFactory.H"
#include "Torseur.H"
#include "Vecteur.H"
#include <math.h>
#include <stdlib.h>


// ----------------------------------------------------------------------------
// Constructeur par defaut
// F.PRADEL - Janv.2000 - Creation
// D.PETIT  - Juin.2000 - Modif
// M.SULAIMAN  - Nov.2015 - Modif
Forme::Forme() :
  m_convex( NULL ),
  m_rayon( 0.0 ),
  Shrinking( 0 )
{
}




// ----------------------------------------------------------------------------
// Constructeur avec initialisation
// G.FERRER - Fevr.2002 - Creation
// M. SULAIMAN - Nov. 2015 - Modification  - Shrinking choice called here
Forme::Forme( Convex *convex_, const Transform &position_ ) :
  m_convex( convex_ )
{
  m_position.composeTransformLeft(position_);
  m_rayon = m_convex->BuildRayon( m_position );
  Shrinking = m_convex->getShrinkingMode();
}




// ----------------------------------------------------------------------------
// Constructeur par copie
// D.PETIT - Juin. 2000 - Creation
Forme::Forme( const Forme& forme ) :
  m_convex( NULL ),
  m_rayon( forme.m_rayon )
{
  m_position.composeTransformLeft(forme.m_position);
  if ( forme.m_convex ) m_convex = forme.m_convex->clone();
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Modif
// destructeur
Forme::~Forme()
{
  delete m_convex;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine la BBox d'une forme
BBox Forme::BoxForme() const
{
  return m_convex->bbox( m_position );
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine l'inertie et l'inertie inverse de la forme.
bool Forme::BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const
{
  bool status = m_convex->BuildInertie( inertie, inertie_1 );

  Matrix m( inertie[0], inertie[1], inertie[2],
	inertie[1], inertie[3], inertie[4],
	inertie[2], inertie[4], inertie[5] );
  Matrix m_1( inertie_1[0], inertie_1[1], inertie_1[2],
	inertie_1[1], inertie_1[3], inertie_1[4],
	inertie_1[2], inertie_1[4], inertie_1[5] );

  Matrix base = m_position.getBasis();
  m   = base * m   * base.transpose();
  m_1 = base * m_1 * base.transpose();

  inertie[0] = m[0][0];
  inertie[1] = m[1][0];
  inertie[2] = m[2][0];
  inertie[3] = m[1][1];
  inertie[4] = m[2][1];
  inertie[5] = m[2][2];

  inertie_1[0] = m_1[0][0];
  inertie_1[1] = m_1[1][0];
  inertie_1[2] = m_1[2][0];
  inertie_1[3] = m_1[1][1];
  inertie_1[4] = m_1[2][1];
  inertie_1[5] = m_1[2][2];

  return status;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
Scalar Forme::DistanteDe( const Forme &voisine ) const
{
  Point A, B;
  int nbIter;
  return closest_points( *m_convex, *voisine.m_convex, m_position,
  	voisine.m_position, A, B, nbIter );
}




// ----------------------------------------------------------------------------
// G.FERRER - Fevr.2000 - Creation
// D.PETIT - Juin. 2000 - Modif
// Valeurs des coordonnees du centre
Point const* Forme::getCentre() const
{
  return m_position.getOrigin();
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin.2000 - Modif
void Forme::getCentre( Scalar *pos ) const
{
  Point const* ori = m_position.getOrigin() ;
  pos[X] =  (*ori)[X];
  pos[Y] =  (*ori)[Y];
  pos[Z] =  (*ori)[Z];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Acces au Convex decrivant la Forme
// G.FERRER - Aout.2004 - Creation
const Convex* Forme::getConvex() const
{
  return m_convex;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Acces a la Transformation du Convex decrivant la Forme
const Transform* Forme::getTransform() const
{
  return &m_position;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Acces a la Transformation du Convex decrivant la Forme
Transform* Forme::getTransform()
{
  return &m_position;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Valeurs de la matrice donnant l'orientation de la forme
Matrix Forme::getOrientation() const
{
  return m_position.getBasis();
}




// ----------------------------------------------------------------------------
// Valeur du rayon externe de la particule
// G.FERRER - Juin.2000 - Creation
Scalar Forme::getRayon() const
{
  return m_rayon;
}




// ----------------------------------------------------------------------------
// Valeur du rayon externe de la particle retrecissante
// M.SULAIMAN - Nov.2015 - Creation
int Forme::getShrinkingMode()
{
  return Shrinking;
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Creation
// Determine le Volume d'une forme
Scalar Forme::getVolume() const
{
  return m_convex->getVolume();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Intersection entre les deux formes ?
// G.FERRER - Aout.2004 - Forme::isContact -> intersect
bool intersect( const Forme &a, const Forme &b )
{
  Vecteur v = *(b.m_position.getOrigin()) - *(a.m_position.getOrigin());
  if ( Norm(v) < EPSILON ) return true;
  return intersect( *a.m_convex, *b.m_convex, a.m_position, b.m_position, v );
}




// ----------------------------------------------------------------------------
// indicatrice du contact entre la forme courante et la forme en argument
// D.PETIT - Juin. 2000 - Creation
bool Forme::isContact( const Forme &voisine ) const
{
  bool contact;

  // Vecteur arbitraire
  Vecteur v = *(voisine.m_position.getOrigin()) - *(m_position.getOrigin());
  if ( Norm(v) < EPSILON ) return true;
  contact = intersect( *m_convex, *(voisine.m_convex),
	m_position, voisine.m_position, v );

  return ( contact );
}




// ----------------------------------------------------------------------------
// Verification de proximite entre deux formes
// G.FERRER - Dec.2000 - Creation
bool Forme::isProche( const Forme &voisine ) const
{
  bool contact;
  BBox box0 = (*this).BoxForme();
  BBox box1 = voisine.BoxForme();
  contact = intersect( box0, box1 );

  return ( contact );
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil 2000 - Creation
void Forme::Rotate( const Quaternion& q )
{
  m_position.rotateOrientation( q );
}




// ----------------------------------------------------------------------------
// D.PETIT - Juil 2000 - Creation
// positionne l'origine de la forme
void Forme::setOrigin( const Scalar *pos )
{
  m_position.setOrigin( pos );
}




// ----------------------------------------------------------------------------
// A.WACHS - Mai 2009 - Creation
// positionne l'origine de la forme
void Forme::setOrigin( const double &gx, const double &gy, const double &gz )
{
  m_position.setOrigin( gx, gy, gz );
}




// ----------------------------------------------------------------------------
// D. RAKOTNIRINA - Mars 2014 - Modification
// positionne l'origine de la forme
void Forme::setOrigin( const Point &pos )
{
  m_position.setOrigin( pos[X], pos[Y], pos[Z] );
}




// ----------------------------------------------------------------------------
// Positionne une forme dans l'espace.
void Forme::setPosition( const Scalar *pos )
{
  m_position.setValue( pos );
}




// ----------------------------------------------------------------------------
// Positionne une forme dans l'espace.
void Forme::setRayon( Scalar r )
{
  m_rayon = r;
}




// ----------------------------------------------------------------------------
// Valeur du rayon externe de la particle retrecissante
// M.SULAIMAN - Nov.2015 - Creation
void Forme::set_shrinking_radius( Scalar R )
{
  m_rayon = R;
  m_convex->set_shrinking_radius( m_rayon );
}




// ----------------------------------------------------------------------------
// Affectation d'une Transformation imposee
void Forme::setTransform( const Transform &transform_ )
{
  m_position = transform_;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Applique la transformation trot � la forme
// D. RAKOTONIRINA - Fev 2014 - Modification
void Forme::composeTransform( const Transform &trot )
{
  m_position.composeTransformLeft( trot );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Compose la transformation de la forme par trot:
// !!! Resultat : d'abord la transformation de la forme puis trot !!!
// D. RAKOTONIRINA - Fev 2014 - Creation
void Forme::composeTransformRight( const Transform &trot )
{
  m_position.composeTransformRight( trot );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Compose la transformation de la forme par une rotation trot par
// rapport au centre de gravite de la forme
// !!! Resultat : d'abord la transformation de la forme puis trot !!!
// !!! Cette composition ne change pas l'origine de la forme !!!
// D. RAKOTONIRINA - Fev 2014 - Creation
void Forme::composeRotationRight( const Transform &trot )
{
  m_position.composeRotationRight( trot );
}




// ----------------------------------------------------------------------------
// D.PETIT - Juin. 2000 - Modif
// operateur creant le vecteur F-A reliant le point A au centre de la forme F.
Vecteur Forme::operator-( const Point& rhs ) const
{
  return Vecteur( *m_position.getOrigin() - rhs );
}




// ----------------------------------------------------------------------------
// Egalite sur les adresses des deux formes
bool Forme::operator==( const Forme &particule ) const
{
  return ( this == &particule );
}



// ----------------------------------------------------------------------------
// Deplacement de la forme d'un vecteur rhs
Forme& Forme::operator+=( const Vecteur& rhs )
{
  m_position.translate( rhs );
  return *this;
}




// ----------------------------------------------------------------------------
// Operateur d'ecriture
// G.FERRER - Juin.2000 - Creation
ostream &operator << ( ostream &fileOut, const Forme &forme )
{
  fileOut << "*Forme\n";
  fileOut << *forme.m_convex;
  fileOut << forme.m_position;
  fileOut << "*END\n";

  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Operateur de lecture
// G.FERRER - Juin.2000 - Creation
istream &operator >> ( istream &fileIn, Forme &forme )
{
  string cle;
  fileIn >> cle;
  if ( forme.m_convex ) fileIn >> *forme.m_convex;
  else forme.m_convex = Convex_BuilderFactory::create( cle, fileIn );

  fileIn >> cle;
  while ( cle != "*END" )
  {
    if ( cle == "%Position&Orientation" ) fileIn >> forme.m_position;
    fileIn >> cle;
  }
  forme.m_rayon = forme.m_convex->BuildRayon( forme.m_position );

  return (fileIn);
}




// ----------------------------------------------------------------------------
// Lecture de l'information de position
// G.FERRER - Janv.2001 - Creation
void Forme::readPosition( istream &fileIn )
{
  string buffer;
  fileIn >> buffer >> m_position;
  m_rayon = m_convex->BuildRayon( m_position );
}




// ----------------------------------------------------------------------------
// Lecture de l'information de position au format de reload 2014
// A.WACHS - Aout 2014 - Creation
void Forme::readPosition2014( istream &fileIn )
{
  m_position.readTransform2014( fileIn );
  m_rayon = m_convex->BuildRayon( m_position );
}




// ----------------------------------------------------------------------------
// Lecture de l'information de position au format de reload 2014 en binaire
// A.WACHS - Aout 2014 - Creation
void Forme::readPosition2014_binary( istream &fileIn )
{
  m_position.readTransform2014_binary( fileIn );
  m_rayon = m_convex->BuildRayon( m_position );
}




// ----------------------------------------------------------------------------
// Lecture de l'information statique
// G.FERRER - Janv.2001 - Creation
void Forme::readStatique( istream &fileIn )
{
  string cle;
  fileIn >> cle;
  if ( m_convex ) fileIn >> *m_convex;
  else m_convex = Convex_BuilderFactory::create( cle, fileIn );
}




// ----------------------------------------------------------------------------
// Ecriture de l'information de position
// G.FERRER - Janv.2001 - Creation
void Forme::writePosition( ostream &fileOut ) const
{
  m_position.writeTransform( fileOut );
}




// ----------------------------------------------------------------------------
// Ecriture de l'information statique
// G.FERRER - Janv.2001 - Creation
void Forme::writeStatique( ostream &statique, Composant const * composant )
{
  statique << *m_convex << "*END\n";
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de l'information de position pour le Fluide
// ATTENTION: valide pour des polygones, polyedres, spheres et cylindres
// circulaires seulement !!
void Forme::writePositionInFluid( ostream &fluid )
{
  Point     pointEnvelop;
  vector<Point> allPoints = m_convex->getEnveloppe();
  vector<Point>::iterator point;

  // m_rayon is not given to stream from here any more
  // it's given from InterfaceFluide2D.cpp//InterfaceFluide3D.cpp line 287
  // fluid << m_rayon << " " << allPoints.size() << "\n";
   fluid << " " << allPoints.size() << "\n";

   // Cas 2D
  if ( Grains_BuilderFactory::getContext() == DIM_2 )
  {
    // Points du polygone
    for (point=allPoints.begin(); point!=allPoints.end(); point++)
    {
      pointEnvelop = m_position(*point);
      fluid << pointEnvelop[X] << " " << pointEnvelop[Y] << "\n";
    }
  }
  // Cas 3D
  else if ( Grains_BuilderFactory::getContext() == DIM_3 )
  {
    // Points du polyedre
    for (point=allPoints.begin(); point!=allPoints.end(); point++)
    {
      pointEnvelop = m_position(*point);
      fluid << pointEnvelop[X] << " " << pointEnvelop[Y] << " "
	<< pointEnvelop[Z] << "\n";
//       fluid << Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
//       		pointEnvelop[X] ) << " "
//       	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
// 		pointEnvelop[Y] ) << " "
// 	<< Grains_Exec::doubleToString( ios::scientific, POSITIONFORMAT,
// 		pointEnvelop[Z] ) << "\n";
    }

    // Faces du polyedre
    vector< vector<int> > const* allFaces  = m_convex->getFaces();
    vector< vector<int> >::const_iterator face;
    if ( allFaces )
    {
      fluid << allFaces->size() << '\n';
      for (face=allFaces->begin(); face!=allFaces->end(); face++)
      {
        vector<int>::const_iterator index;
        fluid << (*face).size() << " ";
        for (index=(*face).begin(); index!=(*face).end(); index++)
          fluid << (*index) << " ";
        fluid << '\n';
      }
    }
    else fluid << "0" << '\n';
  }
  else
  {
    cout << "!!! Warning: Physical dimension undefined (DIM_2 or DIM_3)"
    	<< endl;
    exit(0);
  }
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de l'information de position pour Basilisk (Fluide)
// ATTENTION: valide pour des polygones, polyedres, spheres et cylindres
// circulaires seulement !!
void Forme::writePositionInFluid( BasiliskDataStructure * b )
{
  Point     pointEnvelop;
  vector<Point> allPoints = m_convex->getEnveloppe();
  vector<Point>::iterator point;

  long int r;
  r = allPoints.size();
  b->allPoints  = r;


  // Cas 2D
  if ( Grains_BuilderFactory::getContext() == DIM_2 )
  {
    // alloue un tableau de allPoints.size()x2 dans la structure Basilisk
    b->PointsCoord = (double **)malloc(r*sizeof(double *));
    for (long int i = 0; i < r; i++) {
      b->PointsCoord[i] = (double *)malloc(2 * sizeof(double));
      // printf("Form.cpp 2D case\n");
    }
    long int i = 0;
    // Points du polygone
    for (point=allPoints.begin(); point!=allPoints.end(); point++)
    {
      pointEnvelop = m_position(*point);
      b->PointsCoord[i][0] = pointEnvelop[X];
      b->PointsCoord[i][1] = pointEnvelop[Y];
      i++;
    }
  }
  // Cas 3D
  else if ( Grains_BuilderFactory::getContext() == DIM_3 )
  {
    // alloue un tableau de allPoints.size()x3 dans la structure Basilisk
    b->PointsCoord = (double **)malloc(r*sizeof(double *));
    // printf ("Can debugging pointer  b->PointsCoord = %p\n", (void *) b->PointsCoord);
    for (long int i = 0; i < r; i++) {
      b->PointsCoord[i] = (double *)malloc(3 * sizeof(double));
      // printf ("Can debugging pointer c++ b->PointsCoord[i] = %p\n", (void *) b->PointsCoord[i]);
    }
    long int i = 0;
    // Points du polyedre
    for (point=allPoints.begin(); point!=allPoints.end(); point++)
    {
      pointEnvelop = m_position(*point);
      b->PointsCoord[i][0] = pointEnvelop[X];
      b->PointsCoord[i][1] = pointEnvelop[Y];
      b->PointsCoord[i][2] = pointEnvelop[Z];
      i++;
    }

    // Faces du polyedre
    vector< vector<int> > const* allFaces  = m_convex->getFaces();
    vector< vector<int> >::const_iterator face;
    if ( allFaces )
    {
      // Nombre total de faces
      r = allFaces->size();
      b->allFaces  = r;

      // Alloue un tableau dynamique de r "lignes"
      b->FacesIndex = (long int **)malloc(r*sizeof(long int *));

      // Alloue un vecteur de 1 ligne r collones
      b->numPointsOnFaces = (long int *)malloc(r * sizeof(long int));

      i = 0;
      for (face=allFaces->begin(); face!=allFaces->end(); face++)
      {
        vector<int>::const_iterator index;

	// nombre de points sur la face
	long int rr = (*face).size();

	// qu'on stock dans la Structure Basilisk
	b->numPointsOnFaces[i] = rr;

	// Alloue le tableau dynamique
	b->FacesIndex[i] = (long int *)malloc(rr * sizeof(long int));

	long int j = 0;
	for (index=(*face).begin(); index!=(*face).end(); index++)
	  {
	    b->FacesIndex[i][j] = (*index);
	    j++;
	  }
	i++;
      }
    }
  }
  else
  {
    cout << "!!! Warning: Physical dimension undefined (DIM_2 or DIM_3)"
    	<< endl;
    exit(0);
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de l'objet dans sa configuration courante pour
// post-processing avec GMV
void Forme::GMVoutput( ostream &fileOut ) const
{
  m_convex->GMVoutput( fileOut, m_position );
}




// ----------------------------------------------------------------------------
// CopyTransform (poly�dres)
void Forme::copyTransform( double *vit, int i ) const
{
  m_position.copyTransform( vit, i );
}




// ----------------------------------------------------------------------------
// CopyTransform (poly�dres)
void Forme::copyTransform( double *vit, int i, Vecteur const& vec ) const
{
  m_position.copyTransform( vit, i, vec );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit les points du convexe pour post-processing avec Paraview
void Forme::write_polygonsPts_PARAVIEW( ostream &f, Vecteur const* translation )
	const
{
  m_convex->write_polygonsPts_PARAVIEW( f, m_position, translation );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit le convexe au format STL pour lien avec openFoam
void Forme::write_convex_STL( ostream &f ) const
{
  m_convex->write_convex_STL( f, m_position );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les points du convexe pour post-processing avec Paraview
list<Point> Forme::get_polygonsPts_PARAVIEW( Vecteur const* translation ) const
{
  return ( m_convex->get_polygonsPts_PARAVIEW( m_position, translation ) );
}
