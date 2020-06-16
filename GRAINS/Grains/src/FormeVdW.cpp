#include "MPIWrapperGrains.hh"
#include "FormeVdW.H"
#include "Box.H"
#include "Convex_BuilderFactory.H"
#include "PointContact.hh"
#include "PointC.H"
#include "Grains_Exec.hh"
#include "Particule.H"
#include <fstream>
#include <sstream>
#include "CohContact.H"


// ----------------------------------------------------------------------------
// Constructeur par defaut
FormeVdW::FormeVdW() :
  Forme(),
  m_rayonVdW( 0.0 ),
  m_scaling( NULL ),
  m_positionVdW( NULL ),
  m_positionVdW_computed( false )
{}




// ----------------------------------------------------------------------------
// Constructeur par copie
FormeVdW::FormeVdW( const FormeVdW &forme ) :
  Forme( forme ),
  m_rayonVdW( forme.m_rayonVdW ),
  m_scaling( NULL ),
  m_positionVdW( NULL ),
  m_positionVdW_computed( false )
{
  if ( forme.m_scaling ) m_scaling = new Vecteur( *forme.m_scaling ) ;
  if ( forme.m_positionVdW )
    m_positionVdW = new Transform( *forme.m_positionVdW );
}




// ----------------------------------------------------------------------------
// Constructeur avec initialisation
FormeVdW::FormeVdW( Convex *convex_, const Transform &position_ ) :
  Forme( convex_, position_ ),
  m_rayonVdW( 0.0 ),
  m_scaling( NULL ),
  m_positionVdW( NULL ),
  m_positionVdW_computed( false )
{}




// ----------------------------------------------------------------------------
// Contructeur avec lecture du fichier de persistance
// D. RAKOTONIRINA - Oct. 2014 - Modification
// Prise en compte de la cle = *CompParticule et des types de particules
// elemenataires dans le cas d'une particule composite. Par defaut
// type.empty()	=> Simu avec des particules convexes
FormeVdW::FormeVdW( istream &fileIn, string type )
{
  string cle;
  if ( !type.empty() ) cle = type;
  else fileIn >> cle;

  m_convex = Convex_BuilderFactory::create( cle, fileIn );
  while ( cle != "*END" )
  {
    if ( cle == "*RayonInteraction" ) fileIn >> m_rayonVdW;
    fileIn >> cle;
  }

  m_rayon = m_convex->BuildRayon( m_position );
  Shrinking = m_convex->getShrinkingMode();
  m_scaling = new Vecteur;
  BBox box = m_convex->bbox( TransformIdentity );
  const Vecteur& extent = box.getExtent();

  (*m_scaling)[X] = extent[X] < EPSILON ?
  	1. : ( extent[X] - m_rayonVdW ) / extent[X];
  (*m_scaling)[Y] = extent[Y] < EPSILON ?
  	1. : ( extent[Y] - m_rayonVdW ) / extent[Y];
  (*m_scaling)[Z] = extent[Z] < EPSILON ?
  	1. : ( extent[Z] - m_rayonVdW ) / extent[Z];

  m_positionVdW = new Transform();
  m_positionVdW_computed = false ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur avec decodage du noeud XML
FormeVdW::FormeVdW( DOMNode *root )
{
  // Convex
  DOMNode* forme = ReaderXML::getNode( root, "Convex" );
  m_convex = Convex_BuilderFactory::create( forme );

  m_rayonVdW = ReaderXML::getNodeAttr_Double( forme, "RayonInteraction" );

  // Orientation & Position
  m_position.load( root );

  m_rayon = m_convex->BuildRayon( m_position );
  Shrinking = m_convex->getShrinkingMode();
  m_scaling = new Vecteur;
  BBox box = m_convex->bbox( TransformIdentity );
  const Vecteur& extent = box.getExtent();

  (*m_scaling)[X] = extent[X] < EPSILON ?
        1. : ( extent[X] - m_rayonVdW ) / extent[X];
  (*m_scaling)[Y] = extent[Y] < EPSILON ?
        1. : ( extent[Y] - m_rayonVdW ) / extent[Y];
  (*m_scaling)[Z] = extent[Z] < EPSILON ?
        1. : ( extent[Z] - m_rayonVdW ) / extent[Z];

  m_positionVdW = new Transform();
  m_positionVdW_computed = false ;
}




// ----------------------------------------------------------------------------
// Destructeur
FormeVdW::~FormeVdW()
{
  if ( m_scaling ) delete m_scaling;
  if ( m_positionVdW ) delete m_positionVdW;
}




// ----------------------------------------------------------------------------
// BBox de la forme avec une matrice de transformation modifiee
// On cherche a tenir compte de m_rayonVdW
BBox FormeVdW::BoxForme() const
{
  BBox box = m_convex->bbox( m_position );

  Vecteur extent = box.getExtent(), vec_addVdW( m_rayonVdW );
  extent += vec_addVdW;
  box.setExtent( extent );

  return box;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Grandeurs geometriques caracterisant le contact entre les deux
// formes : point, recouvrement & distance
PointContact FormeVdW::ClosestPoint( FormeVdW &voisine )
  throw(ErreurContact)
{
  const Convex* convexA = m_convex;
  const Convex* convexB = voisine.m_convex;

  // Remarque sur la direction du vecteur de recouvrement
  // Soit A et B les centres des 2 convexes
  // recouvrement = overlap * Vecteur(A vers B)
  // Si contact: overlap < 0. et recouvrement est dirige de B vers A
  // Si pas de contact, la direction importe peu

  // Si deux spheres ou 2 disques 2D, on bascule vers une Intersection optimisee
  if ( convexA->getConvexType() == SPHERE
  	&& convexB->getConvexType() == SPHERE )
    return ClosestPointSPHERE( *this, voisine );
  if ( convexA->getConvexType() == DISQUE2D
   		&& convexB->getConvexType() == DISQUE2D )
    return ClosestPointSPHERE( *this, voisine );

  // De m�me pour une Intersection sphere-Box ou disque2D-Box
  if ( convexA->getConvexType() == SPHERE && convexB->getConvexType() == BOX )
    return ClosestPointSPHEREBOX( *this, voisine );
  if ( convexA->getConvexType() == BOX && convexB->getConvexType() == SPHERE )
    return ClosestPointSPHEREBOX( *this, voisine );
  if ( convexA->getConvexType() == DISQUE2D && convexB->getConvexType() == BOX )
    return ClosestPointSPHEREBOX( *this, voisine );
  if ( convexA->getConvexType() == BOX && convexB->getConvexType() == DISQUE2D )
    return ClosestPointSPHEREBOX( *this, voisine );

  // Cas general
  Vecteur gcagcb = *m_position.getOrigin() - *voisine.m_position.getOrigin();
  if ( Norm(gcagcb) < m_rayon + voisine.m_rayon )
  {
    // Distance entre les deux formes
    Point pointA, pointB;
    int nbIterGJK = 0;
    Transform const* a2w = this->getTransformVdW();
    Transform const* b2w = voisine.getTransformVdW();
    // cout << "a2w=" << *a2w << endl ;
    // cout << "b2w=" << *b2w << endl ;
    Scalar distance = closest_points( *m_convex, *(voisine.m_convex), *a2w,
    	*b2w, pointA, pointB, nbIterGJK );
    if ( distance < EPSILON )
    {
      cout << "ERR FormeVdW::ClosestPoint on Processor "
      	<< (Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
	<< endl;
      throw ErreurContact();
    }

    // DGR : les points A et B sont dans leur repere local
    // Donc on les transforme en WC (world coordinates)
    pointA = (*a2w)( pointA );
    pointB = (*b2w)( pointB );

    // Remarque sur le vecteur ba:
    // pointA est le point de contact du convexe A
    // pointB est le point de contact du convexe B
    // donc pointA - pointB = ba est dirige de B vers A
    Vecteur ba = pointA - pointB;

    // Definition du point de contact en WC
    Point contact = pointA / 2.0 + pointB / 2.0;

    // Calcul du recouvrement reel * n(ij)
    // Si contact: rwA + rwB - distance > 0, donc recouvrement est bien dirige
    // de B vers A
    // Si pas de contact, rwA + rwB - distance < 0 <=> distance - rwA - rwB > 0
    // et la direction du recouvrement n'importe pas
    Vecteur recouvrement = ba / distance;
    recouvrement.round();
    recouvrement *= m_rayonVdW + voisine.m_rayonVdW - distance;

    // Distance relle de recouvrement: distance = distance - rwA - rwB
    // si distance < 0. => contact
    // sinon pas de contact
    distance -= m_rayonVdW + voisine.m_rayonVdW;

    return PointContact( contact, recouvrement, distance, nbIterGJK, 0 );
  }
  else return PointNoContact;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Variante de la methode ci-dessus pour gerer les erreurs lies a des
// contacts non detectes par GJK
PointContact FormeVdW::ClosestPoint_ErreurHandling( const FormeVdW &voisine,
  	const Scalar& factor, const int& id, const int& id_voisine )
  throw(ErreurContact)
{
  // Cas general
  Point pointA, pointB;
  int nbIterGJK = 0;
  Transform a2w = this->getTransformVdW( factor, 0.5 );
  Transform b2w = voisine.getTransformVdW( factor, 0.5 );
  Scalar distance = closest_points( *m_convex, *(voisine.m_convex), a2w, b2w,
	pointA, pointB, nbIterGJK );

  if ( distance < EPSILON )
  {
      cout << "ERR FormeVdW::ClosestPoint_ErreurHandling on Processor "
      	<< (Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
	<< " between components " << id << " and " << id_voisine
	<< endl;
      throw ErreurContact();
  }
  else
  {
    cout << "Handling contact error on Processor "
	<< (Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
	<< " between components " << id << " and " << id_voisine
	<< endl;
  }

  // DGR : les points A et B sont dans leur repere local
  // Donc on les transforme en WC (world coordinates)
  pointA = a2w( pointA );
  pointB = b2w( pointB );

  // Remarque sur le vecteur ba:
  // pointA est le point de contact du convexe A
  // pointB est le point de contact du convexe B
  // donc pointA - pointB = ba est dirige de B vers A
  Vecteur ba = pointA - pointB;
  ba.normalize();

  // Definition du point de contact en WC
  Point contact = pointA / 2.0 + pointB / 2.0;

  Scalar recouvrement_impose = 1. * ( m_rayonVdW
  	+ voisine.m_rayonVdW );

  return PointContact( contact, recouvrement_impose * ba, - recouvrement_impose,
  	nbIterGJK, 0 );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Point d'intersection cas specifique SPHERE-SPHERE
PointContact ClosestPointSPHERE( const FormeVdW &formeA,
	const FormeVdW &formeB )
  throw(ErreurContact)
{
  // Remarque sur la direction du vecteur de recouvrement
  // Soit A et B les centres des 2 convexes
  // recouvrement = overlap * AB
  // Si contact: overlap < 0. et recouvrement est dirige de B vers A
  // Si pas de contact, la direction importe peu

  Point const* pointA  = formeA.getTransform()->getOrigin();
  Point const* pointB  = formeB.getTransform()->getOrigin();
  bool with_cohesion = Grains_Exec::m_withCohesion ;

  Vecteur vecteurAB = *pointB - *pointA;
  Scalar  rayonA    = formeA.getRayon();
  Scalar  rayonB    = formeB.getRayon();

  Scalar  distance  = Norm( vecteurAB ) - ( rayonA + rayonB );
  if( distance>0. )
  {
    if( with_cohesion )
    {
      Point contactPoint  = *pointA + ( rayonA + 0.5 * distance ) *
        vecteurAB / Norm( vecteurAB );
      Vecteur recouvrement = distance * ( vecteurAB / Norm( vecteurAB ) );
      return PointContact( contactPoint, recouvrement, distance, 0, 0 );
    }
    else
      return PointNoContact;
  }
  else
  {
    Scalar rdwA = formeA.getRayonInterAction();
    Scalar rdwB = formeB.getRayonInterAction();
    if ( - distance >= rdwA + rdwB )
    {
      cout << "ERR FormeVdW::ClosestPointSPHERE on Processor "
      	<< (Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
	<< ": " << - distance << " & " << rdwA + rdwB << "\n";
      throw ErreurContact();
    }

    Point contact  = *pointA + ( rayonA + 0.5 * distance ) *
    	vecteurAB / Norm( vecteurAB );
    Vecteur recouvrement = distance * ( vecteurAB / Norm( vecteurAB ) );

    return PointContact( contact, recouvrement, distance, 1, 0 );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Point d'intersection cas specifique SPHERE-BOX
PointContact ClosestPointSPHEREBOX( const FormeVdW &formeA,
	const FormeVdW &formeB )
  throw(ErreurContact)
{
  // Remarque sur la direction du vecteur de recouvrement
  // Soit A et B les centres des 2 convexes
  // recouvrement = overlap * AB
  // Si contact: overlap < 0. et recouvrement est dirige de B vers A
  // Si pas de contact, la direction importe peu

  const Convex* convexA = formeA.getConvex();
  const Convex* convexB = formeB.getConvex();
  Scalar rdwA = formeA.getRayonInterAction();
  Scalar rdwB = formeB.getRayonInterAction();
  Scalar overlap=0.;
  Point contactPoint, contact;

  if ( convexA->getConvexType() == SPHERE
  	|| convexA->getConvexType() == DISQUE2D )
  {
    Box const* convexBoxB = (Box const*)(convexB);
    Point const* pointA  = formeA.getTransform()->getOrigin();
    Scalar rayonA = formeA.getRayon();
    const Transform* transfB = formeB.getTransform();
    Transform w2b;
    w2b.invert( *transfB );
    contactPoint = convexBoxB->IntersectionPointSPHERE( w2b(*pointA), rayonA,
    	overlap );
    if ( overlap < 0. )
    {
      if ( - overlap >=  rdwA + rdwB )
      {
        cout << "ERR FormeVdW::ClosestPointSPHEREBOX on Processor "
      	<< (Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
	<< ": " << - overlap << " & " << rdwA + rdwB << endl;
	throw ErreurContact();
      }
      contact = (*transfB)( contactPoint );
      Vecteur AB = contact - *pointA;
      Vecteur recouvrement( ( overlap / Norm(AB) ) * AB );
      return PointContact( contact, recouvrement, overlap, 1, 0 );
    }
    else return PointNoContact;
  }
  else
  {
    Box const* convexBoxA = (Box const*)(convexA);
    Point const* pointB  = formeB.getTransform()->getOrigin();
    Scalar rayonB = formeB.getRayon();
    const Transform* transfA = formeA.getTransform();
    Transform w2a;
    w2a.invert( *transfA );
    contactPoint = convexBoxA->IntersectionPointSPHERE( w2a(*pointB), rayonB,
    	overlap );
    if ( overlap < 0. )
    {
      if ( - overlap >= rdwA + rdwB )
      {
        cout << "ERR FormeVdW::ClosestPointSPHEREBOX on Processor "
      	<< (Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
	<< ": " << - overlap << " & " << rdwA + rdwB << endl;
	throw ErreurContact();
      }
      contact = (*transfA)( contactPoint );
      Vecteur AB = *pointB - contact;
      Vecteur recouvrement( ( overlap / Norm(AB) ) * AB );
      return PointContact( contact, recouvrement, overlap, 1, 0 );
    }
    else return PointNoContact;
  }
}




// ----------------------------------------------------------------------------
// Valeur du rayon d'approche de la forme
Scalar FormeVdW::getRayonInterAction() const
{
  return m_rayonVdW;
}




// ----------------------------------------------------------------------------
// Valeur du rayon d'approche de la forme
void FormeVdW::setRayonInterAction( Scalar &rayonVdW )
{
  m_rayonVdW = rayonVdW;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Matrice de transformation avec le rayon de Van-der-Waals
Transform const* FormeVdW::getTransformVdW()
{
  if ( !m_positionVdW_computed )
  {
    *m_positionVdW = m_position;
    m_positionVdW->composeWithScaling( (*m_scaling)[X], (*m_scaling)[Y],
  	(*m_scaling)[Z] );
    m_positionVdW_computed = true;
  }

  return ( m_positionVdW );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Matrice de transformation avec le rayon de Van-der-Waals
Transform FormeVdW::getTransformVdW( const Scalar& factor,
	const Scalar& min_scaling ) const
{
  BBox box = m_convex->bbox( TransformIdentity );
  const Vecteur& extent = box.getExtent();

  Scalar scaleX, scaleY, scaleZ;
  scaleX = extent[X] < EPSILON ? 1. :
  	( extent[X] - factor * m_rayonVdW ) / extent[X];
  scaleY = ( extent[Y] < EPSILON ) ? 1. :
  	( extent[Y] - factor * m_rayonVdW ) / extent[Y];
  scaleZ = ( extent[Z]<EPSILON ) ? 1. :
  	( extent[Z] - factor * m_rayonVdW ) / extent[Z];

  if ( min_scaling != 1. )
  {
    scaleX = max( scaleX, min_scaling );
    scaleY = max( scaleY, min_scaling );
    scaleZ = max( scaleZ, min_scaling );
  }

  Transform vdw = m_position;
  vdw.composeWithScaling( scaleX, scaleY, scaleZ );

  return vdw;
}




// ----------------------------------------------------------------------------
// Evaluation d'un contact entre deux formes. ATTENTION: test des BBox
// et non veritable contact geometrique
bool FormeVdW::isProche( const FormeVdW & voisine ) const
{
  bool contact;

  BBox boxA = (*this).BoxForme();
  BBox boxB = voisine.BoxForme();
  const Vecteur& extentA = boxA.getExtent();
  const Vecteur& extentB = boxB.getExtent();
  const Point&   pointA  = boxA.getCenter();
  const Point&   pointB  = boxB.getCenter();
  Vecteur ab = pointA - pointB;

  double x, y, z;
  double vdwA = (*this).m_rayonVdW;
  double vdwB = voisine.m_rayonVdW;
  x = fabs( ab[X] ) - ( extentA[X]
  	+ extentB[X] * ( extentA[X] + vdwA ) / extentA[X]
  	* ( extentB[X] + vdwB ) / extentB[X] );
  y = fabs( ab[Y] ) - ( extentA[Y]
  	+ extentB[Y] * ( extentA[Y] + vdwA ) / extentA[Y]
  	* ( extentB[Y] + vdwB ) / extentB[Y] );
  z = fabs( ab[Z] ) - ( extentA[Z]
  	+ extentB[Z] * ( extentA[Z] + vdwA ) / extentA[Z]
  	* ( extentB[Z] + vdwB ) / extentB[Z] );

  contact = true;
  if ( x>0. || y>0. || z>0. ) contact = false;

  //  contact = intersect((*this).BoxForme(), voisine.BoxForme());
  return contact;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Les deux formes sont-elles en contact ? veritable contact geometrique.
bool FormeVdW::isContact( FormeVdW &voisine )
{
  bool contact = false;

  const Convex* convexA = m_convex;
  const Convex* convexB = voisine.m_convex;

  // Si deux spheres ou 2 disques 2D, on bascule vers une Intersection optimisee
  if ( convexA->getConvexType() == SPHERE
  	&& convexB->getConvexType() == SPHERE )
    return isContactSPHERE( *this,  voisine );
  if ( convexA->getConvexType() == DISQUE2D
  	&& convexB->getConvexType() == DISQUE2D )
    return isContactSPHERE( *this, voisine );

  // De m�me pour une Intersection sphere-Box ou disque2D-Box
  if ( convexA->getConvexType() == SPHERE && convexB->getConvexType() == BOX )
    return isContactSPHEREBOX( *this, voisine );
  if ( convexA->getConvexType() == BOX && convexB->getConvexType() == SPHERE )
    return isContactSPHEREBOX( *this, voisine );
  if ( convexA->getConvexType() == DISQUE2D && convexB->getConvexType() == BOX )
    return isContactSPHEREBOX( *this, voisine );
  if ( convexA->getConvexType() == BOX && convexB->getConvexType() == DISQUE2D )
    return isContactSPHEREBOX( *this, voisine );

  // Cas general
  Point pointA, pointB;
  int nbIterGJK = 0;
  Transform const* a2w = this->getTransformVdW();
  Transform const* b2w = voisine.getTransformVdW();
  Scalar distanceMin = (*this).m_rayonVdW + voisine.m_rayonVdW - EPSILON;
  Scalar distance = closest_points( *m_convex, *(voisine.m_convex), *a2w, *b2w,
	pointA, pointB, nbIterGJK );

  if ( distance < distanceMin ) contact = true;

  return contact;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Y a t il contact, cas specifique SPHERE-SPHERE
bool isContactSPHERE( const FormeVdW &formeA,
		const FormeVdW &formeB )
{
  bool contact = false;

  Point const* pointA  = formeA.getTransform()->getOrigin();
  Point const* pointB  = formeB.getTransform()->getOrigin();

  Vecteur vecteurAB = *pointB - *pointA;
  Scalar  rayonA    = formeA.getRayon();
  Scalar  rayonB    = formeB.getRayon();

  Scalar  distance  = Norm( vecteurAB ) - ( rayonA + rayonB );
  if ( distance < 0. ) contact = true;

  return contact;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Point d'intersection cas specifique SPHERE-BOX
bool isContactSPHEREBOX( const FormeVdW &formeA,
		const FormeVdW &formeB )
{
  bool contact = false;

  const Convex* convexA = formeA.getConvex();
  const Convex* convexB = formeB.getConvex();
  Scalar overlap=0.;
  Point contactPoint;

  if ( convexA->getConvexType() == SPHERE
  	|| convexA->getConvexType() == DISQUE2D )
  {
    Box const* convexBoxB = (Box const*)(convexB);
    Point const* pointA  = formeA.getTransform()->getOrigin();
    Scalar rayonA = formeA.getRayon();
    const Transform* transfB = formeB.getTransform();
    Transform w2b;
    w2b.invert( *transfB );
    contactPoint = convexBoxB->IntersectionPointSPHERE( w2b(*pointA), rayonA,
    	overlap );
    if ( overlap < 0. ) contact = true;
  }
  else
  {
    Box const* convexBoxA = (Box const*)(convexA);
    Point const* pointB  = formeB.getTransform()->getOrigin();
    Scalar rayonB = formeB.getRayon();
    const Transform* transfA = formeA.getTransform();
    Transform w2a;
    w2a.invert( *transfA );
    contactPoint = convexBoxA->IntersectionPointSPHERE( w2a(*pointB), rayonB,
    	overlap );
    if ( overlap < 0. ) contact = true;
  }

  return contact;
}




// ----------------------------------------------------------------------------
// Ecriture de l'information statique
// G.FERRER - Janv.2002 - Creation
// D. RAKOTONIRINA - Sept. 2014 - Modification
// Prise en compte d'une particule composite. Si PARTICULE COMPOSITE, le
// pointeur sur composant N'EST PAS NUL. Il est NUL pour les autres cas
void FormeVdW::writeStatique( ostream &statique , Composant const* composant )
{
  if ( composant )
  {
    string name = composant->getPartName() ;
    statique << "*CompParticule\n" << composant->getVolume()
       <<"\t" << name << "\n";
    statique << "*RayonInteraction\t" << m_rayonVdW << '\n';
    //statique << "*RayonInteraction\t" << "0.\n";
  }
  else
  {
    statique << *m_convex;
    statique << "*RayonInteraction\t" << m_rayonVdW << '\n';
  }
  statique << "*END\n";
}




// ----------------------------------------------------------------------------
// Possibilite d'intersection entre deux formes.
// G.FERRER - Avri.2003 - Creation
bool intersect( const FormeVdW &a, const FormeVdW &b )
{
  return ( intersect( a.BoxForme(), b.BoxForme() ) );
}




// ----------------------------------------------------------------------------
// Initialize a faux le boolean correspondant au calcul de la
// transformation avec scaling par l'epaisseur de croute
void FormeVdW::initializeVdWtransform_to_notComputed()
{
  m_positionVdW_computed = false;
}
