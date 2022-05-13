#include "MPIWrapperGrains.hh"
#include "CompParticule.hh"
#include "Contact_BuilderFactory.hh"
#include "Memento.hh"
#include "Cinematique_BuilderFactory.H"
#include "Grains_BuilderFactory.H"
#include "PointC.H"
#include "ElementParticule.H"
#include "SaveTable.H"
#include "Forme.H"
#include "Grains_Exec.hh"
#include "ContactLaw.hh"
#include <iterator>
#include <algorithm>


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
CompParticule::CompParticule( DOMNode* root,
	const bool &autonumbering, const int &pc, const string &name ):
  Particule( autonumbering ),
  m_name( name ),
  m_compVolume( 0.),
  m_CG( Point(0.,0.,0.) )
{
  // Classe de la particule
  m_ParticuleClasse = pc;

  // Position relative de la particule elemenataire
  Point RelativePosition(0.,0.,0.);

  /* -------------------------------------------------------------- */
  /* Construction de la forme                                       */
  /* Le composite n'a pas de forme au sens strict mais possede une  */
  /* transformation geometrique et un rayon circonscrit             */
  /* -------------------------------------------------------------- */

  m_geoFormeVdw = new FormeVdW( new PointC(), Transform() );

  m_cinematique = Cinematique_BuilderFactory::create(
  	m_geoFormeVdw->getConvex() );

  DOMNode* CompPartMateriau = ReaderXML::getNode( root, "Materiau" );

  // Orientation de la particule composite
  m_geoFormeVdw->getTransform()->load( root );

  m_nomMateriau = ReaderXML::getNodeValue_String( CompPartMateriau );

  Contact_BuilderFactory::defineMaterial( m_nomMateriau, false );

  m_masseVolumique = ReaderXML::getNodeAttr_Double( root, "MasseVolumique");

  //Does it move ?
  if ( ReaderXML::hasNodeAttr_String( root, "Mobilite" ) )
  {
    if (ReaderXML::getNodeAttr_String(root, "Mobilite" ) == "fix" )
      is_mobile = false;
    else if (ReaderXML::getNodeAttr_String(root, "Mobilite" ) == "mobile" )
      is_mobile = true;
    else
    {
      cout << "Mobility not recognized, particle is considered as mobile"<<endl;
      is_mobile = true;
    }
  }
  else is_mobile = true;

  DOMNode* allElPart = ReaderXML::getNode( root, "ElementParticules" );

  DOMNodeList* allElemParticules = ReaderXML::getNodes( allElPart );

  // Nombre de particules elemtaires
  m_nbreElemPart = allElemParticules->getLength();

  // Contient les rayons d'interaction de chaque particule elementaire
  vector<Scalar> EltRayonVdW;
  EltRayonVdW.reserve( m_nbreElemPart );
  for ( size_t j=0;j<m_nbreElemPart; ++j ) EltRayonVdW.push_back( 0. );

  m_elementaryParticules.reserve( m_nbreElemPart );
  ElementParticule* ppp = NULL;
  for ( size_t j=0; j<m_nbreElemPart; ++j )
    m_elementaryParticules.push_back( ppp );

  m_InitialRelativePositions.reserve( m_nbreElemPart );
  for ( size_t j=0; j<m_nbreElemPart; ++j )
    m_InitialRelativePositions.push_back( Vecteur(0.,0.,0.) );
  for ( XMLSize_t i=0; i<m_nbreElemPart; ++i )
  {
    DOMNode* nElemParticule = allElemParticules->item( i );
    // Convention : les particules elementaires sont de la meme classe
    // que le composite
    m_elementaryParticules[i] = new ElementParticule( nElemParticule,
    	pc, this, (int) i );

    // Lecture de la position relative par rapport � 0
    DOMNode* nRelPos  = ReaderXML::getNode( nElemParticule,
	"RelativePosition" );
    RelativePosition[X] = ReaderXML::getNodeAttr_Double( nRelPos, "X" );
    RelativePosition[Y] = ReaderXML::getNodeAttr_Double( nRelPos, "Y" );
    RelativePosition[Z] = ReaderXML::getNodeAttr_Double( nRelPos, "Z" );
    m_elementaryParticules[i]->setPosition( RelativePosition );

    // Affecte le materiau du composite a chaque particule elementaire
    m_elementaryParticules[i]->setMateriau( m_nomMateriau ) ;

    // Affecte la masse volumique du composite a chaque particule elementaire
    m_elementaryParticules[i]->setMasseVolumique( m_masseVolumique ) ;

    m_InitialRelativePositions[i] = RelativePosition ;

    // Lecture des rayons d'interaction de particule elementaire
    DOMNode* forme = ReaderXML::getNode( nElemParticule, "Convex" );

    EltRayonVdW[i] = ReaderXML::getNodeAttr_Double( forme,
	"RayonInteraction" );

    m_geoFormeVdw->setRayonInterAction( *min_element(EltRayonVdW.begin(),
	EltRayonVdW.end()) );
  }


  /* -------------------------------------------------------- */
  /* Sauvegarde des transformations initiales de la particule */
  /* composite et des particules elementaires (insert.xml)    */
  /* -------------------------------------------------------- */
  Matrix ttt;
  ttt.setIdentity();
  m_InitialMatrix.reserve( m_elementaryParticules.size() );
  for ( size_t i=0; i<m_InitialMatrix.size(); ++i )
    m_InitialMatrix.push_back( ttt );
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    m_InitialMatrix.push_back( m_elementaryParticules[i]->getForme()
    ->getTransform()->getBasis() );

  if ( !Grains_Exec::m_isGrainsCompFeatures )
  {
    /* ------------------------------ */
    /* Loading des elements d'inertie */
    /* ------------------------------ */
    DOMNode* nInertie = ReaderXML::getNode( root, "Inertie" );
    DOMNode* nVol     = ReaderXML::getNode( nInertie, "Volume" );
    DOMNode* nCG      = ReaderXML::getNode( nInertie, "GravityCenter" );
    DOMNode* nTen     = ReaderXML::getNode( nInertie, "Tenseur" );

    m_compVolume = ReaderXML::getNodeAttr_Double( nVol, "VOLUME" );

    m_CG[X]  = ReaderXML::getNodeAttr_Double( nCG, "xCG");
    m_CG[Y]  = ReaderXML::getNodeAttr_Double( nCG, "yCG");
    m_CG[Z]  = ReaderXML::getNodeAttr_Double( nCG, "zCG");

    m_masse = m_masseVolumique * m_compVolume;

    for ( size_t i=0; i<6; ++i )
      m_compInertie.push_back( 0. );

    m_compInertie[0] =  ReaderXML::getNodeAttr_Double( nTen, "Ixx" );
    m_compInertie[1] = -ReaderXML::getNodeAttr_Double( nTen, "Ixy" );
    m_compInertie[2] = -ReaderXML::getNodeAttr_Double( nTen, "Ixz" );
    m_compInertie[3] =  ReaderXML::getNodeAttr_Double( nTen, "Iyy" );
    m_compInertie[4] = -ReaderXML::getNodeAttr_Double( nTen, "Iyz" );
    m_compInertie[5] =  ReaderXML::getNodeAttr_Double( nTen, "Izz" );

    /* -------------------------------------------------------- */
    /* Mise a jour des positions relatives apres initialisation */
    /* -------------------------------------------------------- */
    for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
      m_RelativePositions.push_back( m_InitialRelativePositions[i] - m_CG );

    /* ------------------------------------------------------ */
    /* Initialisation de la position de la paticule composite */
    /* et des positions des particules elmentaires            */
    /* ------------------------------------------------------ */
    m_geoFormeVdw->setOrigin(0.,0.,0.);
    for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
      m_elementaryParticules[i]->setPosition( m_RelativePositions[i] );

    // Affecte la rotation du composite a toutes les particules elementaires
    Matrix rota = m_geoFormeVdw->getTransform()->getBasis() ;
    Vecteur trans;

    for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    {
      m_elementaryParticules[i]->getForme()
      ->composeTransform( *(m_geoFormeVdw->getTransform()) );

      trans = ( (rota)*m_RelativePositions[i] ) - m_RelativePositions[i];
      m_elementaryParticules[i]->Translate( trans );
    }

    BuildInertie( m_inertie, m_inertie_1 );
    computeWeight();
    setRayon();
  }
}




// ----------------------------------------------------------------------------
// Constructeur par defaut
CompParticule::CompParticule( const bool &autonumbering ):
  Particule( autonumbering )
{
}




// ----------------------------------------------------------------------------
// Constructeur (MPI)
CompParticule::CompParticule( const int &id_, Particule const* ParticuleRef,
	const double &vx, const double &vy, const double &vz,
	const double &qrotationx, const double &qrotationy,
	const double &qrotationz, const double &qrotations,
	const double &rx, const double &ry, const double &rz,
	const Scalar m[16],
	const ParticuleActivity &activ,
	const int &tag_,
	const int &coordination_number_,
	const bool &updatePosition ) :
  Particule( id_, ParticuleRef,
            vx, vy, vz,
            qrotationx, qrotationy,
	    qrotationz, qrotations,
	    rx, ry, rz, m, activ,
	    tag_, coordination_number_ )
{
  /* Particules elementaires */
  m_elementaryParticules.reserve( ParticuleRef->getNbreElemPart() );
  ElementParticule* ppp = NULL;
  for ( size_t i=0; i<ParticuleRef->getNbreElemPart(); ++i )
    m_elementaryParticules.push_back( ppp );

  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    m_elementaryParticules[i] = new ElementParticule(
      *(ParticuleRef->getElementParticules()[i]), this );

  /* Positions initiales des particules elementaires dans insert.xml */
  m_InitialRelativePositions = ParticuleRef->getInitialRelativePositions();

  /* Positions des particules elementaires par rapport a (0.,0.,0.) */
  /* m_RelativePositions sont les bras de levier                    */
  m_RelativePositions = ParticuleRef->getRelativePositions();

  m_InitialMatrix = ParticuleRef->getInitialMatrix();

  /* Affectation de la cinematique*/
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
  {
    m_elementaryParticules[i]->setVitesseTranslation( Vecteur( vx, vy, vz )
      + ( Vecteur( rx, ry, rz ) ^ ( (m_geoFormeVdw->getTransform()->getBasis())
      * m_RelativePositions[i] ) ) );
    m_elementaryParticules[i]->setVitesseRotation( Vecteur( rx, ry, rz ) );

    m_elementaryParticules[i]->setQuaternionRotation( qrotationx, qrotationy,
      qrotationz, qrotations );
  }

  if ( updatePosition )
    setElementPosition();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur ( SEQUENTIEL PERIODIQUE )
CompParticule::CompParticule( const int &id_, Particule const* ParticuleRef,
	const Vecteur &vtrans,
	const Quaternion &qrot,
	const Vecteur &vrot,
	const Transform &config,
	const ParticuleActivity &activ,
	const int &tag_ ) :
  Particule( id_, ParticuleRef,
	*(ParticuleRef->getVitesseTranslation()),
	*(ParticuleRef->getCinematique()->getRotation()),
	*(ParticuleRef->getVitesseRotation()),
	*(ParticuleRef->getForme()->getTransform()),
	COMPUTE,
	0 )
{
  /* Particules elementaires */
  m_elementaryParticules.reserve( ParticuleRef->getNbreElemPart() );
  ElementParticule* ppp = NULL;
  for ( size_t i=0; i<ParticuleRef->getNbreElemPart(); ++i )
    m_elementaryParticules.push_back( ppp );

  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    m_elementaryParticules[i] = new ElementParticule(
      *(ParticuleRef->getElementParticules()[i]), this );

  /* Positions des particules elementaires par rapport a (0.,0.,0.) */
  /* m_RelativePositions sont les bras de levier                    */
  m_RelativePositions = ParticuleRef->getRelativePositions();

  /* Orientations initiales des particules elementaires */
  m_InitialMatrix = ParticuleRef->getInitialMatrix();

  /* Mise � jour des particules elementaires */
  setElementPosition();

  /* Affectation de la cinematique*/
  Matrix rota = m_geoFormeVdw->getTransform()->getBasis() ;
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
  {
    m_elementaryParticules[i]->setVitesseTranslation( vtrans
      + ( vrot ^ ( rota * m_RelativePositions[i] ) ) );
    m_elementaryParticules[i]->setVitesseRotation( vrot  );

    m_elementaryParticules[i]->setQuaternionRotation( qrot );
  }
}




// ----------------------------------------------------------------------------
// Destructeur
CompParticule::~CompParticule()
{
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    delete m_elementaryParticules[i];
  m_elementaryParticules.clear();
}




// ----------------------------------------------------------------------------
// Constructeur par copie
CompParticule::CompParticule( const CompParticule &copie ) :
  Particule( copie )
{
  m_RelativePositions.reserve( copie.m_RelativePositions.size() );
  for ( size_t i=0; i<copie.m_RelativePositions.size(); ++i )
    m_RelativePositions.push_back( copie.m_RelativePositions[i] );

  m_CG = copie.m_CG;

  m_name = copie.m_name;

  m_coordination_number = copie.m_coordination_number;

  /* !!! L'ideal serait de stocker toutes ces donnees dans la particule
  de reference et non avoir les copies dans chaque instance !!! */
  m_InitialMatrix.reserve( copie.m_InitialMatrix.size() );
  for ( size_t i=0; i<copie.m_InitialMatrix.size(); ++i )
    m_InitialMatrix.push_back( copie.m_InitialMatrix[i] );

  m_InitialRelativePositions.reserve( copie.m_InitialRelativePositions.size() );
  for ( size_t i=0; i<copie.m_InitialRelativePositions.size(); ++i )
    m_InitialRelativePositions.push_back( copie.m_InitialRelativePositions[i] );

  ElementParticule* ppp = NULL ;
  m_elementaryParticules.reserve( copie.m_elementaryParticules.size() );
  for ( size_t i=0; i<copie.m_elementaryParticules.size(); ++i )
  {
    ppp = new ElementParticule( *(copie.m_elementaryParticules[i]),
    	this ) ;
    m_elementaryParticules.push_back( ppp );
  }

  m_compInertie.reserve( copie.m_compInertie.size() );
  for ( size_t i=0; i<copie.m_compInertie.size(); ++i )
    m_compInertie.push_back( copie.m_compInertie[i] );

  m_compVolume = copie.m_compVolume;

  m_nbreElemPart = copie.m_nbreElemPart;

  Composant::setNbComposantsCrees( Composant::getNbComposantsCrees() -
  	int(copie.m_elementaryParticules.size()) );
}




// ----------------------------------------------------------------------------
// Constructeur d'un clone par copie
Particule* CompParticule::createCloneCopy() const
{
  Particule* particule = new CompParticule( *this );
  return particule;
}




// ----------------------------------------------------------------------------
// Determine le rayon circonscrit du composite
void CompParticule::setRayon()
{
  /* R_{comp} = max( abs(GG_i)+R_i ) */

  /* Contient tous les abs(GG_i) + R_i de la particule composite */
  vector<Scalar> GGi_Ri;
  GGi_Ri.reserve( m_elementaryParticules.size() );

  Scalar RayonComposite = 0.;

  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
  {
    GGi_Ri[i] = Norm( m_RelativePositions[i] )
	      + m_elementaryParticules[i]->getRayon();

    if ( GGi_Ri[i] > RayonComposite )
      RayonComposite = GGi_Ri[i];
  }
  m_geoFormeVdw->setRayon( RayonComposite );
}




// ----------------------------------------------------------------------------
// Renvoie le rayon circonscrit du composite
vector<Scalar> CompParticule::getCompFeaturesBox()
{
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    m_elementaryParticules[i]->setPosition( m_InitialRelativePositions[i] );

  /*
  v[0] = x_max, v[1] = x_min
  v[2] = y_max, v[3] = y_min
  v[4] = z_max, v[5] = z_min
  */

  vector<Scalar> allXmax, allYmax, allZmax;
  vector<Scalar> allXmin, allYmin, allZmin;

  allXmax.reserve( m_elementaryParticules.size() );
  allYmax.reserve( m_elementaryParticules.size() );
  allZmax.reserve( m_elementaryParticules.size() );

  allXmin.reserve( m_elementaryParticules.size() );
  allYmin.reserve( m_elementaryParticules.size() );
  allZmin.reserve( m_elementaryParticules.size() );

  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
  {
    allXmax.push_back( m_InitialRelativePositions[i][X] +
	m_elementaryParticules[i]->getRayon() );
    allYmax.push_back( m_InitialRelativePositions[i][Y] +
	m_elementaryParticules[i]->getRayon() );
    allZmax.push_back( m_InitialRelativePositions[i][Z] +
	m_elementaryParticules[i]->getRayon() );

    allXmin.push_back( m_InitialRelativePositions[i][X] -
	m_elementaryParticules[i]->getRayon() );
    allYmin.push_back( m_InitialRelativePositions[i][Y] -
	m_elementaryParticules[i]->getRayon() );
    allZmin.push_back( m_InitialRelativePositions[i][Z] -
	m_elementaryParticules[i]->getRayon() );

  }

  vector<Scalar> box;
  for (int i=0; i<6; ++i )
    box.push_back(0.);

  box[0] = *max_element( allXmax.begin(), allXmax.end() ); // x_max
  box[1] = *min_element( allXmin.begin(), allXmin.end() ); // x_min

  box[2] = *max_element( allYmax.begin(), allYmax.end() ); // y_max
  box[3] = *min_element( allYmin.begin(), allYmin.end() ); // y_min

  box[4] = *max_element( allZmax.begin(), allZmax.end() ); // z_max
  box[5] = *min_element( allZmin.begin(), allZmin.end() ); // z_min


  return( box );
}




// ----------------------------------------------------------------------------
// Renvoie le rayon circonscrit du composite
Scalar CompParticule::getRayon() const
{
  return( m_geoFormeVdw->getRayon() );
}




// ----------------------------------------------------------------------------
// Positionne les particules elementaires dans l'espace
void CompParticule::setElementPosition()
{
  /* ----------------------------------------------------------------------- */
  /* Positionne les particules elementaires en fonction de la transformation */
  /* de la particule "maitresse". En effet, les particules elementaires sont */
  /* positionnees en tenant compte de la position de la particle "maitresse" */
  /* en ajoutant ensuite les bras de levier de chaque particule elementaire. */
  /* La rotation de la particule elementaire est deduite de celle qui est    */
  /* "maitresse" en faisant une composition de rotations:                    */
  /* R(total) = R(element_initial) o R(maitresee)                            */
  /* Les particules elementaires subissent une translation (rigid body       */
  /* motion): ( R(Total) * levier ) - levier                                 */
  /* ----------------------------------------------------------------------- */
  Vecteur trans;
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
  {
    m_elementaryParticules[i]->setPosition( (*m_geoFormeVdw->getCentre()
      + m_RelativePositions[i]) );

    m_elementaryParticules[i]->getForme()->getTransform()->setBasis(
   ( ( m_geoFormeVdw->getTransform()->getBasis() ) * ( m_InitialMatrix[i] ) ) );

    trans  = ( ( ( m_geoFormeVdw->getTransform()->getBasis() ) )
	    * m_RelativePositions[i] ) - m_RelativePositions[i];
    m_elementaryParticules[i]->Translate( trans );
  }
}




// ----------------------------------------------------------------------------
// Positionne la particule composite
void CompParticule::setPosition( const Point &centre )
{
  Particule::setPosition( centre );
  setElementPosition();
}




// ----------------------------------------------------------------------------
// Positionne la particule composite
void CompParticule::setPosition( const Scalar *pos )
{
  Particule::setPosition( pos );
  setElementPosition();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation d'un champ de vitesse de translation
void CompParticule::setVitesseTranslation( const Vecteur &translation )
{
  Particule::setVitesseTranslation( translation );
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    m_elementaryParticules[i]->setVitesseTranslation( *m_cinematique->
	getVitesseTranslation() + ( *m_cinematique->getVitesseRotation()
	^ ( m_geoFormeVdw->getTransform()->getBasis()
	* m_RelativePositions[i] ) ) );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Affectation d'un champ de vitesse de rotation
void CompParticule::setVitesseRotation( const Vecteur &rotation )
{
  Particule::setVitesseRotation( rotation );
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    m_elementaryParticules[i]->setVitesseRotation( rotation );
}




// ----------------------------------------------------------------------------
// Quaternion de rotation dans le vecteur vit en d�butant � la position i
void CompParticule::copyQuaternionRotation( double *vit, int i ) const
{
  Particule::copyQuaternionRotation( vit, i );
}




// ----------------------------------------------------------------------------
// Vitesse de rotation dans le vecteur vit en d�butant � la position i
void CompParticule::copyVitesseRotation( double *vit, int i ) const
{
  Particule::copyVitesseRotation( vit, i );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Vitesse de translation dans le vecteur vit en d�butant � la position i
void CompParticule::copyVitesseTranslation( double *vit, int i ) const
{
  Particule::copyVitesseTranslation( vit, i );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copie de la cinematique au temps t-2dt: vitesse de translation,
// vitese de rotation, variation de QDM translationnelle, variation de QDM
// rotationnalle
void CompParticule::copyCinematiqueNm2( double *vit, int i ) const
{
  Particule::copyCinematiqueNm2( vit, i );
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copie de la transformation dans le vecteur vit
// en d�butant � la position i
void CompParticule::copyTransform( double *vit, int i ) const
{
  Particule::copyTransform( vit, i );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copie de la transformation dans le vecteur vit
// en d�butant � la position i, avec une translation du centre de gravite
void CompParticule::copyTransform( double *vit, int i, Vecteur const& vec ) const
{
  Particule::copyTransform( vit, i, vec );
}




// ----------------------------------------------------------------------------
// Modification du quaternion de rotation.
void CompParticule::setQuaternionRotation( const Scalar &vecteur0,
	const Scalar &vecteur1, const Scalar &vecteur2,
	const Scalar &scalaire )
{
  Particule::setQuaternionRotation( vecteur0, vecteur1, vecteur2, scalaire );
}




// ----------------------------------------------------------------------------
// Modification du quaternion de rotation.
void CompParticule::setQuaternionRotation( const Quaternion &qrot )
{
  Particule::setQuaternionRotation( qrot );
}




// ----------------------------------------------------------------------------
// Renvoie le volume de la particule composite
Scalar CompParticule::getVolume() const
{
  return ( m_compVolume );
}




// ----------------------------------------------------------------------------
// Calcule le poids de la particule composite
void CompParticule::computeWeight()
{
  m_weight = m_masseVolumique * getVolume()
         * ( 1. - Particule::m_fluideMasseVolumique / m_masseVolumique )
	 * Grains_Exec::m_vgravite ;
}




// ----------------------------------------------------------------------------
// Ecrit le convexe pour post-processing avec Paraview
void CompParticule::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    m_elementaryParticules[i]->write_polygonsStr_PARAVIEW(connectivity,
	    offsets, cellstype, firstpoint_globalnumber, last_offset);
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
void CompParticule::write_polygonsPts_PARAVIEW( ostream &f,
	const Transform &transform, Vecteur const* translation ) const
{
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    m_elementaryParticules[i]->write_polygonsPts_PARAVIEW(f,
	transform, translation);
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
int CompParticule::numberOfPoints_PARAVIEW() const
{
  int nbpts = 0 ;
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    nbpts += m_elementaryParticules[i]->getForme()->getConvex()
	    ->numberOfPoints_PARAVIEW();

  return ( nbpts );
}




// ----------------------------------------------------------------------------
// Nombre de polygones elementaires pour post-processing avec Paraview
int CompParticule::numberOfCells_PARAVIEW() const
{
  int nbcells = 0 ;
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    nbcells += m_elementaryParticules[i]->getForme()->getConvex()
	    ->numberOfCells_PARAVIEW();

  return ( nbcells );
}




// ----------------------------------------------------------------------------
// Renvoie les points du convexe pour post-processing avec Paraview (binary)
list<Point> CompParticule::get_polygonsPts_PARAVIEW(
	Vecteur const* translation ) const
{
  list<Point> ppp, ParaviewPoints;
  list<Point>::const_iterator itpp;

  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
  {
    ppp = m_elementaryParticules[i]->get_polygonsPts_PARAVIEW( translation );
    for ( itpp=ppp.begin(); itpp!=ppp.end(); ++itpp )
      ParaviewPoints.push_back( *itpp );
  }

  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
void CompParticule::write_polygonsPts_PARAVIEW( ostream &f,
	Vecteur const* translation )const
{
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    m_elementaryParticules[i]->write_polygonsPts_PARAVIEW( f, translation ) ;
}




// ----------------------------------------------------------------------------
// Applique la transformation trot � chaque particule elementaire
// dans le cas d'une configuration aletoire
void CompParticule::composePosition( const Transform &trot )
{
  m_geoFormeVdw->composeRotationRight( trot );
  CompParticule::setElementPosition();
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils proches ? Variante avec les rayons de VdW.
bool CompParticule::isProcheVdW( const Composant* voisin ) const
{
  bool contact = false;

  for ( size_t i=0; i<m_elementaryParticules.size() && !contact; ++i )
    if ( voisin->isCompParticule() )
      contact = voisin->isProcheVdW( m_elementaryParticules[i] );
    else
      contact = m_elementaryParticules[i]->isProcheVdW( voisin );

  return contact;
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils proches ?
bool CompParticule::isProche( const Composant* voisin ) const
{
  bool contact = false;

  for ( size_t i=0; i<m_elementaryParticules.size() && !contact; ++i )
    if ( voisin->isCompParticule() )
      contact = voisin->isProche( m_elementaryParticules[i] );
    else
      contact = m_elementaryParticules[i]->isProche( voisin );

  return contact;
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils en contact ?
bool CompParticule::isContact( const Composant* voisin ) const
{
  bool contact = false;

  for ( size_t i=0; i<m_elementaryParticules.size() && !contact; ++i )
    if ( voisin->isCompParticule() )
      contact = voisin->isContact( m_elementaryParticules[i] );
    else
      contact = m_elementaryParticules[i]->isContact( voisin );

  return contact;
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils en contact ?
bool CompParticule::isContactVdW( const Composant* voisin ) const
{
  bool contact = false;

  for ( size_t i=0; i<m_elementaryParticules.size() && !contact; ++i )
    if ( voisin->isCompParticule() )
      contact = voisin->isContactVdW( m_elementaryParticules[i] );
    else
      contact = m_elementaryParticules[i]->isContactVdW( voisin );

  return contact;
}




// ----------------------------------------------------------------------------
// Contact entre particule composite et une simple particule ou une particule
// composite.
void CompParticule::InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC )
  throw ( ErreurContact )
{
  list<ContactInfos*>  listContactInfos;
	// int c=0;
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
  {
    if ( voisin->isCompParticule() )
      voisin->SearchContact( m_elementaryParticules[i], dt,
	  temps, LC, listContactInfos );
    else m_elementaryParticules[i]->SearchContact( voisin, dt,
	  temps, LC, listContactInfos )  ;
		// cout << "sub-particle #"<<c<<" and length ContactInfo="<<int(listContactInfos.size()) << endl;
		// c++;
  }

  int nbContact = int(listContactInfos.size());

  for ( list<ContactInfos*>::iterator il=listContactInfos.begin();
      il!=listContactInfos.end(); il++ )
  {
    LC->addToContactsFeatures( temps, (*il)->ContactPoint );

    if ( Contact_BuilderFactory::contactForceModel(
          	(*il)->p0->materiau(), (*il)->p1->materiau() )
      		->computeForces( (*il)->p0, (*il)->p1, (*il)->ContactPoint,
		LC, dt, nbContact ) )
    {
      (*il)->p0->addToCoordinationNumber( 1 );
      (*il)->p1->addToCoordinationNumber( 1 );
    }
    delete *il;
  }
  listContactInfos.clear();
}




// ----------------------------------------------------------------------------
// Contact entre particule composite et une simple particule ou une particule
// composite.
void CompParticule::SearchContact( Composant* voisin, double dt,
      double const& temps, LinkedCell *LC,
      list<ContactInfos*> &listContact )
{
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
      m_elementaryParticules[i]->SearchContact( voisin, dt, temps, LC ,
	listContact );
}




// ----------------------------------------------------------------------------
// Determination du contact avec une particule composite
void CompParticule::InterActionPostProcessing( Composant* voisin,  Scalar dt,
	list<struct PointForcePostProcessing>* listOfContacts )
  throw ( ErreurContact )
{
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    if ( voisin->isCompParticule() )
      voisin->InterActionPostProcessing( m_elementaryParticules[i], dt,
	  listOfContacts );
    else
      m_elementaryParticules[i]->InterActionPostProcessing( voisin, dt,
	  listOfContacts );
}




// ----------------------------------------------------------------------------
// Resolution des equations de la dynamique et double integration pour
// obtenir la nouvelle vitesse et position.
void CompParticule::Deplacer( Scalar temps, double dt )
    throw( ErreurDeplacement )
{
  if (is_mobile)
  Particule::Deplacer( temps, dt );
  setElementPosition();

  /* -------------------------------------------------- */
  /* U_{elem} = U_{comp} + omega ^ ( M(rota) * levier ) */
  /* -------------------------------------------------- */
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
  {
    m_elementaryParticules[i]->setVitesseTranslation( *m_cinematique->
	getVitesseTranslation() + ( *m_cinematique->getVitesseRotation()
	^ ( m_geoFormeVdw->getTransform()->getBasis()
	* m_RelativePositions[i] ) ) );

    m_elementaryParticules[i]->setVitesseRotation( *m_cinematique->
	getVitesseRotation() );
  }
}




// ----------------------------------------------------------------------------
// Determine l'inertie et l'inertie inverse de la particule composite
bool  CompParticule::BuildInertie(Scalar *inertie, Scalar *inertie_1) const
{
  for ( size_t i=0; i<6; ++i )
    inertie[i] = m_compInertie[i] * m_masseVolumique;

  Scalar determinant = inertie[0]*inertie[3]*inertie[5]
    - inertie[0]*inertie[4]*inertie[4]
    - inertie[5]*inertie[1]*inertie[1]
    - inertie[3]*inertie[2]*inertie[2]
    + 2*inertie[1]*inertie[2]*inertie[4];

  inertie_1[0] = (inertie[3]*inertie[5]-inertie[4]*inertie[4])/determinant;
  inertie_1[1] = (inertie[2]*inertie[4]-inertie[1]*inertie[5])/determinant;
  inertie_1[2] = (inertie[1]*inertie[4]-inertie[2]*inertie[3])/determinant;
  inertie_1[3] = (inertie[0]*inertie[5]-inertie[2]*inertie[2])/determinant;
  inertie_1[4] = (inertie[1]*inertie[2]-inertie[0]*inertie[4])/determinant;
  inertie_1[5] = (inertie[0]*inertie[3]-inertie[1]*inertie[1])/determinant;

  return true;
}




// ----------------------------------------------------------------------------
// Initialize a faux le boolean correspondant au calcul de la
// transformation avec scaling par l'epaisseur de croute
void CompParticule::initializeVdWtransform_to_notComputed()
{
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
    m_elementaryParticules[i]->initializeVdWtransform_to_notComputed();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Determine s'il y a intersection entre une Forme et un Composant
// Vrai s'il y a intersection
bool CompParticule::intersectWithForme( const Forme &b )
{
  bool inter = false ;
  for ( size_t i=0; i<m_elementaryParticules.size() && !inter; ++i )
    inter = m_elementaryParticules[i]->intersectWithForme( b );

  return ( inter ) ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie un vecteur de particules elementaires
vector<ElementParticule*> CompParticule::getElementParticules() const
{
  return( m_elementaryParticules );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les positions initiales des particules elementaires a l'insertion
vector<Vecteur> CompParticule::getInitialRelativePositions() const
{
  return ( m_InitialRelativePositions );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les positions de la particule composite apres mise a jour
// du centre de gravite
vector<Vecteur> CompParticule::getRelativePositions() const
{
  return ( m_RelativePositions );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie les matrices de rotations intiales des particules elementaires
vector<Matrix> CompParticule::getInitialMatrix() const
{
  return ( m_InitialMatrix );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre total de particules elementaires
size_t CompParticule::getNbreElemPart() const
{
  return ( m_nbreElemPart );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nom associe a la particule composite
string CompParticule::getPartName() const
{
  return ( m_name );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde du Particule pour Reload
void CompParticule::write( ostream &fileSave, Composant const* composant ) const
{
  Particule::write( fileSave, this );
  fileSave << "\n<ElementParticules>\t"
	   << m_elementaryParticules.size() << "\n";
  for ( size_t i=0; i<m_elementaryParticules.size(); ++i )
  {
    fileSave << "\n<ElementParticule>\n";
    fileSave << m_elementaryParticules[i] << '\t'
	     << m_elementaryParticules[i]->getID() << '\n';
    m_elementaryParticules[i]->getForme()->writeStatique( fileSave );
    writePosition( fileSave, m_RelativePositions[i], m_InitialMatrix[i],
      m_elementaryParticules[i] );
//    Composant::writePosition( fileSave, m_elementaryParticules[i] );
    fileSave << "\n</ElementParticule>\n";
  }
  fileSave << "\n</ElementParticules>\n";
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de l'information de position
void CompParticule::writePosition( ostream &position, const Vecteur
  &relPos, const Matrix &initRot, ElementParticule* particule ) const
{
  Scalar buf = 0. ;
  position << "*Couleur\n"
	   << buf << '\t'
	   << buf << '\t'
	   << buf << '\n';

  position << "%Position&Orientation\n";
  position << "*Position\n";
  position << Grains_Exec::doubleToString( ios::scientific,POSITIONFORMAT,
    relPos[X] ) << " " <<
	      Grains_Exec::doubleToString( ios::scientific,POSITIONFORMAT,
    relPos[Y] ) << " " <<
	      Grains_Exec::doubleToString( ios::scientific,POSITIONFORMAT,
    relPos[Z] ) << endl;
  position << "*Orientation\n";
  initRot.writeMatrix( position );
  position <<  endl;
  position << "*Type\n";
  position << particule->getForme()->getTransform()->getType();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture de la particule pour Reload. Utilise pour les particules
// dans la simulation (actives ou en attentes) pour le format de reload 2014
// format de reload
void CompParticule::read2014( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference )
{
  // Lecture du n� et de la classe de reference
  fileSave >> m_id >> m_ParticuleClasse;

  // Creation de la forme (convex + transform) en utilisant le constructeur
  // de recopie
  m_geoFormeVdw = new FormeVdW(
    *(*ParticuleClassesReference)[m_ParticuleClasse]->getForme() );

  // Materiau, masse, energie et inertie de la particule de class de reference
  m_nomMateriau =
  	(*ParticuleClassesReference)[m_ParticuleClasse]->materiau() ;
  m_masse = (*ParticuleClassesReference)[m_ParticuleClasse]->getMasse() ;
  m_energie = (*ParticuleClassesReference)[m_ParticuleClasse]->getEnergie() ;
  is_mobile = (*ParticuleClassesReference)[m_ParticuleClasse]->getMobilite() ;
  for (size_t i=0;i<6;++i)
  {
    m_inertie[i] =
    	(*ParticuleClassesReference)[m_ParticuleClasse]->getInertie()[i] ;
    m_inertie_1[i] =
    	(*ParticuleClassesReference)[m_ParticuleClasse]->getInertieInverse()[i] ;
  }

  // Lecture du tag
  fileSave >> m_tag;

  // Lecture de la transformation
  m_geoFormeVdw->readPosition2014( fileSave );

  // Lecture de l'activite
  bool actif;
  fileSave >> actif;
  m_activity = ( actif == true ) ? COMPUTE : WAIT;

  // Lecture de la cinematique
  if ( m_cinematique ) delete m_cinematique;
  m_cinematique = Cinematique_BuilderFactory::create(
  	m_geoFormeVdw->getConvex() );
  fileSave >> *m_cinematique;

	// Read the contact memories of the particle (if any)
  readContactMap_2014( fileSave );

  // Calcul de la masse volumique et du poids
  m_compVolume = (*ParticuleClassesReference)[m_ParticuleClasse]->getVolume();
  m_masseVolumique = m_masse / m_compVolume;
  computeWeight();

  // Mise a jour des particules elementaires
  m_nbreElemPart = (*ParticuleClassesReference)[m_ParticuleClasse]
      ->getNbreElemPart();
  /* Construction des particules elementaires */
  m_elementaryParticules.reserve( m_nbreElemPart );
  ElementParticule* ppp = NULL;
  for ( size_t i=0; i<m_nbreElemPart; ++i )
    m_elementaryParticules.push_back( ppp );

  for ( size_t i=0; i<m_nbreElemPart; ++i )
    m_elementaryParticules[i] = new ElementParticule(
      *((*ParticuleClassesReference)[m_ParticuleClasse]
	  ->getElementParticules()[i]), this );

  m_RelativePositions = (*ParticuleClassesReference)[m_ParticuleClasse]
      ->getRelativePositions();
  m_InitialMatrix = (*ParticuleClassesReference)[m_ParticuleClasse]
      ->getInitialMatrix();
  setElementPosition();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture de la particule pour Reload. Utilise pour les particules
// dans la simulation (actives ou en attentes) pour le format de reload 2014
// format de reload
void CompParticule::read2014_binary( istream &fileSave,
	vector<Particule*> const* ParticuleClassesReference,
	bool const& contact_history_storage )
{
  // Lecture du n� et de la classe de reference
	fileSave.read( reinterpret_cast<char*>( &m_id ), sizeof(int) );
	fileSave.read( reinterpret_cast<char*>( &m_ParticuleClasse ), sizeof(int) );

  // Creation de la forme (convex + transform) en utilisant le constructeur
  // de recopie
  m_geoFormeVdw = new FormeVdW(
    *(*ParticuleClassesReference)[m_ParticuleClasse]->getForme() );

  // Materiau, masse, energie et inertie de la particule de class de reference
  m_nomMateriau =
  	(*ParticuleClassesReference)[m_ParticuleClasse]->materiau() ;
  m_masse = (*ParticuleClassesReference)[m_ParticuleClasse]->getMasse() ;
  m_energie = (*ParticuleClassesReference)[m_ParticuleClasse]->getEnergie() ;
  is_mobile = (*ParticuleClassesReference)[m_ParticuleClasse]->getMobilite() ;

  for (size_t i=0;i<6;++i)
  {
    m_inertie[i] =
    	(*ParticuleClassesReference)[m_ParticuleClasse]->getInertie()[i] ;
    m_inertie_1[i] =
    	(*ParticuleClassesReference)[m_ParticuleClasse]->getInertieInverse()[i] ;
  }

  // Lecture du tag
  fileSave.read( reinterpret_cast<char*>( &m_tag ), sizeof(int) );

  // Lecture de la transformation
  m_geoFormeVdw->readPosition2014_binary( fileSave );

  // Lecture de l'activite
  unsigned int iact = 0;
  fileSave.read( reinterpret_cast<char*>( &iact ), sizeof(unsigned int) );
  if ( iact == 1 ) m_activity = COMPUTE;
  else m_activity = WAIT;

  // Lecture de la cinematique
  if ( m_cinematique ) delete m_cinematique;
  m_cinematique = Cinematique_BuilderFactory::create(
  	m_geoFormeVdw->getConvex() );
  m_cinematique->readCineParticule2014_binary( fileSave );

	if (contact_history_storage)
  {
    // Read the contact memories of the particle (if any)
    readContactMap_binary( fileSave );
  }

  // Calcul de la masse volumique et du poids
  m_compVolume = (*ParticuleClassesReference)[m_ParticuleClasse]->getVolume();
  m_masseVolumique = m_masse / m_compVolume;
  computeWeight();

  // Mise a jour des particules elementaires
  m_nbreElemPart = (*ParticuleClassesReference)[m_ParticuleClasse]
      ->getNbreElemPart();
  /* Construction des particules elementaires */
  m_elementaryParticules.reserve( m_nbreElemPart );
  ElementParticule* ppp = NULL;
  for ( size_t i=0; i<m_nbreElemPart; ++i )
    m_elementaryParticules.push_back( ppp );

  for ( size_t i=0; i<m_nbreElemPart; ++i )
    m_elementaryParticules[i] = new ElementParticule(
      *((*ParticuleClassesReference)[m_ParticuleClasse]
	  ->getElementParticules()[i]), this );

  m_RelativePositions = (*ParticuleClassesReference)[m_ParticuleClasse]
      ->getRelativePositions();
  m_InitialMatrix = (*ParticuleClassesReference)[m_ParticuleClasse]
      ->getInitialMatrix();
  setElementPosition();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture de la particule pour Reload. Utilise pour les particules de
// classe de reference dans l'entete du fichier de reload et pour l'ancien
// format de reload
void CompParticule::read( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference )
{
  string buffer, adresse, particuleType, Cle, compCle;

  fileSave >> buffer
	   >> adresse >> m_id;
  SaveTable::create( adresse, this );

  fileSave >> buffer
	   >> m_nomMateriau;

  streampos length = fileSave.tellg();
  fileSave >> compCle >> m_compVolume >> m_name;
  fileSave.seekg( length );

  if ( !ParticuleClassesReference )
  {
    m_geoFormeVdw = new FormeVdW( fileSave );
  }
  else
  {
    string cle;
    fileSave >> cle;
    while ( cle != "*END" ) fileSave >> cle;
  }

  fileSave >> buffer
	   >> m_ParticuleClasse;

  if ( ParticuleClassesReference )
  {
    m_geoFormeVdw = new FormeVdW( new PointC(), Transform() );
  }
  fileSave >> buffer
	   >> m_tag;
  fileSave >> buffer
	   >> m_masse >> m_energie;

  // We check if the mobility is precised in the text reload file then we read
  // it if it exist, otherwhise we put back the input stream to the initial
  // position ( to be conformant with both old and new reload formats (2016))
  length = fileSave.tellg();
  string is_mobilite_precised;
  fileSave >> is_mobilite_precised;
  if ( is_mobilite_precised=="*Mobilite" )
    fileSave >> is_mobile;
  else
    fileSave.seekg( length );

   // We check if the mobility is precised in the text reload file then we read
   // it if it exist, otherwhise we put back the input stream to the initial
   // position ( to be conformant with both old and new reload formats (2016))
  length = fileSave.tellg();
  string is_temperature_evolution_precised;
  fileSave >> is_temperature_evolution_precised;
  if ( is_temperature_evolution_precised=="*TempEvolution" )
    fileSave >> m_temp_evolution;
  else
   fileSave.seekg( length );

  fileSave >> buffer
	   >> m_inertie[0] >> m_inertie[1] >> m_inertie[2]
	   >> m_inertie[3] >> m_inertie[4] >> m_inertie[5];
  fileSave >> buffer
	   >> m_inertie_1[0] >> m_inertie_1[1] >> m_inertie_1[2]
	   >> m_inertie_1[3] >> m_inertie_1[4] >> m_inertie_1[5];

  int buf = 0;
  fileSave >> buffer >> buf >> buf >> buf;
  m_geoFormeVdw->readPosition( fileSave );

  bool actif;
  fileSave >> buffer >> actif;
  // cout << buffer << " " << actif << endl;
  m_activity = ( actif == true ) ? COMPUTE : WAIT;

  if ( m_cinematique ) delete m_cinematique;
  m_cinematique = Cinematique_BuilderFactory::read( fileSave,
  	m_geoFormeVdw->getConvex() );

  m_masseVolumique = m_masse / m_compVolume ;
  computeWeight();

  // Traitement des particules elementaires
  fileSave >> Cle >> m_nbreElemPart;
  ReadElementParticules( fileSave, m_nbreElemPart, m_cinematique );
  setRayon();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Gestion de la lecture des particules elementaires
void CompParticule::ReadElementParticules( istream &fileSave, size_t nbre,
    CineParticule* cinematique )
{
  Matrix ttt;
  ttt.setIdentity();
  Vecteur Nul(0.,0.,0.);

  m_elementaryParticules.reserve( nbre );
  ElementParticule* ppp = NULL;
  for ( size_t i=0; i<nbre; ++i )
    m_elementaryParticules.push_back( ppp );
  m_RelativePositions.reserve( nbre );
  m_InitialMatrix.reserve( nbre );
  m_InitialRelativePositions.reserve( nbre );

  for ( size_t i=0; i<nbre; ++i )
  {
    m_InitialMatrix.push_back( ttt );
    m_RelativePositions.push_back( Nul );
    m_InitialRelativePositions.push_back( Nul );
  }

  string Cle;
	string prev_Cle;
  size_t j = 0;
  while ( Cle != "</ElementParticules>" )
  {
    if ( Cle == "*Cylindre" || Cle == "*Polyhedron"
	|| Cle == "*Sphere" || Cle == "*Box" )
    {
      j++;
      // Creation des particules elementaires
      m_elementaryParticules[j-1] = new ElementParticule( fileSave,
	    Cle, this, std::stoi(prev_Cle) );
    }
		prev_Cle = Cle;
    fileSave >> Cle;
  }
  assert( Cle == "</ElementParticules>" );


  // Mise a jour des particules elementaires
  for ( size_t i=0; i<nbre; ++i )
  {
    m_RelativePositions[i] = *m_elementaryParticules[i]->getPosition();
    m_InitialMatrix[i] = m_elementaryParticules[i]->getForme()
    ->getTransform()->getBasis();
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de sommets ou un code �quivalent
// D. RAKOTONIRINA - Avril. 2015 - Creation
int CompParticule::getNbCorners() const
{
  return 999;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de la position de la particule pour le Fluide
// D. RAKOTONIRINA - Avril. 2015 - Creation
void CompParticule::writePositionInFluid( ostream &fileOut )
{
  setElementPosition();
  Vecteur vitesseT;
  Transform transform;

  /* Donnees des particules elemntaires */
  size_t nbElt = m_elementaryParticules.size();
  size_t nbCyl = nbElt - 1;
  Point pointEnvelop;
  vector<Point> allPointsCyl, allPointsPoly;
  ElementParticule* ppp = NULL;
  ElementParticule* polyhedron = NULL;

  vector<ElementParticule*> AllCylinders(nbCyl, ppp);

  Point pointNul( 0.,0.,0. );
  vector<vector<Point> > RefCylEnv( nbCyl, vector<Point>( 3, pointNul ) );
  vector<Point> RefPolyEnv;
  vector<Point> BottomCenter( nbCyl, pointNul );
  vector<Point> TopCenter( nbCyl, pointNul );
  vector<Scalar> CylRadius( nbCyl, 0.);

  for ( size_t i=0; i<nbElt; ++i )
  {
    if ( m_elementaryParticules[i]->getConvexType() == CYLINDER )
    {
      AllCylinders[i] = m_elementaryParticules[i];
      CylRadius[i] = m_elementaryParticules[i]->getEnveloppe()[1][X];
      RefCylEnv[i] = m_elementaryParticules[i]->getEnveloppe();
      BottomCenter[i] = RefCylEnv[i][0] + m_RelativePositions[i];
      TopCenter[i] = RefCylEnv[i][2] + m_RelativePositions[i];
    }
    else
    {
      polyhedron = m_elementaryParticules[i];
      RefPolyEnv =  m_elementaryParticules[i]->getEnveloppe();
      for ( size_t j=0; j<RefPolyEnv.size(); ++j )
	RefPolyEnv[j] += m_RelativePositions[i];
    }
  }

  size_t polySize = RefPolyEnv.size();
  size_t halfPolySize = size_t( polySize/2 );

  vector<Point> PolyTop( halfPolySize, pointNul ),
    PolyBottom = PolyTop;

  for ( size_t i=0; i<polySize; ++i )
    if ( i<halfPolySize ) PolyBottom[i] = RefPolyEnv[i];
    else PolyTop[i-halfPolySize] = RefPolyEnv[i];


  // Recherche du sommet du polygone qui se trouve sur la surface
  // du disque du cylindre
  Vecteur distance;
  vector<Point> CylPolyPtsBottom( halfPolySize, pointNul ),
    CylPolyPtsTop = CylPolyPtsBottom;
  for ( size_t i=0; i<nbCyl; ++i )
    for ( size_t j=0; j<halfPolySize; ++j )
    {
       distance = - BottomCenter[i] + PolyBottom[j];
       if ( Norm( distance ) <= 1.01*CylRadius[i] )
       {
	  CylPolyPtsBottom[i] = PolyBottom[j];
//	  cout << "BOT = " << CylPolyPtsBottom[i] << endl;
       }
    }

  for ( size_t i=0; i<nbCyl; ++i )
    for ( size_t j=0; j<halfPolySize; ++j )
    {
       distance = - TopCenter[i] + PolyTop[j];
       if ( Norm( distance ) <= 1.01*CylRadius[i] )
       {
	  CylPolyPtsTop[i] = PolyTop[j];
//	  cout << "TOP = " << CylPolyPtsTop[i] << endl;
       }
    }

  // A chaque disque du cylindre sont associes 3 sommets du polygone,
  // un se trouvant sur la surface du cylindre et les deux autres dans les
  // cylindres voisins
  vector<struct tripletPts> TopPts( nbCyl ), BottomPts( nbCyl );
  // *** Trilobe ***
  if ( nbCyl == 3 )
  {
    for ( size_t i=0; i<nbCyl; ++i )
    {
      if ( i == nbCyl-1 )
      {
        BottomPts[i].Center = CylPolyPtsBottom[i];
        BottomPts[i].Poly1  = CylPolyPtsBottom[0];
        BottomPts[i].Poly2  = CylPolyPtsBottom[i-1];

        TopPts[i].Center = CylPolyPtsTop[i];
        TopPts[i].Poly1  = CylPolyPtsTop[0];
        TopPts[i].Poly2  = CylPolyPtsTop[i-1];
      }
      else
      {
        BottomPts[i].Center = CylPolyPtsBottom[i];
        BottomPts[i].Poly1  = CylPolyPtsBottom[i+1];
        BottomPts[i].Poly2  = CylPolyPtsBottom[nbCyl-1-2*i];

        TopPts[i].Center = CylPolyPtsTop[i];
        TopPts[i].Poly1  = CylPolyPtsTop[i+1];
        TopPts[i].Poly2  = CylPolyPtsTop[nbCyl-1-2*i];
      }
    }
  }
  else if ( nbCyl == 4 )
  {
    for ( size_t i=0; i<nbCyl; ++i )
    {
      if ( i == 0 )
      {
        BottomPts[i].Center = CylPolyPtsBottom[i];
        BottomPts[i].Poly1  = CylPolyPtsBottom[i+1];
        BottomPts[i].Poly2  = CylPolyPtsBottom[nbCyl-1];

        TopPts[i].Center = CylPolyPtsTop[i];
        TopPts[i].Poly1  = CylPolyPtsTop[i+1];
        TopPts[i].Poly2  = CylPolyPtsTop[nbCyl-1];
      }
      else if ( i == nbCyl - 1 )
      {
        BottomPts[i].Center = CylPolyPtsBottom[i];
        BottomPts[i].Poly1  = CylPolyPtsBottom[0];
        BottomPts[i].Poly2  = CylPolyPtsBottom[i-1];

        TopPts[i].Center = CylPolyPtsTop[i];
        TopPts[i].Poly1  = CylPolyPtsTop[0];
        TopPts[i].Poly2  = CylPolyPtsTop[i-1];
      }
      else
      {
        BottomPts[i].Center = CylPolyPtsBottom[i];
        BottomPts[i].Poly1  = CylPolyPtsBottom[i+1];
        BottomPts[i].Poly2  = CylPolyPtsBottom[i-1];

        TopPts[i].Center = CylPolyPtsTop[i];
        TopPts[i].Poly1  = CylPolyPtsTop[i+1];
        TopPts[i].Poly2  = CylPolyPtsTop[i-1];
      }
    }
  }

  // Pair of points for each cylinder
  vector<pair<Point,Point> > TopIntersectPts( nbCyl,
    make_pair(pointNul, pointNul) ), BottomIntersectPts = TopIntersectPts;
  vector<double> BottomAngle( nbCyl, 0. ), TopAngle = BottomAngle;

  transform = *polyhedron->getForme()->getTransform();
  for ( size_t i=0; i<nbCyl; ++i )
  {
    intersect2d( BottomPts[i].Center, BottomPts[i].Poly1, BottomCenter[i],
      BottomIntersectPts[i].first, CylRadius[i] );

    intersect2d( TopPts[i].Center, TopPts[i].Poly1, TopCenter[i],
      TopIntersectPts[i].first, CylRadius[i] );

    intersect2d( BottomPts[i].Center, BottomPts[i].Poly2, BottomCenter[i],
      BottomIntersectPts[i].second, CylRadius[i] );

    intersect2d( TopPts[i].Center, TopPts[i].Poly2, TopCenter[i],
      TopIntersectPts[i].second, CylRadius[i] );

    // Calcul de l'angle forme par les interesctions et le centre des bases
    // des cylindres
    BottomAngle[i] = acos( ( BottomIntersectPts[i].first - BottomCenter[i] ) *
	( BottomIntersectPts[i].second - BottomCenter[i] ) /
	( Norm( BottomIntersectPts[i].first - BottomCenter[i] ) *
	Norm( BottomIntersectPts[i].second - BottomCenter[i] ) ) ) ;

    TopAngle[i] = acos( TopIntersectPts[i].first *
	TopIntersectPts[i].second / ( Norm( TopIntersectPts[i].first ) *
	Norm( TopIntersectPts[i].second ) ) ) ;

    BottomIntersectPts[i].first = transform( BottomIntersectPts[i].first );
    BottomIntersectPts[i].second = transform( BottomIntersectPts[i].second );
    TopIntersectPts[i].first = transform( TopIntersectPts[i].first );
    TopIntersectPts[i].second = transform( TopIntersectPts[i].second );
  }


  for ( size_t i=0; i<nbElt; ++i )
    if ( m_elementaryParticules[i]->getConvexType() == CYLINDER )
    {
      transform = *m_elementaryParticules[i]->getForme()->getTransform() ;
      for ( size_t j=0; j<3; ++j )
	/* Une particule (X,Y,Z) apres une autre */
	allPointsCyl.push_back( transform(m_elementaryParticules[i]
	  ->getEnveloppe()[j]
	   ));
    }
    else
    {
      transform = *m_elementaryParticules[i]->getForme()->getTransform() ;
      for ( size_t j=0; j<(m_elementaryParticules[i]->getEnveloppe()).size();
	  ++j )
	  allPointsPoly.push_back( transform( m_elementaryParticules[i]
	    ->getEnveloppe()[j]
	  ) ) ;
    }

  setRayon();

  vector<Point>::iterator point;

  fileOut << getRayon() << " " << allPointsCyl.size() << " "
	  << allPointsPoly.size() << "\n";
  fileOut << nbElt << '\n';

  // Envoie des point d'intersection entre les cylindres et le prisme
  for ( size_t i=0; i<nbCyl; ++i )
  {
    fileOut << TopIntersectPts[i].first << '\t'
	    << TopIntersectPts[i].second << '\t'
	    << BottomIntersectPts[i].first << '\t'
	    << BottomIntersectPts[i].second << '\t'
	    << TopAngle[i] << '\t'
	    << BottomAngle[i] << '\n';
  }

  // Cas 2D
  if ( Grains_BuilderFactory::getContext() == DIM_2 )
  {
    cout << "WARNING!!! Class CompParticule in 2 DIM "
         << "is not implemented!\n"
         << "Need for an assistance! Stop running!" << endl;
    exit(10);
  }
  // Cas 3D
  else if ( Grains_BuilderFactory::getContext() == DIM_3 )
  {
    // Traitement des cylindres d'abord
    // Envoie des points enveloppes
    for (point=allPointsCyl.begin(); point!=allPointsCyl.end(); point++)
    {
      pointEnvelop = (*point);// + *m_geoFormeVdw->getCentre();
      fileOut << pointEnvelop[X] << " " << pointEnvelop[Y] << " "
	<< pointEnvelop[Z] << "\n";
    }
    // Traitement du polyhedre
    for (point=allPointsPoly.begin(); point!=allPointsPoly.end(); point++)
    {
      pointEnvelop = (*point);// + *m_geoFormeVdw->getCentre();
      fileOut << pointEnvelop[X] << " " << pointEnvelop[Y] << " "
	<< pointEnvelop[Z] << "\n";
    }

    // !!! Faces du polyhedre apres le traitement des cylindres !!!
    for ( size_t i=0; i<nbElt; ++i )
    {
      vector< vector<int> > const* allFaces  = m_elementaryParticules[i]
	->getForme()->getConvex()->getFaces();
      vector< vector<int> >::const_iterator face;
      if ( allFaces ) // i.e. Polyhedron
      {
        fileOut << allFaces->size() << '\n';
        for (face=allFaces->begin(); face!=allFaces->end(); face++)
        {
          vector<int>::const_iterator index;
          fileOut << (*face).size() << " ";
          for (index=(*face).begin(); index!=(*face).end(); index++)
            fileOut << (*index) << " ";
          fileOut << '\n';
        }
      }
    }
    for ( size_t i=0; i<nbElt-1; ++i )
      fileOut << "0" << '\n';

  }
  else
  {
    cout << "!!! Warning: Physical dimension undefined (DIM_2 or DIM_3)"
    	<< endl;
    exit(0);
  }
}




/* Return the 2D intersection point of a circle and a segment
-------------------------------------------------------------------*/
void CompParticule::intersect2d( Point &ptA, Point &ptB, Point &ptC,
    Point &ptI, double const &radius )
{
  // !!! The particle is oriented in Y-direction !!!
  vector<double> ptInter(2,0.);
  // Equation of a the segment made by ptA and ptB : y = a*x + b
  // Finding the coefficients a and b
  double aa, bb;
  double x0 = ptC[X], z0 = ptC[Z]; // Center of the circle
  double x1 = ptA[X], z1 = ptA[Z]; // Point inside the circle
  double x2 = ptB[X], z2 = ptB[Z];
  double AA, BB, CC;
  double xx1, xx2, zz1, zz2;
  double dot1, dot2, det;

  double eps = 1.e-9;

  if ( fabs( x1 - x2 ) != 0. )
  {
    if ( fabs(x1) <= eps && eps <= fabs(x2) )
    {
      bb = z1;
      aa = ( z2 - bb ) / x2;
    }
    else if ( fabs(x1) <= eps && fabs(x2) <= eps )
    {
      bb = z1;
      aa = 0.;
    }
    else
    {
      bb = ( z2*x1 - z1*x2 ) / ( x1 - x2 );
      aa = ( z1 - bb ) / x1;
    }
    // Equation of a circle : (x-x0)^2 + (z-z0)^2 = r^2
    // Substituing z = ax + b into (x-x0)^2 + (z-z0)^2 = r^2 gives
    // Ax^2 + Bx + C = 0
    AA = aa*aa + 1.;
    BB = 2.*( aa*bb - aa*z0 - x0 );
    CC = z0*z0 - radius*radius + x0*x0 - 2.*bb*z0 + bb*bb;
    det = BB*BB - 4.*AA*CC;

    // det < 0 and det = 0 are not checked because the line intersects
    // the circle
    // The resulting coordinates (2 solutions)
    xx1 = ( 1./(2.*AA) ) * ( -BB + sqrt( det ) );
    xx2 = ( 1./(2.*AA) ) * ( -BB - sqrt( det ) );
    zz1 = aa*xx1 + bb;
    zz2 = aa*xx2 + bb;
  }
  else
  {
    xx1 = xx2 = x1;
    AA = 1.;
    BB = -2.*z0;
    CC = z0*z0 - radius*radius + (x1-x0)*(x1-x0);
    det = BB*BB - 4.*AA*CC;
    zz1 = ( 1./(2.*AA) ) * ( -BB + sqrt( det ) );
    zz2 = ( 1./(2.*AA) ) * ( -BB - sqrt( det ) );
  }

  // Compute the dot production to check whether the point is in the correct
  // direction or not. If dot > 0, the point lies on the same direction
  // as the vector AB

  dot1 = (xx1-x1) * (x2-x1) + (zz1-z1) * (z2-z1);
  dot2 = (xx2-x1) * (x2-x1) + (zz2-z1) * (z2-z1);

  if ( dot1 > 0. && dot2 < dot1 )
  {
    ptInter[0] = xx1;
    ptInter[1] = zz1;
  }
  else
  {
    ptInter[0] = xx2;
    ptInter[1] = zz2;
  }

  ptI[X] = ptInter[0];
  ptI[Y] = ptC[Y];
  ptI[Z] = ptInter[1];
}
