#include "Voisins.hh"
#include "GrainsGeomTest.H"
#include "Contact_BuilderFactory.hh"
#include "Obstacle_BuilderFactory.H"
#include "ObstacleChargement.H"
#include "Grains_BuilderFactory.H"
#include "PostProcessingWriter.hh"
#include "Paraview_PostProcessingWriter.hh"
#include "PointC.H"
#include "Segment.H"
#include "Transform.H"


//-----------------------------------------------------------------------------
// Constructeur par defaut
GrainsGeomTest::GrainsGeomTest() :
    Grains() 
{}




//-----------------------------------------------------------------------------
// Destructeur
GrainsGeomTest::~GrainsGeomTest()
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction du probleme.
void GrainsGeomTest::Construction( DOMElement* rootElement )
{
  if (m_processorIsActiv)
  { 
    int nbPC;

    // Particules 
    DOMNode* particules = ReaderXML::getNode( rootElement, "Particules" );
    if (particules) 
    {      
      nbPC = int(m_composants.getParticuleClassesReference()->size());
      assert( nbPC <= 1 );
      DOMNodeList* allParticules = ReaderXML::getNodes( rootElement, 
      	"Particule" );
      for (XMLSize_t i=0; i<allParticules->getLength(); i++) 
      {
        DOMNode* nParticule = allParticules->item( i );

        // Remarque: les particules de référence ont un numéro générique -1
        // d'où "false" dans le constructeur pour auto_numbering = false
        Particule* particuleRef = new Particule( nParticule, false,
		nbPC + int(i) );
        m_composants.AjouterClasseParticules( particuleRef ); 
      }
    }

    nbPC = int(m_composants.getParticuleClassesReference()->size());
    assert( nbPC <= 2 );
   
    // Particule 1
    Particule* particule1 = new Particule(
      	*(*m_composants.getParticuleClassesReference())[0] );
    particule1->setActivity( COMPUTE );
    m_composants.Ajouter( particule1 );

    // Particule 2
    Particule* particule2 = new Particule(
      	*(*m_composants.getParticuleClassesReference())[nbPC==2?1:0] );
    particule2->setActivity( COMPUTE );
    m_composants.Ajouter( particule2 );
  
    // Particule 1
    Point position1;
    DOMNode* nPos1 = ReaderXML::getNode( rootElement, "Particule1" );
    string option1 = ReaderXML::getNodeAttr_String( nPos1, "Type" );
    assert( option1 == "Angle" || option1 == "Matrice" ); 
    if ( option1 == "Angle" )
    {  
      // Position
      DOMNode* Centre1 = ReaderXML::getNode( nPos1, "Centre" );     
      position1[X] = ReaderXML::getNodeAttr_Double( Centre1, "X" );
      position1[Y] = ReaderXML::getNodeAttr_Double( Centre1, "Y" );
      position1[Z] = ReaderXML::getNodeAttr_Double( Centre1, "Z" ); 
    
      // Orientation
      Point orientation1;
      DOMNode* nOr1 = ReaderXML::getNode( nPos1, "OrientationAngle" ); 
      if ( nOr1 )
      {
        orientation1[X] = ReaderXML::getNodeAttr_Double( nOr1, "X" );
        orientation1[Y] = ReaderXML::getNodeAttr_Double( nOr1, "Y" );
        orientation1[Z] = ReaderXML::getNodeAttr_Double( nOr1, "Z" );
      }
      
      // Affectation position - orientation     
      m_composants.getParticule(0)->setPosition( position1 );    
      Transform trot1;
      trot1.setBasis(getMatrixRotation( orientation1[X],
    	orientation1[Y], orientation1[Z]) );
      m_composants.getParticule(0)->composePosition( trot1 );    
    }
    else
    {
      // Lecture & affectation
      m_composants.getParticule(0)->getForme()->getTransform()
      	->load( nPos1 );
    } 
         
   
    // Particule 2
    // Position
    Point position2;
    DOMNode* nPos2 = ReaderXML::getNode( rootElement, "Particule2" );
    string option2 = ReaderXML::getNodeAttr_String( nPos2, "Type" );
    assert( option2 == "Angle" || option2 == "Matrice" );     
    if ( option2 == "Angle" )
    {  
      // Position              
      DOMNode* Centre2 = ReaderXML::getNode( nPos2, "Centre" );     
      position2[X] = ReaderXML::getNodeAttr_Double( Centre2, "X" );
      position2[Y] = ReaderXML::getNodeAttr_Double( Centre2, "Y" );
      position2[Z] = ReaderXML::getNodeAttr_Double( Centre2, "Z" ); 
    
      // Orientation
      Point orientation2;
      DOMNode* nOr2 = ReaderXML::getNode( nPos2, "OrientationAngle" ); 
      if (nOr2)
      {
        orientation2[X] = ReaderXML::getNodeAttr_Double( nOr2, "X" );
        orientation2[Y] = ReaderXML::getNodeAttr_Double( nOr2, "Y" );
        orientation2[Z] = ReaderXML::getNodeAttr_Double( nOr2, "Z" );
      }
      
      // Affectation position - orientation     
      m_composants.getParticule(1)->setPosition( position2 );      
      Transform trot2;
      trot2.setBasis(getMatrixRotation( orientation2[X],
    	orientation2[Y], orientation2[Z]) );
      m_composants.getParticule(1)->composePosition( trot2 );       
    }         
    else
    {
      // Lecture & affectation
      m_composants.getParticule(1)->getForme()->getTransform()
      	->load( nPos2 );    
    } 
  }                
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction des forces actives
void GrainsGeomTest::Forces( DOMElement* rootElement )
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction de la simulation.
void GrainsGeomTest::Chargement( DOMElement* rootElement )
{
    // Post-processing writers
    DOMNode* nPostProcessors = ReaderXML::getNode( rootElement, 
    	"PostProcessingWriters" );
    if (nPostProcessors) 
    {
      DOMNode* paraview = ReaderXML::getNode( nPostProcessors, "Paraview" );
      if ( paraview )
      {
	PostProcessingWriter* ppwPara = new Paraview_PostProcessingWriter( 
      		paraview, m_rank, m_nprocs );
        m_composants.addPostProcessingWriter( ppwPara );
      }             
    } 
}




// ----------------------------------------------------------------------------
// Gestion de la simulation
void GrainsGeomTest::Simulation( bool predict, bool isPredictorCorrector,
	bool explicit_added_mass )
{
  if (m_processorIsActiv)
  {  
    Particule *P1 = m_composants.getParticule(0);
    Particule *P2 = m_composants.getParticule(1);
    
    FormeVdW* geoFormeVdw1 = P1->getForme();
    FormeVdW* geoFormeVdw2 = P2->getForme();        
    
    const Convex* convexA = geoFormeVdw1->getConvex();
    const Convex* convexB = geoFormeVdw2->getConvex();
    
    Scalar rayonVdWA = geoFormeVdw1->getRayonInterAction();
    Scalar rayonVdWB = geoFormeVdw2->getRayonInterAction();        
  
    // Cas general
    Vecteur gcagcb = *P1->getPosition() - *P2->getPosition();
    Scalar distance = 1e20, distance1, distance2, distance3;

    // Distance entre les deux formes
    Point pointA, pointB, point__;
    Point pa,pap,pb,pbp;  
    int nbIterGJK = 0;
    Transform a2w = *geoFormeVdw1->getTransformVdW();
    Transform b2w = *geoFormeVdw2->getTransformVdW();
    Transform const* configA = geoFormeVdw1->getTransform();
    Transform const* configB = geoFormeVdw2->getTransform();    
    distance = closest_points( *convexA, *convexB, a2w, b2w, 
    	pointA, pointB, nbIterGJK );
     
    cout << "Is contact = " << P1->isContact( P2 ) << endl;
    cout << "Distance sans VdW = " << 
    	closest_points( *convexA, *convexB, *configA, *configB, 
		point__, point__, nbIterGJK )
	<< endl;
    cout << "Distance avant VdW = " << distance << endl;
    if ( distance < EPSILON )
      cout << "Choc de croute" << endl;
    else
    {    
      distance1 = distance - rayonVdWA - rayonVdWB;
      cout << "Distance1 = " << distance1 << endl; 

      // Pts de contact     
      pap = (*configA)( pointA );
      pbp = (*configB)( pointB );
      pa = a2w( pointA );
      pb = b2w( pointB );      
      Vecteur ba = pa - pb;
      ba.normalize();
      Vecteur pbTopbp = pbp - pb;
      Vecteur paTopap = pap - pa;
      double crouteA = fabs( paTopap * ba );
      double crouteB = fabs( pbTopbp * ba );             
      cout << "Vecteur BA = " << ba;
      cout << "Croute A = " << crouteA << endl;      
      cout << "Croute B = " << crouteB << endl;       
      distance2 = distance - crouteA - crouteB;    
      cout << "Distance2 = " << distance2 << endl;
      
      cout << "Point A = " << pa;
      cout << "Point B = " << pb;      
      Transform InvConfigA,InvConfigB;
      InvConfigA.invert( *configA );
      InvConfigB.invert( *configB );
      Point ShellA, ShellB;
      ShellA = convexA->intersectionToShell( InvConfigA(pa), InvConfigA(pb) );
      ShellB = convexB->intersectionToShell( InvConfigB(pb), InvConfigB(pa) );
      ShellA = (*configA)( ShellA );
      ShellB = (*configB)( ShellB );      
      cout << "Contact point in A = " << pa;
      cout << "Intersection A = " << ShellA;
      cout << "Contact point in B = " << pb;      
      cout << "Intersection B = " << ShellB; 
      
      double crouteAp = Norm( ShellA - pa );
      double crouteBp = Norm( ShellB - pb );      
      cout << "Croute A = " << crouteAp << endl;      
      cout << "Croute B = " << crouteBp << endl;
      distance3 = distance - crouteAp - crouteBp;    
      cout << "Distance3 = " << distance3 << endl;
    }
    
    m_composants.PostProcessing_start( 0., 1.e-6, m_sec, m_fenetres );    
  }   
}




// ----------------------------------------------------------------------------
// Renvoi une transformation de type rotation en fonction des angles
Matrix GrainsGeomTest::getMatrixRotation( const Scalar &angleX,
  	const Scalar &angleY,
	const Scalar &angleZ ) const
{
  Matrix rotation;
  
  Matrix rZ( cos( angleZ * PI / 180. ), - sin( angleZ * PI / 180. ), 0.,
  	sin( angleZ * PI / 180. ), cos( angleZ * PI / 180. ), 0.,
	0., 0., 1. );

  if ( m_dimension == 3 )
  {
    Matrix rX( 1., 0., 0.,
    	0., cos( angleX * PI / 180. ), -sin( angleX * PI / 180. ),
	0., sin( angleX * PI / 180. ), cos( angleX * PI / 180. ) );
    Matrix rY( cos( angleY * PI / 180. ), 0., sin( angleY * PI / 180. ),
    	0., 1., 0.,
    	-sin( angleY * PI / 180. ), 0., cos( angleY * PI / 180. ) );
    Matrix tmp = rY * rZ;
    rotation = rX * tmp;
  }
  else rotation = rZ;
  
  return rotation;
} 
