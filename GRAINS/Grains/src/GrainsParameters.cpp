#include "Voisins.hh"
#include "GrainsParameters.H"
#include "Contact_BuilderFactory.hh"
#include "Obstacle_BuilderFactory.H"
#include "ObstacleChargement.H"
#include "Grains_BuilderFactory.H"
#include "PostProcessingWriter.hh"
#include "Paraview_PostProcessingWriter.hh"
#include "ContactLaw.hh"


//-----------------------------------------------------------------------------
// Constructeur par defaut
GrainsParameters::GrainsParameters() :
    Grains()
  , m_vRelative( 0. ) 
{}




//-----------------------------------------------------------------------------
// Destructeur
GrainsParameters::~GrainsParameters()
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction du probleme.
void GrainsParameters::Construction( DOMElement* rootElement )
{
  assert(rootElement != NULL);
  DOMNode* root = ReaderXML::getNode(rootElement, "Construction");

  unsigned short optio_ = 0x00;
  string restart;

  // Recipient
  DOMNode* recipient = ReaderXML::getNode( root, "Recipient" );
  double lx = ReaderXML::getNodeAttr_Double( recipient, "LX" );
  double ly = ReaderXML::getNodeAttr_Double( recipient, "LY" );
  double lz = ReaderXML::getNodeAttr_Double( recipient, "LZ" );
  DOMNode* origine_recipient = ReaderXML::getNode( root, "Origine" );
  double ox = 0., oy = 0., oz = 0. ;
  if ( origine_recipient )
  {
    ox = ReaderXML::getNodeAttr_Double( origine_recipient, "OX" );
    oy = ReaderXML::getNodeAttr_Double( origine_recipient, "OY" );
    oz = ReaderXML::getNodeAttr_Double( origine_recipient, "OZ" ); 
  }   
  App::setD( lx, ly, lz, ox, oy, oz );
 
  // Decomposition de domaine
  readDomainDecomposition( root, lx - ox, ly - oy, lz - oz ); 
  
  if ( m_processorIsActiv )
  {
    // Particules 
    DOMNode* particules = ReaderXML::getNode( root, "Particules" );
    if ( particules ) 
    {
      optio_ = optio_ | 0x02;
      
      int nbPC = int(m_composants.getParticuleClassesReference()->size());
      DOMNodeList* allParticules = ReaderXML::getNodes( rootElement, 
      	"Particule" );
      for (XMLSize_t i=0; i<allParticules->getLength(); i++) 
      {
        DOMNode* nParticule = allParticules->item( i );
        int nbre = ReaderXML::getNodeAttr_Int( nParticule, "Nombre" );

        // Remarque: les particules de référence ont un numéro générique -1
        // d'où "false" dans le constructeur pour auto_numbering = false
        Particule* particuleRef = new Particule( nParticule, false,
		nbPC + int(i) );
        m_composants.AjouterClasseParticules( particuleRef );
	pair<Particule*,int> ppp( particuleRef, nbre );
	m_newParticules.push_back( ppp );     
      }
    }

    // Obstacles 
    DOMNode* obstacles = ReaderXML::getNode( root, "Obstacles" );
    if ( obstacles ) 
    {
      optio_ = optio_ | 0x04;

      DOMNodeList* allCompObstacles = ReaderXML::getNodes( obstacles );
      for (XMLSize_t i=0; i<allCompObstacles->getLength(); i++) 
      {
        DOMNode* nCompObs = allCompObstacles->item( i );
	Obstacle *obstacle = Obstacle_BuilderFactory::create( nCompObs );
	m_composants.Ajouter( obstacle );   
      }	
//       DOMNode* allObstacles = ReaderXML::getNodeNext( obstacles );
//       Obstacle *obstacle = Obstacle_BuilderFactory::create( allObstacles );
//       m_composants.Ajouter( obstacle );
    }

    // Contact
    DOMNode* contact = ReaderXML::getNode( root, "Contacts" );
    if ( contact )
      Contact_BuilderFactory::define( contact );
    string check_matA, check_matB;
    bool contactLaws_ok = Contact_BuilderFactory::checkContactLawsExist(
    	check_matA, check_matB );
    if ( !contactLaws_ok )
    {
      if ( m_rank == 0 )
        cout << "Pas de loi de contact disponible pour les materiaux : "
	   << check_matA << " & " << check_matB << endl;
       grainsAbort();
    }      
  }
  else optio_ = 0x01 ;

  // Validation & affichage informatif
  if ( optio_ == 0x00 ) {
    if ( m_rank == 0 )
      cout << "ERR : Mise en donnees incomplete sur <Contruction>\n";
    grainsAbort();
  }  
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction des forces actives
void GrainsParameters::Forces( DOMElement* rootElement )
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction de la simulation.
void GrainsParameters::Chargement( DOMElement* rootElement )
{
  if ( m_processorIsActiv )
  {
    assert(rootElement != NULL);
    DOMNode* root = ReaderXML::getNode( rootElement, "Simulation" );

    // Pas de temps
    DOMNode* increment = ReaderXML::getNode( root, "IncreTemps" );
    m_dt = ReaderXML::getNodeAttr_Double( increment, "dt" );

    // Vitesse relative
    DOMNode* nodeVrel = ReaderXML::getNode( root, "VitesseRelative" );
    if ( !nodeVrel )
    {
      cout << "Donner la vitesse relative de contact dans <Simulation> "
      	<< "telle que:" << endl;
      cout << "<VitesseRelative value=\"xxx\"/>" << endl;
      grainsAbort();
    }
    m_vRelative = ReaderXML::getNodeAttr_Double( nodeVrel, "value" );
  }
}




// ----------------------------------------------------------------------------
// Gestion de la simulation
void GrainsParameters::Simulation( bool predict, bool isPredictorCorrector,
	bool explicit_added_mass )
{
  if ( m_processorIsActiv )
  {
    vector<Particule*> const* particuleClasses = 
    	m_composants.getParticuleClassesReference();
    size_t nClasses = particuleClasses->size(), i, j;

    cout << "CONTACTS PARTICULE/PARTICULE" << endl;    
    // Entre particules de la meme classe
    for (i=0;i<nClasses;++i)
    {
      Particule* particule = new Particule( *((*particuleClasses)[i]) );    
      cout << "Bilan du contact Classe " << i << " / Classe " << i << endl;
      Contact_BuilderFactory::contactForceModel( 
       	(*particuleClasses)[i]->materiau(),
      	particule->materiau() )->computeAndWriteEstimates( 
		(*particuleClasses)[i], particule, m_vRelative, cout );    
      delete particule;
    }    

    // Entre particules de classe differente
    for (i=0;i<nClasses;++i)
      for (j=i+1;j<nClasses;++j)
      {
        cout << "Bilan du contact Classe " << i << " / Classe " << j << endl;	
	Contact_BuilderFactory::contactForceModel(
		(*particuleClasses)[i]->materiau(),
      		(*particuleClasses)[j]->materiau() )->computeAndWriteEstimates(
		(*particuleClasses)[i], (*particuleClasses)[j], m_vRelative, 
			cout ); 	
      }

    cout << "CONTACTS PARTICULE/OBSTACLE" << endl;    
    // Particules-obstacles
    list<MonObstacle*> allObs = m_composants.getObstacles()->getObstacles();
    list<MonObstacle*>::iterator il;    
    for (i=0;i<nClasses;++i)
      for (il=allObs.begin();il!=allObs.end();il++)
        if ( (*il)->materiau() != "periode" )
        {
          cout << "Bilan du contact Classe " << i << " / " << (*il)->getName() 
		<< endl;
	  Contact_BuilderFactory::contactForceModel(
	  	(*particuleClasses)[i]->materiau(), (*il)->materiau() )
		->computeAndWriteEstimates( (*particuleClasses)[i], *il,
		m_vRelative, cout ); 
	} 
	
    cout << "DEPLACEMENT PARTICULES" << endl;    
    // Entre particules de la meme classe
    for (i=0;i<nClasses;++i)
    {
      cout << "  Classe " << i << " : v0*dt/croute = " << 
      	m_vRelative * m_dt / (*particuleClasses)[i]->getRayonInteraction() 
	<< endl;
    }   	     
  }
}
