#include "GrainsCoupledWithFluidMPI.hh"
#include "InterfaceFluide_BuilderFactory.hh"
#include "Contact_BuilderFactory.hh"
#include "LinkedCell.H"
#include "AddedMass.H"
#include "Obstacle_BuilderFactory.H"
#include "ObstacleChargement.H"
#include "Grains_Exec.hh"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur par defaut
GrainsCoupledWithFluidMPI::GrainsCoupledWithFluidMPI( Scalar rhoFluide,  
	double grid_size ) :   
  GrainsCoupledWithFluid( rhoFluide, grid_size ),
  GrainsMPI()
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructeur
GrainsCoupledWithFluidMPI::~GrainsCoupledWithFluidMPI()
{
  if ( !m_bd_completed ) BeforeDestructor();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction du probleme.
void GrainsCoupledWithFluidMPI::Construction( DOMElement* rootElement )
{
  assert(rootElement != NULL);
  DOMNode* root = ReaderXML::getNode( rootElement, "Construction" );

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

  // Schema d'integration en temps
  DOMNode* timeIntegrator = ReaderXML::getNode( root, "TimeIntegration" );
  if ( timeIntegrator )
    Grains_Exec::m_TIScheme = ReaderXML::getNodeAttr_String( timeIntegrator, 
    	"Type" );

  if ( m_processorIsActiv )
  {
    // Definition des "applications" accessibles
    // L'algorithme sec est present par defaut.
    m_sec = new LinkedCell();
    m_sec->setName( "LinkedCell" ); 
    m_allApp.push_back( m_sec );    

    // Reload ?
    DOMNode* reload = ReaderXML::getNode( root, "Reload" );
    if ( reload ) 
    {
      // Mode de restart
      string reload_type;
      if ( m_forceReloadSame ) reload_type = "same" ;
      else reload_type = ReaderXML::getNodeAttr_String( reload, "Type" );
      
      if ( reload_type == "new" || reload_type == "same" )
        Grains_Exec::m_ReloadType = reload_type ;

      // Lecture du nom de fichier de simulation precedent pour reload
      // Si le mode est "same", le fichier de reload est le même que
      // le fichier de sortie
      if ( Grains_Exec::m_ReloadType == "new" )
        restart  = ReaderXML::getNodeAttr_String( reload, "Fichier" );
      else
      {
        DOMNode* rootSimu = ReaderXML::getNode( rootElement, "Simulation" ); 
	DOMNode* fileRestartOutput = ReaderXML::getNode( rootSimu, "Fichier" );
        restart = ReaderXML::getNodeValue_String( fileRestartOutput );
	restart = Grains_Exec::restartFileName_AorB( restart, "_RFTable.txt" );
	Grains_Exec::m_reloadFile_suffix = 
		restart.substr( restart.size()-1, 1 );
      }	
      restart = fullResultFileName( restart ); 
      cout << "  Restart du fichier " << restart << endl;
      
      // Extrait le repertoire de reload a partir du fichier principal de
      // restart
      Grains_Exec::m_ReloadDirectory = Grains_Exec::extractRoot( restart );  
            
      string   cle;
      ifstream simulLoad(restart.c_str());
      simulLoad >> cle >> m_temps;
      m_new_reload_format = false ;
      if ( cle == "NEW_RELOAD_FORMAT" ) 
      {
        m_new_reload_format = true ;
	simulLoad >> cle >> m_temps;
      }
      Contact_BuilderFactory::reload( simulLoad );
      m_composants.read( simulLoad, restart, m_new_reload_format );
      Contact_BuilderFactory::set_materialsForObstaclesOnly_reload(
      		m_composants.getParticuleClassesReference() );
      simulLoad >> cle;      
      
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
      
      string reset = ReaderXML::getNodeAttr_String( reload, "Vitesse" );
      m_composants.ResetCinematique( reset );
    }

    // Construction du probleme et affectation des composants
    m_rayon = m_composants.getRayonMax();
    
    DOMNode* rootForces = ReaderXML::getNode( rootElement, "Forces" );
    DOMNode* nLubrication = ReaderXML::getNode( rootForces, "LubricationForce");
    // In case of Lubrication correctio, the linkcell size is increased by half
    // radius (typical lubrication force range) Please remind that in case of
    // micro-scale simulation(GrainsCoupledWithFluid) this value is smaller since 
    // lubrication force is partially resolved.
    // We assume that Grains MPI is always used for meso-scale and 
    // Grains Serial is always used for micro-scale model. 

    if ( nLubrication ) 
    {
      m_rayon = m_rayon + 0.25 * m_rayon; 
      if ( m_rank == 0 ) cout << "   LinkedCell size increased 25 percent"<<
      " due to the lubrication correction " << endl;
      Grains_Exec::m_withlubrication = true;
      double eps_cut;
      if ( ReaderXML::hasNodeAttr_Double( nLubrication,"eps_cut" ) )
        eps_cut = ReaderXML::getNodeAttr_Double( nLubrication,"eps_cut" );  
      else 
        eps_cut = 0.02;
      GrainsCoupledWithFluid::m_lubricationForce =
       new AppFluide_LubricationCorrection( m_gridsize, false, eps_cut );
      m_allApp.push_back( GrainsCoupledWithFluid::m_lubricationForce );      
    }
    else GrainsCoupledWithFluid::m_lubricationForce = NULL;    

    defineLinkedCell( m_rayon );
    
    m_sec->Link( m_composants.getObstacles() );
    
    if ( m_rank == 0 ) 
    {
      cout << "Traitement des contacts particule-obstacle "
    	<< "dans le LinkedCell" << endl;       
      cout << endl << "Schema d'integration en temps = " << 
      	Grains_Exec::m_TIScheme << endl << endl;
    }   
  
    // Nb de particules sur tous les proc
    size_t ninsertedPart = m_wrapper->sum_UNSIGNED_INT(
    	m_composants.nbreParticulesActivesOnProc() );
    m_composants.setNbreParticulesOnAllProc( ninsertedPart ); 
    
    // Postprocess contact energy dissipation rate
    DOMNode* rootSimu = ReaderXML::getNode( rootElement, "Simulation" );
    DOMNode* nContDiss = ReaderXML::getNode( rootSimu,
    "ContactDissipationRate" );
    if (nContDiss) Grains_Exec::m_ContactDissipation = true; 
    
    
  } 
  
  // Nb of particles per class
  int nbPC = int(m_composants.getParticuleClassesReference()->size());
  m_composants.prepare_Polydisperse_vectors( nbPC );
  m_composants.set_NbParticlesPerClass();
  m_composants.compute_SauterMeanDiameter();
  m_composants.compute_ParticleClassesConcentration();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction de la simulation.
void GrainsCoupledWithFluidMPI::Chargement( DOMElement* rootElement )
{
  GrainsCoupledWithFluid::Chargement( rootElement );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction des forces actives
void GrainsCoupledWithFluidMPI::Forces( DOMElement* rootElement )
{
  GrainsCoupledWithFluid::Forces( rootElement );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde de l'etat de simulation */
void GrainsCoupledWithFluidMPI::Save( const string &ext ) const
{}




// ----------------------------------------------------------------------------
// Gestion de la simulation
void GrainsCoupledWithFluidMPI::Simulation( bool predict, 
	bool isPredictorCorrector,
	bool explicit_added_mass )
{
  if ( m_processorIsActiv )
  {
    // Initialisation de la cinematique des obstacles
    // Fait au debut du 1er appel par le fluide uniquement, d'ou l'utilisation
    // d'un compteur statique
    static size_t init_counter = 0 ;
    if ( !init_counter ) 
      m_composants.setCinematiqueObstacleSansDeplacement( m_temps, m_dt );
    ++init_counter;
        
    Scalar time = 0.;
  
    // Predicteur/Correcteur
    if ( predict ) 
    {
      Grains::setMode(true);
      if ( isPredictorCorrector ) saveState();
    } 
    else Grains::setMode(false); 

    // Boucle sur un pas de temps fluide
    while( m_simulTime - time > 0.01 * m_dt ) 
    {
      try {
        time  += m_dt;
        m_temps += m_dt;
	
       if( m_composants.IsShrinking() )  
	         m_composants.ShrinkingRate( m_temps );     
	

        // Initialisation de l'indicateur de calcul
        // de la transformation avec epaiseur de croute a faux 
        m_composants.InitializeVdWState( m_temps, m_dt );
        
        // Initialisation des torseurs de force
        m_composants.InitializeForces( m_temps, m_dt, predict );
 
        // Calcul des forces de contact
        m_sec->CalculerForces( m_temps, m_dt,
            m_composants.getParticulesActives() );
        
        // Calcul des forces de masse ajoutée
        if ( predict && m_explicitAddedMass )
          m_explicitAddedMass->CalculerForces( m_temps, m_dt, 
              m_composants.getParticulesActives() );
        
        // Caclul des forces de trainee hydro en DEM-CFD
        if( app_HydroForce )
          app_HydroForce->CalculerForces( m_temps, m_dt, 
              m_composants.getParticulesActives() );    

        // Caclul du flux de chaleur en DEM-CFD
        if ( app_FluidTemperature )
          app_FluidTemperature->CalculerForces( m_temps, m_dt, 
              m_composants.getParticulesActives() );

        // Deplacement et actualisation des composants
        if ( !b_fixed_particles )
        {
          m_composants.Deplacer( m_temps, m_dt );
          m_composants.Actualiser();
        }
        
        // Caclul de la temperature en DEM-CFD
        if ( app_FluidTemperature )
        {
          m_composants.ComputeTemperature( m_temps, m_dt );
        }

        // Mise à jour du tag des particules intérieures
        m_sec->updateInteriorTag( m_temps,
            m_composants.getParticulesActives(),
            m_composants.getParticulesHalozone(),
            m_wrapper );

        // Création ou Mise à jour des clones
        if ( m_MPIstrategie == "AllgatherGlobal" )
          m_wrapper->UpdateOrCreateClones_AllGatherGlobal( m_temps,
              m_composants.getParticulesClones(),
              m_composants.getParticulesActives(),
              m_composants.getParticulesHalozone(),
              m_composants.getParticuleClassesReference(),
              m_sec );
        else if ( m_MPIstrategie == "AllgatherLocal" )
          m_wrapper->UpdateOrCreateClones_AllGatherLocal( m_temps,
              m_composants.getParticulesClones(),
              m_composants.getParticulesActives(),
              m_composants.getParticulesHalozone(),
              m_composants.getParticuleClassesReference(),
              m_sec );
        else if ( m_MPIstrategie == "SRLocal" )
          m_wrapper->UpdateOrCreateClones_SendRecvLocal( m_temps,
              m_composants.getParticulesClones(),
              m_composants.getParticulesActives(),
              m_composants.getParticulesHalozone(),
              m_composants.getParticuleClassesReference(),
              m_sec );
        else if ( m_MPIstrategie == "SRLocalCommOpt" )
          m_wrapper->UpdateOrCreateClones_SendRecvLocal_GeoLoc( m_temps,
              m_composants.getParticulesClones(),
              m_composants.getParticulesActives(),
              m_composants.getParticulesHalozone(),
              m_composants.getParticuleClassesReference(),
              m_sec );

        // Suppresion des clones hors du Linked Cell (car passes sur un autre
        // processeur)
        m_sec->DestroyOutOfDomainClones( m_temps,
            m_composants.getParticulesClones(),
            m_composants.getParticulesReferencesPeriodiques(),
            m_composants.getParticulesActives() );

        // Mise à jour du tag des particules clones et halozone
        m_sec->updateHalozoneCloneTag( m_temps, 
            m_composants.getParticulesHalozone(),
            m_composants.getParticulesClones(),
            m_wrapper );

        // Actualisation des particules & obstacles dans les cellules
        m_sec->LinkUpdate( m_temps, m_dt, 
            m_composants.getParticulesActives() );
        m_composants.updateGeoLocalisationParticulesHalozone();
      }
      catch (ErreurContact &choc) 
      {
        // Fin de simulation sur choc
        cerr << '\n';
        m_composants.PostProcessingErreurComposants( "ErreurContact",
            choc.getComposants() );
        choc.Message(cerr);
        break;
      } 
      catch (ErreurDeplacement &erreur) 
      {
        // Fin de simulation sur deplacement trop grand
        cerr << '\n';
        erreur.Message(cerr);
        m_composants.PostProcessing( m_temps, m_dt, m_sec, m_rank, 
            m_wrapper->nombre_total_procs(), m_wrapper );
        break;
      } 
      catch (ErreurSimulation &erreur) 
      {
        // Fin de simulation sur erreur
        cerr << '\n';
        erreur.Message(cerr);
        break;
      }
    }
   
    // Caracteristiques globales des contacts pendant la simulation
    Scalar overlap_max_ = AppSec::getOverlapMax();
    Scalar overlap_mean_ = AppSec::getOverlapMean();
    Scalar time_overlapMax_ = AppSec::getTimeOverlapMax(); 
    Scalar nbIterGJK_mean_ = AppSec::getNbIterGJKMean();     
    m_wrapper->ContactsFeatures( overlap_max_, overlap_mean_,
    	time_overlapMax_, nbIterGJK_mean_ );
    AppSec::setContactsFeatures( overlap_max_, overlap_mean_,
    	time_overlapMax_, nbIterGJK_mean_ );
    if ( m_rank == 0 ) 
    {
      cout << "Rayon d'interaction minimal = " << 
    	m_composants.getRayonInteractionMin() << endl;
      cout << "Overlap moyen = " << AppSec::getOverlapMean() << endl;
      cout << "Overlap maximal = " << AppSec::getOverlapMax() << endl;
      cout << "Temps de l'overlap maximal = " << 
	AppSec::getTimeOverlapMax() << endl;     
      cout << "Nb d'iterations moyen de GJK = " << 
    	AppSec::getNbIterGJKMean() << endl; 		
    }    

  }

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Insertion d'une particule dans les algorithmes
bool GrainsCoupledWithFluidMPI::insertParticule( const PullMode& mode )
{
  GrainsCoupledWithFluid::insertParticule( mode );
  return false;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Positionne les particules en attente pour la simulation a partir
// d'un fichier de positions 
void GrainsCoupledWithFluidMPI::setPositionParticulesFichier
	( const PullMode& mode )
{
  GrainsCoupledWithFluid::setPositionParticulesFichier( mode );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Positionne les particules en attente pour la simulation a partir
// d'un bloc structure de positions 
void GrainsCoupledWithFluidMPI::setPositionParticulesBloc
	( const PullMode& mode )
{
  GrainsCoupledWithFluid::setPositionParticulesBloc( mode );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture des donnees pour les simulations MPI
void GrainsCoupledWithFluidMPI::readDomainDecomposition( DOMNode* root,
  	const Scalar& lx, const Scalar& ly, const Scalar& lz )
{
  GrainsMPI::readDomainDecomposition( root, lx, ly, lz );
  Grains_Exec::m_MPI_verbose = 0;
}



 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nom complet du fichier de resultat .resul
string GrainsCoupledWithFluidMPI::fullResultFileName( const string &rootname )
	const
{
  return GrainsMPI::fullResultFileName( rootname );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture des donnees pour les simulations periodiques
void GrainsCoupledWithFluidMPI::readPeriodic( DOMElement* rootElement )
{
  GrainsMPI::readPeriodic( rootElement );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definition du LinkedCell
void GrainsCoupledWithFluidMPI::defineLinkedCell( Scalar const& rayon )
{
  GrainsMPI::defineLinkedCell( rayon );		
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Arret du code en cas de probleme
void GrainsCoupledWithFluidMPI::grainsAbort() const
{
  GrainsMPI::grainsAbort();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modification de la vitesse des particules
void GrainsCoupledWithFluidMPI::UpdateParticulesVelocities(
	const bool &b_set_velocity_nm1_and_diff )
{
  cout << "!!! Warning: GrainsCoupledWithFluidMPI is intended to be coupled "
  	<< "with PeliGRIFF, which uses istringstream instead of text file !!!" 
	<< endl;
  cout << "Thus, UpdateParticulesVelocities(const bool " <<
  	"&b_set_velocity_nm1_and_diff has not been implemented" << endl;	
  grainsAbort();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluidMPI::WriteParticulesInFluid( const string &filename )
	const
{
  cout << "!!! Warning: GrainsCoupledWithFluidMPI is intended to be coupled "
  	<< "with PeliGRIFF, which uses istringstream instead of text file !!!" 
	<< endl;
  cout << "Thus, WriteParticulesInFluid(const string &filename) has not been "
  	<< "implemented" << endl;	
  grainsAbort();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modification de la vitesse des particules
void GrainsCoupledWithFluidMPI::UpdateParticulesVelocities( 
	const vector<vector<double> > &velocities,
	const bool &b_set_velocity_nm1_and_diff )
{
  if ( m_processorIsActiv )
    m_InterfaceFluide->UpdateParticulesVelocities( 
    	*m_composants.getParticulesActives(),
  	m_dt, velocities, b_set_velocity_nm1_and_diff, true );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluidMPI::WriteParticulesInFluid( istringstream &is ) 
	const
{
  if ( m_processorIsActiv )
  {
    vector<Particule*>* allparticules = 
    	m_wrapper->GatherParticules_PostProcessing(
  	*m_composants.getParticulesActives(),
	*m_composants.getParticulesWait(),
	*m_composants.getParticuleClassesReference(),
	m_composants.nbreParticulesOnAllProc() );

    if ( m_rank == 0 ) 
      m_InterfaceFluide->WriteParticulesInFluid( allparticules,
	m_composants.getObstaclesToFluid(), is );

    if ( allparticules )
    { 
      // Destruction du vecteur
      vector<Particule*>::iterator iv;
      for (iv=allparticules->begin();iv!=allparticules->end();iv++)
        delete *iv; 
      allparticules->clear();
      delete allparticules;
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluidMPI::WritePVGCInFluid( const string &filename ) 
	const
{
  cout << "!!! Warning: GrainsCoupledWithFluidMPI is intended to be coupled "
  	<< "with PeliGRIFF, which uses istringstream instead of text file !!!" 
	<< endl;
  cout << "Thus, WritePVGCInFluid(const string &filename) has not been "
  	<< "implemented" << endl;	
  grainsAbort();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluidMPI::WritePVGCInFluid( istringstream &is ) const
{
  if ( m_processorIsActiv )
  {
    vector<Particule*>* allparticules = 
    	m_wrapper->GatherParticules_PostProcessing(
  	*m_composants.getParticulesActives(),
	*m_composants.getParticulesWait(),
	*m_composants.getParticuleClassesReference(),
	m_composants.nbreParticulesOnAllProc() );

    if ( m_rank == 0 ) 
      m_InterfaceFluide->WritePVGCInFluid( allparticules, is );

    if (allparticules)
    { 
      // Destruction du vecteur
      vector<Particule*>::iterator iv;
      for (iv=allparticules->begin();iv!=allparticules->end();iv++)
        delete *iv; 
      allparticules->clear();
      delete allparticules;
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde par defaut de l'etat initial pour post-processing
void GrainsCoupledWithFluidMPI::InitialPostProcessing()
{
  // Par defaut l'etat initial est sauvegarde
  if (m_processorIsActiv)
    m_composants.PostProcessing_start( m_temps, m_dt, m_sec, m_fenetres,
    	m_rank, m_wrapper->nombre_total_procs(), m_wrapper );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde pour post-processing et restart
void GrainsCoupledWithFluidMPI::doPostProcessing()
{  
  if (m_processorIsActiv)
  {
    // Post processing
    m_composants.PostProcessing( m_temps, m_dt, m_sec, m_rank, 
  	m_wrapper->nombre_total_procs(), m_wrapper );
      
    // Sauvegarde du fichier de fin pour Reload
    saveReload( m_temps );
  }  
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Operations a effectuer avant appel du destructeur
void GrainsCoupledWithFluidMPI::BeforeDestructor()
{
  if ( m_processorIsActiv )
  {
    // Avant destruction, sauvegarde de l'etat final
    m_composants.PostProcessing_end();

    // Sauvegarde du fichier de fin pour Reload
//    saveReload( m_temps ); 
  }
  m_bd_completed = true; 
} 
