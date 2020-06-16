#include "Grains_Exec.hh"
#include "SystemState.hh"
#include "EnsComposant.H"
#include "Particule.H"
#include "App.H"
#include "LinkedCell.H"
#include "ObstaclePeriodique.hh"
#include "ParticulePeriodique.hh"


/* Constructeur par defaut 
--------------------------*/
SystemState::SystemState( Scalar const& temps_sauvegarde,
	EnsComposant const& ENScomposants )
{
  list<Particule*> const* particules = ENScomposants.getParticulesActives();
  list<Particule*>::const_iterator particule;
  list<Particule*>::const_iterator clone;
  struct EtatParticule* etat = NULL;
  struct EtatCloneMonoPeriodique* etatClone = NULL;
  pair<ConfigurationMemento*,CineParticuleMemento*> ppp(NULL,NULL);

  // Temps de sauvegarde
  m_temps_sauvegarde = temps_sauvegarde ;
  
  // Etat des particules
  for (particule=particules->begin(); particule!=particules->end(); particule++)
  {
    etat = new EtatParticule;
    etat->PartID = (*particule)->getID();
    ppp = (*particule)->createState();
    etat->memento_config = ppp.first;
    etat->memento_cine = ppp.second;
    etat->PartClasse = (*particule)->getParticuleClasse();
    etat->EtatPeriodicClones = NULL ;
    if ( (*particule)->getNombreClonesPeriodiques() )
    {
      list<Particule*> const* periodicClones = 
      		(*particule)->getPeriodicClones();
      etat->EtatPeriodicClones = new list<struct EtatCloneMonoPeriodique*>;
      for (clone=periodicClones->begin();clone!=periodicClones->end();clone++)
        if ( (*clone)->getNbPeriodes() == 1 )	
        {
          etatClone = new EtatCloneMonoPeriodique;
	  etatClone->memento_config = (*clone)->createConfigState(); 
	  etatClone->pPerClone = *clone ;
          etatClone->contactReferenceObstacle = (*clone)->getObstacle();
	  etat->EtatPeriodicClones->push_back(etatClone);
        }
    }		
    allEtats.push_back(etat);
  }
  
  // Etat des obstacles
  ENScomposants.createStateObstacles( obsStates ); 
}




/* Destructeur 
--------------*/
SystemState::~SystemState()
{
  list<struct EtatParticule*>::iterator il;
  for (il=allEtats.begin();il!=allEtats.end();il++)
  {
    emptyEtatParticule( *il );
    delete *il;
  }
  allEtats.clear();
}




/* Retablit l'etat du systeme tel qu'il a ete sauvegarde
--------------------------------------------------------*/    
void SystemState::RestaureState( Scalar& temps_sauvegarde, 
	EnsComposant& ENScomposants,
	LinkedCell* LC,
	const int &rank)
{
  list<Particule*>* particules = ENScomposants.getParticulesActives();
  list<Particule*>* ParticulesHalozone = ENScomposants.getParticulesHalozone();
  list<Particule*>* particulesClones = ENScomposants.getParticulesClones();
  list<Particule*>* particulesClonesPeriodiques = 
  	ENScomposants.getParticulesClonesPeriodiques(); 
  list<Particule*>::iterator particule,periodicClone;
  list<struct EtatParticule*>::iterator et;
  list<struct EtatCloneMonoPeriodique*>::iterator etPer;
  list<Particule*> toBedestroyed;  
  bool found = false;
  int ID = 0, tag_old = 0, tag_new = 0, nbComp = 0;
  EtatParticule* etat = NULL;
  Particule *pdestroy = NULL, *pRef = NULL, *newPart = NULL, *perclone = NULL;
  vector<Particule*> const* ParticuleClassesReference = 
  	ENScomposants.getParticuleClassesReference();

  // Restauration du temps
  temps_sauvegarde = m_temps_sauvegarde;


  // Initialise l'etat de restauration des clones mono-periodiques
  if ( Grains_Exec::m_periodique == true )  
    for (particule=particulesClonesPeriodiques->begin();
  	particule!=particulesClonesPeriodiques->end(); particule++)
      (*particule)->setRestaureState( false );


  // Mise a jour des particules presentes
  for (particule=particules->begin(); particule!=particules->end(); particule++)
  {
    ID = (*particule)->getID();
      
    // Recherche de la particule dans l'etat sauvegarde
    // Si la particule se trouve dans l'etat sauvegarde => mise a jour    
    found = false; 
    for (et=allEtats.begin();et!=allEtats.end() && !found; )
      if ( (*et)->PartID == ID ) 
      {
        found = true;
	etat = *et;

	// Tag a l'etat present
	tag_new = (*particule)->getTag();
	
	// Suppresion de la cellule actuelle
	LC->remove( *particule );

	// Affectation de la cinematique & configuration de l'etat sauvegarde
	(*particule)->restaureState( etat->memento_config, etat->memento_cine );
	
        // Link avec la cellule
	LC->Link( *particule );
		
	// Tag a l'etat sauvegarde
	tag_old = (*particule)->getTag();		
	
	// Gestion dans les listes de ENScomposants si changement de tag
	if ( tag_new != tag_old )
	{          
	  // Tag present vaut 0	  
	  if ( tag_new == 0 )
	  {
	    if ( tag_old == 1 )
	      ParticulesHalozone->push_back(*particule);
	    else if ( tag_old == 2 )
	      particulesClones->push_back(*particule);
	  }
	  
          // Tag present vaut 1	  
	  if ( tag_new == 1 )
	  {
	    removeParticuleFromList( *ParticulesHalozone, *particule );
	    if ( tag_old == 2 ) particulesClones->push_back(*particule);
	  }
	  
          // Tag present vaut 2	  
	  if ( tag_new == 2 )
	  {
	    removeParticuleFromList( *particulesClones,*particule );
	    if ( tag_old == 1 ) ParticulesHalozone->push_back(*particule);
	  }	  	     	
	}
		
	// Clones mono-periodiques
	if ( etat->EtatPeriodicClones )
	{
	  for (etPer=etat->EtatPeriodicClones->begin();
	  	etPer!=etat->EtatPeriodicClones->end();etPer++)
	  {
            if ( (*particule)->periodicCloneExistence( (*etPer)->pPerClone ) )
	    {
              perclone = (*etPer)->pPerClone ;
	      
	      // Mise à jour de l'etat	    
	      perclone->restaureState( (*etPer)->memento_config,
	      	etat->memento_cine, (*etPer)->contactReferenceObstacle,
		*particule );
	      perclone->setRestaureState( true );
		
              // LinkUpdate avec le LinkedCell
              LC->LinkUpdateActiveParticule( perclone );
	    }
	    else
	    {
              // Creation 
	      perclone = new ParticulePeriodique( *particule,
	    	(*ParticuleClassesReference)[
			(*particule)->getParticuleClasse()],
  		(*etPer)->contactReferenceObstacle,
		*(*etPer)->contactReferenceObstacle->getPeriode(), 1 );	

              // Link avec le LinkedCell
              LC->Link( perclone );
	
              // Ajout dans les differentes listes de EnsComposants
              particulesClonesPeriodiques->push_back(perclone);
	  
	      // Ajout dans la liste des clones de la particule de reference
	      (*particule)->addPeriodicClone( perclone );		    
	    }
	  }	  
	}
	emptyEtatParticule( etat );
	delete etat;
	et = allEtats.erase(et);
      }
      else et++;
      
    // Si la particule ne se trouve pas dans l'etat sauvegarde
    if (!found) toBedestroyed.push_back(*particule);        
  }   

  
  // Destruction des particules n'appartenant pas a l'etat sauvegarde
  for (particule=toBedestroyed.begin();particule!=toBedestroyed.end();
  	particule++)
  {
    pdestroy = *particule;
    
    // Suppression de la particule du LinkedCell
    LC->remove( pdestroy );	
 
    // Suppression des differentes listes
    removeParticuleFromList( *particules, pdestroy );
    tag_new = pdestroy->getTag();
    if ( tag_new == 1 )
      removeParticuleFromList( *ParticulesHalozone, pdestroy );
    else if ( tag_new == 2 )
      removeParticuleFromList( *particulesClones, pdestroy );

    // Destruction de l'objet pointé
    delete pdestroy;    
  }

  
  // Creation des particules appartenant a l'etat sauvegarde
  // et n'existant plus sur ce processeur
  for (et=allEtats.begin();et!=allEtats.end();et++)
  {
    etat = *et;      
    nbComp = Composant::getNbComposantsCrees();

    // Classe de la particule a creer
    pRef = (*ParticuleClassesReference)[etat->PartClasse];
 
    // Creation
    // Le constructeur par copie de Particule incremente par defaut 
    // le nombre de composants crees, or ici on ne rajoute pas de composant
    // donc on le garde a la bonne valeur
    newPart = new Particule(*pRef);
    Composant::setNbComposantsCrees( nbComp );
   
    // Affectation de l'etat
    newPart->setID(etat->PartID);
    newPart->restaureState(etat->memento_config,etat->memento_cine); 
    newPart->setActivity(COMPUTE);   

    // Ajout de la particule dans le  LinkedCell
    LC->Link(newPart);
	
    // Tag a l'etat sauvegarde
    tag_old = newPart->getTag();

    // Ajout dans les differentes listes de EnsComposants
    particules->push_back(newPart); 
    if ( tag_old == 1 ) ParticulesHalozone->push_back(newPart);
    else if ( tag_old == 2 ) particulesClones->push_back(newPart);     
  }


  // Destruction des clones mono-periodiques qui n'existaient pas a l'etat
  // precedent
  if ( Grains_Exec::m_periodique == true )  
    for (particule=particulesClonesPeriodiques->begin();
  	particule!=particulesClonesPeriodiques->end(); )
      if ( (*particule)->getNbPeriodes() == 1 && 
    	!(*particule)->getRestaureState() )      
      {
        pdestroy = *particule ;
	
        // Supprime le clone dans la liste des clones de la particule de ref
        (*particule)->getPeriodicReference()->erasePeriodicClone( pdestroy );

        // Suppression de la particule dans le LinkedCell
	// Avant on fait un link update si d'aventures la particule a changé de
	// cellule sur le dernier pas de temps, car le link update des clones
	// periodiques est effectue en debut de pas de temps et non en fin
        LC->LinkUpdateActiveParticule( pdestroy ); 
	LC->remove( pdestroy );  	       
	
        particule = particulesClonesPeriodiques->erase(particule);
	
        // Destruction de l'objet pointé
        delete pdestroy;   
      }
      else particule++;
    
  
  // Restauration des obstacles
  ENScomposants.restaureStateObstacles( obsStates );
}  	 




/* Detruit un etat
------------------*/
void SystemState::emptyEtatParticule( struct EtatParticule* etat )
{
  delete etat->memento_config;
  delete etat->memento_cine;
  if ( etat->EtatPeriodicClones )
  {
    for (list<struct EtatCloneMonoPeriodique*>::iterator 
  	il=etat->EtatPeriodicClones->begin();
     	il!=etat->EtatPeriodicClones->end();il++)
    {
      delete (*il)->memento_config;
      delete *il;    
    }
    etat->EtatPeriodicClones->clear();
    delete etat->EtatPeriodicClones; 
  }
}
