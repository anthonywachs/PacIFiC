#include "MPIWrapperGrains.hh"
#include "App.H"
#include "Contact_BuilderFactory.hh"
#include "ObstaclePeriodique.hh"
#include "ParticulePeriodique.hh"
#include "Grains_Exec.hh"


// ----------------------------------------------------------------------------
// Collecte sur le processeur master de l'ensemble des clones
// periodiques sur les différents processeurs pour post-processing
list<Particule*>* MPIWrapperGrains::GatherClonesPeriodiques_PostProcessing(
  	const list<Particule*> &particulesClonesPeriodiques,
	vector<Particule*> const& ParticuleClassesReference ) const
{
  list<Particule*>* allClonesPeriodiques=NULL;
  list<Particule*>::const_iterator il;
  int i,j;

  // Taille des messages proportionnel au nombre de particules dans 
  // ParticulesReference
  // ---------------------------------------------------------------
  int nb_part = int(particulesClonesPeriodiques.size());
  for (il=particulesClonesPeriodiques.begin();
  	il!=particulesClonesPeriodiques.end();il++)
    if ((*il)->getTag()==2) nb_part--;

  // Taille des messages à passer
  int *recvcounts = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts[i]=0;
  MPI_Gather( &nb_part, 1, MPI_INT, recvcounts, 1, MPI_INT, 
  	m_rank_masterWorld, m_MPI_COMM_activProc ); 

  // Partition et taille du buffer de réception
  int *displs = new int[m_nprocs];
  displs[0] = 0; 
  for (i=1; i<m_nprocs; ++i) displs[i] = displs[i-1]+recvcounts[i-1];  
  int recvsize = displs[m_nprocs-1]+recvcounts[m_nprocs-1];



  // Communication des entiers: 
  // Ordre par particule: 
  // [numero de particule, classe]
  // -----------------------------
  int NB_INT_PART = 2;
  int *numClass = new int[NB_INT_PART*nb_part];
  for (il=particulesClonesPeriodiques.begin(),i=0;
  	il!=particulesClonesPeriodiques.end();il++)
  {
    if ( (*il)->getTag() == 0 || (*il)->getTag() == 1 )
    {    
      numClass[i] = (*il)->getID();
      numClass[i+1] = (*il)->getParticuleClasse();
      i+=NB_INT_PART;
    }
  }            

  // Partition et taille du buffer de réception
  int *recvcounts_INT = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts_INT[i] = NB_INT_PART * recvcounts[i];
  int *displs_INT = new int[m_nprocs];
  displs_INT[0] = 0; 
  for (i=1; i<m_nprocs; ++i) 
    displs_INT[i] = displs_INT[i-1] + recvcounts_INT[i-1];  
  int recvsize_INT = displs_INT[m_nprocs-1] + recvcounts_INT[m_nprocs-1];
  int *recvbuf_INT = new int[recvsize_INT];

  // Communication des entiers: 2 integer par particule
  MPI_Gatherv( numClass, NB_INT_PART * nb_part, MPI_INT, recvbuf_INT, 
  	recvcounts_INT, displs_INT, MPI_INT, 
	m_rank_masterWorld, m_MPI_COMM_activProc );



  // Communication des doubles: cinématique & configuration
  // Ordre par particule: 
  // [position, vitesse translation, quaternion rotation, vitesse rotation]
  // ----------------------------------------------------------------------
  int NB_DOUBLE_PART = 26;  
  double *features = new double[NB_DOUBLE_PART*nb_part];
  for (il=particulesClonesPeriodiques.begin(),i=0;
  	il!=particulesClonesPeriodiques.end();il++)
  {
    if ( (*il)->getTag() == 0 || (*il)->getTag() == 1 )
    {
      (*il)->copyVitesseTranslation( features, i );
      (*il)->copyQuaternionRotation( features, i+3 );    
      (*il)->copyVitesseRotation( features, i+7 );
      (*il)->copyTransform( features, i+10 ); 
      i += NB_DOUBLE_PART;
    }                  
  }  

  // Partition et taille du buffer de réception
  int *recvcounts_DOUBLE = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) 
    recvcounts_DOUBLE[i] = NB_DOUBLE_PART * recvcounts[i];
  int *displs_DOUBLE = new int[m_nprocs];
  displs_DOUBLE[0] = 0; 
  for (i=1; i<m_nprocs; ++i) 
    displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];  
  int recvsize_DOUBLE = displs_DOUBLE[m_nprocs-1] 
  	+ recvcounts_DOUBLE[m_nprocs-1];
  double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];
      
  // Communication de la cinématique & configuration: 26 double par 
  // particule
  MPI_Gatherv( features, NB_DOUBLE_PART * nb_part, MPI_DOUBLE, recvbuf_DOUBLE, 
  	recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE, 
	m_rank_masterWorld, m_MPI_COMM_activProc ); 


  // Creation des Particules pour Post-processing
  // --------------------------------------------
  if ( m_rank == m_rank_masterWorld )
  {
    allClonesPeriodiques = new list<Particule*>;

    // Clones periodiques sur tous les proc
    for (j=0;j<recvsize;++j)
    {
      // Creation de la particule
      Particule *part_post=new Particule( recvbuf_INT[NB_INT_PART*j],
      		ParticuleClassesReference[recvbuf_INT[NB_INT_PART*j+1]],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+1],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+2],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+5],		
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+9],		
      		&recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],		
		COMPUTE,
		0 ); 		    
      allClonesPeriodiques->push_back(part_post);    
    }
  } 

  delete [] numClass;
  delete [] features;    
  delete [] recvcounts;
  delete [] recvcounts_INT;
  delete [] recvcounts_DOUBLE;    
  delete [] displs;
  delete [] displs_INT;
  delete [] displs_DOUBLE;      
  delete [] recvbuf_INT;
  delete [] recvbuf_DOUBLE;

  return allClonesPeriodiques;
} 




// ----------------------------------------------------------------------------
// Création et mise à jour des clones periodiques 
void MPIWrapperGrains::commParticulesRefPer_AllGatherGlobal_UpdateOrCreate(
	Scalar time,
  	set<Particule*>* particulesReferencesPeriodiques,
  	list<Particule*>* particulesClonesPeriodiques,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC )
{
  set<Particule*>::const_iterator is;
  list<Particule*>::const_iterator il;  
  int i,j;
  bool found=false;
  bool AdamsBashforth = Grains_Exec::m_TIScheme == "SecondOrderAdamsBashforth";
  
  if ( m_isInCommPeriodic )
  {
    // Taille des messages proportionnel au nombre de particules dans 
    // particulesHalozone
    // --------------------------------------------------------------
    int nb_prp = int(particulesReferencesPeriodiques->size());
    for (is=particulesReferencesPeriodiques->begin();
    	is!=particulesReferencesPeriodiques->end();is++)
      if ( (*is)->getTag() == 2 ) nb_prp--;

    // Taille des messages à passer
    int *recvcounts = new int[m_nprocs_per];
    for (i=0; i<m_nprocs_per; ++i) recvcounts[i]=0;
    MPI_Allgather( &nb_prp, 1, MPI_INT, recvcounts, 1, MPI_INT,
    	*m_MPI_COMM_Periodic ); 

    // Partition et taille du buffer de réception
    int *displs = new int[m_nprocs_per];
    displs[0] = 0; 
    for (i=1; i<m_nprocs_per; ++i) displs[i] = displs[i-1]+recvcounts[i-1];  
    int recvsize = displs[m_nprocs_per-1]+recvcounts[m_nprocs_per-1];


    // Communication des entiers: 
    // Ordre par particule: 
    // [numero de particule, classe, rang expéditeur]  
    // ----------------------------------------------
    int NB_INT_PART = 4;
    int *numClassRang = new int[NB_INT_PART*nb_prp];
    for (is=particulesReferencesPeriodiques->begin(),i=0;
    	is!=particulesReferencesPeriodiques->end();
  	is++)
      if ( (*is)->getTag() != 2 ) 
      {
        numClassRang[i] = (*is)->getID();
        numClassRang[i+1] = (*is)->getParticuleClasse();
        numClassRang[i+2] = (*is)->getPeriodicObstaclesID()->front();     
        numClassRang[i+3] = m_rank;
	i+=NB_INT_PART;
      }            

    // Communication des doubles: cinématique & configuration 
    // Ordre par particule: 
    // [position, vitesse translation, quaternion rotation, vitesse rotation]
    // ----------------------------------------------------------------------
    int NB_DOUBLE_PART = AdamsBashforth ? 38 : 26;
    double *features = new double[NB_DOUBLE_PART*nb_prp];
    for (is=particulesReferencesPeriodiques->begin(),i=0;
    	is!=particulesReferencesPeriodiques->end();
  	is++)
      if ( (*is)->getTag() != 2 )     
      {
        (*is)->copyVitesseTranslation( features, i );
        (*is)->copyQuaternionRotation( features, i+3 );    
        (*is)->copyVitesseRotation( features, i+7 );     
        (*is)->copyTransform( features, i+10,
        *(LC->getObstaclePeriodique( (*is)->getPeriodicObstaclesID()->front() )
		->getPeriode()) );
        if ( AdamsBashforth ) (*il)->copyCinematiqueNm2( features, i+26 );
	i += NB_DOUBLE_PART;
      }


    // Partition et taille du buffer de réception
    int *recvcounts_INT = new int[m_nprocs_per];
    for (i=0; i<m_nprocs_per; ++i) 
      recvcounts_INT[i] = NB_INT_PART * recvcounts[i];
    int *displs_INT = new int[m_nprocs_per];
    displs_INT[0] = 0; 
    for (i=1; i<m_nprocs_per; ++i) 
      displs_INT[i] = displs_INT[i-1] + recvcounts_INT[i-1];  
    int recvsize_INT = displs_INT[m_nprocs_per-1] 
    	+ recvcounts_INT[m_nprocs_per-1];
    int *recvbuf_INT = new int[recvsize_INT];

    // Communication des entiers: 4 integer par particule
    MPI_Allgatherv( numClassRang, NB_INT_PART * nb_prp, MPI_INT, recvbuf_INT, 
  	recvcounts_INT, displs_INT, MPI_INT, *m_MPI_COMM_Periodic );

  
    // Partition et taille du buffer de réception
    int *recvcounts_DOUBLE = new int[m_nprocs_per];
    for (i=0; i<m_nprocs_per; ++i) 
      recvcounts_DOUBLE[i] = NB_DOUBLE_PART * recvcounts[i];
    int *displs_DOUBLE = new int[m_nprocs_per];
    displs_DOUBLE[0] = 0; 
    for (i=1; i<m_nprocs_per; ++i) 
      displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];  
    int recvsize_DOUBLE = displs_DOUBLE[m_nprocs_per-1] 
    	+ recvcounts_DOUBLE[m_nprocs_per-1];
    double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];
      
    // Communication de la cinématique & configuration: 26 double par 
    // particule
    MPI_Allgatherv( features, NB_DOUBLE_PART * nb_prp, MPI_DOUBLE, 
    	recvbuf_DOUBLE, recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE, 
	*m_MPI_COMM_Periodic ); 


    // Creation ou maj des clones periodiques
    // --------------------------------------
    list<App*>::iterator app;
    Point gc;
    for (j=0;j<recvsize;++j)
    {
      il = particulesClonesPeriodiques->begin();
      found = false;
      
      // Recherche si le clone existe deja sur ce processeur
      // Si oui => mise à jour
      while( il != particulesClonesPeriodiques->end() && !found )
      {
        if ( recvbuf_INT[NB_INT_PART*j] == (*il)->getPeriodicReferenceID() )
        {
	  (*il)->setPosition( &recvbuf_DOUBLE[NB_DOUBLE_PART*j+10] );
	  Vecteur trans( recvbuf_DOUBLE[NB_DOUBLE_PART*j],
	  	recvbuf_DOUBLE[NB_DOUBLE_PART*j+1],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+2] );
	  (*il)->setVitesseTranslation( trans );
	  (*il)->setQuaternionRotation( recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
	 	recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+5],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+6] );
	  Vecteur rot( recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
	  	recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+9] );
	  (*il)->setVitesseRotation( rot );
	  if ( AdamsBashforth )
	    (*il)->setCinematiqueNm2( &recvbuf_DOUBLE[NB_DOUBLE_PART*j+26] );
	  found = true;		
        }
	else il++;      
      }
      
      // Si non => creation du clone periodique si il appartient à ce processeur
      if (!found)
	if ( LC->isInLinkedCell( recvbuf_DOUBLE[NB_DOUBLE_PART*j+22],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+23],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+24] ) )
	{
          if ( Grains_Exec::m_MPI_verbose )
	  {
	    ostringstream oss;
	    oss << "   t=" << Grains_Exec::doubleToString(time,TIMEFORMAT) 
		<< " Create Periodic Clone                       Ref Id = " 
		<< recvbuf_INT[NB_INT_PART*j]
		<< " Classe = " << recvbuf_INT[NB_INT_PART*j+1] << " " 
		<< recvbuf_DOUBLE[NB_DOUBLE_PART*j+22] << " " 
		<< recvbuf_DOUBLE[NB_DOUBLE_PART*j+23] << " " 
		<< recvbuf_DOUBLE[NB_DOUBLE_PART*j+24]
		<< endl;
	    MPIWrapperGrains::addToMPIString(oss.str());
	  }

          ObstaclePeriodique const* pobs = 
		LC->getObstaclePeriodique( recvbuf_INT[NB_INT_PART*j+2] );	
	  Particule *newClonePer = 
		new ParticulePeriodique( recvbuf_INT[NB_INT_PART*j],
      		(*ParticuleClassesReference)[recvbuf_INT[NB_INT_PART*j+1]],
		pobs,
		*(pobs->getPeriode()),
		1,
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+1],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+2],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+5],		
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+9],		
      		&recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],		
		COMPUTE,
		0 );	
		
	  if ( AdamsBashforth )
	    newClonePer->setCinematiqueNm2(
	    	&recvbuf_DOUBLE[NB_DOUBLE_PART*j+26] );

          // Ajout de la particule dans le LinkedCell
          LC->Link( newClonePer );  

          // Ajout dans les differentes listes de EnsComposants
          particulesClonesPeriodiques->push_back(newClonePer);
	}	  
    } 

    delete [] numClassRang;
    delete [] features;    
    delete [] recvcounts;
    delete [] recvcounts_INT;
    delete [] recvcounts_DOUBLE;    
    delete [] displs;
    delete [] displs_INT;
    delete [] displs_DOUBLE;    
    delete [] recvbuf_INT;
    delete [] recvbuf_DOUBLE; 
  }
  MPI_Barrier( m_MPI_COMM_activProc );  
}




// ----------------------------------------------------------------------------
// AllGather de listes d'entiers avec le comm periodique 
void MPIWrapperGrains::commPer_AllGatherGlobal_listINT( Scalar time,
  	list<int> &IDs )
{
  list<int>::iterator il;  
  int i,j;
  
  if ( m_isInCommPeriodic )
  {
    // Taille des messages proportionnel au nombre d'elements de la liste
    // ------------------------------------------------------------------
    int nbINT = int(IDs.size());

    // Taille des messages à passer
    int *recvcounts = new int[m_nprocs_per];
    for (i=0; i<m_nprocs_per; ++i) recvcounts[i]=0;
    MPI_Allgather( &nbINT, 1, MPI_INT, recvcounts, 1, MPI_INT, 
    	*m_MPI_COMM_Periodic ); 

    // Partition et taille du buffer de réception
    int *displs = new int[m_nprocs_per];
    displs[0] = 0; 
    for (i=1; i<m_nprocs_per; ++i) displs[i] = displs[i-1]+recvcounts[i-1];  
    int recvsize = displs[m_nprocs_per-1]+recvcounts[m_nprocs_per-1];
    int *recvbuf = new int[recvsize];

  
    // Communication des entiers 
    // -------------------------
    // Copie dans un vecteur local
    int *vecIDs = new int[nbINT];
    for (il=IDs.begin(),i=0;il!=IDs.end();il++,++i)
      vecIDs[i] = *il;  
  
    // Communication des entiers
    MPI_Allgatherv( vecIDs, nbINT, MPI_INT, recvbuf, 
  	recvcounts, displs, MPI_INT, *m_MPI_COMM_Periodic );


    // Traitement des infos recues
    // ---------------------------
    IDs.clear();
    for (j=0;j<recvsize;++j) IDs.push_back(recvbuf[j]);

    delete [] vecIDs;
    delete [] recvcounts;
    delete [] displs;
    delete [] recvbuf; 
  }
  
  MPI_Barrier( m_MPI_COMM_activProc );    
} 




// ----------------------------------------------------------------------------
// AllGather de listes d'entiers avec le comm periodique 
void MPIWrapperGrains::commPer_AllGatherGlobal_listINT_unique(Scalar time,
  	list<int> &IDs)
{
  MPIWrapperGrains::commPer_AllGatherGlobal_listINT( time, IDs );
  
  if ( m_isInCommPeriodic )
  {    
    IDs.sort();
    IDs.unique();
  }
  
  MPI_Barrier( m_MPI_COMM_activProc );    
} 




// ----------------------------------------------------------------------------
// Envoi des forces de contact des clones periodiques a leur reference
void MPIWrapperGrains::commParticulesClonesPer_AllGatherGlobal_Forces(
	Scalar time,
  	set<Particule*>* particulesReferencesPeriodiques,
  	list<Particule*>* particulesClonesPeriodiques ) const
{
  set<Particule*>::const_iterator is;
  list<Particule*>::const_iterator il;  
  int i,j;
  bool found=false;

  if ( m_isInCommPeriodic )
  {
    // Taille des messages proportionnel au nombre de particules dans 
    // particulesClonesPeriodiques
    // --------------------------------------------------------------
    int nb_prp = 0;
    for (il=particulesClonesPeriodiques->begin();
    	il!=particulesClonesPeriodiques->end();il++)
      if ( (*il)->getTag() != 2 ) ++nb_prp;    

    // Taille des messages à passer
    int *recvcounts = new int[m_nprocs_per];
    for (i=0; i<m_nprocs_per; ++i) recvcounts[i]=0;
    MPI_Allgather( &nb_prp, 1, MPI_INT, recvcounts, 1, MPI_INT,
    	*m_MPI_COMM_Periodic ); 

    // Partition et taille du buffer de réception
    int *displs = new int[m_nprocs_per];
    displs[0] = 0; 
    for (i=1; i<m_nprocs_per; ++i) displs[i] = displs[i-1] + recvcounts[i-1];  
    int recvsize = displs[m_nprocs_per-1] + recvcounts[m_nprocs_per-1];


    // Communication des entiers: 
    // Ordre par particule: 
    // [numero de particule, classe, rang expéditeur]  
    // ----------------------------------------------
    int NB_INT_PART = 3;
    int *numClassRang = new int[NB_INT_PART*nb_prp];
    for (il=particulesClonesPeriodiques->begin(),i=0;
    	il!=particulesClonesPeriodiques->end();il++)
      if ( (*il)->getTag() != 2 ) 
      {
        numClassRang[i] = (*il)->getPeriodicReferenceID(); 
        numClassRang[i+1] = m_rank;
        numClassRang[i+2] = (*il)->getCoordinationNumber();	
	i += NB_INT_PART;
      }            

    // Communication des doubles: cinématique & configuration 
    // Ordre par particule: 
    // [position, vitesse translation, quaternion rotation, vitesse rotation]
    // ----------------------------------------------------------------------
    int NB_DOUBLE_PART = 6;  
    double *features = new double[NB_DOUBLE_PART*nb_prp];
    for (il=particulesClonesPeriodiques->begin(),i=0;
    	il!=particulesClonesPeriodiques->end();il++)
      if ( (*il)->getTag() != 2 )
      { 
        (*il)->copyForceMoment( features, i );
	i += NB_DOUBLE_PART;
      }


    // Partition et taille du buffer de réception
    int *recvcounts_INT = new int[m_nprocs_per];
    for (i=0; i<m_nprocs_per; ++i) 
      recvcounts_INT[i] = NB_INT_PART * recvcounts[i];
    int *displs_INT = new int[m_nprocs_per];
    displs_INT[0] = 0; 
    for (i=1; i<m_nprocs_per; ++i) 
      displs_INT[i] = displs_INT[i-1] + recvcounts_INT[i-1];  
    int recvsize_INT = displs_INT[m_nprocs_per-1] 
    	+ recvcounts_INT[m_nprocs_per-1];
    int *recvbuf_INT = new int[recvsize_INT];
  
    // Communication des entiers: 2 integer par particule
    MPI_Allgatherv( numClassRang, NB_INT_PART * nb_prp, MPI_INT, recvbuf_INT, 
  	recvcounts_INT, displs_INT, MPI_INT, *m_MPI_COMM_Periodic );

  
    // Partition et taille du buffer de réception
    int *recvcounts_DOUBLE = new int[m_nprocs_per];
    for (i=0; i<m_nprocs_per; ++i) 
      recvcounts_DOUBLE[i] = NB_DOUBLE_PART * recvcounts[i];
    int *displs_DOUBLE = new int[m_nprocs_per];
    displs_DOUBLE[0] = 0; 
    for (i=1; i<m_nprocs_per; ++i) 
      displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];  
    int recvsize_DOUBLE = displs_DOUBLE[m_nprocs_per-1] 
    	+ recvcounts_DOUBLE[m_nprocs_per-1];
    double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];
      
    // Communication des force & moment: 6 double par particule
    MPI_Allgatherv( features, NB_DOUBLE_PART * nb_prp, MPI_DOUBLE, 
    	recvbuf_DOUBLE, recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE, 
	*m_MPI_COMM_Periodic ); 


    // Traitement des infos recues
    // ---------------------------
    for (j=0;j<recvsize;++j)
    {
      found = false;
      is = particulesReferencesPeriodiques->begin();
      
      // Recherche si la reference periodique existe sur ce processeur
      // Si oui => ajout des force & moment
      while( is!=particulesReferencesPeriodiques->end() && !found )
        if ( (*is)->getID() == recvbuf_INT[NB_INT_PART*j] )
	  found = true;
	else is++;

      if ( found )
      {			
        (*is)->addForceMoment( recvbuf_DOUBLE[NB_DOUBLE_PART*j],
	  	recvbuf_DOUBLE[NB_DOUBLE_PART*j+1],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+2],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
	  	recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+5] );
	(*is)->addToCoordinationNumber( recvbuf_INT[NB_INT_PART*j+2] );
      } 
    } 

    delete [] numClassRang;
    delete [] features;    
    delete [] recvcounts;
    delete [] recvcounts_INT;
    delete [] recvcounts_DOUBLE;    
    delete [] displs;
    delete [] displs_INT;
    delete [] displs_DOUBLE;    
    delete [] recvbuf_INT;
    delete [] recvbuf_DOUBLE; 
  }
  
  MPI_Barrier( m_MPI_COMM_activProc );  
}




//-----------------------------------------------------------------------------
// Definition du communicateur pour les simulations periodiques
void MPIWrapperGrains::setCommPeriodic( bool const& hasPeriodicObstacle )
{
  // Communication à tous les proc des numéros de ceux contenant des obstacles
  // periodiques
  int i,numHasPer = hasPeriodicObstacle ? m_rank : -1;
  int *recvcounts = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts[i]=0;
  MPI_Allgather( &numHasPer, 1, MPI_INT, recvcounts, 1, MPI_INT,
  	m_MPI_COMM_activProc );
 
  // Regroupement des numéros de proc contenant des obstacles periodiques
  int j=0;
  m_nprocs_per=0;
  for (i=0; i<m_nprocs; ++i) if ( recvcounts[i] != -1 ) ++m_nprocs_per;
  int *numProcPer = new int[m_nprocs_per];
  for (i=0; i<m_nprocs; ++i) 
    if ( recvcounts[i] != -1 )
    {
      numProcPer[j] = recvcounts[i];
      ++j;
    }
  m_isInCommPeriodic = hasPeriodicObstacle;
  
  // Creation du groupe et du communicateur pour la simulation periodique
  // Ce communicateur regroupe tous les proc contenant des obstacles periodiques
  MPI_Group activ_group;
  MPI_Comm_group( m_MPI_COMM_activProc, &activ_group );
  m_MPI_GROUP_Periodic = new MPI_Group; 
  MPI_Group_incl( activ_group, m_nprocs_per, numProcPer, m_MPI_GROUP_Periodic );
  m_MPI_COMM_Periodic = new MPI_Comm;
  MPI_Comm_create( m_MPI_COMM_activProc, *m_MPI_GROUP_Periodic, 
  	m_MPI_COMM_Periodic );
  MPI_Group_free( &activ_group );

  delete [] recvcounts;
  delete [] numProcPer;
}  
