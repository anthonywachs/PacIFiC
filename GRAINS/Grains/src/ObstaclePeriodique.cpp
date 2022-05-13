#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "ObstaclePeriodique.hh"
#include "Particule.H"
#include "ParticulePeriodique.hh"
#include "CompParticulePeriodique.hh"
#include "SaveTable.H"
#include "App.H"
#include "LinkedCell.H"
#include "Cellule.H"
#include "EnsComposant.H"
#include <assert.h>


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur par defaut
ObstaclePeriodique::ObstaclePeriodique( const string &s ) :
  MonObstacle( s )
{
  Composant::m_nb--;
  m_id = -4;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
ObstaclePeriodique::ObstaclePeriodique( DOMNode* root,
	const string  &name,
	const Vecteur &direction ) :
  MonObstacle(),
  m_vecteur( direction ),
  m_associe( NULL )
{
  m_nom = name;

  // Convex - Position & Orientation
  m_geoFormeVdw = new FormeVdW( root );

  // Materiau inutile pour un obstacle periodique
  m_nomMateriau = "periode";
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructeur
ObstaclePeriodique::~ObstaclePeriodique()
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Deplacement de l'obstacle
list<MonObstacle*> ObstaclePeriodique::Deplacer( double temps, double dt,
	const bool &b_deplaceCine_Comp,
	const bool &b_deplaceF_Comp )
{
  m_deplace = false;
  list<MonObstacle*> obstacleEnDeplacement;
  return obstacleEnDeplacement;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Acces a l'obstacle associe
ObstaclePeriodique* ObstaclePeriodique::getAssocie() const
{
  return m_associe;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Acces au vecteur de translation periodique
Vecteur const* ObstaclePeriodique::getPeriode() const
{
  return &m_vecteur;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Determination du contact avec la Particule
void ObstaclePeriodique::InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC )
  throw (ErreurContact)
{}




// ----------------------------------------------------------------------------
// Determination du contact avec une particule
void ObstaclePeriodique::InterActionPostProcessing( Composant* voisin,
	list<struct PointForcePostProcessing>* listOfContacts )
  throw (ErreurContact)
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reload de l'obstacle & publication dans son referent
void ObstaclePeriodique::reload( Obstacle &obstacle, istream &file )
{
  string buffer, adresse;
  Scalar buf = 0.;

  file >> adresse >> m_id;
  SaveTable::create( adresse, this );
  file>> buffer >> m_nomMateriau;
  m_geoFormeVdw = new FormeVdW( file );
  file >> buffer;
  if ( buffer == "*Couleur" )
    file >> buf >> buf >> buf;
  else if ( buffer == "*ToFluid" )
    file >> m_transferToFluid;
  m_geoFormeVdw->readPosition( file );
  file >> buffer >> m_vecteur >> buffer;
  obstacle.append( this );
  file >> buffer;
  assert(buffer == "</Periodique>");
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajout de la reference a l'ObstaclePeriodique associe
void ObstaclePeriodique::setAssocie( ObstaclePeriodique *obstacle )
{
  m_associe = obstacle;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde De l'Obstacle pour Reload
// D. RAKOTONIRINA - Dec. 2014 - Modification
void ObstaclePeriodique::write( ostream &fileSave, Composant const* composant )
    const
{
  fileSave << "<Periodique>" << '\t' << m_nom << '\n';
  fileSave << this << '\t' << m_id << '\n';
  fileSave << "*Materiau\n" << m_nomMateriau << '\n';
  m_geoFormeVdw->writeStatique( fileSave );
  fileSave << "*ToFluid\n" << m_transferToFluid << '\n';
  m_geoFormeVdw->writePosition( fileSave );
  fileSave << endl;
  fileSave << "<Vecteur> " << m_vecteur << "</Vecteur>" << '\n';
  fileSave << "</Periodique>\n";
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Traitement des particules periodiques: creation et destruction des
// clones periodiques de ces particules
void ObstaclePeriodique::LinkUpdateParticulesPeriodiques(
  	list<Particule*>* particulesClonesPeriodiques,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC )
{
  list<Cellule*>::iterator icell;
  list<Particule*>::iterator particule;
  list<App*>::iterator app;

  for (icell=m_inCells.begin();icell!=m_inCells.end();icell++)
  {
    list<Particule*> particulesInCell = *((*icell)->getParticules());
    for (particule=particulesInCell.begin();particule!=particulesInCell.end();
    	particule++)
    {
      if ( (*particule)->getID() >= 0 )
      {
	// Contact entre la particule et l'obstacle periodique
        if ( (*particule)->isContactVdW( this ) )
        {
          // Si le clone n'existe pas: creation du clone
	  if ( !(*particule)->periodicCloneExistence( m_vecteur ) )
	  {
	    // Creation
	    if ( (*particule)->isCompParticule() )
	    {
	      Particule* perclone = new CompParticulePeriodique( *particule,
		  (*ParticuleClassesReference)[
		  (*particule)->getParticuleClasse()], this, m_vecteur, 1 );
	      // Link avec le LinkedCell
              LC->Link( perclone );

              // Ajout dans les differentes listes de EnsComposants
              particulesClonesPeriodiques->push_back(perclone);

	      // Ajout dans la liste des clones de la particule de reference
	      (*particule)->addPeriodicClone( perclone );
	    }
	    else
	    {
	      Particule* perclone = new ParticulePeriodique( *particule,
		  (*ParticuleClassesReference)[
		  (*particule)->getParticuleClasse()], this, m_vecteur, 1 );
	      // Link avec le LinkedCell
              LC->Link( perclone );

              // Ajout dans les differentes listes de EnsComposants
              particulesClonesPeriodiques->push_back(perclone);

	      // Ajout dans la liste des clones de la particule de reference
	      (*particule)->addPeriodicClone( perclone );
	    }
	  }
        }
        // Pas de contact entre la particule et l'obstacle periodique
        else
        {
	  // La particule possede un clone periodique: il faut le detruire
	  if ( (*particule)->periodicCloneExistence(m_vecteur) )
	  {
	    // Destruction du clone
            Particule *pdestroy =
	    	(*particule)->getPeriodicCloneAndErase( this );

            if ( pdestroy )
	    {
              // Suppression de la particule dans le LinkedCell
              LC->remove( pdestroy );

              // Suppression des differentes listes
              removeParticuleFromList( *particulesClonesPeriodiques, pdestroy );

              // Destruction de l'objet point�
              delete pdestroy;
	    }

	    // Si la particule est hors domaine periodique: translation de
	    // m_vecteur de la particule
	    // et de ses clones dans les autres directions
	    if ( !isIn( *particule, 0. ) )
	    {
	      (*particule)->Translate( m_vecteur );
	      LC->LinkUpdateActiveParticule( *particule );
	      (*particule)->translateAndUpdatePeriodicClones( m_vecteur, LC );
	    }
	  }
        }
      }
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Traitement des particules periodiques
void ObstaclePeriodique::LinkUpdateParticulesPeriodiques_MPI(
  	Scalar time,
  	list<int>& ClonestoDestroy,
  	list<int>& ClonestoParticules,
  	list<int>& PartRefPerHalozone,
  	list<int>& PartRefPerOutDomainHalozone,
  	list<int>& InNotRefPerHalozone,
  	set<Particule*>* particulesReferencesPeriodiques,
	list<Particule*>* ENSparticules,
  	list<Particule*>* particulesHalozone,
	LinkedCell* LC )
{
  list<Cellule*>::iterator icell;
  list<Particule*>::iterator particule;
  Particule *pdestroy = NULL;

  for (icell=m_inCells.begin();icell!=m_inCells.end();icell++)
  {
    list<Particule*> particulesInCell = *((*icell)->getParticules());
    for (particule=particulesInCell.begin();particule!=particulesInCell.end();
    	particule++)
    {
      // Cas des particules actives
      if ( (*particule)->getID() >= 0 && (*particule)->getTag() != 2 )
      {
	// Contact entre la particule et l'obstacle periodique
        if ( (*particule)->isContactVdW( this ) )
        {
          // Si le tag de la particule est 1, infos � transmettre aux autres
	  // proc
	  if ( (*particule)->getTag() == 1 )
	  {
	    PartRefPerHalozone.push_back((*particule)->getID());
	    PartRefPerHalozone.push_back(m_id);
	  }
        }
        // Pas de contact entre la particule et l'obstacle periodique
        else
        {
	  // Si la particule est hors domaine: destruction de la particule
	  if ( !isIn( *particule, time ) )
	  {
	    if ( Grains_Exec::m_MPI_verbose )
	    {
	      ostringstream oss;
	      Point const* gc = (*particule)->getPosition();
              oss << "   t=" << Grains_Exec::doubleToString(time,TIMEFORMAT)
      		<< " Periodic reference out of domain            Id = " <<
      		(*particule)->getID() << "  " <<
		(*gc)[X] << " " << (*gc)[Y] << " " << (*gc)[Z] << endl;
      	      MPIWrapperGrains::addToMPIString( oss.str() );
	    }
	    pdestroy = *particule;

	    // Infos � transmettre aux autres procs si tag = 1
	    if ( (*particule)->getTag() == 1 )
	      PartRefPerOutDomainHalozone.push_back((*particule)->getID());

            // Suppression de la particule dans le LinkedCell
            LC->remove( pdestroy );

            // Suppression des differentes listes
            removeParticuleFromSet( *particulesReferencesPeriodiques,
	    	pdestroy );
            removeParticuleFromList( *ENSparticules, pdestroy );
	    if ( pdestroy->getTag() == 1 )
	      removeParticuleFromList( *particulesHalozone, pdestroy );

	    // Referencer la particule pour changer le statut des clones
	    // de cette particule
	    ClonestoParticules.push_back((*particule)->getID());

            // Destruction de l'objet point�
            delete pdestroy;
	  }
	  // Sinon,
	  // si la particule a des clones: supprimer le numero de l'obstacle
	  // periodique dans ceux repertories dans la particule de reference,
	  // supprimer la particule de la liste des particules de reference
	  // et reference l'ID de la particule pour destruction
	  // sinon ne rien faire
	  else
	  {
	    if ( (*particule)->getNombreClonesPeriodiques() )
	    {
	      if ( Grains_Exec::m_MPI_verbose )
	      {
	        ostringstream oss;
	        Point const* gc = (*particule)->getPosition();
                oss << "   t="
	      	<< Grains_Exec::doubleToString( time, TIMEFORMAT )
      		<< " Particule is not periodic reference anymore Id = " <<
      		(*particule)->getID() << "  " <<
		(*gc)[X] << " " << (*gc)[Y] << " " << (*gc)[Z] << endl;
      	        MPIWrapperGrains::addToMPIString( oss.str() );
              }

	      // supprimer le numero de l'obstacle periodique
	      (*particule)->erasePeriodicObstacleID( m_id );

	      // Supprimer la particule de la liste des particules de reference
	      removeParticuleFromSet( *particulesReferencesPeriodiques,
	      	*particule );

              // Referencer la particule pour detruire les clones
	      // de cette particule
	      ClonestoDestroy.push_back((*particule)->getID());

	      // Infos � transmettre aux autres procs si tag = 1
	      if ( (*particule)->getTag() == 1 )
	      {
	        InNotRefPerHalozone.push_back((*particule)->getID());
	        InNotRefPerHalozone.push_back(m_id);
	      }
	    }
	  }
        }
      }
    }
  }

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie si la particule est hors ou dans le domaine periodique
bool ObstaclePeriodique::isIn( Particule* reference, Scalar time )
{
  bool in = true;

  Vecteur ObstacleParticule = *reference->getPosition() - *getPosition();
  if ( ObstacleParticule * m_vecteur < 0. ) in = false;

  return in;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Indique si un obstacle periodique intersecte une boite
bool ObstaclePeriodique::obstaclePeridiqueIntersect( BBox const& boite ) const
{
  return intersect( m_geoFormeVdw->Forme::BoxForme(), boite );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajout des nouvelles references periodiques
void ObstaclePeriodique::addNewPeriodicReference_MPI(
  	Scalar time,
  	set<Particule*>* particulesReferencesPeriodiques )
{
  list<Cellule*>::iterator icell;
  list<Particule*>::iterator particule;
  list<App*>::iterator app;

  for (icell=m_inCells.begin();icell!=m_inCells.end();icell++)
  {
    list<Particule*> particulesInCell = *((*icell)->getParticules());
    for (particule=particulesInCell.begin();particule!=particulesInCell.end();
    	particule++)
      if ( (*particule)->getID() >= 0 && (*particule)->getTag() != 2 )
      {
	// Contact entre la particule et l'obstacle periodique
        if ( (*particule)->isContactVdW( this ) )
        {
	  // Si c'est une nouvelle reference periodique:
	  // Ajout dans la liste des clones de la particule de reference
	  // et dans celle des references periodiques
	  if ( !(*particule)->periodicCloneExistence( m_id ) )
	  {
            particulesReferencesPeriodiques->insert(*particule);
	    (*particule)->addPeriodicObstacleID( m_id );
	    if ( Grains_Exec::m_MPI_verbose )
	    {
	      ostringstream oss;
	      Point const* gc = (*particule)->getPosition();
              oss << "   t=" << Grains_Exec::doubleToString( time, TIMEFORMAT )
      		<< " New periodic reference                      Id = " <<
      		(*particule)->getID() << "  " <<
		(*gc)[X] << " " << (*gc)[Y] << " " << (*gc)[Z] << endl;
      	      MPIWrapperGrains::addToMPIString( oss.str() );
	    }
	  }
	}
      }
  }
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
int ObstaclePeriodique::numberOfPoints_PARAVIEW() const
{
  return ( m_geoFormeVdw->getConvex()->numberOfPoints_PARAVIEW() );
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
int ObstaclePeriodique::numberOfCells_PARAVIEW() const
{
  return ( m_geoFormeVdw->getConvex()->numberOfCells_PARAVIEW() );
}




// ----------------------------------------------------------------------------
// Nombre de polygones elementaires pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
list<Point> ObstaclePeriodique::get_polygonsPts_PARAVIEW( Vecteur const*
    translation )
const
{
  return ( getForme()->get_polygonsPts_PARAVIEW( translation ) );
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
void ObstaclePeriodique::write_polygonsPts_PARAVIEW( ostream &f,
  	Vecteur const* translation ) const
{
  getForme()->get_polygonsPts_PARAVIEW( translation );
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
void ObstaclePeriodique::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  m_geoFormeVdw->getConvex()->write_polygonsStr_PARAVIEW(connectivity,
	offsets, cellstype, firstpoint_globalnumber, last_offset);
}
