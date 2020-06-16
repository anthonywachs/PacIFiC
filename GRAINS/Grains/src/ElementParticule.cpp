#include "ElementParticule.H"
#include "Memento.hh"
#include "Cinematique_BuilderFactory.H"
#include "Convex_BuilderFactory.H"
#include "FormeVdW.H"
#include "Forme.H"
#include "Convex.H"
#include "Contact_BuilderFactory.hh"
#include "CompParticule.hh"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur par defaut
ElementParticule::ElementParticule() :
  Particule( false )
{
  m_id = -3;
  m_masterComposite = NULL ;
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
ElementParticule::ElementParticule( DOMNode* root, const int &pc,
  CompParticule* masterComposite_, const int &elem_id ) :
  Particule( root, false, pc )
{
  m_masterComposite = masterComposite_ ;
  m_id = -10 - elem_id ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur ( reload )
ElementParticule::ElementParticule( istream &fileSave, string &type,
  CompParticule* masterComposite_, const int &elem_id ) :
  Particule( *masterComposite_ ) /* constructeur par copie */
{
  string buffer, buff;
  m_masterComposite = masterComposite_ ;
  m_id = elem_id ;
  m_geoFormeVdw = new FormeVdW( fileSave, type );
  fileSave >> buffer >> buff >> buff >> buff;
  m_geoFormeVdw->readPosition( fileSave );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur
ElementParticule::ElementParticule( const int &id_,
	Particule const* ParticuleRef,
	const double &vx, const double &vy, const double &vz,
	const double &qrotationx, const double &qrotationy,
	const double &qrotationz, const double &qrotations,
	const double &rx, const double &ry, const double &rz,
	const Scalar m[16],
	const ParticuleActivity &activ,
	const int &tag_,
	const int &coordination_number_ ) :
  Particule( id_, ParticuleRef, vx, vy, vz,
	qrotationx, qrotationy,
	qrotationz, qrotations,
	rx, ry, rz, m, activ, tag_, coordination_number_ )
{
}




// ----------------------------------------------------------------------------
// Constructeur par copie
ElementParticule::ElementParticule( const Particule &copie,
  	CompParticule* masterComposite_ ) :
  Particule( copie )
{
  m_id = copie.getID();
  m_masterComposite = masterComposite_ ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructeur
ElementParticule::~ElementParticule()
{
}




// ----------------------------------------------------------------------------
// Ajout d'une force au torseur des efforts
void ElementParticule::addForce( const Point &point, const Vecteur &force )
{
  m_masterComposite->addForce( point, force );
}




// ----------------------------------------------------------------------------
// Add contact force on each composite particle for postprocessing purposes
void ElementParticule::addContactForcePP( const Vecteur &force )
{
  m_masterComposite->addContactForcePP( force );
  cout << "TOTO" << endl;
}




// ----------------------------------------------------------------------------
// Ajout d'une force s'exercant au centre de gravite au torseur des efforts
void ElementParticule::addBodyForce( const Vecteur &force )
{
  m_masterComposite->addBodyForce( force );
}




// ----------------------------------------------------------------------------
// Ajouter un moment
void ElementParticule::addMoment( const Vecteur &moment )
{
  m_masterComposite->addMoment( moment );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajoute une force et un moment au torseur des efforts exerces sur la
// particule (!!! utile en simulation periodique !!!)
void ElementParticule::addForceMoment( const double &fx, const double &fy,
	const double &fz, const double &mx, const double &my, const double &mz )
{
  Vecteur force( fx, fy, fz );
  m_masterComposite->addBodyForce( force );
  Vecteur torque( mx, my, mz );
  m_masterComposite->addMoment( torque );
}




// ----------------------------------------------------------------------------
// Ajoute un nombre de contacts au nombre de contacts de la particule;
// Utilisation: ajoute a la particule de reference periodique les contacts
// de son clone periodique
void ElementParticule::addToCoordinationNumber( int const& nc )
{
  m_masterComposite->addToCoordinationNumber( nc );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajoute une force et un moment au torseur des efforts exerces sur la
// particule (!!! utile en simulation periodique sequentielle !!!)
Composant* ElementParticule::ReferenceComposant()
{
  return( m_masterComposite );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Acces au type de convex de la particule elementaire
ConvexType ElementParticule::getConvexType() const
{
  return( getForme()->getConvex()->getConvexType() );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Points decrivant l'enveloppe du Convex
vector<Point> ElementParticule::getEnveloppe() const
{
  return( getForme()->getConvex()->getEnveloppe() );
}
