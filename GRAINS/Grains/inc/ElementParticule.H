#ifndef _ElementParticule
#define _ElementParticule

#include "Convex.H"
#include "Particule.H"
#include "ObstaclePeriodique.hh"

class Particule;
class CompParticule;
//class PointContact;

/** @brief Gestion des particules elementaires.

    Constituants elementaires d'une particule composite.

    !!! Remarque : on ne se
    sert jamais du quaternion de rotation d'une particule elementaire mais
    seulement de sa transformation; ce denier n'est donc jamais mis a jour ni en
    sequentiel ni en parallele !!!

    @author D. RAKOTONIRINA - IFP Energies nouvelles - 2014 - Creation */
// ============================================================================
//
class ElementParticule : public Particule
{
public:
  /** @name Constructeur */
  //@{

  /** @brief Constructeur par defaut
  @param autonumbering numerotation automatique ou non */
  ElementParticule( const bool &autonumbering = false );

  /** @brief Constructeur
  @param root Noeud "<Particule>"
  @param pc classe de particules
  @param masterComposite_ composite au niveau superieur
  dont la particule elementaire fait partie */
  ElementParticule( DOMNode* root, const int &pc,
      CompParticule* masterComposite_, const int &elem_id );

  /** @brief Constructeur pour le reload
  @param fileSave Flux de sauvegarde
  @param type particle type
  @param masterComposite_ composite au niveau superieur
  dont la particule elementaire fait partie */
  ElementParticule( istream &fileSave, string &type,
      CompParticule* masterComposite_, const int &elem_id );  

  /** @brief Constructeur avec arguments
  @param id_ numero
  @param ParticuleRef particule de r�f�rence
  @param vx composante x de la vitesse de translation
  @param vy composante y de la vitesse de translation
  @param vz composante z de la vitesse de translation
  @param rx composante x de la vitesse de rotation
  @param ry composante y de la vitesse de rotation
  @param rz composante z de la vitesse de rotation
  @param qrotationx composante x du quaternion de rotation
  @param qrotationy composante y du quaternion de rotation
  @param qrotationz composante z du quaternion de rotation
  @param qrotations scalaire du quaternion de rotation
  @param m position & configuration de la particule
  @param activ activite
  @param tag_ tag
  @param coordination_number_ nombre de contacts de la particule */
  ElementParticule( const int &id_, Particule const* ParticuleRef,
	const double &vx, const double &vy, const double &vz,
	const double &qrotationx, const double &qrotationy,
	const double &qrotationz, const double &qrotations,
	const double &rx, const double &ry, const double &rz,
	const Scalar m[16],
	const ParticuleActivity &activ,
	const int &tag_,
	const int &coordination_number_ = 0 );

  /** @brief Constructeur avec arguments : creer une copie et affecte le
  pointeur sur le master composite
  @param copie particule de reference
  @param masterComposite_ composite au niveau superieur   */
  ElementParticule( const Particule &copie,
  	CompParticule* masterComposite_ );

  /** @brief Destructeur */
  virtual ~ElementParticule();
  //@}


  /**@name Methods */
  //@{
  /** @brief Ajout d'une force au torseur des efforts
  @param force la nouvelle force
  @param point point d'action de la force */
  void addForce( const Point &point, const Vecteur &force );


  /** @brief Add contact force on each composite particle for PP purposes
  @param force la nouvelle force  */
  void addContactForcePP( const Vecteur &force );


  /** @brief Ajoute une force et un moment au torseur des efforts exerces sur la
  particule (utile en simulation periodique)
  @param fx composante x de la force
  @param fy composante y de la force
  @param fz composante z de la force
  @param mx composante x du moment
  @param my composante y du moment
  @param mz composante z du moment */
  void addForceMoment( const double &fx, const double &fy,
  	const double &fz, const double &mx, const double &my,
	const double &mz ) ;

  /** @brief Ajout d'une force s'exercant au centre de gravite au torseur des
  efforts
  @param force la nouvelle force */
  void addBodyForce( const Vecteur &force );

  /** @brief Ajouter un moment
  @param moment La nouvelle force */
  void addMoment( const Vecteur &moment );

  /** Ajoute un nombre de contacts au nombre de contacts de la particule;
  Utilisation: ajoute a la particule de reference periodique les contacts
  de son clone periodique */
  void addToCoordinationNumber( int const& nc );

  /** @brief Est ce que le composant est une particule elementaire ?
  ( Voir CompParticule ) */
  bool isElementParticule() const {return ( true ); }

  /** @brief Gestion des contacts en simmulation PERIODIQUE SEQ:
  Reference du composant lors du calcul des forces de contact
  @return this si contact avec composant autre que CompParticule
  sinon m_masterComposite si contact avec CompParticule */
  Composant* ReferenceComposant();
  //@}


  /**@name Methods Get */
  //@{
  /** @brief Acces au type de convex */
  ConvexType getConvexType() const;

  /** @brief Points decrivant l'enveloppe du Convex.
  @return Un vecteur de points decrivant l'enveloppe du convexe. */
  vector<Point> getEnveloppe() const;

  //@}


private:
  /**@name Parameters */
  //@{
  CompParticule* m_masterComposite; /**< composite au niveau superieur
  	dont la particule elementaire fait partie */
  //@}

  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut interdit,
  une particule periodique n'existe pas sans reference */
  ElementParticule();

  /** @brief Constructeur par Copie.
  Le torseur des forces appliquees est vide.
  Les types de Forme et de Cinematique sont instancies.
  @param copie Particule de reference. */
  ElementParticule( const ElementParticule &copie );
  //@}
};

#endif
