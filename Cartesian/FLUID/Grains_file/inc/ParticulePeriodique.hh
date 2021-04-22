#ifndef _ParticulePeriodique
#define _ParticulePeriodique

#include "Particule.H"
#include "ObstaclePeriodique.hh"


/** @brief Gestion des particules periodiques. 

    Une particule periodique n'a pas d'existence propre, elle est liee a sa
    particule de reference. 
    Elle existe pendant le temps de contact de la particule de reference avec 
    un obstacle periodique. Elle transmet ses contacts a sa reference et se 
    deplace en fonction du deplacement de celle-ci.
    
    Son numero et son tag sont fixes a -2 et 0 respectivement par convention.

    @author GRAINS Project - IFP - 2007 
    D. RAKOTONIRINA - IFP Energies nouvelles - Mars 2014 - Modification */
// ============================================================================
class ParticulePeriodique : virtual public Particule
{
public:
  /** @name Constructeur */
  //@{
  /** @brief Constructeur avec arguments
  @param reference La particule reference de la periodique 
  @param particuleClasseReference classe de la particule
  @param contactReferenceObstacle Obstacle en contact avec la particule de 
  	reference 
  @param trans Vecteur de periodicite par rapport a la particule de reference 
  @param nbper Nombre de periodicite de la particule */
  ParticulePeriodique(Particule* reference,
  	Particule const* particuleClasseReference,
  	ObstaclePeriodique const* contactReferenceObstacle,
	Vecteur const& trans,
	const int &nbper);

  /** @brief Constructeur avec arguments
  @param IDref numero de la particume de reference
  @param particuleClasseReference classe de la particule
  @param contactReferenceObstacle Obstacle en contact avec la particule de 
  	reference 
  @param trans Vecteur de periodicite par rapport a la particule de reference 
  @param nbper Nombre de periodicite de la particule	  
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
  @param tag_ tag */
  ParticulePeriodique(const int &IDref, 
  	Particule const* particuleClasseReference,
  	ObstaclePeriodique const* contactReferenceObstacle,
	Vecteur const& trans,
	const int &nbper,		
	const double &vx, const double &vy, const double &vz,
	const double &qrotationx, const double &qrotationy, 
	const double &qrotationz, const double &qrotations,	 
	const double &rx, const double &ry, const double &rz,	 
	const Scalar m[16],
	const ParticuleActivity &activ, 
	const int &tag_); 
  
  /** @brief Destructeur */
  virtual ~ParticulePeriodique();
  //@}


  /**@name Methods */
  //@{
  /** @brief Mise à jour de la vitesse & position de la particule periodique */
  void updateVitessePositionPeriodique();
  
  /** @brief Ajout des forces & moments de contact à la particule maitre du
  clone 
  @param temps Temps de simulation.
  @param dt Pas de temps */
  void AddForcesFromPeriodicCloneToParticule(Scalar temps, Scalar dt,
  	 bool ContactforceOutput  ) const;  
  //@}


  /** @name Processus de restauration */
  //@{  
  /** @brief Restauration de l'etat 
  @param Pmemento_ configuration    
  @param Cmemento_ cinematique 
  @param contactReferenceObstacle_ pointeur sur l'obstacle periodique de la 
  	particule de reference 
  @param reference_ particule de reference periodique */
  virtual void restaureState(ConfigurationMemento const* Pmemento_,
  	CineParticuleMemento const* Cmemento_,
	ObstaclePeriodique const* contactReferenceObstacle_ = NULL,
	Particule* reference_ = NULL ); 
	
  /** @brief Definit l'etat de restauration du clone 
  @param done_ restaure ou non */
  virtual void setRestaureState( bool const& done_ ) 
  	{ m_restaureState_done = done_; }
  
  /** @brief Recupere l'etat de restauration du clone */
  virtual bool getRestaureState() { return m_restaureState_done; }  	 
  //@}


  /**@name Methods Get */
  //@{
  /** @brief Renvoie l'obstacle en contact avec la reference */
  ObstaclePeriodique const* getObstacle() const;
  
  /** @brief Renvoie le nombre de periodes de la particule periodique: 1 =
  uni-periodicite, 2 = bi-periodicite, 3 = tri-periodicite */
  int getNbPeriodes() const; 
  
  /** @brief Renvoie le vecteur de periodicite */
  Vecteur const* getVecteurPeriodique() const; 
  
  /** @brief Renvoie la particule de reference */
  Particule* getPeriodicReference() const;
  
  /** @brief Renvoie le numero de la particule de reference */
  int getPeriodicReferenceID() const;        
  //@}

protected:

  Vecteur const* m_periodicTranslation; /**< Vecteur de periodicite par rapport 
  	a la particule de reference */
  
private:
  /** @name Parameters */
  //@{
  Particule* m_reference; /**< Particule de reference */
  int m_referenceID; /**< numero de la particule de reference */  
  ObstaclePeriodique const* m_contactReferenceObstacle; /**< Obstacle en 
  	contact avec la particule de reference */
//  Vecteur const* m_periodicTranslation; /**< Vecteur de periodicite par rapport 
//  	a la particule de reference */
  int m_nperiodes; /**< Nombre de periodicite de la particule: 1 =
  	uni-periodicite, 2 = bi-periodicite, 3 = tri-periodicite */
  bool m_restaureState_done; /**< est l'etat du clone a ete restaure ou non ? */
  //@}


  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut interdit, 
  une particule periodique n'existe pas sans reference */
  ParticulePeriodique();
  //@}
};

#endif
