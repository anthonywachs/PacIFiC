// Gestion de la cinematique pour les objets Particule
// ============================================================================
#ifndef _CineParticule
#define _CineParticule

#include "Cinematique.H"
#include "WriterXML.hh"
#include "Vecteur.H"
#include "Quaternion.H"
using namespace solid;

class CineParticuleMemento;
class Particule;
class LeapFrog_2D;
class LeapFrog_3D;
class TimeIntegrator;
class Torseur;
class CineParticule;


istream &operator >> ( istream &fileIn, CineParticule &cinematique );


/** @brief Gestion de la cinematique de Particule.

    @author G.FERRER - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class CineParticule : public Cinematique
{
public:
  /**@name Constructors */
  //@{
  /** @brief Destructeur */
  virtual ~CineParticule();
  //@}


  /**@name Methods Virtual */
  //@{
  /** @brief Classname for BuilderFactory */
  virtual string className() const = 0;

  /** @brief Construction d'un clone de la cinematique.
  @return Le clone construit. */
  virtual CineParticule* clone() const = 0;

  /** @brief Calcul de la variation de quantite de mouvement au temps precedent
  @param torseur Torseur somme des forces appliquees sur la particule 
  @param dt increment de temps. 
  @param particule Particule associee a la cinematique. */
  virtual void CalculerVariationQDM( const Torseur &torseur,
	Scalar dt,const Particule* particule ) = 0;
			   
  /** @brief Calcul du terme de masse ajouté explicite sur les moments
  @param dw difference explicite de vitesse de rotation 
  @param inertie tenseur d'inertie de la particule */
  virtual Vecteur calculIdwExplicite( const Vecteur &dw,
  	const double *inertie ) const = 0;			   
  //@}


  /** @name Methods */
  //@{
  /** @brief Deplacement de la particule 
  @param particule la particule a deplacer
  @param dt pas de temps granulaire */
  double Deplacer( Particule* particule, double dt );

  /** @brief Calcul de la vitesse relative de l'objet.
  @param om Le bras de levier permettant de definir la partie rotation */
  Vecteur Vitesse( const Vecteur &om ) const;
   
  /** @brief Ecriture de l'objet sur le flux de sortie avec un format de
  precision elevee
  @param fileOut Flux de sortie. */
  void writeCineParticule( ostream &fileOut ) const;
  
  /** @brief Ecriture de l'objet sur le flux de sortie avec un format de
  precision elevee avec le format de reload 2014
  @param fileOut Flux de sortie. */
  void writeCineParticule2014( ostream &fileOut ) const; 
  
  /** @brief Ecriture de l'objet en binaire sur le flux de sortie
  avec le format de reload 2014
  @param fileOut Flux de sortie. */
  void writeCineParticule2014_binary( ostream &fileOut );
  
  /** @brief Lecture de l'objet en binaire sur le flux d'entree
  avec le format de reload 2014
  @param StreamIN Flux d'entree. */
  void readCineParticule2014_binary( istream &StreamIN );       
  
  /** @brief Copie de la cinematique au temps t-2dt: vitesse de translation, 
  vitese de rotation, variation de QDM translationnelle, variation de QDM
  rotationnalle
  @param vit vecteur de copie 
  @param i position dans le vecteur vit */
  void copyCinematiqueNm2( double *vit, int i ) const;
  //@}


  /** @name Accessors */
  //@{
  /** @brief Retour du Quaternion de rotation */
  Quaternion const* getRotation() const;

  /** @brief Vitesse de rotation */
  Vecteur const* getVitesseRotation() const;

  /** @brief Vitesse de translation */
  Vecteur const* getVitesseTranslation() const;
  //@}
  

  /** @name Methods Set */
  //@{
  /** @brief  Mise a zero de la cinematique */
  void reset();

  /** @brief Modification du quaternion de rotation.
  @param vecteur0 Composante 0 du vecteur
  @param vecteur1 Composante 1 du vecteur
  @param vecteur2 Composante 2 du vecteur
  @param scalaire Composante du scalaire */
  void setQuaternionRotation( const Scalar &vecteur0, 
	const Scalar &vecteur1, 
	const Scalar &vecteur2, 
	const Scalar &scalaire );
	
  /** @brief Modification du quaternion de rotation.
  @param qrot quaternion de rotation */
  void setQuaternionRotation( const Quaternion &qrot );	

  /** @brief Modification de la vitesse de rotation.
  @param omega vitesse de rotation */
  void setVitesseRotation( const Vecteur &omega );

  /** @brief Modification de la vitesse de translation.
  @param vitesse vitesse de translation */
  void setVitesseTranslation( const Vecteur &vitesse );
  
  /** @brief Cinematique au temps t-2dt: vitesse de translation, 
  vitese de rotation, variation de QDM translationnelle, variation de QDM
  rotationnalle
  @param tab tableau de 12 doubles contenant les valeurs des 4 vecteurs vitesse 
  de translation, vitese de rotation, variation de QDM translationnelle, 
  variation de QDM rotationnalle */
  void setCinematiqueNm2( double const* tab );  
  //@}

  
  /** @name Processus de restauration */
  //@{
  /** @brief Sauvegarde de l'etat */
  void saveState();
  
  /** @brief Cree et renvoie l'etat */
  CineParticuleMemento* createState();  
  
  /** @brief Restauration de l'etat */
  void restaureState();
  
  /** @brief Restauration de l'etat 
  @param memento_ cinematique */
  void restaureState( CineParticuleMemento const* memento_ );  
  //@}


  /**@name Methods Friends */
  //@{
  /** @brief Operateur de lecture.
  @return Flux emetteur
  @param fileIn Flux emetteur
  @param cinematique Objet courant */
  friend istream &operator >> ( istream &fileIn, CineParticule &cinematique );
  
  /** @brief Operateur d'ecriture.
  @return Flux recepteur
  @param fileOut Flux recepteur
  @param cinematique Objet courant. */
  friend ostream &operator << ( ostream &fileOut, 
	const CineParticule &cinematique );
  //@}
    

protected:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut inaccessible. */
  CineParticule();

  /** @brief Constructeur par copie inaccessible
  @param copie l'objet a copier */
  CineParticule( const CineParticule &copie );
  //@}


  /** @name Parameters */
  //@{
  Vecteur m_vitesseT; /**< Vitesse de translation */
  Vecteur m_vitesseR; /**< Vitesse de rotation */  
  Quaternion m_QuaternionVitesseR; /**< Quaternion de vitesse de rotation. 
  	IMPORTANT: ce quaternion correspond a la derivee en temps du quaternion
        de rotation note m_QuaternionRotation dans cette classe, et est egal a 
  	m_QuaternionVitesseR = 0.5 * [0,m_vitesseR] * m_QuaternionRotation, 
  	ou m_vitesseR est la vitesse de rotation de la particule. 
	Pour recuperer m_vitesseR, il faut effectuer l'operation inverse soit 
  	m_vitesseR = 2 * m_QuaternionVitesseR * m_QuaternionRotation^*, 
	ou m_QuaternionRotation^* est le quaternion conjugue de 
	m_QuaternionRotation. 
	Pour plus de details, voir le rapport IFP 56336 de David PETIT. */
  Quaternion m_QuaternionRotation; /**< Quaternion de rotation */ 
  Vecteur m_dUdt; /**< variation de quantite de mouvement translationnelle 
  	au temps precedent */
  Vecteur m_dOmegadt; /**< variation de quantite de mouvement rotationnelle 
  	au temps precedent */
  Vecteur deplacementTranslationnelOverDt; /**< deplacement translationnel sur
  	le pas de temps */
  Vecteur vitesseRMoyenne; /**< vitesse de rotation moyenne 
  	sur le pas de temps */
  Quaternion QuaternionRotationOverDt; /**< quaternion de rotation
  	sur le pas de temps */  				 
  TimeIntegrator* m_timeIntegrationScheme; /**< schema d'integration en temps */
  CineParticuleMemento* m_memento; /**< pour sauver l'état */ 
  //@}
};

#endif
