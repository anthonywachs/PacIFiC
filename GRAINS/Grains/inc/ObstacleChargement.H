#ifndef _ObstacleChargement
#define _ObstacleChargement

#include "Vecteur.H"
using namespace solid;
#include <list>
#include <string>
#include <iostream>
using namespace std;
#include "ReaderXML.hh"


class ObstacleChargement;
bool operator < ( const ObstacleChargement &c0,
	const ObstacleChargement &c1 );
ostream &operator << ( ostream &fileOut, ObstacleChargement &chargement );
istream &operator >> ( istream &fileIn, ObstacleChargement &chargement );


/** @brief Description du chargement sur une obstacle.

    Cet objet permet de stocker les chargements qui seront appliquees sur la 
    obstacle durant la simulation, en particulier les eventuels deplacements.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation */
// ============================================================================
class ObstacleChargement
{
public:
  /**@name Constructeurs */
  //@{
  /** @brief Constructeur par defaut. */
  ObstacleChargement();

  /** @brief Constructeur par copie
  @param copie L'objet de reference. */
  ObstacleChargement( const ObstacleChargement &copie );

  /** @brief Creation d'un chargement d'Obstacle
  @param root Le noeud "<Chargement>" 
  @param dt Pas de temps de la simulation 
  @param rank rang du processus */
  ObstacleChargement( DOMNode* root, double dt, int rank );

  /** @brief Destructeur. */
  ~ObstacleChargement();
  //@}


  /**@name Methods */
  //@{
  /** @brief Nom de l'obstacle sur lequel on applique le deplacement.
  @return L'identificateur */
  string getNom() const;

  /** @brief Temps de chargement actif dans l'espace indique.
  @return Le temps de chargement.
  @param debut Temps de debut
  @param fin Temps de fin */
  Scalar getTime( double debut, double fin ) const;

  /** @brief Le chargement d'obstacle est il actif ?
  @return Vrai si actif
  @param t temps physique
  @param dt Increment de temps */
  bool isActif( Scalar t, Scalar dt ) const;
  
  /** @brief Le chargement d'obstacle est il termine ?
  @return Vrai si termine
  @param t temps physique
  @param dt Increment de temps */
  bool isCompleted( Scalar t, Scalar dt ) const;  

  /** @brief Vitesse de deplacement en translation 
  @param temps Temps de debut
  @param dt Increment de temps */
  Vecteur const* VitesseTranslation( Scalar temps, Scalar dt );

  /** @brief Vitesse de rotation instantanee
  @param temps Temps de debut
  @param dt Increment de temps */
  Vecteur const* VitesseRotation( Scalar temps, Scalar dt );
  
  /** @brief Deplacement en translation sur un pas de temps à un temps donné
  @param temps Temps de debut
  @param dt Increment de temps
  @return Le deplacement demandee. */
  Vecteur DeplacementTranslation( Scalar temps, Scalar dt );  

  /** @brief Deplacement en rotation sur un pas de temps à un temps donné
  @param temps Temps de debut
  @param dt Increment de temps
  @return Le deplacement demandee. */
  Vecteur DeplacementRotation( Scalar temps, Scalar dt ); 

  /** @brief Decalage dans le temps du chargement.
  @param decalage La valeur du decalage. */
  void TempsSlip( Scalar decalage );

  /** @brief Debug
  @param c message de debug */
  void debug( char *c ); 
  
  /** @brief Type de chargement */
  string getType() const; 

  /** @brief Renvoie le signe d'un double
  @param v double a tester */
  Scalar sign( Scalar v ) const;
  //@}


  /**@name Operators */
  //@{
  /** @brief Identite entre deux chargements.
  @return Vrai si le deux chargements ont meme adresse.
  @param chargement Chargement a comparer. */
  bool operator == ( const ObstacleChargement &chargement ) const;
  //@}


  /**@name Methods Friends */
  //@{
  /** @brief Ordre entre deux chargements.
  Un chargement est plus grand si son temps de depart est superieur.
  Cet operator permet un "sort" sur la liste des chargements.
  @return Vrai si c0.tdebut < c1.tdebut */
  friend bool operator < ( const ObstacleChargement &c0,
	const ObstacleChargement &c1 );

  /** @brief Operateur d'ecriture
  @return Flux emetteur
  @param fileOut   Flux recepteur
  @param chargement Objet courant */
  friend ostream &operator << ( ostream &fileOut, 
	ObstacleChargement &chargement );
	
  /** @brief Operateur de lecture
  @return Flux emetteur
  @param fileIn    Flux emetteur
  @param chargement Objet courant */
  friend istream &operator >> ( istream &fileIn, 
	ObstacleChargement &chargement );
  //@}


private:
  /**@name Parameters */
  //@{
  string m_nomObstacle; /**< Nom de l'obstacle associe au chargement */
  string m_type; /**< Type de chargement */
  double m_tdebut; /**< Temps de debut */
  double m_tfin; /**< Temps de fin */
  Vecteur m_vitesse_translation; /**< Vitesse de translation */
  Vecteur m_vitesse_rotation; /**< Vitesse de rotation */
  Scalar m_Sin_amplitude; /**< Amplitude du mouvement sinusoidal */
  Scalar m_Sin_periode; /**< Periode du mouvement sinusoidal */
  Vecteur m_Sin_vitRef; /**< Vecteur de reference du mouvement sinusoidal */  
  Scalar m_freqX, m_freqY, m_freqZ, m_phase;
  Scalar m_ampX, m_ampY, m_ampZ;
  static Vecteur m_prev;
  //@}
};

#endif
