#ifndef _Cellule
#define _Cellule

#include <list>
using namespace std;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Localisation geographique dans le Linkedcell
    @author GRAINS Project - IFP - 2009 */
enum MPIGeoLocalisation
  {
    MPIGEO_NORTH		= 0,
    MPIGEO_NORTH_EAST		= 1,
    MPIGEO_NORTH_WEST		= 2,
    MPIGEO_NORTH_TOP		= 3,
    MPIGEO_NORTH_BOTTOM		= 4,
    MPIGEO_NORTH_EAST_TOP	= 5,
    MPIGEO_NORTH_EAST_BOTTOM	= 6,
    MPIGEO_NORTH_WEST_TOP	= 7,
    MPIGEO_NORTH_WEST_BOTTOM	= 8,
    MPIGEO_SOUTH		= 9,
    MPIGEO_SOUTH_EAST		= 10,
    MPIGEO_SOUTH_WEST		= 11,
    MPIGEO_SOUTH_TOP		= 12,
    MPIGEO_SOUTH_BOTTOM		= 13,
    MPIGEO_SOUTH_EAST_TOP	= 14,
    MPIGEO_SOUTH_EAST_BOTTOM	= 15,
    MPIGEO_SOUTH_WEST_TOP	= 16,
    MPIGEO_SOUTH_WEST_BOTTOM	= 17,
    MPIGEO_EAST			= 18,
    MPIGEO_WEST			= 19,
    MPIGEO_EAST_TOP		= 20,
    MPIGEO_EAST_BOTTOM		= 21,
    MPIGEO_WEST_TOP		= 22,
    MPIGEO_WEST_BOTTOM		= 23,
    MPIGEO_TOP			= 24,
    MPIGEO_BOTTOM		= 25,
    MPIGEO_NONE			= 26
  };


#include "Particule.H"
#include "MonObstacle.H"


/** @brief Cellule de decomposition du recipient (par exemple un reacteur).

    Elle contient 0-N particules, et est utilisee lors du deplacement
    des particules pour determiner les contacts : particule-particule
    et particule-obstacle. Sa forme est obligatoirement cubique.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation */
// ============================================================================
class Cellule
{
public:
  /**@name Constructeurs */
  //@{
  /** @brief Constructeur par defaut */
  Cellule();

  /** @brief Constructeur avec initialisation
  @param id1 Identificateur
  @param x  Position direction X
  @param y  Position direction Y
  @param z  Position direction Z
  @param arete_X Longueur de l'arete des cellules dans la direction X
  @param arete_Y Longueur de l'arete des cellules dans la direction Y
  @param arete_Z Longueur de l'arete des cellules dans la direction Z
  @param xmax_ borne sup�rieur du LinkedCell dans la direction X
  @param ymax_ borne sup�rieur du LinkedCell dans la direction Y
  @param zmax_ borne sup�rieur du LinkedCell dans la direction Z */
  Cellule( int id1, int x, int y, int z,
  	double arete_X, double arete_Y, double arete_Z,
	double xmax_, double ymax_, double zmax_ );

  /** @brief Constructeur avec initialisation
  @param id1 Identificateur
  @param x  Position direction X
  @param y  Position direction Y
  @param z  Position direction Z
  @param arete_X Longueur de l'arete des cellules dans la direction X
  @param arete_Y Longueur de l'arete des cellules dans la direction Y
  @param arete_Z Longueur de l'arete des cellules dans la direction Z
  @param tag_ tag de la cellule
  @param OL origine locale de la grille a laquelle la cellule appartient
  @param xmax_ borne sup�rieur du LinkedCell dans la direction X
  @param ymax_ borne sup�rieur du LinkedCell dans la direction Y
  @param zmax_ borne sup�rieur du LinkedCell dans la direction Z
  @param geoloc_ localisation geographique de la cellule dans le LinkedCell */
  Cellule( int id1, int x, int y, int z, double arete_X, double arete_Y,
	double arete_Z, int tag_, const Point &OL,
	double xmax_, double ymax_, double zmax_,
	MPIGeoLocalisation geoloc_ = MPIGEO_NONE );

  /** @brief Destructeur. */
  ~Cellule();
  //@}


  /**@name Methods Get */
  //@{
  /** @brief Renvoi la liste des cellules du voisinage complet */
  const list<Cellule*>* getVoisinageComplet() const;

  /** @brief Renvoie le nombre d'obstacles dans le voisinage de la cellule */
  int nombreObstacles() const;

  /** @brief Volume global des particules dans la cellule.
  Si une particule appartient a la cellule, son volume est dans la cellule.
  On ne regarde par les intersection entre cellules et particules.
  @return Le volume des particules */
  double getVolumeParticules();

  /** @brief Position du centre de gravite de la cellule.
  @return Le point demande */
  Point const* getCentre() const;

  /** @brief Liste des particules dans la cellule */
  list<Particule*>* getParticules();

  /** @brief Liste des obstacles dans la cellule */
  list<MonObstacle*>* getObstacle();

  /** @brief Localisation geographique */
  MPIGeoLocalisation getGeoLocalisation() const;
  //@}


  /**@name Methods */
  //@{
  /** @brief Ajout de la particule a la cellule.
  Une verification de non appartenance est realisee.
  @param particule La particule a ajouter */
  void add( Particule* particule );

  /** @brief Ajout de la cellule voisine au voisinage pour les contacts
  @param voisine cellule voisine */
  void addVoisineContact( Cellule *voisine );

  /** @brief Ajout de la cellule voisine au voisinage complet
  @param voisine cellule voisine */
  void addVoisine( Cellule *voisine );

  /** @brief Ajout d'un obstacle dans le voisinage de la cellule
  @param obstacle_ l'obstacle */
  void addObstacle( MonObstacle* obstacle_ );

  /** @brief Supression des particules dans la cellule. */
  void clearParticules();

  /** @brief La cellule contient-elle la Particule indique ?
  @return TRUE si le monolithe est geree par la cellule.
  @param monolithe particule */
  bool Contient( Particule* monolithe );

  /** @brief La particule est elle en contact avec un composant present dans
  la cellule. Une verification est effectuee pour eviter une recherche de
  contact de la particule sur elle-meme.
  @return Un contact eventuel.
  @param particule La particule etudiee. */
  bool isContact( const Particule* particule ) const;

  /** @brief La particule est elle en contact avec un composant present dans
  la cellule, variante VdW. Une verification est effectuee pour eviter une
  recherche de contact de la particule sur elle-meme.
  @return Un contact eventuel.m_centre
  @param particule La particule etudiee. */
  bool isContactVdW( const Particule* particule ) const;

  /** @brief La particule est elle en contact avec un obstacle periodique
  @return Un contact eventuel.
  @param particule La particule etudiee. */
  bool isContactVdW_ObstaclePeriodique( const Particule* particule ) const;

  /** @brief La particule est elle proche d'un composant present dans
  la cellule. Une verification est effectuee pour eviter une
  recherche de contact de la particule sur elle-meme.
  @return Un contact eventuel.
  @param particule La particule etudiee. */
  bool isProche( const Particule* particule ) const;

  /** @brief La particule est elle proche d'un composant present dans
  la cellule, variante VdW. Une verification est effectuee pour eviter une
  recherche de contact de la particule sur elle-meme.
  @return Un contact eventuel.
  @param particule La particule etudiee. */
  bool isProcheVdW( const Particule* particule ) const;

  /** @brief La cellule est elle vide (pas de Particule)
  @return VRAI si la cellule est vide */
  bool isEmpty() const;

  /** @brief Mise a jour des particules dans la cellule.
  Les particules hors cellule sont renvoyees pour traitement.
  @param particulesExit Les Particules hors cellule */
  void LinkUpdate( list<Particule*> &particulesExit );

  /** @brief Existe t'il une proximite entre le composant et les composants
  de la cellule.
  A noter, le composant ne se regarde pas lui-meme.
  @return TUE si un risque est detecte
  @param composant Le composant etudie. */
  bool PbProximite( const Composant* composant ) const;

  /** @brief Suppression d'une particule de la cellule.
  @param particule La particule a supprimer. */
  void remove( Particule *particule );

  /** @brief Suppression d'un obstacle de la cellule.
  @param obs L'obstacle a supprimer. */
  void remove( MonObstacle *obs );

  /** @brief Operateur <<
  @param f flux de sortie
  @param C cellule */
  friend ostream& operator << ( ostream &f, const Cellule &C );
  //@}


  /** @name Methods Static */
  //@{
  /** @brief Indices de la cellule associee a la position
  @param position La position demande
  @param id Les indices de la cellule */
  static void GetCellule( const Point &position, int* id );

  /** @brief Definition de l'espace des cellules
  @param nbX Nombre de cellules en X
  @param nbY Nombre de cellules en Y
  @param nbZ Nombre de cellules en Z */
  static void SetEspace( int nbX, int nbY, int nbZ );

  /** @brief Renvoie la geolocalication dans une chaine de caractere
  @param geoloc_ geolocalisation */
  static string getMPIGeoLocalisationName( MPIGeoLocalisation geoloc_ );

  /** @brief Renvoie la geolocalication dans une chaine de caractere
  @param geoloc_ geolocalisation */
  static string getMPIGeoLocalisationName( int geoloc_ );

  /** @brief Renvoie la geolocalication dans une chaine de caractere
  @param geoloc_ geolocalisation */
  static string getMPIGeoLocalisationName_generic( int geoloc_ );
  //@}


  /**@name Operators */
  //@{
  /** @brief  Acces a la cellule par son indice global
  @param xyz indice de la cellule */
  int operator [] ( int xyz ) const;

  /** @brief Comparaison entre deux cellules.
  @return Vrai si identite des deux cellules.
  @param cellule La cellule a comparer. */
  bool operator == ( const Cellule &cellule ) const;
  //@}


  /**@name Class Friend */
  //@{
  friend class LinkedCell; /**< Une cellule est intimement liee a l'algorithme
  	de LinkedCell. */
  //@}


private:
  /**@name Parameters */
  //@{
  int m_number; /**< Numero de la cellule */
  int m_tag;   /**< Tag de la cellule:
  <ul>
    <li> 0=interieur,
    <li> 1=dans la zone de recouvrement d'un autre processeur et
    geometriquement sur ce processeur,
    <li> 2=dans la zone de recouvrement de ce processeur et geometriquement sur
    un autre processeur
  </ul> */
  Point m_centre; /**< coordonnees du centre de la cellule */
  MPIGeoLocalisation m_GeoLocCell; /**< localisation geographique de la cellule
  	dans le LinkedCell, utile pour les cellules taggees 1 uniquement */
  int m_cel[3]; /**< Position de la cellule en x-y-z.
      Cette position est definie par rapport a un caisson exterieur complet,
      meme si la structure reelle est de forme complexe. */
  list<Cellule*> m_voisinesContact; /**< Cellules du voisinage pour les
  	contacts; On utilise : 3 au-dessus, 1 a droite, 9 arriere. */
  list<Cellule*> m_allVoisines; /**< Toutes les cellules du voisinage (cas
  	general 3D: 26 voisines)  */
  list<MonObstacle*> m_obstacles; /**< obstacles dans le voisinage de la
  	cellule */
  list<Particule*> m_particules; /**< Particules contenues dans la cellule */
  static int nbi; /**< Dimension de l'espace des Cellules en X */
  static int nbj; /**< Dimension de l'espace des Cellules en Y */
  static int nbk; /**< Dimension de l'espace des Cellules en Z */
  static double m_arete_X; /**< Longueur de l'arete des Cellules en X */
  static double m_arete_Y; /**< Longueur de l'arete des Cellules en Y */
  static double m_arete_Z; /**< Longueur de l'arete des Cellules en Z */
  static Point m_LC_origine_locale; /**< Origine locale de la grille */
  static double m_xmax; /**< borne sup�rieur du LinkedCell dans la
  	direction X */
  static double m_ymax; /**< borne sup�rieur du LinkedCell dans la
  	direction Y */
  static double m_zmax; /**< borne sup�rieur du LinkedCell dans la
  	direction Z */
  //@}


  /**@name Methods Set */
  //@{
  /** @brief Calcule le centre de la cellule */
  void setCentre();
  //@}
};

#endif
