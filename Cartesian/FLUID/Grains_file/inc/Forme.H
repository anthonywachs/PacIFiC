// F.PRADEL - Janv.2000 - Creation
// G.FERRER - Mars.2000 - Le Contact ne renvoie plus l'objet mais le pointeur
// ============================================================================
#ifndef _Forme
#define _Forme

#include "Erreur.H"
#include "Basic.H"
#include "Transform.H"
#include "Convex.H"
#include "Vecteur.H"
#include "Point.H"
using namespace solid;

class Cellule;
class Torseur;


class Forme;
ostream &operator << ( ostream &fileOut, const Forme &forme );
istream &operator >> ( istream &fileIn, Forme &forme );


/** @brief entite abstraite de geometrie.

    Une forme correspond a une abstraction de geometrie.
    Cette geometrie correspond a un Convexe & sa Position dans l'espace.
    Quelques informations additionnelles sont incluses : rayon externe...

    @author F.PRADEL - Institut Francais du Petrole - 1999 - Creation
    @author G.FERRER - Institut Francais du Petrole - 2000 - Polymorphisme
    @author D.PETIT  - Institut Francais du Petrole - 2000 - Modification
    @author A.WACHS  - Institut Francais du Petrole - 2012 - Modification */
// ============================================================================
class Forme
{
public:
  /**@name Constructeurs */
  //@{
  /** @brief Constructeur par defaut. */
  Forme();

  /** @brief Constructeur avec initialisation, utilise exclusivement par les
  compObstacle dont la forme propre n'a pas de sens
  @param convex_   Le convexe de reference
  @param position_ La transformation de reference. */
  Forme( Convex *convex_, const Transform &position_ );

  /** @brief Constructeur par copie.
  @param forme La forme de reference */
  Forme( const Forme &forme );

  /** @brief Destructeur. */
  virtual ~Forme() ;
  //@}


  /**@name Methods */
  //@{
  /** @brief BBox (boite contenant la forme) d'une forme.
  @return La boite englobant la forme. */
  BBox BoxForme() const;

  /** @brief Determine l'inertie et l'inertie inverse de la forme.
  @return TRUE si l'inertie est calcule.
  @param inertie L'inertie de la forme.
  @param inertie_1 L'inertie inverse de la forme.*/
  bool BuildInertie( Scalar *inertie, Scalar *inertie_1 ) const;

  /** @brief Calcul de la distance entre deux formes.
  @return la distance entre les deux formes.
  @param voisine  La forme voisine. */
  Scalar DistanteDe( const Forme &voisine ) const;

  /** @brief Evaluation d'un contact entre deux formes: veritable contact
  geometrique
  @return TRUE si un contact est present.
  @param voisine La deuxieme forme a etudier. */
  bool isContact( const Forme &voisine ) const;

  /** @brief Evaluation de proximite entre formes.
  La proximite est liee a l'intersection des BBox des deux formes.
  @return TRUE si les deux formes sont a proximitee
  @param voisine La deuxieme forme etudiee. */
  bool isProche( const Forme &voisine ) const;

  /** @brief Rotation de la forme
  @param q Quaternion representant la rotation */
  void Rotate( const Quaternion& q );

  /** @brief Applique la transformation trot � la forme
  @param trot transformation � appliquer */
  void composeTransform( const Transform &trot );

  /** @brief Compose la transformation de la forme par trot:
  !!! Resultat : d'abord la transformation de la forme puis trot !!!
  @param trot transformation � appliquer */
  void composeTransformRight( const Transform &trot );

  /** @brief Compose la transformation de la forme par une rotation trot par
  rapport au centre de gravite de la forme
  !!! Resultat : d'abord la transformation de la forme puis trot !!!
  !!! Cette composition ne change pas l'origine de la forme !!!
  @param trot transformation � appliquer */
  void composeRotationRight( const Transform &trot );

  /** @brief Copie la Transformation du Convex decrivant la Forme
  @param vit tableau de copie
  @param i position de depart dans le tableau de copie */
  void copyTransform( double *vit, int i ) const;

  /** @brief Copie la Transformation du Convex decrivant la Forme translat�e
  d'un vecteur (utile pour les particules periodiques en parallele)
  @param vit tableau de copie
  @param i position de depart dans le tableau de copie
  @param vec vecteur de translation */
  void copyTransform( double *vit, int i, Vecteur const& vec ) const;
  //@}


  /**@name Methods Set */
  //@{
  /** @brief Positionne l'origine de la forme.
  @param pos Definition du point que l'on affecte a l'origine de la forme,
  valeurs en X-Y-Z. */
  void setOrigin( const Scalar *pos );

  /** @brief Positionne l'origine de la forme.
  @param gx coordonn�e x de la position
  @param gy coordonn�e y de la position
  @param gz coordonn�e z de la position */
  void setOrigin( const double &gx, const double &gy, const double &gz );

  /** @brief Positionne l'origine de la forme.
  @param pos Definition du point que l'on affecte a l'origine de la forme,
  valeurs en (X,Y,Z) */
  void setOrigin( const Point &pos );

  /** @brief Positionne une forme dans l'espace.
  @param pos Description de la matrice de base
  (12 coefficients : 3+1+3+1+3+1) puis la position du point */
  void setPosition( const Scalar *pos );

  /** @brief Affectation de la valeur du rayon
  @param r le rayon a affecter*/
  void setRayon( Scalar r );

  /** @brief gets the shrinking choice */
  int getShrinkingMode();

  /** @brief Affectation d'une Transformation imposee */
  void setTransform( const Transform &transform_ );
  //@}


  /**@name Methods Get */
  //@{
  /** @brief Renvoie le point centre de la forme.
  @return centre de la forme */
  Point const* getCentre() const;

  /** @brief Valeurs des coordonnees du centre
  @param pos Les valeurs du centre de la forme. */
  void getCentre( Scalar *pos ) const;

  /** @brief Acces au Convex decrivant la Forme. */
  const Convex* getConvex() const;

  /** @brief Acces a la Transformation du Convex decrivant la Forme. */
  const Transform* getTransform() const;

  /** @brief Acces a la Transformation du Convex decrivant la Forme. */
  Transform* getTransform();

  /** @brief Renvoie la matrice donnant l'orientation de la forme.
  @return La matrice d'orientation de la forme*/
  Matrix getOrientation() const;

  /** @brief Valeur du rayon externe de la forme.
  @return La valeur demandee. */
  Scalar getRayon() const;

  /** @brief sets the shrinking radius
  for shrinking particles */
  void set_shrinking_radius(Scalar CurrentRadius);




  /** @brief Volume de la forme.
  @return Le volume calcule. */
  Scalar getVolume() const;
  //@}


  /**@name Operators */
  //@{
  /** @brief Operateur de creation du vecteur F-A reliant le point A
  au centre de la forme F.
  @return Vecteur joignant le point A au centre de la forme
  @param rhs Le point A */
  Vecteur operator-( const Point& rhs ) const;

  /** @brief Egalite sur les adresses des deux formes.
  @return TRUE si les formes sont identiques (meme adresse)
  @param particule la forme a tester */
  bool operator==( const Forme &particule ) const;

  /** @brief Translate la forme d'un vecteur donne.
  @return La forme actuelle a sa nouvelle position
  @param rhs Vecteur de translation.*/
  Forme& operator+=( const Vecteur& rhs );
  //@}


  /**@name Methods I/O */
  //@{
  /** @brief Lecture de l'information de position.
  @param fileIn Flux emetteur. */
  void readPosition( istream &fileIn );

  /** @brief Lecture de l'information de position au format de reload 2014
  @param fileIn Flux emetteur. */
  void readPosition2014( istream &fileIn );

  /** @brief Lecture de l'information de position au format de reload 2014 en
  binaire
  @param fileIn Flux emetteur. */
  void readPosition2014_binary( istream &fileIn );

  /** @brief Lecture de l'information statique.
  @param statique Flux emetteur. */
  void readStatique( istream &statique );

  /** @brief Ecriture de l'information de position.
  @param position Flux recepteur. */
  void writePosition( ostream &position ) const;

  /** @brief Ecriture de l'information de position pour le Fluide */
  void writePositionInFluid( ostream &fluid );

  /** @brief Ecriture de l'information de position pour Basilisk (fluide) */
  void writePositionInFluid (struct BasiliskDataStructure * b);


  /** @brief Ecriture de l'information statique.
  @param statique Flux recepteur
  @param composant solid object */
  virtual void writeStatique( ostream &statique, Composant const*
      composant = NULL );

  /** @brief Ecriture de l'objet dans sa configuration courante pour
  post-processing avec GMV
  @param fileOut Flux de sortie. */
  void GMVoutput( ostream &fileOut ) const;

  /** @brief Ecrit les points du convexe pour post-processing avec Paraview
  @param f Flux de sortie
  @param translation translation du centre de gravite */
  void write_polygonsPts_PARAVIEW( ostream &f,
  	Vecteur const* translation = NULL ) const;

  /** @brief Ecrit le convexe au format STL pour lien avec openFoam
  @param f Flux de sortie */
  void write_convex_STL( ostream &f ) const;

  /** @brief Renvoie les points du convexe pour post-processing avec Paraview
  @param translation translation du centre de gravite */
  list<Point> get_polygonsPts_PARAVIEW( Vecteur const* translation = NULL )
  	const;
  //@}


  /**@name Methods Friends */
  //@{
  /** @brief Intersection entre deux Forme */
  friend bool intersect( const Forme &a, const Forme &b );

  /** @brief Operateur d'ecriture.
  @return Flux recepteur
  @param fileOut Flux recepteur
  @param forme Objet courant. */
  friend ostream &operator << ( ostream &fileOut,
	const Forme &forme );

  /** @brief Operateur de lecture.
  @return Flux emetteur
  @param fileIn  Flux emetteur
  @param forme Objet courant */
  friend istream &operator >> ( istream &fileIn,
	Forme &forme );
  //@}


  /**@name Class Friend */
  //@{
  /** @brief Autorisation d'acces a la classe de gestion des ecritures. */
  friend class ofstream_Grains;

  /** @brief Autorisation d'acces a la classe de gestion des lectures. */
  friend class ifstream_Grains;
  //@}


protected:
  /**@name Parameters */
  //@{
  Transform m_position; /**< Position et Orientation de l'objet. */
  Convex *m_convex; /**< Forme de Convexe. */
  Scalar m_rayon; /**< Rayon externe de la forme */
  Scalar m_volume; /**< Volume de la forme */
  Scalar m_rayonSphereEquivalente; /**< Rayon de la sphere de meme volume */
  int Shrinking; /**< The shriking choice */
  //@}
};


/** @brief Intersection entre les deux Forme indiquees ? */
bool intersect( const Forme &a, const Forme &b );

#endif
