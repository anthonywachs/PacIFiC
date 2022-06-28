#ifndef _CompParticule
#define _CompParticule

#include "Particule.H"
#include "Forme.H"
#include "ElementParticule.H"
#include "PointContact.hh"


/** @brief Gestion des particules composites : COMPosed PARTICULE
    une particule composite est une particule composee de particules
    convexes elementaires.

    @author D. RAKOTONIRINA - IFP Energies nouvelles - 2014 - Creation */
// ============================================================================
//
class ElementParticule;

class CompParticule : virtual public Particule
{
public:
  /** @name Constructeur */
  //@{
  /** @brief Constructeur
  @param root Noeud "<CompParticule>"
  @param autonumbering numerotation automatique ou non
  @param pc classe de particules
  @param name Nom eventuel de la particule composite */
  CompParticule( DOMNode* root, const bool &autonumbering = true,
  	const int &pc = 0, const string &name = "NonConvex" );

  CompParticule( const bool &autonumbering = true );


  /** @brief Constructeur avec arguments
  @param id_ numero
  @param ParticuleRef particule de reference
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
  @param coordination_number_ nombre de contacts de la particule
  @param updatePosition do we update position or not ? */
  CompParticule( const int &id_, Particule const* ParticuleRef,
	const double &vx, const double &vy, const double &vz,
	const double &qrotationx, const double &qrotationy,
	const double &qrotationz, const double &qrotations,
	const double &rx, const double &ry, const double &rz,
	const Scalar m[16],
	const ParticuleActivity &activ,
	const int &tag_,
	const int &coordination_number_ = 0,
	const bool &updatePosition = false );

  /** @brief Constructeur avec arguments
  @param id_ numero
  @param ParticuleRef particule de reference
  @param vtrans vitesse de translation
  @param vrot vitesse de rotation
  @param qrot quaternion de rotation
  @param config transformation li�e a la particule
  @param activ activite
  @param tag_ tag */
  CompParticule( const int &id_, Particule const* ParticuleRef,
	const Vecteur &vtrans,
	const Quaternion &qrot,
	const Vecteur &vrot,
	const Transform &config,
	const ParticuleActivity &activ,
	const int &tag_ );

  /** @brief Constructeur par Copie.
  Le torseur des forces appliquees est vide.
  Les types de Forme et de Cinematique sont instancies.
  @param copie Particule de reference. */
  CompParticule( const CompParticule &copie );

  /** @brief Destructeur */
  virtual ~CompParticule();
  //@}

  /**@name Methods */
  //@{
  /** @brief Contact entre particule composite et une particule ou un obstacle
  Si le contact est realise, on ajoute le contact au deux composants.
  @exception ErreurContact Si un choc entre croute est detecte
  @param voisin Le composant (monolithe) a etudier
  @param dt Increment de temps
  @param temps temps physique
  @param LC grille de cellules
  @param listContact list of information about contacts */
  virtual void SearchContact( Composant* voisin, double dt,
      double const& temps, LinkedCell *LC,
      list<ContactInfos*> &listContact );

  /** @brief Ecrit les points du convexe pour post-processing avec Paraview
  @param f Flux de sortie
  @param transform transformation courante
  @param translation translation du centre de gravite */
  virtual void write_polygonsPts_PARAVIEW( ostream &f,
  	const Transform &transform,
  	Vecteur const* translation = NULL ) const;

  /** @brief Nombre de points pour post-processing avec Paraview */
  int numberOfPoints_PARAVIEW() const;

  /** @brief Nombre de polygones elementaires pour post-processing avec
  Paraview */
  int numberOfCells_PARAVIEW() const;

  /** @brief Renvoie les points du convexe pour post-processing avec Paraview
  @param translation translation du centre de gravite */
  list<Point> get_polygonsPts_PARAVIEW(
  	Vecteur const* translation = NULL ) const;

  /** @brief Ecrit les points du convexe pour post-processing avec Paraview
  @param f Flux de sortie
  @param translation translation du centre de gravite */
  void write_polygonsPts_PARAVIEW( ostream &f,
  	Vecteur const* translation = NULL ) const;

  /** @brief Applique la transformation trot au composant
  @param trot transformation � appliquer */
  void composePosition( const Transform &trot );

  /** @brief Constructeur d'un clone par copie */
  Particule* createCloneCopy() const ;

  /** @brief Calcul du poids de la particule composite */
  void computeWeight();

  /** Le composant est il une particule composite ?
  @return Vrai si la particule est un composite */
  bool isCompParticule() const {return ( true ); }

  /** @brief Contact entre la particule et une particule.
  Si le contact est realise, on ajoute le contact au deux composants.
  @exception ErreurContact Si un choc entre croute est detecte
  @param voisin Le composant (monolithe) a etudier
  @param dt Increment de temps
  @param temps temps physique
  @param LC grille de cellules */
  virtual void InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC ) throw (ErreurContact);

  /** @brief Contact entre l'obstacle et une particule pour le post processing
  Si le contact est realise, on ajoute le contact au deux composants.
  @param voisin Le composant (monolithe) a etudier
  @param listOfContacts liste de contacts */
  virtual void InterActionPostProcessing( Composant* voisin,  Scalar dt,
	list<struct PointForcePostProcessing>* listOfContacts )
	throw (ErreurContact);

  /** @brief Les deux composants sont ils proches ? La proximite est liee a
  l'intersection des BBox des deux formes. Si les deux composants sont
  identiques (==), le contact n'existe pas.
  Limitation: le composant voisin ne doit PAS etre de type deriv� CompObstacle.
  @return TRUE Si un contact existe
  @param voisin Le deuxieme composant etudie. */
  virtual bool isProche( const Composant* voisin ) const;

  /** @brief Les deux composants sont ils en contact ? Variante avec les rayons
  de VdW. La proximite est liee a l'intersection des BBox des deux formes.
  Si les deux composants sont identiques (==), le contact n'existe pas.
  Limitation: le composant voisin ne doit PAS etre de type deriv� CompObstacle.
  @return TRUE Si un contact existe
  @param voisin Le deuxieme composant etudie. */
  virtual bool isProcheVdW( const Composant* voisin ) const;

  /** @brief Les deux composants sont ils en contact ?
  Si les deux composants sont identiques (==), le contact n'existe pas.
  Limitation: le composant voisin ne doit PAS etre de type deriv� CompObstacle.
  @return TRUE Si un contact existe
  @param voisin Le deuxieme composant etudie. */
  virtual bool isContact( const Composant* voisin ) const;

  /** @brief Les deux composants sont ils en contact ? Variante avec les rayons
  de VdW. Si les deux composants sont identiques (==), le contact n'existe pas.
  Limitation: le composant voisin ne doit PAS etre de type deriv� CompObstacle.
  @return TRUE Si un contact existe
  @param voisin Le deuxieme composant etudie. */
  virtual bool isContactVdW( const Composant* voisin ) const;

  /** @brief Inertie de la particule.
  @param inertie L'inertie de la particule composite
  @param inertie_1 L'inertie inverse de la particule composite */
  bool BuildInertie(Scalar *inertie, Scalar *inertie_1) const;

  /** @brief Resolution des equations de la dynamique et double integration
  pour obtenir la nouvelle vitesse et position.
  <ul>
    <li> v^{n+0.5} & =  &v_i^{n-0.5} + Dt * f_{total}(i)^{n}
    <li> x^{n+1}_i & = & x^{n}_i + Dt * v^{n+0.5}_i
  </ul>
  De plus, cette fonction ne vide pas la liste des torseurs ;
  cette action doit etre realisee par la methode "clearForces()".
  @exception ErreurDeplacement si le deplacement est trop grand.
  @param temps Temps de simulation
  @param dt Pas de temps de resolution */
  void Deplacer( Scalar temps, double dt ) throw(ErreurDeplacement);

  /** @brief Initialize a faux le boolean correspondant au calcul de la
  transformation avec scaling par l'epaisseur de croute */
  void initializeVdWtransform_to_notComputed();

  /** @brief Determine s'il y a intersection entre une Forme et un Composant
  @return Vrai s'il y a intersection */
  bool intersectWithForme( const Forme &b ) ;

  /** @brief Copie du quaternion de rotation dans le vecteur vit en d�butant �
  la position i
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  void copyQuaternionRotation( double *vit, int i ) const;

  /** @brief Copie du vitesse de translation dans le vecteur vit en d�butant �
  la position i
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  void copyVitesseTranslation( double *vit, int i ) const;

  /** @brief Copie du vitesse de rotation dans le vecteur vit en d�butant � la
  position i
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  void copyVitesseRotation( double *vit, int i ) const;

  /** @brief Copie de la cinematique au temps t-2dt: vitesse de translation,
  vitese de rotation, variation de QDM translationnelle, variation de QDM
  rotationnalle
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  void copyCinematiqueNm2(double *vit,int i) const;

  /** @brief Copie de la transformation dans le vecteur vit en d�butant �
  la position i
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  void copyTransform( double *vit, int i ) const;

  /** @brief Copie de la transformation dans le vecteur vit en d�butant �
  la position i, avec une translation du centre de gravite (utile pour les
  particules periodiques en parallele)
  @param vit vecteur de copie
  @param i position dans le vecteur vit
  @param vec vecteur de translation */
  void copyTransform( double *vit, int i, Vecteur const& vec ) const;
  //@}

  /**@name Accessors */
  //@{
  /** @brief Renvoie le nombre total de particules elementaires */
  size_t getNbreElemPart() const;

  /** @brief Renvoie le "volume" de la particule composite
  @return volume*/
  Scalar getVolume() const;

  /** @brief Renvoie le nom associe a la particule composite */
  virtual string getPartName() const;

  /** @brief Renvoie un vecteur de particules elementaires */
  vector<ElementParticule*> getElementParticules() const;

  /** @brief Renvoie les positions initiales des particules elementaires
  a l'insertion */
  vector<Vecteur> getInitialRelativePositions() const;

  /** @brief Renvoie les positions des particules elementaires
  apres mise a jour */
  vector<Vecteur> getRelativePositions() const;

  /** @brief Renvoie les matrices de rotations intiales des particules
  elementaires */
  vector<Matrix> getInitialMatrix() const;

  /** @brief Renvoie le nombre de sommets ou un code �quivalent */
  virtual int getNbCorners() const;
  //@}

  /**@name Methods set */
  //@{

  /** @brief Affectation d'un champ de vitesse de translation.
  @param translation vecteur de translation */
  void setVitesseTranslation( const Vecteur &translation );

  /** @brief Affectation d'un champ de vitesse de rotation.
  @param rotation vecteur de rotation */
  void setVitesseRotation( const Vecteur &rotation );

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

  //@}

  /**@name Methods I/O */
  //@{

  /** @brief Ecrit le convexe pour post-processing avec Paraview
  @param connectivity connectivite des polygones Paraview
  @param offsets decalage dans la connectivite
  @param cellstype type de polygones Paraview
  @param firstpoint_globalnumber numero global du 1er point
  @param last_offset dernier offset utilise pour le convexe precedent */
  void write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const;

  /** @brief Sauvegarde du monolithe pour Reload
  @param fileSave Flux de sauvegarde
  @param composant solid body */
  void write( ostream &fileSave, Composant const* composant = NULL ) const;

  /** @brief Lecture de la particule pour Reload. Utilise pour les particules de
  classe de reference dans l'entete du fichier de reload et pour l'ancien format
  de reload
  @param fileSave Flux de lecture
  @param ParticuleClassesReference particules de r�f�rence
  	pour chaque classe de particules */
  void read( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference = NULL );

  /** @brief Lecture de la particule pour Reload. Utilise pour les particules
  dans la simulation (actives ou en attentes) pour le format de reload 2014
  @param fileSave Flux de lecture
  @param ParticuleClassesReference particules de r�f�rence
  	pour chaque classe de particules */
  void read2014( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference );

  /** @brief Lecture de la particule pour Reload. Utilise pour les particules
  dans la simulation (actives ou en attentes) pour le format de reload 2014
  en binaire
  @param fileSave Flux de lecture
  @param ParticuleClassesReference particules de r�f�rence
  	pour chaque classe de particules */
  void read2014_binary( istream &fileSave, vector<Particule*> const*
  	ParticuleClassesReference );

  /** @brief Ecriture de l'information de position pour le Fluide */
  virtual void writePositionInFluid( ostream &fluid );
  //@}


protected:

  /** @name Methods */
  //@{
  /** @brief Determine le rayon circonscrit du composite */
  void setRayon();

  /** @brief Positionne la particule composite dans l'espace */
  void setPosition( const Point &centre);

  /** @brief Positionne un composant dans l'espace.
  @param pos Description de la matrice de base
  (12 coefficients : 3+1+3+1+3+1) puis la position du point,
  soit un vecteur de 16 doubles */
  void setPosition( const Scalar *pos );

  /** @brief Positionne les particules elementaires dans l'espace */
  void setElementPosition();
  //@}

  /** @name Accessors */
  //@{
  /** @brief Renvoie le rayon circonscrit du composite */
  Scalar getRayon() const;

  /** @brief Renvoie un cube pour le calcul des elements d'inertie
  v[0] = x_max, v[1] = x_min
  v[2] = y_max, v[3] = y_min
  v[4] = z_max, v[5] = z_min */
  vector<Scalar> getCompFeaturesBox() ;
  //@}

  /** @name Parameters */
  //@{
  string m_name; /**< Nom de la particule composite */
  //@}

private:

  /** @name Parameters */
  //@{
  size_t m_nbreElemPart; /**< Nombre de particules elementaires */

  vector<ElementParticule*> m_elementaryParticules; /**< Vecteur comportant l'ensemble
  	des constituants elementaires du composite */

  vector<Vecteur> m_InitialRelativePositions; /**< Vecteur des positions relatives
  apres lecture du fichier insertion */

  vector<Vecteur> m_RelativePositions; /**< Vecteur des positions relatives par
  rapport � 0 apres positionnement en (0.,0.,0.) */

  vector<Matrix> m_InitialMatrix; /**< Vecteur de Matrices de rotation initiale
  de chaque particule elementaire */

  Scalar m_compVolume; /**< Volume approxime de la particule composite */

  Point m_CG; /**< Position du centre de gravite apres approximation */

  vector<Scalar> m_compInertie; /**< Les elements d'inerties A,B,C,D,E,F */
  //@}

  /** @name Methods */
  //@{
  /** @brief Gestion de la lecture des particules elementaires
      @param fileSave Flux de sauvegarde
      @param nbre Nombre de particules elementaires
      @param cinematique Cinematique de la particule composite */
  void  ReadElementParticules( istream &fileSave, size_t nbre,
      CineParticule* cinematique );


  /** @brief Return the 2D intersection point of a circle and a segment
  defined by 2 points A and B
  @param ptA first point inside the circle
  @param ptB second point
  @param ptC center of the circle
  @param ptI intersection point (ouput)
  @param radius of the circle */
  void intersect2d( Point &ptA, Point &ptB, Point &ptC, Point &ptI,
      double const &radius );


  void writePosition( ostream &position, const Vecteur &relPos, const
    Matrix &initRot, ElementParticule* particule ) const;

  struct pairPts
  {
    Point pt1;
    Point pt2;
  };

  struct tripletPts
  {
    Point Center;
    Point Poly1;
    Point Poly2;
  };

  //@}

  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut interdit */
  CompParticule();
  //@}
};

#endif
