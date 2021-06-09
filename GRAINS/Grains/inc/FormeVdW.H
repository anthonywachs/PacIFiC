#ifndef _FormeVdW__
#define _FormeVdW__

#include "Forme.H"
#include "ReaderXML.hh"
#include "WriterXML.hh"
#include "Erreur.H"

class PointContact;


/** @brief Entite de geometrie pour laquelle la recherche de contact utilise 
    un concept de rayon d'interaction similaire a un concept de force de 
    "Van der Waals".

    Dans ce type de forme, le contact n'est pas approche directement mais 
    avec un rayon d'interaction lie a chaque forme. Ceci permet de lisser 
    le risque de saut pour le point de contact pour des convexes de forte
    angularite.

    @author G. FERRER - Institut Francais du Petrole - 2002 - Creation 
    @author D. RAKOTONIRINA - IFP Energies Nouvelles - Sept. 2014
    - Modification*/
// ============================================================================
class FormeVdW : public Forme
{
public:
  /**@name Constructeurs */
  //@{
  /** @brief Constructeur par defaut */
  FormeVdW();

  /** @brief Constructeur par copie. 
  @param forme La forme de reference */
  FormeVdW( const FormeVdW &forme );

  /** @brief Constructeur avec initialisation, utilise exclusivement par les
  compObstacle dont la forme propre n'a pas de sens
  @param convex_   Le convexe de reference
  @param position_ La transformation de reference. */
  FormeVdW( Convex *convex_, const Transform &position_ );

  /** @brief Contructeur avec lecture sur fichier de persistance
  @param fileIn Flux de lecture.
  @param type type de particule elementaire dans le cas d'une particule
  composite */
  FormeVdW( istream &fileIn, string type = ""  );

  /** @brief Constructeur avec decodage du noeud XML
  @param root Le noeud "Obstacle" */
  FormeVdW( DOMNode *root );

  /** Destructeur. */
  virtual ~FormeVdW() ;
  //@}


  /**@name Methods */
  //@{
  /** @brief Renvoie la BBox (boite contenant la forme) de la forme 
  avec prise en compte du rayon de VdW.
  @return La boite englobant la forme. */
  BBox BoxForme() const;

  /** @brief Grandeurs geometriques caracterisant le contact entre les deux 
  formes : point, recouvrement & distance. IMPORTANT: afin de regulariser le
  contact entre 2 formes convexes, on utilise un artefact. Dans un premier
  temps, on realise une compression isotrope de la particule (homothetie) de son
  rayon de Van der Walls (rVdW) grace a la methode getTransformVdW(),
  puis on test le contact, enfin, au lieu d'avoir
  un contact pour d < 0, ou d est la distance minimale entre les 2 convexes, on
  obtient un contact pour d < rVdW(forme1) + rVdW(forme2)  
  @param voisine La forme en voisinage
  @exception Si la penetration depasse le critere VdW
  @return Le conteneur PointContact decrivant le contact */
  PointContact ClosestPoint( FormeVdW &voisine )   
    throw(ErreurContact);
    
  /** @brief Variante de la methode ci-dessus pour gerer les erreurs lies a des
  contacts non detectes par GJK
  @param voisine La forme en voisinage
  @param factor facteur multiplicatif du rayon de Van-der-Waals 
  @exception Si la penetration depasse le critere VdW augmente d'un facteur
  multiplicatif 
  @param id numero du composant 0
  @param id_voisine numero du composant 1 */  
  PointContact ClosestPoint_ErreurHandling( const FormeVdW &voisine,
  	const Scalar& factor, const int& id, const int& id_voisine )   
    throw(ErreurContact);    

  /** @brief Evaluation d'un contact entre deux formes. ATTENTION: test des BBox
  et non veritable contact geometrique, soit probabilite de contact liee a 
  l'intersection des BBox. IMPORTANT: meme artefact que pour la methode 
  ClosestPoint().
  @return TRUE si un contact est present.
  @param voisine La deuxieme forme a etudier. */
  bool isProche( const FormeVdW &voisine ) const; 

  /** @brief Les deux formes sont-elles en contact: veritable contact
  geometrique
  @param voisine La forme voisine
  @return TRUE si le contact est detecte */
  bool isContact( FormeVdW &voisine );
  
  /** @brief Intersection entre 2 formes
  @param a premiere forme
  @param b seconde forme */
  friend bool intersect( const FormeVdW &a, const FormeVdW &b );
  
  /** @brief Initialize a faux le boolean correspondant au calcul de la
  transformation avec scaling par l'epaisseur de croute */
  void initializeVdWtransform_to_notComputed();  
  //@}


  /**@name Methods Get */
  //@{
  /** @brief Valeur du rayon d'approche.
  @return Le rayon d'approche */
  Scalar getRayonInterAction() const;
  
  /** @brief Valeur du rayon d'approche. */
  void setRayonInterAction( Scalar &rayonVdW );

  /** @brief Matrice de transformation avec le rayon de Van-der-Waals.
  Permet d'associer au Convex une matrice de transformation le reduisant
  de la valeur du rayon de Van-der-Waals. */
  Transform const* getTransformVdW();
  
  /** @brief Matrice de transformation avec le rayon de Van-der-Waals.
  Permet d'associer au Convex une matrice de transformation le reduisant
  de la valeur du rayon de Van-der-Waals. 
  @param factor facteur multiplicatif du rayon de Van-der-Waals 
  @param min_scaling le convex ne peut pas etre reduit de plus que min_scaling
  dans chaque direction */
  Transform getTransformVdW( const Scalar& factor, 
  	const Scalar& min_scaling = 1. ) const;  
  //@}
  

  /**@name Methods I/O */
  //@{
  /** @brief Ecriture de l'information statique.
  @param statique Flux recepteur.
  @param composant Pointeur sur le composant pour tester si la particule
  est une particule composite */
  void writeStatique( ostream &statique, Composant const* composant = NULL );
  //@}


protected:
  //@name Parameters */
  //@{  
  Scalar m_rayonVdW; /**< Epaisseur de croute, valeur du recouvrement 
  	autorise */
  Vecteur* m_scaling; /**< Coefficient diagonaux de la matrice de scaling liee a
  	l'epaisseur de croute */
  Transform* m_positionVdW; /** Transformation correspondant a: (1) scaling de 
  	l'epaisseur de croute puis (2) position et orientation de l'objet */
  bool m_positionVdW_computed; /** la transformation m_positionVdW a t elle ete
  	mise a jour au temps courant */			
  //@}


private:

};

/** @brief Point d'intersection cas specifique SPHERE-SPHERE
@param formeA La premiere forme VdW
@param formeB La seconde forme VdW
@exception Si la penetration depasse le critere VdW
@return Le conteneur PointContact decrivant le contact */
PointContact ClosestPointSPHERE( const FormeVdW &formeA,
		const FormeVdW &formeB )   
  throw(ErreurContact);
    
/** @brief Point d'intersection cas specifique SPHERE-BOX
@param formeA La premiere forme VdW
@param formeB La seconde forme VdW
@exception Si la penetration depasse le critere VdW
@return Le conteneur PointContact decrivant le contact */
PointContact ClosestPointSPHEREBOX( const FormeVdW &formeA,
		const FormeVdW &formeB )   
  throw(ErreurContact);
  
/** @brief Y a t il contact, cas specifique SPHERE-SPHERE
@param formeA La premiere forme VdW
@param formeB La seconde forme VdW
@exception Si la penetration depasse le critere VdW
@return Le conteneur PointContact decrivant le contact */
bool isContactSPHERE( const FormeVdW &formeA,
		const FormeVdW &formeB );   
    
/** @brief Y a t il contact, cas specifique SPHERE-BOX
@param formeA La premiere forme VdW
@param formeB La seconde forme VdW
@exception Si la penetration depasse le critere VdW
@return Le conteneur PointContact decrivant le contact */
bool isContactSPHEREBOX( const FormeVdW &formeA,
		const FormeVdW &formeB );   
  
#endif
