// F.Pradel - Janvier 2000 - Creation
// ============================================================================
#ifndef _Quaternion
#define _Quaternion

#include "Vecteur.H"
using namespace solid;
#include <iostream>
using namespace std;
#include "WriterXML.hh"


class Quaternion;
ostream &operator << (ostream &fileOut, const Quaternion &objet);
istream &operator >> (istream &fileIn, Quaternion &objet);


/** @brief Gestion des quaternions.

    il s'agit d'un couple d'un scalaire plus un Vecteur.

    @author F.PRADEL - Institut Francais du Petrole - 2000 - Modification */
// ============================================================================
class Quaternion
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut. */
  Quaternion();

  /** @brief Constructeur avec initialisation.
  @param q Valeur du vecteur x,y,z
  @param d Scalaire du Quaternion */
  Quaternion( const Scalar q, const Scalar d=0.0 );

  /** @brief Construction a partir d'un 
  @param vecteur Le vecteur et avec 
  @param d scalaire par defaut nul */ 
  Quaternion( const Vecteur& vecteur, const Scalar d=0.0 );

  /** @brief Constructeur avec initialisation
  @param x Position en i
  @param y Position en j
  @param z Position en k
  @param d Scalaire par defaut nul */
  Quaternion( Scalar x, Scalar y, Scalar z, Scalar d );

  /** @brief Constructeur par copie.
  @param q Objet a copie */
  Quaternion( const Quaternion& q );

  /** @brief Destructeur. */
  ~Quaternion();
  //@}


  /**@name Methods */
  //@{
  /** @brief Creation du conjugue du quaternion 
  @return Le quaternion conjugue. */
  Quaternion Conjugate() const;

  /** @brief Creation de l'inverse du quaternion.
  @return Le quaternion inverse. */
  Quaternion Inverse() const;
  
  /** @brief Produit d'un vecteur avec un quaternion (le vecteur lhs correspond 
  à un quaternion [0,lhs]) 
  @return Le produit du vecteur et du quaternion
  @param lhs Le vecteur */
  Quaternion multLeftVec( const Vecteur& lhs ) const;
  
  /** @brief Produit de 2 quaternions dont on renvoie la partie vectorielle
  @param rhs quaternion à multiplier */  
  Vecteur multToVecteur( const Quaternion& rhs ) const; 
  
  /** @brief Produit d'un quaternion q par le conjugue du second quaternion q',
  c.a.d. q*conjufue(q'), dont on renvoie la partie vectorielle
  @param rhs quaternion à multiplier */  
  Vecteur multConjugateToVecteur( const Quaternion& rhs ) const; 
  
  /** @brief Ecriture de l'objet sur le flux de sortie avec un format de
  precision elevee
  @param fileOut Flux de sortie. */
  void writeQuaternion( ostream &fileOut ) const;
  
  /** @brief Ecriture de l'objet en binaire sur le flux de sortie
  @param fileOut Flux de sortie. */
  void writeQuaternion_binary( ostream &fileOut ); 

  /** @brief Lecture de l'objet en binaire sur le flux d'entree
  @param StreamIN Flux d'entree. */
  void readQuaternion_binary( istream &StreamIN );              
  //@}


  /**@name Methods Get */
  //@{ 
  /** @brief Acces au scalaire composant le Quaternion
  @return le scalaire demande. */
  Scalar getScalaire() const;

  /** @brief Acces au vecteur composant le Quaternion.
  @return Le vecteur demande. */
  Vecteur const* getVecteur() const;

  /** @brief Prend le quaternion et construit le vecteur rotation associe 
  qui est renvoye 
  @return Vecteur */
  Vecteur getVecteurRotation() const;    
  //@}


  /**@name Methods Set */
  //@{ 
  /** @brief Modification du quaternion.
  @param vecteur  Composante du vecteur.
  @param scalaire Composante du scalaire. */
  void setQuaternion( const Vecteur &vecteur, Scalar scalaire );

  /** @brief Modification du quaternion a partir d'un vecteur rotation
  @param vecteur Composante du vecteur */
  void setQuaternion( const Vecteur &vecteur );

  /** @brief Modification du quaternion.
  @param vecteur0 Composante 0 du vecteur
  @param vecteur1 Composante 1 du vecteur
  @param vecteur2 Composante 2 du vecteur
  @param scalaire Composante du scalaire */
  void setQuaternion( const Scalar &vecteur0, 
	const Scalar &vecteur1, 
	const Scalar &vecteur2, 
	Scalar scalaire );

  /** @brief Modification du quaternion.
  @param scalaire Composante du scalaire */
  void setScalaire( Scalar scalaire );

  /** @brief Modification du quaternion.
  @param vecteur Composante du vecteur. */
  void setVecteur( const Vecteur &vecteur );
  //@}
  

  /**@name Operators */
  //@{
  /** @brief Indexation: (0,1,2) pour les composantes du vecteur et 3 pour le
  scalaire
  @param i index */
  Scalar& operator[]( int i );

  /** @brief Indexation: (0,1,2) pour les composantes du vecteur et 3 pour le
  scalaire
  @param i index */
  Scalar operator[]( int i ) const;

  /** @brief Produit de deux quaternions 
  @return Le produit des deux quaternions
  @param rhs Le deuxieme quaternion */
  Quaternion operator*( const Quaternion& rhs ) const;

  /** @brief Produit d'un quaternion avec un vecteur (le vecteur rhs correspond 
  à un quaternion [0,rhs]) 
  @return Le produit du quaternion et du vecteur
  @param rhs Le vecteur */
  Quaternion operator,( const Vecteur& rhs ) const;

  /** @brief Multiplication par un scalaire 
  @param d scalaire */
  Quaternion operator*( Scalar d );

  /** @brief Operateur - unaire */ 
  Quaternion operator-();

  /** @brief Addition de deux quaternions composante par composante 
  @param rhs quaternion à ajouter */
  Quaternion operator+( const Quaternion& rhs ) const;

  /** @brief Soustraction de deux quaternions composante par composante 
  @param rhs quaternion à soustraire */
  Quaternion operator-( const Quaternion& rhs );

  /** @brief operateur de comparaison pour les valeurs 
  @param rhs quaternion à comparer */
  bool operator == ( const Quaternion& rhs );

  /** @brief operateur de difference pour les valeurs 
  @param rhs quaternion à comparer */
  bool operator != ( const Quaternion& rhs );

  /** @brief operateur d'affectation 
  @param rhs quaternion à affecter */
  Quaternion& operator = ( const Quaternion& rhs );

  /** @brief operateur d'affectation des composantes d'un quaternion a un 
  scalaire 
  @param rhs scalaire */
  Quaternion operator = ( const Scalar rhs );

  /** @brief operateur de multiplication par un scalaire et d'affectation
  @param d scalaire */
  Quaternion& operator*=( Scalar d );

  /** @brief operateur d'addition d'un point et d'affectation
  @param rhs quaternion à ajouter */
  Quaternion& operator+=( const Quaternion& rhs );

  /** @brief operateur de soustraction d'un point et d'affectation
  @param rhs quaternion à soustraire */
  Quaternion& operator-=( const Quaternion& rhs );
  //@}


  /**@name Methods Friend */
  //@{
  /** @brief Operateur d'ecriture 
  @param fileOut flux de sortie
  @param objet le quaternion */ 
  friend ostream &operator << ( ostream &fileOut, const Quaternion &objet );
  
  /** @brief Operateur de lecture 
  @param fileIn flux d'entrée
  @param objet le quaternion */   
  friend istream &operator >> ( istream &fileIn, Quaternion &objet );

  /** @brief Norme du Quaternion.
  @return Valeur de la norme.
  @param qt Le Quaternion etudie. */
  friend Scalar ABS( const Quaternion& qt );

  /** @brief Norme d'un quaternion
  @return la norme du quaternion
  @param qt le quaternion */
  friend Scalar Norm( const Quaternion& qt );

  /** @brief Norme du quaternion au carre
  @return la norme du quaternion au carre
  @param qt le quaternion */
  friend Scalar Norm2( const Quaternion& qt );
  //@}


protected:
  /**@name Parameters */
  //@{  
  Scalar w; /**< Scalaire du quaternion */  
  Vecteur vqt; /**< Composante du vecteur du quaternion */
  //@}
};

/**@name Quaternion : Methods external */
//@{
/** @brief Multiplication d'un scalaire par un quaternion 
@param d scalaire 
@param rhs quaternion à multiplier */
Quaternion operator*( Scalar d, const Quaternion& rhs );

/** @brief Produit d'un vecteur avec un quaternion (le vecteur lhs correspond 
à un quaternion [0,lhs]) 
@return Le produit du vecteur et du quaternion
@param lhs Le vecteur 
@param q Le quaternion */
Quaternion operator,( const Vecteur& lhs, const Quaternion& q );

/** @brief Valeur absolue d'un quaternion 
@param qt le quaternion */
Scalar ABS( const Quaternion& qt );
//@}

#endif

