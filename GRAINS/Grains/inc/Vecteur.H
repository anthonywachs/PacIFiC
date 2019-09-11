#ifndef _Vecteur
#define _Vecteur

#include "Group3.H"

#include "WriterXML.hh"

class Quaternion;


namespace solid
{
  /** @brief Classe de gestion des vecteurs 3D.
      
  @author G.FERRER - Institut Francais du Petrole - 1999 - Creation 
  @author F.PRADEL - Institut Francais du Petrole - 2000 - Modification */
  // ==========================================================================
  class Vecteur : public Group3
  {
  public:
    /** @name Constructeurs */
    //@{
    /** @brief Constructeur par defaut. 
    @param def Valeurs de defaut des composantes du vecteur. */
    Vecteur(Scalar def = 0.);

    /** @brief Constructeur avec initialisation.
    @param x Position en x
    @param y Position en y
    @param z Position en z */
    Vecteur(Scalar x, Scalar y, Scalar z);

    /** @brief Constructeur par copie.
    @param g Objet a copie */
    Vecteur(const Group3 &g);

    /** @brief Destructeur. */
    ~Vecteur();
    //@}


    /** @name Methods */
    //@{
    /** @brief Determine l'axe le plus proche de la 'pointe' du Vecteur.
    @return Entier designant l'axe le plus proche (0->X, 1->Y, 2->Z) */
    int closestAxis() const;

    /** @brief Normalise le Vecteur. */
    void normalize();

    /** @brief Normalise un Vecteur. 
    @return Le vecteur normalise */
    Vecteur normalized() const;

    /** @brief Troncature des decimales de chaque composante a +-EPS */
    void round();

    /** @brief Fonction faisant tourner le vecteur au moyen du quaternion
    @param q Quaternion representant une rotation*/
    void Tourne(const Quaternion& q);

    /** @brief Fonction faisant tourner le vecteur au moyen du quaternion
    @param q Quaternion representant une rotation*/
    Vecteur Tourne(const Quaternion& q) const;
    //@}


    /** @name Operators */
    //@{
    /** @brief Operateur d'affectation 
    @param g2 vecteur à affecter */
    Vecteur& operator = (const Vecteur &g2);

    /** @brief Produit vectoriel de deux vecteurs 
    @param rhs le 2è vecteur */
    Vecteur operator ^ (const Vecteur& rhs) const;        
    //@}


    /** @name Methods Friend */
    //@{
    /** @brief Test si la longueur du vecteur est nul.
    @return Vrai si le vecteur est approximativement nul
    @param v Le vecteur a tester */ 
    friend bool approxZero(const Vecteur& v);

    /** @brief Determine le cosinus de l'angle formee entre deux vecteurs.
    @return Le cos de l'angle forme par deux vecteurs
    @param x Le premier vecteur
    @param y Le deuxieme vecteur */
    friend Scalar cos(const Vecteur& x, const Vecteur& y);

    /** @brief Determine la norme d'un vecteur
    @return le scalaire de la norme
    @param v le vecteur dont on veut la norme */
    friend Scalar Norm(const Vecteur& v);

    /** @brief Determine la norme au carre d'un vecteur 
    @return le scalaire de la norme
    @param v le vecteur dont on veut la norme au carre */
    friend Scalar Norm2(const Vecteur& v);
    //@}
  };

  /** @brief Determine la norme d'un vecteur
  @return le scalaire de la norme
  @param v le vecteur dont on veut la norme */        
  Scalar Norm(const Vecteur& v);

  static Vecteur VecteurNul; /**< Vecteur nul (0.,0.,0.)  */   
} // namespace solid

#endif

