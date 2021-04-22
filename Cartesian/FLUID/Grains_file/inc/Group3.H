// D.PETIT - Juin.2000 - Creation
// ============================================================================
#ifndef _Group3
#define _Group3

#include "Basic.H"

#include <iostream>
using namespace std;


class Group3;
ostream &operator << (ostream &fileOut, const Group3 &objet);
istream &operator >> (istream &fileIn, Group3 &objet);


/** @brief GJK Engine - A Fast and Robust GJK Implementation 
  
    
@author Copyright (C) 1998  Gino van den Bergen. */
// ============================================================================
namespace solid
{
  /** @brief Classe simulant un anneau de R3.
      
  @author D.PETIT - 2000 - Creation */
  // ==========================================================================
  class Group3
  {
  public:
    /**@name Constructors */
    //@{
    /** @brief Constructeur par defaut. 
    @param def valeur des composantes */
    Group3( const Scalar def = 0. );

    /** @brief Constructeur avec initialisation. 
    @param x premier element
    @param y deuxieme element
    @param z troisieme element */
    Group3( const Scalar x, const Scalar y, const Scalar z );

    /** @brief Constructeur par copie.
    @param g l'element a copie. */
    Group3( const Group3 &g );

    /** @brief Destructeur. */
    ~Group3();
    //@}
  

    /**@name Methods Get */
    //@{ 
    /** @brief Renvoi l'element. 
    @return pointeur sur l'element.*/
    const Scalar *getValue() const;
    
    /** @brief Renvoi l'element.
    @return pointeur sur l'element.*/
    Scalar *getValue();    

    /** @brief nombre de composantes d'un element. 
    @return le nombre de composantes.*/
    int size() const;    
    //@}
    
    
    /**@name Methods Set */
    //@{ 
    /** @brief Modifie les composantes d'un element. 
    @param x premier element
    @param y deuxieme element
    @param z troisieme element */
    void setValue( const Scalar x, const Scalar y, const Scalar z );
    
    /** @brief Modifie les composantes d'un element. 
    @param g element a copie. */
    void setValue( const Scalar* g ); 
    
    /** @brief Annule toutes les composantes d'un element */
    void reset();            
    //@}


    /**@name Methods */
    //@{ 
    /** @brief Ecriture de l'objet sur le flux de sortie avec un format de
    precision elevee
    @param fileOut Flux de sortie. */
    void writeGroup3( ostream &fileOut ) const;
    
    /** @brief Ecriture de l'objet en binaire sur le flux de sortie
    @param fileOut Flux de sortie. */
    void writeGroup3_binary( ostream &fileOut ); 
    
    /** @brief Lecture de l'objet en binaire sur le flux d'entree
    @param StreamIN Flux d'entree. */
    void readGroup3_binary( istream &StreamIN );            
    //@}
                 

    /**@name Operators. */
    //@{
    /** @brief Operateur d'indexation. 
    @return La composante i.
    @param i  indice de la composantes de l'element */
    Scalar& operator[]( const int i );
    
    /** @brief Operateur d'indexation. 
    @return la composante i
    @param i  indice de la composantes de l'element */
    const Scalar& operator[]( const int i ) const;

    /** @brief Operateur unaire - 
    @return Objet avec les composantes opposees. */
    Group3 operator - () const;

    /** @brief Produit scalaire de deux objets de type Group3.
    @return Le produit scalaire.
    @param g Le deuxieme objet */
    Scalar operator * ( const Group3 &g ) const;
    
    /** @brief Multiplication par un scalaire
    @return Le nouvel objet
    @param d Le scalaire multipliant les composantes */
    Group3 operator * ( Scalar d ) const;

    /** @brief Division du Group3 par un scalaire.
    @return Le Group3 resultat
    @param d Le scalaire. */
    Group3 operator / ( Scalar d ) const;

    /** @brief Addition de deux Group3 composante par composante
    @return La somme des deux objets.
    @param g2 le deuxieme objet */
    Group3 operator + ( const Group3 &g2 ) const;

    /** @brief Soustraction de deux Group3 composante par composante
    @return La difference des deux objets.
    @param g2 le deuxieme objet */
    Group3 operator - ( const Group3 &g2 ) const;

    /** @brief Comparaison pour les valeurs
    @param g2 le deuxieme objet */
    bool operator == (const Group3 &g2) const;
    
    /** @brief Objets differents par les valeurs 
    @param g2 le deuxieme objet */
    bool operator != ( const Group3 &g2 );

    /** @brief Egalite d'affectation d'un vecteur 
    @param g2 le deuxieme objet */
    Group3& operator = ( const Group3 &g2 );
    
    /** @brief Affectation d'une valeur aux composantes du vecteur.
    @param valeur Valeur a affecter. */
    void operator = ( Scalar valeur );

    /** @brief Multiplication par un scalaire et affectation 
    @param d Le scalaire. */
    Group3& operator *= ( Scalar d );

    /** @brief Division du vecteur courant par un scalaire.
    @return Le vecteur courant.
    @param d Le scalaire. */
    Group3& operator /= ( Scalar d );

    /** @brief Addition d'un vecteur et affectation
    @param g2 le deuxieme objet */
    Group3& operator += ( const Group3 &g2 );

    /** @brief Soustraction d'un vecteur et affectation
    @param g2 le deuxieme objet */
    Group3& operator -= ( const Group3 &g2 );
    //@}


    /**@name Methods Friend */
    //@{
    /** @brief Produit mixte de trois Group3
    @return le resultat du produit mixte
    @param g1 Premier Group3
    @param g2 Deuxieme Group3
    @param g3 Troisieme Group3 */
    friend Scalar triple( const Group3 &g1, 
	const Group3 &g2, const Group3 &g3 );

    /** @brief Multiplication d'un scalaire par un Group3 
    @return Le Group3 = d * g
    @param d Un scalaire
    @param g Un Group3 */
    friend Group3 operator * ( Scalar d, const Group3 &g );

    /** @brief Operateur d'ecriture.
    @return Flux recepteur
    @param fileOut   Flux recepteur
    @param objet Objet courant. */
    friend ostream &operator << ( ostream &fileOut, const Group3 &objet );

    /** @brief Operateur de lecture.
    @return Flux emetteur
    @param fileIn    Flux emetteur
    @param objet Objet courant */
    friend istream &operator >> ( istream &fileIn, Group3 &objet );
    //@}


  protected:
    /**@name Parameters */
    //@{
    Scalar comp[3]; /**< Les trois composantes d'un element. */
    //@}
  };

  /** @brief  Produit triple
  @param g1 groupe 1
  @param g2 groupe 2  
  @param g3 groupe 3 */ 
  Scalar triple( const Group3 &g1, const Group3 &g2, const Group3 &g3 );

  /** @brief Multiplication par un scalaire
  @param d Le scalaire
  @param g groupe */  
  Group3 operator * ( Scalar d, const Group3 &g );

  static Group3 OrigineGroup3; /**< Origine (0.,0.,0.)  */

} // namespace solid

#endif
