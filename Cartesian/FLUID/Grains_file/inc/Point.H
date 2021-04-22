#ifndef _Point
#define _Point

#include "Group3.H"
#include "Vecteur.H"


namespace solid
{
  /** @brief Classe de gestion des points 3D.
      
  @author G.FERRER - Institut Francais du Petrole - 1999 - Creation 
  @author F.PRADEL - Institut Francais du Petrole - 2000 - Modification
  @author D.PETIT  - 2000 - Modification 
  @author A.WACHS  - 2009 - Modification */
  // ==========================================================================
  class Point : public Group3
  {
  public:
    /**@name Constructors */
    //@{
    /** @brief Constructeur par defaut. 
    @param def Valeurs de defaut des composantes du vecteur. */
    Point(Scalar def = 0.);

    /** @brief gets the shrinking choice. */
    int getShrinkingChoice()const;

    /** @brief Constructeur avec initialisation.
    @param x Position en x
    @param y Position en y
    @param z Position en z */
    Point(Scalar x, Scalar y, Scalar z);

    /** @brief Constructeur par copie.
    @param g Objet a copier */
    Point(const Group3 &g);

    /** @brief Destructeur. */
    ~Point();
    //@}


    /**@name Methods */
    //@{
    /** @brief Deplacement de l'objet d'une valeur scalaire.
    Le point se deplace de cette valeur dans les trois directions.
    @param dist Valeur du deplacement */
    void Deplacer(Scalar dist);
    
    /** @brief Deplacement de l'objet à l'aide d'un tableau de 3 valeurs.
    @param dist Valeur des deplacements */
    void Deplacer(const Scalar *dist);
    
    /** @brief Deplacement de l'objet avec les trois valeurs.
    @param distX Deplacement de dist en X
    @param distY Deplacement de dist en Y
    @param distZ Deplacement de dist en Z */
    void Deplacer(Scalar distX, Scalar distY, Scalar distZ);

    /** @brief Distance entre deux points.
    Le signe de la distance this-point est conserve.
    @return Valeur de la distance
    @param point Deuxieme point */    
    Scalar DistanteDe(const Point &point) const;
    
    /** @brief Distance de l'objet a un point.
    @return La distance avec son signe
    @param point Point a etudier */
    Scalar DistanteDe(const Scalar *point) const;
    
    /** @brief Distance de l'objet a un point.
    @return La distance avec son signe
    @param x Position en X
    @param y Position en Y
    @param z Position en Z */
    Scalar DistanteDe(Scalar x, Scalar y, Scalar z) const;
    //@}
  };
  
  static Point OriginePoint; /**< Origine (0.,0.,0.)  */
} // namespace solid

#endif

