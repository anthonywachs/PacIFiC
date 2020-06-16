// G.FERRER - Dece.1999 - Creation
// ============================================================================
#ifndef _Plan
#define _Plan

#include "Basic.H"

#include "Point.H"
#include "Vecteur.H"
using namespace solid;


/** @brief Position d'un point par rapport au plan.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation */
// ============================================================================
enum PositionPlan
{
  /**@name Enumeration */
  //@{
  /// Au dessus du plan
  up,
  /// Contact avec le plan
  in,
  /// Au-dessous du plan
  down,
  /// Sans contact avec le plan
  out
  //@}
};


/** Type de deplacements possibles pour un plan.
    
    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation */
// ============================================================================
enum TypeDeplacement
{
  /**@name Enumeration */
  //@{
  /// Deplacement normal au plan (1 composante)
  normal,
  /// Deplacement dans l'espace (3 composantes)
  espace
  //@}
};


/** @brief Classe de gestion des plans.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation */
// ============================================================================
class Plan
{
public:
  /**@name Constructeurs & Destructeur */
  //@{
  /** @brief Constructeur par defaut. */
  Plan();
  
  /** @brief Constructeur avec initialisation par trois points dans le plan. 
  @param nbre Nombre de points 
  @param xPt Coordonnees x des points
  @param yPt Coordonnees y des points
  @param zPt Coordonnees z des points */
  Plan(int nbre, const Scalar *xPt, const Scalar *yPt, const Scalar *zPt);

  /** @brief Constructeur avec initialisation par trois points dans le plan.
  @param p1 1er point
  @param p2 2eme point
  @param p3 3eme point */
  Plan(const Point& p1, const Point& p2, const Point& p3);

  /** @brief Constructeur avec initialisation.
  Description du plan par ces equations. Ne contenant pas les points 
  de description, le deplacement est impossible.
  @param normale Description de la normale et direction */
  Plan(const Scalar *normale);
  
  /** @brief Constructeur par copie.
  @param plan Objet a copie */
  Plan(const Plan &plan);

  /** @brief Destructeur. */
  ~Plan();
  //@}


  /**@name Methods */
  //@{
  /** @brief Définition du plan à partir de 3 points
  @param p1 point 1
  @param p2 point 2    
  @param p3 point 3 */  
  void  setPlan(const Point& p1, const Point& p2, const Point& p3);
  
  /** @brief Deplacement de l'objet d'une valeur scalaire.
  Le point se deplace de cette valeur dans les trois directions.
  @param dist Valeur du deplacement */
  void Deplacer(Scalar dist);
  
  /** @brief Deplacement de l'objet d'un tenseur pour les trois directions.
  @param dist Valeur des deplacements */
  void Deplacer(const Scalar *dist);
  
  /** @brief Deplacement de l'objet d'un vecteur pour les trois directions.
  @param dist Valeur des deplacements */
  void Deplacer(const Vecteur &dist);
  
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
  @param pt Point a etudier */
  Scalar DistanteDe(const Scalar *pt) const;
  
  /** @brief Distance de l'objet a un point.
  @return La distance avec son signe
  @param x_ Position en X
  @param y_ Position en Y
  @param z_ Position en Z */
  Scalar DistanteDe(Scalar x_, Scalar y_, Scalar z_) const;

  /** @brief Coordonnees d'un point en projection sur le plan.
  @return Le point appartenant au plan.
  @param point Le point dans l'espace. */
  Point Projection(const Point &point) const;
  
  /** @brief Coordonnees d'un point en projection sur le plan.
  @return Le point appartenant au plan.
  @param point Les coordonnes du point dans l'espace. */
  Point Projection(const Scalar *point) const;
  
  /** @brief Coordonnees d'un point en projection sur le plan.
  @return Le point appartenant au plan.
  @param x_ Coordonnes du point dans l'espace. 
  @param y_ Coordonnes du point dans l'espace. 
  @param z_ Coordonnes du point dans l'espace. */
  Point Projection(const Scalar x_,const Scalar y_,const Scalar z_) const;

  /** @brief Validation de la paroi.
  La distance au point doit etre positive par rapport au point indique. 
  Si cela n'est pas le cas, on inverse les sommets 1 et 2.
  @param point Coordonnees du point de validation */
  void Valider(Scalar *point);

  /** @brief Position du point par rapport au plan.
  @return Positionnement
  @param  point Coordonnees du point */
  PositionPlan WhereIs(const Point &point) const;
  
  /** @brief Position du point par rapport au plan.
  @return Positionnement
  @param  pt Coordonnees du point */
  PositionPlan WhereIs(const Scalar *pt) const;

  /** @brief Position du point par rapport au plan.
  @return Positionnement
  @param  x_ coordonnee X
  @param  y_ coordonnee Y
  @param  z_ coordonnee Z */
  PositionPlan WhereIs(Scalar x_, Scalar y_, Scalar z_) const;
  //@}


protected:
  /**@name Parameters */
  //@{
  int nbrePts; /**< Nombre de points */
  Scalar *x; /**< Coordonnees X des points */
  Scalar *y; /**< Coordonnees Y des points */
  Scalar *z; /**< Coordonnees Z des points */
  Scalar a; /**< Description de la normale */
  Scalar b; /**< Description de la normale */
  Scalar c; /**< Description de la normale */
  Scalar d; /**< Decalage du plan par rapport a l'origine : ax + by + cz = d */
  //@}
};

#endif
