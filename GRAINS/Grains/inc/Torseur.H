#ifndef _Torseur
#define _Torseur

#include "Vecteur.H"
#include "Point.H"
using namespace solid;


class Torseur;
ostream &operator << (ostream &fileOut, const Torseur &objet);
istream &operator >> (istream &fileIn, Torseur &objet);


/** Objet decrivant les efforts exerces sur une entite geometrique.
   
    Il indique le point d'application de la force, son intensite et son moment.

    @author F.PRADEL - Institut Francais du Petrole - 1999 - Creation */
// ============================================================================
class Torseur
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  Torseur();

  /** @brief Constructeur par initialisation.
  @param pt Point de reduction
  @param f Vecteur force
  @param m Vecteur de moment */
  Torseur(const Point &pt, const Vecteur &f, const Vecteur &m);

  /** @brief Constructeur par initialisation.
  Le moment reste initialise a 0.0
  @param pt Point de reduction
  @param f Vecteur force */
  Torseur(const Point &pt, const Vecteur &f);

  /** @brief Constructeur par initialisation.
  La force et le moment reste initialise a 0.0
  @param pt Point de reduction */
  Torseur(const Point &pt);

  /** @brief Constructeur par initialisation.
  @param point Point de reduction
  @param load Force
  @param momentum Moment */
  Torseur(double* point, double* load, double* momentum );

  /** @brief Constructeur de copie.
  @param tau Le torseur de reference. */
  Torseur(const Torseur &tau);

  /** @brief Destructeur */
  ~Torseur();
  //@}


  /**@name Methods */
  //@{
  /** @brief Changement du point de reduction du Torseur 
  @param newpoint le nouveau point de reduction */
  void Changer_Point(const Point &newpoint);
  
  /** @brief Ajouter une force qui agit au point de reduction du torseur
  @param force_ La nouvelle force */
  void addForce(const Vecteur &force_); 
  
  /** @brief Ajouter une force qui agit à un point different du point de 
  reduction du torseur
  @param force_ La nouvelle force 
  @param point point d'action de la force */
  void addForce(const Point &point,const Vecteur &force_);   
  
  /** @brief Ajouter une force qui agit au point de reduction du torseur
  @param fx composante x de la force 
  @param fy composante y de la force   
  @param fz composante z de la force */
  void addForce(const double &fx,const double &fy,const double &fz);   
  
  /** @brief Ajouter un moment
  @param moment_ La nouvelle force */
  void addMoment(const Vecteur &moment_);
  
  /** @brief Ajouter un moment
  @param mx composante x du moment
  @param my composante y du moment
  @param mz composante z du moment */
  void addMoment(const double &mx,const double &my,const double &mz);   
  
  /** @brief Copie force & moment du torseur dans le vecteur fm 
  en débutant à la position i 
  @param fm vecteur de copie 
  @param i position dans le vecteur fm */
  void copyForceMoment(double *fm,int i) const;         
  //@}


  /**@name Methods Get */
  //@{
  /** @brief Renvoi la force du Torseur.
  @return La composante de force. */
  Vecteur const* getForce() const;

  /** @brief Renvoi la force du Torseur.
  @return La composante de moment. */
  Vecteur const* getMoment() const;

  /** @brief Renvoi le point de reduction du Torseur.
  @return Le point demande. */
  Point const* getPoint() const;
  //@}  


  /**@name Methods Set */
  //@{
  /** @brief Affectation de la force du torseur.
  @param force_ La nouvelle force */
  void setForce(const Vecteur &force_);

  /** @brief Affectation du point de reduction.
  @param point Le nouveau point de reduction. */
  void setPoint(const Point &point);
  
  /** @brief Affectation de la force, du point de reduction et initialisation
  du moment a zero
  @param point Le nouveau point de reduction. 
  @param force_ La nouvelle force */
  void setToBodyForce(const Point &point,const Vecteur &force_);    
  //@}  


  /**@name Operators */
  //@{
  /** @brief Operateur unaire - qui renvoie l'oppose des vecteurs 
  mais garde le point de reduction intact.
  @return Le torseur oppose au courant. */
  Torseur operator-();

  /** @brief Multiplication du torseur par un scalaire.
  @return Le torseur resultant.
  @param d Le scalaire multipliant la force et le moment. */
  Torseur operator * (double d);

  /** @brief Addition de deux Torseurs au point de reduction du premier.
  @return Le torseur resultant.
  @param rhs Le torseur a ajouter au torseur courant. */
  Torseur operator+(Torseur &rhs);

  /** @brief Soustraction de deux Torseurs au point de reduction du premier.
  @return Le torseur resultant.
  @param rhs Le torseur a soustraire au torseur courant. */
  Torseur operator-(Torseur& rhs);
  
  /** @brief Comparaison entre deux torseurs par les valeurs.
  @return Vrai si les deux torseurs sont egaux.
  @param top2 Le torseur a comparer. */
  bool operator == (Torseur& top2);

  /** @brief Difference entre deux torseurs par les adresses.
  @return Vrai si les deux torseurs sont identiques.
  @param top2 Le torseur a comparer. */
  bool operator != (Torseur& top2);

  /** @brief Egalite membre a membre.
  @return Le torseur courant avec les nouvelles valeurs.
  @param rhs Le torseur de reference. */
  Torseur& operator=(const Torseur & rhs);

  /** @brief Incrementation du torseur courant.
  @return Le torseur courant incremente.
  @param rhs Le torseur a ajouter au torseur courant. */
  Torseur& operator+=(const Torseur& rhs);
  
  /** @brief Incrementation du torseur courant avec un torseur dont on specifie
  le point de reduction pour le calcul du la contribution au moment provenant du
  bras de levier ^ force appliquee. Utile pour les particules periodiques.
  @return Le torseur courant incremente.
  @param rhs Le torseur a ajouter au torseur courant. 
  @param rp_rhs le point de reduction du torseur a ajouter */
  void addWithReductionPoint(const Torseur& rhs,const Point& rp_rhs);  

  /** @brief Decrementation du torseur courant
  @return Le torseur courant decremente.
  @param rhs Le torseur a soustraire au torseur courant. */
  Torseur& operator-=(Torseur& rhs);
  //@}


  /**@name Methods I/O */
  //@{
  /** @brief Ecriture de l'information sous une forme de force.
  Des textes de commentaires sont inseres : *PtReduction...
  @param fileOut Flux recepteur. */
  void write(ostream &fileOut);
  //@}


  /**@name Methods Friend */
  //@{
  /** @brief Operateur d'ecriture.
  @return Flux recepteur
  @param fileOut Flux recepteur
  @param objet Objet courant. */  
  friend ostream &operator << (ostream &fileOut, const Torseur &objet);

  /** @brief Operateur de lecture.
  @return Flux emetteur
  @param fileIn Flux emetteur
  @param objet Objet courant */
  friend istream &operator >> (istream &fileIn, Torseur &objet);
  //@}


private:
  /**@name Parameter */
  //@{
  Point ptReduction; /**< Point de reduction du torseur */
  Vecteur force; /**< Force du torseur */
  Vecteur moment; /**< Moment du torseur */
  //@}
};

#endif
