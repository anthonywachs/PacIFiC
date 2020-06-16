#ifndef _App
#define _App

#include "Particule.H"
#include "MonObstacle.H"
#include "Basic.H"
#include "Erreur.H"
#include <list>
using namespace std;


/** @brief Description des APPlications utilisees dans le probleme grains.

    Une application permet d'evaluer les forces d'interaction agissant
    sur les particules du milieu.

    Package APP.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class App
{
public:
  /**@name Constructors */
  //@{
  /** Constructeur par defaut */
  App();

  /** Destructeur */
  virtual ~App();
  //@}


  /**@name Methodes virtuelles */
  //@{
  /**
    @brief Calcul des forces d'interaction
    @param dt pas de temps granulaire
    @param time temps physique
    @param particules liste de particules
  */
  virtual void CalculerForces( Scalar dt, Scalar time,
    list<Particule*> const* particules ) = 0;

  /**
    @brief Initialization step
    @param dt pas de temps granulaire
    @param time temps physique
    @param particules liste de particules
  */
  virtual void InitializeTemperature( Scalar dt, Scalar time,
    list<Particule*> const* particules,
    const Scalar fluidThermalConductivity_ ) {}

  /**
    @brief Calcul des forces d'interaction pour l'initialisation des
        conditions de cohesion.
    @param dt Pas de temps granulaire
    @param time Temps physique
    @param particules Liste de particules
  */
  virtual void InitializeCohesiveForces( Scalar dt, Scalar time,
  	list<Particule*> const* particules ) {}
  //@}


  /** @name Methods */
  //@{
  /** @brief Affectation du nom
  @param name_ nom de l'App */
  void setName( const string &name_ ) { m_name = name_; }

  /** @brief Renvoi du nom */
  string getName() const { return m_name; }
  //@}


  /** @name Methods Static */
  //@{
  /** @brief Dimension de l'espace de decomposition
  @param xmax coordonn�e max en x de l'espace de decomposition
  @param ymax coordonn�e max en y de l'espace de decomposition
  @param zmax coordonn�e max en z de l'espace de decomposition
  @param ox origine en x de l'espace de decomposition
  @param oy origine en y de l'espace de decomposition
  @param oz origine en z de l'espace de decomposition */
  static void setD( Scalar xmax, Scalar ymax, Scalar zmax,
  	Scalar ox, Scalar oy, Scalar oz );

  /** @brief Dimension locale de l'espace de decomposition
  @param lx_ Les dimensions de l'espace.
  @param ly_ Les dimensions de l'espace.
  @param lz_ Les dimensions de l'espace. */
  static void setDlocale( Scalar lx_, Scalar ly_, Scalar lz_ );

  /** @brief Origine locale
  @param nprocsdir nombre de processeurs dans chaque direction
  @param MPIcoords coordonnees dans la topologie MPI */
  static void setOriginelocale( int const* nprocsdir, int const* MPIcoords );

  /** @brief Get local origin
  @param x coordonnee x de l'origine locale du proc
  @param y coordonnee y de l'origine locale du proc
  @param z coordonnee z de l'origine locale du proc */
  static void getOrigineLocale( double& x, double& y, double& z );

  /** @brief Get global origin
  @param x coordonnee x de l'origine locale du proc
  @param y coordonnee y de l'origine locale du proc
  @param z coordonnee z de l'origine locale du proc */
  static void getOrigineGlobale( double& x, double& y, double& z );

  /** @brief Get sub-domain dimensions
  @param lx X length of the sub-domain
  @param ly Y length of the sub-domain
  @param lz Z length of the sub-domain */
  static void getDimensionsLocales( double& lx, double& ly, double& lz );

  /** @brief Get domain dimensions
  @param lx X length of the domain
  @param ly Y length of the domain
  @param lz Z length of the domain */
  static void getDimesionsGlobales( double& lx, double& ly, double& lz );

  /** @brief Affiche les attributs statiques
  @param f flux de sortie */
  static void affiche_attributs_statiques( ostream &f );

  /** @brief Appartenance du centre de gravite d'une particule au domaine de
  calcul
  @param position position du centre de gravite */
  static bool isInDomain( Point const* position );

  /**
    @brief Appartenance du centre de gravite d'une particule au domaine local
      de calcul
    @param position position du centre de gravite
  */
  static bool isInLocalDomain( Point const* position );

  /**
    @brief Tells if an app corresponding to the argument exists or not.
    @param name_ The app name we are asking for existence.
    @return true or false
  */
  bool isName( const string name_ );



protected:
  /**@name Parameters */
  //@{
  string m_name; /**< Nom de l'App */
  //@}

  /**@name Parameters Static */
  //@{
  static Scalar m_dX; /**< Dimension de l'espace de decomposition X */
  static Scalar m_dY; /**< Dimension de l'espace de decomposition Y */
  static Scalar m_dZ; /**< Dimension de l'espace de decomposition Z */
  static Scalar m_lX; /**< Dimension locale de l'espace de decomposition X */
  static Scalar m_lY; /**< Dimension locale de l'espace de decomposition Y */
  static Scalar m_lZ; /**< Dimension locale de l'espace de decomposition Z */
  static Point m_origine_globale; /**< Origine geometrique de l'espace de
  	decomposition */
  static Point m_origine_locale; /**< Origine geometrique du sous-domaine */
  //@}
};

#endif
