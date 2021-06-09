#ifndef _CylindricalBoxObstacle
#define _CylindricalBoxObstacle

#include "CompObstacle.H"
using namespace std;

#include "ReaderXML.hh"



/**
@brief Class which refers to the composed obstacle of type CylindricalBox
@details It permits to create and manage obstacles inheritating from
  CompObstacle.<BR>
  A %CylindricalBox is composed of a multiple panels assembly.<BR>
  Properties of this assembly (Length, Diameter, Nb_panels, etc) are read in 
    input-file.

    @author M.BERNARD - IFPEN - March 2015 - Creation
*/
// ============================================================================
class CylindricalBoxObstacle : public CompObstacle
{
public:
  /** @name Constructeurs */
  //@{
  /**
    @brief Constructeur par defaut.
    @param s Obstacle name
    @param radius_ Obstacle radius
    @param length_ Obstacle length
  */
  CylindricalBoxObstacle( const string &s = "", const double &radius_ = 0. ,
    const double &length_ = 0. );

  /** @brief Constructeur par decodage d'un noeud XML
  @param root Le noeud "Obstacle" */
  CylindricalBoxObstacle( DOMNode *root );
  
  /** @brief Destructeur.*/
  ~CylindricalBoxObstacle();
  //@}


  /** @name Methods Get */
  //@{
  /** @brief Renvoie le type de l'obstacle */
  string getObstacleType();

  /** @brief Liste des obstacles à transmettre au fluide
  @return The CylindricalBoxObstacle itself */
  list<Obstacle*> getObstaclesToFluid();
  
  /** @brief Points decrivant l'enveloppe
  @return La liste des points sommets */
  vector<Point> getEnveloppe() const;
  //@}

  /** @name Methods */
  //@{
  //@}

  /** @name Methods I/O */
  //@{
  /** @brief Sauvegarde des obstacles pour Reload
  @param fileSave Flux de sauvegarde
  @param composant solid object */
  virtual void write( ostream &fileSave, Composant const* composant = NULL )
      const; 

  /** @brief Reload du groupe d'Obstacle & publication dans son referent
  @param obstacle L'obstacle referent du courant
  @param file Le fichier de persistance */
  void reload( Obstacle &obstacle, istream &file );
  //@}
  


protected:
  /** @name Methods */
  //@{
  /** @brief Generates side panels from informations read in the constructor */
  void generate_sidePanels();

  /** @brief Generates top and/or bottom panels from informations read in 
    the constructor */
  void generate_topAndBotomPanels();
  //@}

  /** @name Parameters */
  //@{  
  double m_length; /**< Height of the CylindricalBox */
  double m_radius; /**< Radius of the CylindricalBox */
  Point m_center; /**< Center of the CylindricalBox */
  double m_thickness; /**< Walls thickness of the CylindricalBox */
  size_t m_nbPanels; /**< Number of panels forming the CylindricalBox */
  double m_overlapMax; /**< Maximum overlap distance, will be set to all 
  	elementary MonObstacle */
  bool m_transfer2Fluid; /**< Obstacle à transferer au fluide dans le cas d'un
  	couplage */
  bool m_bottomPanel; /**< Does the box has a bottom face ? */
  bool m_topPanel; /**< Does the box has a top face ? */
  double m_polygonAngle; /**< Angle of the polygon approching the CylindricalBox
  section. ( \f$ \theta = \frac{2\pi}{nbPanels} ) \f$ */
  //@}
    
};

#endif

