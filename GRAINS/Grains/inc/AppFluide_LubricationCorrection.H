#ifndef _AppFluide_LubricationCorrection
#define _AppFluide_LubricationCorrection

#include "App.H"
#include "LinkedCell.H"
#include "ReaderXML.hh"


/** @brief Force de correction de lubrification

    This force is subjected to the particles which are close to each other. 
    The reason is due to the sub-grid viscouse effects of the film of fluid 
    squeezing out of the particle/particle or particle/wall gap. 

    @author A.ESTEGHAMATIAN  2014 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/** @brief Resistance matrix for particle-particle and particle-wall approach.
Please refer to Dance and Maxy paper for definitions.

@author A.Esteghamatian - IFP Energies nouvelles - 2014 - Creation */
// ============================================================================
struct ResistanceMatrix
{
  Scalar A11_1;
  Scalar A11_2;
  Scalar A22_1;
  Scalar A22_2;
  Scalar B23_1;
  Scalar B23_2;
  Scalar B32_1;
  Scalar B32_2;
  Scalar C23_1;
  Scalar C23_2;
  Scalar C32_1;
  Scalar C32_2;
  Scalar D11_1;
  Scalar D11_2;
  Scalar D22_1;
  Scalar D22_2;
  Scalar A11;
  Scalar A22;
  Scalar B23;
  Scalar B32;
  Scalar C23;
  Scalar C32;
  Scalar D11;
  Scalar D22;
};

class AppFluide_LubricationCorrection : public App
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur 
  @param gridsize taille de maille du maillage fluide 
  @param is_DNS is DNS or meso?  
  @param eps_cut cut off distance numerical const */
  AppFluide_LubricationCorrection( const double &gridsize, bool is_DNS, 
		double eps_cut  );

  /** @brief Destructeur */
  ~AppFluide_LubricationCorrection();
  //@}


  /** @brief set the viscosity of fluid
  @param mu viscosity */
  void set_viscosity( double mu );
  
  
  /** @brief Calculating lubrication force
  @param time temps physique 
  @param dt pas de temps granulaire 
  @param particules liste de particules */
  virtual void CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules ); 
	
	
  /** @brief Calcul des forces lubrification
  @param p0_ Premier Composant (Particule)
  @param p1_ Deuxieme Composant (Particule ou Obstacle)
  @param contactInfos Points & Recouvrement entre les deux Composants
  @param LC grille de cellules
  @param dt Increment de temps du pas de calcul.
 */
  void computeforce( Composant* p0_, 
		     Composant* p1_,
		     LinkedCell *LC,
		     Scalar dt ) ;	

protected:


private:
  /**@name Methods */
  //@{
  /** @brief Fill the resistance matrix 
  @param RM resistance matrix
  @param eps non-dimensional gap
  @param PPorPW 0 if particle-particle and 1 if particle-wall case
  */
  void fill_ResistanceMatrix(struct ResistanceMatrix *RM, Scalar eps,
  	 int PPorPW );
	 
  /**@name Methods */
  //@{
  /** @brief Fill the resistance matrix 
  @param RM resistance matrix
  @param eps non-dimensional gap
  @param PPorPW 0 if particle-particle and 1 if particle-wall case
  */
  void nullify_ResistanceMatrix(struct ResistanceMatrix *RM, Scalar eps,
  	 int PPorPW );
	 
  /**@name Parameters */
  //@{
  list< vector<Composant*> > m_pair; /**< liste de paires de composants sur
  	lesquels la force de lubrication est appliquée */  
  double m_gridsize; /**< taille de maille du maillage fluide */
  double m_eps_cut; /**< dimmless cut-off distance to saturate the lubrication force  */
  bool m_isDNS; /** DEMCFD or DNS? */
  LinkedCell const* m_LC; /**< lien vers la grille de cellules (LinkedCell) */
  double m_viscosity; /**< viscosite du fluide */ 
  //@} 
};

#endif
