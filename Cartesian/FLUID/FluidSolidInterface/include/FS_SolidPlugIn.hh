#ifndef _FS_SOLIDPLUGIN__
#define _FS_SOLIDPLUGIN__

#include <geomVector.hh>
#include <MAC_assertions.hh>
#include <PAC_solvercomputingtime.hh>
#include <string>
#include <sstream>
using std::string;
using std::istringstream ;


/** @brief The Class FS_SolidPlugIn.

Used as a plug-in class to connect a module that handles the Lagrangian tracking
of solid bodies with potential collisions.

@author A. Wachs - Pacific project 2021 */

class FS_SolidPlugIn : public PAC_SolverComputingTime
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor with arguments
      @param insertion_file_ insertion file name
      @param simulation_file_ simulation file name */
      FS_SolidPlugIn( string const& insertion_file_,
        string const& simulation_file_ );

      /** @brief Destructor */
      virtual ~FS_SolidPlugIn();
      //@}


   //-- Virtual methods

      /** @name Virtual methods */
      //@{
      /** @brief Simulation
      @param time_interval fluid time step
      @param predictor if yes, predictor phase, otherwise corrector phase
      @param isPredictorCorrector is the coupling scheme predictor-corrector
      @param contact_force_coef contact forces coefficient
      @param explicit_added_mass whether to treat added mass (and torque) term
        explicitly */
      virtual void Simulation( double const& time_interval,
      	bool const& predictor = true,
        bool const& isPredictorCorrector = false,
        double const& contact_force_coef = 1.,
        bool const& explicit_added_mass = false ) = 0;
	
//       /** @brief Writes solid body features to be sent as a string vector to the 
//       fluid 
//       @param solidbodyfeatures vector of strings containing solid body 
//       features */
//       virtual void getSolidBodyFeatures( vector<string>& solidbodyfeatures ) 
//       	= 0; 
    
      /** @brief Writes solid body features to be sent as a stream to the fluid 
      @param is the stream */
      virtual void getSolidBodyFeatures( istringstream* & is ) = 0;	

      /** @brief Saves results
      @param filename file name
      @param time physical time
      @param cycleNumber cycle number */
      virtual void saveResults( string const& filename, double const& time,
        int const& cycleNumber ) = 0;
	
      /** @brief Transfer hydro force and torque array to solid solver
      @param hydroFT hydro force and torque array */
      virtual void transferHydroFTtoSolid( 
      	vector< vector<double> > const* hydroFT ) const = 0;
	
      /** @brief Check that Paraview writer is activated
      @param solid_resDir particles results directory */
      virtual void checkParaviewPostProcessing( string const& solid_resDir ) 
      	= 0;
	
      /** @brief Set the post-processing translation vector in case of
      projection-translation
      @param tvx x coordinate
      @param tvy y coordinate
      @param tvz z coordinate */
      virtual void setParaviewPostProcessingTranslationVector( 
      	double const& tvx, double const& tvy, double const& tvz ) = 0;
	
      /** @brief Set the initial physical time
      @param time0 initial physical time */
      virtual void setInitialTime( double const& time0 ) = 0;      	
      //@}
    

   protected: //--------------------------------------------------------------

   //-- Attributes  

      /**@name Parameters */
      //@{
      string m_insertion_file; /**< solid insertion file */
      string m_simulation_file; /**< solid simulation file */
      //@}     
    
    
   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      FS_SolidPlugIn();

      /** @brief Copy constructor
      @param copy copied FS_SolidPlugIn object */
      FS_SolidPlugIn( FS_SolidPlugIn const& copy ) ;
      //@}
};

#endif
