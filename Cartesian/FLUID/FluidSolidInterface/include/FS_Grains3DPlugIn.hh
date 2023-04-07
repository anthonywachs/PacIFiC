#ifndef _FS_GRAINS3DPLUGIN__
#define _FS_GRAINS3DPLUGIN__

#include <FS_SolidPlugIn.hh>
class GrainsCoupledWithFluid;
class MAC_Communicator;


/** @brief The Class FS_Grains3DPlugIn.hh.

The Grains3D plug-in to handles the Lagrangian tracking of rigid bodies with 
collisions.

@author A. Wachs - Pacific project 2021 */

class FS_Grains3DPlugIn : public FS_SolidPlugIn
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor with arguments
      @param insertion_file_ insertion file name
      @param simulation_file_ simulation file name
      @param fluid_density fluid density
      @param correct_particle_acceleration particle acceleration is corrected by
      the factor ( 1 - fluid_density / particle_density ) if value is true
      @param b_restart is the run a restart or a new run
      @param grid_size size of the smallest grid cell
      @param is_solidsolver_parallel is Grains3D running in parallel ?
      @param error =0 if the construction is successful */
      FS_Grains3DPlugIn( string const& insertion_file_,
        string const& simulation_file_,
        double const& fluid_density,
	bool const& correct_particle_acceleration,
        bool const& b_restart,
        double const& grid_size,
        bool const& is_solidsolver_parallel,
	int& error );

      /** @brief Destructor */
      ~FS_Grains3DPlugIn();
      //@}


   //-- Methods

      /** @name Methods */
      //@{
      /** @brief Simulation
      @param time_interval fluid time step
      @param predictor if yes, predictor phase, otherwise corrector phase
      @param isPredictorCorrector is the coupling scheme predictor-corrector
      @param contact_force_coef contact forces coefficient
      @param explicit_added_mass whether to treat added mass (and torque) term
        explicitly */
      void Simulation( double const& time_interval,
      	bool const& predictor = true,
        bool const& isPredictorCorrector = false,
        double const& contact_force_coef = 1.,
        bool const& explicit_added_mass = false );

//       /** @brief Write solid body features to be sent as a string vector to the 
//       fluid 
//       @param solidbodyfeatures vector of strings containing solid body 
//       features */
//       void getSolidBodyFeatures( vector<string>& solidbodyfeatures ); 
    
      /** @brief Writes solid body features to be sent as a stream to the fluid 
      @param is the stream */
      void getSolidBodyFeatures( istringstream* & is );	

      /** @brief Saves results
      @param filename file name
      @param time physical time
      @param cycleNumber cycle number */
      void saveResults( string const& filename, double const& time,
        int const& cycleNumber );
	
      /** @brief Transfer hydro force and torque array to solid solver
      @param hydroFT hydro force and torque array */
      void transferHydroFTtoSolid( 
      	vector< vector<double> > const* hydroFT ) const;
	
      /** @brief Check that Paraview writer is activated
      @param solid_resDir particles results directory */
      void checkParaviewPostProcessing( string const& solid_resDir );	
	
      /** @brief Set the post-processing translation vector in case of
      projection-translation
      @param tvx x coordinate
      @param tvy y coordinate
      @param tvz z coordinate */
      void setParaviewPostProcessingTranslationVector( 
      	double const& tvx, double const& tvy, double const& tvz );
      //@}
    	  
  
   protected: //--------------------------------------------------------------


   private: //----------------------------------------------------------------

   //-- Attributes  

      /**@name Parameters */
      //@{
      GrainsCoupledWithFluid* m_Grains3D; /**< pointer to an instance of 
    	Grains3D */
      MAC_Communicator const* m_macCOMM; /**< MPI communicator */       
      size_t m_nb_ranks; /**< Total number of processes */ 
      size_t m_my_rank; /**< Rank of process */ 
      size_t m_is_master; /**< Rank of the master process (set to 0) */ 
      bool m_Grains3D_parallel_mode; /**< whether Grains3D runs in parallel */
      bool m_Grains3D_active_on_this_rank; /**< whether Grains3D 
      	is active on this process */
      //@}


   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      FS_Grains3DPlugIn();    
    
      /** @brief Copy constructor 
      @param copy copied FS_Grains3DPlugIn object */
      FS_Grains3DPlugIn( FS_Grains3DPlugIn const& copy );
      //@}
};

#endif
