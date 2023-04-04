#ifndef _FS_SOLIDPLUGIN_BUILDERFACTORY__
#define _FS_SOLIDPLUGIN_BUILDERFACTORY__

#include <string>
using std::string;
class FS_SolidPlugIn;


/** @brief The Class FS_SolidPlugIn_BuilderFactory.

Solid plug-in builder factory, creates an instance of solid plug-in using the
type of plug-in as the main parameter.

@author A. Wachs - Pacific project 2021 */

class FS_SolidPlugIn_BuilderFactory
{
   public: //-----------------------------------------------------------------

   //-- Static methods

      /** @name Static methods */
      //@{
      /** @brief Creates a solver
      @param type solid plug-in type
      @param insertion_file_ insertion file
      @param simulation_file_ simulation file
      @param fluid_density fluid density
      @param correct_particle_acceleration particle acceleration is corrected by
      the factor ( 1 - fluid_density / particle_density ) if value is true
      @param b_restart is the run a restart or a new run
      @param grid_size size of the smallest grid cell
      @param is_solidsolver_parallel is solid solver running in parallel ?
      @param error =0 if the construction is successful */
      static FS_SolidPlugIn* create( string const& type,
	string const& insertion_file_,
        string const& simulation_file_,
        double const& fluid_density,
	bool const& correct_particle_acceleration,
        bool const& b_restart,
        double const& grid_size,
        bool const& is_solidsolver_parallel,
	int& error );
      //@}	
	
	
   protected: //--------------------------------------------------------------
	

   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      FS_SolidPlugIn_BuilderFactory();

      /** @brief Destructor */
      ~FS_SolidPlugIn_BuilderFactory();
      //@}
};

#endif
