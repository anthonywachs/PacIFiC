#ifndef DS_DirectionSplitting_HH
#define DS_DirectionSplitting_HH

#include <mpi.h>
#include <FV_OneStepIteration.hh>
#include <computingtime.hh>
#include <solvercomputingtime.hh>
#include <DS_NavierStokes.hh>
#include <DS_HeatTransfer.hh>
#include <MAC_DoubleVector.hh>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;

class MAC_Communicator ;
class FV_DiscreteField ;
class FS_SolidPlugIn ;
class DS_AllRigidBodies ;

/** @brief The Class DS_DirectionSplitting.

Server for the intiating the NavierStokes and/or HeatTransfer classes.

@author A. Goyal - Pacific project 2022 */

class DS_DirectionSplitting : public FV_OneStepIteration,
                           public ComputingTime,
                           public SolverComputingTime
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      /** @name Substeps of the step by step progression */
      //@{
      /** @brief Tasks performed at initialization of the algorithm, before
      starting the time stepping loop
      @param t_it time iterator */
      virtual void do_before_time_stepping( FV_TimeIterator const* t_it,
      	std::string const& basename ) ;

      /** @brief Perform one time step
      @param t_it time iterator */
      virtual void do_one_inner_iteration( FV_TimeIterator const* t_it ) ;

      /** @brief Tasks performed at initialization of each time step
      @param t_it time iterator */
      virtual void do_before_inner_iterations_stage(
      	FV_TimeIterator const* t_it );

      /** @brief Tasks performed after of each time step
      @param t_it time iterator */
      virtual void do_after_inner_iterations_stage(
      	FV_TimeIterator const* t_it );

      /** @brief Tasks performed at the end of the time stepping loop */
      virtual void do_after_time_stepping( void );

      /** @brief Save additional data than fields
      @param t_it time iterator
      @param cycleNumber cycle number */
      virtual void do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber  );

      /** @brief Inbuilt function for post processing data
      @param dom domain and fields
      @param exp module explorer */
      virtual void do_more_post_processing( FV_DomainAndFields * dom,
                                           MAC_ModuleExplorer const* exp ) ;
      //@}

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Destructor */
      ~DS_DirectionSplitting( void ) ;

      /** @brief Copy constructor */
      DS_DirectionSplitting( DS_DirectionSplitting const& other ) ;

      /** @brief Operator ==
      @param other the right hand side */
      DS_DirectionSplitting& operator=( DS_DirectionSplitting const& other ) ;

      /** @brief Constructor with arguments
      @param a_owner the MAC-based object
      @param exp to read the data file */
      DS_DirectionSplitting( MAC_Object* a_owner,
      		              FV_DomainAndFields const* dom,
	                       MAC_ModuleExplorer const* exp ) ;

      /** @brief Constructor without argument */
      DS_DirectionSplitting( void ) ;

      /** @brief Create a clone
      @param a_owner the MAC-based object
      @param dom mesh and fields
      @param prms set of parameters
      @param exp to read the data file */
      virtual DS_DirectionSplitting* create_replica(
                     		MAC_Object* a_owner,
                     		FV_DomainAndFields const* dom,
                     		MAC_ModuleExplorer* exp ) const ;
      //@}


   private: //----------------------------------------------------------------

   //-- Class attributes

      static DS_DirectionSplitting const* PROTOTYPE ;

   //-- Attributes

      double rho;
      double mu;
      double kai;
      string AdvectionScheme;
      size_t AdvectionTimeAccuracy;
      size_t space_dimensions;
      bool b_restart;
      bool is_solids;
      bool is_HE, is_NS, is_NSwithHE;
      double RBTemp;

      string insertion_type;
      bool is_stressCal;
      string ViscousStressOrder;
      double surface_cell_scale;
      bool is_surfacestressOUT;
      size_t stressCalFreq;
      bool is_par_motion;
      MAC_DoubleVector* gravity_vector ;

      DS_NavierStokes* FlowSolver ;
      DS_HeatTransfer* HeatSolver ;

      // Grains3D variable
      string solidSolverType;
      FS_SolidPlugIn* solidSolver;
      bool b_solidSolver_parallel;
      string solidSolver_insertionFile;
      string solidSolver_simulationFile;
      istringstream* solidFluid_transferStream;
      DS_AllRigidBodies* allrigidbodies;
      bool b_particles_as_fixed_obstacles;

      double critical_distance_translation;

      MAC_Communicator const* macCOMM;

} ;

#endif
