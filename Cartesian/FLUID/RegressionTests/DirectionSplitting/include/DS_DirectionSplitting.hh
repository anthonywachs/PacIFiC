#ifndef DS_DirectionSplitting_HH
#define DS_DirectionSplitting_HH

#include <mpi.h>
#include <FV_OneStepIteration.hh>
#include <computingtime.hh>
#include <solvercomputingtime.hh>
#include <DS_NavierStokes.hh>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;

class MAC_Communicator ;
class FV_DiscreteField ;

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
      bool b_restart;
      bool is_solids;

      string insertion_type;
      bool is_stressCal;
      string ViscousStressOrder;
      double surface_cell_scale;
      bool is_surfacestressOUT;
      size_t stressCalFreq;
      bool is_par_motion;
      double grid_check_for_solid;

      DS_NavierStokes* FlowSolver ;

} ;

#endif
