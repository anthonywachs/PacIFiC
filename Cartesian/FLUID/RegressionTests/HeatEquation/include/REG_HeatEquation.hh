#ifndef REG_HEATEQUATION_HH
#define REG_HEATEQUATION_HH

#include <FV_OneStepIteration.hh>
#include <geomVector.hh>
#include <computingtime.hh>
#include <solvercomputingtime.hh>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;

class MAC_Communicator ;
class FV_DiscreteField ;
class LA_Vector ;
class REG_HeatEquationSystem ;

/** @brief The Class REG_HeatEquation.

Server for the resolution of the unsteady heat equation by a first order
implicit time integrator and a Finite Volume MAC scheme on rectangular grids.

Equation: dT/dt = ( 1 / Pe ) * lap(T) + bodyterm, where Pe is the Peclet number.

@author A. Wachs - Pacific project 2017 */

class REG_HeatEquation : public FV_OneStepIteration, public ComputingTime,
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
      ~REG_HeatEquation( void ) ;
     
      /** @brief Copy constructor */      
      REG_HeatEquation( REG_HeatEquation const& other ) ;
      
      /** @brief Operator == 
      @param other the right hand side */        
      REG_HeatEquation& operator=( REG_HeatEquation const& other ) ;
      
      /** @brief Constructor with arguments 
      @param a_owner the MAC-based object   
      @param exp to read the data file */                 
      REG_HeatEquation( MAC_Object* a_owner, 
      		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) ;

      /** @brief Constructor without argument */      
      REG_HeatEquation( void ) ;

      /** @brief Create a clone
      @param a_owner the MAC-based object
      @param dom mesh and fields
      @param prms set of parameters      
      @param exp to read the data file */
      virtual REG_HeatEquation* create_replica( 
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const ;
      //@}


   //-- Basic discrete system building
   
      /** @name Basic discrete system building */
      //@{	
      /** @brief Assemble temperature body term */
      void assemble_temperature_bodyterm_rhs(  
      	FV_DiscreteField const* FF,
	LA_Vector* VEC_rhs ) ;	
      
      /** @brief Error compared to analytical solution */ 
      void error_with_analytical_solution ( FV_DiscreteField const* FF,
      	 FV_DiscreteField* FF_ERROR ) ;
	 
      /** @brief Compute L2 norm of solution */ 
      double compute_L2Norm_solution ( FV_DiscreteField const* FF,
      	size_t const& comp ) ;
      //@}

      
   private: //----------------------------------------------------------------
      
   //-- Class attributes

      static REG_HeatEquation const* PROTOTYPE ;

   //-- Attributes

      FV_DiscreteField* TF;
      FV_DiscreteField* TF_STAG;
      FV_DiscreteField* TF_VERTEX;       
      
      FV_DiscreteField* TF_ERROR;
      FV_DiscreteField* TF_STAG_ERROR;
      FV_DiscreteField* TF_VERTEX_ERROR;                  

      REG_HeatEquationSystem* GLOBAL_EQ ;
      REG_HeatEquationSystem* GLOBAL_EQ_STAG ; 
      REG_HeatEquationSystem* GLOBAL_EQ_VERTEX ;            

      size_t nb_procs;
      size_t my_rank;
      size_t is_master;
      size_t dim;
      MAC_Communicator const* pelCOMM;
      
      double peclet ; 
      bool b_bodyterm ;
      bool b_restart ;
      size_t DiffusionTimeAccuracy ;                                     
} ;

#endif
