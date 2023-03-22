#ifndef REG_PROJ_NAVIER_STOKES_HH
#define REG_PROJ_NAVIER_STOKES_HH

#include <FV_OneStepIteration.hh>
#include <geomVector.hh>
#include <PAC_computingtime.hh>
#include <PAC_solvercomputingtime.hh>
#include <REG_ProjectionNavierStokesSystem.hh>
#include <FV_DiscreteField.hh>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;

class MAC_Communicator ;
class MAC_ListIdentity ;


/** @brief The Class REG_ProjectionNavierStokes.

Solver for the incompressible Newtonian NavierStokes problem. 

@author A. Wachs - Pacific project 2018-2019 */

class REG_ProjectionNavierStokes : public FV_OneStepIteration, 
public PAC_ComputingTime, public PAC_SolverComputingTime
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      /** @name Substeps of the step by step progression */
      //@{
      /** @brief Tasks performed at initialization of the algorithm, before
      starting the time stepping loop
      @param t_it time iterator */
      virtual void do_before_time_stepping( FV_TimeIterator const* t_it, 
      	string const& basename ) ;
      
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
	
      // Save other data than FV fields for restart 
      virtual void do_additional_save_for_restart( 
      	FV_TimeIterator const* t_it,
      	size_t const& restartCycleNumber, string const& basename ) ;
      //@}      


   //-- Persistence

      // Extend `list' so that it contains all objects required by the
      // storage and retrieval mechanisms that are not part of the
      // `FV_OneStepIteration::' base class subobject.
      virtual void add_storable_objects( MAC_ListIdentity* list ) ;
            

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Destructor */          
      ~REG_ProjectionNavierStokes( void ) ;
     
      /** @brief Copy constructor */      
      REG_ProjectionNavierStokes( REG_ProjectionNavierStokes const& other ) ;
      
      /** @brief Operator == 
      @param other the right hand side */        
      REG_ProjectionNavierStokes& operator=( 
      	REG_ProjectionNavierStokes const& other ) ;
      
      /** @brief Constructor with arguments 
      @param a_owner the PEL-based object   
      @param dom mesh and fields
      @param exp to read the data file */                 
      REG_ProjectionNavierStokes( MAC_Object* a_owner, 
      		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) ;

      /** @brief Constructor without argument */      
      REG_ProjectionNavierStokes( void ) ;

      /** @brief Create a clone
      @param a_owner the PEL-based object
      @param dom mesh and fields
      @param exp to read the data file */
      virtual REG_ProjectionNavierStokes* create_replica( 
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const ;
      //@}


   //-- Basic discrete system building
   
      /** @name Basic discrete system building */
      //@{
      /** @brief Assemble pressure boundary conditions in
      Navier-Stokes momentum equation vector
      @param VEC_rhs distributed vector */
      void assemble_pressure_DirichletBC_in_momentumEquation( 
      		LA_Vector* VEC_rhs ) ;  
		
      /** @brief Assemble periodic pressure gradient in Navier-Stokes momentum 
      equation vector for dp/dl=1
      @param VEC_rhs distributed vector */
      void assemble_unitary_periodic_pressure_gradient_rhs( 
      	LA_Vector* VEC_rhs ) ;
      //@}


   //-- Solvers

      /** @name Solvers */
      //@{		
      /** @brief Navier&Stokes Uzawa solver
      @param t_it time iterator */      
      void NavierStokes_Projection( FV_TimeIterator const* t_it );
      
      /** @brief Advection-diffusion velocity solver
      @param t_it time iterator */      
      void NavierStokes_AdvectionDiffusion_PredictionStep( 
      	FV_TimeIterator const* t_it );
      
      /** @brief Solver that projects velocity on a divergence free space by
      solving a pressure Poisson problem, and update velocity and pressure
      @param t_it time iterator */      
      void NavierStokes_VelocityPressure_CorrectionStep( 
      	FV_TimeIterator const* t_it );            
      //@}      

      
   //-- Utilities

      /** @name Utilities */
      //@{		            
      /** @ brief Reload other data than FV fields for restart  */
      void do_additional_reload( string const& basename ) ;
      
      /** @ brief Update pressure drop in case of periodic imposed flow rate 
      @param t_it time iterator */
      void update_pressure_drop_imposed_flow_rate( 
      	FV_TimeIterator const* t_it ) ;
	
      /** @ brief Compute and print L2 and Linf norms of div(u) */
      void compute_and_print_divu_norm( void ) ;	            
      //@}       


   private: //----------------------------------------------------------------
      
   //-- Class attributes

      static REG_ProjectionNavierStokes const* PROTOTYPE ;

   //-- Attributes

      FV_DiscreteField* UU;
      FV_DiscreteField* PP;              
      
      REG_ProjectionNavierStokesSystem* GLOBAL_EQ ;      

      // MPI parameters
      size_t nb_ranks;
      size_t my_rank;
      size_t is_master;
      MAC_Communicator const* macCOMM;
      
      // Physical Parameters
      size_t dim;
      double density;
      double viscosity;
      
      // Numerical parameters
      double imposed_CFL ;
      string AdvectionScheme ;
      string resultsDirectory;
      size_t ViscousTimeAccuracy ;
      size_t AdvectionTimeAccuracy ;
      bool b_ExplicitPressureGradient ;
      bool b_HighOrderPressureCorrection ;
      size_t sub_prob_number;
      
      // Restart
      bool b_restart ;

      // Pressure rescaling in case of all non-Dirichlet BCs
      bool b_pressure_rescaling;                  
} ;

#endif
