#ifndef VPH_VISCOPLASTIC_FISTA_HH
#define VPH_VISCOPLASTIC_FISTA_HH

#include <FV_OneStepIteration.hh>
#include <geomVector.hh>
#include <PAC_computingtime.hh>
#include <PAC_solvercomputingtime.hh>
#include <VPH_ViscoplasticSystem.hh>
#include <VPH_HeatTransfer.hh>
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
class MAC_DoubleVector ;


/** @brief The Class VPH_Viscoplastic.

Solver for the incompressible viscoplastic NavierStokes problem with heat
transfer. Solution algorithm for viscoplasticity is either ALG2 or FISTA.

@author A. Wachs - Pacific project 2018-2019 */

class VPH_Viscoplastic : public FV_OneStepIteration, public PAC_ComputingTime,
	public PAC_SolverComputingTime
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
      virtual void do_additional_save_for_restart( FV_TimeIterator const* t_it,
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
      ~VPH_Viscoplastic( void ) ;
     
      /** @brief Copy constructor */      
      VPH_Viscoplastic( VPH_Viscoplastic const& other ) ;
      
      /** @brief Operator == 
      @param other the right hand side */        
      VPH_Viscoplastic& operator=( VPH_Viscoplastic const& other ) ;
      
      /** @brief Constructor with arguments 
      @param a_owner the PEL-based object   
      @param dom mesh and fields
      @param exp to read the data file */                 
      VPH_Viscoplastic( MAC_Object* a_owner, 
      		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) ;

      /** @brief Constructor without argument */      
      VPH_Viscoplastic( void ) ;

      /** @brief Create a clone
      @param a_owner the PEL-based object
      @param dom mesh and fields
      @param exp to read the data file */
      virtual VPH_Viscoplastic* create_replica( 
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const ;
      //@}


   //-- Basic discrete system building
   
      /** @name Basic discrete system building */
      //@{		      	         
      /** @brief Compute the strain rate tensorial field DD
      based on the velocity field and store it in VEC at the matrix level 
      @param VEC distributed vector */
      void compute_strain_rate_tensor_D( LA_Vector* VEC ) ;

      /** @brief Update the strain rate tensor field d using the Lagrange 
      multiplier tensor and the strain rate tensor D(u) 
      @param VEC distributed vector
      @param alpha multiplication factor of the whole term
      @param beta multiplication factor of D(u), beta is a scalar 
      @param first_field_level level of first field in 
      first_field + beta*D(u), second field is always D(u) */
      void update_strain_rate_tensor_d( LA_Vector* VEC, double const& alpha, 
      	double const& beta, size_t const& first_field_level ) ;
	
      /** @brief Update the stress tensor field tau
      using the strain rate tensors D(u) and d and the stress tensor tau_hat
      and store it in VEC at the matrix level 
      @param VEC distributed vector */
      void update_tau( LA_Vector* VEC  ) ;
      
      /** @brief Update the stress tensor field tau_hat
      using the strain stress tensor tau_hat, the current value of the fista
      parameter and the old value of the fista parameter
      and store it in VEC at the matrix level 
      @param VEC distributed vector 
      @param stepsize step size in update formula */
      void update_tau_hat( LA_Vector* VEC, double const& stepsize  ) ;      

      /** @brief Update the Lagrange multiplier field lambda
      using the strain rate tensors D(u) and d 
      and store it in VEC at the matrix level 
      @param VEC distributed vector */
      void update_Lagrange_multiplier( LA_Vector* VEC  ) ;
      
      /** @brief Compute the norm of D-d */
      double compute_norm_Dminusd( void ) ;
      	
      /** @brief Assemble buoyancy in Navier-Stokes momentum equation vector */
      void assemble_buoyancy_rhs( LA_Vector* VEC_rhs ) ;            	
      //@}


   //-- Solvers

      /** @name Solvers */
      //@{		
      /** @brief FISTA viscoplastic solver
      @param t_it time iterator */      
      void FISTA_solver( FV_TimeIterator const* t_it );
      
      /** @brief ALG2 viscoplastic solver
      @param t_it time iterator */      
      void ALG2_solver( FV_TimeIterator const* t_it );      
      //@}      

      
   //-- Utilities

      /** @name Utilities */
      //@{		      
      /** @ brief Reload other data than FV fields for restart  */
      void do_additional_reload( string const& basename ) ; 
      
      /** @brief Return the 2D Euclidian norm of first_field+beta.D(u)
      stored as a vector 
      @param beta multiplication factor of D(u), beta is a scalar 
      @param first_field_level level of first field */
      double Compute_euclidian_norm(
      		size_t i, size_t j, size_t comp,
		FV_SHIFT_TRIPLET shift, double const& beta, 
		size_t const& first_field_level ) ;

      /** @brief Return the 2D Euclidian norm of first_field+beta.D(u)
      stored as a vector 
      @param beta multiplication factor of D(u), beta is a scalar 
      @param first_field_level level of first field */
      double Compute_euclidian_norm(
      		size_t i, size_t j, size_t k, size_t comp,
		FV_SHIFT_TRIPLET shift, double const& beta, 
		size_t const& first_field_level ) ;            
      //@}       


   private: //----------------------------------------------------------------
      
   //-- Class attributes

      static VPH_Viscoplastic const* PROTOTYPE ;

   //-- Attributes

      FV_DiscreteField* UU;
      FV_DiscreteField* PP;
      FV_DiscreteField* DD;                    
      
      VPH_ViscoplasticSystem* GLOBAL_EQ ;      

      // MPI parameters
      size_t nb_ranks;
      size_t my_rank;
      size_t is_master;
      MAC_Communicator const* macCOMM;
      
      // Physical Parameters
      size_t dim;
      double density;
      double viscosity;
      double yield_stress;
      
      // Numerical parameters
      double imposed_CFL ;
      string AdvectionScheme ;
      string resultsDirectory;
      size_t ViscousTimeAccuracy ;
      size_t AdvectionTimeAccuracy ;
      bool b_PreconditionedWithLapP;
      string ViscoplasticSolver;
      double FISTA_1overL_param;
      double VP_tolerance;
      size_t VP_maxiter;
      size_t d_level ;
      size_t D_level ;
      size_t tau_hat_level ;
      size_t tau_level ;
      size_t tau_old_level ;
      double ALG2_aug_param;
      size_t lambda_level ;
      
      // Restart
      bool b_restart ;

      // Pressure rescaling in case of all non-Dirichlet BCs
      bool b_pressure_rescaling;
      
      // Temperature
      FV_DiscreteField const* TT; 
      bool b_with_buoyancy ; 
      double boussinesq_thermal_expansion_coef ;
      double boussinesq_reference_temperature ;
      MAC_DoubleVector* gravity_vector ;
      VPH_HeatTransfer* Solver_Temperature ;  
                              
} ;

#endif
