#ifndef REG_PROJ_NAVIER_STOKES_SYSTEM_HH
#define REG_PROJ_NAVIER_STOKES_SYSTEM_HH

#include <MAC_Object.hh>
#include <utility>
using namespace std;


class MAC_ModuleExplorer ;
class MAC_Timer ;
class MAC_ListIdentity ;
class size_t_vector ;
class intVector ;
class doubleVector;
class LA_Matrix ;
class LA_Vector ;
class LA_SeqVector ;
class LA_Scatter ;
class LA_Solver ;
class LA_StorableVectors ;
class MAC_TimeIterator ;
class FV_SystemNumbering ;
class FV_DiscreteField ;


/** @brief The Class REG_ProjectionNavierStokesSystem.

Generic matrix systems for the resolution of NavierStokes equations 

@author A. Wachs - Pacific project 2018-2019 */

class REG_ProjectionNavierStokesSystem : public MAC_Object
{
   private: //----------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */      
      REG_ProjectionNavierStokesSystem( void ) ;

      /** @brief Destructor */       
      virtual ~REG_ProjectionNavierStokesSystem( void ) ;	

      /** @brief Copy constructor */       
      REG_ProjectionNavierStokesSystem( 
      	REG_ProjectionNavierStokesSystem const& other ) ;

      /** @brief Operator == 
      @param other the right hand side */   
      REG_ProjectionNavierStokesSystem& operator=( 
      	REG_ProjectionNavierStokesSystem const& other ) ;
      
      /** @brief Constructor with arguments 
      @param a_owner the MAC-based object
      @param exp to read the data file 
      @param mac_uu MAC velocity field 
      @param mac_pp MAC pressure field 
      @param NS_Viscous_TimeAccuracy_ accuracy of the viscous term time
      	discretization 
      @param NS_Advection_TimeAccuracy_ accuracy of the advection term time
      	discretization 
      @param b_pressure_rescaling_ rescale the pressure in case of all 
      	non-Dirichlet pressure BCs 
      @param b_ExplicitPressureGradient_ use explicit pressure gradient in
      	advection-diffusion predictor 
      @param b_HighOrderPressureCorrection_ use high order pressure correction 
      */      
      REG_ProjectionNavierStokesSystem ( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
        FV_DiscreteField* mac_uu,
        FV_DiscreteField* mac_pp,
        size_t const& NS_Viscous_TimeAccuracy_,	
	size_t const& NS_Advection_TimeAccuracy_,
	bool const& b_pressure_rescaling_,
	bool const& b_ExplicitPressureGradient_,
	bool const& b_HighOrderPressureCorrection_ );      
      //@}


   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      /** @name Instance delivery and initialization */
      //@{
      /** @brief Create and initialize an instance of 
      REG_ProjectionNavierStokesSystem
      @param a_owner the MAC-based object
      @param exp to read the data file 
      @param mac_uu MAC velocity field 
      @param mac_pp MAC pressure field 
      @param NS_Viscous_TimeAccuracy_ accuracy of the viscous term time
      	discretization 
      @param NS_Advection_TimeAccuracy_ accuracy of the advection term time
      	discretization 
      @param b_pressure_rescaling_ rescale the pressure in case of all 
      	non-Dirichlet pressure BCs 
      @param b_ExplicitPressureGradient_ use explicit pressure gradient in
      	advection-diffusion predictor 
      @param b_HighOrderPressureCorrection_ use high order pressure correction 
      */
      static REG_ProjectionNavierStokesSystem* create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
        FV_DiscreteField* mac_uu,
        FV_DiscreteField* mac_pp,
        size_t const& NS_Viscous_TimeAccuracy_,	
	size_t const& NS_Advection_TimeAccuracy_,
	bool const& b_pressure_rescaling_, 
	bool const& b_ExplicitPressureGradient_,
	bool const& b_HighOrderPressureCorrection_ ) ;


   //-- Access

      /** @name Access */
      //@{      
      /** @brief Return the velocity solution vector */  
      LA_SeqVector const* get_solution_velocity( void ) const ; 

      /** @brief Return the pressure solution vector */  
      LA_SeqVector const* get_solution_pressure( void ) const ;
      
      /** @brief Return the pressure Dirichlet BC vector */  
      LA_Vector* get_pressure_DirichletBC_vector( void ) ;
      
      /** @brief Return the periodic pressure drop vector */  
      LA_Vector* get_unitary_periodic_pressure_drop_vector( void ) ;            
      //@}


   //-- Basic operations on matrices & vectors

      /** @name Basic operations on matrices & vectors */
      //@{
      /** @brief Initialize the velocity unknown vector with field values */
      void initialize_velocity( void );

      /** @brief Initialize the pressure unknown vector with field values */
      void initialize_pressure( void );  

      /** @brief Store velocity vector at previous time step */
      void at_each_time_step( void ) ;

      /** @brief Compute velocity change from one time step to the next one */
      double compute_velocity_change( void ) ;
      
      /** @brief Compute velocity divergence norm */
      double compute_velocity_divergence_norm( void ) ;
      
      /** @brief Nullify velocity advection (inertia) right hand side */
      void nullify_velocity_advection_rhs( void );       
      
      /** @brief Assemble the velocity advection rhs term coef*u*grad(u)
      @param AdvectionScheme advection scheme ("Upwind" or "TVD")
      @param advecting_level level of advecting velocity
      @param coef prefactor
      @param advected_level level of advected velocity */          
      void assemble_velocity_advection( string const& AdvectionScheme,
      	size_t advecting_level, double const& coef, size_t advected_level ) ;
      
      /** @brief Finalize constant matrices */
      void finalize_constant_matrices( void ) ;
      
      /** @brief Compute second order velocity advection-diffusion rhs 
      @param restart if the simulation is a restart
      @param iteration_number iteration number
      @param b_with_advection add advection term 
      @param dpdl periodic pressure gradient source term (is 0 if no periodic
      pressure drop is imposed) */       
      void compute_velocityAdvectionDiffusion_rhs( bool const& b_restart,
      	size_t const& iteration_number, 
	bool const& b_with_advection,
	double const& dpdl ) ;
	
      /** @brief Assemble velocity viscous matrix and rhs
      @param coef_lap laplacian coefficient */
      void assemble_velocity_viscous_matrix_rhs( double const& coef_lap ) ;	
      
      /** @brief Assemble velocity unsteady matrix 
      @param coef_lap mass coefficient */
      void assemble_velocity_unsteady_matrix( double const& coef ) ;
      
      /** @brief Assemble velocity divergence matrix and rhs
      @param coef coefficient (either -1 or 1 ) */
      void assemble_pdivv_matrix_rhs( double const& coef ) ;
      
      /** @brief Assemble pressure laplacian matrix and rhs
      @param coef_lap laplacian coefficient */
      void assemble_pressure_laplacian_matrix_rhs( double const& coef_lap ) ;
	
      /** @brief Correct the pressure laplacian operator in case of all non
      Dirichlet BCs i.e. add coefdiag to one diagonal coefficient to ensure that
      pressure is 0 at this node */ 
      void pressure_laplacian_correction( void ) ;
      
      /** @brief Store u.grad(u)^(n-2) at first sub-iteration when CFL condition
      is not met and 2nd order scheme degenerates to 1st order */
      void store_ugradu_Nm2( size_t const& n_advection_subtimesteps );      	
      //@}


   //-- Solver

      /** @name Solvers */
      //@{
      /** @brief Velocity Diffusion (viscous) solver with velocity advection 
      in the rhs */ 	
      bool VelocityDiffusion_solver( void ) ;       

      /** @brief Velocity-pressure correction solver: (i) solve a pressure
      Poisson problem and (ii) update velocity and pressure 
      @param density fluid density 
      @param viscosity fluid viscosity 
      @param timestep time step magnitude */ 	
      double VelocityPressure_correction_solver( 
      	double const& density, double const& viscosity, 
	double const& timestep ) ;
	
      /** @brief Velocity Advection solver */ 	
      bool VelocityAdvection_solver( void ) ; 	
      //@}

 
   //-- Input / Output methods

      /** @name Input / Output methods */
      //@{               
      //@}  
 
            
   //-- Persistence

      // Add the vectors to be stored for restart
      void add_storable_objects( MAC_ListIdentity* list ) const ;
      
            
   protected: //--------------------------------------------------------


   private: //----------------------------------------------------------      
            
   //-- Initialize matrices & vectors 

      /** @name Initialize matrices & vectors */
      //@{
      /** @brief Create matrices & vectors (without allocating memory)
      @param exp to read the data file */       
      void build_system( MAC_ModuleExplorer const* exp ) ;

      /** @brief Allocate memory for matrices & vectors */  
      void re_initialize( void ) ;
      //@} 
      

   //-- Attributes

      // Fields
      FV_DiscreteField* UU ; 
      FV_DiscreteField* PP ;                              
      
      // Unknowns numbering
      FV_SystemNumbering* UU_NUM ;
      FV_SystemNumbering* PP_NUM ; 

      // Velocity Unsteady Matrix
      LA_Matrix * MAT_A_VelocityUnsteady ;
      LA_Vector * VEC_rhs_A_Velocity ;

      // Velocity Unsteady + Viscous Matrix for Stokes sub-problem
      LA_Matrix * MAT_A_VelocityUnsteadyPlusViscous ;
      LA_Vector * VEC_rhs_A_VelocityViscous ;

      // Velocity Divergence B=p.Div(v) Matrix
      LA_Matrix * MAT_B_VelocityDivergence ;
      LA_Vector * VEC_rhs_B_VelocityDivergence ;

      // Pressure gradient rhs
      LA_Vector * VEC_rhs_Bt_PressureGradient ;

      // Pressure Laplacian preconditioner
      LA_Matrix * MAT_D_PressureLaplacian ;
      LA_Vector * VEC_rhs_D_PressureLaplacian ;

      // Local vectors
      LA_SeqVector * U_LOC ;
      LA_SeqVector * P_LOC ;

      // Unknowns vectors
      LA_Vector * VEC_U;
      LA_Vector * VEC_U_previoustime;
      LA_Vector * VEC_U_timechange;
      LA_Vector * VEC_P;  
      
      // Velocity advection rhs
      LA_Vector * VEC_rhs_VelocityAdvection ;
      LA_Vector * VEC_rhs_VelocityAdvectionDiffusion ;
      LA_Vector * VEC_rhs_VelocityAdvection_Nm2 ;

      // Periodic pressure drop rhs for dp/dl = 1
      LA_Vector * VEC_rhs_A_UnitaryPeriodicPressureGradient ; 

      // Work vectors 
      LA_Vector * VEC_r ;
      LA_Vector * VEC_w ;
      LA_Vector * VEC_mean_pressure ;
      LA_Vector * VEC_unit_pressure ;
      LA_Vector * VEC_inverse_mean_pressure ; 

      // Solvers
      LA_Solver* SOLVER_A_VelocityUnsteadyPlusViscous ;
      LA_Solver* SOLVER_D_PressureLaplacian ;
      LA_Solver* SOLVER_A_VelocityUnsteady ;      
      
      // Storing vectors for reload              
      LA_StorableVectors * VECS_Storage;               
      
      // Time accuracy
      size_t NS_Viscous_TimeAccuracy ;
      size_t NS_Advection_TimeAccuracy ;        
      
      // Pressure rescaling in case of all non-Dirichlet BCs
      bool b_pressure_rescaling; 
      
      // Features of the fractional step projection algorithm
      bool b_ExplicitPressureGradient ;
      bool b_HighOrderPressureCorrection ;  
      
} ; 

#endif
