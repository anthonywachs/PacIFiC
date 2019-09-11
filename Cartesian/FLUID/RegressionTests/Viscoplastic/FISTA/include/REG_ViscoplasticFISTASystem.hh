#ifndef REG_VISCOPLASTIC_FISTA_SYSTEM_HH
#define REG_VISCOPLASTIC_FISTA_SYSTEM_HH

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


/** @brief The Class REG_ViscoplasticFISTASystem.

Generic matrix systems for the resolution of viscoplastic NavierStokes equations 

@author A. Wachs - Pacific project 2018-2019 */

class REG_ViscoplasticFISTASystem : public MAC_Object
{
   private: //----------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */      
      REG_ViscoplasticFISTASystem( void ) ;

      /** @brief Destructor */       
      virtual ~REG_ViscoplasticFISTASystem( void ) ;	

      /** @brief Copy constructor */       
      REG_ViscoplasticFISTASystem( REG_ViscoplasticFISTASystem const& other ) ;

      /** @brief Operator == 
      @param other the right hand side */   
      REG_ViscoplasticFISTASystem& operator=( 
      	REG_ViscoplasticFISTASystem const& other ) ;
      
      /** @brief Constructor with arguments 
      @param a_owner the MAC-based object
      @param exp to read the data file 
      @param mac_uu MAC velocity field 
      @param mac_pp MAC pressure field
      @param mac_dd MAC tensor field 
      @param NS_Viscous_TimeAccuracy_ accuracy of the viscous term time
      	discretization 
      @param NS_Advection_TimeAccuracy_ accuracy of the advection term time
      	discretization 
      @param b_PreconditionedWithLapP_ precondition Uzawa with lap(P)	
      @param b_pressure_rescaling_ rescale the pressure in case of all 
      	non-Dirichlet pressure BCs */      
      REG_ViscoplasticFISTASystem ( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
        FV_DiscreteField* mac_uu,
        FV_DiscreteField* mac_pp,
	FV_DiscreteField* mac_dd,	
        const size_t &NS_Viscous_TimeAccuracy_,	
	const size_t &NS_Advection_TimeAccuracy_,
	const bool& b_PreconditionedWithLapP_, 
	const bool& b_pressure_rescaling_ );      
      //@}


   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      /** @name Instance delivery and initialization */
      //@{
      /** @brief Create and initialize an instance of 
      REG_ViscoplasticFISTASystem
      @param a_owner the MAC-based object
      @param exp to read the data file 
      @param mac_uu MAC velocity field 
      @param mac_pp MAC pressure field
      @param mac_dd MAC tensor field       
      @param NS_Viscous_TimeAccuracy_ accuracy of the viscous term time
      	discretization 
      @param NS_Advection_TimeAccuracy_ accuracy of the advection term time
      	discretization 
      @param b_PreconditionedWithLapP_ precondition Uzawa with lap(P)	
      @param b_pressure_rescaling_ rescale the pressure in case of all 
      	non-Dirichlet pressure BCs */
      static REG_ViscoplasticFISTASystem* create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
        FV_DiscreteField* mac_uu,
        FV_DiscreteField* mac_pp,
	FV_DiscreteField* mac_dd,		
        const size_t &NS_Viscous_TimeAccuracy_,	
	const size_t &NS_Advection_TimeAccuracy_,
	const bool& b_PreconditionedWithLapP_, 
	const bool& b_pressure_rescaling_ ) ;


   //-- Access

      /** @name Access */
      //@{      
      /** @brief Return the convergence criterion for Stokes (div(u)=0) 
      problem */  
      double get_Stokes_convergence_criterion( void ) const ;  
      
      /** @brief Return the velocity solution vector */  
      LA_SeqVector const* get_solution_velocity( void ) const ; 

      /** @brief Return the pressure solution vector */  
      LA_SeqVector const* get_solution_pressure( void ) const ;
      
      /** @brief Return the tensor unknown vector */  
      LA_Vector* get_tensor_unknown_vector( void ) ;
      
      /** @brief Return the strain rate tensor solution vector */  
      LA_SeqVector const* get_solution_strainrate_tensor( void ) const ; 
      
      /** @brief Return the lambda-r*d vector */  
      LA_Vector* get_TauhatMinus1overLd_vector( void ) ;
      
      /** @brief Return the stress tensor solution vector */  
      LA_SeqVector const* get_solution_stress_tensor( void ) const ;
      //@}


   //-- Basic operations on matrices & vectors

      /** @name Basic operations on matrices & vectors */
      //@{
      /** @brief Initialize the velocity unknown vector with field values */
      void initialize_velocity( void );

      /** @brief Initialize the pressure unknown vector with field values */
      void initialize_pressure( void );
      
      /** @brief Initialize the tau_hat rhs vector with field values */
      void initialize_tau_hat( size_t tau_hat_level );  

      /** @brief Store velocity vector at previous time step */
      void at_each_time_step( void ) ;

      /** @brief Compute velocity change from one time step to the next one */
      double compute_velocity_change( void ) ;

      /** @brief Compute velocity change from one iteration to the next one 
      and update previous iteration velocity vector */
      double compute_velocity_change_iteration( void ) ;
      
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
      @param b_with_advection add advection term */       
      void compute_velocityAdvectionDiffusion_rhs( bool const& b_restart,
      	size_t const& iteration_number, 
	bool const& b_with_advection ) ;
	
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
      
      /** @brief Assemble tensor divergence matrix in momentum equation
      @param coef mass prefactor */
      void assemble_tauGradv_tensor_divergence_matrix( double const& coef ) ;
      
      /** @brief Initialize the previous iteration velocity unknown vector 
      i.e. copy current velocity vector to previous iteration velocity vector */
      void initialize_previous_iteration_velocity( void );      

      /** @brief Compute the viscoplastic vector tau_hat-(1/L)*d, i.e.,
      "Stress tensor tau hat minus (1/L) dot strain rate tensor d"
      with field values. Note: in fact what this method does is 
      VEC_TauhatMinus1overLd =- (1/L) * VEC_TENSOR, so this assumes that
      VEC_TauhatMinus1overLd already has tau_hat and VEC_TENSOR already has d 
      @param FISTA_1overL FISTA 1/L parameter */
      void compute_viscoplastic_TauhatMinus1overLd_vector( 
      	double const& FISTA_1overL );
      //@}


   //-- Solver

      /** @name Solvers */
      //@{
      /** @brief Basic Uzawa solver */         
      bool Uzawa_solver( void );
      
      /** @brief Uzawa Navier&Stokes solver
      @param viscosity fluid viscosity
      @param density fluid density      
      @param timestep time step */      
      void Uzawa_NavierStokes_solver( const double &viscosity, 
            const double &density, const double &timestep ); 
      
      /** @brief Preconditioned Uzawa solver for the Stokes problem, where the
      preconditioner is the pressure laplacian 
      @param mupr viscosity 
      @param gamma density/time step ratio */    
      bool PreconditionedStokes_Uzawa_solver( const double &mupr, 
            const double &gamma );           
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
      FV_DiscreteField* DD;
            
      // Unknowns numbering
      FV_SystemNumbering* UU_NUM ;
      FV_SystemNumbering* PP_NUM ; 
      FV_SystemNumbering* DD_NUM ;       

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

      // Work vectors for Uzawa
      LA_Vector * VEC_q ;
      LA_Vector * VEC_t ;
      LA_Vector * VEC_r ;
      LA_Vector * VEC_w ;
      LA_Vector * VEC_x ;
      LA_Vector * VEC_s ;
      LA_Vector * VEC_z ;
      LA_Vector * VEC_mean_pressure ;
      LA_Vector * VEC_unit_pressure ;
      
      // Specific to viscoplastic FISTA
      LA_SeqVector * DD_LOC ;      
      LA_Vector * VEC_TENSOR;
      LA_Vector * VEC_U_previous_iteration;      
      LA_Vector * VEC_rhs_FullViscoplastic ;
      LA_Matrix * MAT_T_StressTensorDivergence ;
      LA_Vector * VEC_TauhatMinus1overLd ;

      // Solvers
      LA_Solver* SOLVER_A_VelocityUnsteadyPlusViscous ;
      LA_Solver* SOLVER_D_PressureLaplacian ;
      LA_Solver* SOLVER_Uzawa ;
      
      // Storing vectors for reload              
      LA_StorableVectors * VECS_Storage;      
         
      // Parameters
      double Uzawa_Stokes_precision;    
      size_t Uzawa_Stokes_maxiter;
      
      // Time accuracy
      size_t NS_Viscous_TimeAccuracy ;
      size_t NS_Advection_TimeAccuracy ;  
      
      // Precondition Uzawa with lap(P)
      bool b_PreconditionedWithLapP;
      
      // Pressure rescaling in case of all non-Dirichlet BCs
      bool b_pressure_rescaling;   
      
} ; 

#endif
