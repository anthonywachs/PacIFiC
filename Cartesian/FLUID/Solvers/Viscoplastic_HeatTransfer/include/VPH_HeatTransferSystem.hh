#ifndef REG_HEAT_TRANSFER_SYSTEM_HH
#define REG_HEAT_TRANSFER_SYSTEM_HH

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


/** @brief The Class VPH_HeatTransferSystem.

Generic matrix systems for the resolution of time-dependent advection diffusion 
heat transfer problem.

@author A. Wachs - Pacific project 2018-2019 */

class VPH_HeatTransferSystem : public MAC_Object
{
   private: //----------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */      
      VPH_HeatTransferSystem( void ) ;

      /** @brief Destructor */       
      virtual ~VPH_HeatTransferSystem( void ) ;	

      /** @brief Copy constructor */       
      VPH_HeatTransferSystem( 
      	VPH_HeatTransferSystem const& other ) ;

      /** @brief Operator == 
      @param other the right hand side */   
      VPH_HeatTransferSystem& operator=( 
      	VPH_HeatTransferSystem const& other ) ;
      
      /** @brief Constructor with arguments 
      @param a_owner the MAC-based object
      @param exp to read the data file 
      @param mac_tt MAC temperature field 
      @param mac_pp MAC velocity field 
      @param HEAT_Diffusion_TimeAccuracy_ accuracy of the diffusion term time
      	discretization 
      @param HEAT_Advection_TimeAccuracy_ accuracy of the advection term time
      	discretization */      
      VPH_HeatTransferSystem ( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
        FV_DiscreteField* mac_tt,
        FV_DiscreteField const* mac_uu,
        size_t const& HEAT_Diffusion_TimeAccuracy_,	
	size_t const& HEAT_Advection_TimeAccuracy_ );      
      //@}


   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      /** @name Instance delivery and initialization */
      //@{
      /** @brief Create and initialize an instance of 
      VPH_HeatTransferSystem
      @param a_owner the MAC-based object
      @param exp to read the data file 
      @param mac_tt MAC temperature field 
      @param mac_pp MAC velocity field 
      @param HEAT_Diffusion_TimeAccuracy_ accuracy of the diffusion term time
      	discretization 
      @param HEAT_Advection_TimeAccuracy_ accuracy of the advection term time
      	discretization */
      static VPH_HeatTransferSystem* create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
        FV_DiscreteField* mac_tt,
        FV_DiscreteField const* mac_uu,
        size_t const& HEAT_Diffusion_TimeAccuracy_,	
	size_t const& HEAT_Advection_TimeAccuracy_ ) ;


   //-- Access

      /** @name Access */
      //@{      
      /** @brief Return the temperature solution vector */  
      LA_SeqVector const* get_solution_temperature( void ) const ;           
      //@}


   //-- Basic operations on matrices & vectors

      /** @name Basic operations on matrices & vectors */
      //@{
      /** @brief Initialize the temperature unknown vector with field values */
      void initialize_temperature( void );

      /** @brief Store temperature vector at previous time step */
      void at_each_time_step( void ) ; 
      
      /** @brief Compute temperature change from one time step to the next one 
      */
      double compute_temperature_change( void ) ;
      
      /** @brief Nullify temperature advection right hand side */
      void nullify_temperature_advection_rhs( void );
      
      /** @brief Assemble the temperature advection rhs term coef*u*grad(T)
      @param AdvectionScheme advection scheme ("Upwind" or "TVD")
      @param advecting_level level of advecting velocity
      @param coef prefactor
      @param advected_level level of advected velocity */          
      void assemble_temperature_advection( string const& AdvectionScheme,
      	size_t advecting_level, double const& coef, size_t advected_level ) ;

      /** @brief Finalize constant matrices */
      void finalize_constant_matrices( void ) ;
      
      /** @brief Compute second order temperature advection-diffusion rhs 
      @param restart if the simulation is a restart
      @param iteration_number iteration number
      @param b_with_advection add advection term */       
      void compute_temperatureAdvectionDiffusion_rhs( bool const& b_restart,
      	size_t const& iteration_number, 
	bool const& b_with_advection ) ;
	
      /** @brief Assemble temperature diffusion matrix and rhs
      @param coef_lap laplacian coefficient */
      void assemble_temperature_diffusion_matrix_rhs( double const& coef_lap ) ;
      
      /** @brief Assemble temperature unsteady matrix 
      @param coef_lap mass coefficient */
      void assemble_temperature_unsteady_matrix( double const& coef ) ;
      
      /** @brief Store u.grad(T)^(n-2) at first sub-iteration when CFL condition
      is not met and 2nd order scheme degenerates to 1st order */
      void store_ugradT_Nm2( size_t const& n_advection_subtimesteps );      
      //@}


   //-- Solver

      /** @name Solvers */
      //@{
      /** @brief Temperature Diffusion solver with velocity advection 
      in the rhs */ 	
      bool TemperatureDiffusion_solver( void ) ;       
	
      /** @brief Temperature Advection solver */ 	
      bool TemperatureAdvection_solver( void ) ;
      //@}

 
   //-- Input / Output methods

      /** @name Input / Output methods */
      //@{               
      //@}  
 
            
   //-- Persistence

      // Add the vectors to be stored for restart
      void add_storable_objects( MAC_ListIdentity* list ) ;
      
            
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
      FV_DiscreteField* TT ; 
      FV_DiscreteField const* UU ;                              
      
      // Unknowns numbering
      FV_SystemNumbering* TT_NUM ;     

      // Temperature Unsteady Matrix
      LA_Matrix * MAT_A_TemperatureUnsteady ;
      LA_Vector * VEC_rhs_A_Temperature ;
      
      // Temperature Unsteady + Viscous Matrix for Stokes sub-problem
      LA_Matrix * MAT_A_TemperatureUnsteadyPlusDiffusion ;
      LA_Vector * VEC_rhs_A_TemperatureDiffusion ;      

      // Local vectors
      LA_SeqVector * T_LOC ;

      // Unknowns vectors
      LA_Vector * VEC_T;
      LA_Vector * VEC_T_previoustime;
      LA_Vector * VEC_T_timechange;
      
      // Temperature advection rhs
      LA_Vector * VEC_rhs_TemperatureAdvection ;
      LA_Vector * VEC_rhs_TemperatureAdvectionDiffusion ;
      LA_Vector * VEC_rhs_TemperatureAdvection_Nm2 ;

      // Solvers
      LA_Solver* SOLVER_A_TemperatureUnsteadyPlusDiffusion ;
      LA_Solver* SOLVER_A_TemperatureUnsteady ;      
      
      // Storing vectors for reload              
      LA_StorableVectors * VECS_Storage;               
      
      // Time accuracy
      size_t HEAT_Diffusion_TimeAccuracy ;
      size_t HEAT_Advection_TimeAccuracy ;        
      
} ; 

#endif
