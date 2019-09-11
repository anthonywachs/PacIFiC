#ifndef REG_HEAT_EQUATION_SYSTEM_HH
#define REG_HEAT_EQUATION_SYSTEM_HH

#include <MAC_Object.hh>
#include <utility>
using namespace std;


class MAC_ModuleExplorer ;
class MAC_Timer ;
class size_t_vector ;
class intVector ;
class doubleVector;
class LA_Matrix ;
class LA_Vector ;
class LA_SeqVector ;
class LA_Scatter ;
class LA_Solver ;
class FV_SystemNumbering ;
class FV_DiscreteField ;
class FV_TimeIterator ;


/** @brief The Class REG_HeatEquationSystem.

Matrix systems for the resolution of the heat equation.

@author A. Wachs - Pacific project 2017 */

class REG_HeatEquationSystem : public MAC_Object
{
   private: //----------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */      
      REG_HeatEquationSystem( void ) ;

      /** @brief Destructor */       
      ~REG_HeatEquationSystem( void ) ;

      /** @brief Copy constructor */       
      REG_HeatEquationSystem( REG_HeatEquationSystem const& other ) ;

      /** @brief Operator == 
      @param other the right hand side */   
      REG_HeatEquationSystem& operator=( REG_HeatEquationSystem const& other ) ;

      /** @brief Constructor with arguments 
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param mac_tf FV temperature field */      
      REG_HeatEquationSystem ( MAC_Object* a_owner,
            MAC_ModuleExplorer const* exp,
            FV_DiscreteField* mac_tf );
      //@}


   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      /** @name Instance delivery and initialization */
      //@{
      /** @brief Create and initialize an instance of REG_HeatEquationSystem
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param mac_tf FV temperature field */
      static REG_HeatEquationSystem* create( 
            MAC_Object* a_owner,
            MAC_ModuleExplorer const* exp,
            FV_DiscreteField* mac_tf ) ;
      //@}


   //-- Access

      /** @name Access */
      //@{
      /** @brief Return the temperature solution vector TF */  
      LA_SeqVector const* get_solution_temperature( void ) const ;     
      //@}

   
   //-- Basic operations on matrices & vectors

      /** @name Basic operations on matrices & vectors */
      //@{
      /** @brief Initialize the temperature unknown vector with field values */
      void initialize_temperature( void );     

      /** @brief Finalize constant matrices 
      @param DiffusionTimeAccuracy 1st or 2nd order treatment of diffusion */
      void finalize_constant_matrices( size_t DiffusionTimeAccuracy ) ;
      
      /** @brief Store temperature vector at previous time step */
      void at_each_time_step( void ) ;
      
      /** @brief Compute temperature change from one time step to the 
      next one */
      double compute_temperature_change( void ); 
      
      /** @brief Assemble temperature unsteady matrix 
      @param coef_lap mass coefficient */
      void assemble_temperature_unsteady_matrix( double const& coef ) ;
      
      /** @brief Assemble temperature diffusion matrix and rhs
      @param coef_lap laplacian coefficient */
      void assemble_temperature_diffusion_matrix_rhs( double const& coef_lap ) ;
      
      /** @brief Return the buoyancy vector */  
      LA_Vector* get_diffrhs_plus_bodyterm_vector( void ) ;       
      //@}

   //-- Solver

      /** @name Solvers */
      //@{
      /** @brief Heat equation solver 
      @param DiffusionTimeAccuracy 1st or 2nd order treatment of diffusion 
      @param iter_number iteration number */ 	
      bool HeatEquation_solver( size_t DiffusionTimeAccuracy,
      	size_t iter_number ) ;
      //@}


   //-- Output methods

      /** @name Output methods */
      //@{          
      //@} 
      

   protected: //--------------------------------------------------------


   private: //----------------------------------------------------------

      /** @name Initialize matrices & vectors */
      //@{
      /** @brief Create matrices & vectors (without allocating memory)
      @param exp to read the data file */       
      void build_system( MAC_ModuleExplorer const* exp ) ;

      /** @brief Allocate memory for matrices & vectors */  
      void re_initialize( void ) ;
      //@}       

      //-- Attributes

      FV_DiscreteField* TF ;     

      // Local vectors
      LA_SeqVector * TF_LOC ;

      // Unknowns vectors
      LA_Vector * VEC_TF ; 
      LA_Vector * VEC_TF_previoustime ;
      LA_Vector * VEC_TF_timechange ;                       

      // Matrices & rhs
      LA_Matrix * MAT_D_TemperatureUnsteadyPlusDiffusion ;
      LA_Matrix * MAT_A_TemperatureUnsteady ;      
      LA_Vector * VEC_rhs_D_TemperatureDiffusionPlusBodyTerm ;
      LA_Vector * VEC_rhs_A_TemperatureUnsteady ;      
      
      // Solvers
      LA_Solver* SOLVER_Temperature ;

      // Unknowns numbering
      FV_SystemNumbering* TF_NUM ;   
} ; 

#endif
