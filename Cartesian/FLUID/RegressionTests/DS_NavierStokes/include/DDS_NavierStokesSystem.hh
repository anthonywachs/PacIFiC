#ifndef REG_HEAT_EQUATION_SYSTEM_HH
#define REG_HEAT_EQUATION_SYSTEM_HH

#include <MAC_Object.hh>
#include <utility>
#include <boolVector.hh>
using namespace std;


class MAC_ModuleExplorer ;
class MAC_Communicator ;
class MAC_Timer ;
class size_t_vector ;
class intVector ;
class doubleVector;
class LA_Matrix ;
class LA_SeqMatrix ;
class LA_Vector ;
class LA_SeqVector ;
class LA_Scatter ;
class LA_Solver ;
class LA_CRSmatrix ;
class FV_SystemNumbering ;
class FV_DiscreteField ;
class FV_TimeIterator ;


/** @brief The Class DDS_NavierStokesSystem.

Matrix systems for the resolution of the heat equation.

@author A. Wachs - Pacific project 2017 */

class DDS_NavierStokesSystem : public MAC_Object
{
   private: //----------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */
      DDS_NavierStokesSystem( void ) ;

      /** @brief Destructor */
      ~DDS_NavierStokesSystem( void ) ;

      /** @brief Copy constructor */
      DDS_NavierStokesSystem( DDS_NavierStokesSystem const& other ) ;

      /** @brief Operator ==
      @param other the right hand side */
      DDS_NavierStokesSystem& operator=( DDS_NavierStokesSystem const& other ) ;

      /** @brief Constructor with arguments
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param mac_UF FV velocity field */
      DDS_NavierStokesSystem ( MAC_Object* a_owner,
            MAC_ModuleExplorer const* exp,
            FV_DiscreteField* mac_UF,
            FV_DiscreteField* mac_PF );
      //@}


   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      /** @name Instance delivery and initialization */
      //@{
      /** @brief Create and initialize an instance of DDS_NavierStokesSystem
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param mac_UF FV velocity field */
      static DDS_NavierStokesSystem* create(
            MAC_Object* a_owner,
            MAC_ModuleExplorer const* exp,
            FV_DiscreteField* mac_UF,
            FV_DiscreteField* mac_PF ) ;
      //@}


   //-- Access

      /** @name Access */
      //@{
      /** @brief Return the velocity solution vector UF */
      LA_SeqVector const* get_solution_velocity( void ) const ;

      /** @brief Return the DS velocity solution vector DS_UF */
      LA_SeqVector const* get_solution_DS_velocity( void ) const ;

      /** @brief Return the DS pressure solution vector DS_PF */
      LA_SeqVector const* get_solution_DS_velocity_P( void ) const ;

      /** @brief Return the DS rhs in x vector */
      LA_SeqVector* get_rhs_DS_velocity_x( size_t const& comp ) ;

      /** @brief Return the DS rhs in y vector */
      LA_SeqVector* get_rhs_DS_velocity_y( size_t const& comp ) ;

      /** @brief Return the DS rhs in z vector */
      LA_SeqVector* get_rhs_DS_velocity_z( size_t const& comp ) ;

      // Domain Decomposition Functions in x for velocity

      /** @brief Return the local rhs vector in x direction */
      LA_SeqVector* get_local_temp_x( size_t const& c ) ;

      /** @brief Return the solution vector for local unknowns in x direction */
      LA_SeqVector* get_local_solution_temp_x( size_t const& c ) ;

      /** @brief Return the Aei*xi vector in x direction */
      LA_SeqVector* get_temp_x( size_t const& c ) ;

      /** @brief Return the solution vector for interface unknowns in x direction */
      LA_SeqVector* get_interface_temp_x( size_t const& c ) ;

      /** @brief Return the domain decomposition matrices in x direction for velocity */

      LA_SeqVector* get_aii_main_diag( size_t const& c, size_t const dir ) ;
      LA_SeqVector* get_aii_super_diag( size_t const& c, size_t const dir ) ;
      LA_SeqVector* get_aii_mod_super_diag( size_t const& c, size_t const dir ) ;
      LA_SeqVector* get_aii_sub_diag( size_t const& c, size_t const dir ) ;

      LA_SeqMatrix* get_aie( size_t const& c, size_t const dir ) ;
      LA_SeqMatrix* get_aei( size_t const& c, size_t const dir ) ;
      LA_SeqMatrix* get_Aee_matrix( size_t const& c, size_t const dir ) ;

      LA_SeqVector* get_product_result( size_t const& c, size_t const dir );
      LA_SeqMatrix* get_Aei_Aii_Aie_product( size_t const& c, size_t const dir ) ;
      LA_SeqVector* get_Aii_Aie_product( size_t const& c, size_t const dir ) ;

       /** @brief Return the schlur complement matrix in x direction */
      LA_SeqMatrix* get_schlur_complement( size_t const& c, size_t const dir ) ; 

       /** @brief Return the Aee matrix in x */
      LA_SeqVector* get_U_vec_u( size_t const& c, size_t const dir ) ;

       /** @brief Return the schlur complement matrix in x direction */
      LA_SeqVector* get_U_vec_v( size_t const& c, size_t const dir ) ;      

      // Domain Decomposition Functions in y for velocity

      /** @brief Return the local rhs vector in y direction */
      LA_SeqVector* get_local_temp_y( size_t const& comp ) ;

      /** @brief Return the solution vector for local unknowns in y direction */
      LA_SeqVector* get_local_solution_temp_y( size_t const& comp ) ;

      /** @brief Return the Aei*xi vector in y direction */
      LA_SeqVector* get_temp_y( size_t const& comp ) ;

      /** @brief Return the solution vector for interface unknowns in y direction */
      LA_SeqVector* get_interface_temp_y( size_t const& comp ) ;

      /** @brief Return the domain decomposition matrices in y direction for velocity */

      /** @brief Return the local rhs vector in z direction */
      LA_SeqVector* get_local_temp_z( size_t const& comp ) ;

      /** @brief Return the solution vector for local unknowns in z direction */
      LA_SeqVector* get_local_solution_temp_z( size_t const& comp ) ;

      /** @brief Return the Aei*xi vector in z direction */
      LA_SeqVector* get_temp_z( size_t const& comp ) ;

      /** @brief Return the solution vector for interface unknowns in z direction */
      LA_SeqVector* get_interface_temp_z( size_t const& comp ) ;

      // Domain Decomposition Functions in x for pressure

      /** @brief Return the local rhs vector in x direction */
      LA_SeqVector* get_local_temp_x_P( void ) ;

      /** @brief Return the solution vector for local unknowns in x direction */
      LA_SeqVector* get_local_solution_temp_x_P( void ) ;

      /** @brief Return the Aei*xi vector in x direction */
      LA_SeqVector* get_temp_x_P( void ) ;

      /** @brief Return the solution vector for interface unknowns in x direction */
      LA_SeqVector* get_interface_temp_x_P( void ) ;

      /** @brief Return the domain decomposition matrices in x direction for pressure */

      LA_SeqVector* get_aii_main_diag_P( size_t const dir ) ;
      LA_SeqVector* get_aii_super_diag_P( size_t const dir ) ;
      LA_SeqVector* get_aii_mod_super_diag_P( size_t const dir ) ;
      LA_SeqVector* get_aii_sub_diag_P( size_t const dir ) ;

      LA_SeqMatrix* get_aie_P( size_t dir ) ;
      LA_SeqMatrix* get_aei_P( size_t dir ) ;
      LA_SeqMatrix* get_Aee_matrix_P( size_t dir ) ;

      LA_SeqMatrix* get_schlur_complement_P( size_t dir ) ;
       /** @brief Return the Aei*inv(Aii)*Aie product in x direction */
      LA_SeqVector* get_product_result_P( size_t const dir );
      LA_SeqMatrix* get_Aei_Aii_Aie_product_P( size_t dir ) ;
      LA_SeqVector* get_Aii_Aie_product_P( size_t dir ) ;

       /** @brief Return the Aee matrix in x */
      LA_SeqVector* get_P_vec_xu( void ) ;

       /** @brief Return the schlur complement matrix in x direction */
      LA_SeqVector* get_P_vec_xv( void ) ;      


      // Domain Decomposition Functions in y for pressure

      /** @brief Return the local rhs vector in y direction */
      LA_SeqVector* get_local_temp_y_P( void ) ;

      /** @brief Return the solution vector for local unknowns in y direction */
      LA_SeqVector* get_local_solution_temp_y_P( void ) ;

      /** @brief Return the Aei*xi vector in y direction */
      LA_SeqVector* get_temp_y_P( void ) ;

      /** @brief Return the solution vector for interface unknowns in y direction */
      LA_SeqVector* get_interface_temp_y_P( void ) ;

       /** @brief Return the Aee matrix in x */
      LA_SeqVector* get_P_vec_yu( void ) ;

       /** @brief Return the schlur complement matrix in x direction */
      LA_SeqVector* get_P_vec_yv( void ) ;      


      // Domain Decomposition Functions in z for pressure

      /** @brief Return the local rhs vector in z direction */
      LA_SeqVector* get_local_temp_z_P( void ) ;

      /** @brief Return the solution vector for local unknowns in z direction */
      LA_SeqVector* get_local_solution_temp_z_P( void ) ;

      /** @brief Return the Aei*xi vector in z direction */
      LA_SeqVector* get_temp_z_P( void ) ;

      /** @brief Return the solution vector for interface unknowns in z direction */
      LA_SeqVector* get_interface_temp_z_P( void ) ;

       /** @brief Return the Aee matrix in x */
      LA_SeqVector* get_P_vec_zu( void ) ;

       /** @brief Return the schlur complement matrix in x direction */
      LA_SeqVector* get_P_vec_zv( void ) ;      

   //-- Basic operations on matrices & vectors

      /** @name Basic operations on matrices & vectors */
      //@{
      /** @brief Initialize the velocity unknown vector with field values */
      void initialize_DS_velocity( void );

      /** @brief Initialize the pressure unknown vector with field values */
      void initialize_DS_pressure( void ); 
      
      /** @brief Initialize the velocity unknown vector with field values */
      void initialize_velocity( void );       

      /** @brief Finalize constant matrices */
      void finalize_constant_matrices( void ) ;

      /** @brief Store velocity vector at previous time step */
      void at_each_time_step( void ) ;

      /** @brief Compute velocity change from one time step to the
      next one */
      double compute_velocity_change( void );

      /** @brief Compute velocity change from one time step to the
      next one with the direction splitting solution method */
      double compute_directionsplitting_velocity_change( void );

      /** @brief Assemble velocity unsteady matrix
      @param coef_lap mass coefficient */
      void assemble_velocity_unsteady_matrix( double const& coef ) ;

      /** @brief Assemble velocity diffusion matrix and rhs
      @param coef_lap laplacian coefficient */
      void assemble_velocity_diffusion_matrix_rhs( double const& coef_lap ) ;

      /** @brief Return the buoyancy vector */
      LA_Vector* get_diffrhs_plus_bodyterm_vector( void ) ;

   //-- Solver

      /** @name Solvers */
      //@{
      /** @brief Heat equation solver */
      bool NavierStokes_solver( void ) ;

      /** @brief Solve the DS splitting problem in x by performing the
      matrix-vector product A_x^-1.Vx and transfer in the distributed vector */
      void DS_NavierStokes_x_solver( size_t const& j, size_t const& k, size_t const& min_i, LA_SeqVector* rhs, LA_SeqVector* interface_rhs , size_t const& comp) ;

      void DS_NavierStokes_x_solver_periodic( size_t const& j, size_t const& k, size_t const& min_i, LA_SeqVector* rhs, LA_SeqVector* interface_rhs , size_t const& comp) ;

      /** @brief Solve the DS splitting problem in y by performing the
      matrix-vector product A_y^-1.Vy and transfer in the distributed vector */
      void DS_NavierStokes_y_solver( size_t const& i, size_t const& k, size_t const& min_j, LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp   ) ;

      void DS_NavierStokes_y_solver_periodic( size_t const& i, size_t const& k, size_t const& min_j, LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp   ) ;

      /** @brief Solve the DS splitting problem in z by performing the
      matrix-vector product A_z^-1.Vz and transfer in the distributed vector */
      void DS_NavierStokes_z_solver( size_t const& i, size_t const& j, size_t const& min_k, LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp  ) ;

      void DS_NavierStokes_z_solver_periodic( size_t const& i, size_t const& j, size_t const& min_k, LA_SeqVector* rhs, LA_SeqVector* interface_rhs, size_t const& comp  ) ;

      /** @brief Solve the domain decomposition problem in x by solving for interface 
      unknowns using schlur complement */
      void DS_NavierStokes_x_interface_unknown_solver( LA_SeqVector* rhs, size_t const& comp ) ;

      /** @brief Solve the domain decomposition problem in x by obtaining the interface
      unknowns and solving for the local unknowns */
      void DS_NavierStokes_x_local_unknown_solver( LA_SeqVector* rhs, size_t const& comp ) ;

      /** @brief Solve the domain decomposition problem in y by solving for interface 
      unknowns using schlur complement */
      void DS_NavierStokes_y_interface_unknown_solver( LA_SeqVector* rhs, size_t const& comp ) ;

      /** @brief Solve the domain decomposition problem in y by obtaining the interface
      unknowns and solving for the local unknowns */
      void DS_NavierStokes_y_local_unknown_solver( LA_SeqVector* rhs, size_t const& comp ) ;

      /** @brief Solve the domain decomposition problem in z by solving for interface 
      unknowns using schlur complement */
      void DS_NavierStokes_z_interface_unknown_solver( LA_SeqVector* rhs, size_t const& comp ) ;

      /** @brief Solve the domain decomposition problem in z by obtaining the interface
      unknowns and solving for the local unknowns */
      void DS_NavierStokes_z_local_unknown_solver( LA_SeqVector* rhs, size_t const& comp ) ;

      /** @brief Solve the DS splitting problem in x by performing the
      matrix-vector product A_x^-1.Vx and transfer in the distributed vector for pressure */
      void DS_NavierStokes_x_solver_P( size_t const& j, size_t const& k, size_t const& min_i, LA_SeqVector* rhs, LA_SeqVector* interface_rhs) ;

      void DS_NavierStokes_x_solver_P_periodic( size_t const& j, size_t const& k, size_t const& min_i, LA_SeqVector* rhs, LA_SeqVector* interface_rhs) ;

      /** @brief Solve the domain decomposition problem in x by solving for interface 
      unknowns using schlur complement for pressure */
      void DS_NavierStokes_x_interface_unknown_solver_P( LA_SeqVector* rhs) ;

      /** @brief Solve the domain decomposition problem in x by obtaining the interface
      unknowns and solving for the local unknowns for pressure */
      void DS_NavierStokes_x_local_unknown_solver_P( LA_SeqVector* rhs ) ;

      /** @brief Solve the DS splitting problem in y by performing the
      matrix-vector product A_y^-1.Vy and transfer in the distributed vector for pressure */
      void DS_NavierStokes_y_solver_P( size_t const& i, size_t const& k, size_t const& min_j, LA_SeqVector* rhs, LA_SeqVector* interface_rhs ) ;

      void DS_NavierStokes_y_solver_P_periodic( size_t const& i, size_t const& k, size_t const& min_j, LA_SeqVector* rhs, LA_SeqVector* interface_rhs ) ;

      /** @brief Solve the domain decomposition problem in y by solving for interface 
      unknowns using schlur complement for pressure */
      void DS_NavierStokes_y_interface_unknown_solver_P( LA_SeqVector* rhs ) ;

      /** @brief Solve the domain decomposition problem in y by obtaining the interface
      unknowns and solving for the local unknowns for pressure */
      void DS_NavierStokes_y_local_unknown_solver_P( LA_SeqVector* rhs ) ;

      /** @brief Solve the DS splitting problem in z by performing the
      matrix-vector product A_z^-1.Vz and transfer in the distributed vector for pressure */
      void DS_NavierStokes_z_solver_P( size_t const& i, size_t const& j, size_t const& min_k, LA_SeqVector* rhs, LA_SeqVector* interface_rhs ) ;

      void DS_NavierStokes_z_solver_P_periodic( size_t const& i, size_t const& j, size_t const& min_k, LA_SeqVector* rhs, LA_SeqVector* interface_rhs ) ;

      /** @brief Solve the domain decomposition problem in z by solving for interface 
      unknowns using schlur complement for pressure */
      void DS_NavierStokes_z_interface_unknown_solver_P( LA_SeqVector* rhs ) ;

      /** @brief Solve the domain decomposition problem in z by obtaining the interface
      unknowns and solving for the local unknowns for pressure */
      void DS_NavierStokes_z_local_unknown_solver_P( LA_SeqVector* rhs ) ;


      //@}

      /** @brief Synchronize the solution vector for velocity*/
      void synchronize_solution_vec( void );

      /** @brief Synchronize the DS solution vector for velocity*/
      void synchronize_DS_solution_vec( void );

      /** @brief Synchronize the DS solution vector for pressure*/
      void synchronize_DS_solution_vec_P( void );

   //-- Output methods

      /** @name Output methods */
      //@{
      /** @brief Display matrices and vectors for debugging purposes */
      void display_debug( void );
      //@}

      /** @brief Compute the product of Aei*inv(Aii)*Aie in x for Velocity*/
      void compute_product_matrix( size_t const& comp, size_t const dir );

      void compute_product_matrix_interior(size_t const& comp, size_t const column, size_t const dir);

      void compute_product_matrix_P( size_t const dir );

      void compute_product_matrix_interior_P(size_t const column, size_t const dir);
      /** @brief Compute the product of Aei*inv(Aii)*Aie in z*/
      double compute_vector_transpose_product(LA_SeqVector* a,LA_SeqVector* b);

   //-- Utilities

      /** @name Utilities */
      //@{
      /** @brief Solve Linear system mat_A*x = rhs using thomas algorithm  */
      static void thomas_algorithm( LA_SeqMatrix* mat_A,LA_SeqVector* rhs) ;

      /** @brief Solve Linear system mat_A*x = rhs with only three vectors of mat_A(x,y,z) using thomas algorithm  */
      static void mod_thomas_algorithm( LA_SeqVector* x,LA_SeqVector* y,LA_SeqVector* z,LA_SeqVector* rhs) ;

      /** @brief Compute the modified super diagonal in Aii_x of velocity for thomas algorithm  */
      void compute_Aii_ref(size_t const& comp, size_t const dir);
      void compute_Aii_ref_P(size_t const dir);

      /** @brief Compute the modified super diagonal in schlur_x of velocity for thomas algorithm  */
      void compute_schlur_x_ref(size_t const& comp);
      /** @brief Compute the modified super diagonal in schlur_y of velocity for thomas algorithm  */
      void compute_schlur_y_ref(size_t const& comp);
      /** @brief Compute the modified super diagonal in schlur_z of velocity for thomas algorithm  */
      void compute_schlur_z_ref(size_t const& comp);

      /** @brief Compute the modified super diagonal in Aii_y of pressure for thomas algorithm  */
      void compute_schlur_x_ref_P(void);
      /** @brief Compute the modified super diagonal in schlur_x of pressure for thomas algorithm  */
      void compute_schlur_y_ref_P(void);
      /** @brief Compute the modified super diagonal in schlur_z of pressure for thomas algorithm  */
      void compute_schlur_z_ref_P(void);

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

      FV_DiscreteField* UF ;
      FV_DiscreteField* PF ; 

      // Local vectors
      LA_SeqVector * UF_LOC ;
      LA_SeqVector * UF_DS_LOC ;
      LA_SeqVector * PF_DS_LOC ;

      // Unknowns vectors
      LA_Vector * VEC_UF ;
      LA_Vector * VEC_UF_previoustime ;
      LA_Vector * VEC_UF_timechange ;
      // Matrices & rhs
      LA_Matrix * MAT_D_velocityUnsteadyPlusDiffusion ;
      LA_Matrix * MAT_A_velocityUnsteady ;
      LA_Vector * VEC_rhs_D_velocityDiffusionPlusBodyTerm ;
      LA_Vector * VEC_rhs_A_velocityUnsteady ;

      // Solvers
      LA_Solver* SOLVER_velocity ;

      // Unknowns numbering
      FV_SystemNumbering* UF_NUM ;
      FV_SystemNumbering* PF_NUM ;

      // Direction splitting matrices
      
      LA_SeqMatrix * MAT_velocityUnsteadyPlusDiffusion_1D_y ;
      LA_SeqMatrix * MAT_velocityUnsteadyPlusDiffusion_1D_y_L ;
      LA_SeqMatrix * MAT_velocityUnsteadyPlusDiffusion_1D_y_U ;
      LA_SeqMatrix * MAT_velocityUnsteadyPlusDiffusion_1D_z ;
      LA_SeqMatrix * MAT_velocityUnsteadyPlusDiffusion_1D_z_L ;
      LA_SeqMatrix * MAT_velocityUnsteadyPlusDiffusion_1D_z_U ;

      // First Order Direction splitting vectors
      LA_SeqVector ** VEC_rhs_velocity_1D_x ;
      LA_SeqVector ** VEC_rhs_velocity_1D_y ;
      LA_SeqVector ** VEC_rhs_velocity_1D_z ;
      LA_SeqVector ** VEC_rhs_velocity_1D_x_P;


      // Domain decomposition matrices in x for velocity
      LA_SeqVector ** Aii_x_main_diagonal;
      LA_SeqVector ** Aii_x_super_diagonal;
      LA_SeqVector ** Aii_x_sub_diagonal;
      LA_SeqVector ** Aii_x_mod_super_diagonal;

      LA_SeqMatrix ** Aie_x ;
      LA_SeqMatrix ** Aei_x ;

      // Aei*inv(Aii)*Aie product matrix in x
      LA_SeqMatrix ** Aei_Aii_Aie_product_x ;

      // inv(Aii)*Aie product vector in x
      LA_SeqVector ** product_result_x ; 

      // Aei*inv(Aii)*Aie product vector (only for one column) in x
      LA_SeqVector ** Aii_Aie_product_x ;

      // Matrices for interface unknowns (only on master proc) in x
      LA_SeqMatrix ** Aee_x;
      LA_SeqMatrix ** schlur_complement_x ;
      LA_SeqMatrix ** schlur_complement_x_ref ;

      // Rhs for local unknowns in x
      LA_SeqVector ** VEC_local_temp_x ;

      // RHS for interface unknowns in x
      LA_SeqVector ** VEC_interface_temp_x ;   

      // Initial Solution vector for local unknowns in x
      LA_SeqVector ** VEC_local_solution_temp_x ;   

      // Vector for sending Aei*xi product in x
      LA_SeqVector ** VEC_temp_x ;      

      // Domain decomposition matrices in y for velocity
      LA_SeqVector ** Aii_y_main_diagonal;
      LA_SeqVector ** Aii_y_super_diagonal;
      LA_SeqVector ** Aii_y_sub_diagonal;
      LA_SeqVector ** Aii_y_mod_super_diagonal;

      LA_SeqMatrix ** Aie_y ;
      LA_SeqMatrix ** Aei_y ;

      // Aei*inv(Aii)*Aie product matrix in y
      LA_SeqMatrix ** Aei_Aii_Aie_product_y ;

      // inv(Aii)*Aie product vector in x
      LA_SeqVector ** product_result_y ; 

      // Aei*inv(Aii)*Aie product vector (only for one column) in y
      LA_SeqVector ** Aii_Aie_product_y ;

      // Matrices for interface unknowns (only on master proc) in y
      LA_SeqMatrix ** Aee_y;
      LA_SeqMatrix ** schlur_complement_y ;
      LA_SeqMatrix ** schlur_complement_y_ref ;

      // Rhs for local unknowns in y
      LA_SeqVector ** VEC_local_temp_y ;

      // RHS for interface unknowns in y
      LA_SeqVector ** VEC_interface_temp_y ;   

      // Initial Solution vector for local unknowns in y
      LA_SeqVector ** VEC_local_solution_temp_y ;   

      // Vector for sending Aei*xi product in y
      LA_SeqVector ** VEC_temp_y ;

      // Domain decomposition matrices in z for velocity
      LA_SeqMatrix ** Aii_z;
      LA_SeqVector ** Aii_z_main_diagonal;
      LA_SeqVector ** Aii_z_super_diagonal;
      LA_SeqVector ** Aii_z_sub_diagonal;
      LA_SeqVector ** Aii_z_mod_super_diagonal;
      LA_SeqMatrix ** Aii_z_ref ;

      LA_SeqMatrix ** Aie_z ;
      LA_SeqMatrix ** Aei_z ;

      // Aei*inv(Aii)*Aie product matrix in z
      LA_SeqMatrix ** Aei_Aii_Aie_product_z ;

      // inv(Aii)*Aie product vector in z
      LA_SeqVector ** product_result_z ; 

      // Aei*inv(Aii)*Aie product vector (only for one column) in z
      LA_SeqVector ** Aii_Aie_product_z ;

      // Matrices for interface unknowns (only on master proc) in z
      LA_SeqMatrix ** Aee_z;
      LA_SeqMatrix ** schlur_complement_z ;
      LA_SeqMatrix ** schlur_complement_z_ref;

      // Rhs for local unknowns in z
      LA_SeqVector ** VEC_local_temp_z;

      // RHS for interface unknowns in z
      LA_SeqVector ** VEC_interface_temp_z ;   

      // Initial Solution vector for local unknowns in z
      LA_SeqVector ** VEC_local_solution_temp_z ;   

      // Vector for sending Aei*xi product in z
      LA_SeqVector ** VEC_temp_z ;

      // Domain decomposition matrices in x for Pressure
      LA_SeqVector * Aii_x_main_diagonal_P;
      LA_SeqVector * Aii_x_super_diagonal_P;
      LA_SeqVector * Aii_x_sub_diagonal_P;
      LA_SeqVector * Aii_x_mod_super_diagonal_P;

      LA_SeqMatrix * Aie_x_P ;
      LA_SeqMatrix * Aei_x_P ;

      // Aei*inv(Aii)*Aie product matrix in x for P
      LA_SeqMatrix * Aei_Aii_Aie_product_x_P ;

      // inv(Aii)*Aie product vector in x for P
      LA_SeqVector * product_result_x_P ; 

      // Aei*inv(Aii)*Aie product vector (only for one column) in x for P
      LA_SeqVector * Aii_Aie_product_x_P ;

      // Matrices for interface unknowns (only on master proc) in x for P
      LA_SeqMatrix * Aee_x_P;
      LA_SeqMatrix * schlur_complement_x_P ;
      LA_SeqMatrix * schlur_complement_x_ref_P ;

      // Rhs for local unknowns in x for P
      LA_SeqVector * VEC_local_temp_x_P ;

      // RHS for interface unknowns in x for P
      LA_SeqVector * VEC_interface_temp_x_P ;   

      // Initial Solution vector for local unknowns in x for P
      LA_SeqVector * VEC_local_solution_temp_x_P ;   

      // Vector for sending Aei*xi product in x for P
      LA_SeqVector * VEC_temp_x_P ;

      // Domain decomposition matrices in y for Pressure
      LA_SeqVector * Aii_y_main_diagonal_P;
      LA_SeqVector * Aii_y_super_diagonal_P;
      LA_SeqVector * Aii_y_sub_diagonal_P;
      LA_SeqVector * Aii_y_mod_super_diagonal_P;

      LA_SeqMatrix * Aie_y_P ;
      LA_SeqMatrix * Aei_y_P ;

      // Aei*inv(Aii)*Aie product matrix in y for P
      LA_SeqMatrix * Aei_Aii_Aie_product_y_P ;

      // inv(Aii)*Aie product vector in y for P
      LA_SeqVector * product_result_y_P ; 

      // Aei*inv(Aii)*Aie product vector (only for one column) in y for P
      LA_SeqVector * Aii_Aie_product_y_P ;

      // Matrices for interface unknowns (only on master proc) in y for P
      LA_SeqMatrix * Aee_y_P;
      LA_SeqMatrix * schlur_complement_y_P ;
      LA_SeqMatrix * schlur_complement_y_ref_P ;

      // Rhs for local unknowns in y for P
      LA_SeqVector * VEC_local_temp_y_P ;

      // RHS for interface unknowns in y for P
      LA_SeqVector * VEC_interface_temp_y_P ;   

      // Initial Solution vector for local unknowns in y for P
      LA_SeqVector * VEC_local_solution_temp_y_P ;   

      // Vector for sending Aei*xi product in y for P
      LA_SeqVector * VEC_temp_y_P ;

      // Domain decomposition matrices in z for Pressure
      LA_SeqVector * Aii_z_main_diagonal_P;
      LA_SeqVector * Aii_z_super_diagonal_P;
      LA_SeqVector * Aii_z_sub_diagonal_P;
      LA_SeqVector * Aii_z_mod_super_diagonal_P;

      LA_SeqMatrix * Aie_z_P ;
      LA_SeqMatrix * Aei_z_P ;

      // Aei*inv(Aii)*Aie product matrix in z for P
      LA_SeqMatrix * Aei_Aii_Aie_product_z_P ;

      // inv(Aii)*Aie product vector in z for P
      LA_SeqVector * product_result_z_P ; 

      // Aei*inv(Aii)*Aie product vector (onlz for one column) in z for P
      LA_SeqVector * Aii_Aie_product_z_P ;

      // Matrices for interface unknowns (onlz on master proc) in z for P
      LA_SeqMatrix * Aee_z_P;
      LA_SeqMatrix * schlur_complement_z_P ;
      LA_SeqMatrix * schlur_complement_z_ref_P ;

      // Rhs for local unknowns in z for P
      LA_SeqVector * VEC_local_temp_z_P ;

      // RHS for interface unknowns in z for P
      LA_SeqVector * VEC_interface_temp_z_P ;   

      // Initial Solution vector for local unknowns in z for P
      LA_SeqVector * VEC_local_solution_temp_z_P ;   

      // Vector for sending Aei*xi product in z for P
      LA_SeqVector * VEC_temp_z_P ;

      // Global velocity solution vectors
      LA_Vector * VEC_DS_UF ;
      LA_Vector * VEC_DS_UF_previoustime ;
      LA_Vector * VEC_DS_UF_timechange ;
      // Global pressure solution vectors
      LA_Vector * VEC_DS_PF ;

      // Two additional vectors in x for 1 proc periodic scenario for velocity
      LA_SeqVector ** U_vec_xu ;
      LA_SeqVector ** U_vec_xv ;

      // Two additional vectors in y for 1 proc periodic scenario for velocity
      LA_SeqVector ** U_vec_yu ;
      LA_SeqVector ** U_vec_yv ;

      // Two additional vectors in z for 1 proc periodic scenario for velocity
      LA_SeqVector ** U_vec_zu ;
      LA_SeqVector ** U_vec_zv ;

      // Two additional vectors in x for 1 proc periodic scenario for pressure
      LA_SeqVector * P_vec_xu ;
      LA_SeqVector * P_vec_xv ;

      // Two additional vectors in y for 1 proc periodic scenario for pressure
      LA_SeqVector * P_vec_yu ;
      LA_SeqVector * P_vec_yv ;

      // Two additional vectors in z for 1 proc periodic scenario for pressure
      LA_SeqVector * P_vec_zu ;
      LA_SeqVector * P_vec_zv ;

      size_t dim;
      MAC_Communicator const* pelCOMM;
      size_t nb_comps;

      /** Processor positions in x,y,z */
      size_t proc_pos_in_i[3];
      /** Number of Processors in x,y,z */
      size_t nb_procs_in_i[3];

      bool is_Uperiodic[3];
      boolVector const* U_periodic_comp;

      bool is_Pperiodic[3];
      boolVector const* P_periodic_comp;
} ;

#endif
