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

/** @brief TDMatrix include all elements of block matrices (ii,ie,ei,ee) */
struct TDMatrix {
   LA_SeqVector ** ii_main;
   LA_SeqVector ** ii_super;
   LA_SeqVector ** ii_sub;
   LA_SeqMatrix ** ie;
   LA_SeqMatrix ** ei;
   LA_SeqMatrix ** ee;
};

/** @brief Product matrix is composed of products elements of block matrices (ii,ie,ei,ee) */
struct ProdMatrix {
   LA_SeqMatrix ** ei_ii_ie;
   LA_SeqVector ** ii_ie;
   LA_SeqVector ** result;
};

/** @brief LocalVector to be used in storing the local values and interface values of DOF */
struct LocalVector {
   LA_SeqVector** local_T;
   LA_SeqVector** local_solution_T;
   LA_SeqVector** T;
   LA_SeqVector** interface_T;
};

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
      /** @brief Return the DS velocity solution vector DS_UF */
      LA_SeqVector const* get_solution_DS_velocity( void ) const ;

      /** @brief Return the DS pressure solution vector DS_PF */
      LA_SeqVector const* get_solution_DS_pressure( void ) const ;

      /** @brief Return the matrix system of spacial discretization */
      TDMatrix* get_A(size_t const field);
      /** @brief Return the Schur complement of spacial discretization */
      TDMatrix* get_Schur(size_t const field);
      /** @brief Return the Schur complement of Schur complement in case of periodic domain */
      TDMatrix* get_DoubleSchur(size_t const field);
      /** @brief Return the product matrix of Schur complement */
      ProdMatrix* get_SchurP(size_t const field);
      /** @brief Return RHS for the Schur complement */
      LocalVector* get_Schur_VEC(size_t const field);
      /** @brief Return the product matrix of spacial discretization */
      ProdMatrix* get_Ap(size_t const field);
      /** @brief Return RHS for the matrix system of spacial discretization */
      LocalVector* get_VEC(size_t const field);

   //-- Basic operations on matrices & vectors

      /** @name Basic operations on matrices & vectors */
      //@{
      /** @brief Initialize the velocity unknown vector with field values */
      void initialize_DS_velocity( void );
      /** @brief Initialize the pressure unknown vector with field values */
      void initialize_DS_pressure( void ); 
      /** @brief Store velocity vector at previous time step */
      void at_each_time_step( void ) ;
      /** @brief Compute velocity change from one time step to the
      next one with the direction splitting solution method */
      double compute_DS_velocity_change( void );

   //-- Solver

      /** @name Solvers */
      //@{
      /** @brief Solve the DS splitting problem in x by performing the
      matrix-vector product A_x^-1.Vx and transfer in the distributed vector */
      void DS_NavierStokes_solver( FV_DiscreteField* FF, size_t const& j, size_t const& k, size_t const& min_i, size_t const& comp, size_t const dir, size_t const field) ;

      //@}

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

      /** @brief Calls interior function for different conditions to compute the product of Aei*inv(Aii)*Aie */
      void compute_product_matrix(struct TDMatrix *arr, struct ProdMatrix *prr, size_t const& comp, size_t const dir, size_t const field );
      /** @brief Compute the product of Aei*inv(Aii)*Aie in any direction for any field*/
      void compute_product_matrix_interior(struct TDMatrix *arr, struct ProdMatrix *prr, size_t const& comp, size_t const column, size_t const dir);

   //-- Utilities

      /** @name Utilities */
      //@{
      /** @brief Solve Linear system mat_A*x = rhs with only three vectors of mat_A(x,y,z) using thomas algorithm  */
      void mod_thomas_algorithm(TDMatrix *arr, LA_SeqVector* rhs, size_t const& comp, size_t const dir);
      /** @brief Compute the modified super diagonal for thomas algorithm  */
      void pre_thomas_treatment( size_t const& comp, size_t const dir, struct TDMatrix *arr);
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
      LA_SeqVector * UF_DS_LOC ;
      LA_SeqVector * PF_DS_LOC ;

      // Global velocity solution vectors
      LA_Vector * VEC_DS_UF ;
      LA_Vector * VEC_DS_UF_previoustime ;
      LA_Vector * VEC_DS_UF_timechange ;
      // Global pressure solution vectors
      LA_Vector * VEC_DS_PF ;
      // Matrices & rhs
      LA_Matrix * MAT_D_velocityUnsteadyPlusDiffusion ;
      LA_Vector * VEC_rhs_D_velocityDiffusionPlusBodyTerm ;

      // Unknowns numbering
      FV_SystemNumbering* UF_NUM ;
      FV_SystemNumbering* PF_NUM ;

      // Direction splitting matrices
      LA_SeqMatrix * MAT_velocityUnsteadyPlusDiffusion_1D ;

      // Spacitial discretization matrices
      struct TDMatrix A[2][3];          // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct ProdMatrix Ap[2][3];       // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct LocalVector VEC[2][3];     // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions

      // Schur complement matrices
      struct TDMatrix Schur[2][3];      // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct ProdMatrix SchurP[2][3];
      struct LocalVector Schur_VEC[2][3];

      // Schur complement of Schur complement
      struct TDMatrix DoubleSchur[2][3];

      size_t dim;
      MAC_Communicator const* pelCOMM;
      size_t nb_comps[2];

      /** Processor positions in x,y,z */
      size_t proc_pos_in_i[3];
      /** Number of Processors in x,y,z */
      size_t nb_procs_in_i[3];

      bool is_periodic[2][3];
      boolVector const* U_periodic_comp;
      boolVector const* P_periodic_comp;
} ;

#endif
