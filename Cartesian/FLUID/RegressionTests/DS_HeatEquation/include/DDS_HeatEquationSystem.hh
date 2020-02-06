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
   LA_SeqVector *** ii_main;
   LA_SeqVector *** ii_super;
   LA_SeqVector *** ii_sub;
   LA_SeqMatrix *** ie;
   LA_SeqMatrix *** ei;
   LA_SeqMatrix *** ee;
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

/** @brief PartInput to be used to store the Input properties of particles in the domian */
struct PartInput {
   LA_SeqMatrix ** coord;               // Coordinates
   LA_SeqVector ** size;                // Size of the sphere
   LA_SeqVector ** temp;         // Temperature of the sphere
};

/** @brief NodeProp to be used to store the nodes properties due to presence of solid particles in the domian */
struct NodeProp {
   LA_SeqVector ** void_frac;               // void_fraction of the node due to particle
   LA_SeqVector ** parID;                   // ID of solid particle on the node
};

/** @brief BoundaryBisec to be used to store the intersection of solids with grids in each direction */
struct BoundaryBisec {
   LA_SeqMatrix ** offset;                  // Direction of intersection relative to node (Column 0 for left and Column 1 for right) 
   LA_SeqMatrix ** value;                   // Value of offset relative to node point
   LA_SeqMatrix ** field;                   // Value of field variable at the intersection
};

/** @brief The Class DDS_HeatEquationSystem.

Matrix systems for the resolution of the heat equation.

@author A. Wachs - Pacific project 2017 */

class DDS_HeatEquationSystem : public MAC_Object
{
   private: //----------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */
      DDS_HeatEquationSystem( void ) ;

      /** @brief Destructor */
      ~DDS_HeatEquationSystem( void ) ;

      /** @brief Copy constructor */
      DDS_HeatEquationSystem( DDS_HeatEquationSystem const& other ) ;

      /** @brief Operator ==
      @param other the right hand side */
      DDS_HeatEquationSystem& operator=( DDS_HeatEquationSystem const& other ) ;

      /** @brief Constructor with arguments
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param mac_tf FV temperature field */
      DDS_HeatEquationSystem ( MAC_Object* a_owner,
            MAC_ModuleExplorer const* exp,
            FV_DiscreteField* mac_tf );
      //@}


   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      /** @name Instance delivery and initialization */
      //@{
      /** @brief Create and initialize an instance of DDS_HeatEquationSystem
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param mac_tf FV temperature field */
      static DDS_HeatEquationSystem* create(
            MAC_Object* a_owner,
            MAC_ModuleExplorer const* exp,
            FV_DiscreteField* mac_tf ) ;
      //@}


   //-- Access

      /** @name Access */
      //@{

      /** @brief Return the DS temperature solution vector DS_TF */
      LA_SeqVector const* get_solution_DS_temperature( void ) const ;

      /** @brief Return the matrix system of spacial discretization */
      TDMatrix* get_A();
      /** @brief Return the product matrix of spacial discretization */
      ProdMatrix* get_Ap();
      /** @brief Return the Schur complement of spacial discretization */
      TDMatrix* get_Schur();
      /** @brief Return the product matrix of Schur complement */
      ProdMatrix* get_SchurP();
      /** @brief Return the Schur complement of Schur complement in case of periodic domain */
      TDMatrix* get_DoubleSchur();
      /** @brief Return RHS for the matrix system of spacial discretization */
      LocalVector* get_VEC();
      /** @brief Return RHS for the Schur complement */
      LocalVector* get_Schur_VEC();
      /** @brief Return information of intersection with solid boundary */
      BoundaryBisec* get_b_intersect();
      /** @brief Return the (presence/absence) of particle vector */
      NodeProp get_node_property();
      /** @brief Return the particle input properties */
      PartInput get_solid();

   //-- Basic operations on matrices & vectors

      /** @name Basic operations on matrices & vectors */
      //@{
      /** @brief Initialize the temperature unknown vector with field values */
      void initialize_temperature( void );

      /** @brief Store temperature vector at previous time step */
      void at_each_time_step( void ) ;

      /** @brief Compute temperature change from one time step to the
      next one */
      double compute_temperature_change( void );

   //-- Solver

      /** @name Solvers */
      //@{

      /** @brief Solve the DS splitting problem in x by performing the
      matrix-vector product A_x^-1.Vx and transfer in the distributed vector */
      void DS_HeatEquation_solver( size_t const& j, size_t const& k, size_t const& min_i, size_t const& comp, size_t const& dir,size_t const& r_index) ;

      //@}

      /** @brief Synchronize the DS solution vector*/
      void synchronize_DS_solution_vec( void );
   //-- Output methods

      /** @name Output methods */
      //@{
      /** @brief Display matrices and vectors for debugging purposes */
      void display_debug( void );
      //@}

      /** @brief Call the interior function for different conditions of procs and periodicity*/
      void compute_product_matrix( struct TDMatrix *arr, struct ProdMatrix *prr, size_t const& comp, size_t const& dir, size_t const& r_index );
      /** @brief Compute the product of Aei*inv(Aii)*Aie in x*/
      void compute_product_matrix_interior( struct TDMatrix *arr, struct ProdMatrix *prr, size_t const& comp, size_t const& column, size_t const& dir, size_t const& r_index);

   //-- Utilities

      /** @name Utilities */
      //@{
      /** @brief Compute the pre-thomas step on the provided matrix arr */
      void pre_thomas_treatment(size_t const& comp, size_t const& dir, struct TDMatrix *arr, size_t const& r_index);
      /** @brief Compute the inverse of 1 tridiagonal matrix  */
      static void mod_thomas_algorithm(TDMatrix *arr, LA_SeqVector* rhs, size_t const& comp, size_t const& dir, size_t const& r_index) ;

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
      LA_SeqVector * TF_DS_LOC ;

      // Matrices & rhs
      LA_Matrix * MAT_D_TemperatureUnsteadyPlusDiffusion ;

      // Unknowns numbering
      FV_SystemNumbering* TF_NUM ;

      // Direction splitting matrices
      LA_SeqMatrix * MAT_TemperatureUnsteadyPlusDiffusion_1D ;

      // Spacitial discretization matrices
      struct TDMatrix A[3];
      struct ProdMatrix Ap[3];
      struct LocalVector VEC[3];

      // Particle structures
      struct PartInput solid;
      struct NodeProp node;
      struct BoundaryBisec b_intersect[3];

      // Schur complement matrices
      struct TDMatrix Schur[3];
      struct ProdMatrix SchurP[3];
      struct LocalVector Schur_VEC[3];

      // Schur complement of Schur complement
      struct TDMatrix DoubleSchur[3];

      // Global Temperature solution vectors
      LA_Vector * VEC_DS_TF ;
      LA_Vector * VEC_DS_TF_previoustime ;
      LA_Vector * VEC_DS_TF_timechange ;

      size_t dim;
      MAC_Communicator const* pelCOMM;
      size_t nb_comps;
      size_t Npart;
      double Rpart;
      bool is_solids;

      size_t proc_pos_in_i[3];
      size_t nb_procs_in_i[3];

      bool is_iperiodic[3];
      boolVector const* periodic_comp;
} ;

#endif
