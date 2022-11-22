#ifndef REG_HEAT_EQUATION_SYSTEM_HH
#define REG_HEAT_EQUATION_SYSTEM_HH

#include <MAC_Object.hh>
#include <utility>
#include <boolVector.hh>
#include <size_t_array2D.hh>
#include <doubleArray2D.hh>
#include <vector>
using namespace std;


class MAC_ModuleExplorer ;
class MAC_Communicator ;
class MAC_Timer ;
class size_t_vector ;
class intVector ;
class doubleVector;
class doubleArray2D;
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

/** For set of variables to pass from NavierStokes to System */
struct NS2System
{
  bool is_solids_ ;
  bool is_stressCal_ ;
};

/** @brief TDMatrix include all elements of block matrices (ii,ie,ei,ee) */
struct TDMatrix {
   LA_SeqVector *** ii_main;
   LA_SeqVector *** ii_super;
   LA_SeqVector *** ii_sub;
   LA_SeqMatrix *** ie;
   LA_SeqMatrix *** ei;
   LA_SeqMatrix *** ee;
};

/** @brief Product matrix is composed of products
elements of block matrices (ii,ie,ei,ee) */
struct ProdMatrix {
   LA_SeqMatrix ** ei_ii_ie;
   LA_SeqVector ** ii_ie;
   LA_SeqVector ** result;
};

/** @brief LocalVector to be used in storing the
local values and interface values of DOF */
struct LocalVector {
   LA_SeqVector** local_T;
   LA_SeqVector** local_solution_T;
   LA_SeqVector** T;
   LA_SeqVector** interface_T;
};

/** @brief NodeProp to be used to store the
nodes properties due to presence of solid particles in the domian */
struct NodeProp {
   LA_SeqVector * void_frac;
   // void_fraction of the node due to particle
   LA_SeqVector * parID;
   // ID of solid particle on the node
   LA_SeqVector * bound_cell;
   // Stores the boundary cell presence in the solids; 1 == 1st boundary cell
};

/** @brief The Class DS_NavierStokesSystem.

Matrix systems for the resolution of the heat equation.

@author A. Goyal - Pacific project 2022 */

class DS_NavierStokesSystem : public MAC_Object
{
   private: //----------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */
      DS_NavierStokesSystem( void ) ;

      /** @brief Destructor */
      ~DS_NavierStokesSystem( void ) ;

      /** @brief Copy constructor */
      DS_NavierStokesSystem( DS_NavierStokesSystem const& other ) ;

      /** @brief Operator ==
      @param other the right hand side */
      DS_NavierStokesSystem& operator=( DS_NavierStokesSystem const& other ) ;

      /** @brief Constructor with arguments
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param mac_UF FV velocity field */
      DS_NavierStokesSystem ( MAC_Object* a_owner,
            MAC_ModuleExplorer const* exp,
            FV_DiscreteField* mac_UF,
            FV_DiscreteField* mac_PF ,
            struct NS2System const& fromNS );
      //@}


   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      /** @name Instance delivery and initialization */
      //@{
      /** @brief Create and initialize an instance of DS_NavierStokesSystem
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param mac_UF FV velocity field */
      static DS_NavierStokesSystem* create(
                                          MAC_Object* a_owner,
                                          MAC_ModuleExplorer const* exp,
                                          FV_DiscreteField* mac_UF,
                                          FV_DiscreteField* mac_PF,
                                          struct NS2System const& transfer );
      //@}


   //-- Access

      /** @name Access */
      //@{
      /** @brief Return the matrix system of spacial discretization */
      TDMatrix* get_A(size_t const& field);
      /** @brief Return the Schur complement of spacial discretization */
      TDMatrix* get_Schur(size_t const& field);
      /** @brief Return the (presence/absence) of particle vector */
      NodeProp get_node_property(size_t const& field, size_t const& time_level);
      /** @brief Return the divergence on pressure node */
      doubleVector* get_node_divergence(size_t const& level);
      /** @brief Return the velocity diffusive terms */
      vector<doubleVector*> get_velocity_diffusion();
      /** @brief Return the velocity advection terms */
      vector<doubleVector*> get_velocity_advection();
      /** @brief Return the velocity field face fractions */
      vector<doubleArray2D*> get_velocity_face_fractions();
      vector<doubleArray2D*> get_velocity_normalRB();
      /** @brief Return the local vector with a vector of row index */
      size_t_array2D* get_row_indexes(size_t const& field
                                  , size_t const& dir
                                  , size_t const& comp);

      /** @brief Return the Schur complement of
      Schur complement in case of periodic domain */
      TDMatrix* get_DoubleSchur(size_t const& field);
      /** @brief Return the product matrix of Schur complement */
      ProdMatrix* get_SchurP(size_t const& field);
      /** @brief Return RHS for the Schur complement */
      LocalVector* get_Schur_VEC(size_t const& field);
      /** @brief Return the product matrix of spacial discretization */
      ProdMatrix* get_Ap(size_t const& field);
      /** @brief Return RHS for the matrix system of spacial discretization */
      LocalVector* get_VEC(size_t const& field);

   //-- Solver

      /** @name Solvers */
      //@{
      /** @brief Solve the DS splitting problem in x by performing the
      matrix-vector product A_x^-1.Vx and transfer in the distributed vector */
      void DS_NavierStokes_solver( FV_DiscreteField* FF
                                 , size_t const& j
                                 , size_t const& k
                                 , size_t const& min_i
                                 , size_t const& comp
                                 , size_t const& dir
                                 , size_t const& r_index
                                 , size_t const& level ) ;

      //@}

   //-- Output methods

      /** @name Output methods */
      //@{
      /** @brief Display matrices and vectors for debugging purposes */
      void display_debug( void );
      //@}

      /** @brief Calls interior function for different conditions
      to compute the product of Aei*inv(Aii)*Aie */
      void compute_product_matrix(struct TDMatrix *arr
                                , struct ProdMatrix *prr
                                , size_t const& comp
                                , size_t const& dir
                                , size_t const& field
                                , size_t const& r_index );
      /** @brief Compute the product of Aei*inv(Aii)*Aie in any
      direction for any field*/
      void compute_product_matrix_interior(struct TDMatrix *arr
                                         , struct ProdMatrix *prr
                                         , size_t const& comp
                                         , size_t const& column
                                         , size_t const& dir
                                         , size_t const& r_index);

   //-- Utilities

      /** @name Utilities */
      //@{
      /** @brief Solve Linear system mat_A*x = rhs with only three
      vectors of mat_A(x,y,z) using thomas algorithm  */
      void mod_thomas_algorithm(TDMatrix *arr
                              , LA_SeqVector* rhs
                              , size_t const& comp
                              , size_t const& dir
                              , size_t const& r_index);
      /** @brief Compute the modified super diagonal for thomas algorithm  */
      void pre_thomas_treatment( size_t const& comp
                               , size_t const& dir
                               , struct TDMatrix *arr
                               , size_t const& r_index);
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

      // Local vector to store row indexes for each
      // field (i.e. 2)
      // and direction (i.e. 3)
      size_t_array2D ** row_index[2][3] ;

      // Unknowns numbering
      FV_SystemNumbering* UF_NUM ;
      FV_SystemNumbering* PF_NUM ;

      // Direction splitting matrices
      LA_SeqMatrix * MAT_velocityUnsteadyPlusDiffusion_1D ;

      // Spacitial discretization matrices
      // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct TDMatrix A[2][3];
      // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct ProdMatrix Ap[2][3];
      // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct LocalVector VEC[2][3];

      // Schur complement matrices
      // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct TDMatrix Schur[2][3];
      struct ProdMatrix SchurP[2][3];
      struct LocalVector Schur_VEC[2][3];

      // Schur complement of Schur complement
      struct TDMatrix DoubleSchur[2][3];

      // Particle structures
      // 2 rows are for fields; 2 columns are for time level (current and last)
      struct NodeProp node[2][2];
      // 0 current timestep, 1 last time step
      vector<doubleVector*> divergence;
      // Local vector to store diffusive terms
      vector<doubleVector*> vel_diffusion;
      // Local vector to store advection terms
      vector<doubleVector*> vel_advection;

      size_t dim;
      MAC_Communicator const* pelCOMM;
      size_t nb_comps[2];

      /** Processor positions in x,y,z */
      size_t proc_pos_in_i[3];
      /** Number of Processors in x,y,z */
      size_t nb_procs_in_i[3];

      bool is_solids;
      bool is_stressCal;
      size_t Npart;
      string level_set_type;
      double Nmax;
      double Rpart,ar;


      bool is_periodic[2][3];
      boolVector const* U_periodic_comp;
      boolVector const* P_periodic_comp;
} ;

#endif
