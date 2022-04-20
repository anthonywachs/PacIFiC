#ifndef DS_HeatTransferSystem_HH
#define DS_HeatTransferSystem_HH

#include <MAC_Object.hh>
#include <utility>
#include <boolVector.hh>
#include <size_t_array2D.hh>
#include <vector>
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


/** @brief The Class DS_HeatTransferSystem.

Matrix systems for the resolution of the heat equation.

@author A. Goyal - Pacific project 2022 */

class DS_HeatTransferSystem : public MAC_Object
{
   private: //----------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */
      DS_HeatTransferSystem( void ) ;

      /** @brief Destructor */
      ~DS_HeatTransferSystem( void ) ;

      /** @brief Copy constructor */
      DS_HeatTransferSystem( DS_HeatTransferSystem const& other ) ;

      /** @brief Operator ==
      @param other the right hand side */
      DS_HeatTransferSystem& operator=( DS_HeatTransferSystem const& other ) ;

      /** @brief Constructor with arguments
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param mac_tf FV temperature field */
      DS_HeatTransferSystem ( MAC_Object* a_owner,
            MAC_ModuleExplorer const* exp,
            FV_DiscreteField* mac_tf);
      //@}


   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      /** @name Instance delivery and initialization */
      //@{
      /** @brief Create and initialize an instance of DS_HeatTransferSystem
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param mac_tf FV temperature field */
      static DS_HeatTransferSystem* create(
            MAC_Object* a_owner,
            MAC_ModuleExplorer const* exp,
            FV_DiscreteField* mac_tf);
      //@}


   //-- Access

      /** @name Access */
      //@{

      /** @brief Return the local vector with a vector of row index */
      size_t_array2D* get_row_indexes(size_t const& field
                                  , size_t const& dir
                                  , size_t const& comp);

      /** @brief Return the temperature diffusive terms */
      vector<doubleVector*> get_temperature_diffusion();

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

   //-- Basic operations on matrices & vectors


   //-- Solver

      /** @name Solvers */
      //@{

      /** @brief Solve the DS splitting problem in x by performing the
      matrix-vector product A_x^-1.Vx and transfer in the distributed vector */
      void DS_HeatEquation_solver( size_t const& j
                                 , size_t const& k
                                 , size_t const& min_i
                                 , size_t const& comp
                                 , size_t const& dir
                                 , size_t const& r_index
                                 , size_t const& level) ;

      //@}

   //-- Output methods

      /** @name Output methods */
      //@{
      /** @brief Display matrices and vectors for debugging purposes */
      void display_debug( void );
      //@}

      /** @brief Call the interior function for different conditions of procs and periodicity*/
      void compute_product_matrix( struct TDMatrix *arr
                                 , struct ProdMatrix *prr
                                 , size_t const& comp
                                 , size_t const& dir
                                 , size_t const& r_index );
      /** @brief Compute the product of Aei*inv(Aii)*Aie in x*/
      void compute_product_matrix_interior( struct TDMatrix *arr
                                          , struct ProdMatrix *prr
                                          , size_t const& comp
                                          , size_t const& column
                                          , size_t const& dir
                                          , size_t const& r_index);

   //-- Utilities

      /** @name Utilities */
      //@{
      /** @brief Compute the pre-thomas step on the provided matrix arr */
      void pre_thomas_treatment(size_t const& comp
                              , size_t const& dir
                              , struct TDMatrix *arr
                              , size_t const& r_index);
      /** @brief Compute the inverse of 1 tridiagonal matrix  */
      static void mod_thomas_algorithm(TDMatrix *arr
                                     , LA_SeqVector* rhs
                                     , size_t const& comp
                                     , size_t const& dir
                                     , size_t const& r_index) ;

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

      // Local vector to store row indexes for each
      // field (i.e. 1)
      // and direction (i.e. 3)
      size_t_array2D ** row_index[1][3] ;

      // Direction splitting matrices
      LA_SeqMatrix * MAT_TemperatureUnsteadyPlusDiffusion_1D ;

      // Spacitial discretization matrices
      struct TDMatrix A[3];
      struct ProdMatrix Ap[3];
      struct LocalVector VEC[3];

      // Schur complement matrices
      struct TDMatrix Schur[3];
      struct ProdMatrix SchurP[3];
      struct LocalVector Schur_VEC[3];

      // Schur complement of Schur complement
      struct TDMatrix DoubleSchur[3];

      // Local vector to store diffusive terms
      vector<doubleVector*> T_diffusion;

      size_t dim;
      size_t nb_comps;

      size_t proc_pos_in_i[3];
      size_t nb_procs_in_i[3];

      bool is_iperiodic[3];
      boolVector const* periodic_comp;
} ;

#endif
