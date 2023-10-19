#ifndef DS_HeatTransfer_HH
#define DS_HeatTransfer_HH

#include <mpi.h>
#include <FV_OneStepIteration.hh>
#include <geomVector.hh>
#include <PAC_computingtime.hh>
#include <PAC_solvercomputingtime.hh>
#include <vector>
#include <string>
#include <boolVector.hh>
#include <size_t_vector.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;

class MAC_Communicator ;
class FV_DiscreteField ;
class LA_Vector ;
class LA_SeqVector ;
class DS_HeatTransferSystem ;
class LA_SeqMatrix ;
class DS_AllRigidBodies ;

/** @brief The Class DS_HeatTransfer.

Server for the resolution of the unsteady heat equation by a first order
implicit time integrator and a Finite Volume MAC scheme on rectangular grids.

Equation: rho*cp*dT/dt + rho*cp*(u.nabla(T))= k * lap(T) + bodyterm

@author A. Goyal - Pacific project 2022 */

struct DS2HE
{
  double rho_ ;
  bool b_restart_ ;
  string AdvectionScheme_ ;
  string ViscousStressOrder_ ;
  size_t AdvectionTimeAccuracy_ ;
  bool is_solids_ ;
  bool is_NSwithHE_ ;
  bool is_stressCal_ ;
  size_t stressCalFreq_;
  bool is_par_motion_ ;
  FV_DomainAndFields const* dom_ ;
  DS_AllRigidBodies* allrigidbodies_ ;
};

/** @brief MPIVar include all vectors required while message passing */
struct MPIVarHT {
   size_t *size;
   double ***send;
   double ***receive;
};

class DS_HeatTransfer : public MAC_Object, public PAC_ComputingTime, 
	public PAC_SolverComputingTime
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      /** @brief Create and initialize an instance of REG_HeatTransfer
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param fromNS structure containing input data from the NS solver */
      static DS_HeatTransfer* create(
		MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp,
                struct DS2HE const& transfer );

      /** @name Substeps of the step by step progression */
      //@{
      /** @brief Tasks performed at initialization of the algorithm, before
      starting the time stepping loop
      @param t_it time iterator */
      virtual void do_before_time_stepping( FV_TimeIterator const* t_it,
      	std::string const& basename ) ;

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
      //@}

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Destructor */
      virtual ~DS_HeatTransfer( void ) ;

      /** @brief Copy constructor */
      DS_HeatTransfer( DS_HeatTransfer const& other ) ;

      /** @brief Operator ==
      @param other the right hand side */
      DS_HeatTransfer& operator=( DS_HeatTransfer const& other ) ;

      /** @brief Constructor with arguments
      @param a_owner the MAC-based object
      @param exp to read the data file */
      DS_HeatTransfer( MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp,
                struct DS2HE const& fromDS );

      /** @brief Constructor without argument */
      DS_HeatTransfer( void ) ;

      //@}


   //-- Basic discrete system building

      /** @name Basic discrete system building */
      //@{
      /** @brief Compute diffusive term of temperature from perivious timestep */
      double compute_un_component ( size_t const& comp
                                  , size_t const& i
                                  , size_t const& j
                                  , size_t const& k
                                  , size_t const& dir
                                  , size_t const& level);
      /** @brief Assemble temperature body term */
      double bodyterm_value ( double const& xC
                            , double const& yC
                            , double const& zC);

      /** @brief Assemble RHS of temperature for first step of Crank_Nicolson time discretization */
      void assemble_DS_un_at_rhs ( FV_TimeIterator const* t_it
                                 , double const& gamma);



      /** @brief Call the functions to assemble temperature and schur complement */
      void assemble_temperature_and_schur( FV_TimeIterator const* t_it) ;

      /** @brief Assemble temperature matrix */
      void assemble_temperature_matrix ( FV_DiscreteField const* FF
                                       , FV_TimeIterator const* t_it);

      /** @brief Assemble schur matrix */
      void assemble_schur_matrix ( FV_DiscreteField const* FF);


      void write_output_field();
      //@}

   //-- Solver

      /** @name Solvers */
      //@{

      /** @brief Assemble local RHS for 1D equation solver */
      double assemble_local_rhs( size_t const& j
                               , size_t const& k
                               , double const& gamma
                               , FV_TimeIterator const* t_it
                               , size_t const& comp
                               , size_t const& dir );

      /** @brief Compute Aei*(Aii)-1*fi required to compute interface unknown */
      void compute_Aei_ui (struct TDMatrix* arr
                         , struct LocalVector* VEC
                         , size_t const& comp
                         , size_t const& dir
                         , size_t const& r_index);

      /** @brief Pack Aei*(Aii)-1*fi and fe for sending to master processor */
      void data_packing ( double const& fe
                        , size_t const& comp
                        , size_t const& dir
                        , size_t const& p);

      /** @brief Unpack the data sent by "data_packing" and compute the interface unknown; and pack ue for sending to slave processor */
      void unpack_compute_ue_pack(size_t const& comp
                                , size_t const& dir
                                , size_t const& p);

      /** @brief Unpack the interface variable sent by master processor to slave processor */
      void unpack_ue(size_t const& comp
                   , double * received_data
                   , size_t const& dir
                   , size_t const& p);

      /** @brief Call the appropriate functions to solve local variable and interface unknown */
      void solve_interface_unknowns(double const& gamma
                                  , FV_TimeIterator const* t_it
                                  , size_t const& comp
                                  , size_t const& dir
                                  , size_t const& level );

      /** @brief Call the appropriate functions to solve any particular direction in the other directions */
      void HeatEquation_DirectionSplittingSolver( FV_TimeIterator const* t_it ) ;

      /** @brief Correct the fluxes and variables on the nodes due to presence of solid objects */
      void nodes_temperature_initialization ( size_t const& level );

      /** @brief Compute advective term based on either Upwind or TVD spacial scheme */
      double compute_adv_component ( size_t const& comp
                                   , size_t const& i
                                   , size_t const& j
                                   , size_t const& k);

      double assemble_advection_Centered( FV_DiscreteField const* AdvectingField
                                        , size_t advecting_level
                                        , double const& coef
                                        , size_t const& i
                                        , size_t const& j
                                        , size_t const& k
                                        , size_t advected_level) const;

      double divergence_of_U ( size_t const& comp
                             , size_t const& i
                             , size_t const& j
                             , size_t const& k
                             , size_t const& level);

      double assemble_advection_TVD( FV_DiscreteField const* AdvectingField
                                   , size_t advecting_level
                                   , double const& coef
                                   , size_t const& i
                                   , size_t const& j
                                   , size_t const& k
                                   , size_t advected_level) const;

      double assemble_advection_Upwind( FV_DiscreteField const* AdvectingField
                                      , size_t advecting_level
                                      , double const& coef
                                      , size_t const& i
                                      , size_t const& j
                                      , size_t const& k
                                      , size_t advected_level) const;

      /** @brief Solve i in j and k; e.g. solve x in y ank z */
      void Solve_i_in_jk (FV_TimeIterator const* t_it
                        , double const& gamma
                        , size_t const& dir_i
                        , size_t const& dir_j
                        , size_t const& dir_k
                        , size_t const& level);

      /** @brief Calculate the temperature change in each iteration */
      double compute_DS_temperature_change( );

      void assemble_temperature_diffusion_terms ( );

      void calculate_row_indexes ( );

      /** @brief Solve interface unknown for all cases */
      void DS_interface_unknown_solver(LA_SeqVector* interface_rhs
                                     , size_t const& comp
                                     , size_t const& dir
                                     , size_t const& r_index);

      void ugradu_initialization (  );
      /** @brief Calculate L2 norm without any analytical solution */
      void output_l2norm ( void ) ;
      //@}


   // Direction splitting communicators

      /** @name Direction splitting communicators */
      /** @brief Create the sub-communicators */
      void create_DS_subcommunicators ( void ) ;

      void processor_splitting ( int const& color
                               , int const& key
                               , size_t const& dir );

      void allocate_mpi_variables (void);
      void deallocate_mpi_variables ( void );


      /** @brief Free the sub-communicators */
      void free_DS_subcommunicators ( void ) ;
      //@}

   private: //----------------------------------------------------------------

   //-- Class attributes

      static DS_HeatTransfer const* PROTOTYPE ;

   //-- Attributes

      FV_DiscreteField* TF;
      FV_DiscreteField const* UF;

      FV_DiscreteField* TF_ERROR;

      FV_DiscreteField* TF_DS_ERROR;

      DS_HeatTransferSystem* GLOBAL_EQ ;

      size_t nb_procs;
      size_t my_rank;
      size_t is_master;
      size_t dim;
      size_t nb_comps;

      MAC_Communicator const* macCOMM;
      MPI_Comm DS_Comm_i[3];

      int rank_in_i[3];
      int nb_ranks_comm_i[3];

      struct MPIVarHT first_pass[3];
      struct MPIVarHT second_pass[3];
      struct MPIVarHT data_for_S[3];

      double rho;
      string AdvectionScheme;
      size_t AdvectionTimeAccuracy;
      string ViscousStressOrder;
      double heat_capacity, thermal_conductivity;
      bool b_bodyterm ;
      bool is_solids;
      bool is_NSwithHE;
      bool is_stressCal;
      size_t stressCalFreq;
      bool b_restart ;
      bool is_iperiodic[3];
      boolVector const* periodic_comp;
      bool is_par_motion;

      // Grains3D variable
      DS_AllRigidBodies* allrigidbodies;
} ;

#endif
