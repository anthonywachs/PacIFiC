#ifndef DS_NavierStokes_HH
#define DS_NavierStokes_HH

#include <mpi.h>
#include <FV_OneStepIteration.hh>
#include <geomVector.hh>
#include <computingtime.hh>
#include <boolVector.hh>
#include <solvercomputingtime.hh>
#include <MAC_DoubleVector.hh>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;

class MAC_Communicator ;
class FV_DiscreteField ;
class LA_Vector ;
class LA_SeqVector ;
class DS_NavierStokesSystem ;
class LA_SeqMatrix ;
class DS_AllRigidBodies ;

/** For set of variables to pass from NavierStokes to System */
struct DS2NS
{
  double rho_ ;
  double mu_ ;
  double kai_ ;
  string AdvectionScheme_ ;
  size_t AdvectionTimeAccuracy_ ;
  bool b_restart_ ;
  bool is_solids_ ;
  bool is_stressCal_;
  string ViscousStressOrder_;
  double surface_cell_scale_;
  bool is_surfacestressOUT_;
  size_t stressCalFreq_;
  bool is_par_motion_;
  FV_DomainAndFields const* dom_ ;
  DS_AllRigidBodies* allrigidbodies_ ;
  double critical_distance_translation_ ;
};

/** @brief The Class DS_NavierStokes.

Server for the resolution of the unsteady momentum equation by a second order
time integrator and a Finite Volume MAC scheme on rectangular grids.

@author A. Goyal - Pacific project 2022 */

/** @brief MPIVar include all vectors required while message passing */
struct MPIVarNS {
   size_t *size;
   double ***send;
   double ***receive;
};

class DS_NavierStokes : public MAC_Object,
                        public ComputingTime,
                        public SolverComputingTime
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      /** @brief Create and initialize an instance of DS_NavierStokes
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param fromDS structure containing input data from the DS framework */
      static DS_NavierStokes* create(
                  MAC_Object* a_owner,
                  MAC_ModuleExplorer const* exp,
                  struct DS2NS const& fromDS );

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
      virtual ~DS_NavierStokes( void ) ;

      /** @brief Copy constructor */
      DS_NavierStokes( DS_NavierStokes const& other ) ;

      /** @brief Operator ==
      @param other the right hand side */
      DS_NavierStokes& operator=( DS_NavierStokes const& other ) ;

      /** @brief Constructor with arguments
      @param a_owner the MAC-based object
      @param exp to read the data file */
      DS_NavierStokes( MAC_Object* a_owner,
                       MAC_ModuleExplorer const* exp,
                       struct DS2NS const& fromNS );

      /** @brief Constructor without argument */
      DS_NavierStokes( void ) ;

      //@}


   //-- Basic discrete system building

      /** @name Basic discrete system building */
      //@{
      /** @brief Call the function to assemble 1D matrices
      for both velocity and pressure field*/
      void assemble_1D_matrices( FV_DiscreteField const* FF
                               , FV_TimeIterator const* t_it ) ;

      /** @brief Assemble 1D matrices for both velocity
      and pressure field in all directions */
      void assemble_field_matrix ( FV_DiscreteField const* FF
                                 , FV_TimeIterator const* t_it );

      /** @brief Assemble 1D schur matrices for both
      velocity and pressure field in all directions */
      void assemble_field_schur_matrix ( FV_DiscreteField const* FF );

      /** @brief Assemble advection term for Upwind spacial scheme */
      double assemble_advection_Upwind( size_t const& advecting_level
                                      , double const& coef
                                      , size_t const& advected_level
                                      , size_t const& i
                                      , size_t const& j
                                      , size_t const& k
                                      , size_t const& component) ;

      /** @brief Assemble advection term for Centered spacial scheme */
      double assemble_advection_Centered( size_t const& advecting_level
                                        , double const& coef
                                        , size_t const& advected_level
                                        , size_t const& i
                                        , size_t const& j
                                        , size_t const& k
                                        , size_t const& component );


      /** @brief Assemble advection term for TVD spacial scheme */
      double assemble_advection_TVD( size_t const& advecting_level
                                   , double const& coef
                                   , size_t const& advected_level
                                   , size_t const& i
                                   , size_t const& j
                                   , size_t const& k
                                   , size_t const& component) const;

      /** @brief Assemble rhs for velocity in any direction */
      double velocity_local_rhs( size_t const& j
                               , size_t const& k
                               , double const&  gamma
                               , FV_TimeIterator const* t_it
                               , size_t const& comp
                               , size_t const& dir );
      /** @brief Compute diffusive term of velocity field from previous timestep */
      double compute_un_component ( size_t const& comp
                                  , size_t const& i
                                  , size_t const& j
                                  , size_t const& k
                                  , size_t const& dir
                                  , size_t const& level);
      /** @brief Compute diffusive term of pressure field from previous timestep */
      double compute_p_component ( size_t const& comp
                                 , size_t const& i
                                 , size_t const& j
                                 , size_t const& k);
      /** @brief Compute advective term based on either Upwind or TVD spacial scheme */
      double compute_adv_component ( size_t const& comp
                                   , size_t const& i
                                   , size_t const& j
                                   , size_t const& k);
      /** @brief Assemble rhs term after calling compute_**_component */
      void assemble_DS_un_at_rhs (FV_TimeIterator const* t_it
                                , double const& gamma);

      void assemble_velocity_diffusion_terms ( );

      void calculate_row_indexes ( FV_DiscreteField const* FF);

      double divergence_of_U( size_t const& i
                            , size_t const& j
                            , size_t const& k
                            , size_t const& component
                            , size_t const& level);

      /** @brief Call functions to assemble rhs for pressure or velocity fields in any direction */
      double assemble_local_rhs( size_t const& j
                               , size_t const& k
                               , double const& gamma
                               , FV_TimeIterator const* t_it
                               , size_t const& comp
                               , size_t const& dir
                               , size_t const& field);

      /** @brief Assemble rhs for pressure in any direction */
      double calculate_velocity_divergence ( size_t const& i
                                           , size_t const& j
                                           , size_t const& k
                                           , size_t const& level
                                           , FV_TimeIterator const* t_it);

      double assemble_velocity_gradients (class doubleVector& grad
                                        , size_t const& i
                                        , size_t const& j
                                        , size_t const& k
                                        , size_t const& level);

      void detect_fresh_cells_and_neighbours();

      void calculate_divergence_weighting (FV_TimeIterator const* t_it);

      double pressure_local_rhs( size_t const& j
                               , size_t const& k
                               , FV_TimeIterator const* t_it
                               , size_t const& dir );

      /** @brief Initialize the velocity on the velocity nodes in MAC grid*/
      void initialize_grid_nodes_on_rigidbody( vector<size_t> const& list );

      void ugradu_initialization ( );

      void correct_pressure_1st_layer_solid (size_t const& level );

      void correct_pressure_2nd_layer_solid (size_t const& level );

      void correct_mean_pressure (size_t const& level );

      /** @brief Solve interface unknowns for
      both fields in any particular direction */
      void solve_interface_unknowns( FV_DiscreteField* FF
                                   , double const& gamma
                                   , FV_TimeIterator const* t_it
                                   , size_t const& comp
                                   , size_t const& dir);
      /** @brief Unpack the interface variable sent
      by master processor to slave processor */
      void unpack_ue(size_t const& comp
                   , double * received_data
                   , size_t const& dir
                   , size_t const& p
                   , size_t const& field);
      /** @brief Unpack the data sent by "data_packing" and compute
      the interface unknown; and pack ue for sending to slave processor */
      void unpack_compute_ue_pack(size_t const& comp
                                , size_t const& dir
                                , size_t const& p
                                , size_t const& field);

      /** @brief Pack Aei*(Aii)-1*fi and fe for sending to master processor */
      void data_packing ( FV_DiscreteField const* FF
                        , size_t const& j
                        , size_t const& k
                        , double const& fe
                        , size_t const& comp
                        , size_t const& dir );
      /** @brief Compute Aei*(Aii)-1*fi required to compute interface unknown */
      void compute_Aei_ui (struct TDMatrix* arr
                         , struct LocalVector* VEC
                         , size_t const& comp
                         , size_t const& dir
                         , size_t const& r_index);

      /** @brief Solve interface unknown for all cases */
      void DS_interface_unknown_solver( LA_SeqVector* rhs
                                      , size_t const& comp
                                      , size_t const& dir
                                      , size_t const& field
                                      , size_t const& r_index ) ;

      /** @brief Solve i in j and k; e.g. solve x in y ank z */
      void Solve_i_in_jk ( FV_DiscreteField* FF
                         , FV_TimeIterator const* t_it
                         , size_t const& dir_i
                         , size_t const& dir_j
                         , size_t const& dir_k
                         , double const& gamma );

      /** Pressure predictor */
      void NS_first_step( FV_TimeIterator const* t_it ) ;

      /** Velocity update */
      void NS_velocity_update( FV_TimeIterator const* t_it ) ;

      /** Pressure update (penalty step) */
      void NS_pressure_update( FV_TimeIterator const* t_it ) ;

      /** Pressure correction */
      void NS_final_step( FV_TimeIterator const* t_it ) ;

      void write_output_field( FV_DiscreteField const* FF) ;

      void output_L2norm_divergence( );

      void output_L2norm_pressure( size_t const& level );

      void output_L2norm_velocity( size_t const& level );
      //@}


   // Direction splitting communicators

      /** @name Direction splitting communicators */
      /** @brief Create the sub-communicators */
      void create_DS_subcommunicators ( void ) ;
      void processor_splitting ( int const& color
                               , int const& key
                               , size_t const& dir );

      void allocate_mpi_variables (FV_DiscreteField const* FF);
      void deallocate_mpi_variables ();

      /** @brief Free the sub-communicators */
      void free_DS_subcommunicators ( void ) ;


      //@}

      //-- Projection-Translation methods

      /** @name Projection-Translation methods */
      //@{
      /** @brief Set translation vector and direction */
      void set_translation_vector();

      /** @brief Build the field projection-translation interpolations */
      void build_links_translation();

      /** @brief Projection of the field on the translated position of the grid
      */
      void fields_projection();

      /** @brief Synchronize the global pressure variable at "level"
      */
      void synchronize_pressure_field ( size_t level );

      /** @brief Synchronize the global velocity variable at "level"
      */
      void synchronize_velocity_field ( size_t level );


   private: //----------------------------------------------------------------

   //-- Class attributes

   //-- Attributes

      FV_DiscreteField* UF;
      FV_DiscreteField* PF;

      DS_NavierStokesSystem* GLOBAL_EQ ;

      size_t nb_procs;
      size_t my_rank;
      size_t is_master;
      size_t dim;
      // 0th element for P and 1st element for U
      size_t nb_comps[2];

      MAC_Communicator const* macCOMM;
      MPI_Comm DS_Comm_i[3];

      int rank_in_i[3];
      int nb_ranks_comm_i[3];

      // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct MPIVarNS first_pass[2][3];
      // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct MPIVarNS second_pass[2][3];
      // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct MPIVarNS data_for_S[2][3];

      double mu;
      double kai;
      string AdvectionScheme;
      size_t AdvectionTimeAccuracy;
      double rho;
      bool b_restart ;
      bool is_solids;

      bool is_stressCal;
      string ViscousStressOrder;
      double surface_cell_scale;
      bool is_surfacestressOUT;
      size_t stressCalFreq;
      bool is_par_motion;

      // Grains3D variable
      DS_AllRigidBodies* allrigidbodies;

      // Grid motion
      bool b_projection_translation;
      bool b_grid_has_been_translated_since_last_output;
      bool b_grid_has_been_translated_at_previous_time;
      double critical_distance_translation;
      geomVector MVQ_translation_vector;
      size_t translation_direction;
      double bottom_coordinate;
      double translated_distance;

      boolVector const* P_periodic_comp;
      boolVector const* U_periodic_comp;
      MAC_DoubleVector* gravity_vector ;
      bool is_periodic[2][3];

} ;

#endif
