#ifndef DDS_NavierStokes_HH
#define DDS_NavierStokes_HH

#include <mpi.h>
#include <FV_OneStepIteration.hh>
#include <geomVector.hh>
#include <computingtime.hh>
#include <boolVector.hh>
#include <solvercomputingtime.hh>
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
class DDS_NavierStokesSystem ;
class LA_SeqMatrix ;

/** @brief The Class DDS_NavierStokes.

Server for the resolution of the unsteady heat equation by a first order
implicit time integrator and a Finite Volume MAC scheme on rectangular grids.

Equation: dT/dt = ( 1 / Pe ) * lap(T) + bodyterm, where Pe is the Peclet number.

@author A. Wachs - Pacific project 2017 */

class DDS_NavierStokes : public FV_OneStepIteration, public ComputingTime,
public SolverComputingTime
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

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
      ~DDS_NavierStokes( void ) ;

      /** @brief Copy constructor */
      DDS_NavierStokes( DDS_NavierStokes const& other ) ;

      /** @brief Operator ==
      @param other the right hand side */
      DDS_NavierStokes& operator=( DDS_NavierStokes const& other ) ;

      /** @brief Constructor with arguments
      @param a_owner the MAC-based object
      @param exp to read the data file */
      DDS_NavierStokes( MAC_Object* a_owner,
      		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) ;

      /** @brief Constructor without argument */
      DDS_NavierStokes( void ) ;

      /** @brief Create a clone
      @param a_owner the MAC-based object
      @param dom mesh and fields
      @param prms set of parameters
      @param exp to read the data file */
      virtual DDS_NavierStokes* create_replica(
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const ;
      //@}


   //-- Basic discrete system building

      void allocate_mpi_vectors_U(void);

      void allocate_mpi_vectors_P(void);

      void deallocate_mpi_vectors_U(void);

      void deallocate_mpi_vectors_P(void);

      /** @name Basic discrete system building */
      //@{
      /** @brief Assemble velocity body term */
      void assemble_velocity_bodyterm_rhs(
      	FV_DiscreteField const* FF,
	LA_Vector* VEC_rhs ) ;

      /** @brief Error compared to analytical solution */
      void error_with_analytical_solution ( FV_DiscreteField const* FF,
      	 FV_DiscreteField* FF_ERROR ) ;

      /** @brief Assemble 1D velocity matrices */
      void assemble_velocity_1D_matrices( FV_TimeIterator const* t_it ) ;

      /** @brief Assemble 1D velocity matrix in x */
      void assemble_velocity_matrix_1D_x (
        FV_DiscreteField const* FF,
        FV_TimeIterator const* t_it,
        double gamma,
        size_t const& comp  );

      /** @brief Assemble 1D velocity matrix in y */
      void assemble_velocity_matrix_1D_y (
        FV_DiscreteField const* FF,
        FV_TimeIterator const* t_it,
        double gamma,
        size_t const& comp  );

      /** @brief Assemble 1D velocity matrix in z */
      void assemble_velocity_matrix_1D_z (
        FV_DiscreteField const* FF,
        FV_TimeIterator const* t_it,
        double gamma,
        size_t const& comp   );

      /** @brief Assemble 1D pressure matrix in x */
      void assemble_pressure_matrix_1D_x (
        FV_DiscreteField const* FF,
        FV_TimeIterator const* t_it);

      /** @brief Assemble 1D pressure matrix in y */
      void assemble_pressure_matrix_1D_y (
        FV_DiscreteField const* FF,
        FV_TimeIterator const* t_it);

      /** @brief Assemble 1D pressure matrix in z */
      void assemble_pressure_matrix_1D_z (
        FV_DiscreteField const* FF,
        FV_TimeIterator const* t_it);

      double assemble_advection_Upwind( size_t advecting_level, double const& coef, 
         size_t advected_level,
         size_t const& i, size_t const& j, size_t const& k, size_t const& component) const;

      double assemble_advection_TVD( size_t advecting_level, double const& coef, size_t advected_level,
         size_t const& i, size_t const& j, size_t const& k, size_t const& component) const;


      //@}

   //-- Solver

      /** @name Solvers */
      //@{
      
      /** @brief Second order Direction splitting with Domain Decomposition of Velocity Update step in Navier Stokes solver */

      /** @brief Assemble rhs for velocity in x */
      double assemble_local_rhs_x( size_t const& j, size_t const& k, double gamma, FV_TimeIterator const* t_it, size_t const& comp  );
      /** @brief Solve local unknowns for velocity in x */
      void solve_x_for_secondorder ( size_t const& j, size_t const& k, double gamma, FV_TimeIterator const* t_it, double * packed_data, size_t const& comp );
      /** @brief Solve interface unknowns for velocity in x */
      void solve_interface_unknowns_x(  double * packed_data, size_t nb_send_data, double gamma,  FV_TimeIterator const* t_it, size_t const& comp );

      /** @brief Assemble rhs for velocity in y */
      double assemble_local_rhs_y( size_t const& i, size_t const& k, double gamma, FV_TimeIterator const* t_it, size_t const& comp   );
      /** @brief Solve local unknowns for velocity in y */
      void solve_y_for_secondorder ( size_t const& i, size_t const& k, double gamma, FV_TimeIterator const* t_it, double * packed_data, size_t const& comp  );
      /** @brief Solve interface unknowns for velocity in y */
      void solve_interface_unknowns_y(  double * packed_data, size_t nb_send_data, double gamma,  FV_TimeIterator const* t_it, size_t const& comp   );

      /** @brief Assemble rhs for velocity in z */
      double assemble_local_rhs_z( size_t const& i, size_t const& j, double gamma, FV_TimeIterator const* t_it , size_t const& comp  );
      /** @brief Solve local unknowns for velocity in z */
      void solve_z_for_secondorder ( size_t const& i, size_t const& j, double gamma, FV_TimeIterator const* t_it, double * packed_data, size_t const& comp  );
      /** @brief Solve interface unknowns for velocity in z */
      void solve_interface_unknowns_z(  double * packed_data, size_t nb_send_data, double gamma,  FV_TimeIterator const* t_it, size_t const& comp   );
     

      /** @brief Second order Direction splitting with Domain Decomposition of Pressure Update step in Navier Stokes solver */

      /** @brief Assemble rhs for pressure in x */
      double pressure_assemble_local_rhs_x( size_t const& j, size_t const& k, FV_TimeIterator const* t_it );
      /** @brief Solve local unknowns for pressure in x */
      void pressure_solve_x_for_secondorder ( size_t const& j, size_t const& k, FV_TimeIterator const* t_it, double * packed_data );
      /** @brief Solve interface unknowns for pressure in x */
      void pressure_solve_interface_unknowns_x(  double * packed_data, size_t nb_send_data, FV_TimeIterator const* t_it );

      /** @brief Assemble rhs for pressure in y */
      double pressure_assemble_local_rhs_y( size_t const& i, size_t const& k, FV_TimeIterator const* t_it   );
      /** @brief Solve local unknowns for velocity in y */
      void pressure_solve_y_for_secondorder ( size_t const& i, size_t const& k, FV_TimeIterator const* t_it, double * packed_data );
      /** @brief Solve interface unknowns for velocity in y */
      void pressure_solve_interface_unknowns_y(  double * packed_data, size_t nb_send_data, FV_TimeIterator const* t_it   );

      /** @brief Assemble rhs for pressure in z */
      double pressure_assemble_local_rhs_z( size_t const& i, size_t const& j, FV_TimeIterator const* t_it );
      /** @brief Solve local unknowns for pressure in z */
      void pressure_solve_z_for_secondorder ( size_t const& i, size_t const& j, FV_TimeIterator const* t_it, double * packed_data );
      /** @brief Solve interface unknowns for pressure in z */
      void pressure_solve_interface_unknowns_z(  double * packed_data, size_t nb_send_data, FV_TimeIterator const* t_it );

      
       /** @brief Navier Stokes solver */

      /** Pressure predictor */

      void NS_first_step( FV_TimeIterator const* t_it ) ;

      /** Velocity update */
      void NS_velocity_update( FV_TimeIterator const* t_it ) ;
      
      /** Pressure update (penalty step) */
      void NS_pressure_update( FV_TimeIterator const* t_it ) ;
      
      /** Pressure correction */
      void NS_final_step( FV_TimeIterator const* t_it ) ;

      void write_pressure_field( FV_TimeIterator const* t_it );

      void write_velocity_field( FV_TimeIterator const* t_it );

      void get_velocity_divergence(void);
      
      void output_L2norm_pressure( size_t level );

      void output_L2norm_velocity( size_t level );

      void compute_CFL( FV_TimeIterator const* t_it, size_t level ) const;
      //@}


   // Direction splitting communicators
     
      /** @name Direction splitting communicators */
      /** @brief Create the sub-communicators */
      void create_DDS_subcommunicators ( void ) ;

      /** @brief Free the sub-communicators */
      void free_DDS_subcommunicators ( void ) ;


      //@}

   private: //----------------------------------------------------------------

   //-- Class attributes

      static DDS_NavierStokes const* PROTOTYPE ;

   //-- Attributes

      FV_DiscreteField* UF;

      FV_DiscreteField* PF;

      DDS_NavierStokesSystem* GLOBAL_EQ ;

      size_t nb_procs;
      size_t my_rank;
      size_t is_master;
      size_t dim;
      size_t nb_comps;

      MAC_Communicator const* pelCOMM;
      MPI_Comm DDS_Comm_x;
      MPI_Comm DDS_Comm_y;
      MPI_Comm DDS_Comm_z;

      int nb_ranks_comm_x;
      int rank_in_x;
      int nb_ranks_comm_y;
      int rank_in_y;
      int nb_ranks_comm_z;
      int rank_in_z;
      
      //MPI_Comm DDS_Comm_z;

      double peclet ;
      double rho;
      double mu;
      double kai;
      string AdvectionScheme;
      size_t AdvectionTimeAccuracy;

      bool b_restart ;

      bool U_is_xperiodic;
      bool U_is_yperiodic;
      bool U_is_zperiodic;
      boolVector const* U_periodic_comp;

      bool P_is_xperiodic;
      bool P_is_yperiodic;
      bool P_is_zperiodic;
      boolVector const* P_periodic_comp;

      // double * schlur_mpi_data_u_x;
      // double * schlur_mpi_data_u_y;
      // double * schlur_mpi_data_u_z;

      double ** mpi_packed_data_U_x;
      double ** mpi_packed_data_U_y;
      double ** mpi_packed_data_U_z;

      // MPI vectors on the master proc
      double *** all_receive_data_U_x;
      double *** all_receive_data_U_y;
      double *** all_receive_data_U_z;

      double *** all_send_data_U_x;
      double *** all_send_data_U_y;
      double *** all_send_data_U_z;

      double * mpi_packed_data_P_x;
      double * mpi_packed_data_P_y;
      double * mpi_packed_data_P_z;

      double ** all_receive_data_P_x;
      double ** all_receive_data_P_y;
      double ** all_receive_data_P_z;

      double ** all_send_data_P_x;
      double ** all_send_data_P_y;
      double ** all_send_data_P_z;
} ;

#endif
