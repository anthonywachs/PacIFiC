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

/** @brief MPIVar include all vectors required while message passing */
struct MPIVar {
   int *size;
   double ***send;
   double ***receive;
};

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

      /** @name Basic discrete system building */
      //@{
      /** @brief Error compared to analytical solution */
      void error_with_analytical_solution ( ) ;

      /** @brief Call the function to assemble 1D matrices for both velocity and pressure field*/
      void assemble_1D_matrices( FV_TimeIterator const* t_it ) ;

      /** @brief Assemble 1D matrices for both velocity and pressure field in all directions */
      double assemble_field_matrix (
        FV_DiscreteField const* FF,
        FV_TimeIterator const* t_it,
        double const& gamma,
        size_t const& comp,
        size_t const& dir, size_t const& field );

      /** @brief Assemble 1D schur matrices for both velocity and pressure field in all directions */
      void assemble_field_schur_matrix (struct TDMatrix *A, size_t const& comp, size_t const& dir, double const& Aee_diagcoef, size_t const& field );

      /** @brief Assemble advection term for Upwind spacial scheme */
      double assemble_advection_Upwind( size_t const& advecting_level, double const& coef, 
         size_t const& advected_level,
         size_t const& i, size_t const& j, size_t const& k, size_t const& component) const;

      /** @brief Assemble advection term for TVD spacial scheme */
      double assemble_advection_TVD( size_t const& advecting_level, double const& coef, size_t const& advected_level,
         size_t const& i, size_t const& j, size_t const& k, size_t const& component) const;

      /** @brief Assemble rhs for velocity in any direction */
      double velocity_local_rhs( size_t const& j, size_t const& k, double const&  gamma, FV_TimeIterator const* t_it, size_t const& comp, size_t const& dir );
      /** @brief Compute diffusive term of velocity field from previous timestep */
      double compute_un_component ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k, size_t const& dir);
      /** @brief Compute diffusive term of pressure field from previous timestep */
      double compute_p_component ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k);
      /** @brief Compute advective term based on either Upwind or TVD spacial scheme */
      double compute_adv_component ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k);
      /** @brief Assemble rhs term after calling compute_**_component */
      void assemble_DS_un_at_rhs (FV_TimeIterator const* t_it, double const& gamma);

      /** @brief Call functions to assemble rhs for pressure or velocity fields in any direction */
      double assemble_local_rhs(size_t const& j,size_t const& k,double const& gamma,FV_TimeIterator const* t_it,size_t const& comp,size_t const& dir,size_t const& field);
      
      /** @brief Assemble rhs for pressure in any direction */
      double pressure_local_rhs( size_t const& j, size_t const& k, FV_TimeIterator const* t_it, size_t const& dir );

      /** @brief Solve interface unknowns for both fields in any particular direction */
      void solve_interface_unknowns( FV_DiscreteField* FF, double const& gamma, FV_TimeIterator const* t_it, size_t const& comp, size_t const& dir, size_t const& field );
      /** @brief Unpack the interface variable sent by master processor to slave processor */
      void unpack_ue(size_t const& comp, double * received_data, size_t const& dir, int const& p, size_t const& field);
      /** @brief Unpack the data sent by "data_packing" and compute the interface unknown; and pack ue for sending to slave processor */
      void unpack_compute_ue_pack(size_t const& comp, size_t const& dir, size_t const& p, size_t const& field);

      /** @brief Pack Aei*(Aii)-1*fi and fe for sending to master processor */
      void data_packing ( FV_DiscreteField const* FF, size_t const& j, size_t const& k, double const& fe, size_t const& comp, size_t const& dir, size_t const& field);
      /** @brief Compute Aei*(Aii)-1*fi required to compute interface unknown */
      void compute_Aei_ui (struct TDMatrix* arr, struct LocalVector* VEC, size_t const& comp, size_t const& dir);

      /** @brief Solve interface unknown for all cases */
      void DS_interface_unknown_solver( LA_SeqVector* rhs, size_t const& comp, size_t const& dir, size_t const& field ) ;
      
      /** @brief Solve i in j and k; e.g. solve x in y ank z */
      void Solve_i_in_jk ( FV_DiscreteField* FF, FV_TimeIterator const* t_it, size_t const& dir_i, size_t const& dir_j, size_t const& dir_k, double const& gamma,size_t const& field );

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

      double get_velocity_divergence(void);
      
      void output_L2norm_pressure( size_t const& level );

      void output_L2norm_velocity( size_t const& level );
      //@}


   // Direction splitting communicators
     
      /** @name Direction splitting communicators */
      /** @brief Create the sub-communicators */
      void create_DDS_subcommunicators ( void ) ;
      void processor_splitting ( int const& color, int const& key, size_t const& dir );

      void allocate_mpi_variables (FV_DiscreteField const* FF, size_t const& field);
      void deallocate_mpi_variables (size_t const& field);

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
      size_t nb_comps[2];               // 0th element for P and 1st element for U

      MAC_Communicator const* pelCOMM;
      MPI_Comm DDS_Comm_i[3];

      int rank_in_i[3];
      int nb_ranks_comm_i[3];

      struct MPIVar first_pass[2][3];           // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
      struct MPIVar second_pass[2][3];          // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions


      double peclet ;
      double rho;
      double mu;
      double kai;
      string AdvectionScheme;
      size_t AdvectionTimeAccuracy;

      bool b_restart ;
      bool is_firstorder ;

      boolVector const* P_periodic_comp;
      boolVector const* U_periodic_comp;
      bool is_periodic[2][3];
} ;

#endif
