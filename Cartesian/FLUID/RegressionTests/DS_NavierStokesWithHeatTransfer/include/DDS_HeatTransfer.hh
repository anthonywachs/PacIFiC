#ifndef DDS_HeatTransfer_HH
#define DDS_HeatTransfer_HH

#include <mpi.h>
#include <FV_OneStepIteration.hh>
#include <geomVector.hh>
#include <computingtime.hh>
#include <solvercomputingtime.hh>
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
class DDS_HeatTransferSystem ;
class LA_SeqMatrix ;

/** @brief The Class DDS_HeatTransfer.

Server for the resolution of the unsteady heat equation by a first order
implicit time integrator and a Finite Volume MAC scheme on rectangular grids.

Equation: rho*cp*dT/dt + rho*cp*(u.nabla(T))= k * lap(T) + bodyterm

@author A. Wachs - Pacific project 2017 */

struct NavierStokes2Temperature
{
  double rho_ ;
  bool b_restart_ ;
  FV_DomainAndFields const* dom_ ;
  FV_DiscreteField const* UF_ ;
  string AdvectionScheme_ ;
  size_t AdvectionTimeAccuracy_ ;
  bool is_solids_ ;
  size_t Npart_ ;
  string solid_filename_ ;
  double loc_thres_ ; 
  string level_set_type_ ;
  bool is_stressCal_ ;
  double Npoints_ ;
  double ar_ ;
  size_t Nrings_ ;
  size_t Pmin_ ;
  size_t pole_loc_ ;
};

/** @brief MPIVar include all vectors required while message passing */
struct MPIVar {
   int *size;
   double ***send;
   double ***receive;
};

class DDS_HeatTransfer : public MAC_Object, public ComputingTime, public SolverComputingTime
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      /** @brief Create and initialize an instance of REG_HeatTransfer
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param fromNS structure containing input data from the NS solver */
      static DDS_HeatTransfer* create(
		MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp,
                struct NavierStokes2Temperature const& fromNS );

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
      virtual ~DDS_HeatTransfer( void ) ;

      /** @brief Copy constructor */
      DDS_HeatTransfer( DDS_HeatTransfer const& other ) ;

      /** @brief Operator ==
      @param other the right hand side */
      DDS_HeatTransfer& operator=( DDS_HeatTransfer const& other ) ;

      /** @brief Constructor with arguments
      @param a_owner the MAC-based object
      @param exp to read the data file */
      DDS_HeatTransfer( MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp,
                struct NavierStokes2Temperature const& fromNS );

      /** @brief Constructor without argument */
      DDS_HeatTransfer( void ) ;

      //@}


   //-- Basic discrete system building

      /** @name Basic discrete system building */
      //@{
      /** @brief Compute diffusive term of temperature from perivious timestep */
      double compute_un_component ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k, size_t const& dir, size_t const& level);
      /** @brief Assemble temperature body term */
      double bodyterm_value ( double const& xC, double const& yC, double const& zC);

      /** @brief Assemble RHS of temperature for first step of Crank_Nicolson time discretization */
      void assemble_DS_un_at_rhs ( FV_TimeIterator const* t_it, double const& gamma);



      /** @brief Call the functions to assemble temperature and schur complement */
      void assemble_temperature_and_schur( FV_TimeIterator const* t_it) ;

      size_t return_row_index ( FV_DiscreteField const* FF, size_t const& comp, size_t const& dir, size_t const& j, size_t const& k );

      /** @brief Returns the node index required in the presence of solids in the system */
      size_t return_node_index ( FV_DiscreteField const* FF, size_t const& comp, size_t const& i, size_t const& j, size_t const& k );


      /** @brief Assemble temperature matrix */
      double assemble_temperature_matrix (
        FV_DiscreteField const* FF,
        FV_TimeIterator const* t_it,
        double const& gamma,
        size_t const& comp, 
        size_t const& dir,
        size_t const& j,
        size_t const& k,
        size_t const& r_index  );

      /** @brief Assemble schur matrix */
      void assemble_schur_matrix (size_t const& comp, size_t const& dir, double const& Aee_diagcoef, size_t const& r_index);


      void write_output_field();
      //@}

   //-- Solver

      /** @name Solvers */
      //@{

      /** @brief Assemble local RHS for 1D equation solver */
      double assemble_local_rhs( size_t const& j, size_t const& k, double const& gamma, FV_TimeIterator const* t_it, size_t const& comp, size_t const& dir );

      /** @brief Compute Aei*(Aii)-1*fi required to compute interface unknown */
      void compute_Aei_ui (struct TDMatrix* arr, struct LocalVector* VEC, size_t const& comp, size_t const& dir, size_t const& r_index);

      /** @brief Pack Aei*(Aii)-1*fi and fe for sending to master processor */ 
      void data_packing ( double const& fe, size_t const& comp, size_t const& dir, size_t const& vec_pos);

      /** @brief Unpack the data sent by "data_packing" and compute the interface unknown; and pack ue for sending to slave processor */ 
      void unpack_compute_ue_pack(size_t const& comp, size_t const& dir, size_t const& p);

      /** @brief Unpack the interface variable sent by master processor to slave processor */ 
      void unpack_ue(size_t const& comp, double * received_data, size_t const& dir, int const& p);

      /** @brief Call the appropriate functions to solve local variable and interface unknown */ 
      void solve_interface_unknowns(double const& gamma,FV_TimeIterator const* t_it, size_t const& comp, size_t const& dir);
      
      /** @brief Call the appropriate functions to solve any particular direction in the other directions */ 
      void HeatEquation_DirectionSplittingSolver( FV_TimeIterator const* t_it ) ;

      /** @brief Calculate the required field parrameters such as void_fractions, intersection pointes with the solid objects */ 
      void node_property_calculation( ) ;

      /** @brief Correct the fluxes and variables on the nodes due to presence of solid objects */ 
      void nodes_temperature_initialization ( size_t const& level );

      /** @brief Returns negative value of the point lies inside the solid, otherwise returns a positive number*/
      double level_set_function (size_t const& m, size_t const& comp, double const& xC, double const& yC, double const& zC, string const& type);

      /** @brief Correct the fluxes and variables on the nodes due to presence of solid objects */ 
      void assemble_intersection_matrix ( size_t const& comp, size_t const& level);                 // Here level:0 -> fluid; 1-> solid
      
      /** @brief Compute advective term based on either Upwind or TVD spacial scheme */
      double compute_adv_component ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k);

      /** @brief Find the intersection using bisection method with the solid interface */ 
      double find_intersection ( size_t const& left, size_t const& right, size_t const& yconst, size_t const& zconst, size_t const& comp, size_t const& dir, size_t const& off, size_t const& level);

      double find_intersection_for_ghost ( double const& xl, double const& xr, double const& yvalue, double const& zvalue, size_t const& id, size_t const& comp, size_t const& dir, double const& dx, size_t const& level, size_t const& off);

      double impose_solid_temperature_for_ghost (size_t const& comp, double const& xg, double const& yg, double const& zg, size_t const& parID );

      double impose_solid_temperature (size_t const& comp, size_t const& dir, size_t const& off, size_t const& i, size_t const& j, size_t const& k, double const& xb, size_t const& parID );

      double ghost_field_estimate_in_box ( FV_DiscreteField* FF, size_t const& comp, size_t const& i0, size_t const& j0, size_t const& k0, double const& x0, double const& y0, double const& z0, double const& dh, size_t const& level, size_t const& parID);

      double assemble_advection_Centered( FV_DiscreteField const* AdvectingField, size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const;
      double assemble_advection_Centered_new( FV_DiscreteField const* AdvectingField, size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const;

      double divergence_of_U ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k);


      double assemble_advection_TVD( FV_DiscreteField const* AdvectingField, size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const;

      double assemble_advection_Upwind( FV_DiscreteField const* AdvectingField, size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const;

      double assemble_advection_Upwind_new( FV_DiscreteField const* AdvectingField, size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const;

      /** @brief Generate the solid particles present in the domain */ 
      void Solids_generation( ) ;

      /** @brief Solve i in j and k; e.g. solve x in y ank z */ 
      void Solve_i_in_jk (FV_TimeIterator const* t_it, double const& gamma, size_t const& dir_i, size_t const& dir_j, size_t const& dir_k );

      /** @brief Solve interface unknown for all cases */ 
      void DS_interface_unknown_solver(LA_SeqVector* interface_rhs, size_t const& comp, size_t const& dir, size_t const& r_index);

      /** @brief Error compared to analytical solution */
      void DS_error_with_analytical_solution ( FV_DiscreteField const* FF,FV_DiscreteField* FF_ERROR ) ;
      
      void ugradu_initialization (  );
      /** @brief Calculate L2 norm without any analytical solution */
      void output_l2norm ( void ) ;
      //@}

      void compute_fluid_particle_interaction( FV_TimeIterator const* t_it, double const& Np);
      void generate_discretization_parameter(class doubleVector& eta, class doubleVector& k, class doubleVector& Rring, size_t const& k0, size_t const& Nring);
      void compute_surface_points(class doubleVector& eta, class doubleVector& k, class doubleVector& Rring, class doubleArray2D& point_coord, class doubleVector& cell_area, size_t const& Nring);
      void compute_temperature_gradient_on_particle(class doubleArray2D& point_coord, class doubleVector& cell_area, class doubleVector& force, size_t const& parID, size_t const& Np );
      double ghost_field_estimate_on_face ( FV_DiscreteField* FF, size_t const& comp, size_t const& i0, size_t const& j0, size_t const& k0, double const& x0, double const& y0, double const& z0, double const& dh, size_t const& face_vec, size_t const& level);


     

   // Direction splitting communicators
     
      /** @name Direction splitting communicators */
      /** @brief Create the sub-communicators */
      void create_DDS_subcommunicators ( void ) ;
      void processor_splitting ( int const& color, int const& key, size_t const& dir );

      void allocate_mpi_variables (void);
      void deallocate_mpi_variables ( void );


      /** @brief Free the sub-communicators */
      void free_DDS_subcommunicators ( void ) ;
      //@}

   private: //----------------------------------------------------------------

   //-- Class attributes

      static DDS_HeatTransfer const* PROTOTYPE ;

   //-- Attributes

      FV_DiscreteField* TF;
      FV_DiscreteField const* UF;

      FV_DiscreteField* TF_ERROR;

      FV_DiscreteField* TF_DS_ERROR;

      DDS_HeatTransferSystem* GLOBAL_EQ ;

      size_t nb_procs;
      size_t my_rank;
      size_t is_master;
      size_t dim;
      size_t nb_comps;
      size_t Npart;

      MAC_Communicator const* pelCOMM;
      MPI_Comm DDS_Comm_i[3];

      int rank_in_i[3];
      int nb_ranks_comm_i[3];

      struct MPIVar first_pass[3];
      struct MPIVar second_pass[3];
      
      double rho, heat_capacity, thermal_conductivity;
      double loc_thres; // Local threshold for the node near the solid interface to be considered inside the solid, i.e. local_CFL = loc_thres*global_CFL
      bool b_bodyterm ;
      bool is_solids;
      bool is_stressCal;
      double Npoints;
      size_t Nrings, Pmin, pole_loc;
      double ar;
      bool b_restart ;
      bool is_iperiodic[3];
      boolVector const* periodic_comp;
      string solid_filename;
      string AdvectionScheme;
      size_t AdvectionTimeAccuracy;
      string level_set_type;
} ;

#endif
