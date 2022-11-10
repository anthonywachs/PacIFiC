#ifndef _DS_ALLRIGIDBODIES__
#define _DS_ALLRIGIDBODIES__

#include <geomVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <doubleArray2D.hh>
#include <boolArray2D.hh>
#include <MAC_Communicator.hh>
#include <MAC_DoubleVector.hh>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using std::vector;
using std::istream ;
using std::ostream ;
using std::istringstream ;
using std::string;
class DS_RigidBody;
class FS_AllRigidBodies;
class FV_DiscreteField;
class FV_Mesh;
class FV_TimeIterator;

/** @brief Interpolation scheme enumeration */
enum scheme_list
{
  quadratic,
  linear01,
  linear12,
  linear0,
  linear1,
  linear2
};

/** @brief The class DS_AllRigidBodies.

The array of all rigid bodies in the Direction Splitting solver.


@author A. Goyal - Pacific project 2022 */

class DS_AllRigidBodies
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_AllRigidBodies();

      /** @brief Constructor with arguments
      @param dimens number of space dimensions
      @param in input stream where features of rigid bodies are read
      @param b_particles_as_fixed_obstacles treat all rigid bodies as fixed
      obstacles
      @param arb_UF Pointer to flow field UF
      @param arb_PF Pointer to flow field PF
      @param arb_scs scale of cell on the rigid body surface as
      compared with the cell of computational grid
      @param arb_macCOMM communicator for MPI communications */
      DS_AllRigidBodies( size_t& dimens
                       , istream& in
                       , bool const& b_particles_as_fixed_obstacles
                       , FV_DiscreteField const* arb_UF
                       , FV_DiscreteField const* arb_PF
                       , double const& arb_rho
                       , MAC_DoubleVector const* arb_gv
                       , double const& arb_scs
                       , MAC_Communicator const* arb_macCOMM
                       , double const& arb_mu
                       , bool const& is_GRAINS
                       , bool const& is_STL
                       , istream& STL_input );

      /** @brief Constructor with arguments
      @param dimens number of space dimensions
      @param in input stream where features of rigid bodies are read
      @param b_particles_as_fixed_obstacles treat all rigid bodies as fixed
      obstacles
      @param arb_UF Pointer to flow field UF
      @param arb_PF Pointer to flow field PF
      @param arb_TF Pointer to flow field TF
      @param arb_scs scale of cell on the rigid body surface as
      compared with the cell of computational grid
      @param arb_macCOMM communicator for MPI communications */
      DS_AllRigidBodies( size_t& dimens
                       , istream& in
                       , bool const& b_particles_as_fixed_obstacles
                       , FV_DiscreteField const* arb_UF
                       , FV_DiscreteField const* arb_PF
                       , FV_DiscreteField const* arb_TF
                       , double const& arb_scs
                       , MAC_Communicator const* arb_macCOMM
                       , double const& arb_mu
                       , double const& arb_RBTemp);

      /** @brief Constructor with arguments
      @param dimens number of space dimensions
      @param in input stream where features of rigid bodies are read
      @param b_particles_as_fixed_obstacles treat all rigid bodies as fixed
      obstacles
      @param arb_TF Pointer to flow field TF
      @param arb_scs scale of cell on the rigid body surface as
      compared with the cell of computational grid
      @param arb_macCOMM communicator for MPI communications */
      DS_AllRigidBodies( size_t& dimens
                       , istream& in
                       , bool const& b_particles_as_fixed_obstacles
                       , FV_DiscreteField const* arb_TF
                       , double const& arb_scs
                       , MAC_Communicator const* arb_macCOMM
                       , double const& arb_mu
                       , double const& arb_RBTemp);

      /** @brief Destructor */
      ~DS_AllRigidBodies();
      //@}


   //-- Get methods

      /**@name Get methods */
      //@{
      /** @brief Returns the total number of rigid bodies */
      size_t get_number_rigid_bodies() const;

      /** @brief Returns the number of particles */
      size_t get_number_particles() const;

      /** @brief Returns minimum coordinate of RB's dir component */
      double get_min_RB_coord( size_t const& dir );

      /** @brief Returns a constant pointer to the FS_AllRigidBodies object that
      contains the vector of all corresponding geometric rigid bodies */
      FS_AllRigidBodies const* get_ptr_FS_AllRigidBodies() const;

      /** @brief Returns a const pointer to the ith DS rigid body */
      DS_RigidBody const* get_ptr_rigid_body( size_t i ) const;

      /** @brief Returns a pointer to the ith DS rigid body */
      DS_RigidBody* get_ptr_rigid_body( size_t i );

      /** @brief Returns the void_fraction on field FF */
      size_t_vector* get_void_fraction_on_grid( FV_DiscreteField const* FF );

      /** @brief Returns the ID of rigid body present on the field FF */
      size_t_vector* get_rigidbodyIDs_on_grid( FV_DiscreteField const* FF );

      /** @brief Returns the ID of rigid body present on the field FF */
      size_t_array2D* get_intersect_vector_on_grid( FV_DiscreteField const* FF );

      /** @brief Returns the intersection distance with the rigid body for FF*/
      doubleArray2D* get_intersect_distance_on_grid( FV_DiscreteField const* FF );

      /** @brief Returns the Dirichlet BC on the near rigid body on FF */
      doubleArray2D* get_intersect_fieldValue_on_grid( FV_DiscreteField const* FF );

      //@}


   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief Updates all rigid bodies
      @param in input stream where features of rigid bodies are read */
      void update( istream& in );
      //@}


   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Writes the geometric attributes of the rigid bodies in a stream
      @param out output stream
      @param indent_width indentation width */
      void display_geometric( ostream& out, size_t const& indent_width ) const;

      /** @brief Writes the attributes of the rigid bodies in a stream
      @param out output stream
      @param indent_width indentation width */
      void display( ostream& out, size_t const& indent_width ) const;

      /** @brief Returns whether a point is inside a rigid body
      @param parID particle ID to check for isIn
      @param pt the point */
      bool isIn( size_t const& parID, geomVector const& pt ) const;

      /** @brief Returns the parID if the point is inside parID, based on isIn
      @param ownID ID of RB owning the pt
      @param pt the point */
      int isIn_any_RB( size_t const& ownID, geomVector const& pt ) const;

      /** @brief Returns the parID if the point is inside parID, based on isIn
      @param ownID ID of RB owning the pt(x,y,z)
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      int isIn_any_RB(  size_t const& ownID,
                        double const& x,
                        double const& y,
                        double const& z ) const;

      /** @brief Returns the parID if the point is inside parID, based on levelset
      @param pt the point */
      int levelset_any_RB( geomVector const& pt ) const;

      /** @brief Returns the parID if the point is inside parID, based on levelset
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      int levelset_any_RB( double const& x,
                           double const& y,
                           double const& z ) const;


      /** @brief Returns the parID if the point is inside parID, based on levelset
      @param ownID ID of RB owning the pt
      @param pt the point */
      int levelset_any_RB( size_t const& ownID, geomVector const& pt ) const;

      /** @brief Returns the parID if the point is inside parID, based on levelset
      @param ownID ID of RB owning the pt(x,y,z)
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      int levelset_any_RB(  size_t const& ownID,
                        double const& x,
                        double const& y,
                        double const& z ) const;

      /** @brief Returns whether a point is inside a rigid body
      @param parID particle ID to check for isIn
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      bool isIn( size_t const& parID,
                 double const& x,
                 double const& y,
                 double const& z ) const;

      /** @brief Returns the level set value of a point from a rigid body
      @param parID particle ID to check for isIn
      @param pt the point */
      double level_set_value( size_t const& parID, geomVector const& pt ) const;

      /** @brief Returns the level set value of a point from a rigid body
      @param parID particle ID to check for isIn
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      double level_set_value( size_t const& parID,
                              double const& x,
                              double const& y,
                              double const& z ) const;

      /** @brief Returns the length of side in fluid
      @param pt1 Point 1
      @param pt2 Point 2 */
      std::tuple<double, geomVector, int, int>
             return_side_fraction(geomVector const& pt1
                                , geomVector const& pt2 );

      /** @brief Returns the flux from RB given the PF node index
      @param i index i
      @param j index j
      @param k index k */
      double calculate_divergence_flux_fromRB ( size_t const& i,
                                                size_t const& j,
                                           		size_t const& k);

      /** @brief Calculate the divergenceredistribution factor for each PF node
      using the normal of RB intersect
      @param i index i
      @param j index j
      @param k index k */
      vector<double> flux_redistribution_factor ( size_t const& i,
                                                  size_t const& j,
                                                  size_t const& k);

      /** @brief Calculates the area of RB section with given vector
      of polygon corners
      @param points vector of points */
      double calculate_area_of_RBplane( vector<geomVector> const& points);
                                      // , geomVector const& normal);

      /** @brief Computes the halo zone for all rigid bodies, required for
      void fraction and intersection calculation on the grid nodes */
      void compute_halo_zones_for_all_rigid_body( );

      /** @brief Computes the void fraction on the grid nodes
      of a given fluid field
      @param FF the fluid field (PF, UF)
      @param is_in_time_iter true if method called in time iteration */
      void compute_void_fraction_on_grid( FV_DiscreteField const* FF
                                        , bool const& is_in_time_iter );

      /** @brief Computes the void fraction on the epsilon grid for PP
      of a given fluid field
      @param FF the fluid field (PP_EPSILON) */
      void compute_void_fraction_on_epsilon_grid( FV_DiscreteField * FF );

      /** @brief Computes the face fraction belonging the each node
      of a given fluid field */
      void compute_cutCell_geometric_parameters_PF();

      /** @brief Compute the velocity face flux from the UF face */
      double divergence_face_flux ( size_t const& p_PF
                                  , size_t const& i
                                  , size_t const& j
                                  , size_t const& k
										    , size_t const& comp);

      /** @brief Computes the intersection of grid nodes of a given fluid field
      with the nearest rigid body of a given fluid field
      @param FF the fluid field (PF, UF)
      @param is_in_time_iter true if method called in time iteration */
      void compute_grid_intersection_with_rigidbody(FV_DiscreteField const* FF
                                                , bool const& is_in_time_iter);

      /** @brief Clear the void fraction and intersection belonging to
      rigid bodies coming from GRAINS
      @param FF the fluid field (PF, UF) */
      void clear_GrainsRB_data_on_grid( FV_DiscreteField const* FF );

      /** @brief Computes the rigid body velocity including the rotation speed
      at a given geometric vector pt
      @param pt a point in space*/
      geomVector rigid_body_velocity( size_t const& parID,
                                          geomVector const& pt );

      /** @brief Computes the rigid body velocity translation velocity
      without the rotation speed
      @param parID rigid body ID */
      geomVector rigid_body_GC_velocity( size_t const& parID);

      /** @brief Returns pointer to gravity center of parID
      @param parID ID of rigid body*/
      geomVector const* get_gravity_centre( size_t const& parID);

      /** @brief Returns the rigid body angular velocity */
      geomVector rigid_body_angular_velocity( size_t const& parID) const;

      /** @brief Return the rigid body temperature, currently a user input
      but can be read from the Grains in future
      @param pt a point in space*/
      geomVector rigid_body_temperature( size_t const& parID,
                                          geomVector const& pt ) const;

      /** @brief Build the variable associated with the rigid bodies
      on the Cartesian computational grid
      @param FF target fluid grid */
      void build_solid_variables_on_fluid_grid( FV_DiscreteField const* FF );

      /** @brief Periodic transformation of a distance in a given direction
      @param delta distance
      @param dir direction  */
      double delta_periodic_transformation( double const& delta
                                          , size_t const& dir);

      /** @brief Periodic transformation of a point in a given direction
      @param x point
      @param dir direction  */
      double periodic_transformation( double const& x
                                    , size_t const& dir);

      /** @brief Intialize the surface variables for all rigid bodies */
      void initialize_surface_variables_for_all_RB( );

      /** @brief Intialize the surface variables on the Cartesian grid */
      void initialize_surface_variables_on_grid();

      /** @brief Compute the surface variables for all rigid bodies */
      void compute_surface_variables_for_all_RB( );

      /** @brief Write the surface variables for all rigid bodies */
      void write_surface_discretization_for_all_RB( );

      /** @brief Return the sum of interpolated field for all
      given list of levels on a point in 2D plane including
      the corrections near the solid interface
      @param FF field to interpolate
      @param comp component of the field to interpolate
      @param point coordinate where the field is interpolated
      @param face_vec Inplane vector of the 2D plane
      @param level vector of field levels to be interpolated   */
      double Bilinear_interpolation ( FV_DiscreteField const* FF
                                    , size_t const& comp
                                    , geomVector const* pt
                                    , size_t_vector const& i0
                                    , size_t_vector const& face_vec
                                    , vector<size_t> const& list);

      /** @brief Return the sum of interpolated field for all
      given list of levels on a point in 2D plane including
      the corrections near the solid interface
      @param FF field to interpolate
      @param comp component of the field to interpolate
      @param point coordinate where the field is interpolated
      @param i0 indexes of the point
      @param interpol_dir direction of the interpolation
      @param sign of the surface normal vector in each direction
      @param level vector of field levels to be interpolated   */
      double Biquadratic_interpolation ( FV_DiscreteField const* FF
                                       , size_t const& comp
                                       , geomVector const* pt
                                       , size_t_vector const& i0
                                       , size_t const& interpol_dir
                                       , int const& sign
                                       , vector<size_t> const& list);

      double Triquadratic_interpolation ( FV_DiscreteField const* FF
                                        , size_t const& comp
                                        , geomVector const* pt
                                        , size_t_vector const& i0
                                        , size_t const& parID
                                        , size_t const& ghost_points_dir
                                        , vector<int> const& sign
                                        , vector<size_t> const& list);

      /** @brief Return the sum of interpolated field for all
      given list of levels on a point in 3D box including
      the corrections near the solid interface
      @param FF field to interpolate
      @param comp component of the field to interpolate
      @param point coordinate where the field is interpolated
      @param level vector of field levels to be interpolated   */
      double Trilinear_interpolation ( FV_DiscreteField const* FF
                                     , size_t const& comp
                                     , geomVector const* pt
                                     , size_t_vector const& i0
                                     , size_t const& parID
                                     , vector<size_t> const& list);

      /** @brief Calculate the first order pressure force and torque on parID
      @param parID rigid body ID */
      void first_order_pressure_stress( size_t const& parID );

      /** @brief Calculate the first order viscous force and torque on parID
      @param parID rigid body ID */
      void first_order_viscous_stress( size_t const& parID );

      /** @brief Calculate the second order viscous force and torque on parID
      @param parID rigid body ID */
      void second_order_viscous_stress( size_t const& parID );

      /** @brief Calculate the first order temperature flux on parID
      @param parID rigid body ID */
      void first_order_temperature_flux( size_t const& parID );

      /** @brief Calculate the second order temperature flux on parID
      @param parID rigid body ID */
      void second_order_temperature_flux( size_t const& parID );

      /** @brief Calculate the pressure force and torque on all rigid bodies */
      void compute_pressure_force_and_torque_for_allRB ();

      /** @brief Calculate the viscous force and torque on all rigid bodies */
      void compute_viscous_force_and_torque_for_allRB (string const& StressOrder);

      /** @brief Calculate the temperature gradient all rigid bodies */
      void compute_temperature_gradient_for_allRB (string const& StressOrder);

      /** @brief Creates the neighbour list for each RB, helps in reducing
      computing time of isIn_any_RB */
      void create_neighbour_list_for_AllRB( );

      /** @brief Return numeric value with the provided field
      @param FF input field */
      int field_num( FV_DiscreteField const* FF );

      /** @brief Solve and update the particle position and velocity
      after solving RB equation of motion explicitly
      @param t_it Time_iterator */
      void solve_RB_equation_of_motion( FV_TimeIterator const* t_it);

      /** @brief Create a list of local rigid bodies */
      void generate_list_of_local_RB( );

      /** @brief Output force and flux summary in a csv file
      @param t_it Time_iterator
      @param b_restart if simulation is restart or not */
      void write_force_and_flux_summary( FV_TimeIterator const* t_it
                                       , bool const& b_restart);

      /** @brief Returns true if a box is within the local domain
      extents
      @param bounds vector with min and max in dir
      @param dir direction to check */
      bool is_bounding_box_in_local_domain( class doubleVector& bounds
                                           , size_t const& dir);
      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{

      //@}

   private: //----------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      size_t m_space_dimension; /**< Space dimension */
      size_t m_npart; /**< number of particles */
      size_t m_nrb; /**< total number of rigid bodies = number of
      	particles + number of obstacles, npart first rigid bodies are always
         particles while ( m_nrb - m_npart ) last rigid bodies are obstacles */
      vector<DS_RigidBody*> m_allDSrigidbodies; /**< the vector of all
    	Direction Splitting rigid bodies */
      FS_AllRigidBodies* m_FSallrigidbodies; /**< the pointer to the
    	FS_AllRigidBodies object that contains the vector of all
    	corresponding geometric rigid bodies */

      // Pointers to the constant fields and primary grid
      FV_DiscreteField const* UF ;
      FV_DiscreteField const* PF ;
      FV_DiscreteField const* TF ;
      FV_Mesh const* MESH ;

      double surface_cell_scale; /**< a variable to store the scale of surface
      cell on the RB as compared with computational grid cell size */

      vector<size_t_vector*> void_fraction; /**< vector of void fraction the
      field grid nodes, 0 in fluid and (parID+1) in the rigid bodies*/

      doubleVector* face_fraction; /**< vector of field grid nodes
      face fraction, only for UF*/

      intVector* face_ownerID; /**< vector of field grid nodes
      face fraction, only for UF*/

      intVector* cell_ownerID; /**< vector of field grid nodes
      face fraction, only for UF*/

      vector<geomVector> face_centroid; /**< centroid of the cut faces
      on UF field, if any*/

      vector<vector<geomVector>> intersect_points; /**< vector of set of
      intersection points on UF field grid faces*/

      vector<geomVector> normalRB; /**< vector of rigid body normal
      intersecting a pressure cell, if any*/

      vector<double> CutRBarea; /**< vector of area of section of RB intersect
      with a pressure cell, if any*/

      // Columns in each variable are (left,right,bottom,top,behind,front)
      vector<size_t_array2D*> intersect_vector;  /**<Direction of intersection*/
      vector<doubleArray2D*> intersect_distance; /**< Value of offset relative
      to node point */
      vector<doubleArray2D*> intersect_fieldValue; /**< Value of field variable
      at the intersection */
      doubleArray2D* pressure_force; /**< Value of force due to
      pressure stress on rigid bodies */
      doubleArray2D* viscous_force; /**< Value of force due to
      viscous stress on rigid bodies */
      doubleVector* temperature_gradient; /**< Value of temperature flux
      on rigid bodies */
      double avg_temperature_gradient; /**< Value of average temperature flux
      on all rigid bodies */
      geomVector avg_pressure_force; /**< Value of average force due to
      pressure stress on all rigid bodies */
      geomVector avg_viscous_force; /**< Value of average force due to
      viscous stress on all rigid bodies */
      doubleArray2D* pressure_torque; /**< Value of torque due to
      pressure stress on rigid bodies */
      doubleArray2D* viscous_torque; /**< Value of torque due to
      viscous stress on rigid bodies */
      geomVector avg_pressure_torque; /**< Value of average torque due to
      pressure stress on all rigid bodies */
      geomVector avg_viscous_torque; /**< Value of average torque due to
      viscous stress on all rigid bodies */
      MAC_Communicator const* m_macCOMM; /**< Variable for communication
      between processors */
      double m_mu; /**< Fluid viscosity */
      double m_rho; /**< fluid density */
      double m_RBTemp; /**< Rigid body temperature */
      MAC_DoubleVector const* gravity_vector ;
      vector<vector<size_t>> neighbour_list;

      vector<size_t> local_RB_list; /**< Stores a list of local RB */

      bool b_GRAINS;
      bool b_STL;

      //@}


    //-- Constructors & Destructor

       /**@name Constructors & Destructor */
       //@{
       /** @brief Copy constructor
       @param copy copied DS_AllRigidBodies object */
       DS_AllRigidBodies( DS_AllRigidBodies const& copy );
       //@}
};

#endif
