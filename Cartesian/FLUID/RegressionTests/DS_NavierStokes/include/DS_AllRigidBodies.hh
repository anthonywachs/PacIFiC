#ifndef _DS_ALLRIGIDBODIES__
#define _DS_ALLRIGIDBODIES__

#include <geomVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <doubleArray2D.hh>
#include <vector>
#include <iostream>
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

/** @brief The class DS_AllRigidBodies.

The array of all rigid bodies in the Direction Splitting solver.


@author A. Wachs - Pacific project 2021 */

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
      @param surface_cell_scale scale of cell on the rigid body surface as
      compared with the cell of computational grid */
      DS_AllRigidBodies( size_t& dimens
                       , istream& in
                       , bool const& b_particles_as_fixed_obstacles
                       , FV_DiscreteField const* arb_UF
                       , FV_DiscreteField const* arb_PF
                       , double const& arb_scs);

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


      /** @brief Computes the hydrodynamic force and torque acting on each rigid
      body and store the values in its corresponding geometric rigid body
      @param PP the pressure field
      @param UU the velocity field */
      void compute_hydro_force_torque( FV_DiscreteField const* PP,
	                                    FV_DiscreteField const* UU );

      /** @brief Computes the halo zone for all rigid bodies, required for
      void fraction and intersection calculation on the grid nodes */
      void compute_halo_zones_for_all_rigid_body( );

      /** @brief Computes the void fraction on the grid nodes
      of a given fluid field */
      void compute_void_fraction_on_grid( );

      /** @brief Computes the intersection of grid nodes of a given fluid field
      with the nearest rigid body of a given fluid field */
      void compute_grid_intersection_with_rigidbody( );

      /** @brief Computes the rigid body velocity including the rotation speed
      at a given geometric vector pt
      @param pt a point in space*/
      geomVector rigid_body_velocity( size_t const& parID,
                                          geomVector const& pt ) const;

      /** @brief Build the variable associated with the rigid bodies
      on the Cartesian computational grid */
      void build_solid_variables_on_grid( );

      /** @brief Intialize the surface variables for each rigid body */
      void initialize_surface_variables_for_all_RB( );

      /** @brief Compute the surface variables for each rigid body */
      void compute_surface_variables_for_all_RB( );
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

      // Pointers to the constant velocity and pressure field
      FV_DiscreteField const* UF ;
      FV_DiscreteField const* PF ;

      double surface_cell_scale; /**< a variable to store the scale of surface
      cell on the RB as compared with computational grid cell size */

      vector<size_t_vector*> void_fraction; /**< vector of void fraction the
      field grid nodes */
      vector<size_t_vector*> rb_ID; /**< vector of rigid body ID on the
      field grid node, if any */

      vector<struct BoundaryBisec*> rb_intersect; /**< 2DArray of intersection
      of field grid node near the rigid body with the rigid */

      // Columns in each variable are (left,right,bottom,top,behind,front)
      vector<size_t_array2D*> intersect_vector;  /**<Direction of intersection*/
      vector<doubleArray2D*> intersect_distance; /**< Value of offset relative
      to node point */
      vector<doubleArray2D*> intersect_fieldValue; /**< Value of field variable
      at the intersection */

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
