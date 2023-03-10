#ifndef _DS_RIGIDBODY__
#define _DS_RIGIDBODY__

#include <geomVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <doubleArray2D.hh>
#include <boolVector.hh>
#include <MAC_assertions.hh>
#include <MAC_Communicator.hh>
#include <vector>
#include <iostream>
#include <map>
#include <iostream>
#include <fstream>
using std::ostream;
using std::vector;
using std::tuple;
class FS_RigidBody;
class FV_DiscreteField;
class FV_Mesh;

/** @brief The class DS_RigidBody.

A moving or stationary rigid body in the Direction Splitting solver.

@author A. Wachs - Pacific project 2021 */

class DS_RigidBody
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_RigidBody();

      /** @brief Constructor with arguments
      @param pgrb pointer to the corresponding geometric rigid body */
      DS_RigidBody( FS_RigidBody* pgrb );

      /** @brief Destructor */
      virtual ~DS_RigidBody();
      //@}


   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief Updates the rigid body features from its corresponding
      geometric rigid body */
      virtual void update() = 0;
      //@}

   //-- Get Methods
      /**@name Get methods */
      //@{
      /** @brief Return the halo zones of the rigid body */
      vector<geomVector*> get_rigid_body_haloZone( ) const;

      /** @brief Returns rigid body velocity including rotation speed at pt
      @param pt the point */
      virtual geomVector get_rigid_body_velocity( geomVector const& pt )
                                                      const = 0;

      /** @brief Returns rigid body angular velocity */
      virtual geomVector get_rigid_body_angular_velocity( ) const = 0;

      /** @brief Returns pointer to the rigid body gravity center */
      virtual geomVector const* get_ptr_to_gravity_centre( ) const = 0;

      /** @brief Returns circumscribed radius */
      virtual double get_circumscribed_radius( ) const = 0;

      /** @brief Returns a tuple of mass and density of RB */
      virtual std::tuple<double,double,double> get_mass_and_density_and_moi() const = 0;

      /** @brief Returns the distance of a point with the rigid body
      with a given ray vector and source
      @param pt the point
      @param direction x, y or z (=0, 1 or 2)
      @param positive true if search in the positive direction of the coordinate
      axis and false otherwise */
      virtual double get_distanceTo( geomVector const& source,
                                     geomVector const& rayDir,
                                     double const& delta ) const = 0;

      /** @brief Return the surface points on the rigid body */
      vector<geomVector*> get_rigid_body_surface_points( ) const;

      /** @brief Return the area of the surface points on the rigid body */
      vector<geomVector*> get_rigid_body_surface_areas( ) const;

      /** @brief Return the normal of the surface points on the rigid body */
      vector<geomVector*> get_rigid_body_surface_normals( ) const;

      //@}

   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Writes the attributes in a stream
      @param out output stream
      @param indent_width indentation width */
      virtual void display( ostream& out, size_t const& indent_width )
      	const = 0;

      /** @brief Compute the halozone of a rigid body
      @param out output stream
      @param indent_width indentation width */
      virtual void compute_rigid_body_halozone( double const& dx ) = 0;

      /** @brief Returns whether a point is inside the rigid body
      @param pt the point */
      virtual bool isIn( geomVector const& pt ) const = 0;

      /** @brief Returns whether a point is inside the rigid body
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      virtual bool isIn( double const& x, double const& y, double const& z )
         const = 0;

      /** @brief Returns the level set value of a point from the rigid body
      @param pt the point */
      virtual double level_set_value( geomVector const& pt ) const = 0;

      /** @brief Returns the level set value of a point from the rigid body
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      virtual double level_set_value( double const& x
                                    , double const& y
                                    , double const& z ) const = 0;

      /** @brief Compute the surface points by discretizing the rigid body
      surface in approximately equal areas (if possible) */
      virtual void compute_surface_points( ) = 0;


      /** @brief Compute number of points on a rigid body
      @param surface_cell_scale scale of surface cell compared with the grid
      @param dx grid size */
      virtual void compute_number_of_surface_variables(
               double const& surface_cell_scale, double const& dx ) = 0;

      /** @brief Write the surface discretization points of all RB in
      respective CSV files */
      void write_surface_discretization( const std::string& file );

      /** @brief Correct surface discretization due to PBC */
      void correct_surface_discretization( FV_Mesh const* MESH );

      /** @brief Translate the surface points already discretized at
      the start of simulation */
      void translate_surface_points( geomVector const& delta);

      /** @brief Initializa the surface variables for a rigid body
      such as points, normals, and area */
      void initialize_surface_variables( );

      /** @brief Update the pressure force on a surface point
      @param i index of the surface point
      @param value the value to assign for pressure force */
      void update_Pforce_on_surface_point( size_t const& i
                                         , geomVector const& value );

      /** @brief Update the viscous force on a surface point
      @param i index of the surface point
      @param value the value to assign for viscous force */
      void update_Vforce_on_surface_point( size_t const& i
                                         , geomVector const& value );

      /** @brief Update the temperature gradient on a surface point
      @param i index of the surface point
      @param value the value to assign for temperature gradient */
      void update_Tgrad_on_surface_point( size_t const& i
                                        , double const& value );

      /** @brief Update the RB position and velocity
      @param pos updated position
      @param vel updated translation velocity */
      virtual void update_RB_position_and_velocity(geomVector const& pos,
                                                   geomVector const& vel,
                                                   geomVector const& ang_vel,
                         vector<geomVector> const& periodic_directions,
                         double const& time_step) = 0;

      /** @brief Update additional parameters of each RB type */
      virtual void update_additional_parameters( ) = 0;

      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      FS_RigidBody* m_geometric_rigid_body; /**< Pointer to the corresponding
    	geometric rigid body */
      vector<geomVector*> m_surface_points; /**< vector of points distributed on
      	the surface of the particle to compute surface integrals */
      vector<geomVector*> m_surface_area; /**< vector of the area associated
         with the points distributed on the surface of the particle */
      vector<geomVector*> m_surface_normal; /**< vector of the normal associated
         with the points distributed on the surface of the particle */
      vector<geomVector> m_surface_Pforce; /**< vector of the pressure force
         on the points distributed on the surface of the particle */
      vector<geomVector> m_surface_Vforce; /**< vector of the viscous force
         on the points distributed on the surface of the particle */
      vector<double> m_surface_Tgrad; /**< vector of the temperature gradient
         on the points distributed on the surface of the particle */
      vector<geomVector*> m_halo_zone; /**< vector of min and max extents
         of rigid body halozone, required for void fraction detection */
      size_t Ntot; /** < Stores the total number of surface points on a rigid
         body */
      //@}


   //-- Methods

      /**@name Methods */
      //@{

      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied DS_RigidBody object */
      DS_RigidBody( DS_RigidBody const& copy );
      //@}
};

#endif
