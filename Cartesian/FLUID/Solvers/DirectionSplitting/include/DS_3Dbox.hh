#ifndef _DS_3DBOX__
#define _DS_3DBOX__

#include <DS_RigidBody.hh>
#include <string>
using std::string;


/** @brief The class DS_3Dbox.

A moving or stationary rigid 3D box of axisymmetric cross-section in the
Direction Splitting solver.

@author A. Wachs - Pacific project 2021 */

class DS_3Dbox: public DS_RigidBody
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_3Dbox();

      /** @brief Constructor with arguments
      @param pgrb pointer to the corresponding geometric rigid body */
      DS_3Dbox( FS_RigidBody* pgrb );

      /** @brief Destructor */
      ~DS_3Dbox();
      //@}


   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief Updates the sphere features from its corresponding
      geometric rigid body */
      void update();
      //@}


   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Writes the attributes in a stream
      @param out output stream
      @param indent_width indentation width */
      void display( ostream& out, size_t const& indent_width ) const;

      /** @brief Computes the min and max extents of the 3D cylinder halozone
      , required for the computation of void fraction */
      void compute_rigid_body_halozone( double const& dx );

      /** @brief Compute the surface points by discretizing the 3D cylinder
      surface in approximately equal areas (if possible) */
      void compute_surface_points( );

      /** @brief Returns whether a point is inside the rigid body
      @param pt the point */
      bool isIn( geomVector const& pt ) const;

      /** @brief Returns whether a point is inside the rigid body
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      bool isIn( double const& x, double const& y, double const& z ) const;

      /** @brief Returns the level set value of a point from the rigid body
      @param pt the point */
      double level_set_value( geomVector const& pt ) const;

      /** @brief Returns the level set value of a point from the rigid body
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      double level_set_value( double const& x
                            , double const& y
                            , double const& z ) const;

      /** @brief Returns the distance of a point with the rigid body
      with a given ray vector and source
      @param pt the point
      @param direction x, y or z (=0, 1 or 2)
      @param positive true if search in the positive direction of the coordinate
      axis and false otherwise */
      double get_distanceTo( geomVector const& source,
                             geomVector const& rayDir,
                             double const& delta ) const;

      /** @brief Returns rigid body velocity including rotation speed at pt
      @param pt the point */
      geomVector get_rigid_body_velocity( geomVector const& pt ) const;

      /** @brief Returns rigid body angular velocity */
      geomVector get_rigid_body_angular_velocity( ) const;

      /** @brief Returns pointer to the rigid body gravity center */
      geomVector const* get_ptr_to_gravity_centre( ) const;

      /** @brief Returns circumscribed radius */
      double get_circumscribed_radius( ) const;

      /** @brief Returns a tuple of mass and density of RB */
      std::tuple<double,double,double> get_mass_and_density_and_moi() const;

      /** @brief Update the RB position and velocity
      @param pos updated position
      @param vel updated translation velocity */
      void update_RB_position_and_velocity(geomVector const& pos,
                                                   geomVector const& vel,
                                                   geomVector const& ang_vel,
                         vector<geomVector> const& periodic_directions,
                         double const& time_step);

      /** @brief Update additional parameters of each RB type */
      void update_additional_parameters( );

      /** @brief Compute number of points on a 3D cylinder
      @param surface_cell_scale scale of surface cell compared with the grid
      @param dx grid size */
      void compute_number_of_surface_variables( double const& surface_cell_scale
                                              , double const& dx );

      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      vector<double> box_min, box_max;
      vector<size_t> Npoints;    /** Surface points on each direction
                                     of 3D box.*/
      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied DS_3Dbox object */
      DS_3Dbox( DS_3Dbox const& copy );
      //@}
};

#endif
