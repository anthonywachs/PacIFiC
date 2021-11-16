#ifndef _FS_RIGIDBODY__
#define _FS_RIGIDBODY__

#include <geomVector.hh>
#include <MAC_assertions.hh>
#include <string>
#include <vector>
#include <iostream>
#include <map>
using std::istream;
using std::ostream;
using std::string;
using std::vector;
using std::tuple;
typedef double Matrix3D[3][3];


/** @brief Geometric shape enumeration */
enum GEOMETRICSHAPE
{
  GEOM_SPHERE,
  GEOM_2DCYLINDER,
  GEOM_3DCYLINDER,
  GEOM_GENERAL_POLYGON,
  GEOM_GENERAL_POLYHEDRON,
  GEOM_SQUARE,
  GEOM_CUBE,
  GEOM_EQUILATERAL_TRIANGLE,
  GEOM_REGULAR_TETRAHEDRON
};


/** @brief The class FS_RigidBody.

A moving or stationary rigid body.

@author A. Wachs - Pacific project 2021 */

class FS_RigidBody
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      FS_RigidBody();

      /** @brief Destructor */
      virtual ~FS_RigidBody();
      //@}


   //-- Get methods

      /**@name Get methods */
      //@{
      /** @brief Returns the shape type */
      GEOMETRICSHAPE get_shape_type() const;

      /** @brief Returns the rigid body type */
      string get_type() const;
      //@}


   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief Updates the rigid body features
      @param in input stream where features of the rigid body are read */
      virtual void update( istream& in ) = 0;

      /** @brief Sets the value of the hydrodynamic force
      @param hf hydrodynamic force */
      void set_hydro_force( geomVector const& hf );

      /** @brief Sets the value of the hydrodynamic torque
      @param ht hydrodynamic torque */
      void set_hydro_torque( geomVector const& ht );

      /** @brief Sets the translational and angular velocity to zero */
      void nullify_velocity();

      /** @brief Changes rigid body type from particle to obstacle */
      void change_from_particle_to_obstacle();
      //@}


   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Writes the attributes in a stream
      @param out output stream
      @param indent_width indentation width */
      virtual void display( ostream& out, size_t const& indent_width )
      	const = 0;

      /** @brief Returns whether a point is inside the rigid body
      @param pt the point */
      virtual bool isIn( geomVector const& pt ) const = 0;

      /** @brief Returns whether a point is inside the rigid body
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      virtual bool isIn( double const& x, double const& y, double const& z )
      	const = 0;

      /** @brief Returns level set value of a point from the rigid body
      @param pt the point */
      virtual double level_set_value( geomVector const& pt ) const = 0;

      /** @brief Returns level set value of a point from the rigid body
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      virtual double level_set_value( double const& x
                                    , double const& y
                                    , double const& z ) const = 0;

      /** @brief Returns whether a line originating from a point intersects the
      rigid body, and if it does the distance from the point to the rigid body
      surface
      @param pt the point
      @param direction x, y or z (=0, 1 or 2)
      @param positive true if search in the positive direction of the coordinate
      axis and false otherwise */
      virtual tuple<bool,double,size_t> distanceTo( geomVector const& pt,
      	size_t const& direction,
      	bool const& positive ) const = 0;
      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      size_t m_space_dimension; /**< Space dimension */
      size_t m_Id; /**< Identification number */
      GEOMETRICSHAPE m_shape_type; /**< Shape type */
      string m_type; /**< type: P=particle, PP=periodic particle, O=obstacle or
      	PO=periodic obstacle */
      geomVector m_gravity_center; /**< Gravity center */
      double m_density; /**< Density */
      double m_volume; /**< Volume */
      double m_mass; /**< Mass */
      double m_circumscribed_radius; /**< Circumscribed radius */
      Matrix3D m_inertia; /**< Moment of inertia tensor */
      Matrix3D m_rotation_matrix; /**< Rotation matrix */
      geomVector m_translational_velocity; /**< Translational velocity */
      geomVector m_angular_velocity; /**< Angular velocity */
      geomVector m_hydro_force; /**< Hydrodynamic force */
      geomVector m_hydro_torque; /**< Hydrodynamic torque */
      double m_temperature; /**< Temperature */
      geomVector m_heat_flux; /**< Heat flux */
      vector<geomVector>* m_periodic_directions; /**< periodic clones whose
      	position is gravity_center + (*periodic_directions)[i] */
      static vector<string> GEOMETRICSHAPE_name;
      //@}


   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Writes the general attributes in a stream
      @param out output stream
      @param indent_width indentation width */
      void display_general( ostream& out, size_t const& indent_width ) const;
      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied FS_RigidBody object */
      FS_RigidBody( FS_RigidBody const& copy );
      //@}
};

#endif
