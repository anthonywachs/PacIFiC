#ifndef _DS_RIGIDBODY__
#define _DS_RIGIDBODY__

#include <geomVector.hh>
#include <MAC_assertions.hh>
#include <vector>
#include <iostream>
#include <map>
using std::ostream;
using std::vector;
using std::tuple;
class FS_RigidBody;
class FV_DiscreteField;


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


      /** @brief Returns whether a line originating from a point intersects the
      rigid body, and if it does the distance from the point to the rigid body
      surface
      @param pt the point
      @param direction x, y or z (=0, 1 or 2)
      @param positive true if search in the positive direction of the coordiante
      axis and false otherwise */
      tuple<bool,double,size_t> distanceTo( geomVector const& pt,
      	size_t const& direction,
      	bool const& positive ) const;

      /** @brief Computes the hydrodynamic force and torque and stores the
      values in the corresponding geometric rigid body
      @param PP the pressure field
      @param UU the velocity field */
      virtual void compute_hydro_force_torque( FV_DiscreteField const* PP,
	FV_DiscreteField const* UU ) = 0;
      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      FS_RigidBody* m_geometric_rigid_body; /**< Pointer to the corresponding
    	geometric rigid body */
      vector<geomVector> m_surface_points; /**< vector of points distributed on
      	the surface of the particle to compute surface integrals */
      //@}


   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Performs numerical integration on the set of points distributed
      on the surface to yield the hydrodynamic force and torque and stores the
      values in the corresponding geometric rigid body
      @param PP the pressure field
      @param UU the velocity field */
      void compute_surface_integrals_hydro_force_torque(
      	FV_DiscreteField const* PP,
	FV_DiscreteField const* UU );
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
