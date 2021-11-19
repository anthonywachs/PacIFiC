#ifndef _FS_SPHERE__
#define _FS_SPHERE__

#include <FS_RigidBody.hh>
#include <iostream>
using std::istream ;


/** @brief Additional geometric parameters for the sphere */
struct FS_Sphere_Additional_Param
{
  double radius; /**< radius of the sphere */
};


/** @brief The class FS_Sphere.

A moving or stationary rigid sphere.

@author A. Wachs - Pacific project 2021 */

class FS_Sphere: public FS_RigidBody
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      FS_Sphere();

      /** @brief Constructor with arguments
      @param in input stream where features of rigid bodies are read
      @param id_ identification number */
      FS_Sphere( istream& in, size_t& id_ );

      /** @brief Destructor */
      ~FS_Sphere();
      //@}


   //-- Get methods

      /**@name Get methods */
      //@{
      /** @brief Returns a constant pointer to the structure containing the
      additional geometric parameters for the sphere */
      struct FS_Sphere_Additional_Param const*
      	get_ptr_FS_Sphere_Additional_Param() const;

      /** @brief Returns a constant pointer to the gravity center of the sphere */
      geomVector const*
         get_ptr_FS_Sphere_gravity_centre() const;

      //@}


   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief Updates the rigid body features
      @param in input stream where features of the rigid body are read */
      void update( istream& in );
      //@}


   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Writes the attributes in a stream
      @param out output stream
      @param indent_width indentation width */
      void display( ostream& out, size_t const& indent_width ) const;

      /** @brief Returns whether a point is inside the sphere
      @param pt the point */
      bool isIn( geomVector const& pt ) const;

      /** @brief Returns whether a point is inside the rigid body
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      bool isIn( double const& x, double const& y, double const& z ) const;

      /** @brief Returns the level set value of a point from a sphere
      @param pt the point */
      double level_set_value( geomVector const& pt ) const;

      /** @brief Returns the level set value of a point from a sphere
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      double level_set_value( double const& x
                            , double const& y
                            , double const& z ) const;

      /** @brief Returns whether a line originating from a point intersects the
      sphere, and if it does the distance from the point to the sphere
      surface
      @param pt the point
      @param direction x, y or z (=0, 1 or 2)
      @param positive true if search in the positive direction of the coordinate
      axis and false otherwise */
      tuple<bool,double,size_t> distanceTo( geomVector const& pt,
      	size_t const& direction,
      	bool const& positive ) const;
      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      struct FS_Sphere_Additional_Param m_agp_sphere; /**< Additional geometric
      	parameters for the sphere */
      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied FS_Sphere object */
      FS_Sphere( FS_Sphere const& copy );
      //@}


   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Sets the rigid body features from an input stream
      @param in input stream where features of the rigid body are read */
      void set( istream& in );
      //@}

};

#endif
