#ifndef _DS_SPHERE__
#define _DS_SPHERE__

#include <DS_RigidBody.hh>
#include <string>
using std::string;


/** @brief The class DS_Sphere.

A moving or stationary rigid sphere in the Direction Splitting solver.

@author A. Wachs - Pacific project 2021 */

class DS_Sphere: public DS_RigidBody
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_Sphere();

      /** @brief Constructor with arguments
      @param pgrb pointer to the corresponding geometric rigid body */
      DS_Sphere( FS_RigidBody* pgrb );

      /** @brief Destructor */
      ~DS_Sphere();
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

      /** @brief Computes the min and max extents of the sphere halozone
      , required for the computation of void fraction */
      void compute_rigid_body_halozone( );

      /** @brief Compute the surface points by discretizing the sphere
      surface in approximately equal areas (if possible) */
      void compute_surface_points( );

      /** @brief Translate the surface points already discretized on sphere
      surface */
      void translate_surface_points( geomVector const& delta);

      /** @brief Compute number of points on a sphere
      @param surface_cell_scale scale of surface cell compared with the grid
      @param dx grid size */
      void compute_number_of_surface_variables( double const& surface_cell_scale
                                              , double const& dx );


      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied DS_Sphere object */
      DS_Sphere( DS_Sphere const& copy );
      //@}
};

#endif
