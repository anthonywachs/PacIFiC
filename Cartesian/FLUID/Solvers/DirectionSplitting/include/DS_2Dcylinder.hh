#ifndef _DS_2DCYLINDER__
#define _DS_2DCYLINDER__

#include <DS_RigidBody.hh>
#include <string>
using std::string;


/** @brief The class DS_2Dcylinder.

A moving or stationary rigid 2Dcylinder in the Direction Splitting solver.

@author A. Wachs - Pacific project 2021 */

class DS_2Dcylinder: public DS_RigidBody
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_2Dcylinder();

      /** @brief Constructor with arguments
      @param pgrb pointer to the corresponding geometric rigid body */
      DS_2Dcylinder( FS_RigidBody* pgrb );

      /** @brief Destructor */
      ~DS_2Dcylinder();
      //@}


   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief Updates the 2Dcylinder features from its corresponding
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

      /** @brief Computes the min and max extents of the 2Dcylinder halozone
      , required for the computation of void fraction */
      void compute_rigid_body_halozone( );

      /** @brief Compute the surface points by discretizing the 2Dcylinder
      surface in approximately equal areas (if possible) */
      void compute_surface_points( );

      /** @brief Compute number of points on a 2Dcylinder
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
      @param copy copied DS_2Dcylinder object */
      DS_2Dcylinder( DS_2Dcylinder const& copy );
      //@}
};

#endif
