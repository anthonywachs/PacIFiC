#ifndef _FS_3DBOX__
#define _FS_3DBOX__

#include <FS_RigidBody.hh>
#include <iostream>
using std::istream ;


/** @brief Additional geometric parameters for the sphere */
struct FS_3Dbox_Additional_Param
{
  vector<geomVector> corners; /**< Corner coordinates of the polyhedron */
  vector<vector<size_t> > facesVec; /**< polygonal faces numbering */
  vector<geomVector> ref_corners; /**< Reference Corner coordinates
                                       of the polyhedron */
  // geomVector* g2; /**< slightly randomly translated gravity center */
  // geomVector coor_min; /**< minimal coordinates of the bounding box */
  // geomVector coor_max; /**< maximal coordinates of the bounding box */
};


/** @brief The class FS_3Dbox.

A moving or stationary rigid 3D cylinder of axisymmetric cross-section.

@author A. Wachs - Pacific project 2021 */

class FS_3Dbox: public FS_RigidBody
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      FS_3Dbox();

      /** @brief Constructor with arguments
      @param in input stream where features of rigid bodies are read
      @param id_ identification number */
      FS_3Dbox( istream& in, size_t& id_ );

      /** @brief Destructor */
      ~FS_3Dbox();
      //@}


   //-- Get methods

      /**@name Get methods */
      //@{
      /** @brief Returns a constant pointer to the structure containing the
      additional geometric parameters for the sphere */
      struct FS_3Dbox_Additional_Param const*
      	get_ptr_FS_3Dbox_Additional_Param() const;

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

      /** @brief Returns the level set value of a point from a cylinder
      @param pt the point */
      double level_set_value( geomVector const& pt ) const;

      /** @brief Returns the level set value of a point from a cylinder
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      double level_set_value( double const& x
                            , double const& y
                            , double const& z ) const;

      /** @brief Returns 3Dbox velocity including rotation speed at pt
      @param pt the point */
      geomVector rigid_body_velocity( geomVector const& pt ) const;

      /** @brief Calculate determinant 4 X 4 for checking
      a point in tetrahedron */
      double calcPointDeterm4by4( const geomVector &pointOne,
              const geomVector &pointTwo, const geomVector &pointThree,
              const geomVector &pointFour ) const;

      /** @brief Check whether a point is inside a tetrahedron */
      bool checkPointInTetrahedron( const geomVector &pointOne,
           const geomVector &pointTwo, const geomVector &pointThree,
           const geomVector &pointFour, const geomVector &pointToCheck ) const;

      /** @brief Calculates the relative
      distance of a point from a tetrahedron */
      double DistOfPointFromTetrahedron( const geomVector &pointOne,
         const geomVector &pointTwo, const geomVector &pointThree,
         const geomVector &pointFour, const geomVector &pointToCheck ) const;

      /** @brief Reverse tranform the corners of the actual 3D box to the
      reference 3D box
      @param pt point to transform */
      void compute_reverseTransformationOfCorners( );
      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      struct FS_3Dbox_Additional_Param m_agp_3dbox; /**< Additional
      	geometric parameters for the 3D cylinder */
      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied FS_3Dbox object */
      FS_3Dbox( FS_3Dbox const& copy );
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
