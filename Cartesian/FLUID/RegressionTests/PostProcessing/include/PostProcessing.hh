#ifndef POST_PROCESSING_HH
#define POST_PROCESSING_HH

#include <FV_DiscreteField.hh>
#include <DS_AllRigidBodies.hh>
#include <DS_RigidBody.hh>
#include <MAC_Communicator.hh>
#include <FV_DomainAndFields.hh>
#include <MAC_DoubleVector.hh>

/** @brief The class PostProcessing

The static function class for data post processing


@author A. Goyal - Pacific project 2022 */

/** @brief The Structure fieldVolumeAverageInBox.
Information to compute the mean value of a field into a defined box */
struct fieldVolumeAverageInBox
{
  FV_DiscreteField const* FF;
  string field_name;
  string box_name;
  MAC_DoubleVector* center;
  MAC_DoubleVector* length;
  bool withPorosity;
};

/** @brief The Structure fieldVolumeAverageAroundRB.
Information to compute the mean value of a field into a defined box */
struct fieldVolumeAverageAroundRB
{
  FV_DiscreteField const* FF;
  string field_name;
  double volumeWidth;
  bool withPorosity;
  size_t kernelType;
};


class PostProcessing
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      PostProcessing();

      /** @brief Constructor with arguments
      @param is_solids RB present or not
      @param allrigidbodies pointer to all RB */
      PostProcessing( bool is_solids_
                    , DS_AllRigidBodies const* allrigidbodies_
                    , size_t const& dim_
                    , MAC_Communicator const* macCOMM_);

      /** @brief Destructor */
      ~PostProcessing( );

      /** @brief Read parameters from prob def for box averaging */
      void prepare_fieldVolumeAverageInBox( MAC_Object* a_owner
                                          , MAC_ModuleExplorer const* exp
                                          , FV_DomainAndFields const* dom);

      /** @brief Read parameters from prob def for averaging around RB*/
      void prepare_fieldVolumeAverageAroundRB( MAC_Object* a_owner
                                          , MAC_ModuleExplorer const* exp
                                          , FV_DomainAndFields const* dom);

      /** @brief Compute average field in the list of boxes
      @param allrigidbodies pointer to all rigid bodies */
      void compute_fieldVolumeAverageInBox(
                                    DS_AllRigidBodies* allrigidbodies);

      /** @brief Compute average field around a RB for all list of boxes
      @param allrigidbodies pointer to all rigid bodies */
      void compute_fieldVolumeAverageAroundRB(
                                    DS_AllRigidBodies* allrigidbodies);

   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{

      //@}

   private: //----------------------------------------------------------------

   //-- Attributes
      /** @name Constructors & Destructor */
      //@{


      /** @brief Copy constructor */
      PostProcessing( PostProcessing const& other ) {}

      //@}

      /**@brief Perform kernel weighting function for averaging calculus
         @param distance between center of the particle and cell considered
         @param radius/characteristic length scale of considered particle
         @param boxWidth of the control volume
         @param kernel_type defined the kernel to be used
       */
       static double kernel( double distance
                           , double radius
                           , double boxWidth
                           , size_t kernel_type );

      /**@brief Return the grid index in a direction of box box_extents
         @param bounds bounds of the box in particular direction
         @param FF field
         @param dir direction
         @param comp field component
       */
       size_t_vector get_local_index_of_extents( class doubleVector& bounds
                                              , FV_DiscreteField const* FF
                                              , size_t const& dir
                                              , size_t const& comp);

      /**@brief Transform a point using PBC
         @param x point
         @param dir transformation direction
       */
       double periodic_transformation( double const& x, size_t const& dir);

   private: //----------------------------------------------------------------

      // Pointers to the constant fields and primary grid
      FV_DiscreteField const* UF ;
      FV_DiscreteField const* PF ;
      FV_DiscreteField const* TF ;
      FV_Mesh const* MESH ;

      bool m_is_solids;
      DS_AllRigidBodies const* m_allrigidbodies;
      size_t m_dim;
      MAC_Communicator const* m_macCOMM; /**< Variable for communication
      between processors */

      list<struct fieldVolumeAverageInBox> m_fieldVolumeAverageInBox_list;

      list<struct fieldVolumeAverageAroundRB> m_fieldVolumeAverageAroundRB_list;


      //@}

};

#endif
