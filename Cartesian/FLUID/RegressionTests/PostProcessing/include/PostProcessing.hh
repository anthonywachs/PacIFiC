#ifndef POST_PROCESSING_HH
#define POST_PROCESSING_HH

#include <FV_DiscreteField.hh>
#include <DS_AllRigidBodies.hh>
#include <MAC_Communicator.hh>
#include <FV_DomainAndFields.hh>

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
  geomVector center;
  geomVector length;
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
      void prepare_fieldVolumeAverageInBox( MAC_ModuleExplorer const* exp
                                          , FV_DomainAndFields const* dom);

      /** @brief Compute average field in the list of boxes */
      void compute_fieldVolumeAverageInBox( );

      /** @brief Compute the L2 norm of a given field
      @param FF field
      @param level field level */
      double compute_mean( FV_DiscreteField const* FF
                         , size_t const& level
                         , size_t const& comp);

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


      //@}

};

#endif
