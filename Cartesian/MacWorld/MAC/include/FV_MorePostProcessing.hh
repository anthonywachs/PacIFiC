#ifndef MAC_APPLICATION_RESTORER_HH
#define MAC_APPLICATION_RESTORER_HH

#include <MAC_Application.hh>
#include <MAC.hh>

#include <string>

class MAC_ObjectReader ;
class FV_StepByStepProgression ;


class FV_MorePostProcessing : public MAC_Application
{
   public: //----------------------------------------------------------------

   //-- Program core execution

      /**
        @brief Run simulation
        @warning Call MAC_Application::run on behalf of the MAC_Application::
            instance whose internal self was restored when `self' was
            initialized.
        @par Content
          Calls FV_StepByStepProgression::do_more_post_processing
      */
      virtual void run( std::string const& inputRestartFileName = 
      	MAC::undefined_string ) ;
      
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

     ~FV_MorePostProcessing( void ) ;
      FV_MorePostProcessing( FV_MorePostProcessing const& other ) ;
      FV_MorePostProcessing& operator=( 
                               FV_MorePostProcessing const& other ) ;

      FV_MorePostProcessing( MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp ) ;

   //-- Plug in

      FV_MorePostProcessing( void ) ;

      virtual FV_MorePostProcessing* create_replica( 
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time ) const ;

   //-- Data deck for restarting

      MAC_Module* create_modified_data_deck_module( void ) ;

   //-- Class Attributes
      
      static FV_MorePostProcessing const* PROTOTYPE ;
      
   //-- Attributes      

      MAC_ObjectReader* READER ;
      FV_StepByStepProgression* APPLI ;
      MAC_ModuleExplorer const* ROOT_EXP ;
} ;

#endif
