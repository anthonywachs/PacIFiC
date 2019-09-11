#ifndef MAC_APPLICATION_RESTORER_HH
#define MAC_APPLICATION_RESTORER_HH

#include <MAC_Application.hh>
#include <MAC.hh>

#include <string>

class MAC_ObjectReader ;

/*
Launchers of  executions of `MAC_Application::' instances, starting from one
of their internal states that was previously saved on disk.

Given a `MAC_Application::' instance whose internal state was saved
on disk using the facilities provided by `MAC_ObjectWriter::',
an object of `MAC_ApplicationRestorer'
   - first creates a `MAC_Application::' instance from the original
     data deck,
   - then modifies its internal state according to one the savings performed
     during the original execution.
Thus, the newly created `MAC_Application::' instance has a state identical
to that of the original `MAC_Application::' instance at its execution point
where disk saving occured (with `MAC_ObjectWriter::'). Calling
`::run' for this new `MAC_Application::' instance will subsequently produce
the same results as continuing the original execution from where it
was saved.
   
PUBLISHED
*/

class MAC_ApplicationRestorer : public MAC_Application
{
   public: //----------------------------------------------------------------

   //-- Program core execution

      // Call `MAC_Application::run' on behalf of the `MAC_Application::'
      // instance whose internal self was restored when `self' was
      // initialized.
      virtual void run( std::string const& inputRestartFileName = 
      	MAC::undefined_string ) ;
      
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

     ~MAC_ApplicationRestorer( void ) ;
      MAC_ApplicationRestorer( MAC_ApplicationRestorer const& other ) ;
      MAC_ApplicationRestorer& operator=( 
                               MAC_ApplicationRestorer const& other ) ;

      // First build the application to be restarted from its original
      // data deck, then restore its internal state.
      MAC_ApplicationRestorer( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp ) ;

   //-- Plug in

      MAC_ApplicationRestorer( void ) ;

      virtual MAC_ApplicationRestorer* create_replica( 
		MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp,
		double const& initial_time ) const ;

   //-- Data deck for restarting

      // Create an instance of the module which has been saved calling
      // `MAC_ObjectWriter::write_data_deck_module()', and possibly
      // modify it according to the data stored in the file of name
      // `appendum_file'
      MAC_Module* create_modified_data_deck_module( 
	MAC_ModuleExplorer const* exp ) ;
					   
      // Input file name for restart
      std::string input_restart_file_name( void ) const;

   //-- Class Attributes

      static MAC_ApplicationRestorer const* PROTOTYPE ;
      
   //-- Attributes      

      MAC_ObjectReader* READER ;
      MAC_Application* APPLI ;
} ;

#endif
