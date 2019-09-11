#ifndef MAC_EXTERNAL_API_HH
#define MAC_EXTERNAL_API_HH

#include <MAC_Object.hh>

class MAC_ObjectRegister ;
class stringVector ;

/*
External applications, performing their specific initialization
and termination.

FRAMEWORK INSTANTIATION
   1. Derive a concrete subclass, say MyExtAPI.
   2. Choose a name for MyAppli, say "my_api", 
      and a priority level, say pl.
   3. Declare all constructors private.
   4. Define an instance to be registered :
      4.1 Implement a default constructor that initializes the
          `MAC_ExternalAPI::' subobject by calling
               `MAC_ExternalAPI( std::string const&, size_t )'
          with "my_api" and pl as argument.
          Example of pseudo-code :
          | MyExtAPI:: MyExtAPI( void ) : MAC_ExternalAPI( "my_api", pl )" ) {}
      4.2 Define and initialize a static instance by calling the default
          constructor.
             declaration (in the header file, eg MyExtAPI.hh) :
             | static MyExtAPI const* SINGLETON ;
             definition (in the implementation file, eg MyExtAPI.cc) :
             | MyExtAPI const* MyExtAPI::SINGLETON = new MyExtAPI() ;'
   5. Implement `::initialize' that does the initialization required
      by the external API at hand.
   6. Implement a private destructor that does the deinitialization required
      by the external API at hand. 

PUBLISHED*/

class MAC_ExternalAPI : public MAC_Object
{

   public: //-----------------------------------------------------------
      
   //-- Management of all registered instances

      // Call `::initialize' for all registered instances in the
      // reverse order of their priority level (instances with
      // higher priority are initialized first).
      static void initialize_all_APIs( int& argc, char**& argv ) ;

      // Call `MAC_Object::destroy' for all registered instances in the
      // order of their priority level (instances with
      // lower priority are terminanted first).
      static void terminate_all_APIs( void ) ;

   protected: //--------------------------------------------------------

   //-- Plug in

      virtual ~MAC_ExternalAPI( void ) ;

      // registration of `self', calling it `a_name' and
      // setting its priority level to `a_priority_level'.
      MAC_ExternalAPI( std::string const& a_name, size_t a_priority_level ) ;

   //-- Current instance management

      // Initialize `self'.
      virtual void initialize( int& argc, char**& argv ) = 0 ;

   private: //----------------------------------------------------------

      MAC_ExternalAPI( void ) ;
      MAC_ExternalAPI( MAC_ExternalAPI const& other ) ;
      MAC_ExternalAPI& operator=( MAC_ExternalAPI const& other ) ;
      
      static MAC_ObjectRegister* plugins_map( void ) ;
      static stringVector& plugins_names( void ) ;

   //-- Attributes

      std::string MY_NAME ;
      size_t MY_PRIORITY ;
} ;

#endif



