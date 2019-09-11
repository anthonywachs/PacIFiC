#ifndef MAC_APPLICATION_HH
#define MAC_APPLICATION_HH

#include <MAC_Object.hh>
#include <MAC.hh>

#include <string>
#include <stringVector.hh>

class MAC_ListIdentity ;
class MAC_ListIterator ;
class MAC_ModuleExplorer ;
class MAC_ObjectRegister ;

/*
Applications, performing their specific tasks.

The program execution consists of five stages :
  1. Initial stage (Big-Bang time) : all statics are initialized and
     the only instance of MAC_Root is created.
  2. The data deck is read and stored in memory.
  3. An instance of a concrete subclass of MAC_Application is created.
  4. Program core execution : the program execution proceeds by performing
     its specific tasks.
  5. Final stage : termination of the only instance of MAC_Root, leading 
     to the termination of all objects belonging to a ownership tree whose
     root node is not the NULL object.
The MAC_Application class provides an interface for executing the 
specific tasks of the above fouth point.

FRAMEWORK INSTANTIATION

   CASE 1 : derivation of a concrete subclass

   1. Derive a concrete subclass, say MyAppli.
   2. Choose a name for MyAppli, say "my_appli".
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the prototype to be registered :
      5.1 Implement a default constructor that initializes the
          `MAC_Application::' subobject by calling
               `MAC_Application( std::string const& )'
          with "my_appli" as argument.
          Example of pseudo-code :
          | MyAppli:: MyAppli( void ) : MAC_Application( "my_appli" ) {}
      5.2 Define and initialize a static instance by calling the default
          constructor.
             declaration (in the header file, eg MyAppli.hh) :
             | static MyAppli const* PROTOTYPE ;
             definition (in the implementation file, eg MyAppli.cc) :
             | MyAppli const* MyAppli::PROTOTYPE = new MyAppli() ;'
   6. Implement a private constructor that initializes the 
      `MAC_Application::' subobject by calling
                 `MAC_Application( MAC_Object*, MAC_ModuleExplorer const* )'
      or
                 `MAC_Application( MAC_Object*, stringVector& )'
      Example of pseudo-code :
      | MyAppli:: MyAppli( MAC_Object* a_owner,
      |                    MAC_ModuleExplorer const* exp )
      |    : MAC_Application( a_owner, exp ), ...
      | { ... }
   7. Implement the `::create_replica' method that allocates an object
      of type `MyAppli' initialized using the private constructor described
      above, and subsequently return a pointer to that object.
      Example of pseudo-code :
      | MyAppli* MyAppli::create_replica( MAC_Object* a_owner,
      |                                   MAC_ModuleExplorer const* exp ) const
      | {
      |    MAC_LABEL( "MyAppli::create_replica" ) ;
      |    MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;
      |    MyAppli* result = new MyAppli( a_owner, exp ) ;
      |    MAC_CHECK( create_replica_POST( result, a_owner, exp ) ;
      |    return result ;
      | }
   8. Implement the `::run' method

   CASE 2 : derivation of an abstract subclass

   1. Derive an abstract subclass, say MyAppli.
   2. Implement a protected virtual destructor.
   3. Implement a protected constructor that initializes the
      `MAC_Application::' subobject by calling
               `MAC_Application( std::string const& )'
      Example of pseudo-code :
      | MyAppli:: MyAppli( std::string const& name ) 
      |    : MAC_Application( name ) {}
      This constructor is devoted to be used by the concrete subclasses 
      of MyAppli for the registration of their prototype.
   4. Implement a protected constructor that initializes the 
      `MAC_Application::' subobject by calling
                 `MAC_Application( MAC_Object*, MAC_ModuleExplorer const* )'.
      Example of pseudo-code :
      | MyAppli:: MyAppli( MAC_Object* a_owner,
      |                    MAC_ModuleExplorer const* exp )
      |    : MAC_Application( a_owner, exp ), ...
      | { ... }
      This constructor is devoted to be used to initialize the MyAppli
      base class subobject when creating objects of concrete subclasses
      of MyAppli (such creations are performed in the `create_replica::'
      method whose implementation is deferred into those concrete subclasses).
*/

class MAC_Application : public MAC_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance of `MAC_Application::' according
      // to the data attainable by `exp'.
      static MAC_Application* make( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time = 0. ) ;

//       Create application from args NOT SUPPORTED for now
//       // Create and return an instance according to the command-line
//       // arguments gathered in `args'. 
//       static MAC_Application* make( MAC_Object* a_owner,
//                                     stringVector& args ) ;

   //-- Program core execution(0.2)

      // Perform the specific tasks of the application (Called by main()).
      virtual void run( std::string const& inputRestartFileName = 
      	MAC::undefined_string ) = 0 ;

   //-- Persistence

      void register_storable_objects( void ) ;

      void write_storable_objects( double const& time = 0 ) const ;

      void restore_registered_objects( MAC_ObjectReader* ret ) const ;
      
   //-- Reload & follow
   
      static bool is_follow( void ); 
      
      static bool is_reload( void );         

   protected: //--------------------------------------------------------

   //-- Plug in

      virtual ~MAC_Application( void ) ;

      // Registration of an instance as `name'.
      MAC_Application( std::string const& name ) ;

      // In the constructor called by `::create_replica' or 
      // `::create_replica_from_args' : initialization the base class subobject
      // (`exp' can possibly be 0 ).
      MAC_Application( MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp,
		bool const& is_follow_ = false,
		bool const& is_reload_ = false ) ;

      virtual MAC_Application* create_replica( 
		MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp,
		double const& initial_time ) const = 0 ;

      // IMPLEMENTATION : raise a fatal error.
      // Create application from args NOT SUPPORTED for now
//       virtual MAC_Application* create_replica_from_args( 
//                                    MAC_Object* a_owner,
//                                    stringVector& args ) const ;

      bool is_a_prototype( void ) const ;

   //-- Command line(1010.0)

      void notify_error_in_arguments( void ) const ;

      virtual void print_usage( void ) const ;
      virtual void print_options( void ) const ;
      virtual void print_operands( void ) const ;
      virtual void print_exit_status( void ) const ;

      std::string usage_title( std::string const& name ) const ;
      std::string options_title( void ) const ;
      std::string operands_title( void ) const ;
      std::string exit_status_title( void ) const ;

   //-- Preconditions, Postconditions, Invariant    

      bool create_replica_PRE( MAC_Object* a_owner,
                               MAC_ModuleExplorer const* exp ) const ;
      bool create_replica_POST( MAC_Application const* result,
				MAC_Object* a_owner,
				MAC_ModuleExplorer const* exp ) const ;

//       Create application from args NOT SUPPORTED for now
//       bool create_replica_from_args_POST( MAC_Application const* result,
// 		         		  MAC_Object* a_owner,
// 				          stringVector& args ) const ;

      virtual bool invariant( void ) const ;

   //-- Persistence
      
      // Extend `list' (with the `MAC_ListIdentity::extend' method) so that it
      // contains all objects required by the storage and retrieval mechanisms.
      // IMPLEMENTATION : do nothing, i.e. leave `list' unchanged.
      virtual void add_storable_objects( MAC_ListIdentity* list ) ;

      // name of the module containing the data related to
      // the storage mechanism (of persistence)
      // IMPLEMENTATION : "MAC_ObjectWriter"
      virtual std::string const& object_writer_module_name( void ) const ;
      
      // Restart cycle number
      size_t restart_cycle_number( void ) const;
      
      // Current output restart file name
      std::string output_restart_file_name( void ) const;
      
      // Force writing restart files at last time
      bool force_write_restart_files_at_last_time( void ) const;
      
      // Input file name for restart
      virtual std::string input_restart_file_name( void ) const;
      
      // Swap restart file names in case of last_two_cycles and input 
      // restart file name is similar to first output restart file name
      void swap_restart_file_names( std::string const& inputfilename );       

   private: //----------------------------------------------------------

      MAC_Application( void ) ;
      MAC_Application( MAC_Application const& other ) ;
      MAC_Application& operator=( MAC_Application const& other ) ;

      void print_usage_then_exit( int exit_status = 0 ) const ;

      void initialize_objects_storage( MAC_ModuleExplorer const* exp ) ;

      static MAC_ObjectRegister* plugins_map( void ) ;

   //-- Attributes
      
      bool const IS_PROTO ;

      MAC_ObjectWriter* SAVER ;

      // List of the persistent objects :
      MAC_ListIdentity* persistent_objects ;
      mutable MAC_ListIterator* persistent_objects_it ;
      
      static bool B_FOLLOW;
      static bool B_RELOAD;
} ;

#endif



