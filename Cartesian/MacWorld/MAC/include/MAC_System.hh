#ifndef MAC_SYSTEM_HH
#define MAC_SYSTEM_HH

#include <MAC_Object.hh>

#include <iosfwd>
class stringVector ;

/*
interfaces to operating system features that are available on
all operating systems
*/

class MAC_System : public MAC_Object
{
   public: //-----------------------------------------------------------------
      
      // host name
      static std::string const& host_name( void ) ;
      
      // process id
      static int process_id( void ) ;

      // character used to form path name description on native system
      static char path_name_separator( void ) ;

      // Does `str' matches `pattern' description ?
      static bool matches(  std::string const& str,
                            std::string const & pattern ) ;
      
      // search for `filename' entry in `directory'
      static void find( std::string const& directory,
                        std::string const & filename,
                        stringVector & result,
                        bool recurse=true,
                        std::string const & prefix="" ) ;

      // execute command `cmd' in `directory' and save output in `output_file'
      static int run( stringVector const& cmd, 
                      std::string const& directory,
                      std::string const& output_file ) ;
      
      // build directory
      static bool mkdir( std::string const& directory ) ;

      static bool changedir( std::string const& directory ) ;
      
      static bool copy( std::string const& src, std::string const& dest ) ;

      // remove the file called `filename'
      static bool erase( std::string const& filename ) ;

      // working directory
      static std::string const& working_directory( void ) ;

      // kernel name
      static std::string const& sysname( void ) ;

      // currently amount of used TIME
      static double user_time( void ) ;

      // absolute TIME elapsed from 1/1/1970 in s
      static double epoch_time( void ) ;

      // delay execution of the calling process for (at least) 
      // `msec' milliseconds
      static void sleep( size_t msec ) ;

      // currently amount of used memory
      static size_t used_memory( void ) ;
      
      // Name of compiler used to build MAC library.
      static std::string const& compiler_name( void ) ;

      // Stripped non-directory suffix from file name `path_and_name'.
      // Directory separator is given by `separator'.
      static std::string dirname( std::string const& path_and_name,
                                  char separator = path_name_separator() ) ;
      
      // Return `path_and_name' with any leading directory components removed.
      // Directory separator is given by `separator'.
      static std::string basename( std::string const& path_and_name,
                                   char separator = path_name_separator() ) ;

      // Return absolute path corresponding to relative or absolute `filename'.
      static std::string absolute_path( std::string const& filename ) ;
      
      // Active exception treatment.
      static void exception_trapping( void ) ;

      // Is big endian encoding format supported on the current machine ?
      static bool big_endian_encoding( void ) ;

      // new recommended block size for a block of memory that needs
      // to be made larger (`current_size' is the actual size of the block
      // and `expected_size' is the new needed size)
      static size_t new_block_size( size_t current_size, 
                                    size_t expected_size ) ;

      // file access on reading
      static bool can_read( std::string const& filename ) ;
      
      // file access on writing
      static bool can_write( std::string const& filename ) ;

      // close file descriptor
      static int close_file_descriptor(int fildes) ;

      // open and possibly create a file or device
      static int open_file_descriptor(std::string filename, int flags ) ;

      // process management
      class Process 
      {
         public :
            
            Process( stringVector const& args ) ;
            
            size_t id( void ) const ;
            
            static void wait_for_child_processes( size_t nb, Process ** child ) ;
            
            ~Process(void) ;
            
            void set_name( std::string a_name ) ;
            
            std::string const& name( void ) const ;
            
         private : 
            
            bool running ;
            std::string my_name ;
            size_t pid ;
            void * sdef ;
      } ;
      
      // exit current process
      static void exit( int exit_status ) ;

      // environment variable content (empty string if not defined)
      static std::string getenv( std::string const& str ) ;

   protected: //--------------------------------------------------------------
      
   private: //----------------------------------------------------------------

      MAC_System( void ) ;
     ~MAC_System( void ) ;
      MAC_System( MAC_System const& other ) ;
      MAC_System& operator=( MAC_System const& other ) ;

      static std::string const& string( std::string const& cmd ) ;
      
      static MAC_System const* SINGLETON ;
      
} ;

#endif
