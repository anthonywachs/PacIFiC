#include <MAC_System.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <stringVector.hh>
#include <size_t_vector.hh>

#include <dir.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <dir.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <new>


#ifdef _WIN32
#include <direct.h>
#include <io.h>
#include <process.h>
#include <tchar.h>
#include <Windows.h>
#ifdef max
#undef max
#endif
#define SYSCALL(X) _##X
#else
#include <sys/utsname.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#define SYSCALL(X) ::X
#endif

#include <time.h>

#define MAC_SIGNAL_HANDLING 0

#ifdef _WIN32

//no_doc----------------------------------------------------------------------
static unsigned int old_exponent_format = 
     _set_output_format(_TWO_DIGIT_EXPONENT);
//no_doc----------------------------------------------------------------------

//no_doc----------------------------------------------------------------------
MAC_System::Process::Process( stringVector const& args )
//no_doc----------------------------------------------------------------------
  : running(true)
{
   MAC_ASSERT( args.size() > 0 ) ;

   STARTUPINFO si ;
   std::string line ;
   for( size_t i=0 ; i<args.size() ; i++ )
      line += args(i) + " " ;
         
   LPTSTR szCmdline=_tcsdup(line.c_str()) ;

   ZeroMemory( &si, sizeof(si) );
   si.cb = sizeof(si);

   PROCESS_INFORMATION *pi = new PROCESS_INFORMATION() ;
   ZeroMemory( pi, sizeof(PROCESS_INFORMATION) );
   
   // Start the child process. 
   if( !CreateProcess( NULL,   // No module name (use command line)
                       szCmdline,      // Command line
                       NULL,           // Process handle not inheritable
                       NULL,           // Thread handle not inheritable
                       FALSE,          // Set handle inheritance to FALSE
                       0,              // No creation flags
                       NULL,           // Use parent's environment block
                       NULL,           // Use parent's starting directory 
                       &si,            // Pointer to STARTUPINFO structure
                       pi )            // Pointer to PROCESS_INFORMATION structure
      ) 
   {
      MAC_Error::object()->raise_plain(
         "MAC_CoupledApplications : unable to start " + line ) ;
   }
   sdef = (void*) pi ;
   pid = pi->dwProcessId ;
}

//no_doc----------------------------------------------------------------------
size_t MAC_System:: Process::id( void ) const
//no_doc----------------------------------------------------------------------
{   
   return pid ;
}

//no_doc----------------------------------------------------------------------
void MAC_System:: Process::wait_for_child_processes( size_t nb, Process ** child )
//no_doc----------------------------------------------------------------------
{
   HANDLE * handles = new HANDLE[nb] ;
   
   for( size_t i=0 ; i<nb ; i++ )
   {
      size_t n = 0 ;
      size_t_vector index(nb) ;
      for( size_t p=0 ; p<nb ; p++ )
      {
         if( child[p]->running )
         {
            index(n) = p ;
            PROCESS_INFORMATION *pi =  ( PROCESS_INFORMATION* )child[p]->sdef ;
            handles[n++] = pi->hProcess ;
         }
      }
      int idx = WaitForMultipleObjects( n, handles, false, INFINITE );
      idx -= WAIT_OBJECT_0 ;
      MAC_ASSERT( idx>=0 && idx< n ) ;

      int lidx = index(idx) ;
      MAC_ASSERT( child[lidx]->running ) ;

      DWORD lpExitCode = 0 ;
      PROCESS_INFORMATION *pi = ( PROCESS_INFORMATION* )child[lidx]->sdef ;
      MAC_ASSERT( GetExitCodeProcess( pi->hProcess, &lpExitCode ) != 0 ) ;
	  int ret = lpExitCode ;
      MAC::out() << "Exit code " << ret 
                 << " returned by " << child[lidx]->name() << std::endl ;
      if( ret!=0 ) 
      {
        MAC::out() << "Exit code " << ret << std::endl ;
        MAC_Error::object()->raise_plain(
            child[lidx]->name() + " terminated unsuccessfuly" ) ;
      }
      child[lidx]->running = false ;
   } 
   delete [] handles ;

}

   
//no_doc----------------------------------------------------------------------
MAC_System:: Process::~Process(void)
//no_doc----------------------------------------------------------------------
{
   int ret ;
   PROCESS_INFORMATION *pi = ( PROCESS_INFORMATION* )sdef ;
   if( running ) TerminateProcess( pi->hProcess, ret ) ;
   delete pi ;
}

#else

//no_doc----------------------------------------------------------------------
MAC_System:: Process::Process( stringVector const& args )
//no_doc----------------------------------------------------------------------
      : running(true)
{
   MAC_ASSERT( args.size() > 0 ) ;

   pid = fork() ;
   if( pid == 0 ) 
   {
      char ** argv = new char * [args.size()+1] ;
      for( size_t i=0 ; i<args.size() ; i++ )
      {
         const char * str = args(i).c_str() ;
         size_t n =  args(i).length() +1 ;
               
         argv[i] = new char [ n ] ;
         for( size_t c=0 ; c<n ; c++ )
            argv[i][c] = str[c] ;
         
      }
      argv[args.size()] = 0 ;
            
      int ret = execv( argv[0],
                       argv ) ;
      if( ret!=0 ) 
      {
         MAC_Error::object()->raise_plain(
            "Unable to start " + args(0) ) ;
      }
   }
}

//no_doc----------------------------------------------------------------------
size_t MAC_System:: Process::id( void ) const
//no_doc----------------------------------------------------------------------
{   
   return pid ;
}

//no_doc----------------------------------------------------------------------
void MAC_System:: Process::wait_for_child_processes( size_t nb, Process ** child )
//no_doc----------------------------------------------------------------------
{
   for( size_t i=0 ; i<nb ; i++ )
   {
      int status ;
      size_t a_id = wait(&status) ;
      size_t idx = nb + 1 ;
      for( size_t j=0 ; j<nb ; j++ )
      {
         if( child[j]->pid == a_id ) idx=j ;
      }
      MAC_ASSERT( idx < nb ) ;
      child[idx]->running = false ;
      bool ok = WIFEXITED(status) && WEXITSTATUS(status)==0 ;
      MAC::out() << "Status " << ( ok ? "successful" : "unsuccessful" )
                 << " returned by " << child[idx]->name() << std::endl ;
      if( !ok ) 
      {
         MAC_Error::object()->raise_plain(
            "Process terminated unsuccessfuly" ) ;
      }
   } 
}


//no_doc----------------------------------------------------------------------
MAC_System:: Process::~Process(void)
//no_doc----------------------------------------------------------------------
{
   if( running ) kill( pid, 9 ) ;
}
      
#endif

//no_doc----------------------------------------------------------------------
void
MAC_System::Process::set_name( std::string a_name )
//no_doc----------------------------------------------------------------------
{
   my_name = a_name ;
}

//no_doc----------------------------------------------------------------------
std::string const& 
MAC_System::Process::name( void ) const
//no_doc----------------------------------------------------------------------
{
   return my_name ;
}

//-------------------------------------------------------------------------
std::string const&
MAC_System:: host_name( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::host_name" ) ;
   static std::string result ;
   if( result.empty() )
   {
#ifndef _WIN32
      char host[256] ;
      gethostname( host, 256 ) ;
      result = host ;
#else
      // Cette fonction existe dans .net : System::Dns:gethostname
      result = "unknown host" ;
#endif   
   }
   
   return( result ) ;
}

//-------------------------------------------------------------------------
int
MAC_System:: process_id( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::process_id" ) ;
   static int result = -1 ;
   if( result==-1 )
   {
       result = SYSCALL(getpid)() ;
   }
   
   return result ;
}

//-------------------------------------------------------------------------
char
MAC_System:: path_name_separator( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::path_name_separator" ) ;
#ifndef _WIN32 
   return '/' ;
#else
   return '\\' ;
#endif
}

//-------------------------------------------------------------------------
bool 
MAC_System:: matches( std::string const& str,
                      std::string const & pattern )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::matches" ) ;
   
   bool result = ( pattern=="*" ) || ( str == pattern ) ;
   if( !result ) 
   {
      size_t idx = pattern.find( "*" ) ;
      if( idx==0 && str.length()>=pattern.length()-1 ) 
      {
         std::string res = pattern.substr( 1, pattern.length()-1 ) ;
         result = str.substr( str.length()-res.length(), str.length()-1 )== res ;
      }
      else if ( idx==pattern.length()-1 )
      {
         std::string res = pattern.substr( 0, pattern.length()-1 ) ;
         result = str.find( res ) == 0 ;
      }
   }
//   MAC::out() << str << "<" << result << ">" << pattern << std::endl ;
   
   return result ;
   
}
   
//-------------------------------------------------------------------------
void 
MAC_System:: find( std::string const& directory,
                   std::string const& filename,
                   stringVector& result, 
                   bool recurse,
                   std::string const& prefix )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::find" ) ;
   
   static int level = 0 ;
   level++ ;
   if( level >= 1024 ) 
   {
      MAC::out() << "Maximal recursion level reached" << std::endl <<
         " Directory contains probably cyclic path : " << std::endl <<
         prefix << std::endl ;
   }
   else
   {
      std::string complete = directory ;
      if( !prefix.empty() ) complete = prefix + path_name_separator() + complete ;
      DIR* dir = opendir( complete.c_str() ) ;
      if( dir!=0 ) 
      {
         struct dirent* item ;
         while( (item=readdir(dir)) ) 
         {
            std::string f = item->d_name ;
            if( f !="." && f!=".." )
            {
               struct stat st ;
               std::string absolute_path = complete+ path_name_separator() + f ;
               
               int ret = stat( absolute_path.c_str(), &st ) ;
               
               int mode = st.st_mode ;
               if( ret == 0 ) 
               {
                  bool is_regular = ( mode & S_IFREG) ;
                  bool is_directory = ( mode & S_IFDIR ) ;
                  
                  if( matches( f, filename ) && is_regular ) 
                  {
                     result.extend( absolute_path ) ;
                  }
                  else if( recurse && is_directory )
                  {
                     find( f, filename, result, recurse, complete ) ;   
                  }
               }
            }
         }
         closedir(dir) ;
      }
   }
   level-- ;
   if( level == 0 )
   {
      result.sort() ;
   }
}

//-------------------------------------------------------------------------
int 
MAC_System:: run( stringVector const& cmd, 
                  std::string const& directory,
                  std::string const& output_file )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: mac_run" ) ;
   std::string const here = MAC_System::working_directory() ;
   MAC_System::changedir( directory ) ;
   
//    if( fork()==0 )
//    {
//       char *const argv = new char* [cmd.size()+1] ;
//       for( size_t i=0 ; i<cmd.size() ; i++ )
//       {
//          argv[i] = new char [ cmd(i).length() ] ;
//          strncpy( argv[i], cmd(i).c_str(), cmd(i).length() ) ;
//       }
//       argv[cmd.size()]=0 ;
      
//       execv( argv[0], argv ) ;
//    }
   std::string line ;
   for( size_t i=0 ; i<cmd.size() ; i++ )
   {
      line += MAC_System::string( cmd(i) )+" " ;
   }
   if( !output_file.empty() )
      line += " > " + MAC_System::string( output_file ) ;

#ifdef _WIN32
   line = "\"" + line + "\"" ;
#endif
   
   int result = system( line.c_str() ) ;
   
   MAC_System::changedir( here ) ;
   return result ;
   
}

//-------------------------------------------------------------------------
bool
MAC_System:: mkdir( std::string const& directory )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: mkdir" ) ;
   DIR * dir = opendir( directory.c_str() ) ;
   int sys_res = 0 ;
   if( dir!=0 ) 
   {
      closedir(dir) ;
   }
   else
   {
      std::string line ;
      line = "mkdir " ;
#ifndef _WIN32
      line += "-p " ;
#endif
      line += MAC_System::string( directory ) ;
      sys_res = system( line.c_str() ) ;
   }

   bool result = ( sys_res == 0 ) ;
   return result ;
}

//-------------------------------------------------------------------------
bool 
MAC_System:: changedir( std::string const& directory )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: changedir" ) ;

   int sys_res = SYSCALL(chdir)( directory.c_str() ) ;
   bool result = ( sys_res == 0 ) ;

   return( result ) ;
}


//-------------------------------------------------------------------------
bool
MAC_System:: copy( std::string const& src, std::string const& dest )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: copy" ) ;
   
   std::string line ;
#ifndef _WIN32
   line = "cp " ;
#else
   line = "copy " ;
#endif
   line += MAC_System::string( src ) ;
   line += " " ;
   line += MAC_System::string( dest ) ;
   int sys_res = system( line.c_str() ) ;
   

   bool result = ( sys_res == 0 ) ;
   return result ;
}

//-------------------------------------------------------------------------
bool
MAC_System:: erase( std::string const& filename )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: erase" ) ;

   return( remove( filename.c_str() ) == 0 ) ;
}
        
//-------------------------------------------------------------------------
std::string const&
MAC_System:: working_directory( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::working_directory" ) ;
   static std::string result  ;
   
   char tmp[2048] ;

   MAC_ASSERT( SYSCALL(getcwd)(tmp,2048)!=0 ) ;

   result = tmp ;
   
   return result ;
}

//-------------------------------------------------------------------------
std::string const&
MAC_System:: sysname( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::sysname" ) ;
   static std::string result  ;
   if( result.empty() )
   {
#ifdef _WIN32
      result = "windows" ;
#else
      struct utsname u;
      uname(&u) ;
      result = u.sysname ;
#endif
   }
   
   return result ;   
}

//-------------------------------------------------------------------------
double
MAC_System:: user_time( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::user_time" ) ;
   static double result = -1.0 ;
#ifndef _WIN32
   static struct rusage r;
   static double last_time = 0.0 ;
   
   getrusage(RUSAGE_SELF, &r);
   result = r.ru_utime.tv_sec + 1.e-6 * r.ru_utime.tv_usec ;
   result = MAC::max( last_time, result ) ;
   last_time = result ;
#else
   result = 0.001 * GetTickCount() ;
#endif
   return( result ) ;   
}

//-------------------------------------------------------------------------
double
MAC_System:: epoch_time( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::epoch_time" ) ;
   
   return( (double)time(0) ) ;   
}


//-------------------------------------------------------------------------
void
MAC_System:: sleep( size_t msec )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::sleep" ) ;
   
#ifndef _WIN32
   usleep( 1000*msec ) ;   
#else
   Sleep( msec ) ;
#endif
      
}

     
//-------------------------------------------------------------------------
size_t
MAC_System:: used_memory( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: used_memory" ) ;
   size_t result = 0 ;
   
#ifndef _WIN32
   if( sysname()!="Linux" )
   {
      static struct rusage r;
   
      getrusage( RUSAGE_SELF, &r );
      result = r.ru_idrss * getpagesize() ;
   }
   else
   {
      std::ostringstream os ;
      
      std::string word ;
      os << "/proc/" << process_id() << "/status" ;
      
      std::ifstream in( os.str().c_str() ) ;
      if( !in )
      {
         std::cout << "MAC_System : Unable to open "
                   << os.str() << std::endl ;
      }
      else
      {
         while(!in.eof())
         {
            in >> word ;
            if( word == "VmSize:" )
            {
               in >> result ;
               in >> word ;
               if( !( word == "kB" ) )
               {
                  MAC_Error::object()->raise_plain( "Unit is "+word ) ;
               }
               result *= 1000 ;
               break ;
            }
         }
         in.close() ;
      }
   }
#endif
   return( result ) ;
}

//------------------------------------------------------------------------
std::string const&
MAC_System:: compiler_name( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::compiler_name" ) ;
   static std::string const result =
#ifdef __SUNPRO_CC
   "CC" ;
#endif
#ifdef __GNUC__
   "gcc" ;
#endif
#ifdef __ICC
   "icc" ;
#endif
#ifdef __DECCXX
   "cxx" ;
#endif
#ifdef __XLC121__
   "xlC" ;
#endif
#ifdef __PGI
   "pgCC" ;
#endif
#ifdef _WIN32
   "VC++" ;
#endif
   return( result ) ;
}

//----------------------------------------------------------------------
std::string
MAC_System:: basename( std::string const& path_and_name,
                       char separator )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: basename" ) ;
   size_t idx = path_and_name.find_last_of( separator ) ;
   std::string result ;
   if( idx==0 ) 
   {
      result = separator ;
   }
   else if( idx==path_and_name.length()-1 )
   {
      result = path_and_name.substr( 0, path_and_name.length()-1 ) ;
   }
   else if( idx<path_and_name.length() )
   {
      result = path_and_name.substr( idx+1, path_and_name.length()-idx-1 ) ;
   }
   else
   {
      result = path_and_name ;
   }
   return( result ) ;
}



//----------------------------------------------------------------------
std::string
MAC_System:: dirname( std::string const& path_and_name,
                      char separator )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: dirname" ) ;
   size_t idx = path_and_name.find_last_of( separator ) ;
   std::string result ;
   if( idx==0 ) 
   {
      result = separator ;
   }
   else if( idx>0 && idx<path_and_name.length() )
   {
      result = path_and_name.substr( 0, idx ) ;
   }
   else
   {
      result = "." ;
   }
   return( result ) ;
}

//---------------------------------------------------------------------------
std::string
MAC_System:: absolute_path( std::string const& filename ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: absolute_path" ) ;
   MAC_CHECK_PRE( filename.length() > 0 ) ;

   std::string result = filename ;
   char sep = path_name_separator() ;
   
   bool relative = ! ( filename.at(0)==sep || sep == '\\' &&
                       ( filename.length()>1 && filename.at(1)==':' ) ) ;
   
   if( relative ) // relative path
      result = working_directory() + sep + result ;

   stringVector vec( result, sep ) ;
   for( size_t i=1 ; i<vec.size() ; i++ )
   {
      if( vec(i).empty() ) vec.remove_at(i) ;
   }
   
   bool replacement = true ;
   while( replacement )
   {
      replacement = false ;
      for( size_t i=0 ; i<vec.size() ; i++ )
      {
         if( vec(i)=="." )
         {
            vec.remove_at( i ) ;
            replacement = true ;
            break ;
         }
         if( vec(i)==".." && i>0 )
         {
            vec.remove_at( i ) ; vec.remove_at( i-1 ) ;
            replacement = true ;
            break ;
         }
         
      }
   }
   result = vec(0) ;
   for( size_t i=1 ; i<vec.size() ; i++ ) 
      result = result + sep + vec(i) ;

   return result ;
}

#if MAC_SIGNAL_HANDLING==0

extern "C" 
{
   void mac_signal_handler( int sig )
   {
      std::string mess ;
      switch (sig) 
      {
         case SIGSEGV : mess = "Segmentation violation\n Illegal memory access\n (a null pointer may have been dereferenced)" ; 
            break ;
#ifndef _WIN32
         case SIGBUS : mess = "Bus Error\n Illegal memory access" ; 
            break ;
#endif
         case SIGFPE : mess = "Floating Point Exception, probably divide by zero" ; 
            break ;
         case SIGINT : mess = "Interrupted by user" ; 
            break ;
            
         default : mess = "Unrecognized signal" ;
      }
      MAC_Error::object()->raise_internal(
         "Pelicans signal handler :\n" + mess ) ;
   }
} 

void newHandler( void )
{
   MAC_Error::object()->raise_internal( "unable to satisfy memory request" ) ;
}

//------------------------------------------------------------------------
void MAC_System::exception_trapping( void )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::exception_trapping" ) ;
   
   std::set_new_handler( newHandler ) ;
#ifndef _WIN32
   signal( SIGBUS, mac_signal_handler) ;
#endif
   signal( SIGFPE, mac_signal_handler) ;
   signal( SIGSEGV, mac_signal_handler) ;
   if( MAC_Exec::communicator()->nb_ranks()==1 )
   {
      // WARNING:
      //   SIGINT is used by MPI to kill slave processes
      //   error if the signal is trapped
      signal( SIGINT, mac_signal_handler) ;
   }  
}

#endif

#if MAC_SIGNAL_HANDLING==1
extern "C" 
{
   typedef void (*MAC_sighandler_t)(int);
   #define NBSIGNALMAX NSIG
   MAC_sighandler_t signals[NBSIGNALMAX] ;
   
   void mac_signal_handler( int sig )
   {
      static int prem = 0 ;
      prem++ ;
      std::ostringstream mess ;
      bool end_execution = false ;
      switch( sig ) 
      {
         case SIGSEGV : 
            mess << "Segmentation violation\n Illegal memory access\n (a null pointer may have been dereferenced)" ;
            end_execution = true ;
            break ;
#ifndef _WIN32
         case SIGBUS : 
            mess << "Bus Error\n Illegal memory access" ; 
            end_execution = true ;
            break ;
#endif
         case SIGFPE : 
            mess << "Floating Point Exception, probably divide by zero" ; 
            end_execution = true ;
            break ;
         case SIGINT : 
            mess << "Interrupted by user" ; 
            end_execution = true ;
            break ;
         default : 
            mess << " signal received " << sig ;
      }
      MAC_sighandler_t handler = signals[sig] ;
      if( end_execution )
      {
         if( prem == 1 ) 
            MAC_Error::object()->raise_internal(
               "Pelicans signal handler :\n" + mess.str() ) ;
      }
      else
      {
         MAC_Error::object()->trace( mess.str() ) ;
         if( handler != 0 ) handler( sig ) ;
      }
   }
} 

void newHandler( void )
{
   MAC_Error::object()->raise_internal( "unable to satisfy memory request" ) ;
}

//no_doc------------------------------------------------------------------
void MAC_System::exception_trapping( void )
//no_doc------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System::exception_trapping" ) ;
   std::set_new_handler( newHandler ) ;

   // Interrupt all signals
   for( size_t i=0 ; i<NBSIGNALMAX ; i++ )
   {
#ifndef _WIN32
      if( i!=SIGCHLD ) //PETSc
#endif
         signals[i] = signal( i, mac_signal_handler ) ;
   }
   
#ifndef _WIN32
   signal( SIGBUS, mac_signal_handler) ;
#endif
}

#endif

//------------------------------------------------------------------------
size_t 
MAC_System:: new_block_size( size_t current_size, size_t expected_size )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: new_block_size" ) ;
   MAC_CHECK_PRE( expected_size > current_size ) ;

   size_t result = expected_size ;

   if( expected_size > 4 )
   {
      size_t dd = 2 * current_size ;
      if( dd > expected_size ) result = dd ;
   }
    
   MAC_CHECK_POST( result >= expected_size ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_System:: big_endian_encoding( void )
//----------------------------------------------------------------------
{
   bool result ;
   int integer = 1 ;
   char * c = (char*) &integer ;
   result = *c!=(char)1 ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_System:: can_read( std::string const& filename )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: can_read" ) ;
   struct stat st ;
   int ret = stat( filename.c_str(), &st ) ;
   bool result = ret==0 ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_System:: can_write( std::string const& filename )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: can_write" ) ;
   
   std::string dummy = dirname( filename ) ;
   dummy += path_name_separator() ;
   dummy += "_dummy_file_" ;
    
   std::ofstream out( dummy.c_str() ) ;
//   bool result = out!=0 ; // Not accepted from gcc-9.x.x
   bool result = out.is_open() ;   
   if( result ) 
   {
      out.close() ;
      MAC_System::erase( dummy ) ;
   }
   return( result );
}

//----------------------------------------------------------------------
std::string const&
MAC_System:: string( std::string const& cmd )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: string" ) ;

   static std::string result ;
   static char const s = '\"' ; // For blanks or special characters in paths
   bool alpha = true ;
   
   for( size_t i=0 ; i<cmd.length() ; i++ )
      alpha = alpha && isalpha( cmd[i] ) ;
   
   if( cmd[0] == s || alpha )
   {
      result = cmd ;
   } 
   else
   {
      result = s+cmd+s ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
int
MAC_System:: open_file_descriptor( std::string filename, int flags )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: open_file_descriptor" ) ;
   
   return( SYSCALL(open)( filename.c_str(), flags, S_IREAD ) );
}

//----------------------------------------------------------------------
int
MAC_System:: close_file_descriptor( int fildes )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: close_file_descriptor" ) ;
   

   return( SYSCALL(close)( fildes ) );
}

//----------------------------------------------------------------------
void
MAC_System:: exit( int exit_status )
//----------------------------------------------------------------------
{
   SYSCALL(exit)(exit_status) ;
}

//----------------------------------------------------------------------
std::string
MAC_System:: getenv( std::string const& str )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_System:: getenv" ) ;
   
   std::string result = "" ;
   
   char * res = ::getenv(str.c_str()) ;
   if( res!=0 ) result = res ;
   
   return result ;
}



