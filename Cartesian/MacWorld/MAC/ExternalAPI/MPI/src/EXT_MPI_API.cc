#include <EXT_MPI_API.hh>

#include <MAC_Bool.hh>
#include <MAC_Exec.hh>
#include <MAC_Error.hh>
#include <MAC_String.hh>
#include <MAC_System.hh>
#include <MAC_Variable.hh>
#include <MAC_assertions.hh>

// Both stdio.h and the MPI C++ interface use SEEK_SET, SEEK_CUR, SEEK_END.
// This is really a bug in the MPI-2 standard.
// A possibility would be to undefine the 3 names SEEK_SET, SEEK_CUR, SEEK_END
//    #undef SEEK_SET
//    #undef SEEK_CUR
//    #undef SEEK_END
// Our solution is to define MPICH_IGNORE_CXX_SEEK which works at least 
// with MPICH2
#define MPICH_IGNORE_CXX_SEEK 1
#include <mpi.h>

#include <unistd.h>

#include <iostream>
#include <sstream>

using std::ostringstream ;
using std::endl ;

EXT_MPI_API* EXT_MPI_API:: SINGLETON = new EXT_MPI_API() ;

struct EXT_MPI_API_ERROR
{
   static void n0( std::string const& func ) ;
} ;


//----------------------------------------------------------------------
EXT_MPI_API:: EXT_MPI_API( void )
//----------------------------------------------------------------------
   : MAC_ExternalAPI( "EXT_MPI_API", 8 )
{   
}




//----------------------------------------------------------------------
EXT_MPI_API:: ~EXT_MPI_API( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_MPI_API:: ~EXT_MPI_API" ) ;
   int mpierr ;

   int size ;
   mpierr = MPI_Comm_size( MPI_COMM_WORLD, &size ) ;
   if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Comm_size" ) ;

   int err = MAC_Exec::exit_code() ;
   if( size>1 && err!=0 )
   {
      int rank ;

      mpierr = MPI_Comm_rank( MPI_COMM_WORLD, &rank ) ;
      if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Comm_rank" ) ;

      std::cout << "Parallel execution interrupted by processor #" << rank 
                << "." << std::endl ;

      mpierr = MPI_Abort( MPI_COMM_WORLD, err ) ;
      if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Abort" ) ;
   }

   mpierr = MPI_Finalize() ;
   if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Finalize" ) ;
}




//----------------------------------------------------------------------
void
EXT_MPI_API:: initialize( int& argc, char **& argv )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_MPI_API:: initialize" ) ;
   
   int mpierr ;

   std::string const cwd = MAC_System::working_directory() ;

   mpierr = MPI_Init( &argc, &argv ) ;
   if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Init" ) ;

   int size ;
   mpierr = MPI_Comm_size( MPI_COMM_WORLD, &size ) ;
   if( mpierr != MPI_SUCCESS ) EXT_MPI_API_ERROR::n0( "MPI_Comm_size" ) ;

   if( size == 1 )
   {
      MAC_System::changedir( cwd ) ;
   }
  
   MAC_Exec::add_variable_to_execution_context(
                      MAC_Variable::object( "BS_with_MPI" ),
                      MAC_Bool::create( 0, true ) ) ;
   
#ifndef MPIRUN
# error \
Macro MPIRUN must be set when compiling Pelicans (opt: -DMPIRUN=<value>).
#else
   if( !MAC_System::can_read( MPIRUN ) )
   {
      ostringstream mesg ;
      mesg << "*** EXT_MPI_API:" << endl ;
      mesg << "    unable to read the file called" << endl ;
      mesg << "       \"" << MPIRUN << "\""  << endl << endl ; 
      mesg << "This file was specified by the macro MPIRUN when" << endl ;
      mesg << "compiling MAC, eg via the compiler option:" << endl ;
      mesg << "  -DMPIRUN=" << MPIRUN << endl ;
      mesg << "or in the extra-makefile with an instruction such as:" << endl ;
      mesg << "  MPIRUN=" << MPIRUN << endl << endl ;
      mesg << "The macro MPIRUN should contain the full path of the" << endl ;
      mesg << "command used to run an MPI application" ;
      MAC_Error::object()->raise_plain( mesg.str() ) ;
   }
   MAC_Exec::add_variable_to_execution_context(
                      MAC_Variable::object( "SS_MPI_RUN" ),
                      MAC_String::create( 0, MPIRUN ) ) ;
#undef MPIRUN
#endif
}




//internal---------------------------------------------------------------
void
EXT_MPI_API_ERROR:: n0( std::string const& func )
//internal---------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** EXT_MPI_API:" << endl ;
   mesg << "    call to " << func << " failed" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}
