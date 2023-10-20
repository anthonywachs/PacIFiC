#ifndef EXT_PETSc_API_HH
#define EXT_PETSc_API_HH

#include <MAC_ExternalAPI.hh>
#include <MAC_Timer.hh>
#include <MAC.hh>
#include <mpi.h>
class MAC_ModuleExplorer ;
// Some defs to allow PETSc tracelogs.
// #define PETSC_USE_DEBUG 1
// !!! IMPORTANT !!!
// the above line (#define PETSC_USE_DEBUG 1) needs to be commented when petsc 
// is compiled with --with-debugging=0, i.e., in optimized mode.
// Interestingly, this does not cause any trouble petsc-3.0.0, but it does 
// with version 3.2. Even #define PETSC_USE_DEBUG 0 does not work.
// Not very clear why ... the error we get is the following, at the linking
// stage of the code that use libmac0.so :
// libmac0.so: undefined reference to `petscstack'
// collect2: ld a retourné 1 code d'état d'exécution
// The problem is documented on the web by googling "undefined reference to
// `petscstack'" 
// !!! END IMPORTANT !!!
#define PETSC_USE_LOG 1
#define PETSC_USE_STACK 1
#include <mpi.h>

// From version 3.6, extern "C" causes odd compiling problems related to c++
// complex numbers. It turned out that extern "C" is not even needed.
//extern "C"
//{
#include <petscao.h>
#include "petscmat.h"
#include "petscksp.h"
#include "petscpc.h"
#include "petscversion.h"
//}

/*
PETSc applications, performing their specific initialization
and termination.

PUBLISHED
*/

#define PETSC_IMPLEMENTATION 3
class EXT_PETScAPI : public MAC_ExternalAPI
{

   public: //-----------------------------------------------------------

      static bool parse_options( MAC_ModuleExplorer const* exp,
                                 bool verbose ) ;
      static void going_to_do( char const* action ) ;
      static void verify( char const* action, int result ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      EXT_PETScAPI( void ) ;
     ~EXT_PETScAPI( void ) ;
      EXT_PETScAPI( EXT_PETScAPI const& other ) ;
      EXT_PETScAPI& operator=( EXT_PETScAPI const& other ) ;

      static MAC_Timer* timer ;
   //-- Current instance management

      virtual void initialize( int& argc, char **& argv ) ;

   //-- Class attributes

      static EXT_PETScAPI* SINGLETON ;
} ;

#define PETSc_do(X) { EXT_PETScAPI::going_to_do( #X ) ; { MAC_Marker pspy( #X ) ; EXT_PETScAPI::verify( #X, X ) ; } }

#endif



