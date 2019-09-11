#ifndef EXT_MPI_API_HH
#define EXT_MPI_API_HH

#include <MAC_ExternalAPI.hh>

/*
MPI applications, performing their specific initialization
and termination.

PUBLISHED
*/

class EXT_MPI_API : public MAC_ExternalAPI
{

   public: //-----------------------------------------------------------
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~EXT_MPI_API( void ) ;
      EXT_MPI_API( void ) ;
      EXT_MPI_API( EXT_MPI_API const& other ) ;
      EXT_MPI_API& operator=( EXT_MPI_API const& other ) ;

   //-- Current instance management
      
      virtual void initialize( int& argc, char **& argv ) ;

   //-- Class attributes

      static EXT_MPI_API* SINGLETON ;
} ;

#endif



