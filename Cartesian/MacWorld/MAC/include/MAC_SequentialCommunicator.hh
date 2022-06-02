#ifndef MAC_SEQUENTIAL_COMMUNICATOR_HH
#define MAC_SEQUENTIAL_COMMUNICATOR_HH

#include <MAC_Communicator.hh>

/*
Sequential communicator.

PUBLISHED
*/

class MAC_SequentialCommunicator : public MAC_Communicator
{
   public: //-----------------------------------------------------------

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      MAC_SequentialCommunicator( void ) ;
     ~MAC_SequentialCommunicator( void ) ;
      MAC_SequentialCommunicator( MAC_SequentialCommunicator const& other ) ;
      MAC_SequentialCommunicator& operator=(
                                  MAC_SequentialCommunicator const& other ) ;
      
   //-- Characteristics

      // IMPLEMENTATION : 1
      virtual size_t nb_ranks( void ) const ;

      // IMPLEMENTATION : 0
      virtual size_t rank( void ) const ;
      
   //-- Point-to-point blocking communication
      
      virtual void send( size_t dest, int const* value, int nb ) const ;
      virtual void receive( size_t src, int* value, int nb ) const ;

      virtual void send( size_t dest, long long int const* value, int nb ) 
      	const ;
      virtual void receive( size_t src, long long int* value, int nb ) const ;
      
      virtual void send( size_t dest, double const* value, int nb ) const ;
      virtual void receive( size_t src, double* value, int nb ) const ;
      
      virtual void send( size_t dest, char const* value, int nb ) const ;
      virtual void receive( size_t src, char* value, int nb ) const ;
      
   //-- Point-to-point non-blocking communication(501.)
      
      virtual void* Isend( size_t dest, int const* value, int nb ) const ;
      virtual void* Ireceive( size_t src, int* value, int nb ) const ;
      
      virtual void* Isend( size_t dest, double const* value, int nb ) const ;
      virtual void* Ireceive( size_t src, double* value, int nb ) const ;
      
   //-- Collective communication: broadcast
      
      virtual void broadcast( int* value, size_t nb, size_t root ) const ;
      
      virtual void broadcast( double* value, size_t nb, size_t root ) const ;
      
      virtual void broadcast( char* value, size_t nb, size_t root ) const ;
      
   //-- Collective communication: gather
      
      virtual void gather( double const* value, size_t nb, 
                           double* result, size_t root ) const ;
      
      virtual void gather( int const* value, size_t nb, 
                           int* result, size_t root ) const ;

      virtual void gather_v( double const* values,
                                 size_t nb,
                                 double* result,
                                 intVector const& partition,
                                 intVector const& start,
				 size_t root) const ;
      
   //-- Collective communication: gather to all
      
      virtual void all_gather( int const* value, size_t nb, 
                               int* result ) const ;
      
      virtual void all_gather_v( double const* values,
                                 size_t nb,
                                 double* result,
                                 intVector const& partition,
                                 intVector const& start ) const ;
	
      virtual void all_gather( double const* value, size_t nb, 
      	double* result  ) const ;	
      
   //-- Collective communication: all to all scatter/gather
      
      virtual void all_to_all( int const* value, size_t nb,
                               int* result ) const ;

   //-- Collective communication: reductions
      
      virtual bool same_value_everywhere( double val ) const ;

      virtual bool same_value_everywhere( int val ) const ;
      
   //-- Collective communication: synchronization
      
      virtual void barrier( void ) const ;
      
   //-- Class attributes
      
      static MAC_SequentialCommunicator* SINGLETON ;
} ;

#endif
