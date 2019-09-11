#ifndef EXT_MPI_COMMUNICATOR_HH
#define EXT_MPI_COMMUNICATOR_HH

#include <MAC_Communicator.hh>

/*
MPI communicator.

PUBLISHED
*/

class EXT_MPIcommunicator : public MAC_Communicator
{
   public: //-----------------------------------------------------------

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      EXT_MPIcommunicator( void ) ;
     ~EXT_MPIcommunicator( void ) ;
      EXT_MPIcommunicator( EXT_MPIcommunicator const& other ) ;
      EXT_MPIcommunicator& operator=( EXT_MPIcommunicator const& other ) ;
      
      enum { TAG_INT, TAG_DOUBLE, TAG_CHAR } ;

   //-- Characteristics
         
      virtual size_t nb_ranks( void ) const ;
         
      virtual size_t rank( void ) const ;
         
   //-- Point-to-point blocking communication

      virtual void send( size_t dest, int const* value, int nb ) const ;
      virtual void receive( size_t src, int* value, int nb ) const ;
      
      virtual void send( size_t dest, double const* value, int nb ) const ;
      virtual void receive( size_t src, double* value, int nb ) const ;

      virtual void send( size_t dest, char const* value, int nb ) const ;
      virtual void receive( size_t src, char* value, int nb ) const ;
      
   //-- Point-to-point non-blocking communication
      
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
				 size_t root ) const ; 			   
      
   //-- Collective communication: gather to all
      
      virtual void all_gather( int const* value, size_t nb, 
                               int* result ) const ;

      virtual void all_gather( double const* value, size_t nb, 
      	double* result ) const ;
      
      virtual void all_gather_v( double const* values,
                                 size_t nb,
                                 double* result,
                                 intVector const& partition,
                                 intVector const& start ) const ;

   //-- Collective communication: all to all scatter/gather
      
      virtual void all_to_all( int const* value, size_t nb, 
                               int* result ) const ;
      
   //-- Collective communication: reductions
         
      virtual bool boolean_and( bool value ) const ;

      virtual bool boolean_or( bool value ) const ;
      
      virtual double sum( double value ) const ;

      virtual void sum_vector( double* values, int nb ) const ;
      
      virtual double min( double value ) const ;

      virtual double max( double value ) const ;
      
      virtual bool same_value_everywhere( double val ) const ;

      virtual bool same_value_everywhere( int val ) const ;
      
      void reduce_vector( doubleVector& vec, size_t root = 0) const ;
      
      void reduce_vector_max( doubleVector& vec, size_t root = 0 ) const ;
      
      virtual size_t sum( size_t value ) const ;
      
      virtual unsigned long long int sum( unsigned long long int value ) 
      	const ; 
      
      virtual size_t min( size_t value ) const ;

      virtual size_t max( size_t value ) const ;            
      
   //-- Collective communication: synchronization
      
      virtual void barrier( void ) const ;
      
      virtual void wait( void* request ) const ;
      
   //-- Class attributes
      
      static EXT_MPIcommunicator* SINGLETON ;

   //-- Attributes
      
      mutable int count ;
} ;

#endif
