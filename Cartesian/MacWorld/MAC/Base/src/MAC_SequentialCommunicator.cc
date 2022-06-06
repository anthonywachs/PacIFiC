#include <MAC_SequentialCommunicator.hh>

#include <MAC_assertions.hh>


MAC_SequentialCommunicator*
MAC_SequentialCommunicator:: SINGLETON = new MAC_SequentialCommunicator() ;
   
//----------------------------------------------------------------------
MAC_SequentialCommunicator:: MAC_SequentialCommunicator( void )
//----------------------------------------------------------------------
   : MAC_Communicator( "MAC_SequentialCommunicator" )
{
}

//----------------------------------------------------------------------
MAC_SequentialCommunicator:: ~MAC_SequentialCommunicator( void )
//----------------------------------------------------------------------
{
   SINGLETON = 0 ;
}

//----------------------------------------------------------------------
size_t
MAC_SequentialCommunicator:: nb_ranks( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: nb_ranks" ) ;

   static size_t result = 1 ;
   
   MAC_CHECK_POST( nb_ranks_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
MAC_SequentialCommunicator:: rank( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: rank" ) ;

   static size_t result = 0 ;

   MAC_CHECK_POST( rank_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: send( size_t dest, int const* value, 
                                   int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: send( int const* )" ) ;
   MAC_CHECK( send_PRE( dest, value, nb ) ) ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: send( size_t dest, long long int const* value, 
                                   int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: send( long long int const* )" ) ;
   MAC_CHECK( send_PRE( dest, value, nb ) ) ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: send( size_t dest, double const* value, 
                                   int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: send( double const* )" ) ;
   MAC_CHECK( send_PRE( dest, value, nb ) ) ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: send( size_t dest, char const* value, 
                                   int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: send( char const* )" ) ;
   MAC_CHECK( send_PRE( dest, value, nb ) ) ;
}

//----------------------------------------------------------------------
void*
MAC_SequentialCommunicator:: Isend( size_t dest, int const* value,
                                    int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: Isend( int const*)" ) ;
   MAC_CHECK( send_PRE( dest, value, nb ) ) ;
   return( 0 ) ;
}

//----------------------------------------------------------------------
void*
MAC_SequentialCommunicator:: Isend( size_t dest, double const* value,
                                    int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: Isend( double const* )" ) ;
   MAC_CHECK( send_PRE( dest, value, nb ) ) ;
   return( 0 ) ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: receive( size_t src, int* value,
                                      int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: receive( int* )" ) ;
   MAC_CHECK( receive_PRE( src, value, nb ) ) ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: receive( size_t src, long long int* value,
                                      int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: receive( long long int* )" ) ;
   MAC_CHECK( receive_PRE( src, value, nb ) ) ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: receive( size_t src, double* value, 
                                      int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: receive( double* )" ) ;
   MAC_CHECK( receive_PRE( src, value, nb ) ) ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: receive( size_t src, char* value, 
                                      int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: receive( char* )" ) ;
   MAC_CHECK( receive_PRE( src, value, nb ) ) ;
}

//----------------------------------------------------------------------
void*
MAC_SequentialCommunicator:: Ireceive( size_t src, int* value,
                                       int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: Ireceive( int* )" ) ;
   MAC_CHECK( receive_PRE( src, value, nb ) ) ;
   return( 0 ) ;
}

//----------------------------------------------------------------------
void*
MAC_SequentialCommunicator:: Ireceive( size_t src, double* value,
                                       int nb  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: Ireceive( double* )" ) ;
   MAC_CHECK( receive_PRE( src, value, nb ) ) ;
   return( 0 ) ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: all_gather_v( double const* values,
                                           size_t nb,
                                           double* result,
                                           intVector const& partition,
                                           intVector const& start ) const 
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<nb ; ++i ) result[i] = values[i] ;   
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: all_gather(
                       int const* value, size_t nb, int* result  ) const
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<nb ; ++i ) result[i] = value[i] ;   
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: all_to_all( int const* value,
                                         size_t nb,
                                         int* result  ) const
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<nb ; ++i ) result[i] = value[i] ;   
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: gather( double const* value, size_t nb, 
                                     double* result, size_t root ) const
//----------------------------------------------------------------------
{
   for( size_t i=0 ; i<nb ; ++i ) result[i] = value[i] ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: gather( int const* value, size_t nb, 
                                     int* result, size_t root ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: gather(int)" ) ;

   for( size_t i=0 ; i<nb ; ++i ) result[i] = value[i] ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: broadcast(
                              int* value,size_t nb, size_t root ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: broadcast(
                            double* value,size_t nb, size_t root ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: broadcast(
                              char* value,size_t nb, size_t root ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: barrier( void ) const
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
bool
MAC_SequentialCommunicator:: same_value_everywhere( double val ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_SequentialCommunicator:: same_value_everywhere( int val ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: all_gather( double const* value, size_t nb, 
	double* result  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: all_gather(double)" ) ;

   for(size_t i=0 ; i<nb ; i++ )
      result[i] = value[i] ;
   
}

//----------------------------------------------------------------------
void
MAC_SequentialCommunicator:: gather_v( double const* values,
                                           size_t nb,
                                           double* result,
                                           intVector const& partition,
                                           intVector const& start,
					   size_t root) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SequentialCommunicator:: gather_v(int)" ) ;

   for(size_t i=0 ; i<nb ; i++ )
      result[i] = values[i] ;
   
}
