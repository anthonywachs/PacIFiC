#include <MAC_Timer.hh>

#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_System.hh>
#include <MAC_assertions.hh>

#include <iostream>

//-------------------------------------------------------------------------
MAC_Timer* MAC_Timer:: create( MAC_Object* a_owner )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Timer:: create" ) ;

   MAC_Timer* result = new MAC_Timer( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
MAC_Timer:: MAC_Timer( MAC_Object* a_owner )
//-------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , CUMUL_TIME( 0. )
   , CURRENT_TIME( MAC::bad_double() )
   , RUNNING( false )
   , CUMUL_ELAPSED_TIME( 0. )
   , CURRENT_ELAPSED_TIME( MAC::bad_double() )
{
   MAC_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
MAC_Timer:: ~MAC_Timer( void )
//-------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   if( is_running() )
   {
      stop() ;
   }
}

//-------------------------------------------------------------------------
void
MAC_Timer:: start( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Timer:: start" ) ;
   MAC_CHECK_PRE( !is_running() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   CURRENT_TIME = MAC_System::user_time() ;
   CURRENT_ELAPSED_TIME = MAC_System::epoch_time() ;
   RUNNING = true ;

   MAC_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
void
MAC_Timer:: stop( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Timer:: stop" ) ;
   MAC_CHECK_PRE( is_running() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   CUMUL_TIME += ( MAC_System::user_time() - CURRENT_TIME ) ;
   CURRENT_TIME = MAC::bad_double() ;
   CUMUL_ELAPSED_TIME += ( MAC_System::epoch_time() - CURRENT_ELAPSED_TIME ) ;
   CURRENT_ELAPSED_TIME = MAC::bad_double() ;
   RUNNING = false ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( !is_running() ) ;
}

//-------------------------------------------------------------------------
void
MAC_Timer:: reset( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Timer:: reset" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   if( is_running() )
   {
      stop() ;
   }
   CUMUL_TIME = 0. ;
   CUMUL_ELAPSED_TIME = 0. ;
   
   MAC_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
double
MAC_Timer:: time( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Timer:: time" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   double result = CUMUL_TIME ;
   if( is_running() )
   {
      result += ( MAC_System::user_time() - CURRENT_TIME ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( result>=0.0 ) ;
   return( result );
}
   
//-------------------------------------------------------------------------
double
MAC_Timer:: elapsed_time( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Timer:: elapsed_time" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   double result = CUMUL_ELAPSED_TIME ;
   if( is_running() )
   {
      result += ( MAC_System::epoch_time() - CURRENT_ELAPSED_TIME ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( result>=0.0 ) ;
   return( result );
}
   
//-------------------------------------------------------------------------
bool
MAC_Timer:: is_running( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Timer:: is_running" ) ;
   MAC_CHECK_INV( invariant() ) ;
   return( RUNNING ) ;
}

//-------------------------------------------------------------------------
void
MAC_Timer:: print( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   static size_t size = MAC_Exec::communicator()->nb_ranks() ;
   
   print_time( time(), os, indent_width ) ;
   if( size>1 ) 
   {
      os << " (" ;
      print_time( elapsed_time(), os, 0 ) ;
      os << ")" ;
   }
}

//-------------------------------------------------------------------------
void
MAC_Timer:: print_time( double a_time, std::ostream& os, size_t indent_width ) 
//-------------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space ;
   double t = a_time ;
   
   int const hh = (int) ( t/3600. ) ;
   int const mm = (int) ( (t-3600.*hh)/60. ) ;
   if( hh==0 && mm==0 )
   {
      os << t << " s" ;
   }
   else
   {
      if( hh > 0 )
      {
         os << hh << ":" ;
      }
      if( mm==0 )
      {
         os << "00:" ;
      }
      else if( mm<10 )
      {
         os << "0" << mm << ":" ;
      }
      else
      {
         os << mm << ":" ;
      }
      double const ss = ( (int) (10.*(t-3600.*hh-60.*mm)) )/10. ;
      if( ss==0 )
      {
         os << "00" ;
      }
      else if( ss<10 )
      {
         os << "0" << ss ;
      }
      else
      {
         os << ss ;
      }
      os << " h:m:s" ;
   }
}

//-------------------------------------------------------------------------
bool
MAC_Timer:: invariant( void ) const
//-------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   MAC_ASSERT( !( RUNNING && CURRENT_TIME==MAC::bad_double() ) ) ;
   MAC_ASSERT( !( RUNNING && CURRENT_ELAPSED_TIME==MAC::bad_double() ) ) ;
   return( true ) ;
}
