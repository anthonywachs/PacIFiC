#include <MAC_MemoryTracer.hh>

#include <MAC_assertions.hh>
#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Root.hh>
#include <MAC_System.hh>

#include <iostream>
#include <fstream>
#include <sstream>

struct MAC_MemoryTracer_ERROR
{
   static void n0( std::string const& file_name ) ;
   static void n1( void ) ;
} ;

//----------------------------------------------------------------------
MAC_MemoryTracer*
MAC_MemoryTracer:: object( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MemoryTracer:: object" ) ;

   static MAC_MemoryTracer* result = new MAC_MemoryTracer() ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == MAC_Root::object() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_MemoryTracer:: MAC_MemoryTracer( void )
//----------------------------------------------------------------------
   : MAC_Object( MAC_Root::object() )
   , MEMORY_TRACE( false )
   , OFILENAME( "" )
   , INDENT( "" )
   , EVENTS( 0 )
   , MEM0( 0 )
   , OBJ0( 0 )
{
}

//----------------------------------------------------------------------
MAC_MemoryTracer:: ~MAC_MemoryTracer( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
MAC_MemoryTracer:: enable_memory_trace( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MemoryTracer:: enable_memory_trace" ) ;
   MAC_CHECK_PRE( !memory_trace_enabled() ) ;

   init_trace_file() ;
   MAC::out() << "*** MAC_MemoryTracer: memory trace enabled" << std::endl
              << "    trace_file: " << OFILENAME  << std::endl ;
   
   MEMORY_TRACE = true ;
   
   MAC_CHECK_POST( memory_trace_enabled() ) ;
}

//----------------------------------------------------------------------
void
MAC_MemoryTracer:: disable_memory_trace( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MemoryTracer:: disable_memory_trace" ) ;
   MAC_CHECK_PRE( memory_trace_enabled() ) ;

   MAC::out() << "*** MAC_MemoryTracer: memory trace disabled" << std::endl ;
   
   MEMORY_TRACE = false ;
   
   MAC_CHECK_POST( !memory_trace_enabled() ) ;
}

//----------------------------------------------------------------------
bool
MAC_MemoryTracer:: memory_trace_enabled( void ) const
//----------------------------------------------------------------------
{
   return( MEMORY_TRACE ) ;
}

//----------------------------------------------------------------------
size_t
MAC_MemoryTracer:: used_memory( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MemoryTracer:: used_memory" ) ;
   return( MAC_System::used_memory() ) ;
}

//----------------------------------------------------------------------
void
MAC_MemoryTracer:: display_memory( std::ostream& os, size_t memory )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MemoryTracer:: display_memory" ) ;
//   MAC_CHECK_PRE( os ) ; // Not accepted from gcc-9.x.x
   MAC_ASSERT( os.good() ) ;

   static size_t const mo = 1024*1024 ;
   static size_t const go = 1024*1024*1024 ;

   if( memory > go )
   {
      os << ( (double) memory )/go << " Go" ;
   }
   else if( memory > mo )
   {
      os << ( (double) memory )/mo << " Mo" ;
   }
   else
   {
      os << memory << " octets" ;
   }
}

//----------------------------------------------------------------------
void
MAC_MemoryTracer:: start_event( std::string const& label )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MemoryTracer:: start_event" ) ;
   MAC_CHECK_PRE( !label.empty() ) ;

   if( MEMORY_TRACE )
   {
      std::ofstream file( OFILENAME.c_str(), std::ios::out | std::ios::app ) ;
      if( !file )
      {
         MAC_MemoryTracer_ERROR:: n0( OFILENAME ) ;
      }
      size_t const mem = used_memory() ;
      size_t const nb_objs = MAC_Object::GetNumberOf_MAC_objects() ;
      file << INDENT << "### Start: " << label << " (memory: " ;
      display_memory( file, mem ) ;
      file << ", objects: " << nb_objs << ")" << std::endl ;
      file << INDENT << "#" << std::endl ;
      EVENTS.append( label ) ;
      MEM0.append( mem ) ;
      OBJ0.append( nb_objs ) ;
      INDENT += "   " ;
   }
}

//----------------------------------------------------------------------
void
MAC_MemoryTracer:: stop_event( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MemoryTracer:: stop_event" ) ;
   
   if( MEMORY_TRACE )
   {
      if( INDENT.size() < 3 ) MAC_MemoryTracer_ERROR:: n1() ;
      INDENT.erase( INDENT.length()-3 ) ;
      std::ofstream file( OFILENAME.c_str(), std::ios::out | std::ios::app ) ;
      if( !file )
      {
         MAC_MemoryTracer_ERROR:: n0( OFILENAME ) ;
      }
      size_t const mem = used_memory() ;
      size_t const nb_objs = MAC_Object::GetNumberOf_MAC_objects() ;
      size_t const last = EVENTS.size()-1 ;
      file << INDENT << "#          diff memory: " ;
      if( mem>=MEM0(last) )
      {
         display_memory( file, mem-MEM0(last) ) ;
      }
      else
      {
         file << "-" ;
         display_memory( file, MEM0(last)-mem ) ;
      }
      file << ", diff objects: " ;
      if( nb_objs>=OBJ0(last) )
      {
         file << nb_objs-OBJ0(last) ;
      }
      else
      {
         file << "-" << OBJ0(last)-nb_objs ;
      }
      file << std::endl ;
      file << INDENT << "### Stop:  " ;
      file << EVENTS(last) << " (memory: " ;
      display_memory( file, mem ) ;
      file << ", objects: " << nb_objs << ")" << std::endl ;
      file << std::endl ;
      EVENTS.resize( last ) ;
      MEM0.resize( last ) ;
      OBJ0.resize( last ) ;
   }
}

//----------------------------------------------------------------------
void
MAC_MemoryTracer:: trace( std::string const& a_message )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MemoryTracer:: trace" ) ;
   MAC_CHECK_PRE( !a_message.empty() ) ;

   if( MEMORY_TRACE )
   {
      std::ofstream file( OFILENAME.c_str(), std::ios::out | std::ios::app ) ;
      if( !file )
      {
         MAC_MemoryTracer_ERROR:: n0( OFILENAME ) ;
      }
      file << INDENT << a_message << std::endl ;
      file << std::endl ;
      file.close() ;
   }
}

//----------------------------------------------------------------------
std::string const&
MAC_MemoryTracer:: message( std::string const& label ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MemoryTracer:: message" ) ;
   static std::string result ;
   std::ostringstream msg ;
   msg << INDENT << "# " << label << " (memory_used: " ;
   display_memory( msg, used_memory() ) ;
   msg << ", "
       << MAC_Object::GetNumberOf_MAC_objects() << " objects)" ;
   result = msg.str() ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
MAC_MemoryTracer:: indent( void ) const
//----------------------------------------------------------------------
{
   return( INDENT ) ;
}

//----------------------------------------------------------------------
void
MAC_MemoryTracer:: init_trace_file( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MemoryTracer:: init_trace_file" ) ;

   if( OFILENAME.empty() )
   {
      std::stringstream m ;
      m << "memory" ;
      MAC_Communicator const* com = MAC_Exec::communicator() ;
      if( com->nb_ranks() > 1 )
      {
         m << "#" << com->rank() ;
      }
      m << ".txt" ;
      OFILENAME = m.str() ;
      std::ofstream file( OFILENAME.c_str(),
                          std::ios::out | std::ios::trunc ) ;
      if( !file )
      {
         MAC_MemoryTracer_ERROR:: n0( OFILENAME ) ;
      }
      else
      {
         std::string const s = std::string( 60, '#' ) ;
         file << s << std::endl ;
         file << "#" << std::endl ;
         file << "# MAC_MemoryTracer generated file" << std::endl ;
         file << "#" << std::endl ;
         file << s << std::endl ;
         file << std::endl ;
      }
      file.close() ;
   }
}

//internal--------------------------------------------------------------
void
MAC_MemoryTracer_ERROR:: n0( std::string const& file_name )
//internal--------------------------------------------------------------
{
   MAC_Error::object()->raise_plain(
      "*** MAC_MemoryTracer error:\n"
      "    Unable to open file \""+file_name+"\" for writing" ) ;
}

//internal--------------------------------------------------------------
void
MAC_MemoryTracer_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   MAC_Error::object()->raise_plain(
      "*** MAC_MemoryTracer error:\n"
      "    attempt to decrease indentation below the zero limit" ) ;
}
