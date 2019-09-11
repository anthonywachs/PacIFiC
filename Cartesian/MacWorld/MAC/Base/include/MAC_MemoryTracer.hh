#ifndef MAC_MEMORY_TRACER_HH
#define MAC_MEMORY_TRACER_HH

#include <MAC_Object.hh>

#include <size_t_vector.hh>
#include <stringVector.hh>

#include <string> // size_t

/*
Server devoted to memory trace.

PUBLISHED
*/

class MAC_MemoryTracer : public MAC_Object
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_MemoryTracer* object( void ) ;
      
   //-- Enable/disable trace(5000.)

      void enable_memory_trace( void ) ;
      void disable_memory_trace( void ) ;
      
      bool memory_trace_enabled( void ) const ;
      
   //-- Memory used(5001.)
      
      static size_t used_memory( void ) ;
      static void display_memory( std::ostream& os, size_t memory ) ;

   //-- Memory trace(5002.)

      void start_event( std::string const& label ) ;
      void stop_event( void ) ;

      void trace( std::string const& a_message ) ;
      std::string const& indent( void ) const ;
      std::string const& message( std::string const& label ) const ;
      
   protected: //------------------------------------------------------------

   private: //---------------------------------------------------------------

      MAC_MemoryTracer( void ) ;
      MAC_MemoryTracer( MAC_MemoryTracer const& other ) ;
      MAC_MemoryTracer& operator=( MAC_MemoryTracer const& other ) ;

     ~MAC_MemoryTracer( void ) ;
      
      void init_trace_file( void ) ;
      
  //--- Class attributes

      bool MEMORY_TRACE ;
      std::string OFILENAME ;
      
      std::string INDENT ;
      stringVector EVENTS ;
      size_t_vector MEM0 ;
      size_t_vector OBJ0 ;
} ;

#endif
