#ifndef MAC_OBJECT_READER_HH
#define MAC_OBJECT_READER_HH

#include <MAC_Object.hh>
#include <fstream>
#include <stack>
#include <string>
#include <iosfwd>

class MAC_Data ;
class MAC_Module ;
class MAC_ModuleExplorer ;
class MAC_ModuleIterator ;

/*
Servers used to retrieve objects stored with
associated `MAC_ObjectWriter::' instances.

PUBLISHED
*/

class MAC_ObjectReader : public MAC_Object
{
   public: //----------------------------------------------------------------

   //-- Instance delivery and initialization
      
      static MAC_ObjectReader* create( MAC_Object* a_owner,
                                       MAC_ModuleExplorer const* exp ) ;

   //-- Cycles

      MAC_Module* header_module( void ) const ;

      size_t nb_cycles( void ) const ;

      void seek_cycle( size_t cycle_number ) ;
      
      void print( std::ostream& os, size_t indent_width ) const;

      void close_cycle( void ) ;

      bool positioned_in_a_valid_cycle( void ) const ;

   //-- Object retrieval

      void start_object_retrieval( std::string const& class_name ) ;

      size_t current_object_number( void ) const ;

      bool has_entry( std::string const& keyword ) const ;
      MAC_Data const* data_of_entry( std::string const& keyword ) const ;

      void end_object_retrieval( void ) ; 
      
      std::string next_object_name_in_current_module( void ) const;
      
      std::string next_object_class_in_current_module( void ) const;      
      
   //-- Input file
   
      std::string input_file_name( void ) const; 
      
   //-- Restart time from "RFTable.txt"
   
      double get_initial_time( void ) const;        
    
   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

      MAC_ObjectReader( void ) ;
     ~MAC_ObjectReader( void ) ;
      MAC_ObjectReader( MAC_ObjectReader const& other ) ;
      MAC_ObjectReader& operator=( MAC_ObjectReader const& other ) ;

      MAC_ObjectReader( MAC_Object* a_owner,
                        MAC_ModuleExplorer const* exp ) ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      MAC_Module const* ROOT_MOD ;
      MAC_Module* HEADER_MOD ;      
      std::string IFILE_NAME ;
      std::string HEADERFILE ; 
      size_t NB_CYCLES ;
      size_t LAST_CYCLE ;
      int iOBJECT ;
      std::stack< MAC_Module const* > MODS ;
      std::stack< MAC_ModuleIterator* > MOD_ITS ;
      double initial_time ; 
} ;

#endif
