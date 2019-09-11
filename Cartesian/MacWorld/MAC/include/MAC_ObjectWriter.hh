#ifndef MAC_OBJECT_WRITER_HH
#define MAC_OBJECT_WRITER_HH

#include <MAC_Object.hh>

#include <stack>
#include <string>

class MAC_Data ;
class MAC_ModuleExplorer ;

/*
Servers used to store objects so that they can be retrieved with
associated `MAC_ObjectReader::' instances.

Objects are stored in files according to some options
specified in the Hierarchical Data Structure provided 
at creation of `self'. That data structure is briefly described below.

The entry of keyword "output_format" defines the format of the
saving files. There are two possibilities:
   - "text": human readable but not exact (truncated values)
   - "hybrid": parts of the data remain readable but 
               double or integer values are stored in binary format 
               (file with "bin" extension)
                
Several saving strategies are available :

   - all the saved cycles are stored in the same file

     example :
     
        MODULE MAC_ObjectWriter
           type = "all_cycles_in_one_file"
           file_name = join( getcwd(), "saving.mac" )
           output_format = "hybrid"
        END MODULE MAC_ObjectWriter

        A text file named "saving.mac" is created to store all the cycles
        (a companion binary file named "saving.mac.bin" is also created 
        to store the double and integer values).

   - each saved cycle is stored in a separate file (one cycle per file)

     example :
     
        MODULE MAC_ObjectWriter
           type = "cycles_in_separate_files"
           files_basename = join( getcwd(), "saving" )
           output_format = "hybrid"
        END MODULE MAC_ObjectWriter

        A sequence of text files named "saving.00001.mac", 
        "saving.00002.mac",... is created to store respectively 
        the first cycle, the second cycle,...
        (a sequence of companion binary files named "saving.00001.mac.bin", 
        "saving.00002.mac.bin",... is also created to stored
         the double and integer values).

   - only the last two cycles are stored

     example :
     
        MODULE MAC_ObjectWriter
           type = "last_two_cycles"
           file_name_0 = join( getcwd(), "saving_0.mac" )
           file_name_1 = join( getcwd(), "saving_1.mac" )
           output_format = "hybrid"
        END MODULE MAC_ObjectWriter

        The text files "saving_0.mac" and "saving_1.mac" are created to store
        the last two cycles ; the time of last modification of these files 
        identifies that of the more recent saving.
        (the companion binary files named "saving_0.mac.bin" and 
        "saving_1.mac.bin" are also created to store with the double and 
        integer values).


See `MAC_ObjectReader::' and `MAC_ApplicationRestorer::' for restoring
objects stored with `MAC_ObjectWriter::' objects.

PUBLISHED
*/

class MAC_ObjectWriter : public MAC_Object
{
   public: //----------------------------------------------------------------

   //-- Instance delivery and initialization
      
      static MAC_ObjectWriter* create( MAC_Object* a_owner,
                                       MAC_ModuleExplorer const* exp,
                                       MAC_ModuleExplorer const* header_exp ) ;

   //-- Cycles(1.)

      // Start a new cycle.
      void start_cycle( double const& time ) ;

      // Terminate the current cycle. `::finalize_object' must have been
      // called as many times as `::start_new_object'. If not, a fatal error 
      // is raised.
      void terminate_cycle( void ) ;

      // Is there a cycle that is started and not terminated ?
      bool has_an_opened_cycle( void ) const ;

      // cycle number
      size_t cycle_number( void ) const ;

   //-- Object storing(2.)

      // Notify that the storage of a new object is starting so that all 
      // subsequent calls to `::add_entry' are relative that object,
      // until `::finalize_object' or `::start_new_object' are called.
      void start_new_object( std::string const& class_name ) ;

      // nonzero number associated to the object being currently stored if any,
      // 0 otherwize
      size_t current_object_number( void ) const ;

      void add_entry( std::string const& keyword, MAC_Data* data ) ;

      // Notify that the storage of the current object is completed.
      void finalize_object( void ) ;
      
      // Current file name
      std::string get_current_file_name( void ) const;
      
      // Force writing at last time
      bool force_write_at_last_time( void ) const;
      
      // Swap restart file names in case of last_two_cycles and  
      // OFILE_NAME0 is equal to inputfilename
      void swap_restart_file_names( std::string const& inputfilename );

   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

      // all_cycles :      save all the cycles in one file
      // per_one_cycle :   save all the cycles, but one per file
      // last_two_cycles : save only the two last cycles
      enum MAC_ObjectWriterType 
      {
         all_cycles,
         per_one_cycle,
         last_two_cycles
      } ;
         

      MAC_ObjectWriter( void ) ;
     ~MAC_ObjectWriter( void ) ;
      MAC_ObjectWriter( MAC_ObjectWriter const& other ) ;
      MAC_ObjectWriter& operator=( MAC_ObjectWriter const& other ) ;

      MAC_ObjectWriter( MAC_Object* a_owner,
                        MAC_ObjectWriterType const writer_type,
                        MAC_ModuleExplorer const* exp,
                        MAC_ModuleExplorer const* header_exp ) ;

   //-- Internals

      void initialize_saving_file( void ) ;
      void set_file_name( void ) ;
      void write_header( void ) const ;
      void write_communicator( void ) const ;
      void write_rftable_file( double const& time );

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      MAC_ObjectWriterType const TYPE ;

      std::string const OFILE_FORMAT ;
      MAC_ModuleExplorer const* const HEADER_EXP ;

      // File names :
      std::string OFILE_NAME ;
      std::string OFILE_NAME0 ;
      std::string OFILE_NAME1 ;
      std::string RFTABLE ;
      std::string HEADERFILE ;
      
      bool B_WRITE_LAST_ITER;
      
      size_t iCYCLE ;
      size_t NB_OBJECTS ;
      std::stack< MAC_Module* > MODS ;
} ;

#endif
