#ifndef FV_MATLAB_POSTPROCESSING_WRITER_HH
#define FV_MATLAB_POSTPROCESSING_WRITER_HH


#include <FV_PostProcessingWriter.hh>
#include <MAC_Object.hh>
#include <FV_DiscreteField.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_IntVector.hh>

#include <string>
#include <vector>
#include <stringVector.hh>
#include <list>
using std::list ;
using std::string ;

class size_t_vector ;
class LA_Vector ;
class MAC_Communicator ;
class MAC_Module ;
class MAC_ModuleExplorer ;
class MAC_Vector ;
class FV_Mesh ;
class FV_DiscreteField ;
class FV_TimeIterator ;


/** @brief The Class MatlabPostProcessingWriter.

Write mesh and fields data in the Matlab format

@author A. Hammouti - Particulate flow project 2011-2013 */
class FV_MatlabPostProcessingWriter  : public FV_PostProcessingWriter
{
   public: //-----------------------------------------------------------
 
      void write_cycle( FV_TimeIterator const* t_it,
      			size_t cycle_number ) ;
			      

   //-- Data clearing

      void clearResultFiles( void ) ;
      
      size_t getPreviousCycleNumber( void ) ;
      
      void readTimeFile( FV_TimeIterator const* t_it,
      			size_t& cycle_number ) ; 	


   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Constructor without argument */   
      FV_MatlabPostProcessingWriter( void ) ;

      /** @brief Constructor with argument
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param com MPI communicator
      @param a_fields list of FV_DiscreteField to be output
      @param a_primary_mesh primary grid
      @param a_binary true if files written in binary mode */
      FV_MatlabPostProcessingWriter( MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp ,
                MAC_Communicator const* com,
       		list< FV_DiscreteField const* > a_fields,
       		FV_Mesh const* a_primary_mesh,
        	bool a_binary );  
      
      /** @brief Destructor */      
      ~FV_MatlabPostProcessingWriter( void ) ;
    
      /** @brief Copy constructor */  
      FV_MatlabPostProcessingWriter(
          FV_MatlabPostProcessingWriter const& other ) ;
      
      /** @brief Operator == 
      @param other the right hand side */   
      FV_MatlabPostProcessingWriter& operator
     	    =( FV_MatlabPostProcessingWriter const& other ) ;
     
      /**
        @brief Constructor called by
            FV_MatlabPostProcessingWriter::create_replica
        @param a_owner the MAC-based object
      */
      FV_MatlabPostProcessingWriter( MAC_Object* a_owner ) ;
      //@}


   //-- Plug in

      /** @name Plug in */
      //@{
      /** @brief Registration of an instance
      @param a_name instance name */
      FV_MatlabPostProcessingWriter( std::string const& a_name ) ;

      /** @brief Create and initialize an instance of FV_PostProcessingWriter
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param com MPI communicator
      @param a_fields list of FV_DiscreteField to be output
      @param a_primary_mesh primary grid
      @param a_binary true if files written in binary mode */
      FV_PostProcessingWriter* create_replica (MAC_Object* a_owner, 
	       	 MAC_ModuleExplorer const* exp,
       		 MAC_Communicator const* com,
       		 list< FV_DiscreteField const* > a_fields,
       		 FV_Mesh const* a_primary_mesh,
        	 bool a_binary ) const; 
      //@}

       
   //-- Class attributes

      static FV_MatlabPostProcessingWriter const* PROTOTYPE ;

   // -- Writing process main steps
      
      void write_time_file( FV_TimeIterator const* t_it,
                           std::string const& vtr_filename ) ;
			         
      void output_domain( void ) ;
      
      void output_proc_domain( void ) ;

      void build_module( MAC_Module* matlab,
                      bool parallel )  ;

      std::string output_file_name( size_t nb,
                                    bool parallel,
                                    size_t rank ) ;

   //-- Data writing

      void write_module( MAC_ModuleExplorer* matlab,
                      std::ofstream& file,
                      size_t level,
                      bool parallel ) ;

      void start_output( size_t size, size_t number ) ;

      void header_output( void ) ;

      void write_double( std::ofstream& file, double val ) ;

      void write_int( std::ofstream& file, int val ) ;

      size_t store_int( int val ) ;
      
      size_t store_double( double val ) ;

      void check_allocated( size_t size ) ;

   //-- Attributes

      MAC_ModuleExplorer const* EXP ;
      MAC_Communicator const* COM ;

      list< FV_DiscreteField const* > FVFIELDS ;
      FV_Mesh const* PRIMARY_GRID ;
      
      std::string RES_DIRECTORY ;
      std::string BASE_FILENAME ;
      std::string TIME_FILENAME ;
      std::string TIME_FILENAME_PVD; //JL Pierson 2016
      stringVector TIME_STRINGS ;
      size_t CYCLE_NUMBER ;
      std::string OUTPUT_DOMAIN_FILENAME ;
      std::string OUTPUT_PROC_GRID_FILENAME ;
      bool BINARY ;
      char * BUFFER ;
      size_t ALLOCATED ;
      size_t OFFSET ;
      size_t CURRENT_LENGTH ;
      size_t NB_SPACE_DIMENSION ;
      size_t_vector WHOLE_EXTENT;
      intVector EXTENT;
      vector< intVector >* EXTENT_AllProc; 
      doubleVector BORDER; 
      size_t NXPROCS ;
      size_t NYPROCS ;
      size_t NZPROCS ;
      size_t NXPROCS_LOC ;
      size_t NYPROCS_LOC ;
      size_t NZPROCS_LOC ;
      size_t MY_COL ;
      size_t MY_ROW ;
      size_t MY_PLN ;
  
      size_t_vector const* p_GLOBAL_MAX_INDEX ;
      size_t_vector const* p_GLOBAL_MIN_INDEX ;      
      size_t_vector const* p_LOCAL_MAX_INDEX ;
      size_t_vector const* p_LOCAL_MIN_INDEX ;


} ;

#endif
