#ifndef FV_PARAVIEW_POSTPROCESSING_WRITER_HH
#define FV_PARAVIEW_POSTPROCESSING_WRITER_HH


#include <FV_PostProcessingWriter.hh>
#include <FV_DiscreteField.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_IntVector.hh>


#include <string>
#include <vector>
#include <stringVector.hh>
#include <list>
using std::list ;
using std::string ;

class MAC_Communicator ;
class MAC_Module ;
class MAC_ModuleExplorer ;
class MAC_Vector ;
class FV_Mesh ;
class FV_DiscreteField ;
class FV_TimeIterator ;


/** @brief The Class ParaviewPostProcessingWriter.

Write mesh and fields data in the Paraview vtk format

@author G. Vinay - Particulate flow project 2010-2012 */

class FV_ParaviewPostProcessingWriter  : public FV_PostProcessingWriter
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
      FV_ParaviewPostProcessingWriter( void ) ;

      /** @brief Constructor with argument
      @param a_owner the MAC-based object
      @param exp to read the data file
      @param com MPI communicator
      @param a_fields list of FV_DiscreteField to be output
      @param a_primary_mesh primary grid
      @param a_binary true if files written in binary mode */
      FV_ParaviewPostProcessingWriter( MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp ,
                MAC_Communicator const* com,
       		list< FV_DiscreteField const* > a_fields,
       		FV_Mesh const* a_primary_mesh,
        	bool a_binary );       
      
      /** @brief Destructor */      
      ~FV_ParaviewPostProcessingWriter( void ) ;
    
      /** @brief Copy constructor */  
      FV_ParaviewPostProcessingWriter(
          FV_ParaviewPostProcessingWriter const& other ) ;
      
      /** @brief Operator == 
      @param other the right hand side */   
      FV_ParaviewPostProcessingWriter& operator
     	    =( FV_ParaviewPostProcessingWriter const& other ) ;
     
      /**
        @brief Constructor called by
            FV_ParaviewPostProcessingWriter::create_replica
        @param a_owner the MAC-based object
      */
      FV_ParaviewPostProcessingWriter( MAC_Object* a_owner ) ;
      //@}

   
   //-- Plug in

      /** @name Plug in */
      //@{
      /** @brief Registration of an instance
      @param a_name instance name */
      FV_ParaviewPostProcessingWriter( std::string const& a_name ) ;

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
        	 bool a_binary ) const ; 
       //@}

       
     //-- Class attributes

     static FV_ParaviewPostProcessingWriter const* PROTOTYPE ;

     //-- Writing process main steps
      
      void write_pvd_file( FV_TimeIterator const* t_it,
                           std::string const& vtr_filename ) ;

      void build_vtr( MAC_Module* vtk,
                      bool parallel )  ;
		      
      void write_grid( MAC_Module* base,
                       bool parallel )  ;

      std::string output_file_name( size_t nb,
                                    bool parallel,
                                    size_t rank ) ;

   //-- Data writing

      void write_vtk( MAC_ModuleExplorer* vtk,
                      std::ofstream& file,
                      size_t level,
                      bool parallel ) ;

      void start_output( size_t size, size_t number ) ;

      void write_double( std::ofstream& file, double val ) ;

      void write_int( std::ofstream& file, int val ) ;

      size_t store_int( int val ) ;

      void flush( std::ofstream& file )         ;

      void reserve_double( size_t size ) ;

      void check_allocated( size_t size ) ;

      void compress_segment( size_t seg ) ;
            

   //-- Attributes

      MAC_ModuleExplorer const* EXP ;
      MAC_Communicator const* COM ;

      list< FV_DiscreteField const* > FVFIELDS ;
      FV_Mesh const* PRIMARY_GRID ;
      
      std::string RES_DIRECTORY ;
      std::string BASE_FILENAME ;
      std::string PVD_FILENAME ;
      stringVector PVD_STRINGS ;
      size_t CYCLE_NUMBER ;
      bool BINARY ;
      char * BUFFER ;
      size_t ALLOCATED ;
      size_t OFFSET ;
      size_t CURRENT_LENGTH ;
      size_t NB_SPACE_DIMENSION ;
      intVector WHOLE_EXTENT;
      intVector EXTENT;
      vector< intVector >* EXTENT_AllProc; 
      
      size_t_vector const* p_GLOBAL_MAX_INDEX ;
      size_t_vector const* p_GLOBAL_MIN_INDEX ;      
      size_t_vector const* p_LOCAL_MAX_INDEX ;
      size_t_vector const* p_LOCAL_MIN_INDEX ;
      
} ;

#endif
