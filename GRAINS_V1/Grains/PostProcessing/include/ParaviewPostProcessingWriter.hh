#ifndef _PARAVIEWPOSTPROCESSINGWRITER_HH_
#define _PARAVIEWPOSTPROCESSINGWRITER_HH_

#include "PostProcessingWriter.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using std::ofstream;
using std::ostringstream;
using std::vector;


/** @brief The class ParaviewPostProcessingWriter.

    Writes data in files for post-processing with Paraview.

    @author A.WACHS - Institut Francais du Petrole - 2010 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
//=============================================================================
class ParaviewPostProcessingWriter : public PostProcessingWriter
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with XML node, rank and number of processes 
    as input parameters
    @param dn XML node
    @param rank_ process rank 
    @param nbranks_ number of processes 
    @param verbose outputs writer features if true */
    ParaviewPostProcessingWriter( DOMNode* dn, int const& rank_, 
    	int const& nbranks_, bool const& verbose = true );
  
    /** @brief Constructor with input parameters 
    @param rank_ process rank 
    @param nbranks_ number of processes 
    @param name_ files name
    @param root_ file root name
    @param isBinary whether to write in binary
    @param verbose outputs writer features if true */
    ParaviewPostProcessingWriter( int const& rank_,
  	int const& nbranks_,
  	const string &name_,
	const string &root_,
  	const bool &isBinary, 
	bool const& verbose = true );  

    /** @brief Destructor */
    virtual ~ParaviewPostProcessingWriter();
    //@}


    /** @name Methods */
    //@{
    /** @brief Initializes the post-processing writer
    @param time physical time
    @param dt time step magnitude
    @param particles active particles
    @param pwait inactive particles
    @param periodic_clones periodic particles
    @param referenceParticles reference particles
    @param obstacle obstacles 
    @param LC linked-cell grid
    @param insert_windows insertion windows */
    void PostProcessing_start( double const& time,
  	double const& dt,
  	list<Particle*> const* particles,
	list<Particle*> const* pwait,
	list<Particle*> const* periodic_clones,	
	vector<Particle*> const* referenceParticles,
	Obstacle* obstacle,
	LinkedCell const* LC,
	vector<Window> const& insert_windows );

    /** @brief Writes data
    @param time physical time
    @param dt time step magnitude
    @param particles active particles
    @param pwait inactive particles
    @param periodic_clones periodic particles
    @param referenceParticles reference particles
    @param obstacle obstacles 
    @param LC linked-cell grid */
    void PostProcessing( double const& time,
  	double const& dt,
  	list<Particle*> const* particles,
	list<Particle*> const* pwait,
	list<Particle*> const* periodic_clones,		
	vector<Particle*> const* referenceParticles,
	Obstacle* obstacle,
	LinkedCell const* LC );

    /** @brief Finalizes writing data */
    void PostProcessing_end();
    
    /** @brief Gets the post-processing writer type */
    string getPostProcessingWriterType() const;
	
    /** @brief Sets the initial cycle number
    @param cycle0 initial cycle number */
    void setInitialCycleNumber( const int& cycle0 );

    /** @brief Writes components involved in a displacement or a contact error
    @param filename file root name
    @param errcomposants list of the 2 components invovled */
    void writeErreurComponentsPostProcessing( string const& filename,
  	list<Component*> const& errcomposants );  

    /** @brief Writes particles data
    @param particles active particles
    @param partFilename output file name
    @param forceForAllTag writes particle data regardless of the particle tag 
    @param processwrites whether this process writes data */
    void writeParticlesPostProcessing_Paraview(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag = false,
	bool const& processwrites = true );
    //@}
  
  
    /** @name Methods to write in binary (undocumented) */
    //@{
    void start_output_binary( int size, int number ) ;
    void write_double_binary( double val ) ;
    void write_int_binary( int val ) ;
    int store_int_binary( int val ) ;
    void flush_binary( std::ofstream& file, string const& calling ) ;
    void check_allocated_binary( int size ) ;
    void compress_segment_binary( int seg, string const& calling ) ;
    //@}  


  private:
    string m_ParaviewFilename_dir; /**< output directory name */
    string m_ParaviewFilename; /**< output file name */  
    ostringstream m_Paraview_saveObstacles_pvd; /**< obstacles output stream */	
    vector<ostringstream*> m_Paraview_saveParticles_pvd; /**< particles output
  	stream */
    ostringstream m_Paraview_savePeriodicCloneParticles_pvd; /**< periodic 
    	particles output stream */	
    ostringstream m_Paraview_saveParticleVelocityVectors_pvd; /**< particle 
    	translational velocity output stream */
    ostringstream m_Paraview_saveContactForceVectors_pvd; /**< contact force 
    	output stream */		
    ostringstream m_Paraview_saveContactForceChains_pvd; /**< force chain 
    	network output stream */		
    int m_ParaviewCycleNumber; /**< cycle number */
    bool m_binary; /**< whether to write data in binary */
    bool m_postProcessObstacle; /**< whether to write obstacles data */
    bool m_initialCycleNumber_forced; /**< whether to force the initial cycle
    	number to be equal to a given value */
    bool m_network; /**< whether to write force chain network data */
    bool m_mpiio_singlefile; /**< whether to use MPI I/O writing routines to
    	output data in a single file */
    bool m_pertype; /**< whether to write particle files per type */
    char * BUFFER ;
    int ALLOCATED ;
    int OFFSET ;
    int CURRENT_LENGTH ;
    list<string> empty_string_list; 

  
    /** @name Methods */
    //@{  
    /** @brief Writes data at one physical time
    @param time physical time
    @param dt time step magnitude
    @param particles active particles
    @param periodic_clones periodic particles
    @param referenceParticles reference particles
    @param obstacle obstacles 
    @param LC linked-cell grid */
    void one_output( double const& time,
  	double const& dt,
	list<Particle*> const* particles,
	list<Particle*> const* periodic_clones,
	vector<Particle*> const* referenceParticles,
  	Obstacle* obstacle,
	LinkedCell const* LC ); 

    /** @brief Writes obstacles data
    @param allObstacles list of simple obstcales 
    @param obsFilename output file name */
    void writeObstaclesPostProcessing_Paraview(
  	list<SimpleObstacle*> const& allObstacles,
  	string const& obsFilename );
	
    /** @brief Writes spherical particles data in a vector form containing the
    center of mass coordinates of each particle
    @param particles active particles
    @param partFilename output file name
    @param forceForAllTag writes particle data regardless of the particle tag 
    @param processwrites whether this process writes data */
    void writeSpheresPostProcessing_Paraview(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag = false,
	bool const& processwrites = true );	
	
    /** @brief Writes particle translational and angular velocity vectors
    @param particles active particles
    @param partFilename output file name 
    @param processwrites whether this process writes data */
    void writeParticleVelocityVectorsPostProcessing_Paraview(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& processwrites = true );

    /** @brief Writes contact force vectors 
    @param particles active particles   
    @param LC linked-cell gruid
    @param partFilename output file name
    @param time physical time 
    @param processwrites whether this process writes data */
    void writeContactForceVectorsPostProcessing_Paraview(
  	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& partFilename, double const& time,
	bool const& processwrites = true );
				
    /** @brief Writes force chain network data
    @param particles active particles   
    @param LC linked-cell gruid
    @param filename output file name
    @param time physical time */
    void writeContactForceChains_Paraview(
  	list<Particle*> const* particles,
  	LinkedCell const* LC, const string &filename, double const& time );

    /** @brief Writes a pvtp file 
    @param filename output file name 
    @param pointVector list of vectorial field names at data points
    @param pointScalar list of scalar field names at data points
    @param cellScalar list of scalar field names in cells */
    void writePVTP_Paraview( const string &filename,
  	list<string> const* pointVector,
	list<string> const* pointScalar,
	list<string> const* cellScalar );

    /** @brief Writes a pvtu file 
    @param filename output file name 
    @param pointVector list of vectorial field names at data points
    @param pointScalar list of scalar field names at data points
    @param cellScalar list of scalar field names in cells */
    void writePVTU_Paraview( const string &filename,
  	list<string> const* pointVector,
	list<string> const* pointScalar,
	list<string> const* cellScalar );
	
    /** @brief Writes linked-cell grid data
    @param LC linked-cell grid
    @param partFilename output file name 
    @param local use local coordinates if true, global if false */
    void writeLinkedCellPostProcessing_Paraview(
  	LinkedCell const* LC, string const& partFilename,
	bool const& local );

    /** @brief Writes insertion windows data
    @param insert_windows insertion windows
    @param partFilename output file name */
    void writeInsertionPostProcessing_Paraview(
   	vector<Window> const& insert_windows,
  	string const& partFilename );
	
    /** @brief Writes periodic boundary data
    @param LC linked-cell grid
    @param partFilename output file name  */
    void writePeriodicBoundaryPostProcessing_Paraview(
  	LinkedCell const* LC, string const& partFilename );	
		
    /** @brief Updates obstacles indicator
    @param time physical time
    @param dt time step magnitude
    @param obstacle obstacles */
    void updateObstaclesIndicator( double const& time,
  	double const& dt, Obstacle* obstacle );
  
    /** @brief Delete all result files */
    void clearResultFiles() const;
  
    /** @brief Reads a pvd file and transfers its content to an output stream
    @param filename input file name
    @param ossflux output stream corresponding to the input file */
    void readPVDFile( string const& filename, ostringstream& ossflux );
  
    /** @brief Gets the last cycle number in case of reload */
    int getPreviousCycleNumber() const;
  
    /** @brief Writes a long ostringstream character by character (in case the
    classical fileOUT << oss.str() crashes)
    @param fileOUT output file name
    @param oss the output stream to be written */  
    void writeBigOSS( ofstream& fileOUT, ostringstream const& oss );
	
    /** @brief Writes particles data in a single MPI file in text mode with 
    MPI I/O routines
    @param particles active particles
    @param partFilename output file name
    @param forceForAllTag writes particle data regardless of the particle tag 
    @param processwrites whether this process writes data */
    void writeParticlesPostProcessing_Paraview_MPIIO_text(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag = false,
	bool const& processwrites = true ); 
	
    /** @brief Writes particles data in a single MPI file in binary mode with 
    MPI I/O routines
    @param particles active particles
    @param partFilename output file name
    @param forceForAllTag writes particle data regardless of the particle tag 
    @param processwrites whether this process writes data */
    void writeParticlesPostProcessing_Paraview_MPIIO_binary(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag = false,
	bool const& processwrites = true ); 
	
    /** @brief Writes spherical particles data in a vector form containing the
    center of mass coordinates of each particle in a single MPI file in text 
    mode with MPI I/O routines
    @param particles active particles
    @param partFilename output file name
    @param forceForAllTag writes particle data regardless of the particle tag 
    @param processwrites whether this process writes data */
    void writeSpheresPostProcessing_Paraview_MPIIO_text(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag = false,
	bool const& processwrites = true ); 
	
    /** @brief Writes spherical particles data in a vector form containing the
    center of mass coordinates of each particle in a single MPI file in binary 
    mode with MPI I/O routines
    @param particles active particles
    @param partFilename output file name
    @param forceForAllTag writes particle data regardless of the particle tag 
    @param processwrites whether this process writes data */
    void writeSpheresPostProcessing_Paraview_MPIIO_binary(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag = false,
	bool const& processwrites = true );
	
    /** @brief Writes particle translational and angular velocity vectors in a 
    single MPI file in text mode with MPI I/O routines
    @param particles active particles
    @param partFilename output file name 
    @param processwrites whether this process writes data */
    void writeParticleVelocityVectorsPostProcessing_Paraview_MPIIO_text(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag = false,
	bool const& processwrites = true );
	
    /** @brief Writes particle translational and angular velocity vectors in a 
    single MPI file in binary mode with MPI I/O routines
    @param particles active particles
    @param partFilename output file name 
    @param processwrites whether this process writes data */
    void writeParticleVelocityVectorsPostProcessing_Paraview_MPIIO_binary(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag = false,
	bool const& processwrites = true );
	
    /** @brief Writes contact force vectors in a single MPI file in text mode 
    with MPI I/O routines 
    @param particles active particles   
    @param LC linked-cell gruid
    @param partFilename output file name
    @param time physical time 
    @param processwrites whether this process writes data */
    void writeContactForceVectorsPostProcessing_Paraview_MPIIO_text(
  	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& partFilename, double const& time,
	bool const& processwrites = true );
	
    /** @brief Writes contact force vectors in a single MPI file in binary mode 
    with MPI I/O routines 
    @param particles active particles   
    @param LC linked-cell gruid
    @param partFilename output file name
    @param time physical time 
    @param processwrites whether this process writes data */
    void writeContactForceVectorsPostProcessing_Paraview_MPIIO_binary(
  	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& partFilename, double const& time,
	bool const& processwrites = true );		
    //@}
      

  protected:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    ParaviewPostProcessingWriter();
  //@}  
};

#endif
  
