#ifndef _RAWDATAPOSTPROCESSINGWRITER_HH_
#define _RAWDATAPOSTPROCESSINGWRITER_HH_

#include "PostProcessingWriter.hh"
#include <iostream>
#include <fstream>
using std::ofstream;


/** @brief The class RawDataPostProcessingWriter

    Writes particle data as arrays of raw data in files

    @author M.BERNARD - IFPEN - 2012 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
//=============================================================================
class RawDataPostProcessingWriter : public PostProcessingWriter
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
    RawDataPostProcessingWriter( DOMNode* dn, int const& rank_, 
    	int const& nbranks_, bool const& verbose = true );

    /** @brief Destructor */
    virtual ~RawDataPostProcessingWriter();
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
    virtual void PostProcessing_start( double const& time,
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
    virtual void PostProcessing( double const& time,
  	double const& dt,
  	list<Particle*> const* particles,
	list<Particle*> const* pwait,
	list<Particle*> const* periodic_clones,		
	vector<Particle*> const* referenceParticles,
	Obstacle* obstacle,
	LinkedCell const* LC );

    /** @brief Finalizes writing data */
    virtual void PostProcessing_end();
  
    /** @brief Gets the post-processing writer type */
    virtual string getPostProcessingWriterType() const;
    //@}
    
    
  private:
    /** @name Methods */
    //@{
    /** @brief Writes data in parallel mode at one physical time
    @param time physical time
    @param nb_total_part total number of particles
    @param data_Global vector containing particle data  */
    void one_output_MPI( double const& time, size_t& nb_total_part,
  	vector< vector<double> > const* data_Global );
	
    /** @brief Writes data in serial mode at one physical time
    @param time physical time
    @param nb_total_part total number of particles     
    @param particles active particles */
    void one_output_Standard( double const& time, size_t& nb_total_part,
  	list<Particle*> const* particles );
	
    /** @brief Delete all result files */
    void clearResultFiles() const ; 	  
  
    /** @brief Creates output files and open streams
    @param mode File opening mode (here : ios::app) */
    void prepareResultFiles( ios_base::openmode mode ) ;
    //@}
  

  private:
    /**@name Parameters */
    //@{  
    ofstream m_gc_coordinates_x; /**< center of mass x-coordinate output 
    	stream */
    ofstream m_gc_coordinates_y; /**< center of mass y-coordinate output 
    	stream */ 
    ofstream m_gc_coordinates_z; /**< center of mass z-coordinate output 
    	stream */
    ofstream m_translational_velocity_x; /**< x-translational velocity output 
    	stream */ 
    ofstream m_translational_velocity_y; /**< y-translational velocity output 
    	stream */    
    ofstream m_translational_velocity_z; /**< z-translational velocity output 
    	stream */
    ofstream m_angular_velocity_x; /**< x-angular velocity output stream */ 
    ofstream m_angular_velocity_y; /**< y-angular velocity output stream */    
    ofstream m_angular_velocity_z; /**< z-angular velocity output stream */
    ofstream m_coordination_number; /**< coordination number output stream */
    ofstream m_particle_class; /**< particle class output stream */
    string m_filerootname; /**< files root name */
    bool m_binary; /**< whether to write data in binary */
    int m_ndigits; /**< number of digits after the decimal in the scientific
    	format in text mode writing */ 
       


  protected:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    RawDataPostProcessingWriter();  
    //@}  
};

#endif
  
