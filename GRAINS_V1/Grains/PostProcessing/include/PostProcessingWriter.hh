#ifndef _POSTPROCESSINGWRITER_HH_
#define _POSTPROCESSINGWRITER_HH_

#include "Basic.hh"
#include <list>
#include <string>
using std::list;
using std::string;
#include "ReaderXML.hh"
#include "LinkedCell.hh"

class Particle;
class Obstacle;
class Component;


/** @brief The class PostProcessingWriter.

    Writes results in files for post-processing by an external software.

    @author A.WACHS - Institut Francais du Petrole - 2010 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
//=============================================================================
class PostProcessingWriter
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with rank and number of processes as input parameters
    @param rank_ process rank 
    @param nbranks_ number of processes */
    PostProcessingWriter( int const& rank_, int const& nbranks_ );

    /** @brief Constructor with XML node, rank and number of processes 
    as input parameters
    @param dn XML node
    @param rank_ process rank 
    @param nbranks_ number of processes */
    PostProcessingWriter( DOMNode* dn, int const& rank_, int const& nbranks_);

    /** @brief Destructor */
    virtual ~PostProcessingWriter();
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
	vector<Window> const& insert_windows ) = 0;

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
	LinkedCell const* LC ) = 0;

    /** @brief Finalizes writing data */
    virtual void PostProcessing_end() = 0;
  
    /** @brief Sets the initial cycle number
    @param cycle0 initial cycle number */
    virtual void setInitialCycleNumber( int const& cycle0 );
  
    /** @brief Gets the post-processing writer type */
    virtual string getPostProcessingWriterType() const = 0;
  
    /** @brief Writes components involved in a displacement or a contact error
    @param filename file root name
    @param errcomposants list of the 2 components invovled */
    virtual void writeErreurComponentsPostProcessing( string const& filename,
  	list<Component*> const& errcomposants );  
    //@}

  
    /** @name Static Methods */
    //@{
    /** @brief Allocates post-processing windows and initializes all processes
    to true
    @param nbRank number of processes */
    static void allocate_PostProcessingWindow( int const& nbRank );
  
    /** @brief Sets whether process of rank rank_ writes data
    @param rank_ process rank
    @param bPPWindow whether this process writes data */
    static void set_PostProcessingWindow( int const& rank_,
  	bool const& bPPWindow );

    /** @brief Returns the vector of processes invovled in writing data */
    static vector<bool> get_PostProcessingWindow();
    //@}
  

  protected:
    /**@name Parameters */
    //@{  
    int m_rank; /**< process rank */
    int m_nprocs; /**< number of processes */
    //@}

    /**@name Static Parameters */
    //@{  
    static vector<bool> m_bPPWindow;
    //@}

    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    PostProcessingWriter();
    //@}
  
};

#endif
  
