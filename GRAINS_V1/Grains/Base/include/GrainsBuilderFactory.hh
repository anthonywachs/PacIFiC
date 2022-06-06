#ifndef _GRAINSBUILDERFACTORY_HH_
#define _GRAINSBUILDERFACTORY_HH_

#include "Grains.hh"
#include "GrainsCoupledWithFluid.hh"
#include "ReaderXML.hh"


/** @brief Space dimension of the simulation */
enum EAPPLI 
{
  DIM_3,
  DIM_2,
  UNDEFINED
};


/** @brief The class GrainsBuilderFactory.

    Creates the appropriate Grains application depending on options.
    
    @author A.WACHS - Institut Francais du Petrole - 2009 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
//=============================================================================
class GrainsBuilderFactory
{
  public:
    /**@name Static methods */
    //@{
    /** @brief Creates and returns a standard Grains application
    @param root XML root ("<Grains3D>" or "<Graind2D>") */
    static Grains* create( DOMElement* root );
    
    /** @brief Creates and returns a GrainsCoupledWithFluid application
    @param root XML root ("<Grains3D>" or "<Graind2D>")
    @param fluid_density_ fluid density */
    static GrainsCoupledWithFluid* createCoupledWithFluid( DOMElement* root,
    	double fluid_density_ );      

    /** @brief Returns the space dimension of the simulation */
    static EAPPLI getContext();
    
    /** @brief Adds the path to the dtd files using the GRAINS_HOME variable to
    a copy of the input file. Returns the name of this copy
    @param filename input file name 
    @param rank rank of process 
    @param nprocs total number of processes */
    static string init( string const& filename, int const& rank, 
  	int const& nprocs );  
    //@}


  private:
    /** @name Parameters */
    //@{
    static EAPPLI m_context; /**< Space dimension of the simulation */
    //@}


    /**@name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    GrainsBuilderFactory() {}

    /** @brief Destructor (forbidden) */
    ~GrainsBuilderFactory() {}
    //@}
};

#endif

