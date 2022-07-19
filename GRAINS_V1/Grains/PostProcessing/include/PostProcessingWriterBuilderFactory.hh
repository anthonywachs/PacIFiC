#ifndef _POSTPROCESSINGWRITERBUILDERFACTORY_HH_
#define _POSTPROCESSINGWRITERBUILDERFACTORY_HH_

#include "ReaderXML.hh"
#include <string>
using std::string;
class PostProcessingWriter;


/** @brief The class PostProcessingWriterBuilderFactory.

    Creates the appropriate post-processing writer depending on options.

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class PostProcessingWriterBuilderFactory
{
  public:
    /** @name Methods Static */
    //@{
    /** @brief Creates a post-processing writer from an XML node, the process
    rank and the number of processus
    @param nPPW XMl noeud
    @param rank_ process rank 
    @param nbranks_ number of processes */
    static PostProcessingWriter* create( DOMNode* nPPW,
  	int const& rank_, int const& nbranks_ );
    //@}


  private:
    /**@name Constructors & Destructor */
    //@{
    /** @brief Default constructor (forbidden) */
    PostProcessingWriterBuilderFactory() {}

    /** @brief Destructor (forbidden) */
    ~PostProcessingWriterBuilderFactory() {}
    //@}
};

#endif
