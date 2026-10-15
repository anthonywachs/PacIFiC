#ifndef _POSTPROCESSINGWRITERFACTORY_HH_
#define _POSTPROCESSINGWRITERFACTORY_HH_

#include "PostProcessingWriter.hh"
#include "ReaderXML.hh"
#include <string>

// =================================================================================================
/** @brief The class PostProcessingWriterFactory.

    Creates the appropriate post-processing writer depending on options.

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring
    @author A.Yazdani - 2024 - Porting to GPU */
// =================================================================================================
template <typename T>
class PostProcessingWriterFactory
{
private:
    /**@name Constructors & Destructor */
    //@{
    /** @brief Default constructor (forbidden) */

    PostProcessingWriterFactory();

    /** @brief Destructor (forbidden) */

    ~PostProcessingWriterFactory();
    //@}

public:
    /** @name Methods */
    //@{
    /** @brief Creates a post-processing writer from an XML node
        @param nPPW XML node */

    static PostProcessingWriter<T>* create(DOMNode* nPPW);
    //@}
};

#endif
