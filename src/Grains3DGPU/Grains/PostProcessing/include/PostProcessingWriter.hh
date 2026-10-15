#ifndef _POSTPROCESSINGWRITER_HH_
#define _POSTPROCESSINGWRITER_HH_

#include <filesystem>
#include <regex>

#include "Basic.hh"
#include "ComponentManager.hh"
#include "GrainsMemBuffer.hh"
#include "ReaderXML.hh"

// PostProcessingWriter types
enum PostProcessingWriterType
{
    PARAVIEW,
    RAW
};

// =================================================================================================
/** @brief The class PostProcessingWriter.

    Writes results in files for post-processing by an external software.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class PostProcessingWriter
{
protected:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */

    PostProcessingWriter();
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Destructor */

    virtual ~PostProcessingWriter();
    //@}

    /** @name Get methods */
    //@{

    virtual PostProcessingWriterType getPostProcessingWriterType() const = 0;
    //@}

    /** @name Methods */
    //@{
    /** @brief Removes post-processing files already in the directory
     @param directory directory of the files to be removed
     @param patterns Regex patterns of the files to be removed */

    void clearPostProcessingFiles(const std::filesystem::path&   directory,
                                  const std::vector<std::regex>& patterns) const;

    /** @brief Initializes the post-processing writer */

    virtual void PostProcessing_start() = 0;

    /** @brief Writes post-processing data
        @param rb Arrays of rigid bodies
        @param cm Component manager
        @param currentTime Current simulation time */

    virtual void PostProcessing(const GrainsMemBuffer<RigidBody<T>*>&       rb,
                                const std::unique_ptr<ComponentManager<T>>& cm,
                                const T                                     currentTime)
        = 0;

    /** @brief Finalizes writing data */

    virtual void PostProcessing_end() = 0;
    //@}
};

#endif
