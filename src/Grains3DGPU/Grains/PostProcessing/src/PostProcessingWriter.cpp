#include "PostProcessingWriter.hh"
#include "GrainsUtils.hh"

// ------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
PostProcessingWriter<T>::PostProcessingWriter()
{
}

// ------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
PostProcessingWriter<T>::~PostProcessingWriter()
{
}

// ------------------------------------------------------------------------------------------------
// Removes post-processing files already in the directory
template <typename T>
void PostProcessingWriter<T>::clearPostProcessingFiles(
    const std::filesystem::path& directory, const std::vector<std::regex>& patterns) const
{
    for(const auto& entry : std::filesystem::directory_iterator(directory))
    {
        if(entry.is_regular_file())
        {
            const std::string filename = entry.path().filename().string();

            bool match = false;
            for(const auto& pattern : patterns)
            {
                if(std::regex_match(filename, pattern))
                {
                    match = true;
                    break;
                }
            }

            if(match)
            {
                std::error_code ec;
                std::filesystem::remove(entry.path(), ec);
                if(ec)
                    std::cerr << "Failed to remove: " << filename << " (" << ec.message() << ")\n";
            }
        }
    }
    Gout("Post-processing files removed!");
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class PostProcessingWriter<float>;
template class PostProcessingWriter<double>;
