#include "PostProcessingWriterFactory.hh"
#include "GrainsUtils.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "RawDataPostProcessingWriter.hh"

// -------------------------------------------------------------------------------------------------
// Creates a post-processing writer from an XML node
template <typename T>
PostProcessingWriter<T>* PostProcessingWriterFactory<T>::create(DOMNode* nPPW)
{
    std::string              PPWName = ReaderXML::getNodeName(nPPW);
    PostProcessingWriter<T>* ppw     = NULL;

    if(PPWName == "RawData")
        ppw = new RawDataPostProcessingWriter<T>(nPPW);
    else if(PPWName == "Paraview")
        ppw = new ParaviewPostProcessingWriter<T>(nPPW);
    else
        GAbort("Unknown postprocessing writer. Aborting Grains!");

    return (ppw);
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class PostProcessingWriterFactory<float>;
template class PostProcessingWriterFactory<double>;