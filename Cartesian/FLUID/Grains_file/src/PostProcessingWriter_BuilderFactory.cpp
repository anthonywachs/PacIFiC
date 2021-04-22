#include "PostProcessingWriter_BuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "Paraview_PostProcessingWriter.hh"
#include "CompFeatures_PostProcessingWriter.hh"
#include "Text_PostProcessingWriter.hh"
#include "GMV_PostProcessingWriter.hh"
#include "STL_PostProcessingWriter.hh"
#include "Matlab_PostProcessingWriter.hh"
#include "DLMFD_PostProcessingWriter.hh"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation de de l'integrateur en temps
PostProcessingWriter* PostProcessingWriter_BuilderFactory::create(DOMNode* nPPW,
  	int const& rank_,int const& nbranks_)
{
  string PPWName = ReaderXML::getNodeName(nPPW);
  PostProcessingWriter* ppw = NULL;

  if ( PPWName == "Paraview" )
    ppw = new Paraview_PostProcessingWriter( nPPW, rank_, nbranks_ );
  else if( PPWName == "CompFeatures" )
    ppw = new CompFeatures_PostProcessingWriter( nPPW, rank_, nbranks_ );
  else if ( PPWName == "PositionVitesse" ||
            PPWName == "PV_Parallel" ||
            PPWName == "Text" )
    ppw = new Text_PostProcessingWriter( nPPW, rank_, nbranks_ );
  else if ( PPWName == "GMV" )
    ppw = new GMV_PostProcessingWriter( nPPW, rank_, nbranks_ );
  else if ( PPWName == "STL" )
    ppw = new STL_PostProcessingWriter( nPPW, rank_, nbranks_ );
  else if ( PPWName == "Matlab" )
    ppw = new Matlab_PostProcessingWriter( nPPW, rank_, nbranks_ );
  else if ( PPWName == "DLMFD_Serial" )
    ppw = new DLMFD_PostProcessingWriter( nPPW, rank_, nbranks_ );
  else  
  {
    cout << "Unknown post-processing writer " << PPWName << " in:" << endl;
    cout << "<PostProcessingWriters>" << endl;
    cout << "<" << PPWName << "  .../>" << endl;    
    cout << "</PostProcessingWriters>" << endl;
  }
  
  return ppw;    
}
