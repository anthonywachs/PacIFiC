#include "PostProcessingWriter.hh"

vector<bool> PostProcessingWriter::m_bPPWindow;


//-----------------------------------------------------------------------------
// Constructeur par defaut
PostProcessingWriter::PostProcessingWriter()
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Will this sub-domain make outputs ?
// M.BERNARD - Fev.2015
void PostProcessingWriter::allocate_PostProcessingWindow( const int& nbRank )
{
  PostProcessingWriter::m_bPPWindow.resize( nbRank, true );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Will this sub-domain make outputs ?
// M.BERNARD - Fev.2015
void PostProcessingWriter::set_PostProcessingWindow( const int& rank, 
	const bool& bPPWindow )
{
  PostProcessingWriter::m_bPPWindow[rank] = bPPWindow;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Does this sub-domain make outputs ?
// M.BERNARD - Fev.2015
vector<bool> PostProcessingWriter::get_PostProcessingWindow( )
{
  return PostProcessingWriter::m_bPPWindow;
}
