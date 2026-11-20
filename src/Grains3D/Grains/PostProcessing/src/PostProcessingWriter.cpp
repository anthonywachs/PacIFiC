#include "PostProcessingWriter.hh"

vector<bool> PostProcessingWriter::m_bPPWindow;


// ----------------------------------------------------------------------------
// Default constructor
PostProcessingWriter::PostProcessingWriter()
{}




// ----------------------------------------------------------------------------
// Constructor with rank and number of processes as input parameters
PostProcessingWriter::PostProcessingWriter( int const& rank_, 
	int const& nbranks_ )
  : m_rank( rank_ )
  , m_nprocs( nbranks_ )
  , m_externallyControlled( false )
{}




// ----------------------------------------------------------------------------
// Constructor with XML node, rank and number of processes as input parameters
PostProcessingWriter::PostProcessingWriter( DOMNode* dn, int const& rank_, 
	int const& nbranks_ )
  : m_rank( rank_ )
  , m_nprocs( nbranks_ )
  , m_externallyControlled( false )  
{}




// ----------------------------------------------------------------------------
// Destructor
PostProcessingWriter::~PostProcessingWriter()
{}




// ----------------------------------------------------------------------------
// Sets the initial cycle number
void PostProcessingWriter::setInitialCycleNumber( int const& cycle0, 
	double const& time )
{}




// ----------------------------------------------------------------------------
// Writes components involved in a motion or a contact error
void PostProcessingWriter::writeErreurComponents_Paraview( 
	string const& filename,
  	list<Component*> const& errcomposants )
{}
  



// ----------------------------------------------------------------------------
// Allocates post-processing windows and initializes all processes to true
void PostProcessingWriter::allocate_PostProcessingWindow( int const& nbRank )
{
  PostProcessingWriter::m_bPPWindow.resize( nbRank, true );
}




// ----------------------------------------------------------------------------
// Sets whether process of rank rank_ writes data
void PostProcessingWriter::set_PostProcessingWindow( int const& rank_,
  	bool const& bPPWindow )
{
  PostProcessingWriter::m_bPPWindow[rank_] = bPPWindow;
}




// ----------------------------------------------------------------------------
// Returns the vector of processes invovled in writing data
vector<bool> PostProcessingWriter::get_PostProcessingWindow( )
{
  return ( PostProcessingWriter::m_bPPWindow );
}




// ----------------------------------------------------------------------------
// Sets the externally controlled boolean
void PostProcessingWriter::setExternallyControlled( bool const& ext )
{
  m_externallyControlled = ext;
}




// ----------------------------------------------------------------------------
// Returns whether the post-processing writer is externally controlled
bool PostProcessingWriter::getExternallyControlled() const
{
  return ( m_externallyControlled );
} 
