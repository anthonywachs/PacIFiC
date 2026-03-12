#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "Group3.hh"


namespace solid
{
  size_t Group3::m_sizeofGroup3 = 3 * sizeof( double );




  // --------------------------------------------------------------------------
  // Output operator
  ostream& operator << ( ostream& fileOut, Group3 const& g )
  {
    fileOut << g[X] << " " << g[Y] << " " << g[Z];
    return ( fileOut );
  }




  // --------------------------------------------------------------------------
  // Input operator
  istream& operator >> ( istream& fileIn, Group3& g )
  {
    fileIn >> g[X] >> g[Y] >> g[Z];
    return ( fileIn );
  }




  // --------------------------------------------------------------------------
  // Writes the object with a high precision format given by
  // FORMAT16DIGITS defined in GrainsExec.hh
  void Group3::writeGroup3( ostream& fileOut ) const
  {
    fileOut << GrainsExec::doubleToString( ios::scientific, FORMAT16DIGITS,
  	m_comp[X] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT16DIGITS,
  	m_comp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT16DIGITS,
  	m_comp[Z] );
  }




  // --------------------------------------------------------------------------
  // Writes the object in binary format
  void Group3::writeGroup3_binary( ostream& fileOut )
  {
    fileOut.write( reinterpret_cast<char*>( &m_comp[X] ), sizeof(double) );
    fileOut.write( reinterpret_cast<char*>( &m_comp[Y] ), sizeof(double) );
    fileOut.write( reinterpret_cast<char*>( &m_comp[Z] ), sizeof(double) );
  }




  // --------------------------------------------------------------------------
  // Reads the object in binary format
  void Group3::readGroup3_binary( istream& StreamIN )
  {
    StreamIN.read( reinterpret_cast<char*>( &m_comp[X] ), sizeof(double) );
    StreamIN.read( reinterpret_cast<char*>( &m_comp[Y] ), sizeof(double) );
    StreamIN.read( reinterpret_cast<char*>( &m_comp[Z] ), sizeof(double) );
  }
}
